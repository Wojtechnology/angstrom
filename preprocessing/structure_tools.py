"""Structure utilities for the Runs N' Poses memorisation visualiser.

Everything here works on RDKit molecules (ligands) and gemmi structures (receptors).
Predicted complexes are superposed onto the ground-truth binding site so that every
method's pose can be overlaid in one frame in the viewer.
"""
from __future__ import annotations

import io
import math
from dataclasses import dataclass, field

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolTransforms
from rdkit.Chem import rdDetermineBonds  # noqa: F401  (import registers functionality)

VDW = {"H": 1.1, "C": 1.7, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.8, "S": 1.8, "Cl": 1.75,
       "Br": 1.85, "I": 1.98, "B": 1.92, "Se": 1.9, "Si": 2.1}
CLASH_RATIO = 0.75  # same threshold PoseBusters uses for minimum protein distance
POCKET_CUTOFF = 8.0  # Å from any ligand heavy atom -> pocket residue (for superposition)
CONTACT_CUTOFF = 4.5  # Å for the residue list shown in the viewer


# --------------------------------------------------------------------------- receptors
@dataclass
class Receptor:
    """Protein heavy atoms of one complex plus per-chain CA tables."""
    structure: gemmi.Structure
    chains: dict = field(default_factory=dict)  # chain name -> {"seq": str, "ca": (n,3), "resnums": [..]}

    @classmethod
    def from_file(cls, path: str, polymer_only: bool = True) -> "Receptor":
        st = gemmi.read_structure(str(path))
        st.setup_entities()
        st.remove_hydrogens()
        st.remove_waters()
        rec = cls(structure=st)
        model = st[0]
        for ch in model:
            seq, ca, resnums = [], [], []
            for res in ch:
                if res.het_flag == "H" and polymer_only:
                    continue
                info = gemmi.find_tabulated_residue(res.name)
                if info is None or not info.is_amino_acid():
                    continue
                atom = res.find_atom("CA", "*")
                if atom is None:
                    continue
                seq.append(info.one_letter_code.upper() or "X")
                ca.append([atom.pos.x, atom.pos.y, atom.pos.z])
                resnums.append(res.seqid.num)
            if seq:
                rec.chains[ch.name] = {"seq": "".join(seq), "ca": np.array(ca), "resnums": resnums}
        return rec

    def heavy_atoms(self):
        """Yield (element, xyz, chain, resname, resnum, atom_name) for protein heavy atoms."""
        for ch in self.structure[0]:
            for res in ch:
                info = gemmi.find_tabulated_residue(res.name)
                if info is None or not info.is_amino_acid():
                    continue
                for at in res:
                    if at.element.name == "H":
                        continue
                    yield at.element.name, np.array([at.pos.x, at.pos.y, at.pos.z]), ch.name, res.name, res.seqid.num, at.name

    def transform(self, R: np.ndarray, t: np.ndarray):
        for ch in self.structure[0]:
            for res in ch:
                for at in res:
                    p = np.array([at.pos.x, at.pos.y, at.pos.z]) @ R.T + t
                    at.pos = gemmi.Position(*p)
        for c in self.chains.values():
            c["ca"] = c["ca"] @ R.T + t

    def to_pdb(self, protein_only: bool = True) -> str:
        st = self.structure.clone()
        if protein_only:
            for ch in st[0]:
                for i in range(len(ch) - 1, -1, -1):
                    info = gemmi.find_tabulated_residue(ch[i].name)
                    if info is None or not info.is_amino_acid():
                        del ch[i]
        st.remove_empty_chains()
        st.shorten_chain_names()
        return st.make_pdb_string()


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Return R, t minimising |P@R.T + t - Q|."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - pc @ R.T


def align_seqs(a: str, b: str):
    """Global sequence alignment via gemmi; returns list of (i, j) matched index pairs."""
    res = gemmi.align_string_sequences(list(a), list(b), [])
    pairs = []
    i = j = 0
    for op in res.cigar_str().replace("M", "M ").replace("I", "I ").replace("D", "D ").split():
        n, kind = int(op[:-1]), op[-1]
        for _ in range(n):
            if kind == "M":
                pairs.append((i, j)); i += 1; j += 1
            elif kind == "I":
                i += 1
            else:
                j += 1
    return pairs, res.identity() if hasattr(res, "identity") else None


def pocket_residue_indices(rec: Receptor, lig_xyz: np.ndarray, cutoff: float):
    """Per chain, indices of CA table entries whose residue has any heavy atom within cutoff."""
    out = {}
    for ch in rec.structure[0]:
        if ch.name not in rec.chains:
            continue
        resnums = rec.chains[ch.name]["resnums"]
        idx_of = {r: k for k, r in enumerate(resnums)}
        hits = set()
        for res in ch:
            if res.seqid.num not in idx_of:
                continue
            for at in res:
                p = np.array([at.pos.x, at.pos.y, at.pos.z])
                if np.min(np.linalg.norm(lig_xyz - p, axis=1)) <= cutoff:
                    hits.add(idx_of[res.seqid.num]); break
        out[ch.name] = sorted(hits)
    return out


def superpose_on_pocket(pred: Receptor, gt: Receptor, gt_lig_xyz: np.ndarray):
    """Superpose predicted receptor onto the GT binding site (CA atoms within POCKET_CUTOFF).

    Chains are matched by sequence alignment. Returns dict with pocket/global CA RMSD and the
    chain mapping; the predicted receptor is transformed in place.
    """
    pocket = pocket_residue_indices(gt, gt_lig_xyz, POCKET_CUTOFF)
    P, Q, Pall, Qall = [], [], [], []
    mapping = {}
    used = set()
    for gname, g in gt.chains.items():
        best = None
        for pname, p in pred.chains.items():
            if pname in used:
                continue
            pairs, _ = align_seqs(g["seq"], p["seq"])
            ident = sum(g["seq"][i] == p["seq"][j] for i, j in pairs) / max(1, len(g["seq"]))
            if best is None or ident > best[0]:
                best = (ident, pname, pairs)
        if best is None or best[0] < 0.3:
            continue
        ident, pname, pairs = best
        used.add(pname)
        mapping[gname] = pname
        pk = set(pocket.get(gname, []))
        for i, j in pairs:
            Qall.append(g["ca"][i]); Pall.append(pred.chains[pname]["ca"][j])
            if i in pk:
                Q.append(g["ca"][i]); P.append(pred.chains[pname]["ca"][j])
    if len(P) < 3:
        P, Q = Pall, Qall
    if len(P) < 3:
        raise RuntimeError("could not match enough CA atoms for superposition")
    P, Q = np.array(P), np.array(Q)
    R, t = kabsch(P, Q)
    pred.transform(R, t)
    pocket_rmsd = float(np.sqrt(np.mean(np.sum((P @ R.T + t - Q) ** 2, axis=1))))
    Pall, Qall = np.array(Pall), np.array(Qall)
    global_rmsd = float(np.sqrt(np.mean(np.sum((Pall @ R.T + t - Qall) ** 2, axis=1))))
    return {"pocket_ca_rmsd": pocket_rmsd, "global_ca_rmsd_pocket_aligned": global_rmsd,
            "n_pocket_ca": int(len(P)), "chain_map": mapping, "R": R, "t": t}


# --------------------------------------------------------------------------- ligands
def ligand_from_structure(path: str, chain_name: str, template: Chem.Mol) -> Chem.Mol:
    """Extract ligand heavy atoms of one chain from a predicted CIF/PDB and assign bond orders
    from the ground-truth template molecule (same chemical species)."""
    st = gemmi.read_structure(str(path))
    st.remove_hydrogens()
    model = st[0]
    xyz, elems, names = [], [], []
    for ch in model:
        if ch.name != chain_name:
            continue
        for res in ch:
            for at in res:
                if at.element.name in ("H", "D"):
                    continue
                xyz.append([at.pos.x, at.pos.y, at.pos.z])
                elems.append(at.element.name)
                names.append(at.name)
    if not xyz:
        raise RuntimeError(f"no ligand atoms in chain {chain_name}")
    xyz = np.array(xyz)
    # Build a bare molecule with connectivity from distances, then take bond orders from template.
    raw = Chem.RWMol()
    conf = Chem.Conformer(len(xyz))
    for k, (e, p) in enumerate(zip(elems, xyz)):
        a = Chem.Atom(e.capitalize() if len(e) > 1 else e)
        a.SetNoImplicit(True)
        raw.AddAtom(a)
        conf.SetAtomPosition(k, p.tolist())
    raw.AddConformer(conf)
    mol = raw.GetMol()
    rdDetermineBonds.DetermineConnectivity(mol, useVdw=True)
    for a in mol.GetAtoms():
        a.SetNoImplicit(False)
    tmpl = Chem.RemoveHs(template)
    if tmpl.GetNumAtoms() != mol.GetNumAtoms():
        raise RuntimeError(f"atom count mismatch template={tmpl.GetNumAtoms()} model={mol.GetNumAtoms()}")
    mol = AllChem.AssignBondOrdersFromTemplate(tmpl, mol)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    return mol


def transform_mol(mol: Chem.Mol, R: np.ndarray, t: np.ndarray):
    conf = mol.GetConformer()
    pos = conf.GetPositions() @ R.T + t
    for i, p in enumerate(pos):
        conf.SetAtomPosition(i, p.tolist())


def symmetry_rmsd(pred: Chem.Mol, ref: Chem.Mol) -> float:
    """Symmetry-corrected heavy-atom RMSD without re-alignment."""
    return float(rdMolAlign.CalcRMS(Chem.RemoveHs(pred), Chem.RemoveHs(ref)))


# --------------------------------------------------------------------------- diagnostics
def _bonded_distance_matrix(mol: Chem.Mol):
    return Chem.GetDistanceMatrix(mol)


def ligand_diagnostics(pred: Chem.Mol, ref: Chem.Mol, protein_atoms: list):
    """Atom-level localisation of geometry problems in a predicted ligand pose.

    Bond lengths/angles are compared with the ground-truth pose of the very same molecule,
    which is an experimentally refined geometry, so deviations are strain the model invented.
    """
    pred = Chem.RemoveHs(pred); ref = Chem.RemoveHs(ref)
    match = ref.GetSubstructMatch(pred, useChirality=False)
    if len(match) != pred.GetNumAtoms():
        match = tuple(range(pred.GetNumAtoms()))
    pc, rc = pred.GetConformer(), ref.GetConformer()
    ppos = pc.GetPositions()

    bonds = []
    for b in pred.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        lp = rdMolTransforms.GetBondLength(pc, i, j)
        lr = rdMolTransforms.GetBondLength(rc, match[i], match[j])
        ratio = lp / lr if lr > 0 else 1.0
        bonds.append({"atoms": [i, j], "pred": round(lp, 3), "ref": round(lr, 3), "ratio": round(ratio, 3),
                      "flag": bool(ratio < 0.85 or ratio > 1.15)})

    angles = []
    for a in pred.GetAtoms():
        nbrs = [n.GetIdx() for n in a.GetNeighbors()]
        for x in range(len(nbrs)):
            for y in range(x + 1, len(nbrs)):
                i, j, k = nbrs[x], a.GetIdx(), nbrs[y]
                ap = rdMolTransforms.GetAngleDeg(pc, i, j, k)
                ar = rdMolTransforms.GetAngleDeg(rc, match[i], match[j], match[k])
                dev = ap - ar
                angles.append({"atoms": [i, j, k], "pred": round(ap, 1), "ref": round(ar, 1),
                               "dev": round(dev, 1), "flag": bool(abs(dev) > 20)})

    # intra-ligand clashes: atoms >= 4 bonds apart closer than CLASH_RATIO * vdW sum
    topo = _bonded_distance_matrix(pred)
    intra = []
    n = pred.GetNumAtoms()
    for i in range(n):
        for j in range(i + 1, n):
            if topo[i, j] < 4:
                continue
            d = float(np.linalg.norm(ppos[i] - ppos[j]))
            lim = CLASH_RATIO * (VDW.get(pred.GetAtomWithIdx(i).GetSymbol(), 1.7) + VDW.get(pred.GetAtomWithIdx(j).GetSymbol(), 1.7))
            if d < lim:
                intra.append({"atoms": [i, j], "dist": round(d, 2), "limit": round(lim, 2)})

    # protein-ligand clashes and contacts
    clashes, contacts = [], {}
    if protein_atoms:
        P = np.array([p[1] for p in protein_atoms])
        for i in range(n):
            d = np.linalg.norm(P - ppos[i], axis=1)
            near = np.where(d < CONTACT_CUTOFF)[0]
            ri = VDW.get(pred.GetAtomWithIdx(i).GetSymbol(), 1.7)
            for k in near:
                el, _, ch, resn, resi, an = protein_atoms[k]
                key = f"{ch}:{resn}{resi}"
                contacts[key] = min(contacts.get(key, 99.0), float(d[k]))
                lim = CLASH_RATIO * (ri + VDW.get(el.capitalize() if len(el) > 1 else el, 1.7))
                if d[k] < lim:
                    clashes.append({"atom": i, "protein": {"chain": ch, "resname": resn, "resnum": int(resi), "atom": an},
                                    "dist": round(float(d[k]), 2), "limit": round(lim, 2)})

    # stereo: compare CIP labels of tetrahedral centres with the ground truth
    Chem.AssignStereochemistryFrom3D(pred); Chem.AssignStereochemistryFrom3D(ref)
    stereo = []
    for i in range(n):
        a, r = pred.GetAtomWithIdx(i), ref.GetAtomWithIdx(match[i])
        if a.GetAtomicNum() != 6:  # PoseBusters only judges carbon centres; P/S labels are ambiguous
            continue
        cp, cr = a.GetPropsAsDict().get("_CIPCode"), r.GetPropsAsDict().get("_CIPCode")
        if cp or cr:
            stereo.append({"atom": i, "pred": cp, "ref": cr, "flag": bool(cp != cr)})

    # aromatic ring flatness
    rings = []
    ri = pred.GetRingInfo()
    for ring in ri.AtomRings():
        if not all(pred.GetAtomWithIdx(k).GetIsAromatic() for k in ring):
            continue
        pts = ppos[list(ring)]
        c = pts.mean(0)
        _, _, vt = np.linalg.svd(pts - c)
        dev = float(np.max(np.abs((pts - c) @ vt[2])))
        rings.append({"atoms": list(ring), "max_dev": round(dev, 3), "flag": bool(dev > 0.25)})

    contact_list = [{"residue": k, "min_dist": round(v, 2)} for k, v in sorted(contacts.items(), key=lambda kv: kv[1])]
    flagged = set()
    for b in bonds:
        if b["flag"]:
            flagged.update(b["atoms"])
    for a in angles:
        if a["flag"]:
            flagged.add(a["atoms"][1])
    for c in intra:
        flagged.update(c["atoms"])
    for c in clashes:
        flagged.add(c["atom"])
    for s in stereo:
        if s["flag"]:
            flagged.add(s["atom"])
    for r in rings:
        if r["flag"]:
            flagged.update(r["atoms"])
    return {"bonds": bonds, "angles": [a for a in angles if a["flag"]], "n_angles": len(angles),
            "intra_clashes": intra, "protein_clashes": clashes, "stereo": stereo, "rings": rings,
            "contacts": contact_list, "flagged_atoms": sorted(flagged),
            "summary": {"bad_bonds": sum(b["flag"] for b in bonds), "bad_angles": sum(a["flag"] for a in angles),
                        "intra_clashes": len(intra), "protein_clashes": len(clashes),
                        "stereo_mismatches": sum(s["flag"] for s in stereo),
                        "nonplanar_rings": sum(r["flag"] for r in rings)}}


# --------------------------------------------------------------------------- minimisation
def _mmff(mol):
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    if props is None:
        raise RuntimeError("MMFF parameters unavailable")
    return props


def _energy(mol, props, pos=None):
    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=0)
    return float(ff.CalcEnergy(pos.ravel().tolist()) if pos is not None else ff.CalcEnergy())


def minimise_with_trajectory(mol: Chem.Mol, max_frames: int = 30, its_per_frame: int = 15, n_global: int = 30, seed: int = 42):
    """MMFF94s minimisation of a ligand pose recording a trajectory.

    Returns (summary dict, list of heavy-atom coordinate frames, mol with hydrogens at the local minimum).
    Strain is reported relative to (a) the pose after a restrained relaxation of hydrogens and tiny
    heavy-atom deviations and (b) the nearest local minimum and (c) the global minimum from a
    conformer search. Displacement per heavy atom localises where the pose had to move.
    """
    mh = Chem.AddHs(mol, addCoords=True)
    props = _mmff(mh)
    heavy = [a.GetIdx() for a in mh.GetAtoms() if a.GetAtomicNum() > 1]
    start = mh.GetConformer().GetPositions()

    # (a) relax hydrogens with heavy atoms restrained close to the predicted positions
    ff = AllChem.MMFFGetMoleculeForceField(mh, props)
    for i in heavy:
        ff.MMFFAddPositionConstraint(i, 0.05, 200.0)
    ff.Minimize(maxIts=500)
    e_pose = _energy(mh, props)
    pos_pose = mh.GetConformer().GetPositions()

    # (b) free local minimisation, recording frames
    ff = AllChem.MMFFGetMoleculeForceField(mh, props)
    frames = [pos_pose[heavy].copy()]
    energies = [e_pose]
    converged = False
    for _ in range(max_frames - 1):
        rc = ff.Minimize(maxIts=its_per_frame)
        pos = np.array(ff.Positions()).reshape(-1, 3)
        frames.append(pos[heavy].copy())
        energies.append(float(ff.CalcEnergy()))
        if rc == 0:
            converged = True
            break
    ff.Minimize(maxIts=2000)  # finish to the local minimum
    pos_local = np.array(ff.Positions()).reshape(-1, 3)
    if not converged:
        frames.append(pos_local[heavy].copy()); energies.append(float(ff.CalcEnergy()))
    conf = mh.GetConformer()
    for i, p in enumerate(pos_local):
        conf.SetAtomPosition(i, p.tolist())
    e_local = float(ff.CalcEnergy())

    # (c) global minimum from a conformer search
    mg = Chem.Mol(mh)
    mg.RemoveAllConformers()
    cids = AllChem.EmbedMultipleConfs(mg, numConfs=n_global, randomSeed=seed)
    e_global = e_local
    if len(cids):
        res = AllChem.MMFFOptimizeMoleculeConfs(mg, mmffVariant="MMFF94s", maxIters=2000)
        e_global = min(min(e for _, e in res), e_local)

    disp = np.linalg.norm(pos_local[heavy] - start[heavy], axis=1)
    rmsd_drift = float(np.sqrt(np.mean(disp ** 2)))
    summary = {
        "e_pose": round(e_pose, 2), "e_local": round(e_local, 2), "e_global": round(e_global, 2),
        "strain_local": round(e_pose - e_local, 2), "strain_global": round(e_pose - e_global, 2),
        "rmsd_drift": round(rmsd_drift, 3), "max_atom_displacement": round(float(disp.max()), 3),
        "atom_displacement": [round(float(d), 3) for d in disp],
        "n_frames": len(frames), "converged": converged, "energies": [round(e, 2) for e in energies],
        "frame_rmsd": [round(float(np.sqrt(np.mean(np.sum((f - frames[0]) ** 2, axis=1)))), 3) for f in frames],
    }
    return summary, frames, mh


def frames_to_sdf(mol_heavy: Chem.Mol, frames: list, energies: list) -> str:
    """Multi-record SDF with one heavy-atom conformer per minimisation frame."""
    m = Chem.RemoveHs(mol_heavy)
    buf = io.StringIO()
    w = Chem.SDWriter(buf)
    for k, (f, e) in enumerate(zip(frames, energies)):
        conf = m.GetConformer()
        for i, p in enumerate(f):
            conf.SetAtomPosition(i, p.tolist())
        m.SetProp("_Name", f"frame {k}")
        m.SetProp("energy", str(e))
        w.write(m)
    w.close()
    return buf.getvalue()


def mol_to_sdf(mol: Chem.Mol, name: str = "") -> str:
    m = Chem.Mol(mol)
    if name:
        m.SetProp("_Name", name)
    return Chem.MolToMolBlock(m) + "$$$$\n"
