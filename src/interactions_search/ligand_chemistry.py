"""Optional chemical templates; PDB atom order and coordinates remain authoritative."""
import warnings
from pathlib import Path

from rdkit import Chem


def load_reference(*, smiles=None, sdf=None):
    if smiles is not None and sdf is not None:
        raise ValueError("Usá solo una referencia: SMILES o SDF.")
    if smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
    elif sdf is not None:
        path = Path(sdf)
        if not path.is_file():
            raise ValueError(f"No existe el SDF de referencia: {path}")
        records = list(Chem.SDMolSupplier(str(path), removeHs=False))
        if len(records) != 1:
            raise ValueError("El SDF de referencia debe contener exactamente una molécula.")
        mol = records[0]
    else:
        return None
    if mol is None or mol.GetNumHeavyAtoms() == 0:
        raise ValueError("Referencia química inválida o vacía.")
    return Chem.RemoveHs(mol)


def apply_reference(mol, reference):
    """Transfer chemistry through a full heavy-atom graph match, retaining explicit H.

    Equivalent graph matches prefer compatible explicit H and existing multiple
    bonds/charges in the PDB, with first-match tie breaking. No reference
    conformer is used. A different heavy-atom graph or incompatible explicit H
    raises ValueError instead of silently falling back to PDB inference.
    """
    heavy_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    heavy = Chem.RWMol(mol)
    for idx in reversed(range(mol.GetNumAtoms())):
        if mol.GetAtomWithIdx(idx).GetAtomicNum() == 1:
            heavy.RemoveAtom(idx)
    heavy = heavy.GetMol()
    if (heavy.GetNumAtoms() != reference.GetNumAtoms()
            or heavy.GetNumBonds() != reference.GetNumBonds()):
        raise ValueError("La referencia no coincide con el grafo de átomos pesados del PDB.")
    try:
        def skeleton(source):
            graph = Chem.Mol(source)
            for atom in graph.GetAtoms():
                atom.SetFormalCharge(0)
                atom.SetIsAromatic(False)
                atom.SetNumExplicitHs(0)
                atom.SetNoImplicit(False)
            for bond in graph.GetBonds():
                bond.SetBondType(Chem.BondType.SINGLE)
                bond.SetIsAromatic(False)
            graph.UpdatePropertyCache(strict=False)
            return graph

        matches = skeleton(heavy).GetSubstructMatches(skeleton(reference), uniquify=False,
                                                     maxMatches=1001)
        if not matches:
            raise ValueError("El grafo de la referencia no coincide con el del PDB.")
        if len(matches) > 1000:
            raise ValueError("Demasiadas correspondencias del grafo (>1000); referencia ambigua.")

        def compatibility(match):
            consistent_h = True
            evidence = 0
            for atom in reference.GetAtoms():
                target = mol.GetAtomWithIdx(heavy_indices[match[atom.GetIdx()]])
                explicit_h = sum(n.GetAtomicNum() == 1 for n in target.GetNeighbors())
                consistent_h &= explicit_h <= atom.GetTotalNumHs()
                if target.GetFormalCharge() != 0:
                    evidence += target.GetFormalCharge() == atom.GetFormalCharge()
            for bond in reference.GetBonds():
                target = heavy.GetBondBetweenAtoms(match[bond.GetBeginAtomIdx()],
                                                   match[bond.GetEndAtomIdx()])
                if target.GetBondTypeAsDouble() > 1:
                    evidence += target.GetBondType() == bond.GetBondType()
            return consistent_h, evidence

        best_match = max(matches, key=compatibility)
        if len(matches) > 1:
            warnings.warn('La referencia tiene varias correspondencias con el PDB; '
                          'se prioriza compatibilidad con H, enlaces múltiples y cargas '
                          'del PDB (primera en empates). Revisá grupos simétricos/protonación.',
                          UserWarning, stacklevel=2)
        mapping = [heavy_indices[i] for i in best_match]
        result = Chem.Mol(mol)
        for src in reference.GetAtoms():
            dst = result.GetAtomWithIdx(mapping[src.GetIdx()])
            explicit_h = sum(n.GetAtomicNum() == 1 for n in dst.GetNeighbors())
            total_h = src.GetTotalNumHs()
            if explicit_h > total_h:
                raise ValueError("Los hidrógenos explícitos del PDB contradicen la referencia.")
            dst.SetFormalCharge(src.GetFormalCharge())
            dst.SetIsAromatic(src.GetIsAromatic())
            dst.SetNumExplicitHs(total_h - explicit_h)
            dst.SetNoImplicit(True)
        for src in reference.GetBonds():
            dst = result.GetBondBetweenAtoms(mapping[src.GetBeginAtomIdx()],
                                            mapping[src.GetEndAtomIdx()])
            dst.SetBondType(src.GetBondType())
            dst.SetIsAromatic(src.GetIsAromatic())
        Chem.SanitizeMol(result)
    except (ValueError, RuntimeError) as exc:
        raise ValueError(f"No se pudo aplicar la referencia química al PDB: {exc}") from exc
    return result
