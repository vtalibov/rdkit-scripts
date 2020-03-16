from rdkit import Chem

def all_bond_rings(mol):
    ring_info = mol.GetRingInfo()
    rings = [set(r) for r in ring_info.BondRings()]
    return rings
