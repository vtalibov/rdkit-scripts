import pandas as pd
from rdkit import Chem


def txt_to_list(file):
    '''Reads .txt file, returns a list of strings.'''
    list_out = []
    for line in open(file):
        if line[-1] == '\n':
            list_out.append(line.replace('\n',""))
        else:
            list_out.append(line)
    return list_out


def prepare_smcm_en_values(file):
    '''Reads multiple files with SMART patterns
    and corresponding SMCM scores, returns pandas
    dataframe with rdkit.Mol objects and SMCM scores.'''
    smcm_scores = pd.DataFrame()
    for fin in txt_to_list(file):
        smcm_scores = pd.concat([smcm_scores, pd.read_csv(fin, header=None)],
                                ignore_index=True)
    smcm_scores[0] = [Chem.MolFromSmarts(entry) for entry in smcm_scores[0]]
    return smcm_scores


def en_score(mol, file='SMCM_values.txt'):
    '''Calculates elecronegativity-based part of SMCM'''
    smcm_scores = prepare_smcm_en_values(file)
    score = 0
    for i in range(len(smcm_scores)):
        if mol.HasSubstructMatch(smcm_scores[0][i]):
            score = score + len(mol.GetSubstructMatches(smcm_scores[0][i])) * \
                    smcm_scores[1][i]
    return score


def ring_atoms(mol):
    '''Returns sets of ring system atoms.'''
    rings = mol.GetRingInfo()
    return [set(a) for a in rings.AtomRings()]


def rings_shared_atoms(mol):
    '''Returns sets of atoms, shared between two ring system.'''
    shared_atoms = []
    atoms = ring_atoms(mol)
    for ring1 in atoms:
        for ring2 in atoms:
            shared = ring1.intersection(ring2)
            if ring1 != ring2 and len(shared) != 0:
                if shared not in shared_atoms:
                    shared_atoms.append(shared)
        del atoms[atoms.index(ring1)]
    return shared_atoms

# Two ring_score functions. The first (commented) calculates the score
# additively. Not sure if it is correct way.

# def ring_score(mol):
#     '''Perfomns ring assignment and returns
#     ring systems-dependent part of SMCM.'''
#     score = 0
#     rings = ring_atoms(mol)
#     for ring in rings:
#         two = [3, 4, 7, 8, 9]
#         if len(ring) in two:
#             score += 2
#         else:
#             score += 1
#     shared_atoms = rings_shared_atoms(mol)
#     for shared in shared_atoms:
#         if len(shared) == 1:    # spiro
#             score += 3
#         elif len(shared) == 2:  # fused
#             score += 2
#         elif len(shared) > 2:   # bridge head
#             score += 4
#     return score


def ring_score(mol):
    '''Perfomns ring assignment and returns
    ring systems-dependent part of SMCM.'''
    rings = ring_atoms(mol)
    shared_atoms = rings_shared_atoms(mol)
    simple_rings = [5, 6]   # penta-and hexaatomic ring systems
    score = 0
    for ring in rings:
        junction_coef = 1
        if len(ring) in simple_rings:
            ring_coef = 1
        else:
            ring_coef = 2
        # only the first instance of junction is significant.
        for shared in shared_atoms:
            if shared.issubset(ring):
                if len(shared) == 1:    # Spiro
                    junction_coef = 3
                elif len(shared) == 2:
                    junction_coef = 2   # Fused
                elif len(shared) > 2:
                    junction_coef = 4   # Bridged
        score += ring_coef * junction_coef
    return score


def smcm_score(mol, file='SMCM_values.txt'):
    '''Takes RDKit mol object and returns Synthethic
    Molecular Complexity Metric (SMCM) descriptor,
    similar to described by Allu and Oprea in
    J. Chem. Inf. Model, 2005, 45, 1237-1243.'''
    return en_score(mol, file)+ring_score(mol)
