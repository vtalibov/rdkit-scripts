from rdkit import Chem
import pandas as pd


def txt_to_list(file):
    '''Reads .txt file, returns a list of strings.'''
    list_out = []
    for line in open(file):
        if line[-1] == '\n':
            list_out.append(line.replace('\n',""))
        else:
            list_out.append(line)
    return list_out
 
 
def ring_atoms(mol):
    rings = mol.GetRingInfo()
    return [set(a) for a in rings.AtomRings()]
 
 
def rings_shared_atoms(mol):
    shared_atoms = []
    atoms = ring_atoms(mol)
    for set1 in atoms:
        for set2 in atoms:
            intersect = set1.intersection(set2)
            if set1 != set2 and len(intersect) != 0:
                if intersect not in shared_atoms:
                    shared_atoms.append(intersect)
    return shared_atoms


def prepare_SMCM_EN_values(file):
    '''Reads multiple files with SMART patterns and corresponding SMCM scores, returns pandas dataframe with rdkit.Mol objects and SMCM scores.'''
    SMCM_scores = pd.DataFrame()
    for f in txt_to_list(file):
        SMCM_scores = pd.concat([SMCM_scores, pd.read_csv(f, header=None)], ignore_index = True)
    SMCM_scores[0] = [Chem.MolFromSmarts(entry) for entry in SMCM_scores[0]]
    return SMCM_scores
 
 
def EN_score(mol, file = 'SMCM_values.txt'):
    '''Calculates elecronegativity-based part of SMCM'''
    SMCM_scores = prepare_SMCM_EN_values(file)
    score = 0
    for i in range(len(SMCM_scores)):
        if mol.HasSubstructMatch(SMCM_scores[0][i]):
            score = score + len(mol.GetSubstructMatches(SMCM_scores[0][i]))*SMCM_scores[1][i]
    return score
 
 
def ring_score(mol):
    score = 0
    rings = ring_atoms(mol)
    for ring in rings:
        two = [3,4,7,8,9]
        if len(ring) in two:
            score = score + 2
        else:
            score = score + 1
    
    shared_atoms = rings_shared_atoms(mol)
    for set in shared_atoms:
        if len(set) == 1:
            score = score + 3
        elif len(set) == 2:
            score = score + 2
        elif len(set) > 2:
            score = score + 4
            
    return score
 
 
def SMCM_score(mol, file='SMCM_values.txt'):
    return EN_score(mol, file)+ring_score(mol)
