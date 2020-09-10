from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D


def embed_confs(mol, seed=10, numConfs=27):
    ''' Embeds multiple conformations to a molecule,
    returns a list of conformation IDs. If the molecule
    is given as a SMILES string, rdkit.mol object is
    created and hydrated.
    '''
    if type(mol) is str:
        mol = Chem.AddHs(Chem.MolFromSmiles(mol), addCoords=True)

    conformers = AllChem.EmbedMultipleConfs(mol, pruneRmsThresh=0.5, numConfs=numConfs, randomSeed=seed)
    return list(conformers)


def mmff_opt_confs(mol, seed=10, numConfs=27):
    ''' Performes MMFF'94 optiomization of embed conformations,
    returns convergence and absolute energy value.
    '''
    embed_confs(mol, seed, numConfs)
    energies = Chem.rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol)
    return mol, energies


def converged_conformers(mol, energies, cutoff=10):
    '''Purges unconverged FF-optimized conformers,
    returns mol and a table with conformer ID and its
     relative minimized energy.
    (rdchem.Mol, [(confID, energy),...])
    '''
    conf_n_energies = list()
    for i in range(len(energies)):
        if energies[i][0] == 0 and (energies[i][1] - min(energies)[1]) <= cutoff:
            conf_n_energies.append((i, (energies[i][1] - min(energies)[1])))
        else:
            mol.RemoveConformer(i)
    return mol, conf_n_energies


def average_PMIscore(mol, cutoff=10):
    if type(mol) is str:
        mol = Chem.AddHs(Chem.MolFromSmiles(mol))
    mol = mmff_opt_confs(mol)
    mol = converged_conformers(mol[0], mol[1], cutoff)
    NPR1 = list()
    NPR2 = list()
    for i in range(mol[0].GetNumConformers()):
        mol[0].ClearComputedProps() # RDKit version 2017.09.1, not needed in newer versions
        NPR1.append(Descriptors3D.NPR1(mol[0], mol[1][i][0]))
        NPR2.append(Descriptors3D.NPR2(mol[0], mol[1][i][0]))
    return (sum(NPR1)+sum(NPR2))/len(NPR1)


if __name__ == '__main__':
    import sys
    print(average_PMIscore(sys.argv[1]))
