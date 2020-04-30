import rdkit
from rdkit import DataStructs
from rdkit.Chem import AllChem, rdMolAlign
from pymol import cmd
from rdkit import Chem, RDConfig

# Take a pymol pdb object and then
'''
template = Chem.MolFromSmiles('c1nccc2n1ccc2')

from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
smiles = 'CCC'
m = Chem.MolFromSmiles(smiles)

# Read in and define the reference ligand
#ref_ligand = next(Chem.readfile('pdb', ref_input))

# Documents
#mol = AllChem.MolFromPDBFile()

fig = Draw.MolToMPL(template)
#fig = Draw.MolToMPL(smiles)
plt.title('2D Interactions')
plt.axis('off')
fig.savefig('mol.jpeg', bbox_inches='tight')
'''

def AtomsWithin(ligand, protein, cutoff):

    # load in ligand and protein
    cmd.load(ligand)
    cmd.load(protein)

    #Saves the selected atoms within x radius per atom
    selected_name = ligand + '_' + protein + '.pdb'
    cmd.save(selected_name, protein + ' within ' + cutoff + ' of ' + ligand)

    selection_nums = []
    stored.residues = []
    selected_resns = []

    cmd.iterate(selector.process(selected_name), 'stored.residues.append(resv)')
    selection_nums = stored.residues

    stored.residues = []
    # Letters

    # Selected name minus the pdb
    cmd.iterate(selector.process(selected_name), 'stored.residues.append(resn)')
    selected_resns = stored.residues

    x = 0
    Pair_List = [0]*len(selected_resns)

    while x < len(selected_resns):
        Pair_List[x] = [selection_nums[x], selected_resns[x]]
        print(selection_nums[x] + " " + selected_resns[x])


#AtomsWithin(ligand, protein, cutoff)
