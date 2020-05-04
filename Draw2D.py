import rdkit
from rdkit import DataStructs
from rdkit.Chem import AllChem, rdMolAlign
#from pymol import cmd
from rdkit import Chem, RDConfig

from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve, S, solveset

# Gets the max for use in the residue mapping
def GetMaxCoords(pos):
    max_val = 0
    max_c = 0.5

    # Cycle through to find the max
    for atom in pos:
        x = atom[0]
        y = atom[1]

        if abs(x) > max_val:
            max_val = x

        elif abs(y) > max_val:
            max_val = y

        else:
            pass

    if max_val < 0:
        e_val = 0.50
    else:
        e_val = 0.25

    max_c = ((1 / max_val) + e_val) * 0.98
    print("max: " + str(max_c))
    return max_c
##############################################################################################################
# Get the x and y positions

def GetAtomPos(pos, AtomNum):

    CurrNum = 1
    Coords = []

    # Cycle through to find the correct atomic number
    for atom in pos:
        if CurrNum == AtomNum:
            Coords = atom
            CurrNum += 1
        else:
            CurrNum += 1
    return Coords[0], Coords[1]
##############################################################################################################

# Gets the max for use in the residue mapping
def GetResCoords(pos, AtomNum, max_val):

    CurrNum = 1
    Coords = []

    # Cycle through to find the correct atomic number
    for atom in pos:
        if CurrNum == AtomNum:
            Coords = atom
            CurrNum += 1
        else:
            CurrNum += 1

    # Variables
    if Coords[0] < 0:
        x2 = abs(Coords[0]) + 0.25
    else:
        x2 = abs(Coords[0]) + 0.50


    if Coords[1] < 0:
        y2 = abs(Coords[1] + 0.25)
    else:
        y2 = abs(Coords[1])

    r = max_val

    # y1 = mx1+b
    print(x2)
    print(y2)
    x1, y1 = symbols('x1 y1')
    eq1 = Eq(x1**2 + y1**2 - r**2)
    eq2 = Eq((y1-y2)/(x1-x2)*x2 + 0.4)


    x1, y1 = solve((eq1, eq2), x1, y1)

    print("x1: " + str(x1) + " y1: " + str(y1))
    # r^2 = x2^2 + y2^21
    return x1,y1



##############################################################################################################

# Finds the residues within the atoms of the ligand
# Out of those residues, a list is then created [[Res, Num, Chain], corresponding ligand atom]
def AtomsWithin(ligand, protein, cutoff):

    # load in ligand and protein
    cmd.load(ligand)
    cmd.load(protein)

    #Saves the selected atoms within x radius of the ligand
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

    # Makes a pair list of resns and their corresponding numbers
    while x < len(selected_resns):
        Pair_List[x] = [selection_nums[x], selected_resns[x]]
        print(selection_nums[x] + " " + selected_resns[x])


def MoreProcessing(ResList):

    # For each of the residues in the list, go through and find the atoms that correspond with the resns, then take
    # an average of the coordinates (x and y)
    residue_number = 0
    for residue in ResList:
        # Residue is [Resn, number]


        ligand_atoms = []
        cmd.select('obj', ligand + ' within ' + cutoff + ' of ' + protein +' and resn ' + residue[0] + ' and resv ' + residue[1])
        stored.residues = []
        cmd.iterate(selector.process('obj'), 'stored.residues.append(resv)')
        ligand_atoms = stored.residues

        x_tally = 0
        y_tally = 0
        total_num = 0
        for atom_num in ligand_atoms:
            x, y = GetAtomPos(moleule, atom_num)
            x_tally += x
            y_tally += y
            total_num += 1

        # Checks what the corresponding x, y values are, and adds them to a list with the resn string and num

        ResList[residue_number].append([x_tally/total_num, y_tally/total_num])
        residue_number += 1




mol = Chem.MolFromSmiles('c1nccc2n1ccc2')

print("get res coordinates")
AllChem.EmbedMolecule(mol, AllChem.ETKDG())

# Whole position array
pos = mol.GetConformer(0).GetPositions()

# Gets the maximum coordinate number (adjusted)
max_num = GetMaxCoords(pos)

# Gets the x and y values of a residue with
x1, y1 = GetResCoords(pos, 3, max_num)

#GetResCoords(mol, 2, max_val)

#AtomsWithin(ligand, protein, cutoff)

#def MainFunction():
m = Chem.MolFromSmiles('c1nccc2n1ccc2')

# AllChem.Compute2DCoords(m)
AllChem.EmbedMolecule(m, AllChem.ETKDG())
fig = Draw.MolToMPL(m)
pos = m.GetConformer(0).GetPositions()

print(pos)

# Turn off after
plt.grid(True)

bbox_props = dict(boxstyle="circle,pad=0.3", ec="b", lw=2)
# max_c = ((1/2.05) + 0.25)*0.98
import numpy as np

# Add for each of the residues that needs to happen
x1_graph = (complex(x1[0]).real)
y1_graph = (complex(y1[0]).imag)

text_add = plt.text(x1_graph , y1_graph, "Met", size=15, bbox=bbox_props)
coord_min = 0
coord_max = 0

fig.savefig('mol.jpeg', bbox_inches='tight', va='center', ha='center')

# Prints the block
print(Chem.MolToMolBlock(m))



