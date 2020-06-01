import rdkit
from rdkit import DataStructs
from rdkit.Chem import AllChem, rdMolAlign
#from pymol import cmd
from rdkit import Chem, RDConfig

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw.mplCanvas import Canvas

import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve, S, solveset
global posi, max_val

# Gets the max for use in the residue mapping
def GetMaxCoords(pos, canvas):
    max_val = 0
   # max_c = 0.5

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
    x, y = canvas.rescalePt([max_val, max_val])
    print("max: ")
    print(x)
    return x
##############################################################################################################
# Get the x and y positions

def GetAtomPos(pos, AtomNum):

    CurrNum = 1
    Coords = []
    print("Atom Num")
    print(AtomNum)
    # Cycle through to find the correct atomic number
    for atom in pos:
        if CurrNum == AtomNum:
            Coords = atom
            CurrNum += 1
        else:
            CurrNum += 1
    return Coords[0], Coords[1]
##############################################################################################################

def GetResCoords(x2, y2):
    global max_val
    CurrNum = 1
    Coords = []
 #   max_val = max_val[0]
    # Variables

    r = max_val

    # y1 = mx1+b
    print("x2")
    print(x2)
    print("y2")
    print(y2)
    print("r")
    print(r)

    x1, y1 = symbols('x1 y1')
    eq1 = Eq(x1**2 + y1**2 - r**2, 0)
    eq2 = Eq((y1-y2)/(x1-x2)*x2 + 0.1, 0)

    x1, y1 = solve((eq1, eq2), x1, y1)

   # print("x1: " + str(x1) + " y1: " + str(y1))
    # r^2 = x2^2 + y2^21
   # print("get res coords")
   # print("X1: " + str(x1))
   # print("y1: " + str(y1))

    x_r = complex(x1[0]).real
    y_r = complex(y1[0]).imag

    return [x_r, y_r]



##############################################################################################################

# Finds the residues within the atoms of the ligand
# Out of those residues, a list is then created [[Res, Num, Chain], corresponding ligand atom]
def AtomsWithin(ligand, protein, cutoff):

    # load in ligand and protein
    cmd.load(ligand)
    cmd.load(protein)

    ligand_name = ligand.rsplit('.', 1)[0]
    protein_name = protein.rsplit('.', 1)[0]

    #Saves the selected atoms within x radius of the ligand
    selected_name = ligand_name + '_' + protein_name + '.pdb'
    cmd.save(selected_name, protein_name + ' within ' + str(cutoff) + ' of ' + ligand_name)
    cmd.load(selected_name)
    selected_name = selected_name.rsplit('.', 1)[0]

    stored.residues = []
    cmd.iterate(selector.process(selected_name), 'stored.residues.append(resi)')
    selection_nums = stored.residues
    print(selection_nums)

    stored.residues = []
    # Selected name minus the pdb
    cmd.iterate(selector.process(selected_name), 'stored.residues.append(chain)')
    selected_chains = stored.residues
    print(selected_chains)

    stored.residues = []
    # Selected name minus the pdb
    cmd.iterate(selector.process(selected_name), 'stored.residues.append(resn)')
    selected_resns = stored.residues
    print(selected_resns)

    x = 0
    Pair_List = []
    # Makes a pair list of resns and their corresponding numbers
    while x < len(selected_resns):
        item = [selection_nums[x], selected_resns[x], selected_chains[x]]
        if x > 0:
            if item not in Pair_List:
                Pair_List.append(item)

        else:
            Pair_List.append(item)

      #  print(selection_nums[x] + " " + selected_resns[x])

        x += 1

    for item in Pair_List:
        print("pair list")
        print(item)

    return Pair_List, protein_name, ligand_name

# Then we go through the resn list to find the ligand atom numbers
# Protein and ligand are already loaded
def LigandAtomsforResn(ResList, protein, ligand, pos, cutoff, canvas):
    # For each of the residues in the list, go through and find the atoms that correspond with the resns, then take
    # the distance a
    residue_number = 0
    for residue in ResList:
        # Residue is [Resn, number, chain, [[x1, y1], [x2, y2], [x3, y3]]]

        ligand_atoms = []
        select_string = ligand + ' within ' + cutoff + ' of ' + protein +' and resn ' + residue[1] + ' and resi ' + str(residue[0]) + ' and chain ' + residue[2]
        print(select_string)
        cmd.select('obj', select_string)
        stored.residues = []
        cmd.iterate(selector.process('obj'), 'stored.residues.append(index)')
        ligand_atoms = stored.residues
        print("resi")
        print(stored.residues)

        global posi
        # Appends the stores atoms, but we want a list of those positions
        PosList = []
        for i in stored.residues:
                x = posi[i-1][0]
                y = posi[i-1][1]
          #      print(x)
          #      print(y)
                pt0 = canvas.rescalePt([x, y])
                pt1 = GetResCoords(x, y)
                print(pt0)
                print(pt1)
                PosList.append([pt0, pt1])


        ResList[residue_number].append(PosList)
        residue_number += 1

    return ResList

'''
        x_tally = 0
        y_tally = 0
        total_num = 0
        for atom_num in ligand_atoms:
            x, y = GetAtomPos(pos, atom_num)
            x_tally += x
            y_tally += y
            total_num += 1

        # Checks what the corresponding x, y values are, and adds them to a list with the resn string and num

        ResList[residue_number].append([x_tally/total_num, y_tally/total_num])
        residue_number += 1
        '''


def MainFunc():
    mol = Chem.MolFromSmiles('C1=CC=NC(=C1)C2=NC=C(C=C2)CC(C(=O)O)N')

    print("get res coordinates")
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    # Whole position array
    pos = mol.GetConformer(0).GetPositions()
    print("pos")
    print(pos)
    global posi
    global max_val
    posi = pos
    size = (300, 300)


    canvas = Canvas(size, name='SmilesString', imageType='png')

    max_val = GetMaxCoords(pos, canvas)



   # p0 = (0, size[1] // 2)
   # p1 = (size[0] - 0, size[1] // 2)


    #canvas.addCanvasLine(p0, p1, lineWidth=2, color=(0, 0, 0))
    #
    #canvas.addCanvasText('text', [0.5 ,0.5], font='font.sans-serif', color=(0, 0, 0))

    Pair_List, protein_name, ligand_name = AtomsWithin('lig1.pdb', 'protein.pdb', '3.0')

    ResList = LigandAtomsforResn(Pair_List, protein_name, ligand_name, pos, '3.0', canvas)

    AllChem.Compute2DCoords(mol)
    fig = Draw.MolToMPL(mol, canvas=canvas)


    for item in ResList:
        print(item)
        x1_graph, y1_graph = canvas.rescalePt(item[3][0][0])
      #  x1_graph, y1_graph = item[3][0][0], item[3][1][0]

        print("x1 graph")
        print(x1_graph)

        print("y1 graph")
        print(y1_graph)

        bbox_props = dict(boxstyle="circle,pad=0.3", ec="b", lw=2)
        # Add for each of the residues that needs to happen
    #    x_r = complex(x1_graph[0]).real
    #    y_r = complex(y1_graph[0]).imag

     #   p0 = (x1_graph, size[1] // 2)
     #   p1 = (size[0] - y1_graph, size[1] // 2)

        plt.text(x1_graph, y1_graph, item[1] + ' ' + item[0], size=10, bbox=bbox_props)

    #x = canvas.rescalePt([0.4756, 1.353])
    #print(x[0])
    #print(x[1])
    #plt.text(x, x,'Test', size=10, bbox=bbox_props)

    # print(Draw.calcAtomGaussians(mol))
    #plot([x1, x2], [y1, y2], color='k', linestyle='-', linewidth=2)
    #plt.plot([0, 0], [0.5, 0.5], 'k-', lw=2)
    fig.savefig('Test123.jpeg', bbox_inches='tight', va='center', ha='center')
    # Now for each resn in Pair_List need to find the corresponding atoms

    print("Finished processing.")

MainFunc()