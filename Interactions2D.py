import pymol
import os
from pymol import cmd
from pymol import stored
from pymol import selector

#import os
import numpy as np
import glob
import sys
import oddt
import oddt.interactions as interactions
import numpy as np
from rdkit import Chem
from scipy.spatial import distance
global Calpha, Cbeta

""" Note that the elementlist.csv file is from https://sciencenotes.org/list-elements-atomic-number/ 
"""

class ProtLigInfo:
  def __init__(self, LigAtomType, ProtAtomType, ProtResName,
               ProtResNum, Dist):
    self.LigAtomType = LigAtomType
    self.ProtAtomType = ProtAtomType
    self.ProtResName = ProtResName
    self.ProtResNum = ProtResNum
    self.Dist = Dist

# Take the two molecules, protein and ligand and use each of those oddt modules

def FindCACB(proteinpath):
    cmd.load(proteinpath)
    # cmd.select("crystal and n. CB")
    pname = os.path.basename(proteinpath)
    namenoext = pname.split('.')[0]
   # cmd.select(pname + " and n. CA", 'sele')
    stored.idA =[]
    cmd.iterate(selector.process(namenoext + " and n. CA"), 'stored.idA.append(ID)')

    stored.idB = []
    cmd.iterate(selector.process(namenoext + " and n. CB"), 'stored.idB.append(ID)')
    return stored.idA, stored.idB


# Protein molecule, and a list of ligands
def InteractionCheck(proteinpath, Listoflig, cur_dir):
    global Calpha, Cbeta
    Calpha, Cbeta = FindCACB(proteinpath)
    os.chdir(os.path.dirname(proteinpath))

    protein = next(oddt.toolkit.readfile('pdb', proteinpath))
    protein.protein = True

    for ligand_object in Listoflig:
        ligandname = ligand_object.PoseNameExt

        #Ligand_Name = ligandname.split('.')[0]

        ResReport = ligand_object.PoseName + "_ResidueReport.csv"
        path = os.path.join(cur_dir, 'Fingerprint', ResReport)

        file = open(path, 'w')
        file.write("Ligand interactions with protein residues\n")
        file.close()

        # Read in and define the reference ligand
        ligand = next(oddt.toolkit.readfile('pdb', ligandname))

        # Hydrophobic interactions
        p_hydroph, l_hydroph = interactions.hydrophobic_contacts(protein, ligand)
        InteractionsFile(p_hydroph, l_hydroph, path, 'hydrophobic')

        # h bonds
        p_hbonds, l_hbonds, strict = interactions.hbonds(protein, ligand)
        InteractionsFile(p_hbonds, l_hbonds, path, 'hydrogen bond')

        # halogens
        p_halogen, l_halogen, strict = interactions.halogenbonds(protein, ligand)
        InteractionsFile(p_halogen, l_halogen, path, 'halogen bond')

        # pistacking bonds
        pi_interactions = interactions.pi_stacking(protein, ligand)
        InteractionsFile(pi_interactions[0], pi_interactions[2], path, 'pi stacking')

        # salt bridges
        p_salt_bridges, l_salt_bridges = interactions.salt_bridges(protein, ligand)
        InteractionsFile(p_salt_bridges, l_salt_bridges, path, 'salt bridge')

        # pi_cation
        p_pi_cation, l_pi_cation, strict = interactions.pi_cation(protein, ligand)
        InteractionsFile(p_pi_cation, l_pi_cation, path, 'pi cation')

        # acceptor_metal bonds
        p_acceptor_metal_a, acceptor_metal_a, strict = interactions.acceptor_metal(protein, ligand)
        InteractionsFile(p_acceptor_metal_a, acceptor_metal_a, path, 'acceptor metal')

        # pi_metal bonds

        p_pi_metal, l_pi_metal, strict = interactions.pi_metal(protein, ligand)
        InteractionsFile(p_pi_metal, l_pi_metal, path, 'pi metal')


def InteractionsFile(protein, ligand, FilePath, Interaction_Name):
    # If it is not an empty array
    global Calpha, Cbeta
    if len(protein) != 0:

        OneLine = 0
        for x in range(len(protein)):

            try:
                p_coords = protein[x]['coords']
                OneLine = 1

            except ValueError:
                # Key is not present
                pass

        if OneLine:
            file = open(FilePath, 'a')
            file.write("Ligand " + Interaction_Name + " interactions with protein residues\n")
            file.write('Ligand Atom Numbering, Residue Atom Numbering,' +
                       'Resn, Resi, Ligand to Protein Distance (Angstrom)\n')
        ListOfDist = []


        # Put everything into a list; use an object so that it is neater
        for y in range(len(protein)):

            try:
              #  print('coords are here')
                p_coords = protein[y]['coords']
                l_coords = ligand[y]['coords']

                # want the atomic number, and then convert to symbol later.
                # For AA, want isbeta, isalpha

                l_atomicnum = ligand[y]['atomicnum']
                p_atomicnum = protein[y]['atomicnum']

             #   if protein[y]['isalpha']:
             #       p_atomsym = ANumtoASym(p_atomicnum, os.path.dirname(__file__)) + 'A'
             #   elif protein[y]['isbeta']:
            #        p_atomsym = ANumtoASym(p_atomicnum, os.path.dirname(__file__)) + 'B'
             #   else:
             #       p_atomsym = ANumtoASym(p_atomicnum, os.path.dirname(__file__))

                p_atomsym = ANumtoASym(p_atomicnum, os.path.dirname(__file__))

                p_id = protein[y]['id']
                if p_id in Calpha:
                    p_atomsym = p_atomsym + 'A'
                elif p_id in Cbeta:
                    p_atomsym = p_atomsym + 'B'

                l_atomsym = ANumtoASym(l_atomicnum, os.path.dirname(__file__))

                l_type = str(l_atomsym) + str(ligand[y]['id'])
                p_type = str(p_atomsym) + str(protein[y]['id'])

                dst = distance.euclidean(p_coords, l_coords)
                item = ProtLigInfo(l_type, p_type, protein[y]['resname'], str(protein[y]['resnum']), dst)
                ListOfDist.append(item)

            except ValueError:
                # Key is not present
                pass


        # Sort the list according to the distance
        ListOfDist.sort(key=lambda z: z.Dist, reverse=False)

        if len(ListOfDist) != 0:
            # Make a list of the items that are already included
            # Resn, Resi list
            ResnResi = []
            r0 = ListOfDist[0]

            file.write(r0.LigAtomType + ', ' + r0.ProtAtomType + ', ' +
                       r0.ProtResName + ', ' + r0.ProtResNum + ', ' +
                       str(round(r0.Dist, 2)) + '\n')
            ResnResi.append([r0.ProtResName, r0.ProtResNum])
           # print("first printed")
           # print(r0.ProtResName + " " + r0.ProtResNum)

            for residue in ListOfDist:
                InResnResiList = 0
                for x in ResnResi:
                    if ((residue.ProtResName == x[0]) and (residue.ProtResNum == x[1])):
                        InResnResiList = 1

                if InResnResiList == 0:
                    file.write(residue.LigAtomType + ', ' + residue.ProtAtomType + ', ' +
                               residue.ProtResName + ', ' + residue.ProtResNum + ', ' +
                               str(round(residue.Dist, 2)) + '\n')
                    ResnResi.append([residue.ProtResName, residue.ProtResNum])

            file.write('\n')
            file.close()

def ANumtoASym(ANum, dirname):
    filename = os.path.join(dirname, 'elementlist.csv')
    file = open(filename, 'r')
    lines = file.readlines()
    for line in lines:
        symlist = line.split(',')
        if symlist[0] == str(ANum):
          #  print(symlist[0])
          #  print(str(ANum))
            file.close()
            return symlist[1]
    file.close()
    return 'X'

#proteinpath = '/home/justine/protein.pdb'
#ligandpaths = ['ligand1.pdb', 'ligand2.pdb', 'ligand3.pdb']
#InteractionCheck(proteinpath, ligandpaths)

