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
from rdkit.Chem.rdmolfiles import MolToMolFile
global proteinpath, molpath
from oddt.toolkits.rdk import Outputfile
import shutil

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


def FindPymolType(AtomID):
    global proteinpath
    pbase = os.path.basename(proteinpath)
    pname = pbase.split('.')[0]
  #  print(molpath)
    cmd.load(proteinpath)
    stored.name = []
    cmd.iterate(selector.process(pname + " and id " + str(AtomID)), 'stored.name.append(name)')
    cmd.delete(pname)

    if stored.name == []:
        return ""
    else:
        return stored.name[0]


# Protein molecule, and a list of ligands
def InteractionCheck(ppath, Listoflig, cur_dir):
    global proteinpath
    proteinpath = ppath
    os.chdir(os.path.dirname(proteinpath))
#    pname = os.path.basename(proteinpath)

    # protein = next(oddt.toolkit.readfile('pdb', proteinpath, removeHs=False, cleanupSubstructures=False, sanitize=False))
    try:
        protein = next(oddt.toolkit.readfile('pdb', proteinpath, removeHs=False))
        protein.protein = True
    except Exception as e:

        print("Input structure could not be split into protein and ligand. Please check ligand identifier.")
        f2 = open(os.path.join(os.path.basename(proteinpath), 'ErrorLog.txt'), 'w')
        f2.write(str(e))
        f2.close()


    for ligand_object in Listoflig:
        ligandname = ligand_object.PoseNameExt

        ResReport = ligand_object.PoseName + "_ResidueReport.csv"
        path = os.path.join(cur_dir, 'Fingerprint', ResReport)

        file = open(path, 'w')
        file.write("Ligand interactions with protein residues\n")
        file.close()

        # Read in and define the reference ligand
        ligand = next(oddt.toolkit.readfile('pdb', ligandname, removeHs=False))

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
    global proteinpath
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
                p_coords = protein[y]['coords']
                l_coords = ligand[y]['coords']

                l_atomicnum = ligand[y]['atomicnum']
                p_atomicnum = protein[y]['atomicnum']

                p_atomsym = ANumtoASym(p_atomicnum, os.path.dirname(__file__))

                # To adjust for the 0 versus 1 indexing
                p_atomsym = FindPymolType(protein[y]['id']+1)

                l_atomsym = ANumtoASym(l_atomicnum, os.path.dirname(__file__))

                l_type = str(l_atomsym) + ' ' + str(ligand[y]['id']+1)
                p_type = str(p_atomsym) + ' ' + str(protein[y]['id']+1)

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
            file.close()
            return symlist[1]
    file.close()
    return 'X'
