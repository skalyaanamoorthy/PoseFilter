import pymol
import os
from pymol import cmd
#import fpkit.similarity as fps
#import fpkit.filters as filters
#import pandas as pd
#from .oddt_fingerprint import SimpleInteractionFingerprint, SPLIF, similarity_SPLIF, InteractionFingerprint, dice
#import os
import numpy as np
import glob
import sys
import numpy as np
import oddt
import oddt.fingerprints as fp
from .General import GenerateRotList
from .General import CreateHeatMap
from .General import GeneralSimCheck, File_write
from .Interactions2D import InteractionCheck
from .General import RMSInfo
from .Preprocessing import PDBInfo_Wrapper

'''

oddt and rdkit modules were used to create the SPLIF and Interaction fingerprints. These modules were also used to
examine the types of protein-ligand interactions that were then printed to .csv files.
Here are the interaction functions, some of which were used: https://oddt.readthedocs.io/en/latest/_modules/oddt/interactions.html
Here are the fingerprint functions: https://oddt.readthedocs.io/en/latest/_modules/oddt/fingerprints.html

'''

#######################################################################################################

# Takes the np array and then writes to file
# Takes in a LigandList, Array to write to the file, and type of fingerprint
def Fingerprint_write(LigList, FArr, Type, PDB_code):

    f2 = open(Type + '_' + PDB_code + '_' + '_Fingerprint.csv', 'w')
    f2.write('Fingerprint values'+'\n')
    f2.write(' ,')

    for file in LigList:
        f2.write(file + ', ')

    f2.write('\n')
    row_t = 0
    for row in FArr:
        f2.write(LigList[row_t] + ', ')
        for item in row:
            f2.write(str(item) + ', ')
        f2.write('\n')
        row_t = row_t + 1

######################################################################################################################

# Take the protein name, list of ligands and then the type of fingerprint required (SimpleInteraction, Interaction, SPLIF)

# Type is "SPLIF" or "Simple_Interaction"
def Fingerprint(proteinName, Listoflig, Type):
  #  print("List of Lig: ")
  #  print(Listoflig)

    # Setting up the array
    ligsize = len(Listoflig)
    All_Fingerprint = np.zeros((ligsize, ligsize))

    cur_row = 0
    for ref in Listoflig:
      #  print("Type: ")
      #  print(Type)
        # Divide by the type of fingerprint
        if Type == "SPLIF":
            # Fill up the row with the fingerprint data
            All_Fingerprint[cur_row, :] = SPLIF_Fingerprint(ref, Listoflig, proteinName)

        else:
            All_Fingerprint[cur_row, :] = Simple_Interaction_Fingerprint(ref, Listoflig, proteinName)
        # Puts list into the corresponding row of the np array

        cur_row += 1

    # Ideally the create heat map would be in the wrapper
    return Listoflig, All_Fingerprint

######################################################################################################################

# Convert the pdbqt files to pdb using pymol and then use the toolkit

# Take the .pdb ligand files (List of ligands) and the .pdb protein file
# for each, make a fingerprint and output
def SPLIF_Fingerprint(ref_input, Listoflig, proteinpath):
    F_Scores = [0]*len(Listoflig)

    #protein = next(oddt.toolkit.readfile('pdb', proteinpath, removeHs=False, cleanupSubstructures=False, sanitize=False))
    protein = next(oddt.toolkit.readfile('pdb', proteinpath, removeHs=False))
    protein.protein = True

    # Read in and define the reference ligand
    #ref_ligand = next(oddt.toolkit.readfile('pdb', ref_input, removeHs=False, cleanupSubstructures=False, sanitize=False))
    ref_ligand = next(oddt.toolkit.readfile('pdb', ref_input, removeHs=False))
    ref = fp.SPLIF(ref_ligand, protein)

    # Loop through each ligand in the list
    count = 0
   # print(Listoflig)
    for ligandpath in Listoflig:
        #ligand = next(oddt.toolkit.readfile('pdb', ligandpath, removeHs=False, cleanupSubstructures=False, sanitize=False))
        ligand = next(oddt.toolkit.readfile('pdb', ligandpath, removeHs=False))
        fp_query = fp.SPLIF(ligand, protein)
        # similarity score for current query
        F_Scores[count] = fp.similarity_SPLIF(ref, fp_query, rmsd_cutoff=3.)
        count = count + 1
    return F_Scores


######################################################################################################################

# Take the .pdb ligand files (List of ligands) and the .pdb protein file
def Simple_Interaction_Fingerprint(ref_input, Listoflig, proteinpath):
    F_Scores = [0]*len(Listoflig)

    # Read in protein
    #protein = next(oddt.toolkit.readfile('pdb', proteinpath, removeHs=False, cleanupSubstructures=False, sanitize=False))
    protein = next(oddt.toolkit.readfile('pdb', proteinpath, removeHs=False))
    protein.protein = True

    # Read in and define the reference ligand
    #ref_ligand = next(oddt.toolkit.readfile('pdb', ref_input, removeHs=False, cleanupSubstructures=False, sanitize=False))
    ref_ligand = next(oddt.toolkit.readfile('pdb', ref_input, removeHs=False))
    ref = fp.SimpleInteractionFingerprint(ref_ligand, protein)

    # Loop through each ligand in the list
    count = 0
    for ligandpath in Listoflig:
        #ligand = next(oddt.toolkit.readfile('pdb', ligandpath, removeHs=False, cleanupSubstructures=False, sanitize=False))
        ligand = next(oddt.toolkit.readfile('pdb', ligandpath, removeHs=False))
        fp_query = fp.SimpleInteractionFingerprint(ligand, protein)

        # similarity score for current query
        cur_score = fp.dice(ref, fp_query)
        F_Scores[count] = cur_score
        count = count + 1
    # Returns a list of the fingerprint scores
    return F_Scores

######################################################################################################################

# Type is a list, goes through each in the list
def Fingerprint_Wrapper(files, Type, PDB_code, SI_cutoff, SPLIF_cutoff, TextInteraction, cur_dir, pname):
    ligand_path = os.path.join(cur_dir, "PDBLigand")
    protein_path = os.path.join(cur_dir, "PDBLigand", pname)

    # For each of the values in the list
    for fprint in Type:
        os.chdir(ligand_path)

        # Define these variables
        Listoflig, All_Fingerprint = Fingerprint(pname, files, fprint)


        if fprint == "SPLIF":
            cutoff = SPLIF_cutoff

        else:
            cutoff = SI_cutoff

        PoseObjects = []

        for index in range(len(files)):
            # PoseName, PoseNameExt, RotNum
            item = RMSInfo(files[index].split('.')[0], files[index], "")
            PoseObjects.append(item)

        os.chdir(cur_dir)
        GeneralSimCheck(PoseObjects, All_Fingerprint, cur_dir, "Fingerprint", float(cutoff), fprint)
        # Write to a .csv file

        # 2D resn, and corresponding atoms
        if TextInteraction == 1:
            x = 0
           # interaction_pathway = os.path.join(working_dir, 'Fingerprint')
            os.chdir(cur_dir)
            InteractionCheck(protein_path, PoseObjects, cur_dir)
        os.chdir(cur_dir)
        File_write(files, All_Fingerprint, "Fingerprint", fprint, PDB_code, cur_dir)
        CreateHeatMap("Fingerprint", All_Fingerprint, files, fprint, "mako", PDB_code, cur_dir)

    print("Fingerprint finished processing.")

######################################################################################################################

