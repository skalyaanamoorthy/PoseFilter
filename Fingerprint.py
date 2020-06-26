import pymol
import os
from pymol import cmd
#import fpkit.similarity as fps
#import fpkit.filters as filters
#import pandas as pd
import oddt
import oddt.fingerprints as fp
#import os
import numpy as np
import glob
import sys
import numpy as np
from .General import GenerateRotList
from .General import CreateHeatMap
from .General import FilterFiles, GeneralSimCheck, File_write, DirSearch, natural_sort
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
    # Setting up the array
    ligsize = len(Listoflig)
    All_Fingerprint = np.zeros((ligsize, ligsize))

    cur_row = 0
    for ref in Listoflig:

        # Divide by the type of fingerprint
        if Type == "SPLIF":
            # Fill up the row with the fingerprint data
            All_Fingerprint[cur_row, :] = SPLIF_Fingerprint(ref, Listoflig, proteinName)

        elif Type == "Interaction":
            All_Fingerprint[cur_row, :] = Interaction_Fingerprint(ref, Listoflig, proteinName)

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

    # Read in protein
    protein = next(oddt.toolkit.readfile('pdb', proteinpath))
    protein.protein = True

    # Read in and define the reference ligand
    ref_ligand = next(oddt.toolkit.readfile('pdb', ref_input))
    ref = fp.SPLIF(ref_ligand, protein)

    # Loop through each ligand in the list
    count = 0
   # print(Listoflig)
    for ligandpath in Listoflig:
        ligand = next(oddt.toolkit.readfile('pdb', ligandpath))
        fp_query = fp.SPLIF(ligand, protein)

        # similarity score for current query
        F_Scores[count] = fp.similarity_SPLIF(ref, fp_query, rmsd_cutoff=3.)
        count = count + 1
    return F_Scores


######################################################################################################################

# Convert the pdbqt files to pdb using pymol and then use the toolkit

# Take the .pdb ligand files (List of ligands) and the .pdb protein file
def Interaction_Fingerprint(ref_input, Listoflig, proteinpath):
    F_Scores = [0]*len(Listoflig)

    # Read in protein
    protein = next(oddt.toolkit.readfile('pdb', proteinpath))
    protein.protein = True

    # Read in and define the reference ligand
    ref_ligand = next(oddt.toolkit.readfile('pdb', ref_input))
    ref = fp.InteractionFingerprint(ref_ligand, protein)

    # Loop through each ligand in the list
    count = 0
    for ligandpath in Listoflig:
        ligand = next(oddt.toolkit.readfile('pdb', ligandpath))
        fp_query = fp.InteractionFingerprint(ligand, protein)

        # similarity score for current query
        cur_score = fp.dice(ref, fp_query)
        F_Scores[count] = cur_score
        count = count + 1
    return F_Scores

######################################################################################################################

# Takes in a list with some extra files that may not be ligands, only outputs the files that are of the same format
def ListtoLig(all_files):
    # Check the files for duplicates

    ListofDup = []
    for item in all_files:
        res = ''.join([i for i in item if not i.isdigit()])
        if ListofDup is None:
            # add element and num

            # add to list of duplicates
            ListofDup.append([res, 0, [item]])

        else:
            # Check each item
            d_count = 0
            added_item = 0
            for y in ListofDup:
                if y[0] == res:
                    # Add one to that entry
                    ListofDup[d_count][1] = y[1] + 1
                    ListofDup[d_count][2].append(item)
                    added_item = 1
                d_count += 1

            if added_item == 0:
                # add entry to listofdup
                ListofDup.append([res, 0, [item]])

    max_item = ""
    max_num = 0
    max_dup = 0
    # Go through each of the items from the list and
    dup_index = 0

    for dup in ListofDup:
        if dup[1] > max_num:
            max_item = dup[0]
            max_num = dup[1]
            max_dup = dup_index
        dup_index += 1

    # returns a list of the ligands, including their extension
    return(ListofDup[max_dup][2])

######################################################################################################################

# Take the .pdb ligand files (List of ligands) and the .pdb protein file
def Simple_Interaction_Fingerprint(ref_input, Listoflig, proteinpath):
    F_Scores = [0]*len(Listoflig)

    # Read in protein
    protein = next(oddt.toolkit.readfile('pdb', proteinpath))
    protein.protein = True

    # Read in and define the reference ligand
    ref_ligand = next(oddt.toolkit.readfile('pdb', ref_input))
    ref = fp.SimpleInteractionFingerprint(ref_ligand, protein)

    # Loop through each ligand in the list
    count = 0
    for ligandpath in Listoflig:
        ligand = next(oddt.toolkit.readfile('pdb', ligandpath))
        fp_query = fp.SimpleInteractionFingerprint(ligand, protein)

        # similarity score for current query
        cur_score = fp.dice(ref, fp_query)
        F_Scores[count] = cur_score
        count = count + 1
    # Returns a list of the fingerprint scores
    return F_Scores

######################################################################################################################

def PDBQTtoPDB(savename, mol):
    cmd.load(mol)
    save_as = savename + ".pdb"
    mol_s = mol.rsplit(".", 1)[0]
    cmd.save(save_as, mol_s)
    # returns the name of the newly made pdb file
    return save_as

######################################################################################################################

def ComplextoLigand(ComplexId, ResId):

    # for each file
    files = DirSearch(ComplexId)
    Ligands = []

    for complex in files:
        cmd.load(complex)

        # Ligand name, file extension
        l_name, file_ext = complex.rsplit('.', 1)

        cmd.save(l_name + '_l.pdb', l_name + ' and resn ' + ResId)
        Ligands.append(l_name + '_l.pdb')

    # Protein file
    protein = 'protein_l.pdb'
    # Ligand name, file extension
    ligand_name, file_ext = files[0].rsplit('.', 1)
    cmd.save('protein_l.pdb', ligand_name + ' and not resn ' + ResId)

    return Ligands, protein


######################################################################################################################

# Type is a list, goes through each in the list
def Fingerprint_Wrapper(info, Type, PDB_code, SI_cutoff, SPLIF_cutoff, TextInteraction):
    working_dir = ""
    # Tab 1, just the ligands and the protein file

    # Files to clean up
    Clean_up = []

    if len(info) == 2:
        pfile, LigandID = info
        working_dir = os.path.dirname(pfile)
        os.chdir(working_dir)
        protein_name = pfile

        # Ligand files from the DirSearch, either .pdb or other
        Ligand_files = DirSearch(LigandID)
        Ligand_Names_total = list(Ligand_files)
      #  for x in Ligand_files:
       #     print("FP")
        #    print("Ligand files:" + x)

        # Ligand name, file extension
        ligand_name, file_ext = Ligand_files[0].rsplit('.', 1)

        PDBLigands = []
        if file_ext is not "pdb":
            for lig in Ligand_files:
                ligand_name, file_ext = lig.rsplit('.', 1)
                PDBQTtoPDB(ligand_name + '_l', lig)
                PDBLigands.append(ligand_name + '_l.pdb')

            PDBQTtoPDB('protein_l', os.path.basename(pfile))
            protein_name = 'protein_l.pdb'

            # Add to clean up list
            Clean_up = PDBLigands

            Ligand_files = PDBLigands

        # need to check if the ligand files are .pdb

        #Ligand_files.sort()
        lf = natural_sort(Ligand_files)
        Ligand_files = lf

        LigFiles = natural_sort(Ligand_files)
        Ligand_files = LigFiles

    else:
        # Tab 2

        pdir, ComplexID, LigResID = info
        working_dir = pdir
        os.chdir(pdir)

        # Complex to ligand, by using the residue identifier
        Ligand_files, protein_name = ComplextoLigand(ComplexID, LigResID)
        Ligand_Names_total = list(Ligand_files)

        # Add to clean up list
        Clean_up = Ligand_files

        lf = natural_sort(Ligand_files)
        Ligand_files = lf


    # For each of the values in the list
    for fprint in Type:
        # Define these variables
        Listoflig, All_Fingerprint = Fingerprint(protein_name, Ligand_files, fprint)

      #  SI_cutoff,SPLIF_cutoff

        if fprint is "SPLIF":
            cutoff = SPLIF_cutoff

        else:
            cutoff = SI_cutoff

        RotNumList = ["0"] * len(Listoflig)
        PoseObjects = []

        # Rotlist (minus the .pdb), Rotstruct
        FPList = []
        for el in Listoflig:
            el_name, ext = el.rsplit('.', 1)
            FPList.append(el_name)

        RotationLabel = []
        for rotLab in Ligand_Names_total:
            RotationLabel.append(rotLab.replace('_l', ''))

        RotationLabel = natural_sort(RotationLabel)
        for index in range(len(FPList)):
            item = RMSInfo(RotationLabel[index], "", FPList[index])
            PoseObjects.append(item)

        GeneralSimCheck(PoseObjects, All_Fingerprint, working_dir, "Fingerprint", float(cutoff), fprint)
        # Write to a .csv file

        # 2D resn, and corresponding atoms
        if TextInteraction == 1:
            x = 0
           # interaction_pathway = os.path.join(working_dir, 'Fingerprint')
            os.chdir(working_dir)
            InteractionCheck(protein_name, PoseObjects)

        File_write(natural_sort(RotationLabel), All_Fingerprint, "Fingerprint", fprint, PDB_code, working_dir)
        CreateHeatMap("Fingerprint", All_Fingerprint, natural_sort(RotationLabel), fprint, "rocket", PDB_code, working_dir)
    # In the case that there were complexes, we want to clean up the files that we made
  #  os.chdir(working_dir)
    for lig in Clean_up:
        os.remove(lig)

    os.remove(protein_name)
    print("Fingerprint finished processing.")

######################################################################################################################
