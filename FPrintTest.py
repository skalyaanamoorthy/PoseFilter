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

global Testing
Testing = 0
# Will be used as a switch- when it is on then we have lots of output to test, off we do not have any

#######################################################################################################

# Takes the np array and then writes to file
# Takes in a LigandList, Array to write to the file, and type of fingerprint
def Fingerprint_write(LigList, FArr, Type, PDB_code):

# 'Fingerprint_SPLIF.csv'

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

# Type is "SPLIF", "Interaction", or "Simple_Interaction"
def Fingerprint(proteinName, Listoflig, Type):
    # Setting up the array
    ligsize = len(Listoflig)
    All_Fingerprint = np.zeros((ligsize, ligsize))
    print("prot name: ")
    print(proteinName)

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
    print(All_Fingerprint)

    # Ideally the create heat map would be in the wrapper
    return Listoflig, All_Fingerprint

# Give it a string to indicate the file to write to as well as the type of fingerprint
#    Fingerprint_write(Listoflig, All_Fingerprint, FingerprintTypeString)

######################################################################################################################

# Convert the pdbqt files to pdb using pymol and then use the toolkit

# Take the .pdb ligand files (List of ligands) and the .pdb protein file
# for each, make a fingerprint and output
def SPLIF_Fingerprint(ref_input, Listoflig, proteinpath):
    F_Scores = [0]*len(Listoflig)
    # Read in protein
    print("protein path: " + proteinpath)
    print("list of lig: ")
    print(Listoflig)
 #   proteinpath = r"C:\Users\Justine\PycharmProjects\Oligomer_script\Vina_docking\Dimers\Docking\1FX9\1fx9.pdbqt"
    protein = next(oddt.toolkit.readfile('pdb', proteinpath))
    protein.protein = True

    # Read in and define the reference ligand
    ref_ligand = next(oddt.toolkit.readfile('pdb', ref_input))
    ref = fp.SPLIF(ref_ligand, protein)
   # print(ref)

    # Loop through each ligand in the list
    count = 0
   # print(Listoflig)
    for ligandpath in Listoflig:
        ligand = next(oddt.toolkit.readfile('pdb', ligandpath))
        fp_query = fp.SPLIF(ligand, protein)
      #  print(fp_query)

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
    print(Listoflig)
    for ligandpath in Listoflig:
        ligand = next(oddt.toolkit.readfile('pdb', ligandpath))
        fp_query = fp.InteractionFingerprint(ligand, protein)

        # similarity score for current query
        cur_score = fp.dice(ref, fp_query)
       # tan = fp.tanimoto(fp_query, ref)
        F_Scores[count] = cur_score
        print(cur_score)
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
    print(Listoflig)
    for ligandpath in Listoflig:
        ligand = next(oddt.toolkit.readfile('pdb', ligandpath))
        fp_query = fp.SimpleInteractionFingerprint(ligand, protein)

        # similarity score for current query
        cur_score = fp.dice(ref, fp_query)
       # tan = fp.tanimoto(fp_query, ref)
        F_Scores[count] = cur_score
        print(cur_score)
        count = count + 1
  #  print(F_Scores)
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
# Fingerprint(proteinName, Listoflig, Type)
# For the actual fingerprint input we need to take in the .pdb/.pdbqt protein file and then output the fingerprint
# based on that
# Take the function from .olig to
'''
# Takes in a pdb or a pdbqt file
def ProteintoLigListComplex(pfile):
    # Changes the directory
    working_dir = os.path.dirname(pfile)
    os.chdir(working_dir)

    # pfile just want the extention as well as the name + extension

    protein_name = os.path.basename(pfile)

    # Protein name, file extension
    protein_name, file_ext = protein_name.rsplit('.', 1)
    print("Protein true name: " + protein_name)

    # Saving as a .pdb file
    print("protein save 1: " + protein_name)
    print("protein save 2: " + protein_name + file_ext)
    PDBQTtoPDB(protein_name, protein_name + "." + file_ext)

    protein_name = protein_name + '.pdb'
    print("Protein and pdb again: " + protein_name)
    all_files = []

    query = glob.glob('*.' + file_ext)

    caught = 0
    # len = 0
    # Deletes the extra files
    for x in query:

        # This is the protein
        if protein_name == x:
            pass
        elif "rotation" in x:
            pass
        else:
            #     len += 1
            if caught == 0:
                energy_file = x
            caught = 1
            all_files.append(x)

    LigandList = ListtoLig(all_files)
    print(LigandList)
    Saved_Complexes = [""] * len(LigandList)
    # Goes through each ligand to make a complex
    SNum = 0
    testlig = LigandList[0]



    for ligand in LigandList:
        print("ligand: " + ligand)
        lig_name = ligand.rsplit('.', 1)[0]
        # #MakeComplex(dir, "Justlig", "lig2.pdbqt", "protein.pdbqt")

        # Commenting out for now
        print("Lig_name: " + lig_name)
        print("Ligand: " + ligand)
        PDBQTtoPDB(lig_name, ligand)

        #   print("loading complex")
        #   cmd.load(lig_name + '.pdb')

        # These are the names of the complexes that will be iterated through the program
        Saved_Complexes[SNum] = lig_name + ".pdb"
        SNum += 1
    print("Saved complexes: ")
    print(Saved_Complexes)
    print(protein_name)
    return protein_name, Saved_Complexes
'''
######################################################################################################################
# Type is a list, goes through each in the list
def Fingerprint_Wrapper(pfile, Type, PDB_code, SI_cutoff, I_cutoff, SPLIF_cutoff):

    # Reinitialize
    cmd.reinitialize()

    # Make sure the dir is correct
    working_dir = os.path.dirname(pfile)
    os.chdir(working_dir)

    # Takes in the protein file as well as the type of operation that will be performed
   # proteinName, Saved_Complexes = ProteintoLigListComplex(pfile)

    Saved_Complexes = ProteintoLigList(pfile)
    Saved_Complexes.sort()

    # First, load the first file as a reference
    cmd.load(Saved_Complexes[0], 'obj1')

    # Assigns the min and max res, as well as the chain to use: A, B, C
    res_min, res_max, toAlph = PDBInfo_Wrapper(Saved_Complexes[0])
    print(res_min)
    print(res_max)
    olig_num = len(toAlph)
    for x in toAlph:
        print(x)
    print("olig_num: " + str(olig_num))
    if Testing:
        print(olig_num)

    RotList, RotStruct, RotNumList = GenerateRotList(FilterFiles(Saved_Complexes))

    FP_List = []
    for x in range(len(RotNumList)):
       # print('FP List (RotList): ' + RotList[x] + '.pdb')
      #  print('UNK: ' + Saved_Complexes[x] + '_UNK' + str(x) + '.pdb')
        #print('Rot Struct: ' + RotStruct[x])
        FP_List.append(RotList[x] + '.pdb')

# need the protein .pdb file

    # pfile just want the extention as well as the name + extension

    protein_name = os.path.basename(pfile)

    # Protein name, file extension
    protein_name, file_ext = protein_name.rsplit('.', 1)

    # Saving as a .pdb file
    PDBQTtoPDB(protein_name, protein_name + "." + file_ext)
    # Protein name for fingeprint portion
    protein_name = protein_name + '.pdb'

    # For each of the values in the list
    for fprint in Type:
        # Define these variables
        Listoflig, All_Fingerprint = Fingerprint(protein_name, FP_List, fprint)

      #  SI_cutoff, I_cutoff, SPLIF_cutoff

        if fprint is "SPLIF":
            cutoff = SPLIF_cutoff
        elif fprint is "Interaction":
            cutoff = I_cutoff
        else:
            cutoff = SI_cutoff


        GeneralSimCheck(RotList, RotStruct, RotNumList, All_Fingerprint, working_dir, "Fingerprint", float(cutoff), fprint)
        # Write to a .csv file


        File_write(Listoflig, All_Fingerprint, "Fingerprint", fprint, PDB_code, working_dir)
        CreateHeatMap("Fingerprint", All_Fingerprint, Listoflig, fprint, "rocket", PDB_code, working_dir)

        # Make sure the dir is correct
        working_dir = os.path.dirname(pfile)
        os.chdir(working_dir)

    print("Fingerprint finished processing.")



######################################################################################################################

#Fingerprint_Wrapper('protein.pdb', 'SPLIF')
#ProteintoLigListComplex('/home/justine/PycharmProjects/PoseFilter/Test/Trimer/protein.pdbqt')
#Fingerprint_Wrapper('/home/justine/PycharmProjects/PoseFilter/Test/Trimer/protein.pdbqt', "SPLIF")

# Make sure the dir is correct
working_dir = os.path.dirname('/home/justine/Documents/DockingPart2/newprotfile.pdb')
os.chdir(working_dir)
#LigList = ['ligand1_complex_UNK0.pdb', 'ligand1_complex_UNK1.pdb', 'ligand1_complex_UNK2.pdb']
#x =SPLIF_Fingerprint('ligand1_complex_UNK0.pdb', LigList, protein)
#print(x)

#x =SPLIF_Fingerprint('ligand1_complex_UNK0.pdb', LigList, protein)
#print(x)

#x =Simple_Interaction_Fingerprint('ligand1_complex_UNK0.pdb', LigList, protein)
#print(x)

LigList = ['pose1.pdb', 'pose2.pdb', 'pose3.pdb']
protein = '/home/justine/Documents/DockingPart2/newprotfile.pdb'

x =SPLIF_Fingerprint('pose1.pdb', LigList, protein)
print(x)

x =Interaction_Fingerprint('pose1.pdb', LigList, protein)
print(x)

x =Simple_Interaction_Fingerprint('pose1.pdb', LigList, protein)
print(x)