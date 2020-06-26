import pymol
from pymol import cmd
from pymol import stored
from pymol import selector
import os
import re
import fileinput
import shutil
import glob
import numpy as np
import os
import sys
from .General import CreateHeatMap
from statistics import mode
from collections import Counter
from pymol import cmd
from .Preprocessing import PDBInfo_Wrapper
from .General import GenerateRotList, LigandtoComplex, FilterFiles, GeneralSimCheck, File_write, DirSearch, natural_sort
from .Fingerprint import PDBQTtoPDB
from .General import RMSInfo
global olig_num, working_dir, Testing
Testing = 0
# Will be used as a switch- when it is on then we have lots of output to test, off we do not have any


######################################################################################################################

# Takes an optional PDB code and uses that for labelling purposes
def Olig_MainLoop(energy_files,OriginalLigands, PDB_code, PDB_len, cutoff, UNK_var, alpha):

    global olig_num, working_dir
  #  print("energy files: ")
  #  print(energy_files)

    # First, load the first file as a reference
    cmd.load(energy_files[0], 'obj1')
  #  print(energy_files[0])

    fdsf = FilterFiles(energy_files)
    myobjects = cmd.get_object_list()
 #   for obj in myobjects:
  #      print("obj: " + obj)

    RotList, RotStruct, RotNumList = GenerateRotList(myobjects[1:], UNK_var, alpha)


    RMSItemList = []
    for y in range(len(RotList)):
        RMSItem = RMSInfo(OriginalLigands[y], RotNumList[y], RotList[y])
        RMSItemList.append(RMSItem)

    # Delete obj1
    cmd.delete('obj1')

    All_RMS = np.zeros((PDB_len, PDB_len))
    row_n = 0
    # For the designated lowest rotation, we cycle through
    for obj_x in RotList:
        # Calculates RMSList values and puts into list
        RMSList = RMS_Calc(obj_x, RotList)

        # Puts list into the corresponding row of the np array
        All_RMS[row_n, :] = RMSList

        row_n = row_n + 1

    # Check if the number is a float/int or something
    # RotStruct not needed

    GeneralSimCheck(RMSItemList, All_RMS, working_dir, "RMS", float(cutoff), "")

    RotationLabel = []
    for rotLab in RotStruct:
        RotationLabel.append(rotLab.replace('complex_', ''))

    File_write(RotationLabel, All_RMS, "RMS", "", PDB_code, working_dir)
    # "YlGnBu"  "Greens_r"
    CreateHeatMap("", All_RMS, RotationLabel, "RMS", "rocket_r", PDB_code, working_dir)

    # Remove the final rotation files at the end (others were previously removed)
    for item in RMSItemList:
        os.remove(item.ObjName + '.pdb')
        Rot_beginning = item.ObjName.rsplit('_', 1)[0]
        os.remove(Rot_beginning + '_' + item.RotNum + 'rot.pdb')


######################################################################################################################

# For a given object, figure out RMS values for it compared to the other objects, and then return it
def RMS_Calc(Ref_obj, OList):
    RMSList = [0]*len(OList)
    count = 0

    for obj in OList:
        rms_val = cmd.rms_cur(Ref_obj, obj)
        RMSList[count] = str(rms_val)
        count += 1

    # RMS values of ref_obj compared to all other objects
    return RMSList

###############################################################################################################

# Changing to incorporate selecting a protein.pdb/pdbqt file first
# Instead let's have an info parameter that's a list:
# Tab 1: [pfile, ligID]
# Tab 2: [dir, complexid, ligresid]
@cmd.extend
def OligWrapper(info, PDB_code, cutoff, alpha):
    # pfile is instead info list
    global olig_num, working_dir
    cmd.reinitialize()
    OriginalLigands = []
  #  UNK_var = ""

    InfoProtein = ""
    IsComplex = 0

    # Tab 1, protein file and ligand list
    if len(info) == 2:
        pfile, LigandID = info
        working_dir = os.path.dirname(pfile)
        os.chdir(working_dir)

        # Define just the protein file
        InfoProtein = os.path.basename(pfile)
        #     file_ext = protein_name.rsplit('.', 1)[1]

        # Creates ligand and protein complexes
        Saved_Complexes, LigandFile, OriginalLigands = LigandtoComplex(pfile, LigandID)

        OL = natural_sort(list(OriginalLigands))
        OriginalLigands = OL

        # Just want the names of the ligand files in order:
        sc = natural_sort(Saved_Complexes)
        ToDelete = list(sc)
        Saved_Complexes = sc

        cmd.load(LigandFile)

        # Just get the ligand name and not the extension
        LigandName = LigandFile.rsplit('.', 1)[0]

        # Check the residue ID
        stored.residues = []
        cmd.iterate(selector.process(LigandName),'stored.residues.append(resn)')
        UNK_var = stored.residues[0]
        cmd.delete(LigandName)

    else:
        # Tab 2
        pdir, ComplexID, LigResID = info
        os.chdir(pdir)
        working_dir = pdir

        # sets the residue ID
        UNK_var = LigResID
        Saved_Complexes = natural_sort(DirSearch(ComplexID))

        ToDelete = []
        # Check if .pdb

        OriginalLigands = natural_sort(list(Saved_Complexes))
        for lig in OriginalLigands:
            print('lig')
            print(lig)

     #   if ('pdb' not in Complex_ext):

        for compl in Saved_Complexes:
            Complex_N, Complex_ext = compl.split('.')
            PDBQTtoPDB(Complex_N+'xyzl', compl)
            ToDelete.append(Complex_N + 'xyzl.pdb')

            # set the saved complexes as those ones
            Saved_Complexes = list(ToDelete)

        for lig2 in Saved_Complexes:
            print("saved complex")
            print(lig2)

        # Need to prepare just the protein for processing
        cmd.load(Saved_Complexes[0])
        InfoProtein = 'OnlyProtLig.pdb'
        LName = Saved_Complexes[0].rsplit('.', 1)[0]
        cmd.save('OnlyProtLig.pdb', LName + ' and not resn ' + LigResID)

    # Preprocessing, need to make sure everything is a complex

    # First, reinitialize the PyMOL window
    cmd.reinitialize()
    energy_file = Saved_Complexes[0]
   # print("energy file: " + energy_file)
    PDB_len = len(Saved_Complexes)

    # Assigns the min and max res, as well as the chain to use: A, B, C
    # Should be using the protein for the PDBInfo
    res_min, res_max, toAlph = PDBInfo_Wrapper(InfoProtein)
    olig_num = len(toAlph)

    Olig_MainLoop(Saved_Complexes, OriginalLigands, PDB_code, PDB_len, cutoff, UNK_var, alpha)

#    else:
    if len(info) == 3:
        os.remove('OnlyProtLig.pdb')

    for complex in ToDelete:
    #    print("remove: " + complex)
        os.remove(complex)

        # Remove protein name


    print("RMS finished processing.")
cmd.extend('OligWrapper', OligWrapper)
