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

global olig_num, working_dir, Testing
Testing = 0
# Will be used as a switch- when it is on then we have lots of output to test, off we do not have any

######################################################################################################################

# Takes an optional PDB code and uses that for labelling purposes
def Olig_MainLoop(energy_files, PDB_code, PDB_len, cutoff, UNK_var, alpha):

    global olig_num, working_dir
    print("energy files: ")
    print(energy_files)

    # First, load the first file as a reference
    cmd.load(energy_files[0], 'obj1')
    print(energy_files[0])

    # Resets the olig num here
  #  ChainCheck('obj1')
  #  print(olig_num)
   # files =
 #   energy_files = files
    # Gets all objects of this type

    fdsf = FilterFiles(energy_files)
    myobjects = cmd.get_object_list()
    for obj in myobjects:
        print("obj: " + obj)

    RotList, RotStruct, RotNumList = GenerateRotList(myobjects[1:], UNK_var, alpha)
    # Delete obj1
    cmd.delete('obj1')

    All_RMS = np.zeros((PDB_len, PDB_len))
    row_n = 0
    # For the designated lowest rotation, we cycle through
    for obj_x in RotList:
        print("RotList: " + obj_x)
        # Calculates RMSList values and puts into list
        RMSList = RMS_Calc(obj_x, RotList)

        # Puts list into the corresponding row of the np array
        All_RMS[row_n, :] = RMSList

        # Clean up the extra UNK files now
        #os.remove(obj_x + '.pdb')

        row_n = row_n + 1

    # Print RMS values to .CSV
    # Replace with general analysis writing thing
   # RMS_analysis(All_RMS, RotList, PDB_code)
#(LigList, Total_Array, Type, TypeInfo, PDB_code, working_dir)
    # Analyze the symmetry


    # Check if the number is a float/int or something
    GeneralSimCheck(RotList, RotStruct, RotNumList, All_RMS, working_dir, "RMS", float(cutoff), "")

    RotationLabel = []
    for rotLab in RotStruct:
        RotationLabel.append(rotLab.replace('complex_', ''))

    File_write(RotationLabel, All_RMS, "RMS", "", PDB_code, working_dir)
    # "YlGnBu"  "Greens_r"
    CreateHeatMap("", All_RMS, RotationLabel, "RMS", "rocket_r", PDB_code, working_dir)


######################################################################################################################
######################################################################################################################
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

######################################################################################################################

#def InputFileDir(pfile, PDBCode):
#    cmd.reinitialize()
#    Saved_Complexes = ProteintoLigList(pfile)
#    Saved_Complexes.sort()
#    OligWrapper(Saved_Complexes, PDBCode)

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
  #  UNK_var = ""

    InfoProtein = ""


    # Tab 1, protein file and ligand list
    if len(info) == 2:
        pfile, LigandID = info
        working_dir = os.path.dirname(pfile)
        os.chdir(working_dir)

        # Define just the protein file
        InfoProtein = os.path.basename(pfile)
        #     file_ext = protein_name.rsplit('.', 1)[1]

        # Creates ligand and protein complexes
        Saved_Complexes, LigandFile = LigandtoComplex(pfile, LigandID)

        # Just want the names of the ligand files in order:
        Ligand_File_Names = DirSearch(LigandID)

        # sort
      #  Saved_Complexes.sort()
        sc = natural_sort(Saved_Complexes)
        Saved_Complexes = sc

        cmd.load(LigandFile)

        # Just get the ligand name and not the extension
        LigandName = LigandFile.rsplit('.', 1)[0]

        # Check the residue ID
        stored.residues = []
        #cmd.iterate(selector.process('pose10'),'stored.residues.append(resn)')
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
        Saved_Complexes = DirSearch(ComplexID)
        Saved_Complexes.sort()

        # Need to prepare just the protein for processing
        cmd.load(Saved_Complexes[0])
        InfoProtein = 'OnlyProtLig.pdb'
        LName = Saved_Complexes[0].rsplit('.', 1)[0]
        cmd.save('OnlyProtLig.pdb', LName + ' and not resn ' + LigResID)

    # Preprocessing, need to make sure everything is a complex

    # First, reinitialize the PyMOL window
    cmd.reinitialize()
    energy_file = Saved_Complexes[0]
    print("energy file: " + energy_file)
    PDB_len = len(Saved_Complexes)

    # Assigns the min and max res, as well as the chain to use: A, B, C
    # Should be using the protein for the PDBInfo
    res_min, res_max, toAlph = PDBInfo_Wrapper(InfoProtein)
    olig_num = len(toAlph)
    for x in toAlph:
        print(x)
    print("olig_num: " + str(olig_num))
    if Testing:
        print(olig_num)

    Olig_MainLoop(Saved_Complexes, PDB_code, PDB_len, cutoff, UNK_var, alpha)

    if len(info) == 2:
        # Clean up files afterwards
        print("complex clean up files")
        for complex in Saved_Complexes:
            print("remove: " + complex)
            os.remove(complex)
       # print("remove: " + protein_name)
       # os.remove(protein_name)

  #  del(RotList)
    print("RMS finished processing.")
  #  cmd.reinitialize()
cmd.extend('OligWrapper', OligWrapper)
# Test cases

#InputFileDir(r"C:\Users\Justine\PycharmProjects\Oligomer_script\Vina_docking\Dimers\Docking\1FX9\1fx9.pdbqt")
#InputFileDir(r"C:\Users\Justine\PycharmProjects\Oligomer_script\Vina_docking\Tetramers\3J5P\protein.pdbqt")
# Tetramer
#OligWrapper(2, r"C:\Users\Justine\PycharmProjects\Oligomer_script\Vina_docking\Tetramers\Oligomer analysis", "test")
#OligWrapper(4, "", "5va1")

# Dimer
#OligWrapper(2, r"C:\Users\Justine\PycharmProjects\Oligomer_script\Vina_docking\Dimers\Oligomer_Analysis\1FX9", "")
#OligWrapper(2, r"C:\Users\Justine\PycharmProjects\Oligomer_script\Vina_docking\Dimers\1FX9\Testing", "")

#Trimer
#OligWrapper(3, r"C:\Users\Justine\PycharmProjects\Oligomer_script\Vina_docking\Trimers\Oligomer_analysis\5EIL", "")
#OligWrapper("tri", "./1rer_Cisa_01_dummy_0.pdb", "1rer")

#OligWrapper(sys.argv[0], sys.argv[1], sys.argv[2])
