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
import sys
from .General import CreateHeatMap
from statistics import mode
from collections import Counter
from pymol import cmd
from .Preprocessing import PDBInfo_Wrapper
from .General import GenerateRotList, ProteintoLigList, FilterFiles, GeneralSimCheck, File_write

global olig_num, working_dir, Testing
Testing = 0
# Will be used as a switch- when it is on then we have lots of output to test, off we do not have any

######################################################################################################################

# Takes an optional PDB code and uses that for labelling purposes
def Olig_MainLoop(energy_files, PDB_code, PDB_len, cutoff):

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

    RotList, RotStruct, RotNumList = GenerateRotList(myobjects[1:])
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

    File_write(RotList, All_RMS, "RMS", "", PDB_code, working_dir)
    # "YlGnBu"  "Greens_r"
    CreateHeatMap("", All_RMS, RotList, "RMS", "rocket_r", PDB_code, working_dir)


######################################################################################################################
'''
# Takes the np array and then writes to file
def RMS_analysis(All_RMS, RotList, PDB_code):
    # store all of the rms values

    f2 = open('RMS_' + PDB_code + '.csv', 'w')
    f2.write('Root mean squared values'+'\n')
    f2.write(' ,')

    for file in RotList:
        f2.write(file + ', ')

    f2.write('\n')
    row_t = 0
    for row in All_RMS:
        f2.write(RotList[row_t] + ', ')
        for item in row:
            f2.write(str(item) + ', ')
        f2.write('\n')
        row_t = row_t + 1
'''
######################################################################################################################
'''
# Checks for similarities from a row. Saves the files to corresponding: Similar/Unique directories
def SimCheck(RotList, RotStruct, RotNumList, All_RMS, working_dir):
    print("RotStruct len: " + str(len(RotStruct)))
    print("RotList len: " + str(len(RotList)))
    cutoff = 2

    for x in RotStruct:
        print("RotStruct: " + str(x))

    # Puts the Rotation structure into this new variable
    Rot_type = [RotStruct, [0]*len(RotStruct)]

    sim_file = open('Similar.csv', 'w')
    sim_file.write('Reference, Structure, RMS, Rotation #\n')

    # The first row is the name, the second is the type: 's' for similar, nothing for unique
    # Uses the global numpy variable to check for similarities and then sorts the files.
    x = len(RotList)-1
    while x > 0:
        y = x
        while y >= 0:
            if x == y:
                pass
            else:
                if Testing:
                    print("y val: " + str(y))
                    print("x val: " + str(x))
                if All_RMS[x][y] < cutoff:
                    # Structure
                    Rot_type[1][y] = 'u'
                   # print(Rot_type[0][y] + " is similar")
                    # Need to make one unique out of the set
                    # Reference
                    Rot_type[1][x] = 's'

                    # Could now write to a file to indicate which files are similar
                    sim_file.write(str(Rot_type[0][x]) + ', ' + str(Rot_type[0][y]) + ', ' +
                                   str(All_RMS[x][y]) + ', ' + RotNumList[y] + '\n')

                # Not similar; pass for now
                else:
                    pass

            y = y - 1
        x = x - 1
    # Use this to create two new lists, sList, uList
  #  simfinder = RotStruct
    # Check to see if the directories exist
    if not os.path.exists('Similar'):
        os.makedirs('Similar')

    if not os.path.exists('Unique'):
        os.makedirs('Unique')

    # Go through the structures again and add to the proper directories
    for w in range(len(RotStruct)):
        # Name
        rot = str(Rot_type[0][w])

        if Rot_type[1][w] == 's':
            # Move from dir to dir

            # Want a new path that has the directory similar on it
            path_1 = os.path.join(working_dir, rot)
            print(path_1)
            path_2 = os.path.join(working_dir, 'Similar', rot)
            print(path_2)
            shutil.move(path_1, path_2)

        else:
            path_1 = os.path.join(working_dir, rot)
            print(path_1)
            path_2 = os.path.join(working_dir, 'Unique', rot)
            print(path_2)
            shutil.move(path_1, path_2)

    sim_file.close()
'''
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
def OligWrapper(pfile, PDB_code, cutoff):
    global olig_num, working_dir

    working_dir = os.path.dirname(pfile)
    os.chdir(working_dir)

    cmd.reinitialize()
    Saved_Complexes = ProteintoLigList(pfile)
    Saved_Complexes.sort()

    PDBfiles = Saved_Complexes

    # First, reinitialize the PyMOL window
    cmd.reinitialize()
    energy_file = PDBfiles[0]
    print("energy file: " + energy_file)
    PDB_len = len(PDBfiles)

    # Assigns the min and max res, as well as the chain to use: A, B, C
    res_min, res_max, toAlph = PDBInfo_Wrapper(energy_file)
    olig_num = len(toAlph)
    for x in toAlph:
        print(x)
    print("olig_num: " + str(olig_num))
    if Testing:
        print(olig_num)

    Olig_MainLoop(PDBfiles, PDB_code, PDB_len, cutoff)

  #  del(RotList)
    print("RMS finished processing.")
  #  cmd.reinitialize()

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
