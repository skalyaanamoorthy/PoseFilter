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
from collections import Counter
from pymol import cmd
from .Preprocessing import PDBInfo_Wrapper
from .General import GenerateRotList, GeneralSimCheck, File_write, DirSearch, natural_sort
from .General import RMSInfo
global olig_num, working_dir, Testing
Testing = 0
# Will be used as a switch- when it is on then we have lots of output to test, off we do not have any

######################################################################################################################
# Takes an optional PDB code and uses that for labelling purposes
# Should be able to get the files from the directory main_dir/PDBComplex

def CreateRMS(files, PDB_code, cutoff, UNK_var, alpha, nonidentical, cur_dir, pname):
    print(files)
    global olig_num, working_dir
    working_dir = cur_dir
    protein_path = os.path.join(cur_dir, "PDBLigand", pname)
   # cmd.do('set retain_order,1')
    res_min, res_max, toAlph = PDBInfo_Wrapper(protein_path, nonidentical)
    olig_num = len(toAlph)
  #  print("res min : " + str(res_min))
  #  print("res max : " + str(res_max))
  #  print("to alpha ")
  #  print(toAlph)


    # First, load the first file as a reference
    cmd.reinitialize()
    cmd.do('set retain_order,1')
  #  cmd.do('set pdb_retain_ids,1')

    os.chdir(os.path.join(cur_dir, "PDBComplex"))

    cmd.load(files[0], 'obj1')
    files_no_ext = []
    for file in files:
        file_name, ext = file.split('.')
        cmd.load(file)
        files_no_ext.append(file_name)

    RotList, RotStruct, RotNumList = GenerateRotList(files_no_ext, UNK_var, alpha)

    RMSItemList = []
    for y in range(len(RotList)):
        # PoseName, PoseNameExt, RotNum
        f_name, f_ext = files[y].split('.')
        RMSItem = RMSInfo(f_name, files[y], RotNumList[y])
        RMSItemList.append(RMSItem)

    # Delete obj1
   # cmd.delete('obj1')

    file_len = len(files)
    All_RMS = np.zeros((file_len, file_len))
    row_n = 0
    # For the designated lowest rotation, we cycle through
    for obj_x in RotList:
        # Calculates RMSList values and puts into list
        RMSList = RMS_Calc(obj_x, RotList)

        # Puts list into the corresponding row of the np array
        All_RMS[row_n, :] = RMSList

        row_n = row_n + 1

    # Change directory back to the main one, not the complex directory.
    print(cur_dir)
    os.chdir(cur_dir)

    GeneralSimCheck(RMSItemList, All_RMS, cur_dir, "RMS", float(cutoff), "")

    File_write(RotStruct, All_RMS, "RMS", "", PDB_code, cur_dir)
    # "YlGnBu"  "Greens_r"

    os.chdir(cur_dir)
    CreateHeatMap("", All_RMS, RotStruct, "RMS", "rocket_r", PDB_code, cur_dir)

    # Remove the final rotation files at the end (others were previously removed)
    os.chdir(os.path.join(cur_dir, "PDBComplex"))
    for item in RMSItemList:
     #   os.remove(item.PoseName + '.pdb')
      #  Rot_beginning = item.PoseName.rsplit('_', 1)[0]
        os.remove(item.PoseName + '_' + item.RotNum + 'rot.pdb')
        os.remove(item.PoseName + '_UNK' + item.RotNum + '.pdb')
    print("RMS Finished Processing.")

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

  #  print("RMS finished processing.")
