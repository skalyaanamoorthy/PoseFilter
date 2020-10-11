import pymol
from pymol import stored
from pymol import cmd
from pymol import selector
#import centroid
import os
import re
import fileinput
import shutil
import glob
import numpy as np

import sys
from statistics import mode
from collections import Counter
import matplotlib.pyplot as plt
from pylab import savefig
import seaborn as sns

#global olig
global Testing
# Fix these.
global chain, r_min, r_max, olig, UNK, r_alpha
Testing = 0
# Will be used as a switch- when it is on then we have lots of output to test, off we do not have any

class RMSInfo:
  def __init__(self, PoseName, PoseNameExt, RotNum):
    #RotStruct is pose1_complex_0rot.pdb
    #RotList is pose1_complex_UNK0
    self.PoseName = PoseName
    self.PoseNameExt = PoseNameExt
    self.RotNum = RotNum

######################################################################################################################
#CreateHeatMap("Fingerprint", All_Fingerprint, Listoflig, fprint, "Greens", PDB_code)
def CreateHeatMap(variety, mol_data, Ligands, type, map_type, PDB_code, working_dir):
    path = os.path.join(working_dir, variety, type)
    os.chdir(path)

    mol_data_round = np.around(mol_data, decimals=2)

    sns.set()
    mask = np.zeros_like(mol_data_round)
    mask[np.triu_indices_from(mask,1)] = True
    with sns.axes_style("white"):
        ax = sns.heatmap(mol_data_round, cmap=map_type, xticklabels= Ligands, yticklabels= Ligands, cbar=False, mask=mask, square=True, annot=True, annot_kws={'size':8})
        Map_title = type + ' ' + PDB_code
        ax.set(title = Map_title)
        figure = ax.get_figure()
       # figure.tight_layout()
        figure.savefig(type + ".png", bbox_inches = "tight")
    plt.clf()
    os.chdir(working_dir)
  #  cmd.reinitialize()

######################################################################################################################

def olig_rot(i):
    # olig_num is important here
    global olig
    iList = [0]*olig
    ref_List = [0]*olig
    RMSArr = [0]*olig

    # Create the initial lists, based on the oligomer's size.
    for n in range(olig):
        iList[n] = n
        ref_List[n] = n

    # Calls the rotation m number of times
    mod_add = 0
    for m in range(olig):
        # First, goes through to modify iList; adds the m value and then takes mod olig
        for x in range(len(iList)):
            iList[x] = (iList[x] + mod_add) % olig

        mod_add = 1

         # Rotation is called
        RMSArr[m] = rotation(olig_string(i, iList, ref_List), i, m)
    #    print(RMSArr[m])

    # The number of rotations directly correlate to the olig_num value (tetra = 4 - 1 -> 3, tri -> 2, di -> 1)
 #   print(i)

    return RMSArr

######################################################################################################################
# This function generates a pair_fit string that can then be called multiple times, depending on the type of olig
# Takes in i, and Pair_List, which is a list of numbers that should be added to the expression olig # times

def olig_string(i, i_List, ref_List):
    global chain, olig, r_min, r_max, r_alpha
    # if dimer x = 0, 1
    # if trimer x = 0, 1, 2
   # print("olig string nums: ")
 #   print("olig: " + str(olig))
 #   print("res min: " + str(r_min))
 #   print("res max: " + str(r_max))
 #   print("chain: ")
 #   print(chain)
    #r_max = "60"

    if r_alpha:
        fit_string = "/CA"
    else:
        fit_string = "/N+C+O+CA"

    total_string = "pair_fit /"
    for x in range(olig):

        sub_string = i + "//" + chain[i_List[x]] + "/" + str(r_min) + "-" + str(r_max) +\
                     fit_string + ", /obj1//" + chain[ref_List[x]] + "/" + str(r_min) + "-" + str(r_max) + fit_string
        total_string = total_string + sub_string

        if x < (olig - 1):
            total_string = total_string + ', /'

    return total_string

######################################################################################################################

# Does rotation, makes the rot structures
def rotation(rot_str, i, num):
    global UNK
 #   cmd.do('set retain_order,1')
    cmd.do(rot_str)

    # Save the rotation and the ligand
    cmd.save(i + '_' + str(num) + 'rot.pdb', i)

    # Save just the ligand
    cmd.save(i + '_UNK' + str(num) + '.pdb', i + ' and resn ' + UNK)
    #  print("UNK Rotation: " + i + '_UNK' + str(num) + '.pdb')

    # Calculates the RMS between the UNK and the obj1 UNK reference
    rms_val = cmd.rms_cur('obj1 and resn ' + UNK, i + ' and resn ' + UNK)

    return rms_val

#####################################################################################################################
# Generate rotation list from the already labelled 'obj1' structure
def GenerateRotList(objects, UNK_var, alpha):
    # Get the variables from the Preprocessing file
    #from .Olig import olig_num, toAlph, res_min, res_max
    from .Preprocessing import minofres, maxofres, ListChains, o_num
    global r_min, r_max, chain, olig, UNK, r_alpha
    r_alpha = alpha
    olig = o_num
    r_min = minofres
    r_max = maxofres
    chain = ListChains
    UNK = UNK_var

    # For each object
    RotStruct = []
    RotList = []
    RotNumList = []
    for i in objects:
        # Can probably remove the tname stuff
        i = i.strip()

        # Does the rotation
        # olig_rot needs to give back info about files to delete, rms

        # Performs the rotation on one object, gives back the rmsArr which is a list of all of the rms values
        # we need to know for later which is the lowest and corresponds to which structure

        rmsArr = olig_rot(i)

        # Removes the necessary files, calculates the lowest RMS
        New_rot, x_val, Rot_S = distComp(rmsArr, i)
        RotStruct.append(Rot_S)
        RotList.append(New_rot)
        RotNumList.append(x_val)

    return RotList, RotStruct, RotNumList

######################################################################################################################
# distComp functions choose the closest rotation file, and then deletes the files that are further away
# (those files are not of concern)

# Takes the RMS 1D Array: outputs the Lowest_RMS value (float), x_val (string), and Rotation file name (string)
def distComp(RMSArr, i):

    global olig
  #  cmd.do('set retain_order,1')
    # Find the min index and then loop through to delete the other files
    min_index = RMSArr.index(min(RMSArr))
 #   print("Min index: " + str(min_index))
    for x in range(olig):
        # Remove all of the protein + ligand complexes (rotations); maybe would not even
        # Save those to begin with (rotation function)
        #os.remove(i + '_' + str(x) + 'rotation.pdb')

        # If we have reached the index, then load! Remove all of the other UNK values
        if x == min_index:
            # Loads into PyMOL
            cmd.load(i + '_UNK' + str(x) + '.pdb')

            Lowest_RMS = i + '_UNK' + str(x)
            x_val = str(x)
            Rot_S = i + '_' + str(x) + 'rot.pdb'

        else:
            os.remove(i + '_UNK' + str(x) + '.pdb')
            os.remove(i + '_' + str(x) + 'rot.pdb')

    # Returns a string of a name that has the lowest RMS
    return Lowest_RMS, x_val, Rot_S

######################################################################################################################

# Since this is for the olig function we want to make complexes with the ligand and protein files
def MakeComplex(savename, ligand, protein):
  #  os.chdir(dir)
   # cmd.do('set retain_order,1')
    cmd.load(ligand)
    cmd.load(protein)
    save_as = savename + ".pdb"
    ligand_s = ligand.rsplit(".", 1)[0]
    protein_s = protein.rsplit(".", 1)[0]
    cmd.save(save_as, ligand_s + ' or ' + protein_s)
  #  cmd.load(save_as)
    # returns the name of the newly made pdb file
    cmd.delete(ligand_s)
    cmd.delete(protein_s)
    return save_as

######################################################################################################################
# Takes a list of files and makes sure that they are not including saved words such as rotation, UNK, etc.
# Loads the appropriate files
def FilterFiles(files):
    # Load files; disregard any for now that would have been output from before
  #  cmd.do('set retain_order,1')
    final_files = []
    for x in files:
        if "rot" in x:
            pass
        elif "UNK" in x:
            pass
        elif "COM" in x:
            pass
        else:
          #  print(x)
            final_files.append(x.strip('.pdb'))
            cmd.load(x)
    return final_files

######################################################################################################################

# Checks for similarities from a row. Saves the files to corresponding: Similar/Unique directories
# Type is RMS or FP

def GeneralSimCheck(RMSItemList, Total_Array, working_dir, Dir_Type, cutoff, Add_Type):
    length = len(RMSItemList)
    RotListPDB = []
    for item in RMSItemList:
        RotListPDB.append(item.PoseNameExt)

    #RotListPDB
    #Rot_type = [RotListPDB, [0]*length]

    Similar_path, Unique_path = CreateDirs(Dir_Type, Add_Type)

    # In the fingerprint case, similar path is unique one
    similar_file_path = os.path.join(Similar_path, 'Similar.csv')

    sim_file = open(similar_file_path, 'w')
    sim_file.write('Reference, Structure, ' + Dir_Type + ', Rotation #\n')

    # The first row is the name, the second is the type: 's' for similar, nothing for unique
    # Uses the global numpy variable to check for similarities and then sorts the files.
    s_array = []
    u_array = []
    set_list = []
    our_dir = ""
  #  print(Total_Array)
   # x = length-1
    x = 0
    while x < length:
        y = x + 1
        cur_set = []
        while y < length:

            if Dir_Type == "RMS":
                Is_similar = (Total_Array[y][x] <= cutoff)
                our_dir = os.path.join(working_dir, "PDBComplex")
            else:
                Is_similar = (Total_Array[y][x] >= cutoff)
                our_dir = os.path.join(working_dir, "PDBLigand")
              #  print("similar val:")
              #  print(Total_Array[y][x])
              #  print(x)
              #  print(y)
              #  print(RMSItemList[x].PoseName)
              #  print(RMSItemList[y].PoseName)

            if Is_similar:
              #  print(RMSItemList[x].PoseName + "similar")
              #  print(RMSItemList[y].PoseName + "similar")
                # Append the reference structure to the similar list

                if RMSItemList[x] not in cur_set:
                    cur_set.append(RMSItemList[x])

                if RMSItemList[y] not in cur_set:
                    cur_set.append(RMSItemList[y])


                PoseName_x = RMSItemList[x].PoseName
                PoseName_y = RMSItemList[y].PoseName

                # Could now write to a file to indicate which files are similar
                sim_file.write(PoseName_x + ', ' + PoseName_y + ', ' +
                               str(Total_Array[x][y]) + ', ' + RMSItemList[y].RotNum + '\n')

            y = y + 1
        set_list.append(cur_set)
        x = x + 1

    for s_set in set_list:
       # print("new set")
        u_set = 0
        for i in s_set:
            if (i not in s_array) and not any(j in u_array for j in s_set):
         #       print("add to uarray (1):")
         #       print(i.PoseName)
                u_array.append(i)
                s_array.append(i)
                u_set = 1
            elif (u_set == 0) and (i not in s_array) and (i not in u_array):
                u_set = 1
                u_array.append(i)
                s_array.append(i)
            elif i not in s_array:
                s_array.append(i)

    # For each set in the master set (1,2,3) (3,4)
    for item in RMSItemList:
        if (item not in s_array) and (item not in u_array):
            u_array.append(item)
        #    print("add to uarray: " + item.PoseName)
        elif (item not in s_array):
            u_array.append(item)
        #    print("add to uarray: " + item.PoseName)

    # Use this to create two new lists, sList, uList

    if Dir_Type == "RMS":
        rot_add = '_rot'
    else:
        rot_add = ''

    # Go through the structures again and add to the proper directories
   # print(s_array)
   # print(u_array)
    for s_item in s_array:

        # Want a new path that has the directory similar on it
        path_1 = os.path.join(our_dir, s_item.PoseNameExt)
        path_2 = os.path.join(Similar_path, s_item.PoseNameExt)
        shutil.copy(path_1, path_2)

        path_3 = os.path.join(Similar_path, s_item.PoseNameExt)
        os.rename(path_2, path_3)

    for u_item in u_array:
        path_1 = os.path.join(our_dir, u_item.PoseNameExt)
        path_2 = os.path.join(Unique_path, u_item.PoseNameExt)
        shutil.copy(path_1, path_2)

        # Path to copy to
        path_3 = os.path.join(Unique_path, u_item.PoseNameExt)
        os.rename(path_2, path_3)

    sim_file.close()

###################################################################################################
# Checks dirs and creates the appropriate folders for moving the RMS/Fingeprint stuff into
def CreateDirs(Dir_Type, Add_Type):
    # In the case of RMS
    if Add_Type == "":
     #   print("Just make similar/unique in DIR")
        if not os.path.exists(Dir_Type):
            os.makedirs(Dir_Type + '/' + 'Similar')
            os.makedirs(Dir_Type + '/' + 'Unique')
        else:
            # Overwrite
            pass
          #  os.makedirs('Similar')
          #  os.makedirs('Unique')
        Similar_path = os.path.join(Dir_Type, 'Similar')
        Unique_path = os.path.join(Dir_Type, 'Unique')

    # In the case of Fingerprint
    else:
     #   print("fingerprint special case")
        if not os.path.exists(Dir_Type):
            os.makedirs(Dir_Type + '/' + Add_Type + '/' + 'Similar')
            os.makedirs(Dir_Type + '/' + Add_Type + '/' + 'Unique')
        else:
            # We have the main fingerprint path, but might not have the Add_Type one
          #  Add_Type_Path = os.path.join(Dir_Type)
            if not os.path.exists(Dir_Type + '/' + Add_Type):
                # Make Type folder and the similar and unique paths
                os.makedirs(Dir_Type + '/' + Add_Type + '/' + 'Similar')
                os.makedirs(Dir_Type + '/' + Add_Type + '/' + 'Unique')
            else:
                pass

        Similar_path = os.path.join(Dir_Type, Add_Type, 'Similar')
        Unique_path = os.path.join(Dir_Type, Add_Type, 'Unique')

    return Similar_path, Unique_path

###################################################################################################
# Takes the np array and then writes to file
# Takes in a LigandList, Array to write to the file, and type of fingerprint
# Needs to make a folder of "Type" if that does not exist" -> in the case of the fingerprint file make a directory of Type info too
def File_write(LigList, Total_Array, Type, TypeInfo, PDB_code, working_dir):
    path = os.path.join(working_dir, Type, TypeInfo)
    os.chdir(path)
    if TypeInfo != "":
        TypeInfo = TypeInfo + '_'

    if PDB_code != "":
        PDB_code = PDB_code + '_'

    FileName = TypeInfo + PDB_code + Type + '.csv'
    f2 = open(FileName, 'w')
    f2.write(Type + ' values'+'\n')
    f2.write(' ,')

    for file in LigList:
        f2.write(file + ', ')

    f2.write('\n')
    row_t = 0
    for row in Total_Array:
        f2.write(LigList[row_t] + ', ')
        for item in row:
            # Round to two decimal places
            f2.write(str(round(item, 2)) + ', ')
        f2.write('\n')
        row_t = row_t + 1
   # Make sure the dir is correct

    os.chdir(working_dir)

#####################################################################################################

def DirSearch(keyword):
 #   file_ext = filewithext.rsplit('.', 1)[1]
    all_files = []
    query = glob.glob("*" + keyword + "*")
   # print("*" + keyword + "*")
   # print(query)

    for x in query:
     #   print("query:")
     #   print(query)
        if "rot" in x:
            pass
        elif "UNK" in x:
            pass
        elif keyword in x:
            all_files.append(x)
        else:
            pass

    return all_files

##############################################################################################################

import re

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

#os.chdir('/home/justine/PycharmProjects/PoseFilter/Examples/Trimer_Example')
#objects = 'pose1.pdb, pose2.pdb, pose3.pdb'
#GenerateRotList(objects, "UNK", 0)