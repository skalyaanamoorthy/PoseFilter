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

######################################################################################################################
# Takes a parent directory and cycles through each of the folders in that directory.

######################################################################################################################
#CreateHeatMap("Fingerprint", All_Fingerprint, Listoflig, fprint, "Greens", PDB_code)
def CreateHeatMap(variety, mol_data, Ligands, type, map_type, PDB_code, working_dir):
    path = os.path.join(working_dir, variety, type)
    os.chdir(path)

    sns.set()
    mask = np.zeros_like(mol_data)
    mask[np.triu_indices_from(mask)] = True
    with sns.axes_style("white"):
        ax = sns.heatmap(mol_data, cmap=map_type, xticklabels= Ligands, yticklabels= Ligands, cbar=False, mask=mask, square=True, annot=True, annot_kws={'size':8})
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
        print(RMSArr[m])

    # The number of rotations directly correlate to the olig_num value (tetra = 4 - 1 -> 3, tri -> 2, di -> 1)
    print(i)

    return RMSArr

######################################################################################################################
# This function generates a pair_fit string that can then be called multiple times, depending on the type of olig
# Takes in i, and Pair_List, which is a list of numbers that should be added to the expression olig # times

def olig_string(i, i_List, ref_List):
    global chain, olig, r_min, r_max, r_alpha
    # if dimer x = 0, 1
    # if trimer x = 0, 1, 2
    print("olig string nums: ")
    print("olig: " + str(olig))
    print("res min: " + str(r_min))
    print("res max: " + str(r_max))
    print("chain: ")
    print(chain)
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
    cmd.do('set retain_order,1')
    cmd.do(rot_str)

    # Save the rotation and the ligand
    cmd.save(i + '_' + str(num) + 'rot.pdb', i)

    # Save just the ligand
    cmd.save(i + '_UNK' + str(num) + '.pdb', i + ' and resn ' + UNK)
    print("UNK Rotation: " + i + '_UNK' + str(num) + '.pdb')

    # Calculates the RMS between the UNK and the obj1 UNK reference
    rms_val = cmd.rms_cur('obj1 and resn ' + UNK, i + ' and resn ' + UNK)
    print("Rotation " + str(num) + " rms val: " + str(rms_val))

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
        print("RMS ARRAY")
        for x in rmsArr:
            print(x)

        # Removes the necessary files, calculates the lowest RMS
        New_rot, x_val, Rot_S = distComp(rmsArr, i)
        RotStruct.append(Rot_S)
        RotList.append(New_rot)
        RotNumList.append(x_val)

    return RotList, RotStruct, RotNumList

######################################################################################################################
# distComp functions choose the closest rotation file, and then deletes the files that are further away
# (those files are not of concern)

# Takes the RMSArray and then outputs the lowest rms value; ONE value
def distComp(RMSArr, i):
    global olig
    # Find the min index and then loop through to delete the other files
    min_index = RMSArr.index(min(RMSArr))
    print("Min index: " + str(min_index))
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
            pass
        # for now
          #  os.remove(i + '_UNK' + str(x) + '.pdb')

    # Returns a string of a name that has the lowest RMS
    return Lowest_RMS, x_val, Rot_S

####################################################################################
# Takes in a list with some extra files that may not be ligands, only outputs the files that are of the same format
'''
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
'''


######################################################################################################################

# Since this is for the olig function we want to make complexes with the ligand and protein files
def MakeComplex(savename, ligand, protein):
  #  os.chdir(dir)
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
    final_files = []
    for x in files:
        if "rot" in x:
            pass
        elif "UNK" in x:
            pass
        elif "COM" in x:
            pass
        else:
            print(x)
            final_files.append(x.strip('.pdb'))
            cmd.load(x)
    return final_files

######################################################################################################################

# Checks for similarities from a row. Saves the files to corresponding: Similar/Unique directories
# Type is RMS or FP
def GeneralSimCheck(RotList, RotStruct, RotNumList, Total_Array, working_dir, Dir_Type, cutoff, Add_Type):
    print("RotStruct len: " + str(len(RotStruct)))
    print("RotList len: " + str(len(RotList)))
   # cutoff = 2

    for x in RotStruct:
        print("RotStruct: " + str(x))
    # Puts the Rotation structure into this new variable
   # Rot_type = [RotStruct, [0]*len(RotStruct)] Puts the UNK structure into this new variable

    RotListPDB = []
    for y in RotList:
        RotListPDB.append(y + '.pdb')

    Rot_type = [RotListPDB, [0]*len(RotStruct)]

    Similar_path, Unique_path = CreateDirs(Dir_Type, Add_Type)

    similar_file_path = os.path.join(Similar_path, 'Similar.csv')

    sim_file = open(similar_file_path, 'w')
    sim_file.write('Reference, Structure, ' + Dir_Type + ', Rotation #\n')

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

                if Dir_Type == "RMS":
                    Is_similar = (Total_Array[x][y] <= cutoff)
                else:
                    Is_similar = (Total_Array[x][y] >= cutoff)

                if Is_similar:
               # if Total_Array[x][y] <= cutoff:
                    # Structure
                    Rot_type[1][y] = 'u'
                   # print(Rot_type[0][y] + " is similar")
                    # Need to make one unique out of the set
                    # Reference
                    Rot_type[1][x] = 's'

                    # Could now write to a file to indicate which files are similar
                    sim_file.write(str(Rot_type[0][x]) + ', ' + str(Rot_type[0][y]) + ', ' +
                                   str(Total_Array[x][y]) + ', ' + RotNumList[y] + '\n')

                # Not similar; pass for now
                else:
                    pass

            y = y - 1
        x = x - 1
    # Use this to create two new lists, sList, uList

    # Go through the structures again and add to the proper directories
    for w in range(len(RotStruct)):
        # Name
        rot = str(Rot_type[0][w])
        print("rot:" + rot)

        if Rot_type[1][w] == 's':
            # Move from dir to dir

            # Want a new path that has the directory similar on it
            path_1 = os.path.join(working_dir, rot)
            print(path_1)
            path_2 = os.path.join(Similar_path, rot)
            print(path_2)
            shutil.copy(path_1, path_2)

        else:
            path_1 = os.path.join(working_dir, rot)
            print(path_1)
            path_2 = os.path.join(Unique_path, rot)
            print(path_2)
            shutil.copy(path_1, path_2)

    sim_file.close()

###################################################################################################
# Checks dirs and creates the appropriate folders for moving the RMS/Fingeprint stuff into
def CreateDirs(Dir_Type, Add_Type):
    # In the case of RMS
    if Add_Type == "":
        print("Just make similar/unique in DIR")
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
        print("fingerprint special case")
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
               # os.makedirs(Dir_Type + '/' + Add_Type + '/' + 'Similar')
               # os.makedirs(Dir_Type + '/' + Add_Type + '/' + 'Unique')
        # Switch these
        Unique_path = os.path.join(Dir_Type, Add_Type, 'Similar')
        Similar_path= os.path.join(Dir_Type, Add_Type, 'Unique')

    return Similar_path, Unique_path

###################################################################################################
# Takes the np array and then writes to file
# Takes in a LigandList, Array to write to the file, and type of fingerprint
# Needs to make a folder of "Type" if that does not exist" -> in the case of the fingerprint file make a directory of Type info too
def File_write(LigList, Total_Array, Type, TypeInfo, PDB_code, working_dir):
    path = os.path.join(working_dir, Type, TypeInfo)
    os.chdir(path)
    if TypeInfo is not "":
        TypeInfo = TypeInfo + '_'

    if PDB_code is not "":
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
            f2.write(str(item) + ', ')
        f2.write('\n')
        row_t = row_t + 1
   # Make sure the dir is correct

    os.chdir(working_dir)

#####################################################################################################
# Checks if individual ligands or a protein + ligand complex
# proteintoliglist is a function that gets the ligand list from the directory
# In general, want to have the directory, specified by pfile or by the dir itself
# Want a list of ligand files or complex files in the directory
# Need shared lists of the complexes, or just the ligands so we don't need to regenerate stuff
# IComplexList = []
# ILigandList = []
# Depending on the user input, one of these will be created.

#def LigandtoComplex(pfile, LigandId):
    # Change dir
#    files = DirSearch(LigandId)
    # for each of those ligands, make a ligand + protein complex


def DirSearch(keyword):
 #   file_ext = filewithext.rsplit('.', 1)[1]
    all_files = []
    query = glob.glob("*" + keyword + "*")
    print("*" + keyword + "*")
    print(query)

    for x in query:
        print("query:")
        print(query)
        if "rot" in x:
            pass
        elif "UNK" in x:
            pass
        elif keyword in x:
            all_files.append(x)
        else:
            pass

    return all_files


# Returns the appropriate ligand list for the circumstances
def ListGenerate(AnalysisType, TabInput):
    # Tab 1 means protein and listoflig
    pass
  #  if TabInput == "Ligand":
        # Take the protein file as ext
   #     protein_name = os.path.basename(pfile)
   #     file_ext = protein_name.rsplit('.', 1)[1]

  #  else:

######################################################################################################################
# Directory must have one protein file (specified) and multiple ligand files
# Takes the protein .pdb/.pdbqt file and then finds the appropriate ligand files of that extension
# outputs the ligand files in a list format.

#Proteintoliglist
def LigandtoComplex(pfile, keyword):
    working_dir = os.path.dirname(pfile)
    os.chdir(working_dir)
    print(working_dir)

    protein_name = os.path.basename(pfile)
    file_ext = protein_name.rsplit('.', 1)[1]
    all_files = DirSearch(keyword)

    #   LigandList = ListtoLig(all_files)
    #  print(LigandList)
    Saved_Complexes = [""] * len(all_files)
    # Goes through each ligand to make a complex
    SNum = 0

    # Make the protein + ligand complexes
    for ligand in all_files:
        print("ligand: " + ligand)
        lig_name = ligand.rsplit('.', 1)[0]
        print("Make complex: ")
        MakeComplex(lig_name + '_complex', ligand, protein_name)

        # These are the names of the complexes that will be iterated through the program
        Saved_Complexes[SNum] = lig_name + '_complex' + ".pdb"
        SNum += 1
    print("Saved complexes: ")
    # Saved complexes are protein + ligand .pdb named
    return Saved_Complexes, all_files[0]

##############################################################################################################

import re

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)