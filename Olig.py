from __future__ import print_function
import pymol
from pymol import cmd
from pymol import stored
from pymol import selector
#import centroid
import os
import re
import fileinput
import shutil
import glob
import numpy as np
import sys
from .Fingerprint import CreateHeatMap
#from ResTest import PDB_Info
from statistics import mode
#from center_of_mass import com
#import center_of_mass.py
from collections import Counter



# All_RMS is a np olig_num x olig_num array; RMS_rot is a list of the olig rotations that are lowest RMS value
global All_RMS, res_min, res_max, olig_num, RotList, RotStruct, RotNumList

# Will be used as a switch- when it is on then we have lots of output to test, off we do not have any
global Testing
Testing = 0

RotList = []
RotNumList = []
RotStruct = []
olig_num = 0
global working_dir, toAlph
#toAlph = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']



from pymol import cmd

# we have a .pdb file that was chosen. Need to get the .pdb ID from it, as well as the Res_min and Res_max numbers
# the chain order, and the chain labels

# returns minval, maxval, (List of chainstrings)

def PDB_Info(file):

    # Need to load the file so we can do things with it
    cmd.load(file, 'obj')
    # The file will be a docked protein .pdb
    chains = ChainOrder(file)
    ChainNum = len(chains)

    # Contains sublists, each with [ChainNum, [ListofChains]], the list of chains which corresponds to the chainNum
    ListofChains = []

    # Chain order; go through .pdb and get the order
    # For each chain, check the number of res; keep a list of the chains
    stored.residues = []

    # Keeps a list of the identical chains
    for chain in chains:
        label = 'chain ' + chain
        #print(label)

        # Res_min, res_max
        stored.residues = []
        cmd.iterate(selector.process(label), 'stored.residues.append(resv)')
        val = float(stored.residues[-1] - stored.residues[0])
        minval = stored.residues[0]
        maxval = stored.residues[-1]
       # print("val: " + str(val))

        if len(ListofChains) == 0:
            # Add the first chain to the list
            ListofChains.append([[minval, maxval], [label]])

        # The list is not empty, but the chain num corresponds
        else:
            # Loop through the ListofChains list
            x = 0
            Added = 0
            while x in range(len(ListofChains)):
                chaindiff = float(ListofChains[x][0][1]) - float(ListofChains[x][0][0])
                if chaindiff == val:
                    # Add val to the list
                    ListofChains[x][1].append(label)
                    Added = 1
                    break

                x += 1

                # Had not been added yet, need to add a new entry
                if Added == 0:
                    ListofChains.append([[minval, maxval], [label]])

    max = 0
    minval = 0
    maxval = 0
   # print(len(ListofChains))
    chainstrings = []
    for item in ListofChains:
        if Testing:
            print("item 1: " + str(item[1]))
            print("item 0: " + str(item[0]))

        if len(item[1]) > max:
            minval = item[0][0]
            maxval = item[0][1]
            max = maxval-minval

            if Testing:
                print("min: " + str(minval) + " max: " + str(maxval))

            chainstrings = item[1]

    if Testing:
        print("max " + str(max))
        print(str(maxval))
        print(str(minval))
    newchains =[]
    for x in chainstrings:
        newchains.append(x.split()[1])

    if Testing:
        print(chainstrings)
        print(newchains)

    # newchains before
    return minval, maxval, chains

# The chainstrings list has all of the chains that are part of the oligomer,
# IN the order that they should be in!!
# Max has the difference -> save the min an max as well though to return.
   # for i in chainstrings:
      #  print("chainstrings: " + i)

######################################################################################################################
# Gets the chain order Prints a list of chains in the correct order
def ChainOrder(energy_file):
    atomList = []
    f = open(energy_file, 'r')
    lines = f.readlines()
    for line in lines:
        #print(line)
        if len(line.split()) >= 6 and len(line.split()[2]) <= 3:
            if "ATOM" in line.split()[0]:
              #  print(line)
                atom = line.split()[4]
                atom = atom[0]
                if atom not in atomList:
                #    print(line)
                    atomList.append(atom)
        else:
            pass
  #  for x in atomList:
  #      print("atomlist: " + str(x))
    return atomList
#   return res_min, res_max


# Should make this based on the olig_size
#All_RMS = []

# User inputs the type of oligomer
# Function takes in if di, tri, tetramer, as well as the lowest energy file (name or .pdb, we can strip)

######################################################################################################################

def olig_rot(i):
    # olig_num is important here
    global olig_num
    iList = [0]*olig_num
    ref_List = [0]*olig_num
    RMSArr = [0]*olig_num

    # Create the initial lists, based on the oligomer's size.
    for n in range(olig_num):
        iList[n] = n
        ref_List[n] = n

    # Calls the rotation m number of times
    mod_add = 0
    for m in range(olig_num):
        # First, goes through to modify iList; adds the m value and then takes mod olig_num
        for x in range(len(iList)):
            iList[x] = (iList[x] + mod_add) % olig_num

        mod_add = 1

         # Rotation is called

        #printing
       # for p in range(len(iList)):
        #    print('iList el: ' + str(iList[p]))
         #   print('refList el: ' + str(ref_List[p]))

        RMSArr[m] = rotation(olig_string(i, iList, ref_List), i, m)
        print(RMSArr[m])

    # The number of rotations directly correlate to the olig_num value (tetra = 4 - 1 -> 3, tri -> 2, di -> 1)
    print(i)
   # for x in RMSArr:
     #   print("RMS Arr: " + str(x))

    return RMSArr


######################################################################################################################
# distComp functions choose the closest rotation file, and then deletes the files that are further away
# (those files are not of concern)

# Takes the RMSArray and then outputs the lowest rms value
def distComp(RMSArr, i):
    global RotStruct
    # Find the min index and then loop through to delete the other files
    min_index = RMSArr.index(min(RMSArr))

    for x in range(olig_num):

        # If we have reached the index, then load!
        if x == min_index:
          #  cmd.save(i + '_' + str(x) + 'rotation.pdb', i + ' and not resn PSD')
          #  cmd.save(i + '_UNK' + str(x) + '.pdb', i + ' and resn UNK')
            cmd.load(i + '_UNK' + str(x) + '.pdb', i + '_UNK' + str(x) + '.pdb', )

            Lowest_RMS = i + '_UNK' + str(x)
            x_val = str(x)
            RotStruct.append(i + '_' + str(x) + 'rotation.pdb')

          #  Rtn_Str.append(Lowest_RMS + '.pdb')

        else:
           # os.remove(i + '_' + str(x) + 'rotation.pdb')

           # REMOVE later!! Saves the rotation
          #  cmd.save(i + '_' + str(x) + 'rotation.pdb', i + ' and not resn PSD')
            os.remove(i + '_UNK' + str(x) + '.pdb')

    # Returns a string of a name that has the lowest RMS
    return Lowest_RMS, x_val

######################################################################################################################

# Takes a number input and outputs the appropriate letter
def NumtoAlph(num):
    # Take these from the chain...
    Letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
   # Letters = ['B', 'A', 'C', 'D' ]
    return Letters[num]

######################################################################################################################
# This function generates a pair_fit string that can then be called multiple times, depending on the type of olig
# Takes in i, and Pair_List, which is a list of numbers that should be added to the expression olig # times

def olig_string(i, i_List, ref_List):
    global toAlph
    # if dimer x = 0, 1
    # if trimer x = 0, 1, 2
    total_string = "pair_fit /"
    for x in range(olig_num):

        sub_string = i + "//" + toAlph[i_List[x]] + "/" + str(res_min) + "-" + str(res_max) +\
                     "/N+C+O+CA, /obj1//" + toAlph[ref_List[x]] + "/" + str(res_min) + "-" + str(res_max) + "/N+C+O+CA"
        total_string = total_string + sub_string

        if x < (olig_num - 1):
            total_string = total_string + ', /'

    return total_string


######################################################################################################################

# Does rotation, makes the rot structures
def rotation(rot_str, i, num):
    cmd.do('set retain_order,1')
    cmd.do(rot_str)
    # cmd.save(i + '_' + str(num) + 'rotation.pdb', i + ' and not(resn UNK or resn PSD)')
    cmd.save(i + '_' + str(num) + 'rotation.pdb', i + ' and not resn PSD')
  # cmd.save(i + '_COM' + str(num) + '.pdb', i + ' and resn PSD')
    cmd.save(i + '_UNK' + str(num) + '.pdb', i + ' and resn UNK')

    # dist_val = cmd.distance('d' + num, i + ' and resn PSD', 'ref_lig')

    rms_val = cmd.rms_cur('obj1 and resn UNK', i + ' and resn UNK')
    print("rms val: " + str(rms_val))

    # All of the RMS_vals will be compared later
    return rms_val


#####################################################################################################################

# Takes an optional PDB code and uses that for labelling purposes
def Olig_MainLoop(energy_files, PDB_code):

    global All_RMS, RotInfo, olig_num
    print("energy files: ")
    print(energy_files)

    # First, load the first file as a reference
    cmd.load(energy_files[0], 'obj1')
    print(energy_files[0])

    # Resets the olig num here
  #  ChainCheck('obj1')
  #  print(olig_num)

    # Load files; disregard any for now that would have been output from before
    for x in energy_files:
        if "rotation" in x:
            pass
        elif "UNK" in x:
            pass
        elif "COM" in x:
            pass
        else:
            print(x)
            cmd.load(x)

    # Gets all objects of this type
    myobjects = cmd.get_object_list()

    # Cycle through the objects in this list from position 1 to the end.
    global RotList, RotStruct, RotNumList

    obj_num = 0

    # For each object
    for i in myobjects[2:]:

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
        New_rot, x_val = distComp(rmsArr, i)
        RotList.append(New_rot)
        RotNumList.append(x_val)

    row_n = 0
    # For the designated lowest rotation, we cycle through
    for obj_x in RotList:

        # Calculates RMSList values and puts into list
        RMSList = RMS_Calc(obj_x, RotList)

        # Puts list into the corresponding row of the np array
        All_RMS[row_n, :] = RMSList

        # Clean up the extra UNK files now
        os.remove(obj_x + '.pdb')

        row_n = row_n + 1

    # Now need to analyze the RMS values and then print them.
    RMS_analysis()
    CreateHeatMap(All_RMS, RotList, "RMS", "Greens_r")

    SimCheck()

######################################################################################################################

# Takes the np array and then writes to file
def RMS_analysis():
    # store all of the rms values
    global All_RMS, RotList, RotNumList

    f2 = open('RMS.csv', 'w')
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

######################################################################################################################

# Checks for similarities from a row. Saves the files to corresponding: Similar/Unique directories
def SimCheck():
    global RotList, RotStruct, RotNumList, All_RMS

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
                if All_RMS[x][y] < 2:
                    # Structure
                    Rot_type[1][y] = 'u'
                    print(Rot_type[0][y] + "is similar")
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
    simfinder = RotStruct
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




######################################################################################################################
# Directory must have one protein file (specified) and multiple ligand files
# Takes the protein .pdb/.pdbqt file and then finds the appropriate ligand files of that extension

def ProteintoLigList(pfile):
    global working_dir
    working_dir = os.path.dirname(pfile)
    os.chdir(working_dir)
    print(working_dir)

    # pfile just want the extention as well as the name + extension

    protein_name = os.path.basename(pfile)
    file_ext = protein_name.rsplit('.', 1)[1]
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
    MakeComplex(working_dir, testlig.rsplit('.', 1)[0], testlig, protein_name)

    for ligand in LigandList:
        print("ligand: " + ligand)
        lig_name = ligand.rsplit('.', 1)[0]
        # #MakeComplex(dir, "Justlig", "lig2.pdbqt", "protein.pdbqt")
        MakeComplex(working_dir, lig_name, ligand, protein_name)
        #   print("loading complex")
        #   cmd.load(lig_name + '.pdb')

        # These are the names of the complexes that will be iterated through the program
        Saved_Complexes[SNum] = lig_name + ".pdb"
        SNum += 1
    print("Saved complexes: ")
    return Saved_Complexes


######################################################################################################################

def InputFileDir(pfile, PDBCode):
    cmd.reinitialize()
    Saved_Complexes = ProteintoLigList(pfile)
    Saved_Complexes.sort()
    OligWrapper(Saved_Complexes, PDBCode)
   # return Saved_Complexes

####################################################################################

# Since this is for the olig function we want to make complexes with the ligand and protein files
def MakeComplex(dir, savename, ligand, protein):
    os.chdir(dir)
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


####################################################################################
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

# Changing to incorporate selecting a protein.pdb/pdbqt file first
def OligWrapper(PDBfiles, PDB_code):

    global RotList, RotNumList, RotStruct, res_min, res_max, toAlph, olig_num
    RotList = []
    RotNumList = []
    RotStruct = []
    # First, reinitialize the PyMOL window
    cmd.reinitialize()

   # energyfiles = InputFileDir(protein_file)



    # Call the function to get the files from the directory.. but now it would be based on the protein.pdb/pdbqt file

    # Not super necessary but gets the first energy file

    # Need to redefine some things
    #print("splitting ext")
   # print()

    energy_file = PDBfiles[0] #.rsplit(".", 1)[0]
    print("energy file: " + energy_file)
    PDB_len = len(PDBfiles)

    # Initialize the global numpy array
    global All_RMS
    All_RMS = np.zeros((PDB_len, PDB_len))


   # print("energy file: " + energy_file)

    # Assigns the min and max res, as well as the chain to use
    res_min, res_max, toAlph = PDB_Info(energy_file)
    print(res_min)
    print(res_max)
    olig_num = len(toAlph)
    if Testing:
        print(olig_num)

    Olig_MainLoop(PDBfiles, PDB_code)

    del(All_RMS)
  #  del(RotList)
    print("Finished processing.")
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
