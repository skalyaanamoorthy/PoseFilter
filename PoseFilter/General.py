
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

        RMSArr[m] = rotation(olig_string(i, iList, ref_List), i, m)
        print(RMSArr[m])

    # The number of rotations directly correlate to the olig_num value (tetra = 4 - 1 -> 3, tri -> 2, di -> 1)
    print(i)

    return RMSArr


######################################################################################################################
# This function generates a pair_fit string that can then be called multiple times, depending on the type of olig
# Takes in i, and Pair_List, which is a list of numbers that should be added to the expression olig # times

def olig_string(i, i_List, ref_List):
    global toAlph, olig_num, res_min, res_max
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

    rms_val = cmd.rms_cur('obj1 and resn UNK', i + ' and resn UNK')
    print("rms val: " + str(rms_val))

    # All of the RMS_vals will be compared later
    return rms_val

#####################################################################################################################
# Generate rotation list from the already labelled 'obj1' structure
def GenerateRotList(objects):
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

# Takes the RMSArray and then outputs the lowest rms value, ONE value
def distComp(RMSArr, i):
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
            Rot_S = i + '_' + str(x) + 'rotation.pdb'

          #  Rtn_Str.append(Lowest_RMS + '.pdb')

        else:
           # os.remove(i + '_' + str(x) + 'rotation.pdb')

           # REMOVE later!! Saves the rotation
          #  cmd.save(i + '_' + str(x) + 'rotation.pdb', i + ' and not resn PSD')
            os.remove(i + '_UNK' + str(x) + '.pdb')

    # Returns a string of a name that has the lowest RMS
    return Lowest_RMS, x_val, Rot_S

######################################################################################################################