import pymol
from pymol import cmd
from pymol import stored
from pymol import selector
import copy
import os
import re
import fileinput
import shutil
import glob
import numpy as np
import sys
from statistics import mode
from collections import Counter
global Testing
Testing = 0
# Will be used as a switch- when it is on then we have lots of output to test, off we do not have any
global o_num, maxofres, minofres, ListChains
ListChains = []
maxofres = -1
minofres = -1
# need to define these here and call them from here too.

######################################################################################################################
# Gets the chain order Prints a list of chains in the correct order
def ChainOrder(energy_file):
    atomList = []
    f = open(energy_file, 'r')
    lines = f.readlines()
    for line in lines:
        if len(line.split()) >= 6 and len(line.split()[2]) <= 3:
            if "ATOM" in line.split()[0]:
                atom = line.split()[4]
                atom = atom[0]
                if atom not in atomList:
                    atomList.append(atom)
        else:
            pass
    return atomList


def CheckIdentical(totalchains):
   # cmd.reinitialize()
    chainresns = []
    for chain in totalchains:
        label = 'chain ' + chain

        stored.residues = []
        cmd.iterate(selector.process(label), 'stored.residues.append(resn)')
        chainresns.append(stored.residues)

    dupchain = []
    for x in chainresns:
        dupchain.append(chainresns.count(x))
       # print(chainresns)
        print(chainresns.count(x))

    max_index = dupchain.index(max(dupchain))
    finalchains = []
    finalchains.append(totalchains[max_index])
    indexc = 0
    for x in range(len(chainresns)):
        if x != max_index and chainresns[x] == chainresns[max_index]:
            finalchains.append(totalchains[x])
    print(finalchains)
    return finalchains


######################################################################################################################
# Protein file is chosen. Need to get the Res_min and Res_max numbers, chain labels and chain order
# returns minval, maxval, (List of chainstrings)

def PDB_Info(file, nonidentical):
    # Need to load the file so we can do things with it
  #  cmd.do('set retain_order,1')
    cmd.load(file, 'obj')
    # The file will be a docked protein .pdb
    chains = ChainOrder(file)
    ChainNum = len(chains)

    # Contains sublists, each with [ChainNum, [ListofChains]], the list of chains which corresponds to the chainNum
    ListofChains = []

    # Chain order; go through .pdb and get the order
    # For each chain, check the number of res; keep a list of the chains
    stored.residues = []
    print(chains)
    ichains = CheckIdentical(chains)

    if nonidentical:
        chains = ichains

    print("identicalchains")
    print(ichains)

    # Keeps a list of the identical chains
    for chain in chains:
        label = 'chain ' + chain
        print(label + 'chain')

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

    chainstrings = []
    for item in ListofChains:

        if len(item[1]) > max:
            minval = item[0][0]
            maxval = item[0][1]
            max = maxval-minval

            if Testing:
                print("min: " + str(minval) + " max: " + str(maxval))

            chainstrings = item[1]

    newchains =[]
    for x in chainstrings:
        newchains.append(x.split()[1])

    cmd.delete('obj')
    global o_num, minofres, maxofres, ListChains
    o_num = len(chains)
    minofres = minval
    maxofres = maxval
    ListChains = chains
    return minval, maxval, chains

# Does not redo the operation; saves the past operation
def PDBInfo_Wrapper(file, nonidentical):
    global ListChains, minofres, maxofres
    minval, maxval, chains = PDB_Info(file, nonidentical)
    minofres = minval
    maxofres = maxval
    ListChains = chains
    return minval, maxval, chains

#infoc = (PDB_Info('/home/justine/PycharmProjects/PoseFilter/Examples/6ewp/6ewp.pdbqt', 1))
#print(infoc)
#CheckIdentical(infoc[2])