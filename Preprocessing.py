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

######################################################################################################################
# Protein file is chosen. Need to get the Res_min and Res_max numbers, chain labels and chain order
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

        # Res_min, res_max
        stored.residues = []
        cmd.iterate(selector.process(label), 'stored.residues.append(resv)')
        print("stored residues")
        print(stored.residues)
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

    cmd.delete('obj')
    global o_num, minofres, maxofres, ListChains
    o_num = len(chains)
    minofres = minval
    maxofres = maxval
    ListChains = chains
    return minval, maxval, chains

# Does not redo the operation; saves the past operation
def PDBInfo_Wrapper(file):
    global ListChains, minofres, maxofres
 #   if minofres > -1:
 #       return minofres, maxofres, ListChains
    minval, maxval, chains = PDB_Info(file)
    minofres = minval
    maxofres = maxval
    ListChains = chains
    return minval, maxval, chains
