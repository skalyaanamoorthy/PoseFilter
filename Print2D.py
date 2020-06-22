import os
from pymol import cmd
from pymol import stored
from pymol import selector
# Gets the max for use in the residue mapping
##############################################################################################################


class ResInfo:
  def __init__(self, Res, Num, Chain, indices, elems, Lig_Ind, Lig_Elem, Lig_Dist, Res_ind, Res_el):
    self.Res = Res
    self.Num = Num
    self.Chain = Chain
    self.indices = indices
    self.elems = elems
    self.Lig_Ind = Lig_Ind
    self.Lig_Elem = Lig_Elem
    self.Lig_Dist = Lig_Dist
    self.Res_ind = Res_ind
    self.Res_el = Res_el

# Finds the residues within the atoms of the ligand
# Out of those residues, a list is then created [[Res, Num, Chain], corresponding ligand atom]
def AtomsWithin(ligand, protein, cutoff):

    # load in ligand and protein
    cmd.load(ligand)
    cmd.load(protein)

    ligand_name = ligand.rsplit('.', 1)[0]
    protein_name = protein.rsplit('.', 1)[0]

    #Saves the selected atoms within x radius of the ligand
    selected_name = ligand_name + '_' + protein_name + '.pdb'
    cmd.save(selected_name, protein_name + ' within ' + str(cutoff) + ' of ' + ligand_name)
    cmd.load(selected_name)
    selected_name = selected_name.rsplit('.', 1)[0]

 #   print("resi")
    stored.residues = []
    cmd.iterate(selector.process(selected_name), 'stored.residues.append(resi)')
    selection_nums = stored.residues

 #   print("chain")
    stored.residues = []
    # Selected name minus the pdb
    cmd.iterate(selector.process(selected_name), 'stored.residues.append(chain)')
    selected_chains = stored.residues

  #  print("resn")
    stored.residues = []
    # Selected name minus the pdb
    cmd.iterate(selector.process(selected_name), 'stored.residues.append(resn)')
    selected_resns = stored.residues

    x = 0
    Resn_List = []
    Resi_tally = []
    # Makes a pair list of resns and their corresponding numbers
 #   print("before while")
    while x < len(selected_resns):

   #     print("in while")
        stored.residues = []
        # Selected name minus the pdb
        resn = selected_resns[x]
        resi = selection_nums[x]
        chain = selected_chains[x]

        index_sel = protein_name + ' and resn ' + resn + ' and resi ' + resi + ' and chain ' + chain
        cmd.iterate(selector.process(index_sel), 'stored.residues.append(index)')

        Resn_indices = stored.residues
  #      print("resn indices")
   #     print(Resn_indices)

        stored.residues = []
        cmd.iterate(selector.process(index_sel), 'stored.residues.append(elem)')

        Resn_Atoms = stored.residues

    #    print(Resn_Atoms)

        # indices of the selected values
        # Res, Num, Chain, Index, Atom, LigIndex, LigAtom, LigDist)
        item = ResInfo(resn, resi, chain, Resn_indices, Resn_Atoms, 0, 0, 0, 0, 0)

        if item.Num not in Resi_tally:
            Resn_List.append(item)

        x += 1

      #  print("Item Num: " + item.Num)
      #  print("ResiTally: " + x)
        Resi_tally.append(item.Num)

    os.remove(selected_name + '.pdb')
    cmd.delete(selected_name)

    return Resn_List, protein_name, ligand_name

##############################################################################################################

# Then we go through the resn list to find the ligand atom numbers
# Protein and ligand are already loaded

# For a residue, we find the ligand distance for each of the residue atoms
# Pick the smallest distance
# Then out of the possible ligand atoms, we choose the smallest distance
def LigandAtomsforResn(ResList, protein, ligand, cutoff):
    # For each of the residues in the list, go through and find the atoms that correspond with the resns, add to list

    Resn_lig_dist = 1000
    Resn_lig_index = 1000

    for residue in ResList:

        ligand_indices = []
        select_string = ligand + ' within ' + cutoff + ' of ' + protein +' and resn ' + residue.Res + ' and resi ' + str(residue.Num) + ' and chain ' + residue.Chain
        cmd.select('obj', select_string)

        # Ligand indices
        stored.residues = []
        cmd.iterate(selector.process('obj'), 'stored.residues.append(index)')
        ligand_indices = stored.residues

        # Ligand atoms
        stored.residues = []
        cmd.iterate(selector.process('obj'), 'stored.residues.append(elem)')
        ligand_elem = stored.residues

        # ligand atoms are the ones we want to check against.
        dist_counter = 1000
        res_x = 0
        # For each of the Resn atoms, loop through that with each of the ligand atoms
        # inner loop is the ligand atoms cycling through; we keep track of the distance that is the smallest
        for Resn_index in residue.indices:
            lig_x = 0
            for ligand_atom in ligand_indices:

                # For each of these we want a distance measured
                lig_selection = ligand + ' and index ' + str(ligand_atom)
                prot_selection = protein + ' and resn ' + residue.Res + ' and resi ' + str(residue.Num) + ' and chain ' + residue.Chain + ' and index ' + str(Resn_index)

                dist = cmd.get_distance(lig_selection, prot_selection)
        #        print("distance: ")
        #        print(dist)

                if (dist < float(cutoff)) and dist < dist_counter:
                    residue.Lig_Dist = dist
                    residue.Lig_Ind = ligand_atom
                    residue.Lig_Elem = ligand_elem[lig_x]
                    residue.Res_ind = Resn_index
                    residue.Res_el = residue.elems[res_x]
                dist_counter = dist
                lig_x += 1
            res_x += 1

    return ResList

##############################################################################################################
# Writes to a .csv file
# Atoms within "cutoff" of Protein Residues

def FileWrite(ResList, cutoff, LigName):
    cwd = os.getcwd()
  #  os.chdir(FilePath)

    # Sort the list according to the distance
    ResList.sort(key=lambda x: x.Lig_Dist, reverse=False)

    LigName = LigName.rsplit('.', 1)[0]
    FilePath = os.path.join(cwd, 'Fingerprint', LigName + 'ResnInteractions.csv')

    file = open(FilePath, 'w')
    file.write("Ligands within " + str(cutoff) + " of protein residues\n")
    file.write('Ligand Atom #, Ligand Atom, Residue Atom, Resn, Resi, Chain, Residue Atom #, Ligand to Protein Distance\n')

    # For each of the Resn
    for Resn in ResList:
        file.write(str(Resn.Lig_Ind) + ', ' + Resn.Lig_Elem + ', ' + Resn.Res_el + ', ' + Resn.Res + ', '
                   + str(Resn.Num) + ' ,' + Resn.Chain + ' , ' + str(Resn.Res_ind) + ', ' + str(Resn.Lig_Dist))
        file.write('\n')

##############################################################################################################
def Print2D(Ligand, LigandName, Protein, cutoff):
    #cutoff = 2.5
#    print("Atomswithin:")
    Pair_List, protein_name, ligand_name = AtomsWithin(Ligand, Protein, str(cutoff))

 #   print("LigandAtomsforresn")
    ResList = LigandAtomsforResn(Pair_List, protein_name, ligand_name, str(cutoff))

  #  print("Filewrite")
    FileWrite(ResList, str(cutoff), LigandName)

    cmd.delete(ligand_name)
    cmd.delete('obj')
    cmd.delete(protein_name)
