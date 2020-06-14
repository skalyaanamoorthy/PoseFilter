import os
from pymol import cmd
from pymol import stored
from pymol import selector
# Gets the max for use in the residue mapping
##############################################################################################################


class ResInfo:
  def __init__(self, Res, Num, Chain, Index, LigIndex, LigDist):
    self.Res = Res
    self.Num = Num
    self.Chain = Chain
    self.Index = Index
    self.LigIndex = LigIndex
    self.LigDist = LigDist

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

    stored.residues = []
    cmd.iterate(selector.process(selected_name), 'stored.residues.append(resi)')
    selection_nums = stored.residues

    stored.residues = []
    # Selected name minus the pdb
    cmd.iterate(selector.process(selected_name), 'stored.residues.append(chain)')
    selected_chains = stored.residues

    stored.residues = []
    # Selected name minus the pdb
    cmd.iterate(selector.process(selected_name), 'stored.residues.append(resn)')
    selected_resns = stored.residues

    x = 0
    Resn_List = []
    Resi_tally = []
    # Makes a pair list of resns and their corresponding numbers
    while x < len(selected_resns):

        stored.residues = []
        # Selected name minus the pdb
        resn = selected_resns[x]
        resi = selection_nums[x]
        chain = selected_chains[x]

        index_sel = protein_name + ' and resn ' + resn + ' and resi ' + resi + ' and chain ' + chain
        cmd.iterate(selector.process(index_sel), 'stored.residues.append(index)')

        # indices of the selected values
        item = ResInfo(resn, resi, chain, stored.residues, 0, 0)

        if item.Num not in Resi_tally:
            Resn_List.append(item)

        Resi_tally[x] = item.Num
        x += 1

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
    residue_number = 0
    Resn_lig_dist = 1000
    Resn_lig_index = 1000

    for residue in ResList:

        ligand_atoms = []
        select_string = ligand + ' within ' + cutoff + ' of ' + protein +' and resn ' + residue.Res + ' and resi ' + str(residue.Num) + ' and chain ' + residue.Chain
      #  print(select_string)
        cmd.select('obj', select_string)
        # Ligand atoms
        stored.residues = []

        cmd.iterate(selector.process('obj'), 'stored.residues.append(index)')
        ligand_atoms = stored.residues
      #  print("resi")
      #  print(stored.residues)

        # ligand atoms are the ones we want to check against.

        Lig_dist = 100000
        Lig_index = 100000

        for ligand_atom in ligand_atoms:
        # Now need to cycle through each of the protein resn atoms
            for ind in residue.Index:
                # For each of these we want a distance measured
                 # 'ligand1 and index 5', 'prot and resn glu and resi 5 and chain A and index 77'
                lig_selection = ligand + ' and index ' + str(ligand_atom)
                prot_selection = protein + ' and resn ' + residue.Res + ' and resi ' + str(residue.Num) + ' and chain ' + residue.Chain + ' and index ' + str(ind)
          #      print("prot selection: " + prot_selection)
           #     print("lig selection: " + lig_selection)
                dist = cmd.get_distance(lig_selection, prot_selection)

                if dist < Lig_dist:
                    Lig_dist = dist
                    Lig_index = ligand_atom

            if Lig_dist < Resn_lig_dist:
                residue.LigDist = Lig_dist
                residue.LigIndex = Lig_index

        residue_number += 1

    return ResList

##############################################################################################################
# Writes to a .csv file
# Atoms within "cutoff" of Protein Residues

def FileWrite(ResList, cutoff, LigName):
    cwd = os.getcwd()
  #  os.chdir(FilePath)
    LigName = LigName.rsplit('.', 1)[0]
    FilePath = os.path.join(cwd, 'Fingerprint', LigName + 'ResnInteractions.csv')

    file = open(FilePath, 'w')
    file.write("Ligands within " + str(cutoff) + " of protein residues\n")
    file.write('Resn, Resi, Chain, nearest ligand atom, distance\n')

    # For each of the Resn
    for Resn in ResList:
        file.write(Resn.Res + ', ' + str(Resn.LigIndex) + ', ' + Resn.Chain + ', ' + str(Resn.LigIndex) + ', ' + str(Resn.LigDist))
        file.write('\n')

##############################################################################################################
def Print2D(Ligand, LigandName, Protein, cutoff):
    #cutoff = 2.5
    Pair_List, protein_name, ligand_name = AtomsWithin(Ligand, Protein, str(cutoff))

    ResList = LigandAtomsforResn(Pair_List, protein_name, ligand_name, str(cutoff))

    FileWrite(ResList, str(cutoff), LigandName)

    cmd.delete(ligand_name)
    cmd.delete('obj')
    cmd.delete(protein_name)
    print("Print2D finished processing.")

#Print2D("lig1.pdb", "ligand", "protein.pdb", 2.5)