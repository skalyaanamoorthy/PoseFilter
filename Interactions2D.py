import pymol
import os
from pymol import cmd
import oddt
import oddt.interactions as interactions
#import os
import numpy as np
import glob
import sys
import numpy as np
from scipy.spatial import distance

#  file.write(ligand[x]['atomtype'] + ', ' + protein[x]['atomtype'] + ', ' + protein[x]['resname'] + ', ' + str(protein[x]['resnum']) + ', ' + str(dst) + '\n')
class ProtLigInfo:
  def __init__(self, LigAtomType, ProtAtomType, ProtResName, ProtResNum, Dist):
    self.LigAtomType = LigAtomType
    self.ProtAtomType = ProtAtomType
    self.ProtResName = ProtResName
    self.ProtResNum = ProtResNum
    self.Dist = Dist

# Take the two molecules, protein and ligand and use each of those oddt modules

# Protein molecule, and a list of ligands
def InteractionCheck(proteinpath, Listoflig):

    dirname = os.path.dirname(proteinpath)
   # os.chdir(dirname)

    protein = next(oddt.toolkit.readfile('pdb', proteinpath))
    protein.protein = True

    for ligand_object in Listoflig:
        ligandname = ligand_object.ObjName + '.pdb'

        #Ligand_Name = ligandname.split('.')[0]
        ResReport = ligand_object.OrigPoseName.split('.')[0] + "_ResidueReport.csv"
        path = os.path.join(dirname, 'Fingerprint', ResReport)

        file = open(path, 'w')
        file.write("Ligand interactions with protein residues\n")
        file.close()

        # Read in and define the reference ligand
        ligand = next(oddt.toolkit.readfile('pdb', ligandname))

        # Hydrophobic interactions
        p_hydroph, l_hydroph = interactions.hydrophobic_contacts(protein, ligand)
        InteractionsFile(p_hydroph, l_hydroph, path, 'hydrophobic')

        # h bonds
        p_hbonds, l_hbonds, strict = interactions.hbonds(protein, ligand)
        InteractionsFile(p_hbonds, l_hbonds, path, 'hydrogen bond')

        # halogens
        p_halogen, l_halogen, strict = interactions.halogenbonds(protein, ligand)
        InteractionsFile(p_halogen, l_halogen, path, 'halogen bond')

        # pistacking bonds
        pi_interactions = interactions.pi_stacking(protein, ligand)
        InteractionsFile(pi_interactions[0], pi_interactions[2], path, 'pi stacking')

        # salt bridges
        p_salt_bridges, l_salt_bridges = interactions.salt_bridges(protein, ligand)
        InteractionsFile(p_salt_bridges, l_salt_bridges, path, 'salt bridge')

        # pi_cation
        p_pi_cation, l_pi_cation, strict = interactions.pi_cation(protein, ligand)
        InteractionsFile(p_pi_cation, l_pi_cation, path, 'pi cation')

        # acceptor_metal bonds
        p_acceptor_metal, acceptor_metal, strict = interactions.acceptor_metal(protein, ligand)
        InteractionsFile(p_acceptor_metal, acceptor_metal, path, 'acceptor metal')

        # pi_metal bonds
        p_pi_metal, l_pi_metal, strict = interactions.pi_metal(protein, ligand)
        InteractionsFile(p_pi_metal, l_pi_metal, path, 'pi metal')


def InteractionsFile(protein, ligand, FilePath, Interaction_Name):
    # If it is not an empty array
    if len(protein) != 0:
        file = open(FilePath, 'a')
        file.write("Ligand " + Interaction_Name + " interactions with protein residues\n")
        file.write('Ligand Atom Type, Residue Atom Type, Resn, Resi, Ligand to Protein Distance\n')
        ListOfDist = []


        # Put everything into a list; use an object so that it is neater
        for y in range(len(protein)):
            p_coords = protein[y]['coords']
            l_coords = ligand[y]['coords']
            dst = distance.euclidean(p_coords, l_coords)
            item = ProtLigInfo(ligand[y]['atomtype'], protein[y]['atomtype'], protein[y]['resname'], str(protein[y]['resnum']), dst)
            ListOfDist.append(item)

        # Sort the list according to the distance
        ListOfDist.sort(key=lambda z: z.Dist, reverse=False)

        # Make a list of the items that are already included
        # Resn, Resi list
        ResnResi = []
        for residue in ListOfDist:
            ResnResi_item = [residue.ProtResName, residue.ProtResNum]
            for x in ResnResi:
                if ((residue.ProtResName == x[0]) and (str(residue.ProtResNum) == x[1])):
                    pass
                else:
                    file.write(residue.LigAtomType + ', ' + residue.ProtAtomType + ', ' + residue.ProtResName + ', ' + residue.ProtResNum + ', ' + str(residue.Dist) + '\n')
                ResnResi.append(ResnResi_item)
        file.write('\n')
        file.close()


#proteinpath = '/home/justine/protein.pdb'
#ligandpaths = ['ligand1.pdb', 'ligand2.pdb', 'ligand3.pdb']
#InteractionCheck(proteinpath, ligandpaths)

