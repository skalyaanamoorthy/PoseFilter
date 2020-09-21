# Take the input files of any extension, and then return them in two folders
# One folder will be PDBLig and one will be PDBComplex

# The file list may consist of pdb or other types of files
# type is either "ligand" or 0 or "complex" or 1"
import os
from shutil import copyfile
#from pymol import cmd


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



# Since this is for the olig function we want to make complexes with the ligand and protein files
def MakeComplex(ligand, protein):
    cmd.load(ligand)
    cmd.load(protein)
    ligand_s = ligand.rsplit(".", 1)[0]
    protein_s = protein.rsplit(".", 1)[0]
    cmd.save(ligand, ligand_s + ' or ' + protein_s)
    cmd.delete(ligand_s)
    cmd.delete(protein_s)

######################################################################################################################

# Also specify the directory that will be used for the saving
def toPDB(FullMolPath):
    cmd.load(FullMolPath)
    basename = os.path.basename(FullMolPath)
    name, ext = basename.split('.')[0] + '.pdb'
    cmd.save(name + '.pdb', name)
    # returns the name of the newly made pdb file


def ComplextoLigand(ComplexId, ResId):

    # for each file
    files = DirSearch(ComplexId)
    newLigands = []


    for complex in files:
        cmd.load(complex)

        # Ligand name, file extension
        l_name, file_ext = complex.rsplit('.', 1)

        cmd.save(l_name + '_l.pdb', l_name + ' and resn ' + ResId)
        newLigands.append(l_name + '_l.pdb')

    # Protein file
    protein = 'protein_l.pdb'
    # Ligand name, file extension
    ligand_name, file_ext = files[0].rsplit('.', 1)
    cmd.save('protein_l.pdb', ligand_name + ' and not resn ' + ResId)

    # returns the oldligands(complex), new ligands, and protein
    return files, newLigands, protein

######################################################################################################################
# Directory must have one protein file (specified) and multiple ligand files
# Takes the protein .pdb/.pdbqt file and then finds the appropriate ligand files of that extension
# outputs the ligand files in a list format.

#Proteintoliglist
def LigandtoComplex(pfile, keyword):
    working_dir = os.path.dirname(pfile)
    os.chdir(working_dir)
  #  print(working_dir)

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
      #  print("ligand: " + ligand)
        lig_name = ligand.rsplit('.', 1)[0]
      #  print("Make complex: ")
        MakeComplex(lig_name + '_complex', ligand, protein_name)

        # These are the names of the complexes that will be iterated through the program
        Saved_Complexes[SNum] = lig_name + '_complex' + ".pdb"
        SNum += 1
   # print("Saved complexes: ")
    # Saved complexes are protein + ligand .pdb named
    return Saved_Complexes, all_files[0], all_files

def InputFileSort (filelist, type, maindir, pdir):
   # working_dir = os.getcwd()
   print("inputfilesort")

    Ligand_path = os.path.join(maindir, 'PDBLigand')
    Complex_path = os.path.join(maindir, 'PDBComplex')

    if not os.path.exists(Ligand_path):
            os.makedirs(Ligand_path)

    if not os.path.exists(Complex_path):
            os.makedirs(Complex_path)

    print("type")

    if type == "ligand":
        LigandFileSort(filelist, maindir, Ligand_path, Complex_path, pdir)
    else:
        pass
     #   ComplexFileSort(filelist, maindir, Complex_path, pdir)



def LigandFileSort(filelist, maindir, Ligand_path, Complex_path, pdir):
    # Check the first entry of the list for the type.
    # If the type is pdb then just copy over
    print("ligand file sort")

    for file in filelist:
        print("file!")
        print(file)
        # Splits on the '.' into the name and ext
        file_name, file_ext = file.split('.')

        # copy file from source to destination
        New_ligpath = os.path.join(Ligand_path, file_name + '.' + file_ext)
        old_ligpath = os.path.join(maindir, file_name + '.' + file_ext)

        copyfile(old_ligpath, New_ligpath)


    if filelist[0].split('.')[1] == 'pdb':
        pdbfiles = filelist

    else:
        # Change directories to the ligand files and convert all to pdb.
        os.chdir(Ligand_path)
        pdbfiles = []
        for file in filelist:
            # Splits on the '.' into the name and ext
            file_name, file_ext = file.split('.')

            pdbfile = toPDB(file_name)
            pdbfiles.append(pdbfile)

        for nonpdb in filelist:
            os.remove(nonpdb)

    # Now deal with changing to complexes.
    # Copy over the files to the complex directory
    for ligand in pdbfiles:
        # new complex path
        new_complex_path = os.path.join(Complex_path, ligand)

        # copy file
        copyfile(New_ligpath, new_complex_path)

        # modify to complex
        MakeComplex(new_complex_path, pdir)

filelist = ['pose1', 'pose2', 'pose3']
InputFileSort(filelist, "ligand", "/home/justine/PycharmProjects/Trimer_Example", "/home/justine/PycharmProjects/Trimer_Example/protfix.pdb")