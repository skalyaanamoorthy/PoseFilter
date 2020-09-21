# Take the input files of any extension, and then return them in two folders
# One folder will be PDBLig and one will be PDBComplex

# The file list may consist of pdb or other types of files
# type is either "ligand" or 0 or "complex" or 1"
import os
from pymol import stored
from pymol import selector
from shutil import copyfile
from pymol import cmd
import glob

def DirSearch(keyword, crystalstruct):
    all_files = []
    query = glob.glob("*" + keyword + "*")
    cquery = glob.glob("*" + crystalstruct + "*")

    for x in query:
        if "rot" in x:
            pass
        elif "UNK" in x:
            pass
        elif keyword in x:
            all_files.append(x)
        else:
            pass

    if crystalstruct != "":
        for x in cquery:
            if crystalstruct in x:
                all_files.append(x)
            else:
                pass

    return all_files


# Takes the ligand file name and the protein
def MakeComplex(ligand, proteinpath):
    cmd.load(ligand)
    cmd.load(proteinpath)
    protein_name = os.path.basename(proteinpath)
    ligand_s = ligand.split(".")[0]
    protein_s = protein_name.split(".")[0]
    cmd.save(ligand, ligand_s + ' or ' + protein_s)
    cmd.delete(ligand_s)
    cmd.delete(protein_s)

######################################################################################################################

# Also specify the directory that will be used for the saving
def toPDB(FullMolPath):
    cmd.load(FullMolPath)
    basename = os.path.basename(FullMolPath)
    name, ext = basename.split('.')
    cmd.save(name + '.pdb', name)
    # returns the name of the newly made pdb file
    return name + '.pdb'


######################################################################################################################

def InputFileSort (filelist, type, maindir, pdir, ResId, pname):
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
        UNK_var, pdbfiles = LigandFileSort(filelist, maindir, Ligand_path, Complex_path, pdir)
        return UNK_var, pdbfiles

    else:
        pdbfiles = ComplexFileSort(filelist, maindir, Ligand_path, Complex_path, ResId, pname)
        return ResId, pdbfiles


def LigandFileSort(filelist, maindir, Ligand_path, Complex_path, ppath):
    # Check the first entry of the list for the type.
    # If the type is pdb then just copy over
    print("ligand file sort")

    for file in filelist:
        print("file")
        print(file)
        # Splits on the '.' into the name and ext
        file_name, file_ext = file.split('.')

        # copy file from source to destination
        New_ligpath = os.path.join(Ligand_path, file_name + '.' + file_ext)
        old_ligpath = os.path.join(maindir, file_name + '.' + file_ext)

        copyfile(old_ligpath, New_ligpath)

    New_ppath = os.path.join(Ligand_path, os.path.basename(ppath))
    copyfile(ppath, New_ppath)

    # Check if pdb, if it isn't, make pdb and delete other file
    print("pfile base name")

    pfilewext = os.path.basename(ppath)
    print(pfilewext)

    pfile, pext = pfilewext.split('.')
    if pext != 'pdb':
        os.chdir(Ligand_path)
        toPDB(ppath)
        os.remove(os.path.basename(ppath))

    if filelist[0].split('.')[1] == 'pdb':
        pdbfiles = filelist

    else:
        # Change directories to the ligand files and convert all to pdb.
        os.chdir(Ligand_path)
        pdbfiles = []
        for file in filelist:
            # Splits on the '.' into the name and ext
            file_name, file_ext = file.split('.')

            pdbfile = toPDB(file)
            pdbfiles.append(pdbfile)

        for nonpdb in filelist:
            os.remove(nonpdb)

    os.chdir(Ligand_path)
    stored.residues = []
    cmd.load(pdbfiles[0])

    file_name = pdbfiles[0].split('.')[0]

    cmd.iterate(selector.process(file_name), 'stored.residues.append(resn)')
    UNK_var = stored.residues[0]
    cmd.delete(file_name)

    # Now deal with changing to complexes.
    # Copy over the files to the complex directory
    # switch to the complex directory

    os.chdir(Complex_path)
   # print(Complex_path)
   # print("dealing with complex files")
   # print(pdbfiles)
    for ligand in pdbfiles:

        # new complex path
        new_complex_path = os.path.join(Complex_path, ligand)
        new_ligand_path = os.path.join(Ligand_path, ligand)
      #  print("new complex path: " + new_complex_path)

        # copy file
        copyfile(new_ligand_path, new_complex_path)
     #   print("copy file done. " + ligand)
        # modify to complex
        MakeComplex(ligand, ppath)

    # Change dir back
    os.chdir(maindir)
    return UNK_var, pdbfiles


def ComplexFileSort(filelist, maindir, Ligand_path, Complex_path, ResId, pname):
    # Check the first entry of the list for the type.
    # If the type is pdb then just copy over
    print("complex file sort")

    for file in filelist:
        print("file!")
        print(file)
        # Splits on the '.' into the name and ext
        file_name, file_ext = file.split('.')

        # copy file from source to destination
        New_comppath = os.path.join(Complex_path, file_name + '.' + file_ext)
        old_comppath = os.path.join(maindir, file_name + '.' + file_ext)

        copyfile(old_comppath, New_comppath)


    if filelist[0].split('.')[1] == 'pdb':
        pdbfiles = filelist

    else:
        # Change directories to the complex files and convert all to pdb.
        os.chdir(Complex_path)
        pdbfiles = []
        for file in filelist:
            # Splits on the '.' into the name and ext
            file_name, file_ext = file.split('.')

            pdbfile = toPDB(file)
            pdbfiles.append(pdbfile)

        for nonpdb in filelist:
            os.remove(nonpdb)

    complex_name, ext = pdbfiles[0].split('.')

    # Protein file preparation
    #  print("saved the protein file!!!")
    cmd.load(pdbfiles[0])

    os.chdir(Ligand_path)

    cmd.save(pname, complex_name + ' and not resn ' + ResId)
    cmd.delete(complex_name)

    # Now deal with changing to complexes.

    for complex in pdbfiles:
        # copy the file over
        new_complex_path = os.path.join(Complex_path, complex)
        new_ligand_path = os.path.join(Ligand_path, complex)

        # copy file
        copyfile(new_complex_path, new_ligand_path)

        # name, file extension
        complex_name, file_ext = complex.split('.')
        cmd.load(complex)

        cmd.save(complex_name + '.pdb', complex_name + ' and resn ' + ResId)
        cmd.delete(complex_name)

    # Change dir back
    os.chdir(maindir)
    return pdbfiles


#filelist = ['1fx9_01.mol2', 'pose2.mol2', 'pose3.mol2']
#os.chdir("/home/justine/PycharmProjects/Complex")

#print(InputFileSort(DirSearch("vina", ""), "complex", "/home/justine/PycharmProjects/PoseFilter/Examples/Dimer_Example", "", "UNK", "protein1.pdb"))

#cur_dir = "/home/justine/PycharmProjects/PoseFilter/Examples/Trimer_Example"
#os.chdir(cur_dir)

#InputFileSort((DirSearch("pose", ""), "ligand", cur_dir, "", ""))
