# PoseFilter User Guide 

PoseFilter is a PyMOL plugin and assists in the analysis of docked ligands through identification of unique
oligomeric poses by utilizing RMSD and interaction fingerprint analysis methods.  

## Installing Open-Source PyMOL
### Windows Installation

1. Anaconda 3 can be installed from the website: `https://www.anaconda.com/products/individual`. During GUI installation,
if anaconda is added to the path (through the checkbox), allowing for the `conda` keyword to be used in the command
prompt, step 2 may be skipped.
2. To use conda commands in the Windows command prompt, pathways to Anaconda must be added to the path variable. Search
Windows for "edit the system environment variables," under "system properties" one can choose "Environment Variables"
select "Path" then "Edit..." and add the following paths.
Add the Anaconda3 folder, and Anaconda3\Scripts to the path variable, which should look similar to:
`C:\tools\Anaconda3\Scripts`
`C:\tools\Anaconda3`
3. Close and then open the command prompt again, selecting to "run as administrator." Typing "conda" into the command line
will give conformation that Anaconda3 was installed correctly and the paths were added. A menu with commands should appear.
4. In this command prompt, type the following:
`conda create -n pymol python=3.7 rdkit matplotlib seaborn pip numpy conda-forge::pmw conda-forge::oddt=0.7 tpeulen::pymol-open-source scikit-learn  `
5. Activate the anaconda environment.
`activate pymol`
6. Open source PyMOL can now be opened through the activated anaconda directory. Type pymol into the command line to 
open the open source program.
7. If the command line is closed, `activate pymol` needs to be typed in again before `pymol` to open the program.
After activating the environment, the paths can be added to the environment variables, as Anaconda3 was. Type `path` into the command
prompt, and then copy and paste the links into the path environment variable. Once this is done, `pymol` can be used 
in the command prompt to open the program.


### Linux Installation
1. Ensure that Anaconda3 is installed, if not then install it through the following link:
https://docs.anaconda.com/anaconda/install/linux/
`conda create -n pymol python=3.7 rdkit matplotlib seaborn pip numpy conda-forge::pmw conda-forge::oddt=0.7 tpeulen::pymol-open-source scikit-learn`
2. Activate the anaconda environment:
`source activate pymol`
3. Launch pymol by typing `pymol` in the terminal. If an error occurs, try to install some additional requirements.
Ensure to run as root: https://pymolwiki.org/index.php/Linux_Install.


### MacOS Installation
1. Ensure that Anaconda3 is installed, if not then install it through the following link: https://docs.anaconda.com/anaconda/install/mac-os/
2. Use the following command to create an environment named pymol and to install the proper packages, using anaconda:
`conda create -n pymol python=3.7 rdkit matplotlib seaborn pip numpy conda-forge::pmw conda-forge::oddt=0.7 tpeulen::pymol-open-source scikit-learn `
3. Activate the anaconda environment:
`source activate pymol`
4. Launch pymol by typing `pymol` in the terminal. 


## Installing the PoseFilter Plugin 

   1. Ensure that open-source PyMOL is installed in the conda environment (as outlined above).
   2. The PoseFilter files can be obtained from this repository. Download and extract the files.
   3. In PyMOL, click on "Plugin" at the top bar, then "Plugin Manager." Select the tab "Install New Plugin" and then
     under "Install from local file" select "choose file...". Find the downloaded GitHub folder and select the
      "\_\_init\_\_.py" file from the PoseFilter folder, then "Open" and press "Okay" twice. 

## PoseFilter Use 

### GUI Input (RMS/ Interaction Fingerprint) 

In the PyMOL window, the PoseFilter plugin can be opened from the ‘Plugins’ tab at the top of the page.
There are currently two types of input options are available and are given as two separate tabs in the
PoseFilter GUI panel.
Please use one of these tab options at a time to choose your input files.

#### Tab 1 Description
![Input Type 1](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/GUITab1.jpg)

The first tab (File input - type 1) gives the option to process a
folder that contain protein and ligand poses as separate files.
1. Ensure that the only the protein structure and ligand pose files are available in .pdb format in this directory. 
2. Choose the receptor/protein structure file using the “Choose a protein file” option.
3. Provide a filename identifier that corresponds to the ligand poses in the “ligand filename identifier” box.
Note this should match the filename of the ligand files in the folder and not the protein. For example, if the ligand poses
are available in files named as: pose1.pdb, pose2.pdb, pose3.pdb and the protein file is named as protein.pdb,
the keyword “pose” should be used in the “Ligand filename identifier” box.  
4. Users can provide an output suffix for labelling purposes. This could be the protein PDB ID, for example.
This is optional and the box can be left blank.

#### Tab 2 Description
![Input Type 2](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/GUITab2.jpg)

The second tab (File input - type 2) allows the input of protein + ligand complex files.

1. A directory containing only complex files should be selected by the user. 
2. A keyword may be chosen for the files. In the case of protein + ligand complex files labelled as complex1.pdb,
complex2.pdb and complex3.pdb an appropriate keyword to use is “complex.” 
3. The ligand residue identifier refers to the three-letter code for the ligand residue name as available in the
complex file. Eg., "LIG" or "UNK."
4. Users can provide a suffix name for labelling purposes. This could be the protein PDB ID, for example.
This is optional and the box can be left blank. 

## Parameters 

### RMS Analysis 

1. The RMS similarity cutoff (> 0) is the value used for sorting the poses as ‘Similar’ and ‘Unique’. The default
auto-filled value is 2.0 Å, and this value can be changed by the user.  
2. The “align only α carbon” box can be checked to just align this atom for every rotation. By default, the
protein backbone atoms are used for fitting the receptor for every rotation. This option is recommended when
flexible/multiple receptor docking was used.
3. There is a checkbox called "nonidentical chains." If this box is checked, the program will find which chains are
 identical (number and order of amino acids), and perform the rotations based on this. This program makes the
 assumption that chains are identical because otherwise one or two amino acid substitutions would result in chains
 labelled
 as "nonidentical." The rotation choice is left up to the user. 
3. The “Calculate pose similarity with RMS” button will run the RMSD-based calculations. The outputs, RMS.csv file
and RMS.png (the heatmap file) will be placed in a sub-folder named ‘RMS’. Pose files are initially organized into
 PDBLigand and PDBComplex folders, for organizational purposes. Then pose files are sorted and moved to ‘Unique’ and
 ‘Similar’ folders based on the user defined cut-off values. There is also a ‘Similar.csv’ file that is generated in
 the ‘Similar’ folder, which shows the similar poses’ relationship to each other. 
4. Depending on the number of folders available, the RMS analysis will take a minute or so, but more time with
increasing poses and a larger n-mer. The RMS analysis is finished once “RMS analysis complete” is displayed in the
PyMOL command window. 

### Interaction Fingerprint Analysis 

1. The Fingerprint similarity cutoff (0-1) is the value used for sorting the files as 'Similar' and 'Unique'.
The default auto-filled value is 0.5, but this value can be changed by the user.  
2. Checkboxes are available to select the type of interaction fingerprints analysis, 'Simple Interaction' and 'SPLIF'.
Both or one can be selected from these options. If none are selected, the ‘Simple Interaction’ fingerprint will be generated. 
3. The “Protein-ligand interaction output” checkbox can be selected, to give a text-based interaction summary, which
includes protein and ligand data. 
4. The “calculate pose similarity with interaction fingerprints” button will run the calculations, and generate the
analysis output, which is placed in a sub-folder named ‘Fingerprint’. 
5. 'Fingerprint.csv' and 'Fingerprint.png' files are generated. Pose files are sorted and moved to 'Unique' and
'Similar' folders based on the user defined cut-off value.  As mentioned before, ‘Similar.csv’ file will be
generated in the ‘Similar’ folder, which shows the similar poses’ relationship to each other. 
6. Once the button has been pressed the calculation will take a few minutes to complete. The "Fingerprint analysis
complete” will be displayed in the PyMOL command window once this process has finished. 

## Example Files

Files in the example folder were prepared using AutoDock Vina, Schrodinger Maestro was used for preprocessing for files
  contained in the "Additional Examples" folder. These files can be run with the PoseFilter plugin.
  
  Run with tab option 1: Trimer and tetramer examples; files except 1JF4 in "Additional Examples" folder
  
  Run with tab option 2: Dimer; 1JF4 (from "Additional Example" folder), which can also be run with the "crystal.pdb"
  file, containing the pose from the PDB database crystal structure.
Other 


## Command Line Input 

### RMS 
1. After following the installation instructions, PoseFilter can be used through the command line as well. The following
commands should be typed into the command line after a new PyMOL session is opened.

     `from pmg_tk.startup.PoseFilter import LigandRMSProcess`
     `from pmg_tk.startup.PoseFilter import ComplexRMSProcess`

    
2. The command line can be used for the same result as the GUI first tab option, where a folder is processed that
contains ligand files and a corresponding protein file. 

3. The commands depend on the input type. 

#### Input Type 1 (Ligand)

Type the following into the command line, filling in the appropriate variables.
`LigandRMSProcess(pfile, keyword, label, RMS_Cutoff, alpha, nonidentical)` 
Ex: LigandRMSProcess('/home/.../PoseFilter/Examples/6ewp/6ewp.pdbqt', 'pose', '', 2.0, 0, 0)

pfile: the full pathway for the protein file, which is in the same directory as the ligand pose files.

keyword: keyword matching the ligand files (for pose1.pdb, pose2.pdb, pose3.pdb), keyword would be `pose`.

label: PDB code, or any kind of label to label the RMS/Fingerprint files.

RMS_cutoff: The RMS similarity cutoff (> 0) is the value used for sorting the poses as ‘Similar’ and ‘Unique’. The default
auto-filled value is 2.0 Å, and this value can be changed by the user.  

alpha: boolean value (0 or 1), 1 indicating that only the α carbon will be aligned for the rotations.

nonidentical: boolean value (0 or 1), 1 indicating that the chains are nonidentical, where a check will be made, and
rotations will only occur for exactly identical chains (residue count and type being equal).

#### Input Type 2 (Complex)

Type the following into the command line, filling in the appropriate variables.
`ComplexRMSProcess(folder_dir, keyword, label, ResId, crystal_struct, RMS_Cutoff, alpha, nonidentical)`
Ex: ComplexRMSProcess('/home/.../PoseFilter/Examples/Dimer_Example', 'pose', '', 'UNK', '', 2.0, 0,0)

folder_dir: The directory that contains the protein complexes.

keyword: keyword matching the complex files (for complex1.pdb, complex2.pdb, complex3.pdb), keyword would be `complex`.

label: PDB code, or any kind of label to label the RMS/Fingerprint files.

ResId: The ligand residue identifier refers to the three-letter code for the ligand residue name as available in the
complex file. Eg., "LIG" or "UNK."

crystal_struct: Full pathway of a crystal structure (protein with ligand), which will be included into the analysis.

RMS_cutoff: The RMS similarity cutoff (> 0) is the value used for sorting the poses as ‘Similar’ and ‘Unique’. The default
auto-filled value is 2.0 Å, and this value can be changed by the user.  

alpha: Boolean value (0 or 1), 1 indicating that only the α carbon will be aligned for the rotations.

nonidentical: Boolean value (0 or 1), 1 indicating that the chains are nonidentical, where a check will be made, and
rotations will only occur for exactly identical chains (residue count and type being equal).

### Interaction Fingerprint 
After following the installation instructions, PoseFilter can be used through the command line as well. The following
commands should be typed into the command line after a new PyMOL session is opened.

     `from pmg_tk.startup.PoseFilter import LigandFP`
     `from pmg_tk.startup.PoseFilter import ComplexFP`

#### Input Type 1 (Ligand)
`LigandFP(pfile, keyword, label, FPList, FP_SI, FP_SPLIF, TextInteraction)`
Ex: LigandFP('/home/.../PoseFilter/Examples/6ewp/6ewp.pdbqt', 'pose', '', ['SPLIF'], 0.5, 0.5, 1)

pfile: The full pathway for the protein file, which is in the same directory as the ligand pose files.

keyword: keyword matching the ligand files (for pose1.pdb, pose2.pdb, pose3.pdb), keyword would be `pose`.

label: PDB code, or any kind of label to label the RMS/Fingerprint files.

FPList: A list of the fingerprint types that will be performed, for example: `['SInteraction', 'SPLIF']` to perform both. Or
`[SPLIF]` to perform just SPLIF.

FP_SI: Simple interaction fingerprint cutoff value. The Fingerprint similarity cutoff (0-1) is the value used for
sorting the files as 'Similar' and 'Unique'.

FP_SPLIF: SPLIF fingerprint cutoff value (0-1).
TextInteraction: bool value (0 or 1). If 1, then generates CSV files containing interaction information of various types.


#### Input Type 2 (Complex)
`ComplexFP(folder_dir, keyword, label, ResId, crystal_struct, FPList, FP_SI, FP_SPLIF, TextInteraction)`
Ex: ComplexFP('/home/.../PoseFilter/Examples/Dimer_Example', 'pose', '', 'UNK', '', ['SPLIF'], 0.5, 0.5, 1)

folder_dir: The directory that contains the protein complexes.

keyword: keyword matching the complex files (for complex1.pdb, complex2.pdb, complex3.pdb), keyword would be `complex`.

label: PDB code, or any kind of label to label the RMS/Fingerprint files.

ResId: The ligand residue identifier refers to the three-letter code for the ligand residue name as available in the
complex file. Eg., "LIG" or "UNK."

crystal_struct: Full pathway of a crystal structure (protein with ligand), which will be included into the analysis.

FPList: A list of the fingerprint types that will be performed, for example: `['SInteraction', 'SPLIF']` to perform both. Or
`[SPLIF]` to perform just SPLIF.

FP_SI: Simple interaction fingerprint cutoff value. The Fingerprint similarity cutoff (0-1) is the value used for
sorting the files as 'Similar' and 'Unique'.

FP_SPLIF: SPLIF fingerprint cutoff value (0-1).

TextInteraction: bool value (0 or 1). If 1, then generates CSV files containing interaction information of various types.
 