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
will give confirmation that Anaconda3 was installed correctly and the paths were added. A menu with commands should appear.
4. In this command prompt, type the following:
`conda create -n pymol python=3.7 rdkit=2020.03.3.0 matplotlib seaborn pip numpy conda-forge::pmw tpeulen::pymol-open-source scikit-learn git `
5. Activate the anaconda environment.
`conda activate pymol`
6. Install oddt.
`pip install git+git://github.com/oddt/oddt`

7. Open source PyMOL can now be opened through the activated anaconda directory. Type pymol into the command line to 
open the open source program.
8. If the command line is closed, `conda activate pymol` needs to be typed in again before `pymol` to open the program.
After activating the environment, the paths can be added to the environment variables, as Anaconda3 was. Type `path` into the command
prompt, and then copy and paste the links into the path environment variable. Once this is done, `pymol` can be used 
in the command prompt to open the program.
9. If there are any package errors, oddt can be downgraded to the corresponding build:
`pip install git+git://github.com/oddt/oddt.git@88a5481e0a74348df5f2a9b27132148a3d5b94c1`


### Linux Installation
1. Ensure that Anaconda3 is installed, if not then install it through the following link:
https://docs.anaconda.com/anaconda/install/linux/
2. Open a terminal in the Anaconda3 directory named "envs." Use the following command to create an environment
named pymol that contains the proper packages:
`conda create -n pymol python=3.7 rdkit=2020.03.3 matplotlib seaborn pip numpy conda-forge::pmw tpeulen::pymol-open-source scikit-learn git`
2. Activate the anaconda environment.
`source activate pymol`
3. Install oddt.
`pip install git+git://github.com/oddt/oddt`
4. Launch pymol by typing `pymol` in the terminal. If an error occurs, try to install some additional requirements.
Ensure to run as root: https://pymolwiki.org/index.php/Linux_Install. This is done using an activated environment.
If the terminal is closed, pymol needs to be reactivated in order to be opened again.


### MacOS Installation
Homebrew should be used to install the required packages.

1. Install open-source-pymol using homebrew.

`brew install brewsci/bio/pymol`

2. Need to check version of python that homebrew is using. For the following steps, ensure that the appropriate python path is being referenced.
3. Install additional packages:
`python3 -m pip install matplotlib`
`python3 -m pip install seaborn`
`python3 -m pip install git+git://github.com/oddt/oddt`
4. Install rdkit using homebrew:
`brew tap rdkit/rdkit`
`brew install rdkit`
5. Launch pymol by typing `pymol` in the terminal. This can be done using an activated environment. If the terminal is
closed, pymol needs to be reactivated in order to be opened.
6. If there are any package errors, oddt can be downgraded to the corresponding build:
`pip install git+git://github.com/oddt/oddt.git@88a5481e0a74348df5f2a9b27132148a3d5b94c1`
7. An alternative is to build using conda and the linux instructions.
Troubleshooting: ensure that dependencies are installed though the wiki: https://pymolwiki.org/index.php/MAC_Install.
If no module named or initialized failed errors occur, try installing through: https://github.com/schrodinger/pymol-open-source
and running `python setup.py install`.

## Installing the PoseFilter Plugin 

   1. Ensure that open-source PyMOL is installed in the conda environment (as outlined above).
   2. The PoseFilter files can be obtained from this repository. Download and extract the files. For simplicity, ensure
      that the extracted folder is named `PoseFilter`
   3. In PyMOL, click on "Plugin" at the top bar, then "Plugin Manager." Select the tab "Install New Plugin" and then
     under "Install from local file" select "choose file...". Find the downloaded GitHub folder and select the
      "\_\_init\_\_.py" file from the PoseFilter folder, then "Open" and press "Okay" twice. 

## PoseFilter Use 

### GUI Input (RMS/Interaction Fingerprints) 

In the PyMOL window, the PoseFilter plugin can be opened from the ‘Plugins’ tab at the top of the page.
There are currently two types of input options are available and are given as two separate tabs in the
PoseFilter GUI panel.
Please use one of these tab options at a time to choose your input files.

#### Tab 1 Description
![Input Type 1](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/GUITab1.jpg)

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
![Input Type 2](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/GUITab2.jpg)

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