#PoseFilter User Guide 

Summary: PoseFilter is a PyMOL plugin and assists in the analysis of docked ligands through identification of unique oligomeric poses by utilizing RMSD and interaction fingerprint analysis methods.  

##Installing the PoseFilter Plugin 

   1. First, open-source PyMOL must be installed. Windows: https://pymolwiki.org/index.php/Windows_Install#Open-Source_PyMOL Linux: https://pymolwiki.org/index.php/Linux_Install  Mac: https://pymolwiki.org/index.php/MAC_Install 

   2. Next, Python modules ‘RDKit’ , ‘Open Drug Discovery Toolkit’, ‘Matplotlib’ (https://matplotlib.org/) and ‘Seaborn’ (https://github.com/mwaskom/seaborn) are needed. These can be installed in the command line through use of anaconda3. Ensure that anaconda3 has the most recent version of these modules. 

    conda install -c rdkit rdkit  

    conda install -c oddt oddt  

    conda install -c anaconda matplotlib 

    conda install -c anaconda seaborn 

   If not, pip (https://pip.pypa.io/en/stable/reference/pip_install/) can be installed and used to retrieve the most recent version of these modules from GitHub. 

   The PoseFilter files can be obtained from this repository. Download and extract the files. Choose the installation file of “_init_.py” from the PoseFilter folder, and press “install”. 

##PoseFilter Use 

###GUI Input (RMS/ Interaction Fingerprint) 

In the PyMOL window, the PoseFilter plugin can be opened from the ‘Plugins’ tab at the top of the page.
There are currently two types of input options are available and are given as two separate tabs in the PoseFilter GUI panel.
Please use one of these tab options at a time to choose your input files.

####Tab 1 Description
The first tab (File input - type 1) gives the option to process a
folder that contain protein and ligand poses are separate files.
1. Ensure that the only the protein structure and ligand pose files are available in .pdb format in this directory. 
2. Choose the receptor/protein structure file using the “Choose a protein file” option.
3. Provide a filename identifier that correspond to the ligand poses in the “ligand filename identifier” box. Note this should match the filename of the ligand files in the folder and not the protein. For example, if the ligand poses are available in files named as: pose1.pdb, pose2.pdb, pose3.pdb and the protein file is named as protein.pdb, the keyword “pose” should be used in the “Ligand filename identifier” box.  
4. Users can provide an output suffix for labelling purposes. This could be the protein PDB ID, for example. This is optional and the box can be left blank. 

Note: This option is useful if you have performed rigid receptor docking.  

####Tab 2 Description 
The second tab (File input - type 2) allows the input of protein + ligand complex files.

1. A directory containing only complex files should be selected by the user. 
2. A keyword may be chosen for the files. In the case of protein + ligand complex files labelled as complex1.pdb, complex2.pdb and complex3.pdb an appropriate keyword to use is “complex.” 
3. The ligand residue identifier refers to the three-letter code for the ligand residue name as available in the complex file. Eg., “LIG” or “UNK” 
4.Users can provide a suffix name for labelling purposes. This could be the protein PDB ID, for example. This is optional and the box can be left blank. 

##Parameters 

###RMS Analysis 

1. The RMS similarity cutoff (> 0) is the value used for sorting the poses as ‘Similar’ and ‘Unique’. The default auto-filled value is 2.0 Å, and this value can be changed by the user.  
2. The “align only α carbon” box can be checked to just align this atom for every rotation. By default, the protein backbone atoms are used for fitting the receptor for every rotation. This option is recommended when flexible/multiple receptor docking was used. 
3. The “Calculate pose similarity with RMS” button will run the RMSD-based calculations. The outputs, RMS.csv file and RMS.png (the heatmap file) will be placed in a sub-folder named ‘RMS’. Pose files are sorted and moved to ‘Unique’ and ‘Similar’ folders based on the user defined cut-off values. There is also a ‘Similar.csv’ file that is generated in the ‘Similar’ folder, which shows the similar poses’ relationship to each other. 
4. Depending on the number of folders available, the RMS analysis will take a minute or so, but more time with increasing poses and a larger n-mer. The RMS analysis is finished once “RMS analysis complete” is displayed in the PyMOL command window. 

###Interaction Fingerprint Analysis 

1. The Fingerprint similarity cutoff (0-1) is the value used for sorting the files as 'Similar' and 'Unique'. The default auto-filled value is 0.5, but this value can be changed by the user.  
2. Checkboxes are available to select the type of interaction fingerprints analysis, 'Simple Interaction' and 'SPLIF'. Both or one can be selected from these options. If none are selected, the ‘Simple Interaction’ fingerprint will be generated. 
3. The “Protein-ligand interaction output” checkbox can be selected, to give a text-based interaction summary, which includes protein and ligand data. 
4. The “calculate pose similarity with interaction fingerprints” button will run the calculations, and generate the analysis output, which is placed in a sub-folder named ‘Fingerprint’. 
5. 'Fingerprint.csv' and 'Fingerprint.png' files are generated. Pose files are sorted and moved to 'Unique' and 'Similar' folders based on the user defined cut-off value.  As mentioned before, ‘Similar.csv’ file will be generated in the ‘Similar’ folder, which shows the similar poses’ relationship to each other. 
6. Once the button has been pressed the calculation will take a few minutes to complete. The "Fingerprint analysis complete” will be displayed in the PyMOL command window once this process has finished. 

##Command Line Input 

###RMS 
After following the installation instructions, PoseFilter can be used through the command line as well. 
1. Type the following into the PyMOL command line:

     `from pmg_tk.startup.PoseFilter import OligWrapper`
    
2. The command line can be used for the same result as the GUI first tab option, where a folder is processed that contains ligand files and a corresponding protein file. 

3. The format is as shown below, where the first field is filled with the protein pathway, the second contains the keyword to recognize the ligand files from the directory. The next field is the output suffix, and then the RMS cutoff value. The last field is a 0 or 1, where 1 represents aligning α helices, and 0 represents aligning with the four-atom option. 

    `OligWrapper(['protein.pdb','ligand keyword'],'Output suffix',2.0,0)` 

    Similarly, the command line can be used to run the RMS analysis in the case of a directory containing protein + ligand complex files. 

    In this case, the fields to be filled out are the protein directory pathway, the complex keyword (to identify the complex files in the directory), then the residue ID (the string that the ligand residues are labelled as), the output suffix, RMS cutoff value. The last field is a 0 or 1, where 1 represents aligning α helices, and 0 represents aligning with the four-atom option. 

    `OligWrapper(['protein.pdbqt', 'complex keyword', 'residue ID'],'Output suffix',2.0,0)`

###Interaction Fingerprint  
1. The following must be typed into the pyMOL command line after the PoseFilter installation has been completed:
`from pmg_tk.startup.PoseFilter import Fingerprint_Wrapper`
2. The following command can be used to run the interaction fingerprint analysis. The arguments are the protein directory, ligand keyword, then either one or both ‘SPLIF’ and ‘SInteraction’ in the array, the output suffix, then the SPLIF and Simple interaction cutoff values. The last value is a Boolean value that contains a 1 (to include the 2D interaction output), and a 0 otherwise.  
`Fingerprint_Wrapper(['protein.pdbqt','ligand keyword'],['SPLIF', 'SInteraction'],'Output suffix', 0.5, 0.5, 1)`
3. Similarly, the command line can be used to run the Interaction Fingerprint analysis in the case of a directory containing protein + ligand complex files. 
4. In this case, the fields to be filled out are the protein directory pathway, the complex keyword (to identify the complex files in the directory), then the residue ID (the string that the ligand residues are labelled as), the output suffix, then the SPLIF and Simple interaction cutoff values, as well as the 0 or 1 to indicate whether the 2D interaction output is included. 
`Fingerprint_Wrapper(['protein.pdb', 'complex keyword', 'residue ID'],['SPLIF', 'SInteraction'],'Output suffix', 0.5, 0.5, 1)`

 