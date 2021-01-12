# PoseFilter Worked Examples

Please refer to the user guide for explanation of PoseFilter options. Described is usage of these options
with example files (found in the Examples folder).

PoseFilter is a PyMOL plugin and assists in the analysis of docked ligands through identification of unique
oligomeric poses by utilizing RMSD and interaction fingerprint analysis methods.

## Description

Files in the example folder were prepared using AutoDock Vina, Schrodinger Maestro was used for preprocessing for files
  contained in the "Additional Examples" folder. These files can be run with the PoseFilter plugin.
  
  Run with tab option 1: In the Examples folder; Dimer_Example, Trimer_Example, Tetramer_Example and Heteromer_Example folders
  
  Run with tab option 2: Monomer_Example (from Examples folder), which can also be run with the "crystal.pdb"
  file, containing the pose from the PDB database crystal structure.

### Dimer Example
#### Dimer Input
##### Dimer Tab 1

![Dimer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Dimer_Tab1.jpg)

##### Dimer Tab 2

![Dimer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Dimer_Tab2.jpg)

#### Dimer Output

##### RMS output
![Dimer RMS Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Dimer_RMS.jpg)

##### SPLIF Fingerprint output
![Dimer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Dimer_SPLIF.jpg)

##### Simple Interaction Fingerprint output
![Dimer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Dimer_SInteraction.jpg)



### Trimer Example
#### Trimer Input
##### Trimer Tab 1
![Trimer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Trimer_Tab1.jpg)

##### Trimer Tab 2
![Trimer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Trimer_Tab2.jpg)

#### Trimer Output

##### RMS output
![Trimer RMS Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Trimer_RMS.jpg)

##### SPLIF Fingerprint output
![Trimer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Trimer_SPLIF.jpg)

##### Simple Interaction Fingerprint output
![Trimer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Trimer_SInteraction.jpg)



### Tetramer Example
#### Tetramer Input
##### Tetramer Tab 1
![Tetramer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Tetramer_Tab1.jpg)

##### Tetramer Tab 2
![Tetramer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Tetramer_Tab2.jpg)

#### Tetramer Output

##### RMS output
![Tetramer RMS Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Tetramer_RMS.jpg)

##### SPLIF Fingerprint output
![Tetramer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Tetramer_SPLIF.jpg)

##### Simple Interaction Fingerprint output
![Tetramer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Tetramer_SInteraction.jpg)



### Monomer Example
#### Monomer Input
##### Monomer Tab 1
![Monomer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Monomer_Tab1.jpg)

##### Monomer Tab 2
![Monomer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Monomer_Tab2.jpg)

#### Monomer Output

##### RMS output
![Monomer RMS Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Monomer_RMS.jpg)

##### SPLIF Fingerprint output
![Monomer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Monomer_SPLIF.jpg)

##### Simple Interaction Fingerprint output
![Monomer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Monomer_SInteraction.jpg)



### Heteromer Example
#### Heteromer Input
##### Heteromer Tab 1
![Heteromer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Heteromer_Tab1.jpg)

##### Heteromer Tab 2
![Heteromer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Heteromer_Tab2.jpg)

#### Heteromer Output

##### RMS output
![Heteromer RMS Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Heteromer_RMS.jpg)

##### SPLIF Fingerprint output
![Heteromer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Heteromer_SPLIF.jpg)

##### Simple Interaction Fingerprint output
![Heteromer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Heteromer_SInteraction.jpg)



## Command Line Input 

### RMS 
1. After following the installation instructions, PoseFilter can be used through the command line as well. The following
commands should be typed into the command line after a new PyMOL session is opened. The name (in this case PoseFilter)
should correspond to the folder. If it is 'PoseFilter-master,' then replace that keyword instead of PoseFilter in the following
commands.
   
    `from pmg_tk.startup.PoseFilter import LigandRMSProcess`
   
    `from pmg_tk.startup.PoseFilter import ComplexRMSProcess`

    
2. The command line can be used for the same result as the GUI first tab option, where a folder is processed that
contains ligand files and a corresponding protein file. 

3. The commands depend on the input type. 

#### Input Type 1 (Ligand)

pfile: the full pathway for the protein file, which is in the same directory as the ligand pose files.

keyword: keyword matching the ligand files (for pose1.pdb, pose2.pdb, pose3.pdb), keyword would be `pose`.

label: PDB code, or any kind of label to label the RMS/Fingerprint files.

RMS_cutoff: The RMS similarity cutoff (> 0) is the value used for sorting the poses as ‘Similar’ and ‘Unique’. The default
auto-filled value is 2.0 Å, and this value can be changed by the user.  

alpha: boolean value (0 or 1), 1 indicating that only the α carbon will be aligned for the rotations.

nonidentical: boolean value (0 or 1), 1 indicating that the chains are nonidentical, where a check will be made, and
rotations will only occur for exactly identical chains (residue count and type being equal).

##### Input Example
Type the following into the command line, filling in the appropriate variables. Shown below is a sammple.
`LigandRMSProcess(pfile, keyword, label, RMS_Cutoff, alpha, nonidentical)`

LigandRMSProcess('/home/.../PoseFilter/Dimer_Example/Tab1_input/1FX9.pdbqt', 'pose', '', 2.0, 0, 0)

#### Input Type 2 (Complex)

Type the following into the command line, filling in the appropriate variables.
`ComplexRMSProcess(folder_dir, keyword, label, ResId, crystal_struct, RMS_Cutoff, alpha, nonidentical)`
Ex: ComplexRMSProcess('/home/.../PoseFilter/Dimer_Example/Tab2_input', 'pose', '', 'MJI', '', 2.0, 0,0)

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

     from pmg_tk.startup.PoseFilter import LigandFP
     from pmg_tk.startup.PoseFilter import ComplexFP

#### Input Type 1 (Ligand)
`LigandFP(pfile, keyword, label, FPList, FP_SI, FP_SPLIF, TextInteraction)`
Ex: LigandFP('/home/.../PoseFilter/Dimer_Example/Tab1_input/1FX9.pdbqt', 'pose', '', ['SPLIF'], 0.5, 0.5, 1)

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
Ex: ComplexFP('/home/.../PoseFilter/Dimer_Example/Tab2_input', 'pose', '', 'UNK', '', ['SPLIF'], 0.5, 0.5, 1)

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
 