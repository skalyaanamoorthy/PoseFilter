# PoseFilter Worked Examples

Please refer to the user guide for explanation of PoseFilter options. Described is usage of these options
with example files (found in the "Examples" folder).

PoseFilter is a PyMOL plugin and assists in the analysis of docked ligands through identification of unique
oligomeric poses by utilizing RMSD and interaction fingerprint analysis methods.

## Description

Files in the example folder were prepared using AutoDock Vina, Schrodinger Maestro was used for preprocessing for files
  contained in the "Additional Examples" folder. These files can be run with the PoseFilter plugin.

Each example type can be run with tab 1 or tab 2 options. The fingerprint assumes that the ligand and protein can be separated
from each other: either pre-separated through ligand and protein input, or by the ligand identifier keyword. In addition, a crystal
structure is included with the Monomer example as well, and can be selected for the tab 2 input. Command line inputs are listed with each example.
Detailed command line example instructions are listed at the end of the document.

To use PoseFilter can be used through the command line as well, type the following commands into the command line after
a new PyMOL session is opened. The name (in this case PoseFilter) should correspond to the folder. The package name can be
checked in "Plugin" -> "Plugin Manager" If it is not "PoseFilter" please use the appropriate keyword  in the following commands.
   
`from pmg_tk.startup.PoseFilter import LigandRMSProcess`

`from pmg_tk.startup.PoseFilter import ComplexRMSProcess`


### Dimer Example
#### Dimer Input
##### Dimer Tab 1 GUI

![Dimer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Dimer_Tab1.jpg)

##### Dimer Tab 1 Command Line

RMS:

`LigandRMSProcess('/home/.../PoseFilter/Dimer_Example/Tab1_input/1FX9.pdbqt', 'pose', 'Dimer', 2.0, 0, 0)`

Fingerprint:

`LigandFP('/home/.../PoseFilter/Dimer_Example/Tab1_input/1FX9.pdbqt', 'pose', 'Dimer', ['SInteraction', 'SPLIF'], 0.5, 0.5, 1)`

##### Dimer Tab 2 GUI

![Dimer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Dimer_Tab2.jpg)

##### Dimer Tab 2 Command Line

RMS:

`ComplexRMSProcess('/home/.../PoseFilter/Dimer_Example/Tab2_input', 'pose', 'Dimer', 'MJI', '', 2.0, 0,0)`

Fingerprint:

`ComplexFP('/home/.../PoseFilter/Dimer_Example/Tab2_input', 'pose', 'Dimer', 'MJI', '', ['SInteraction', 'SPLIF'], 0.5, 0.5, 1)`

#### Dimer Output

##### RMS output
![Dimer RMS Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Dimer_RMS.jpg)

##### SPLIF Fingerprint output
![Dimer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Dimer_SPLIF.jpg)

##### Simple Interaction Fingerprint output
![Dimer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Dimer_SInteraction.jpg)



### Trimer Example
#### Trimer Input
##### Trimer Tab 1 GUI
![Trimer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Trimer_Tab1.jpg)

##### Trimer Tab 1 Command Line

RMS:

`LigandRMSProcess('/home/.../PoseFilter/Trimer_Example/Tab1_input/5EIL.pdb', 'pose', 'Trimer', 2.0, 0, 0)`

Fingerprint:

`LigandFP('/home/.../PoseFilter/Trimer_Example/Tab1_input/5EIL.pdb', 'pose', 'Trimer', ['SInteraction', 'SPLIF'], 0.5, 0.5, 1)`

##### Trimer Tab 2 GUI
![Trimer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Trimer_Tab2.jpg)

##### Trimer Tab 2 Command Line

RMS:

`ComplexRMSProcess('/home/.../PoseFilter/Trimer_Example/Tab2_input', 'pose', 'Trimer', 'UNK', '', 2.0, 0,0)`

Fingerprint:

`ComplexFP('/home/.../PoseFilter/Trimer_Example/Tab2_input', 'pose', 'Trimer', 'UNK', '', ['SInteraction', 'SPLIF'], 0.5, 0.5, 1)`

#### Trimer Output

##### RMS output
![Trimer RMS Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Trimer_RMS.jpg)

##### SPLIF Fingerprint output
![Trimer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Trimer_SPLIF.jpg)

##### Simple Interaction Fingerprint output
![Trimer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Trimer_SInteraction.jpg)



### Tetramer Example
#### Tetramer Input
##### Tetramer Tab 1 GUI
![Tetramer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Tetramer_Tab1.jpg)

##### Tetramer Tab 1 Command Line

RMS:

`LigandRMSProcess('/home/.../PoseFilter/Tetramer_Example/Tab1_input/5kmh.pdbqt', 'pose', 'Tetramer', 2.0, 0, 0)`

Fingerprint:

`LigandFP('/home/.../PoseFilter/Tetramer_Example/Tab1_input/5kmh.pdbqt', 'pose', 'Tetramer', ['SInteraction', 'SPLIF'], 0.5, 0.5, 1)`

##### Tetramer Tab 2 GUI
![Tetramer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Tetramer_Tab2.jpg)

##### Tetramer Tab 2 Command Line

RMS:

`ComplexRMSProcess('/home/.../PoseFilter/Tetramer_Example/Tab2_input', 'pose', 'Tetramer', '6U8', '', 2.0, 0,0)`

Fingerprint:

`ComplexFP('/home/.../PoseFilter/Tetramer_Example/Tab2_input', 'pose', 'Tetramer', '6U8', '', ['SInteraction', 'SPLIF'], 0.5, 0.5, 1)`

#### Tetramer Output

##### RMS output
![Tetramer RMS Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Tetramer_RMS.jpg)

##### SPLIF Fingerprint output
![Tetramer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Tetramer_SPLIF.jpg)

##### Simple Interaction Fingerprint output
![Tetramer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Tetramer_SInteraction.jpg)



### Monomer Example
#### Monomer Input
##### Monomer Tab 1 GUI
![Monomer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Monomer_Tab1.jpg)

##### Monomer Tab 1 Command Line

RMS:

`LigandRMSProcess('/home/.../PoseFilter/Monomer_Example/Tab1_input/protein1.pdb', 'pose', 'Monomer', 2.0, 0, 0)`

Fingerprint:

`LigandFP('/home/.../PoseFilter/Monomer_Example/Tab1_input/protein1.pdb', 'pose', '', ['SInteraction', 'SPLIF'], 0.5, 0.5, 1)`

##### Monomer Tab 2 GUI
![Monomer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Monomer_Tab2.jpg)

##### Monomer Tab 2 Command Line

RMS:

`ComplexRMSProcess('/home/.../PoseFilter/Monomer_Example/Tab2_input', 'pose', 'Monomer', 'SHH', '', 2.0, 0,0)`

Fingerprint:

`ComplexFP('/home/.../PoseFilter/Monomer_Example/Tab2_input', 'pose', 'Monomer', 'SHH', '', ['SInteraction', 'SPLIF'], 0.5, 0.5, 1)`

#### Monomer Output

##### RMS output
![Monomer RMS Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Monomer_RMS.jpg)

##### SPLIF Fingerprint output
![Monomer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Monomer_SPLIF.jpg)

##### Simple Interaction Fingerprint output
![Monomer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Monomer_SInteraction.jpg)



### Heteromer Example
#### Heteromer Input
##### Heteromer Tab 1 GUI
![Heteromer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Heteromer_Tab1.jpg)

##### Heteromer Tab 1 Command Line

RMS:

`LigandRMSProcess('/home/.../PoseFilter/Heteromer_Example/Tab1_input/1TOQ.pdbqt', 'pose', 'Heteromer', 2.0, 0, 0)`

Fingerprint:

`LigandFP('/home/.../PoseFilter/Heteromer_Example/Tab1_input/1TOQ.pdbqt', 'pose', 'Heteromer', ['SInteraction', 'SPLIF'], 0.5, 0.5, 1)`

##### Heteromer Tab 2 GUI
![Heteromer Input](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Heteromer_Tab2.jpg)

##### Heteromer Tab 2 Command Line

RMS:

`ComplexRMSProcess('/home/.../PoseFilter/Dimer_Example/Tab2_input', 'pose', 'Heteromer', 'MJI', '', 2.0, 0,0)`

Fingerprint:

`ComplexFP('/home/.../PoseFilter/Dimer_Example/Tab2_input', 'pose', 'Heteromer', 'UNK', '', ['SInteraction', 'SPLIF'], 0.5, 0.5, 1)`

#### Heteromer Output

##### RMS output
![Heteromer RMS Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Heteromer_RMS.jpg)

##### SPLIF Fingerprint output
![Heteromer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Heteromer_SPLIF.jpg)

##### Simple Interaction Fingerprint output
![Heteromer SPLIF FP Output](https://github.com/skalyaanamoorthy/PoseFilter/blob/master/Snapshots/Heteromer_SInteraction.jpg)



 