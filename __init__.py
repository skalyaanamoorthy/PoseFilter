from __future__ import absolute_import
from __future__ import print_function
from .RMS import CreateRMS
from .Fingerprint import Fingerprint_Wrapper
from .General import natural_sort
from .Input_Organization import InputFileSort, DirSearch
from PyQt5.uic import loadUi
import sys
from pymol.Qt import QtWidgets
from pymol.Qt.utils import getSaveFileNameWithExt
from PyQt5.QtWidgets import QDialog, QFileDialog, QTabWidget
from PyQt5.QtCore import *
from pymol import cmd
from .Input_Organization import LigandRMSProcess, ComplexRMSProcess, LigandFP, ComplexFP
import os
import glob

"""

Followed direction from the Pymol Plugins tutorial the GitHub source code:
https://pymolwiki.org/index.php/Plugins_Tutorial
https://github.com/Pymol-Scripts/pymol2-demo-plugin


"""

def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('PoseFilter', run_plugin_gui)

   # cmd.extend(LigandRMSProcess, "LigandRMSProcess")

dialog = None

def run_plugin_gui():
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()


def make_dialog():

    # create a new Window
    dialog = QDialog()
    uifile = os.path.join(os.path.dirname(__file__), 'olig_gui_tabs_edit.ui')
    global form
    form = loadUi(uifile, dialog)

    # Set the initial values of the cutoff boxes
    form.RMSCutoff.setText("2.0")
    form.FP_SICutoff.setText("0.5")
    form.FP_SPLIFCutoff.setText("0.5")

    def FileGeneration():
        index = QTabWidget.currentIndex(form.Mytab)
        if index == 0:
            pfile = form.file_select.text()
            keyword = form.Ligandkeyword.text()
            folder_dir = os.path.dirname(pfile)
            os.chdir(folder_dir)
            all_files = DirSearch(keyword, "") # including the crystal structure
            pname = os.path.basename(pfile).split('.')[0] + '.pdb'
            UNK, pdb_files = InputFileSort(all_files, "ligand", folder_dir, pfile, "", pname)
            return pdb_files, UNK, folder_dir, pname

        else:
            folder_dir = form.dir_select.text()
            complex = form.complexid.text()
            resInput = form.resInput.text()
            crystal_struct = form.crystal_structure.text()
            if crystal_struct != "":
                crystal_keyword = os.path.basename(crystal_struct).split('.')[0]
            else:
                crystal_keyword = ""
            os.chdir(folder_dir)
            all_files = DirSearch(complex, crystal_keyword)
            pname = "protein1.pdb"
            UNK, pdb_files = InputFileSort(all_files, "complex", folder_dir, "", resInput, pname)
            return pdb_files, UNK, folder_dir, pname

    # callback for the "Browse" button
    def browse_filename():
        QFilename = QFileDialog.getOpenFileName(None, "Choose a protein file...")
        # Change the text in the form
        form.file_select.setText(QFilename[0])

    def browse_crystalstructure():
        QFilename = QFileDialog.getOpenFileName(None, "Choose a crystal structure...")
        # Change the text in the crystal structure portion of the form
        form.crystal_structure.setText(QFilename[0])

    # getting a directory
    def get_dir():
        filedir = str(QFileDialog.getExistingDirectory(None, "Select Directory"))

        if filedir:
            form.dir_select.setText(filedir)


    def fingerprint():
        cmd.reinitialize()
        cmd.do('set retain_order,1')
        ErrorGenerated = 0

        AnyChecked = 0
        IText = []

        if form.SPLIF.isChecked():
            IText.append("SPLIF")
            AnyChecked = 1

        if form.Simple_Interaction.isChecked() or AnyChecked == 0:
            IText.append("SInteraction")

        TextInteraction = 0
        if form.Prox_check.isChecked():
            TextInteraction = 1

        FP_SI = form.FP_SICutoff.text()
        FP_SPLIF = form.FP_SPLIFCutoff.text()

        try:
            files, UNK, cur_dir, pname = FileGeneration()
        except:
            print("Error: please check file inputs.")
            ErrorGenerated = 1

        if ErrorGenerated == 0:
            files = natural_sort(files)
            try:
                Fingerprint_Wrapper(files, IText, form.PDBCODE.text(), FP_SI, FP_SPLIF, TextInteraction, cur_dir, pname)
            except:
                print("Error with fingerprint generation. Please check input files, required packages (specifically ODDT and rdkit builds"
                      "specified in installation), then try again.")

        form.dir_select.setText("")
        form.file_select.setText("")
        form.PDBCODE.setText("")
        form.Ligandkeyword.setText("")
        form.file_select.setText("")
        form.crystal_structure.setText("")
        form.complexid.setText("")


    def run():
        global dialog
        cmd.reinitialize()
        cmd.do('set retain_order,1')
        ErrorGenerated = 0
       # cmd.do('set pdb_retain_ids,1')
       # InfoArray = TabInfo()
        alpha = 0
        if form.Alignalpha.isChecked():
            alpha = 1
        nonidentical = 0
        if form.Nonchains.isChecked():
            nonidentical = 1
        try:
            files, UNK, cur_dir, pname = FileGeneration()
        except:
            print("Error: please check input files.")
            ErrorGenerated = 1

        if ErrorGenerated == 0:
            files = natural_sort(files)
            CreateRMS(files, form.PDBCODE.text(), form.RMSCutoff.text(), UNK, alpha, nonidentical, cur_dir, pname)

    form.Button_browse.clicked.connect(browse_filename)
    form.dirButton.clicked.connect(get_dir)
  #  form.dirButton.clicked.connect(get_dir)
    form.crystal_structure_b.clicked.connect(browse_crystalstructure)

    # Oligomer portion is initiated
    form.Docking_analysis.clicked.connect(run)
    form.Fingerprint_button.clicked.connect(fingerprint)
    form.close_button.clicked.connect(dialog.close)

    return dialog