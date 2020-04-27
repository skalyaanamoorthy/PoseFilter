from __future__ import absolute_import
from __future__ import print_function
from .Olig import OligWrapper
from .Fingerprint import Fingerprint_Wrapper
from PyQt5.uic import loadUi
import sys
from pymol.Qt import QtWidgets
from pymol.Qt.utils import getSaveFileNameWithExt
from PyQt5.QtWidgets import QDialog, QFileDialog
from PyQt5.QtCore import *
from pymol import cmd
import os
import glob

file_text = ""

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Pose Filter', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    '''
    Open our custom dialog
    '''
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()


def make_dialog():
    # entry point to PyMOL's API

    # create a new Window
    dialog = QDialog()

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'olig_gui_fingerprint_new.ui')
    global form
    form = loadUi(uifile, dialog)

    # Set the initial values of the cutoff boxes
    form.RMSCutoff.setText("2.0")
    form.FP_SICutoff.setText("0.5")
   # form.FP_IntCutoff.setText("0.6")
    form.FP_SPLIFCutoff.setText("0.5")

    # callback for the "Browse" button
    def browse_filename():
        global file_text
        QFilename = QFileDialog.getOpenFileName(None, "Choose a protein file...")
        file_text = QFilename[0]
        # Change the text in the form
        form.file_select.setText(QFilename[0])
        print(file_text)

    # getting a directory
    def get_dir():
        global file_text
        filedir = str(QFileDialog.getExistingDirectory(None, "Select Directory"))
        # self.dlg.download_path.setText(filename)
        #form.file_select.text = filedir
        #filename = getSaveFileNameWithExt(dialog, 'Save As...', filter='PDB File (*.PDB)')
        if filedir:
            form.file_select.setText(filedir)
            file_text = filedir


    def fingerprint():
        global file_text
        cmd.reinitialize()

        print('Run fingerprint portion.')
     #   print("Protein file: " + file_text)
        AnyChecked = 0
        IText = []

        if form.SPLIF.isChecked():
            IText.append("SPLIF")
            AnyChecked = 1
       # if form.Interaction.isChecked():
       #     IText.append("Interaction")
       #     AnyChecked = 1
        if form.Simple_Interaction.isChecked() or AnyChecked == 0:
            IText.append("SInteraction")
        print(IText)

        FP_SI = form.FP_SICutoff.text()
       # FP_I =form.FP_IntCutoff.text()
        FP_SPLIF = form.FP_SPLIFCutoff.text()

        Fingerprint_Wrapper(file_text, IText, form.PDBCODE.text(), FP_SI, FP_SPLIF)

    def run():
        global file_text, dialog
        cmd.reinitialize()
        print(file_text)
        print("Running the oligomer script.")

        OligWrapper(file_text, form.PDBCODE.text(), form.RMSCutoff.text())

    form.Button_browse.clicked.connect(browse_filename)
    form.Docking_analysis.clicked.connect(run)
    form.Fingerprint_button.clicked.connect(fingerprint)
    form.close_button.clicked.connect(dialog.close)
    file_text = form.file_select.text()
    return dialog