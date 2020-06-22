from __future__ import absolute_import
from __future__ import print_function
from .RMS import OligWrapper
from .Fingerprint import Fingerprint_Wrapper
from PyQt5.uic import loadUi
import sys
from pymol.Qt import QtWidgets
from pymol.Qt.utils import getSaveFileNameWithExt
from PyQt5.QtWidgets import QDialog, QFileDialog, QTabWidget
from PyQt5.QtCore import *
from pymol import cmd
import os
import glob

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Pose Filter', run_plugin_gui)

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

    # create a new Window
    dialog = QDialog()
    uifile = os.path.join(os.path.dirname(__file__), 'olig_gui_tabs_edit.ui')
    global form
    form = loadUi(uifile, dialog)

    # Set the initial values of the cutoff boxes
    form.RMSCutoff.setText("2.0")
    form.FP_SICutoff.setText("0.5")
   # form.FP_IntCutoff.setText("0.6")
    form.FP_SPLIFCutoff.setText("0.5")
    form.Prox_val.setText("2.5")

    def TabInfo():
        index = QTabWidget.currentIndex(form.Mytab)
        if index == 0:
            print("Tab1")
            pfile = form.file_select.text()
            print(pfile)
            keyword = form.Ligandkeyword.text()
            print(keyword)
            InfoArray = [pfile, keyword]
        else:
            print("Tab2")
            cur_dir = form.dir_select.text()
            complex = form.complexid.text()
            resInput = form.resInput.text()
            print(cur_dir)
            print(complex)
            print(resInput)
            InfoArray = [cur_dir, complex, resInput]

        return InfoArray

    # callback for the "Browse" button
    def browse_filename():
        QFilename = QFileDialog.getOpenFileName(None, "Choose a protein file...")
        # Change the text in the form
        form.file_select.setText(QFilename[0])

    # getting a directory
    def get_dir():
        filedir = str(QFileDialog.getExistingDirectory(None, "Select Directory"))

        if filedir:
            form.dir_select.setText(filedir)
            print(filedir)


    def fingerprint():
        cmd.reinitialize()
     #   print('Run fingerprint portion.')

        InfoArray = TabInfo()

        AnyChecked = 0
        IText = []

        if form.SPLIF.isChecked():
            IText.append("SPLIF")
            AnyChecked = 1

        if form.Simple_Interaction.isChecked() or AnyChecked == 0:
            IText.append("SInteraction")
        print(IText)

        TextInteraction = [0, form.Prox_val.text()]
        if form.Prox_check.isChecked():
            TextInteraction[0] = 1

        FP_SI = form.FP_SICutoff.text()
        FP_SPLIF = form.FP_SPLIFCutoff.text()

        Fingerprint_Wrapper(InfoArray, IText, form.PDBCODE.text(), FP_SI, FP_SPLIF, TextInteraction)

    def run():
        global dialog
        cmd.reinitialize()
        InfoArray = TabInfo()
        alpha = 0
        if form.Alignalpha.isChecked():
            alpha = 1

        print("Running the oligomer script.")
        OligWrapper(InfoArray, form.PDBCODE.text(), form.RMSCutoff.text(), alpha)

    form.Button_browse.clicked.connect(browse_filename)
    form.dirButton.clicked.connect(get_dir)

    # Oligomer portion is initiated
    form.Docking_analysis.clicked.connect(run)
    form.Fingerprint_button.clicked.connect(fingerprint)
    form.close_button.clicked.connect(dialog.close)

    return dialog