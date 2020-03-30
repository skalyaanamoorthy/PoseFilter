
from __future__ import absolute_import
from __future__ import print_function

from .Olig import InputFileDir

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
#import Olig_Link

#olig = 4
file_text = ""


# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

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
    uifile = os.path.join(os.path.dirname(__file__), 'olig_gui_fingerprint.ui')
    global form
    form = loadUi(uifile, dialog)


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

            #return filename

    # callback for the "Browse" button
    def browse_filename():
        global file_text
        QFilename = QFileDialog.getOpenFileName(None, "Choose a protein file...")
        file_text = QFilename[0]

    def fingerprint():
        Test()
        pass

    def run():
        global file_text, dialog

        print(file_text)
        print("Running the oligomer script.")
        #InputFileDir(file_text)
        InputFileDir(file_text)
       # complexes = FindLigands(file_text)
        #OligWrapper(complexes, file_text)

    # Implement fingerprint
  #  def fingerprint():
   #     pass
       # global file_text, dialog
     #   FingerprintWrapper()

    #global file_text

   # form.Button_browse.clicked.connect(get_dir)



    form.Button_browse.clicked.connect(browse_filename)
    form.Docking_analysis.clicked.connect(run)
    form.Fingerprint_button.clicked.connect(fingerprint)
    form.close_button.clicked.connect(dialog.close)
    file_text = form.file_select.text()


    return dialog