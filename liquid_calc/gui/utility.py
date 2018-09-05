# -*- coding: utf-8 -*-
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018, Benedict J Heinen"
__email__ = "benedict.heinen@gmail.com"

from PyQt5.QtGui import QDoubleValidator, QIntValidator
from PyQt5.QtWidgets import QFileDialog, QLineEdit, QStyledItemDelegate


def get_filename(io='open', caption='Load Data File', directory=None):
    if io == 'open':
        file_name = QFileDialog.getOpenFileName(caption=caption, directory=directory, filter='All Files (*);;Chi Files (*.chi);;Data Files (*.dat);;xy Files (*.xy)')
        file_name = file_name[0]
    elif io == 'save':
        file_name = QFileDialog.getSaveFileName(caption=caption, directory=directory, filter='All Files (*);;Chi Files (*.chi);;Data Files (*.dat);;xy Files (*.xy)')
        file_name = file_name[0]
        if not file_name.lower().endswith(('.dat', '.chi', '.xy')):
            file_name += '.dat'
    else:
        raise ValueError('Bad argument')
    
    if file_name:
        return file_name
    else:
        return None


class ValidatedItemDelegate(QStyledItemDelegate):

    def __init__(self, parent):
        super(ValidatedItemDelegate, self).__init__(parent)

    def createEditor(self, parent, option, index):
    
        if not index.isValid():
            return 0
        
        if index.column() == 2:
            _editor = QLineEdit(parent)
            _editor.setFrame(False)
            _validator = QIntValidator()
            _editor.setValidator(_validator)
            
        elif index.column() == 3:
            _editor = QLineEdit(parent)
            _editor.setFrame(False)
            _validator = QDoubleValidator()
            _editor.setValidator(_validator)
        else:
            return 0
     
        return _editor


