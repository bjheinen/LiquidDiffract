# -*- coding: utf-8 -*-
"""
GUI frontend for python implementation of Eggert method.
"""
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018, Benedict J Heinen"
__email__ = "benedict.heinen@gmail.com"

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QVBoxLayout, QLabel, QLineEdit, QCheckBox, \
                            QComboBox, QGroupBox, QGridLayout, QFrame, \
                            QDialogButtonBox, QMessageBox
from PyQt5.QtGui import QDoubleValidator, QIntValidator, QDialog
import numpy as np
#import os

#from . import bkg_ui
#from . import optim_ui
#from . import results_ui
from core.core import __name__, __version__



    
class AboutDialog(QDialog):
    def __init__(self, parent):
        QDialog.__init__(self, parent)
    
######################
#  def accept(self):
#    for k in self.values:
#      sec, opt=k.split('/')
#      widget, l = self.values[k]
#      config.setValue(sec, opt, l(widget))
#    QtGui.QDialog.accept(self)
        
#def launch_edit_dialog():
#        dialog = EditDialog()
#        dialog_result = dialog.exec_()
#
#        if dialog_result == QtWidgets.QDialog.Accepted:
#            model.RestActionsM.update_title(
#                mc_global.active_rest_action_id_it,
#                dialog.rest_action_title_qle.text()
#            )
#            model.RestActionsM.update_rest_action_image_path(
#                mc_global.active_rest_action_id_it,
#                dialog.temporary_image_file_path_str
#            )
#        else:
#            pass
#
#        return dialog_result 