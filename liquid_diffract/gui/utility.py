# -*- coding: utf-8 -*-
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018, Benedict J Heinen"
__email__ = "benedict.heinen@gmail.com"

import os
import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIntValidator, QDoubleValidator, QDialog, QPixmap
from PyQt5.QtWidgets import QFileDialog, QStyledItemDelegate, \
                            QMessageBox, QFrame, QGroupBox, \
                            QVBoxLayout, QGridLayout, QDialogButtonBox, \
                            QLabel, QLineEdit, QCheckBox, QComboBox, QTextBrowser
from core.core import __name__, __version__


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
            _validator = QIntValidator()
            _editor.setValidator(_validator)
        else:
            return 0
     
        return _editor


class PreferencesDialog(QDialog):
    def __init__(self, preferences):
        super(PreferencesDialog, self).__init__()
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setWindowTitle
        self.title = 'Additional Preferences | ' + __name__ + 'v' + __version__
        self.setWindowTitle(self.title)
        self.resize(300, 500)
     
        self.vlayout = QVBoxLayout()
        self.vlayout.setContentsMargins(5, 3, 5, 7)
        self.vlayout.setSpacing(10)

        self.app_settings_gb = AppSettingsGroupBox(preferences)
        self.data_settings_gb = DataSettingsGroupBox(preferences)
        self.refine_settings_gb = RefineSettingsGroupBox(preferences)
        
        self.button_box = QDialogButtonBox()
        self.button_box.addButton('&Cancel', QDialogButtonBox.RejectRole)
        self.button_box.addButton('&Apply', QDialogButtonBox.AcceptRole)

        
        self.vlayout.addWidget(self.app_settings_gb)           
        self.vlayout.addWidget(self.data_settings_gb)
        self.vlayout.addWidget(self.refine_settings_gb)
        self.vlayout.addWidget(self.button_box)
        self.setLayout(self.vlayout)
        
        self.button_box.accepted.connect(self.accept_preferences)
        self.button_box.rejected.connect(self.rejected)
        
        #self.accepted.connect(self.accept_preferences)
        self.rejected.connect(self.close)

    def accept_preferences(self):

        # Log mode
        if self.app_settings_gb.log_mode_input.currentText() == 'Append':
            _append_log_mode = 1
        else:
            _append_log_mode = 0

        try:
            _window_length = np.int(self.data_settings_gb.window_length_input.text())

            _poly_order = np.int(self.data_settings_gb.poly_order_input.text())
            _op_method = self.refine_settings_gb.op_method_input.currentText()
        
            _disp = np.int(self.refine_settings_gb.disp_check.isChecked())
            _maxiter =  np.int(self.refine_settings_gb.maxiter_input.text())
            _maxfun = np.int(self.refine_settings_gb.maxfun_input.text())
            _ftol = np.float(self.refine_settings_gb.ftol_input.text())
            _gtol = np.float(self.refine_settings_gb.gtol_input.text())

        except ValueError:
            _message = ['Missing Values!', 'Please ensures all values are set properly']
            self.error_msg = ErrorMessageBox(_message)
            self.error_msg.show()
            return
        
        if not _window_length%2:
            #replace with a warning messagebox
            _message = ['Error!', 'Please ensure window_length is positive odd integer!']
            self.error_msg = ErrorMessageBox(_message)
            self.error_msg.show()
            return
        
        if _poly_order >= _window_length:
            _message = ['Error!', 'Please ensure poly_order < window_length!']
            self.error_msg = ErrorMessageBox(_message)
            self.error_msg.show()
            return


        # Set minimisation options dictionary
        _minimisation_options = {'disp': _disp,
                                'maxiter': _maxiter,
                                'ftol': _ftol
                                }
        # Add additional optiosn for l_bfgs_b method                        }
        if _op_method == 'L-BFGS-B':
            _minimisation_options['maxfun'] = _maxfun
            _minimisation_options['gtol'] = _gtol
        
        # Set preferences dictionary to return
        self._preferences = {'append_log_mode': _append_log_mode,
                             'window_length': _window_length,
                             'poly_order': _poly_order,
                             'op_method': _op_method,
                             'minimisation_options': _minimisation_options}
        #print(_preferences)

        # handle for ValueError if nothing entered
        self.done(1)
    
    def get_preferences(self):
        return self._preferences


class AppSettingsGroupBox(QGroupBox):
    
    def __init__(self, preferences):
        super(AppSettingsGroupBox, self).__init__()
        self.setTitle('Application Options')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        
        self.create_widgets()
        self.set_data(preferences)
        self.style_widgets()
        self.create_layout()

    def create_widgets(self):
        
        self.log_mode_label = QLabel('Log mode: ')
        self.log_mode_input = QComboBox() 
        self.log_mode_input.insertItem(0, 'Append')
        self.log_mode_input.insertItem(1, 'Overwrite')
        
        
    def set_data(self, preferences):
        self.log_mode_input.setCurrentIndex(1 - preferences['append_log_mode'])

    def style_widgets(self):  
        
        self.log_mode_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.log_mode_input.setMaximumWidth(95)

    def create_layout(self):
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(20, 10, 20, 7)
        self.main_layout.setSpacing(25)


        self.grid_layout = QGridLayout()
        self.grid_layout.setSpacing(15)
        self.grid_layout.addWidget(self.log_mode_label, 0, 0)
        self.grid_layout.addWidget(self.log_mode_input, 0, 1)        
        
        self.main_layout.addLayout(self.grid_layout)
        self.setLayout(self.main_layout)  
     
class DataSettingsGroupBox(QGroupBox):
    
    def __init__(self, preferences):
        super(DataSettingsGroupBox, self).__init__()
        self.setTitle('Data Treatment Options')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        
        self.create_widgets()
        self.set_data(preferences)
        self.style_widgets()
        self.create_layout()

    def create_widgets(self):
        
        self.smoothing_label = QLabel('Savitsky-golay filter parameters (Data smoothing): ')
        self.window_length_label = QLabel('Window size')
        self.window_length_input = QLineEdit()
        self.poly_order_label = QLabel('Poly order')
        self.poly_order_input = QLineEdit()
        
        self.hline = QFrame()
        self.hline.setFrameShape(QFrame.HLine)
        self.hline.setFrameShadow(QFrame.Sunken)
        self.hline.setObjectName("hline")
    
    def set_data(self, preferences):
        self.window_length_input.setText(np.str(preferences['window_length']))
        self.poly_order_input.setText(np.str(preferences['poly_order']))

    def style_widgets(self):  
        
        self.smoothing_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
    
        self.window_length_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.window_length_input.setAlignment(Qt.AlignRight)
        self.window_length_input.setValidator(QIntValidator())
        self.window_length_input.setMaximumWidth(70)
    
        self.poly_order_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.poly_order_input.setAlignment(Qt.AlignRight)
        self.poly_order_input.setValidator(QIntValidator())
        self.poly_order_input.setMaximumWidth(70)    
    

    def create_layout(self):
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(20, 10, 20, 7)
        self.main_layout.setSpacing(25)


        self.grid_layout = QGridLayout()
        self.grid_layout.setSpacing(15)
        
        self.grid_layout.addWidget(self.smoothing_label, 0, 0)
        self.grid_layout.addWidget(self.window_length_label, 1, 0)
        self.grid_layout.addWidget(self.window_length_input, 1, 1)
        self.grid_layout.addWidget(self.poly_order_label, 2, 0)
        self.grid_layout.addWidget(self.poly_order_input, 2, 1)
        
        self.main_layout.addLayout(self.grid_layout)
        self.setLayout(self.main_layout) 



class RefineSettingsGroupBox(QGroupBox):
    
    def __init__(self, preferences):
        super(RefineSettingsGroupBox, self).__init__()
        self.setTitle('Optimisation Options')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        
        self.create_widgets()
        self.set_data(preferences)
        self.style_widgets()
        self.create_layout()
        self.create_signals()


    def create_widgets(self):
        
        self.op_method_label = QLabel('Minimisation Algorithm: ')
        self.op_method_input = QComboBox()
        # Use limited-memory BFGS code for optimising rho
        # See http://users.iems.northwestern.edu/~nocedal/lbfgsb.html for details
        self.op_method_input.insertItem(0, 'L-BFGS-B')
        self.op_method_input.insertItem(1, 'SLSQP')
        #self.op_method_input.insertItem(2, 'COBYLA')
        self.op_method_input.insertItem(2, 'Nelder-Mead')
        #self.op_method_input.insertItem(3, 'Basin Hopping')
        #self.op_method_input.insertItem(4, 'Brute')
        
        
        self.disp_label = QLabel('Output convergence information? ')
        self.disp_check = QCheckBox()
        
        self.maxiter_label = QLabel('Maximum number of iterations: ')
        self.maxiter_input = QLineEdit()     

        # L-BFGS-B specific options
        self.maxfun_label = QLabel('Maximum number of function evaluations: ')
        self.maxfun_input = QLineEdit()
        
        self.ftol_label = QLabel('Function convergence limit (ftol): ')
        self.ftol_input = QLineEdit()   
        
        # gtol not in slsqp
        self.gtol_label = QLabel('Gradient convergence limit (gtol): ')
        self.gtol_input = QLineEdit()
        
        #self.xatol for nelder-mead
        #bfgsb > all of them
        #slsqp > grey out gtol, maxfun
        self.hline = QFrame()
        self.hline.setFrameShape(QFrame.HLine)
        self.hline.setFrameShadow(QFrame.Sunken)
        self.hline.setObjectName("hline")

    def set_data(self, preferences):
        self.op_method_input.setCurrentText(preferences['op_method'])
        _min_options = preferences['minimisation_options']
        self.disp_check.setChecked(_min_options['disp'])
        self.maxiter_input.setText(np.str(_min_options['maxiter']))
        self.ftol_input.setText(np.str(_min_options['ftol']))
        
        # Handle missing maxfun/gtol if op_method!='L-BFGS-B'        
        try:
            self.maxfun_input.setText(np.str(_min_options['maxfun']))
            self.gtol_input.setText(np.str(_min_options['gtol']))
        except KeyError:
            pass

    def style_widgets(self):  
        
        self.op_method_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        
        self.disp_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

             
        self.maxiter_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.maxiter_input.setAlignment(Qt.AlignRight)
        self.maxiter_input.setValidator(QIntValidator())
        self.maxiter_input.setMaximumWidth(70)

        self.maxfun_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.maxfun_input.setAlignment(Qt.AlignRight)
        self.maxfun_input.setValidator(QIntValidator())
        self.maxfun_input.setMaximumWidth(70)

        self.ftol_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.ftol_input.setAlignment(Qt.AlignRight)
        self.ftol_input.setValidator(QDoubleValidator())
        self.ftol_input.setMaximumWidth(70)
        
        self.gtol_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.gtol_input.setAlignment(Qt.AlignRight)
        self.gtol_input.setValidator(QDoubleValidator())
        self.gtol_input.setMaximumWidth(70)
        
        if self.op_method_input.currentText() == 'SLSQP':
            self.maxfun_label.setEnabled(False)
            self.maxfun_input.setEnabled(False)
            self.gtol_label.setEnabled(False)
            self.gtol_input.setEnabled(False)


    def create_layout(self):
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(20, 10, 20, 7)
        self.main_layout.setSpacing(25)


        self.grid_layout = QGridLayout()
        self.grid_layout.setSpacing(15)
        self.grid_layout.addWidget(self.op_method_label, 0, 0)
        self.grid_layout.addWidget(self.op_method_input, 0, 1)        
        self.grid_layout.addWidget(self.hline, 1, 0)
        self.grid_layout.addWidget(self.disp_label, 2, 0)
        self.grid_layout.addWidget(self.disp_check, 2, 1)
        self.grid_layout.addWidget(self.maxiter_label, 3, 0)
        self.grid_layout.addWidget(self.maxiter_input, 3, 1)
        self.grid_layout.addWidget(self.maxfun_label, 4, 0)
        self.grid_layout.addWidget(self.maxfun_input, 4, 1)
        self.grid_layout.addWidget(self.ftol_label, 5, 0)
        self.grid_layout.addWidget(self.ftol_input, 5, 1)
        self.grid_layout.addWidget(self.gtol_label, 6, 0)
        self.grid_layout.addWidget(self.gtol_input, 6, 1)        

        self.main_layout.addLayout(self.grid_layout)

        self.setLayout(self.main_layout)  
     
   
    def create_signals(self):
        self.op_method_input.currentIndexChanged.connect(self.op_method_changed)
        

    def op_method_changed(self):
        if self.op_method_input.currentText() == 'L-BFGS-B':
            self.maxfun_label.setEnabled(True)
            self.maxfun_input.setEnabled(True)
            self.gtol_label.setEnabled(True)
            self.gtol_input.setEnabled(True)
            # Set default values when the method is changed
            self.disp_check.setChecked(False)
            self.maxiter_input.setText('15000')
            self.maxfun_input.setText('15000')
            self.ftol_input.setText('2.22e-8')
            self.gtol_input.setText('1e-10')
            
        elif self.op_method_input.currentText() == 'SLSQP':
            self.maxfun_label.setEnabled(False)
            self.maxfun_input.setEnabled(False)
            self.gtol_label.setEnabled(False)
            self.gtol_input.setEnabled(False)
            # Set default values when the method is changed
            self.disp_check.setChecked(False)
            self.maxiter_input.setText('200')
            self.ftol_input.setText('1e-6')

        elif self.op_method_input.currentText() == 'Nelder-Mead':
            pass
        else:
            pass

class ErrorMessageBox(QMessageBox):
    def __init__(self, _message):
        super(ErrorMessageBox, self).__init__()
        self.setIcon(QMessageBox.Warning)
        self.setStandardButtons(QMessageBox.Ok)
        self.setText(_message[0])
        self.setInformativeText((_message[1]))
        self.setWindowTitle(__name__ + ' v' + __version__)
        self.adjustSize()
        
##### About Dialog from here
        
        
class AboutDialog(QDialog):
    def __init__(self):
        super(AboutDialog, self).__init__()
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setWindowTitle
        self.title = __name__ + 'v' + __version__
        self.setWindowTitle(self.title)
        self.resize(590, 555)
     
        self.vlayout = QVBoxLayout()
        self.vlayout.setContentsMargins(5, 15, 5, 7)
        self.vlayout.setSpacing(10)

        self.logo_box = QLabel()
        self.logo_path = os.path.join(os.path.abspath(os.getcwd()), 'data', 'icons', 'logo.png') 
        self.logo = QPixmap(self.logo_path)
        self.logo_box.setPixmap(self.logo)
        self.logo_box.setAlignment(Qt.AlignCenter)

        self.text_display = QTextBrowser()
        self.text_display.setReadOnly(True)
        self.text_display.setFrameStyle(QFrame.NoFrame)
        self.text_display.setStyleSheet("* { background-color: rgba(0, 0, 0, 0); }")
        self.text_display.setOpenExternalLinks(True)
        #self.text_display.setTextInteractionFlags(Qt.LinksAccessibleByMouse)
        
                                 
                                    
        _text = (
            '<p class="western" align="center"><strong>' + __name__ + '</strong></p>'
            '<p class="western" align="center">v' + __version__ + '</p>'
            '<p class="western" align="center">&nbsp;</p>'
            '<p class="western" align="center"><em>LiquidDiffract is a Python implementation of the iterative procedure of Eggert et al. (2002) to obtain information on macroscopic bulk properties (density) and local atomic arrangement (pair distribution function, g(r)) from XRD data of liquids and amorphous solids.</em></p>'
            '<p class="western" align="center">&nbsp;</p>'
            '<p class="western" align="center"><a class="western" href="https://github.com/bjheinen/LiquidDiffract"><span style="color: #000080;"><span lang="zxx"><u>https://github.com/bjheinen/LiquidDiffract</u></span></span></a></p>'
            '<p class="western" align="center">&nbsp;</p>'
            '<p class="western" align="center"><span style="font-size: small;">Copyright &copy; 2019 &ndash; Benedict J Heinen</span></p>'
            '<p class="western" align="center"><span style="font-size: small;">This program comes with absolutely no warranty or guarantee.<br>'
            'See the </span><a class="western" href="https://www.gnu.org/licenses/gpl.html"><span style="color: #000080;"><span style="font-size: small;"><span lang="zxx"><u>GNU General Public Licence, version 3 or later</u></span></span></span></a><span style="font-size: small;"> for details.</span></p>'
                )
        self.text_display.textCursor().insertHtml(_text)


        self.vlayout.addWidget(self.logo_box)
        self.vlayout.addWidget(self.text_display)           

              
        self.setLayout(self.vlayout)
        
   

#
#
#    @Slot()
#    def about(self):
#        """About Spyder"""
#        versions = get_versions()
#        # Show Mercurial revision for development version
#        revlink = ''
#        if versions['revision']:
#            rev = versions['revision']
#            revlink = " (<a href='https://github.com/spyder-ide/spyder/"\
#                      "commit/%s'>Commit: %s</a>)" % (rev, rev)
#        QMessageBox.about(self,
#            _("About %s") % "Spyder",
#            """<b>Spyder %s</b> %s
#            <br>The Scientific Python Development Environment
#            <br>Copyright &copy; The Spyder Project Contributors
#            <br>Licensed under the terms of the MIT License
#            <p>Created by Pierre Raybaut.
#            <br>Developed and maintained by the
#            <a href="%s/blob/master/AUTHORS">Spyder Project Contributors</a>.
#            <br>Many thanks to all the Spyder beta testers and regular users.
#            <p>For help with Spyder errors and crashes, please read our
#            <a href="%s">Troubleshooting page</a>, and for bug reports and
#            feature requests, visit our <a href="%s">Github website</a>.
#            For project discussion, see our <a href="%s">Google Group</a>.
#            <p>This project is part of a larger effort to promote and
#            facilitate the use of Python for scientific and engineering
#            software development. The popular Python distributions
#            <a href="https://www.anaconda.com/download/">Anaconda</a> and
#            <a href="https://winpython.github.io/">WinPython</a>
#            also contribute to this plan.
#            <p>Python %s %dbits, Qt %s, %s %s on %s
#            <p><small>Most of the icons for the Spyder 2 theme come from the Crystal
#            Project (&copy; 2006-2007 Everaldo Coelho). Other icons for that
#            theme come from <a href="http://p.yusukekamiyamane.com/"> Yusuke
#            Kamiyamane</a> (all rights reserved) and from
#            <a href="http://www.oxygen-icons.org/">
#            The Oxygen icon theme</a></small>.
#            """
#            % (versions['spyder'], revlink, __project_url__, __trouble_url__,
#               __project_url__, __forum_url__, versions['python'],
#               versions['bitness'], versions['qt'], versions['qt_api'],
#               versions['qt_api_ver'], versions['system']))