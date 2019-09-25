# -*- coding: utf-8 -*-
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018, Benedict J Heinen"
__email__ = "benedict.heinen@gmail.com"


from PyQt5.QtCore import Qt, pyqtSignal                                       
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QWidget, QFrame, QGridLayout, QVBoxLayout, \
                            QHBoxLayout, QGroupBox, QPushButton, QLineEdit, \
                            QDoubleSpinBox, QLabel, QScrollArea, QMessageBox, \
                            QCheckBox, QSplitter
import numpy as np
from scipy.optimize import minimize
import os
# Local relative imports
from LiquidDiffract.gui import plot_widgets
from LiquidDiffract.gui import utility
from LiquidDiffract.core import data_manip
from LiquidDiffract.core.core import __name__, __version__

class BkgUI(QWidget):
    
    # plots_changed is handled in main_widget.MainContainer and allows the 
    # next tab to be updated
    plots_changed = pyqtSignal()
    # file_name_changed connects to main_widget.MainContainer.update_filename
    # this method is used to set the working directory based on the data_file
    # that is loaded
    file_name_changed = pyqtSignal()
    
    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QHBoxLayout(self)
        self.layout.setSpacing(0)
        
        # Make Config Widget
        self.bkg_config_widget = BkgConfigWidget()
        self.config_scroll_area = QScrollArea()
        self.config_scroll_area.setFrameShape(QFrame.NoFrame)
        self.config_scroll_area.setWidget(self.bkg_config_widget)
        self.config_scroll_area.setWidgetResizable(True)
        
        # Make vertical line separator
        self.vline = QFrame()
        self.vline.setFrameShape(QFrame.VLine)
        self.vline.setFrameShadow(QFrame.Sunken)
        self.vline.setObjectName("vline")
        
        # Make Plot Widget
        self.bkg_plot_widget = plot_widgets.BkgPlotWidget()
        self.plot_scroll_area = QScrollArea()
        self.plot_scroll_area.setWidget(self.bkg_plot_widget)
        self.plot_scroll_area.setWidgetResizable(True)
        self.plot_scroll_area.setFrameShape(QFrame.NoFrame)
        
        
        self.hsplitter = QSplitter(Qt.Horizontal)
        self.hsplitter.addWidget(self.config_scroll_area)
        self.hsplitter.addWidget(self.plot_scroll_area)
        self.hsplitter.setStretchFactor(0, 2)
        self.hsplitter.setStretchFactor(1, 5)
        self.layout.addWidget(self.hsplitter)
        
        
        #self.layout.addWidget(self.config_scroll_area)        
        #self.layout.addWidget(self.vline)
        #self.layout.addWidget(QWidget())
        #self.layout.addWidget(self.plot_scroll_area)

        #self.layout.setStretch(0,1)
        #self.layout.setStretch(1,0)
        #self.layout.setStretch(2,5)

        self.setLayout(self.layout)

        self.data_file = os.path.abspath(os.getcwd())
        self.file_name_changed.emit()
        
        __a = np.asarray([])
        self.data = {'data_raw_x': __a, 'data_raw_y':   __a,
                     'data_x':     __a, 'data_y':       __a,
                     'bkg_raw_x':  __a, 'bkg_raw_y':    __a,
                     'bkg_x':      __a, 'bkg_y':        __a,
                     'bkg_y_sc':   __a, 'bkg_raw_y_sc': __a,
                     'cor_x':      __a, 'cor_y':        __a
                     }
        self.bkg_file = None
        self.data_file = None
        
        self.create_signals()

    def create_signals(self):
        self.bkg_config_widget.data_files_gb.load_data_btn.clicked.connect(self.load_data)
        self.bkg_config_widget.data_files_gb.load_bkg_btn.clicked.connect(self.load_bkg)
        self.bkg_config_widget.data_files_gb.plot_raw_check.stateChanged.connect(self.plot_data)
        self.bkg_config_widget.bkg_subtract_gb.bkg_sub_btn.clicked.connect(self.sub_bkg)
        self.bkg_config_widget.bkg_subtract_gb.scale_sb.valueChanged.connect(self.plot_data)
        self.bkg_config_widget.bkg_subtract_gb.toggled.connect(self.sub_bkg)
        self.bkg_config_widget.bkg_subtract_gb.auto_sc_btn.clicked.connect(self.auto_scale_bkg)
        
    def load_data(self):
        __file_name = utility.get_filename(io='open')
        if __file_name:
            self.data_file = __file_name
        else:
            return
        try:
            self.data['data_raw_x'], self.data['data_raw_y'] = np.loadtxt(self.data_file, unpack=True)
            self.data['data_x'], self.data['data_y'] = data_manip.rebin_data(self.data['data_raw_x'], self.data['data_raw_y'])
        except ValueError:
            self.data_file = None
            self.load_file_error()
            return
        self.bkg_config_widget.data_files_gb.data_filename_lbl.setText(self.data_file.split('/')[-1])
        self.file_name_changed.emit()
        # Delete any processed data when loading new file            
        self.data['cor_x'] = np.asarray([])
        self.data['cor_y'] = np.asarray([])
        self.plots_changed.emit()    
        if not self.bkg_config_widget.bkg_subtract_gb.isChecked():
            self.sub_bkg()
            self.plots_changed.emit()
        self.plot_data()

    def load_bkg(self):
        __file_name = utility.get_filename(io='open', caption='Load Background File')
        if __file_name:
            self.bkg_file = __file_name
        else:
            return
        try:
            self.data['bkg_raw_x'], self.data['bkg_raw_y'] = np.loadtxt(self.bkg_file, unpack=True)
        except ValueError:
            self.bkg_file = None
            self.load_file_error()
            return
        self.bkg_config_widget.data_files_gb.bkg_filename_lbl.setText(self.bkg_file.split('/')[-1])
        self.data['bkg_x'], self.data['bkg_y'] = data_manip.rebin_data(self.data['bkg_raw_x'], self.data['bkg_raw_y'])
        self.plot_data()

    def plot_data(self):
        if self.data['bkg_y'].size:
            # First scale bkg data
            _bkg_scaling = self.bkg_config_widget.bkg_subtract_gb.scale_sb.value()
            self.data['bkg_raw_y_sc'] = self.data['bkg_raw_y'] * _bkg_scaling
            self.data['bkg_y_sc'] = self.data['bkg_y'] * _bkg_scaling
        if self.data['cor_x'].size:
            # Only re-subtract if already subtract button clicked
            if self.bkg_config_widget.bkg_subtract_gb.isChecked() and self.bkg_file != None:
                # Subtract background from raw data if option chosen
                self.data['cor_y'] = self.data['data_y'] - self.data['bkg_y_sc']
            else:
                self.data['cor_y'] = self.data['data_y']
        if self.bkg_config_widget.data_files_gb.plot_raw_check.isChecked():
            _plot_raw = 1
        else:
            _plot_raw = 0
        self.bkg_plot_widget.update_plots(self.data, _plot_raw)
        # emit signal that data has changed to be picked up by tab2
        self.plots_changed.emit()

    def sub_bkg(self):
        # Cor x = data x
        self.data['cor_x'] = self.data['data_x']
        self.plot_data()
        
    def auto_scale_bkg(self):
        if self.bkg_file == None:
            self.missing_bkg_file_error()
            return
        bkg_scaling = minimize(data_manip.bkg_scaling_residual, 1, 
                               args=(self.data['data_y'], self.data['bkg_y']), 
                               method='nelder-mead', 
                               options={'xtol': 1e-8, 'disp': False})
        bkg_scaling = bkg_scaling.x
        self.bkg_config_widget.bkg_subtract_gb.scale_sb.setValue(bkg_scaling)

    def load_file_error(self):
        message = ['Error loading file!', 'Unable to load file.\nPlease check filename is correct and make sure header lines are commented (#)']
        self.warning_message(message)
        
    def missing_bkg_file_error(self):
        message = ['No background file!', 'Please load a background file']
        self.warning_message(message)
        
    def warning_message(self, message):
        self.error_msg = QMessageBox()
        self.error_msg.setIcon(QMessageBox.Warning)
        self.error_msg.setStandardButtons(QMessageBox.Ok)
        self.error_msg.setText(message[0])
        self.error_msg.setInformativeText((message[1]))
        self.error_msg.setWindowTitle(__name__ + ' v' + __version__)
        self.error_msg.adjustSize()
        self.error_msg.show()


class BkgConfigWidget(QWidget):

    def __init__(self):
        super(BkgConfigWidget, self).__init__()

        self.vlayout = QVBoxLayout()
        self.vlayout.setContentsMargins(0, 0, 5, 0)
        self.vlayout.setSpacing(10)
        self.data_files_gb = DataFilesGroupBox()
        self.bkg_subtract_gb = BkgSubtractGroupBox()
        self.data_conv_gb = DataConvertGroupBox()
        self.vlayout.addWidget(self.data_files_gb, 1)           
        self.vlayout.addWidget(self.bkg_subtract_gb, 1)
        self.vlayout.addWidget(self.data_conv_gb, 1)
        self.vlayout.addWidget(QWidget(), 3)      
        self.setLayout(self.vlayout)
        

class DataFilesGroupBox(QGroupBox):

    def __init__(self, *args):
        super(DataFilesGroupBox, self).__init__(*args)
        self.setTitle('Data Files')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        
        self.grid_layout = QGridLayout()
        self.grid_layout.setContentsMargins(25, 10, 25, 7)
        self.grid_layout.setSpacing(5)

        self.load_data_btn = QPushButton("Load Data")
        self.data_filename_lbl = QLabel("None")
        self.data_filename_lbl.setAlignment(Qt.AlignCenter)
        self.load_bkg_btn = QPushButton("Load Background")
        self.bkg_filename_lbl = QLabel("None")
        self.bkg_filename_lbl.setAlignment(Qt.AlignCenter)
        
        self.data_lbl_frame = QScrollArea()
        self.data_lbl_frame.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.data_lbl_frame.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.data_lbl_frame.setWidgetResizable(True)
        self.data_lbl_frame.setFrameShape(QFrame.NoFrame)

        self.bkg_lbl_frame = QScrollArea()        
        self.bkg_lbl_frame.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.bkg_lbl_frame.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)        
        self.bkg_lbl_frame.setWidgetResizable(True)
        self.bkg_lbl_frame.setFrameShape(QFrame.NoFrame)
        
        
        self.data_lbl_frame.setWidget(self.data_filename_lbl)
        self.bkg_lbl_frame.setWidget(self.bkg_filename_lbl)
        
        
        self.plot_raw_check = QCheckBox()
        self.plot_raw_lbl = QLabel('Plot raw (unbinned) data?')
        self.plot_raw_check.setChecked(False)
        
        self.grid_layout.addWidget(self.load_data_btn, 0, 0)
        self.grid_layout.addWidget(self.data_lbl_frame, 0, 1)
        self.grid_layout.addWidget(self.load_bkg_btn, 1, 0)
        self.grid_layout.addWidget(self.bkg_lbl_frame, 1, 1)
        self.grid_layout.addWidget(self.plot_raw_lbl, 2, 0)
        self.grid_layout.addWidget(self.plot_raw_check, 2, 1)

        self.setLayout(self.grid_layout)


class BkgSubtractGroupBox(QGroupBox):

    def __init__(self, *args):
        super(BkgSubtractGroupBox, self).__init__(*args)
        self.setTitle('Subtract Bkg?')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        self.setCheckable(True)
        
        self.create_widgets()
        self.style_widgets()
        self.create_layout()
        self.create_signals()
        
    def create_widgets(self):
    
        self.scale_lbl = QLabel('Bkg Scaling: ')      
        self.scale_sb = QDoubleSpinBox()        
        self.scale_step = QLineEdit('0.01')      
        self.auto_sc_btn = QPushButton('Auto-scale Bkg')
        self.bkg_sub_btn = QPushButton('Subtract Background')

    def style_widgets(self):
        
        self.scale_lbl.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        
        self.scale_sb.setValue(1.0)
        self.scale_sb.setSingleStep(0.01)
        self.scale_sb.setDecimals(3)
        self.scale_sb.setMinimumWidth(80)
        self.scale_sb.setAlignment(Qt.AlignRight)
        
        self.scale_step.setMaximumWidth(60)
        self.scale_step.setValidator(QDoubleValidator())
        self.scale_step.setAlignment(Qt.AlignRight)
        
    def create_layout(self): 
        self.outer_layout = QVBoxLayout()
        self.outer_layout.setContentsMargins(25, 10, 25, 7)
        
        self.inner_layout = QHBoxLayout()
        self.inner_layout.setContentsMargins(0, 0, 0, 0)   
        self.inner_layout.setSpacing(8)
        self.inner_layout.addWidget(self.scale_lbl, 0)
        self.inner_layout.addWidget(self.scale_sb, 1)
        self.inner_layout.addWidget(self.scale_step, 2)

        self.outer_layout.addLayout(self.inner_layout)
        self.outer_layout.addWidget(self.auto_sc_btn)
        self.outer_layout.addWidget(self.bkg_sub_btn)

        self.setLayout(self.outer_layout)

    def create_signals(self):
        self.scale_step.editingFinished.connect(self.scale_step_changed)

    def scale_step_changed(self):
        self.scale_sb.setSingleStep(float(str(self.scale_step.text())))
        

class DataConvertGroupBox(QGroupBox):

    def __init__(self, *args):
        super(DataConvertGroupBox, self).__init__(*args)
        self.setTitle('Convert 2-theta to Q-space')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        self.setCheckable(True)
        self.setChecked(False)

        self.data_file = None
        self.two_theta_data = None
        
        self.load_conv_data_btn = QPushButton('2 theta Data')
        self.data_filename_lbl = QLabel('None')
        
        self.lambda_lbl = QLabel('Lambda: ')
        #self.lambda_lbl.setAlignment(Qt.AlignRight)
        self.lambda_input = QLineEdit('')
        self.lambda_input.setValidator(QDoubleValidator())
        self.lambda_input.setMaximumWidth(60)
        #self.lambda_input.setAlignment(Qt.AlignLeft)
        
        self.data_filename_lbl.setAlignment(Qt.AlignCenter)
        self.save_conv_data_btn = QPushButton('Convert')
        
        self.outer_layout = QVBoxLayout()
        self.outer_layout.setContentsMargins(25, 25, 25, 25)
        self.outer_layout.setSpacing(10)

        self.inner_layout = QHBoxLayout()
        self.inner_layout.setContentsMargins(0,0,0,0)
        self.inner_layout.setSpacing(10)
        
        self.outer_layout.addWidget(self.load_conv_data_btn)
        self.outer_layout.addWidget(self.data_filename_lbl)
        
        self.inner_layout.addWidget(self.lambda_lbl)
        self.inner_layout.addWidget(self.lambda_input)

        self.outer_layout.addLayout(self.inner_layout)
        self.outer_layout.addWidget(self.save_conv_data_btn)
        
        self.setLayout(self.outer_layout)
        
        self.create_signals()


    def create_signals(self):
        self.load_conv_data_btn.clicked.connect(self.load_two_theta)
        self.save_conv_data_btn.clicked.connect(self.save_q_space)
        
    def load_two_theta(self):

        __file_name = utility.get_filename(io='open', caption='Load 2-theta data file')
        if not __file_name:
            return
        self.data_file = __file_name
        try:
            self.two_theta_data = np.loadtxt(self.data_file, unpack=False)
        except ValueError:
            _message = ['Error loading file!', 'Unable to load file.\nPlease check filename is correct and make sure header lines are commented (#)']
            self.error_msg = QMessageBox()
            self.error_msg.setIcon(QMessageBox.Warning)
            self.error_msg.setStandardButtons(QMessageBox.Ok)
            self.error_msg.setText(_message[0])
            self.error_msg.setInformativeText((_message[1]))
            self.error_msg.setWindowTitle(__name__ + ' v' + __version__)
            self.error_msg.adjustSize()
            self.error_msg.show()
            self.data_file = None
            self.two_theta_data = None
            return
        self.data_filename_lbl.setText(self.data_file.split('/')[-1])
        __split_path = os.path.splitext(self.data_file)
        self.default_fname = __split_path[0] + '_qspace' + __split_path[1]        

        
    def save_q_space(self):
        if self.data_file is None or self.two_theta_data is None:
            return
        try:
            __lambda = np.float(self.lambda_input.text())
        # Return if no wavelength set
        except ValueError:
            _message = ['No Wavelength Set!', 'Please set wavelength value']
            self.error_msg = QMessageBox()
            self.error_msg.setIcon(QMessageBox.Warning)
            self.error_msg.setStandardButtons(QMessageBox.Ok)
            self.error_msg.setText(_message[0])
            self.error_msg.setInformativeText((_message[1]))
            self.error_msg.setWindowTitle(__name__ + ' v' + __version__)
            self.error_msg.adjustSize()
            self.error_msg.show()
            return
        try:
            # Get out filename
            self.convert_filename = utility.get_filename(io='save', caption='Save Q-space Data', directory=self.default_fname)
        except NameError:
            return
        if not self.convert_filename:
            return
        # Convert 2theta data
        __q_data = data_manip.convert_two_theta(self.two_theta_data[0], __lambda)
        __out_data = np.column_stack((__q_data, self.two_theta_data[1]))
        np.savetxt(self.convert_filename, __out_data)
        # Message user success
        self.success_msg = QMessageBox()
        self.success_msg.setIcon(QMessageBox.Information)
        self.success_msg.setStandardButtons(QMessageBox.Ok)
        self.success_msg.setText('Data converted to Q-Space!')
        self.success_msg.setInformativeText(('Lambda: ' + str(__lambda) + '\nConverted data: ' + self.convert_filename))
        self.success_msg.setWindowTitle(__name__ + ' v' + __version__)
        self.success_msg.adjustSize()
        self.success_msg.show()
        # Clear variables
        del __lambda, __q_data, __out_data, self.convert_filename
