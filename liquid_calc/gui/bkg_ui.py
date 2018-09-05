# -*- coding: utf-8 -*-
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018, Benedict J Heinen"
__email__ = "benedict.heinen@gmail.com"


from PyQt5.QtCore import Qt, pyqtSignal                                       
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QWidget, QFrame, QGridLayout, QVBoxLayout, \
                            QHBoxLayout, QGroupBox, QPushButton, QLineEdit, \
                            QDoubleSpinBox, QLabel, QScrollArea
import numpy as np
# Local relative imports
from . import plot_widgets
from . import utility
from core import data_manip

class BkgUI(QWidget):
    
    plots_changed = pyqtSignal()
    
    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QHBoxLayout(self)
        self.layout.setSpacing(0)
        
        # Make Config Widget
        self.bkg_config_widget = BkgConfigWidget()
        
        # Make vertical line separator
        self.vline = QFrame()
        self.vline.setFrameShape(QFrame.VLine)
        self.vline.setFrameShadow(QFrame.Sunken)
        self.vline.setObjectName("vline")
        
        # Make Plot Widget
        self.bkg_plot_widget = plot_widgets.BkgPlotWidget()
        
        self.layout.addWidget(self.bkg_config_widget)        
        self.layout.addWidget(self.vline)
        #self.layout.addWidget(QWidget())
        self.layout.addWidget(self.bkg_plot_widget)

        self.layout.setStretch(0,1)
        self.layout.setStretch(1,0)
        self.layout.setStretch(2,5)

        self.setLayout(self.layout)

        self.create_signals()

        self.data = {'data_x': np.asarray([]), 'data_y': np.asarray([]),
                     'bkg_x':  np.asarray([]), 'bkg_y':  np.asarray([]),
                     'bkg_y_sc': np.asarray([]),
                     'cor_x':  np.asarray([]), 'cor_y':  np.asarray([])
                     }
        self.bkg_file = None
        self.data_file = None

    def create_signals(self):
        self.bkg_config_widget.data_files_gb.load_data_btn.clicked.connect(self.load_data)
        self.bkg_config_widget.data_files_gb.load_bkg_btn.clicked.connect(self.load_bkg)
        self.bkg_config_widget.bkg_subtract_gb.bkg_sub_btn.clicked.connect(self.sub_bkg)
        self.bkg_config_widget.bkg_subtract_gb.scale_sb.valueChanged.connect(self.plot_data)
        self.bkg_config_widget.bkg_subtract_gb.toggled.connect(self.sub_bkg)
        
    def load_data(self):
        __file_name = utility.get_filename(io='open')
        if __file_name:
            self.data_file = __file_name
        else:
            return
        self.bkg_config_widget.data_files_gb.data_filename_lbl.setText(self.data_file.split('/')[-1])
        try:
            self.data['data_x'], self.data['data_y'] = np.loadtxt(self.data_file, unpack=True)
        except ValueError as e:
            print('Please check header lines in data file')
        self.data['data_x'], self.data['data_y'] = data_manip.rebin_data(self.data['data_x'], self.data['data_y'])
        self.plot_data()

    def load_bkg(self):
        __file_name = utility.get_filename(io='open', caption='Load Background File')
        if __file_name:
            self.bkg_file = __file_name
        else:
            return
        self.bkg_config_widget.data_files_gb.bkg_filename_lbl.setText(self.bkg_file.split('/')[-1])
        try:
            self.data['bkg_x'], self.data['bkg_y'] = np.loadtxt(self.bkg_file, unpack=True)
        except ValueError as e:
            print('Please check header lines in data file')
        self.data['bkg_x'], self.data['bkg_y'] = data_manip.rebin_data(self.data['bkg_x'], self.data['bkg_y'])
        self.plot_data()

    def plot_data(self):
        if self.data['bkg_y'].size:
            # First scale bkg data
            _bkg_scaling = self.bkg_config_widget.bkg_subtract_gb.scale_sb.value()
            self.data['bkg_y_sc'] = self.data['bkg_y'] * _bkg_scaling
        if self.data['cor_x'].size:
            # Only re-subtract if already subtract button clicked
            if self.bkg_config_widget.bkg_subtract_gb.isChecked() and self.bkg_file != None:
                # Subtract background from raw data if option chosen
                self.data['cor_y'] = self.data['data_y'] - self.data['bkg_y_sc']
            else:
                self.data['cor_y'] = self.data['data_y']
            
        self.bkg_plot_widget.update_plots(self.data)
        # emit signal that data has changed to be picked up by tab2
        self.plots_changed.emit()

    def sub_bkg(self):
        # Cor x = data x        
        self.data['cor_x'] = self.data['data_x']
        self.plot_data()


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
        
        self.grid_layout.addWidget(self.load_data_btn, 0, 0)
        self.grid_layout.addWidget(self.data_lbl_frame, 0, 1)
        self.grid_layout.addWidget(self.load_bkg_btn, 1, 0)
        self.grid_layout.addWidget(self.bkg_lbl_frame, 1, 1)

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

        
        self.load_conv_data_btn = QPushButton("2 theta Data")
        self.data_filename_lbl = QLabel("Filename: ")
        self.data_filename_lbl.setAlignment(Qt.AlignCenter)
        self.save_conv_data_btn = QPushButton("Convert")
        
        self.outer_layout = QVBoxLayout()
        self.outer_layout.setContentsMargins(25, 5, 25, 25)
        self.outer_layout.setSpacing(1)

        self.inner_layout = QHBoxLayout()
        self.inner_layout.setContentsMargins(0,0,0,0)
        self.inner_layout.setSpacing(5)
        
        self.inner_layout.addWidget(self.load_conv_data_btn)
        self.inner_layout.addWidget(self.data_filename_lbl)
        self.outer_layout.addLayout(self.inner_layout)
        self.outer_layout.addWidget(self.save_conv_data_btn)
        
        self.setLayout(self.outer_layout)
