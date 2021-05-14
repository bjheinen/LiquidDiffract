# -*- coding: utf-8 -*-
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2021, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

# What imports do I need??

import os
import numpy as np
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QDoubleValidator, QFont
from PyQt5.QtWidgets import QWidget, QFrame, QGridLayout, QVBoxLayout, \
                            QHBoxLayout, QGroupBox, QPushButton, QRadioButton, \
                            QLineEdit, QDoubleSpinBox, QLabel, QScrollArea, \
                            QMessageBox, QCheckBox, QSplitter, QToolTip
# LiquidDiffract imports
from LiquidDiffract.gui import plot_widgets
from LiquidDiffract.gui import utility
from LiquidDiffract.core import data_utils
import LiquidDiffract.core.core as core
from LiquidDiffract.version import __appname__, __version__


class StructureUI(QWidget):

    # plots_changed is handled in main_widget.MainContainer and allows the
    # next tab to be updated
    #copy###plots_changed = pyqtSignal()
    # file_name_changed connects to main_widget.MainContainer.update_filename
    # this method is used to set the working directory based on the data_file
    # that is loaded
    #copy###file_name_changed = pyqtSignal()

    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QHBoxLayout(self)
        self.layout.setSpacing(0)
        
        # Make config widget
        self.structure_config_widget = StructureConfigWidget()
        self.config_scroll_area = QScrollArea()
        self.config_scroll_area.setFrameShape(QFrame.NoFrame)
        self.config_scroll_area.setWidget(self.structure_config_widget)
        self.config_scroll_area.setWidgetResizable(True)
        
        
        # Make vertical line separator
        self.vline = QFrame()
        self.vline.setFrameShape(QFrame.VLine)
        self.vline.setFrameShadow(QFrame.Sunken)
        self.vline.setObjectName('vline')
        
        
        # Make plot widget
        self.structure_plot_widget = plot_widgets.StructurePlotWidget()
        self.plot_scroll_area = QScrollArea()
        self.plot_scroll_area.setWidget(self.structure_plot_widget)
        self.plot_scroll_area.setWidgetResizable(True)
        self.plot_scroll_area.setFrameShape(QFrame.NoFrame)
        
        # Create resizeable splitter
        self.hsplitter = QSplitter(Qt.Horizontal)
        self.hsplitter.addWidget(self.config_scroll_area)
        self.hsplitter.addWidget(self.plot_scroll_area)
        self.hsplitter.setStretchFactor(0, 2)
        self.hsplitter.setStretchFactor(1, 5)

        self.layout.addWidget(self.hsplitter)
        
        self.setLayout(self.layout)
        
        # Initialise some data?
        

class StructureConfigWidget(QWidget):
    
    def __init__(self):
        super(StructureConfigWidget, self).__init__()
        
        self.vlayout = QVBoxLayout()
        self.vlayout.setContentsMargins(0, 0, 5, 0)
        self.vlayout.setSpacing(10)
        
        self.plot_view_gb = PlotViewGroupBox()
        self.monatomic_gb = MonatomicGroupBox()
        self.polyatomic_gb = PolyatomicGroupBox()

        self.vlayout.addWidget(self.plot_view_gb, 1)        
        self.vlayout.addWidget(self.monatomic_gb, 1)
        self.vlayout.addWidget(self.polyatomic_gb, 1)

        # Dummy widget to scale vertically
        self.vlayout.addWidget(QWidget(), 3)
        
        self.setLayout(self.vlayout)
        
        # creat signals so gb toggles?? self.create_signals()

class PlotViewGroupBox(QGroupBox):
    
    def __init__(self, *args):
        super(PlotViewGroupBox, self).__init__(*args)
        self.setTitle('Plot View')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')

        self.create_widgets()
        self.style_widgets()
        self.create_layout()
        self.create_signals()
        
    def create_widgets(self):
        
        self.plot_view_label = QLabel('Select function to view: ')
        self.rdf_btn = QRadioButton('RDF(r)')
        self.tr_btn = QRadioButton('T(r)')
        
    def style_widgets(self):
        
        self.rdf_btn.setToolTip('RDF(r) = 4πρr<sup>2</sup>g(r)')
        self.tr_btn.setToolTip('T(r) = RDF(r) / r')
        
        
        self.rdf_btn.setChecked(True)
        
        
        
    def create_layout(self):
        
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(20, 1, 20, 5)
        self.main_layout.setSpacing(1)
        
        self.hbtn_layout = QHBoxLayout()
        self.hbtn_layout.setContentsMargins(0, 0, 0, 0)
        self.hbtn_layout.addWidget(self.rdf_btn)
        self.hbtn_layout.addWidget(self.tr_btn)
        
        self.main_layout.addWidget(self.plot_view_label)
        self.main_layout.addLayout(self.hbtn_layout)
        self.setLayout(self.main_layout)
        
    def create_signals(self):
        pass


class MonatomicGroupBox(QGroupBox):
    
    def __init__(self, *args):
        super(MonatomicGroupBox, self).__init__(*args)
        self.setTitle('Integrate N-1 (Monatomic)')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        self.setCheckable(True)
                           
        self.create_widgets()
        self.style_widgets()               
        self.create_layout()
        self.create_signals()
        
    def create_widgets(self):
        
        self.r0_label = QLabel('r<sub>0</sub>')
        self.r0_input = QDoubleSpinBox()
        
        self.rpmax_label = QLabel('r\'<sub>max</sub>')
        self.rpmax_input = QDoubleSpinBox()
        
        self.rmax_label = QLabel('r<sub>max</sub>')
        self.rmax_input = QDoubleSpinBox()
        
        self.rmin_label = QLabel('r<sub>min</sub>')
        self.rmin_input = QDoubleSpinBox()
        
        self.find_limits_btn = QPushButton('Auto-Refine Int. Limits')
        
        self.calc_N_btn = QPushButton('Calc Coordination Number')
                
        self.Na_output_label = QLabel('N<sub>a</sub> : ')
        self.Na_output = QLineEdit()
        
        self.Nb_output_label = QLabel('N<sub>b</sub> : ')
        self.Nb_output = QLineEdit()
        
        self.Nc_output_label = QLabel('N<sub>c</sub> : ')
        self.Nc_output = QLineEdit()
    
    def style_widgets(self):
        
        self.r0_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.rpmax_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.rmax_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.rmin_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        
        self.Na_output_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.Nb_output_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.Nc_output_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        
        # May need to change max widths
        self.r0_input.setMaximumWidth(82)
        self.rpmax_input.setMaximumWidth(82)
        self.rmax_input.setMaximumWidth(82)
        self.rmin_input.setMaximumWidth(82)
        
        self.Na_output.setMaximumWidth(100)
        self.Nb_output.setMaximumWidth(100)
        self.Nc_output.setMaximumWidth(100)
        
        # Set results as read-only
        self.Na_output.setReadOnly(True)
        self.Nb_output.setReadOnly(True)
        self.Nc_output.setReadOnly(True)
        
        # Set input spinbox settings
        # Might need to set minimum width also to stop scale miss
        # Does dblsb need a validator? - check and the ndelete this comment
        self.r0_input.setSingleStep(0.01)
        self.r0_input.setDecimals(3)
        self.r0_input.setAlignment(Qt.AlignRight)
        self.rpmax_input.setSingleStep(0.01)
        self.rpmax_input.setDecimals(3)
        self.rpmax_input.setAlignment(Qt.AlignRight)
        self.rmax_input.setSingleStep(0.01)
        self.rmax_input.setDecimals(3)
        self.rmax_input.setAlignment(Qt.AlignRight)
        self.rmin_input.setSingleStep(0.01)
        self.rmin_input.setDecimals(3)
        self.rmin_input.setAlignment(Qt.AlignRight)        
        
        # Set tooltips
        self.r0_label.setToolTip('Leading edge of first peak in RDF(r)')
        self.r0_input.setToolTip('Leading edge of first peak in RDF(r)')
        
        self.rmin_input.setToolTip('First minimum after first peak in RDF(r)')
        self.rmin_label.setToolTip('First minimum after first peak in RDF(r)')
        
        self.rpmax_label.setToolTip('Centre of first peak in T(r)')
        self.rpmax_input.setToolTip('Centre of first peak in T(r)')
        
        self.rmax_label.setToolTip('Centre of first peak in RDF(r)')
        self.rmax_input.setToolTip('Centre of first peak in RDF(r)')
        
        self.Na_output_label.setToolTip('N<sub>a</sub> = 2 ∫ 4πρr[rg(r)]<sub>sym</sub> dr | r = r<sub>0</sub> to r = r\'<sub>max</sub>')
        self.Na_output.setToolTip('N<sub>a</sub> = 2 ∫ 4πρr[rg(r)]<sub>sym</sub> dr | r = r<sub>0</sub> to r = r\'<sub>max</sub>')
        
        self.Nb_output_label.setToolTip('N<sub>b</sub> = 2 ∫ 4πρ[r<sup>2</sup>g(r)]<sub>sym</sub> dr | r = r<sub>0</sub> to r = r<sub>max</sub>')
        self.Nb_output.setToolTip('N<sub>b</sub> = 2 ∫ 4πρ[r<sup>2</sup>g(r)]<sub>sym</sub> dr | r = r<sub>0</sub> to r = r<sub>max</sub>')
        
        self.Nc_output_label.setToolTip('N<sub>c</sub> = ∫ 4πρr<sup>2</sup>g(r) dr | r = r<sub>0</sub> to r = r<sub>min</sub>')
        self.test = self.Nc_output.setToolTip('N<sub>c</sub> = ∫ 4πρr<sup>2</sup>g(r) dr | r = r<sub>0</sub> to r = r<sub>min</sub>')
        
        self.setStyleSheet(" QToolTip{ font: 11pt;}")
        # Set ToolTip style sheets separately
        self.Na_output_label.setStyleSheet(" QToolTip{ min-width: 480px; font: 18pt;}")
        self.Na_output.setStyleSheet(" QToolTip{ min-width: 480px; font: 18pt;}")
        self.Nb_output_label.setStyleSheet(" QToolTip{ min-width: 475px; font: 18pt;}")
        self.Nb_output.setStyleSheet(" QToolTip{ min-width: 475px; font: 18pt;}")
        self.Nc_output_label.setStyleSheet(" QToolTip{ min-width: 400px; font: 18pt;}")
        self.Nc_output.setStyleSheet(" QToolTip{ min-width: 400px; font: 18pt;}")  

    
    def create_layout(self):
        self.outer_layout = QVBoxLayout()
        self.outer_layout.setContentsMargins(25, 10, 25, 7)
        self.outer_layout.setSpacing(30)

        self.int_limits_layout = QVBoxLayout()
        self.int_limits_layout.setContentsMargins(0, 0, 0, 0)
        self.int_limits_layout.setSpacing(15)
        
        self.int_results_layout = QVBoxLayout()
        self.int_results_layout.setContentsMargins(0, 0, 0, 0)
        self.int_results_layout.setSpacing(15)        


        # Grid for integration limits
        self.int_limits_grid = QGridLayout()
        self.int_limits_grid.setContentsMargins(0, 0, 0, 0)
        self.int_limits_grid.setSpacing(15)
        self.int_limits_grid.addWidget(self.r0_label, 0, 0)
        self.int_limits_grid.addWidget(self.r0_input, 0, 1)
        self.int_limits_grid.addWidget(self.rmin_label, 0, 2)
        self.int_limits_grid.addWidget(self.rmin_input, 0, 3)
        self.int_limits_grid.addWidget(self.rpmax_label, 1, 0)
        self.int_limits_grid.addWidget(self.rpmax_input, 1, 1)
        self.int_limits_grid.addWidget(self.rmax_label, 1, 2)
        self.int_limits_grid.addWidget(self.rmax_input, 1, 3)
        
        # Grid for N1 results
        self.int_results_grid = QGridLayout()
        self.int_results_grid.setContentsMargins(0, 0, 0, 0)
        self.int_results_grid.setSpacing(15)
        self.int_results_grid.addWidget(self.Na_output_label, 0, 0)
        self.int_results_grid.addWidget(self.Na_output, 0, 1)
        self.int_results_grid.addWidget(self.Nb_output_label, 1, 0)
        self.int_results_grid.addWidget(self.Nb_output, 1, 1)
        self.int_results_grid.addWidget(self.Nc_output_label, 2, 0)
        self.int_results_grid.addWidget(self.Nc_output, 2, 1)
        
        
        self.int_limits_layout.addWidget(self.find_limits_btn)
        self.int_limits_layout.addLayout(self.int_limits_grid)
        
        self.int_results_layout.addWidget(self.calc_N_btn)
        self.int_results_layout.addLayout(self.int_results_grid)


        self.outer_layout.addLayout(self.int_limits_layout)
        self.outer_layout.addLayout(self.int_results_layout)

        
        self.setLayout(self.outer_layout)
            
    def create_signals(self):
        pass
                       
    
class PolyatomicGroupBox(QGroupBox):
    
    def __init__(self, *args):
        super(PolyatomicGroupBox, self).__init__(*args)
        self.setTitle('Fit Gaussians (Polyatomic)')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        self.setCheckable(True)
                           
        self.create_widgets()
        self.style_widgets()               
        self.create_layout()
        self.create_signals()
        
    def create_widgets(self):
        
        self.r0_label = QLabel('r<sub>0</sub>')
        self.r0_input = QDoubleSpinBox()
        
        self.rpmax_label = QLabel('r\'<sub>max</sub>')
        self.rpmax_input = QDoubleSpinBox()
        
        self.rmax_label = QLabel('r<sub>max</sub>')
        self.rmax_input = QDoubleSpinBox()
        
        self.rmin_label = QLabel('r<sub>min</sub>')
        self.rmin_input = QDoubleSpinBox()
        
        self.find_limits_btn = QPushButton('Auto-Refine Int. Limits')
        
        self.calc_N_btn = QPushButton('Calc Coordination Number')
                
        self.Na_output_label = QLabel('N<sub>a</sub> : ')
        self.Na_output = QLineEdit()
        
        self.Nb_output_label = QLabel('N<sub>b</sub> : ')
        self.Nb_output = QLineEdit()
        
        self.Nc_output_label = QLabel('N<sub>c</sub> : ')
        self.Nc_output = QLineEdit()
    
    def style_widgets(self):
        
        self.r0_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.rpmax_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.rmax_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.rmin_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        
        self.Na_output_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.Nb_output_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.Nc_output_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        
        # May need to change max widths
        self.r0_input.setMaximumWidth(82)
        self.rpmax_input.setMaximumWidth(82)
        self.rmax_input.setMaximumWidth(82)
        self.rmin_input.setMaximumWidth(82)
        
        self.Na_output.setMaximumWidth(100)
        self.Nb_output.setMaximumWidth(100)
        self.Nc_output.setMaximumWidth(100)
        
        # Set results as read-only
        self.Na_output.setReadOnly(True)
        self.Nb_output.setReadOnly(True)
        self.Nc_output.setReadOnly(True)
        
        # Set input spinbox settings
        # Might need to set minimum width also to stop scale miss
        # Does dblsb need a validator? - check and the ndelete this comment
        self.r0_input.setSingleStep(0.01)
        self.r0_input.setDecimals(3)
        self.r0_input.setAlignment(Qt.AlignRight)
        self.rpmax_input.setSingleStep(0.01)
        self.rpmax_input.setDecimals(3)
        self.rpmax_input.setAlignment(Qt.AlignRight)
        self.rmax_input.setSingleStep(0.01)
        self.rmax_input.setDecimals(3)
        self.rmax_input.setAlignment(Qt.AlignRight)
        self.rmin_input.setSingleStep(0.01)
        self.rmin_input.setDecimals(3)
        self.rmin_input.setAlignment(Qt.AlignRight)        
        
        # Set tooltips
        self.r0_label.setToolTip('Leading edge of first peak in RDF(r)')
        self.r0_input.setToolTip('Leading edge of first peak in RDF(r)')
        
        self.rmin_input.setToolTip('First minimum after first peak in RDF(r)')
        self.rmin_label.setToolTip('First minimum after first peak in RDF(r)')
        
        self.rpmax_label.setToolTip('Centre of first peak in T(r)')
        self.rpmax_input.setToolTip('Centre of first peak in T(r)')
        
        self.rmax_label.setToolTip('Centre of first peak in RDF(r)')
        self.rmax_input.setToolTip('Centre of first peak in RDF(r)')
        
        self.Na_output_label.setToolTip('N<sub>a</sub> = 2 ∫ 4πρr[rg(r)]<sub>sym</sub> dr | r = r<sub>0</sub> to r = r\'<sub>max</sub>')
        self.Na_output.setToolTip('N<sub>a</sub> = 2 ∫ 4πρr[rg(r)]<sub>sym</sub> dr | r = r<sub>0</sub> to r = r\'<sub>max</sub>')
        
        self.Nb_output_label.setToolTip('N<sub>b</sub> = 2 ∫ 4πρ[r<sup>2</sup>g(r)]<sub>sym</sub> dr | r = r<sub>0</sub> to r = r<sub>max</sub>')
        self.Nb_output.setToolTip('N<sub>b</sub> = 2 ∫ 4πρ[r<sup>2</sup>g(r)]<sub>sym</sub> dr | r = r<sub>0</sub> to r = r<sub>max</sub>')
        
        self.Nc_output_label.setToolTip('N<sub>c</sub> = ∫ 4πρr<sup>2</sup>g(r) dr | r = r<sub>0</sub> to r = r<sub>min</sub>')
        self.test = self.Nc_output.setToolTip('N<sub>c</sub> = ∫ 4πρr<sup>2</sup>g(r) dr | r = r<sub>0</sub> to r = r<sub>min</sub>')
        
        self.setStyleSheet(" QToolTip{ font: 11pt;}")
        # Set ToolTip style sheets separately
        self.Na_output_label.setStyleSheet(" QToolTip{ min-width: 480px; font: 18pt;}")
        self.Na_output.setStyleSheet(" QToolTip{ min-width: 480px; font: 18pt;}")
        self.Nb_output_label.setStyleSheet(" QToolTip{ min-width: 475px; font: 18pt;}")
        self.Nb_output.setStyleSheet(" QToolTip{ min-width: 475px; font: 18pt;}")
        self.Nc_output_label.setStyleSheet(" QToolTip{ min-width: 400px; font: 18pt;}")
        self.Nc_output.setStyleSheet(" QToolTip{ min-width: 400px; font: 18pt;}")  

    
    def create_layout(self):
        self.outer_layout = QVBoxLayout()
        self.outer_layout.setContentsMargins(25, 10, 25, 7)
        self.outer_layout.setSpacing(30)

        self.int_limits_layout = QVBoxLayout()
        self.int_limits_layout.setContentsMargins(0, 0, 0, 0)
        self.int_limits_layout.setSpacing(15)
        
        self.int_results_layout = QVBoxLayout()
        self.int_results_layout.setContentsMargins(0, 0, 0, 0)
        self.int_results_layout.setSpacing(15)        


        # Grid for integration limits
        self.int_limits_grid = QGridLayout()
        self.int_limits_grid.setContentsMargins(0, 0, 0, 0)
        self.int_limits_grid.setSpacing(15)
        self.int_limits_grid.addWidget(self.r0_label, 0, 0)
        self.int_limits_grid.addWidget(self.r0_input, 0, 1)
        self.int_limits_grid.addWidget(self.rmin_label, 0, 2)
        self.int_limits_grid.addWidget(self.rmin_input, 0, 3)
        self.int_limits_grid.addWidget(self.rpmax_label, 1, 0)
        self.int_limits_grid.addWidget(self.rpmax_input, 1, 1)
        self.int_limits_grid.addWidget(self.rmax_label, 1, 2)
        self.int_limits_grid.addWidget(self.rmax_input, 1, 3)
        
        # Grid for N1 results
        self.int_results_grid = QGridLayout()
        self.int_results_grid.setContentsMargins(0, 0, 0, 0)
        self.int_results_grid.setSpacing(15)
        self.int_results_grid.addWidget(self.Na_output_label, 0, 0)
        self.int_results_grid.addWidget(self.Na_output, 0, 1)
        self.int_results_grid.addWidget(self.Nb_output_label, 1, 0)
        self.int_results_grid.addWidget(self.Nb_output, 1, 1)
        self.int_results_grid.addWidget(self.Nc_output_label, 2, 0)
        self.int_results_grid.addWidget(self.Nc_output, 2, 1)
        
        
        self.int_limits_layout.addWidget(self.find_limits_btn)
        self.int_limits_layout.addLayout(self.int_limits_grid)
        
        self.int_results_layout.addWidget(self.calc_N_btn)
        self.int_results_layout.addLayout(self.int_results_grid)


        self.outer_layout.addLayout(self.int_limits_layout)
        self.outer_layout.addLayout(self.int_results_layout)

        
        self.setLayout(self.outer_layout)
            
    def create_signals(self):
        pass
