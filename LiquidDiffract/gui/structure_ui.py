# -*- coding: utf-8 -*-
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2021, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

# What imports do I need??

import os
import numpy as np
from PyQt5.QtCore import Qt, pyqtSignal, QSignalMapper, QObject
from PyQt5.QtGui import QDoubleValidator, QFont
from PyQt5.QtWidgets import QWidget, QFrame, QGridLayout, QVBoxLayout, \
                            QHBoxLayout, QGroupBox, QPushButton, QRadioButton, \
                            QLineEdit, QDoubleSpinBox, QLabel, QScrollArea, \
                            QMessageBox, QCheckBox, QSplitter, QToolTip, QComboBox, QStyle
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
        self.clear_data()

        self.create_signals()

    def create_signals(self):

        self.structure_config_widget.plot_view_gb.rdf_btn.toggled.connect(self.toggle_plot_view)

        self.structure_config_widget.monatomic_gb.r0_input.valueChanged.connect(self.update_plot_data)
        self.structure_config_widget.monatomic_gb.rpmax_input.valueChanged.connect(self.update_plot_data)
        self.structure_config_widget.monatomic_gb.rmax_input.valueChanged.connect(self.update_plot_data)
        self.structure_config_widget.monatomic_gb.rmin_input.valueChanged.connect(self.update_plot_data)

        self.int_limit_signalMapper = QSignalMapper()
        self.int_limit_signalMapper.mapped[QObject].connect(self.update_int_limits)

    def create_int_limit_signals(self):

        try:
            self.int_limit_signalMapper.setMapping(self.structure_plot_widget.r0_line_rdf, self.structure_plot_widget.r0_line_rdf)
            self.int_limit_signalMapper.setMapping(self.structure_plot_widget.r0_line_tr, self.structure_plot_widget.r0_line_tr)
            self.int_limit_signalMapper.setMapping(self.structure_plot_widget.rpmax_line_rdf, self.structure_plot_widget.rpmax_line_rdf)
            self.int_limit_signalMapper.setMapping(self.structure_plot_widget.rpmax_line_tr, self.structure_plot_widget.rpmax_line_tr)
            self.int_limit_signalMapper.setMapping(self.structure_plot_widget.rmax_line_rdf, self.structure_plot_widget.rmax_line_rdf)
            self.int_limit_signalMapper.setMapping(self.structure_plot_widget.rmax_line_tr, self.structure_plot_widget.rmax_line_tr)
            self.int_limit_signalMapper.setMapping(self.structure_plot_widget.rmin_line_rdf, self.structure_plot_widget.rmin_line_rdf)
            self.int_limit_signalMapper.setMapping(self.structure_plot_widget.rmin_line_tr, self.structure_plot_widget.rmin_line_tr)

            self.structure_plot_widget.r0_line_rdf.sigDragged.connect(self.int_limit_signalMapper.map)
            self.structure_plot_widget.r0_line_tr.sigDragged.connect(self.int_limit_signalMapper.map)
            self.structure_plot_widget.rpmax_line_rdf.sigDragged.connect(self.int_limit_signalMapper.map)
            self.structure_plot_widget.rpmax_line_tr.sigDragged.connect(self.int_limit_signalMapper.map)
            self.structure_plot_widget.rmax_line_rdf.sigDragged.connect(self.int_limit_signalMapper.map)
            self.structure_plot_widget.rmax_line_tr.sigDragged.connect(self.int_limit_signalMapper.map)
            self.structure_plot_widget.rmin_line_rdf.sigDragged.connect(self.int_limit_signalMapper.map)
            self.structure_plot_widget.rmin_line_tr.sigDragged.connect(self.int_limit_signalMapper.map)

        except AttributeError:
            return


    def update_int_limits(self, _sender):

        _new_pos = _sender.value()
        _limit = _sender.name()
        if _limit == 'r0':
            self.structure_config_widget.monatomic_gb.r0_input.setValue(_new_pos)
        elif _limit == 'rpmax':
            self.structure_config_widget.monatomic_gb.rpmax_input.setValue(_new_pos)
        elif _limit == 'rmax':
            self.structure_config_widget.monatomic_gb.rmax_input.setValue(_new_pos)
        elif _limit == 'rmin':
            self.structure_config_widget.monatomic_gb.rmin_input.setValue(_new_pos)

        else:
            raise ValueError


    def clear_data(self):

        self.data = {'rdf_x': np.asarray([]), 'rdf_y': np.asarray([]),
                     'tr_x': np.asarray([]), 'tr_y': np.asarray([]),
                     'fr_x': np.asarray([]), 'fr_y': np.asarray([]),
                     'rho': None, 'composition': None, 'sq_x': np.asarray([])}

    def set_atoms(self):

        _atom_list = self.data['composition'].keys()
        if self.structure_config_widget.polyatomic_gb.atom_list == _atom_list:
            pass

        else:
            self.structure_config_widget.polyatomic_gb.atom_list = _atom_list
            for _peak in self.structure_config_widget.polyatomic_gb.peak_dict.values():
                _peak.populate_atom_list(_atom_list)


    def plot_data(self):

        self.update_plot_data()
        self.structure_plot_widget.update_plot_windows(self.data)

    def update_plot_data(self):
        self.data['r0'] = self.structure_config_widget.monatomic_gb.r0_input.value()
        self.data['rpmax'] = self.structure_config_widget.monatomic_gb.rpmax_input.value()
        self.data['rmax'] = self.structure_config_widget.monatomic_gb.rmax_input.value()
        self.data['rmin'] = self.structure_config_widget.monatomic_gb.rmin_input.value()

        self.structure_plot_widget.update_plots(self.data)

        self.create_int_limit_signals()

    def toggle_plot_view(self):

        if self.structure_config_widget.plot_view_gb.rdf_btn.isChecked():
            self.structure_plot_widget.pg_layout_widget_tr.setVisible(False)
            self.structure_plot_widget.pg_layout_widget_rdf.setVisible(True)


        elif self.structure_config_widget.plot_view_gb.tr_btn.isChecked():
            self.structure_plot_widget.pg_layout_widget_rdf.setVisible(False)
            self.structure_plot_widget.pg_layout_widget_tr.setVisible(True)
        else:
            pass


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
        self.dr_btn = QRadioButton('D(r)')

    def style_widgets(self):

        self.rdf_btn.setToolTip('RDF(r) = 4πρr<sup>2</sup>g(r)')
        self.tr_btn.setToolTip('T(r) = RDF(r) / r')
        self.dr_btn.setToolTip('D(r) = 4πρr[g(r) - 1]')
        #self.residuals_btn.setToolTip('Plot residuals from fit?')

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

        self.atom_list = []
        self.init_peaks()

        self.create_signals()

    def create_widgets(self):

        self.fit_limit_label = QLabel('Region to fit: ')
        self.min_limit_label = QLabel('min: ')
        self.max_limit_label = QLabel('max: ')

        self.min_limit_input = QLineEdit()
        self.max_limit_input = QLineEdit()

        self.add_peak_btn = QPushButton('Add Peak')

    def style_widgets(self):

        self.fit_limit_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.min_limit_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.max_limit_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)

        self.min_limit_input.setMaximumWidth(82)
        self.max_limit_input.setMaximumWidth(82)

        # Set validators
        self.min_limit_input.setValidator(QDoubleValidator())
        self.max_limit_input.setValidator(QDoubleValidator())

    def create_layout(self):
        self.outer_layout = QVBoxLayout()
        self.outer_layout.setContentsMargins(25, 10, 25, 7)
        self.outer_layout.setSpacing(30)

        self.fit_settings_grid = QGridLayout()
        self.fit_settings_grid.setContentsMargins(0, 0, 0, 0)
        self.fit_settings_grid.setSpacing(15)
        
        self.fit_settings_grid.addWidget(self.fit_limit_label, 0, 0)
        self.fit_settings_grid.addWidget(self.min_limit_label, 0, 1)
        self.fit_settings_grid.addWidget(self.min_limit_input, 0, 2)
        self.fit_settings_grid.addWidget(self.max_limit_label, 1, 1)
        self.fit_settings_grid.addWidget(self.max_limit_input, 1, 2)

        self.outer_layout.addLayout(self.fit_settings_grid)
        self.outer_layout.addWidget(self.add_peak_btn)

        self.setLayout(self.outer_layout)

    def init_peaks(self):

        self.peak_count = 0
        self.peak_dict = {}

    def create_signals(self):

        self.add_peak_btn.clicked.connect(self.add_peak)

        self.del_peak_signalMapper = QSignalMapper()
        self.del_peak_signalMapper.mapped[int].connect(self.del_peak)

        self.rename_peak_signalMapper = QSignalMapper()
        self.rename_peak_signalMapper.mapped[int].connect(self.rename_peak)

    def add_peak(self):

        if len(self.peak_dict) == 0:
            self.init_peaks()

        _new_peak = GaussianPeakGroupBox()
        _new_peak._unique_idx = self.peak_count
        _new_peak.title.setText(_new_peak.title.text() + ' #' + str(_new_peak._unique_idx+1) + ' ')
        _new_peak.populate_atom_list(self.atom_list)
        self.peak_dict[_new_peak._unique_idx] = _new_peak
        self.outer_layout.addWidget(_new_peak)

        # Link peak delete button to signal mapper
        self.del_peak_signalMapper.setMapping(self.peak_dict[_new_peak._unique_idx].del_peak_btn, self.peak_dict[_new_peak._unique_idx]._unique_idx)
        self.peak_dict[_new_peak._unique_idx].del_peak_btn.clicked.connect(self.del_peak_signalMapper.map)

        # Link atom inputs to name of peak
        self.rename_peak_signalMapper.setMapping(self.peak_dict[_new_peak._unique_idx].alpha_input, self.peak_dict[_new_peak._unique_idx]._unique_idx)
        self.peak_dict[_new_peak._unique_idx].alpha_input.currentTextChanged.connect(self.rename_peak_signalMapper.map)

        self.rename_peak_signalMapper.setMapping(self.peak_dict[_new_peak._unique_idx].beta_input, self.peak_dict[_new_peak._unique_idx]._unique_idx)
        self.peak_dict[_new_peak._unique_idx].beta_input.currentTextChanged.connect(self.rename_peak_signalMapper.map)

        self.peak_count += 1

        del _new_peak

    def del_peak(self, _peak_index):

        # Remove from layout
        self.outer_layout.removeWidget(self.peak_dict[_peak_index])
        # Delete C/C++ object
        self.peak_dict[_peak_index].deleteLater()
        # Delete Py object
        self.peak_dict[_peak_index] = None
        # Remove peak entry from dict
        del self.peak_dict[_peak_index]
        del _peak_index

    def rename_peak(self, _peak_index):

        _target_peak = self.peak_dict[_peak_index]
        _alpha = _target_peak.alpha_input.currentText()
        _beta = _target_peak.beta_input.currentText()
        _current_title = _target_peak.title.text().split('(')[0]
        if (_alpha == '') & (_beta == ''):
            _target_peak.title.setText(_current_title)
        else:
            _title = _current_title + '(' + _alpha + '-' + _beta + ')'
            _target_peak.title.setText(_title)


class GaussianPeakGroupBox(QFrame):

    def __init__(self, *args):
        super(GaussianPeakGroupBox, self).__init__(*args)

        self.create_widgets()
        self.style_widgets()
        self.create_layout()
        #self.create_signals()

    def create_widgets(self):

        self.title = QLabel('Gaussian Peak')

        self.hline = QFrame()
        self.hline.setFrameShape(QFrame.HLine)
        self.hline.setFrameShadow(QFrame.Sunken)
        self.hline.setObjectName('hline')

        self.del_peak_btn = QPushButton()

        self.alpha_label = QLabel('alpha')
        self.beta_label = QLabel('beta')

        self.alpha_input = QComboBox()
        self.alpha_input.setCurrentIndex(-1)

        self.beta_input = QComboBox()
        self.beta_input.setCurrentIndex(-1)

        self.refine_label = QLabel('Refine?')

        self.N_label = QLabel('N: ')
        self.N_input = QLineEdit()
        self.N_refine = QCheckBox()

        self.r_label = QLabel('r: ')
        self.r_input = QLineEdit()
        self.r_refine = QCheckBox()

        self.s_label = QLabel('sigma: ')
        self.s_input = QLineEdit()
        self.s_refine = QCheckBox()

    def style_widgets(self):

        self.title.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        self.del_peak_btn.setIcon(self.style().standardIcon(getattr(QStyle, 'SP_TitleBarCloseButton')))
        self.del_peak_btn.setMaximumWidth(25)
        self.del_peak_btn.setFlat(True)

        self.alpha_label.setAlignment(Qt.AlignVCenter | Qt.AlignCenter)
        self.beta_label.setAlignment(Qt.AlignVCenter | Qt.AlignCenter)

        self.refine_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        self.N_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.r_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.s_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)

        self.N_input.setMaximumWidth(82)
        self.r_input.setMaximumWidth(82)
        self.s_input.setMaximumWidth(82)

        # Set validators
        self.N_input.setValidator(QDoubleValidator())
        self.r_input.setValidator(QDoubleValidator())
        self.s_input.setValidator(QDoubleValidator())

    def create_layout(self):

        self.vlayout = QVBoxLayout()
        self.vlayout.setContentsMargins(0, 5, 0, 0)
        self.vlayout.setSpacing(15)

        self.title_layout = QVBoxLayout()
        self.title_layout.setContentsMargins(0, 0, 0, 0)
        self.title_layout.setSpacing(5)
        self.title_inner_layout = QHBoxLayout()
        self.title_inner_layout.setContentsMargins(0, 0, 0, 0)
        self.title_inner_layout.setSpacing(0)

        self.atoms_grid_layout = QGridLayout()
        self.atoms_grid_layout.setContentsMargins(0, 0, 0, 0)
        self.atoms_grid_layout.setSpacing(10)

        self.title_inner_layout.addWidget(self.title)
        self.title_inner_layout.addWidget(self.del_peak_btn)
        self.title_layout.addLayout(self.title_inner_layout)
        self.title_layout.addWidget(self.hline)

        self.atoms_grid_layout.addWidget(self.alpha_label, 0, 0)
        self.atoms_grid_layout.addWidget(self.beta_label, 0, 1)
        self.atoms_grid_layout.addWidget(self.alpha_input, 1, 0)
        self.atoms_grid_layout.addWidget(self.beta_input, 1, 1)

        self.params_grid_layout = QGridLayout()
        self.params_grid_layout.setContentsMargins(0, 0, 0, 0)
        self.params_grid_layout.setSpacing(10)

        self.params_grid_layout.addWidget(self.refine_label, 0, 2)
        self.params_grid_layout.addWidget(self.N_label, 1, 0)
        self.params_grid_layout.addWidget(self.N_input, 1, 1)
        self.params_grid_layout.addWidget(self.N_refine, 1, 2)
        self.params_grid_layout.addWidget(self.r_label, 2, 0)
        self.params_grid_layout.addWidget(self.r_input, 2, 1)
        self.params_grid_layout.addWidget(self.r_refine, 2, 2)
        self.params_grid_layout.addWidget(self.s_label, 3, 0)
        self.params_grid_layout.addWidget(self.s_input, 3, 1)
        self.params_grid_layout.addWidget(self.s_refine, 3, 2)

        self.vlayout.addLayout(self.title_layout)
        self.vlayout.addLayout(self.atoms_grid_layout)
        self.vlayout.addLayout(self.params_grid_layout)

        self.setLayout(self.vlayout)

    def populate_atom_list(self, _atom_list):
        self.alpha_input.clear()
        self.beta_input.clear()
        self.alpha_input.addItems(_atom_list)
        self.beta_input.addItems(_atom_list)
        self.alpha_input.setCurrentIndex(-1)
        self.beta_input.setCurrentIndex(-1)