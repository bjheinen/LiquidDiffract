# -*- coding: utf-8 -*-
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2018-2024, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

import os
import numpy as np
from scipy.optimize import minimize
from qtpy.QtCore import Qt, Signal
from qtpy.QtGui import QDoubleValidator
from qtpy.QtWidgets import QWidget, QFrame, QGridLayout, QVBoxLayout, \
                            QHBoxLayout, QGroupBox, QPushButton, QLineEdit, \
                            QDoubleSpinBox, QLabel, QScrollArea, QMessageBox, \
                            QCheckBox, QSplitter, QButtonGroup, QRadioButton, \
                            QAbstractSpinBox
# LiquidDiffract imports
from LiquidDiffract.gui import plot_widgets
from LiquidDiffract.gui import utility
from LiquidDiffract.core import data_utils
from LiquidDiffract.core.core import calc_self_shielding
from LiquidDiffract.version import __appname__, __version__


class BkgUI(QWidget):

    # plots_changed is handled in main_widget.MainContainer and allows the
    # next tab to be updated
    plots_changed = Signal()
    # file_name_changed connects to main_widget.MainContainer.update_filename
    # this method is used to set the working directory based on the data_file
    # that is loaded
    file_name_changed = Signal()

    def __init__(self, parent):
        super(BkgUI, self).__init__(parent)
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

        self.setLayout(self.layout)

        self.data_file = os.path.abspath(os.getcwd())
        self.file_name_changed.emit()

        __a = np.asarray([])
        self.data = {'data_raw_x': __a, 'data_raw_y':   __a,
                     'data_x':     __a, 'data_y':       __a,
                     'bkg_raw_x':  __a, 'bkg_raw_y':    __a,
                     'bkg_x':      __a, 'bkg_y':        __a,
                     'bkg_y_sc':   __a, 'bkg_raw_y_sc': __a,
                     'cor_x':      __a, 'cor_y':        __a,
                     'data_correction': 0
                     }
        self.bkg_file = None
        self.data_file = None

        self.create_signals()

    def create_signals(self):
        self.bkg_config_widget.data_files_gb.load_data_btn.clicked.connect(self.load_data)
        self.bkg_config_widget.data_files_gb.load_bkg_btn.clicked.connect(self.load_bkg)
        self.bkg_config_widget.data_files_gb.plot_raw_check.stateChanged.connect(self.plot_data)
        self.bkg_config_widget.data_files_gb.plot_log_check.stateChanged.connect(self.plot_data)
        self.bkg_config_widget.bkg_subtract_gb.bkg_sub_btn.clicked.connect(self.sub_bkg)
        self.bkg_config_widget.bkg_subtract_gb.scale_sb.valueChanged.connect(self.plot_data)
        self.bkg_config_widget.bkg_subtract_gb.toggled.connect(self.sub_bkg)
        self.bkg_config_widget.bkg_subtract_gb.auto_sc_btn.clicked.connect(self.auto_scale_bkg)
        self.bkg_config_widget.data_files_gb.dq_input.editingFinished.connect(self.dq_changed)
        self.bkg_config_widget.data_corrections_gb.zero_shift_check.stateChanged.connect(self.plot_data)
        self.bkg_config_widget.data_corrections_gb.zero_shift_button_group.buttonClicked.connect(self.plot_data)
        self.bkg_config_widget.data_corrections_gb.shift_sb.valueChanged.connect(self.plot_data)
        self.bkg_config_widget.data_corrections_gb.attenuation_slab_check.stateChanged.connect(self.plot_data)
        self.bkg_config_widget.data_corrections_gb.wavelength_input.textChanged.connect(self.plot_data)
        self.bkg_config_widget.data_corrections_gb.attenuation_factor_input.textChanged.connect(self.plot_data)
        self.bkg_config_widget.data_corrections_gb.thickness_input.textChanged.connect(self.plot_data)
        self.bkg_config_widget.data_corrections_gb.rotation_input.valueChanged.connect(self.plot_data)
        self.bkg_config_widget.data_corrections_gb.inspect_btn.clicked.connect(self.inspect_attenuation_factor)

    def rebin_data(self, bkg='check', x_lim=None, suppress_plots=False):
        '''
        Helper function to rebin background/data arrays when necessary.
        Because the S(Q) array is padded before fourier transform operations
        in core.core the size should be checked here. Checking takes place
        before rebinning as this could hang on very large files.
        '''
        try:
            _dx = np.float64(self.bkg_config_widget.data_files_gb.dq_input.text())
        except ValueError:
            _dx = 0
        if bkg == 1:
            if self.data['bkg_raw_x'][-1] > ((2**self.fft_N / 2) * _dx):
                raise RuntimeWarning('Dataset size exceeds 2**N!')
            self.data['bkg_x'], self.data['bkg_y'] = data_utils.rebin_data(self.data['bkg_raw_x'], self.data['bkg_raw_y'], dx=_dx, x_lim=x_lim)
            if suppress_plots == False:
                self.plot_data()

        elif bkg == 0:
            if self.data['data_raw_x'][-1] > ((2**self.fft_N / 2) * _dx):
                 raise RuntimeWarning('Dataset size exceeds 2**N!')
            self.data['data_x'], self.data['data_y'] = data_utils.rebin_data(self.data['data_raw_x'], self.data['data_raw_y'], dx=_dx, x_lim=x_lim)
            # Delete any processed data when loading new file
            self.data['cor_x'] = np.asarray([])
            self.data['cor_y'] = np.asarray([])
            if not self.bkg_config_widget.bkg_subtract_gb.isChecked():
                self.sub_bkg()
            self.plots_changed.emit()
            if suppress_plots == False:
                self.plot_data()

        elif bkg == 'check':
            if self.data['data_raw_x'].size and self.data['bkg_raw_x'].size:
                _x_lim = (0, min(self.data['bkg_raw_x'][-1], self.data['data_raw_x'][-1]))
            else:
                _x_lim = None
            if self.data['bkg_raw_x'].size and self.bkg_file != None:
                 self.rebin_data(bkg=1, x_lim=_x_lim, suppress_plots=True)
            if self.data['data_raw_x'].size and self.data_file != None:
                 self.rebin_data(bkg=0, x_lim=_x_lim, suppress_plots=True)
            if suppress_plots == False:
                self.plot_data()

    def load_data(self):
        __file_name = utility.get_filename(io='open')
        if __file_name:
            self.data_file = __file_name
        else:
            return
        try:
            self.data['data_raw_x'], self.data['data_raw_y'] = np.loadtxt(self.data_file, unpack=True)
            # Convert to angstroms if data_units == 1 (nm)
            if self.data_units:
                self.data['data_raw_x'] /= 10.0
        except ValueError:
            try:
                _header_len = self.check_file_header(self.data_file)
            except UnicodeDecodeError:
                self.data_file = None
                self.load_file_error()
                return
            if _header_len:
                self.data['data_raw_x'], self.data['data_raw_y'] = np.loadtxt(self.data_file, unpack=True, skiprows=_header_len)
                # Convert to angstroms if data_units == 1 (nm)
                if self.data_units:
                    self.data['data_raw_x'] /= 10.0
            else:
                self.data_file = None
                return
        try:
            self.rebin_data(bkg=0)
        except RuntimeWarning:
             self.oversize_file_error()
             self.data_file = None
             return
        self.bkg_config_widget.data_files_gb.data_filename_lbl.setText(self.data_file.split('/')[-1])
        print(f'Data file: {self.data_file}')
        self.file_name_changed.emit()

    def load_bkg(self):
        __file_name = utility.get_filename(io='open', caption='Load Background File')
        if __file_name:
            self.bkg_file = __file_name
        else:
            return
        try:
            self.data['bkg_raw_x'], self.data['bkg_raw_y'] = np.loadtxt(self.bkg_file, unpack=True)
            if self.data_units:
                self.data['bkg_raw_x'] /= 10.0
        except ValueError:
            try:
                _header_len = self.check_file_header(self.bkg_file)
            except UnicodeDecodeError:
                self.bkg_file = None
                self.load_file_error()
                return
            if _header_len:
                self.data['bkg_raw_x'], self.data['bkg_raw_y'] = np.loadtxt(self.bkg_file, unpack=True, skiprows=_header_len)
                # Convert to angstroms if data_units == 1 (nm)
                if self.data_units:
                    self.data['bkg_raw_x'] /= 10.0
            else:
                self.bkg_file = None
                return
        try:
            self.rebin_data(bkg=1)
        except RuntimeWarning:
            self.oversize_file_error()
            self.bkg_file = None
            return
        self.bkg_config_widget.data_files_gb.bkg_filename_lbl.setText(self.bkg_file.split('/')[-1])
        print(f'Background File: {self.bkg_file}')

    def check_file_header(self, _fname):
        _check_dialog = utility.CheckFileDialog(_fname)
        _check_dialog.resize(int(self.screen_size[0]*0.6), int(self.screen_size[1]*0.5))
        if _check_dialog.exec() == utility.CheckFileDialog.Accepted:
            return _check_dialog.get_header_len()
        else:
            return None

    def plot_data(self):
        if self.data['bkg_y'].size:
            # First scale bkg data
            _bkg_scaling = self.bkg_config_widget.bkg_subtract_gb.scale_sb.value()
            self.data['bkg_raw_y_sc'] = self.data['bkg_raw_y'] * _bkg_scaling
            self.data['bkg_y_sc'] = self.data['bkg_y'] * _bkg_scaling
        if self.data['cor_x'].size:
            # Only re-subtract if already subtract button clicked
            if self.bkg_config_widget.bkg_subtract_gb.isChecked() and self.bkg_file is not None:
                # Subtract background from raw data if option chosen
                try:
                    self.data['cor_y'] = self.data['data_y'] - self.data['bkg_y_sc']
                except ValueError:
                    self.rebin_data(suppress_plots=True)
                    # Check if data_x and bkg_x now the same
                    if np.allclose(self.data['data_x'], self.data['bkg_x']):
                        self.data['cor_x'] = self.data['data_x']
                        self.plot_data()
                        return
                    else:
                        self.bkg_match_error()
                        return
            else:
                self.data['cor_y'] = self.data['data_y']
        _plot_raw = self.bkg_config_widget.data_files_gb.plot_raw_check.isChecked()
        _plot_log = self.bkg_config_widget.data_files_gb.plot_log_check.isChecked()
        self.data_corrections()
        self.bkg_plot_widget.update_plots(self.data, _plot_raw, _plot_log)
        # emit signal that data has changed to be picked up by second tab (optim_ui)
        self.plots_changed.emit()

    def sub_bkg(self):
        # Cor x = data x
        self.data['cor_x'] = self.data['data_x']
        self.plot_data()

    def auto_scale_bkg(self):
        if self.bkg_file is None:
            self.missing_bkg_file_error()
            return
        if len(self.data['data_y']) != len(self.data['bkg_y']):
            self.rebin_data(suppress_plots=True)
            if len(self.data['data_y']) != len(self.data['bkg_y']):
                self.bkg_match_error()
                return
        bkg_scaling = minimize(data_utils.bkg_scaling_residual, 1,
                               args=(self.data['data_y'], self.data['bkg_y']),
                               method='nelder-mead',
                               options={'xtol': 1e-8, 'disp': False})
        bkg_scaling = bkg_scaling.x
        self.bkg_config_widget.bkg_subtract_gb.scale_sb.setValue(bkg_scaling)

    def data_corrections(self):
        self.data['shift_correction'] = 0
        self.data['attenuation_correction'] = 1.0
        # Skip if no data present
        if not self.data['cor_x'].size:
            return
        if self.bkg_config_widget.data_corrections_gb.zero_shift_check.isChecked():
            _shift_correction = self.zero_shift_data()
        else:
            _shift_correction = 0
        self.data['shift_correction'] = _shift_correction
        if self.bkg_config_widget.data_corrections_gb.attenuation_slab_check.isChecked():
            _attenuation_correction = self.attenuation_correction_slab()
        else:
            _attenuation_correction = 1.0
        self.data['attenuation_correction'] = _attenuation_correction

    def zero_shift_data(self):
        if self.bkg_config_widget.data_corrections_gb.zero_first_val_btn.isChecked():
            _shift = self.data['cor_y'][0]
        elif self.bkg_config_widget.data_corrections_gb.zero_min_val_btn.isChecked():
            _shift = np.min(self.data['cor_y'])
        elif self.bkg_config_widget.data_corrections_gb.custom_shift_btn.isChecked():
            _shift = -1 * self.bkg_config_widget.data_corrections_gb.shift_sb.value()
        else:
            raise NotImplementedError()
            return
        self.data['cor_y'] = (self.data['cor_y'] - _shift) + 1e-22
        return _shift * -1

    def attenuation_correction_slab(self):
        # Get attenuation info
        try:
            _wavelength = float(self.bkg_config_widget.data_corrections_gb.wavelength_input.text())
            _attenuation_factor = float(self.bkg_config_widget.data_corrections_gb.attenuation_factor_input.text())
            _thickness = float(self.bkg_config_widget.data_corrections_gb.thickness_input.text())
            _rotation = self.bkg_config_widget.data_corrections_gb.rotation_input.value()
        except ValueError:
            return 1.0
            # do i need an error message?
        # Convert attenuation factor to m^-1
        _attenuation_factor *= 100.0
        # Convert thickness to metres
        _thickness /= 1000.0
        # Convert rotation from normal (beta), to rotation from beam direction (alpha)
        _alpha = 90 - _rotation
        # core.calc_self_shielding(Q, mu, alpha, thickness, wavelength)
        _attenuation = calc_self_shielding(self.data['cor_x'], _attenuation_factor, _alpha, _thickness, _wavelength)
        # I_measured = I_actual * A_s,s --> I_actual = I_measured / A_s,s (or * 1/A_s,s)
        self.data['cor_y'] = self.data['cor_y'] / _attenuation
        return 1 / _attenuation

    def inspect_attenuation_factor(self):
        # Return if no data
        if not self.data['cor_x'].size or not isinstance(self.data['attenuation_correction'], np.ndarray):
            return
        # Grab two-theta data
        try:
            _wavelength = float(self.bkg_config_widget.data_corrections_gb.wavelength_input.text())
        except ValueError:
            return
        _two_theta = data_utils.convert_q_space(self.data['cor_x'], _wavelength)
        _attenuation_data = self.data['cor_x'], _two_theta, 1/self.data['attenuation_correction']
        _inspect_dialog = utility.AttenuationInspectDialog(_attenuation_data)
        _inspect_dialog.resize(int(self.screen_size[0]*0.7), int(self.screen_size[1]*0.7))
        _inspect_dialog.exec()

    def dq_changed(self):
        try:
            self.rebin_data()
        except RuntimeWarning:
             self.dq_error()
             self.bkg_config_widget.data_files_gb.dq_input.setText('0.02')
             self.dq_changed()

    def data_units_changed(self):
        if self.data_units == 1:
            self.data['data_raw_x'] /= 10.0
            self.data['bkg_raw_x'] /= 10.0
        elif self.data_units == 0:
            self.data['data_raw_x'] *= 10.0
            self.data['bkg_raw_x'] *= 10.0
        else:
            raise ValueError
        self.rebin_data(bkg='check')

    def load_file_error(self):
        message = ['Error loading file!', 'Unable to load file.\nPlease check filename,\nensure header lines are commented (#),\nand data is in Q-space']
        self.warning_message(message)

    def oversize_file_error(self):
        message = ['Error loading file!', 'Re-binned array size\nis too large! Check the data-units preferences setting, increase the Q-step or increase length of Fourier transform array (N)']
        self.warning_message(message)

    def missing_bkg_file_error(self):
        message = ['No background file!', 'Please load a background file']
        self.warning_message(message)
    
    def bkg_match_error(self):
        message = ['Error subtracting background!', 'Data and background\n do not match!']
        self.warning_message(message)

    def dq_error(self):
        message = ['Error setting Q-step!', 'Re-binned array size\nis too large! Increase the Q-step or length of Fourier transform array (N)']
        self.warning_message(message)

    def warning_message(self, message):
        self.error_msg = QMessageBox()
        self.error_msg.setIcon(QMessageBox.Warning)
        self.error_msg.setStandardButtons(QMessageBox.Ok)
        self.error_msg.setText(message[0])
        self.error_msg.setInformativeText((message[1]))
        self.error_msg.setWindowTitle(__appname__ + ' v' + __version__)
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
        self.data_corrections_gb = DataCorrectionsGroupBox()
        self.vlayout.addWidget(self.data_files_gb, 1)
        self.vlayout.addWidget(self.bkg_subtract_gb, 1)
        self.vlayout.addWidget(self.data_conv_gb, 1)
        self.vlayout.addWidget(self.data_corrections_gb, 1)
        self.vlayout.addWidget(QWidget(), 3)
        self.setLayout(self.vlayout)


class DataFilesGroupBox(QGroupBox):

    def __init__(self, *args):
        super(DataFilesGroupBox, self).__init__(*args)
        self.setTitle('Data Files')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        
        self.create_widgets()
        self.style_widgets()
        self.create_layout()

    def create_widgets(self):
        self.load_data_btn = QPushButton("Load Data")
        self.data_filename_lbl = QLabel("None")

        self.load_bkg_btn = QPushButton("Load Background")
        self.bkg_filename_lbl = QLabel("None")

        self.data_lbl_frame = QScrollArea()
        self.bkg_lbl_frame = QScrollArea()

        self.plot_raw_lbl = QLabel('Plot raw (unbinned) data?')
        self.plot_raw_check = QCheckBox()
        self.plot_raw_check.setChecked(False)

        self.plot_log_label = QLabel('Plot log I(Q)?')
        self.plot_log_check = QCheckBox()
        self.plot_log_check.setChecked = False

        self.dq_label = QLabel('Re-binned resolution (Q-step):')
        self.dq_input = QLineEdit('0.02')

    def style_widgets(self):
        self.data_filename_lbl.setAlignment(Qt.AlignCenter)
        self.bkg_filename_lbl.setAlignment(Qt.AlignCenter)

        self.data_lbl_frame.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.data_lbl_frame.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.data_lbl_frame.setWidgetResizable(True)
        self.data_lbl_frame.setFrameShape(QFrame.NoFrame)

        self.bkg_lbl_frame.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.bkg_lbl_frame.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.bkg_lbl_frame.setWidgetResizable(True)
        self.bkg_lbl_frame.setFrameShape(QFrame.NoFrame)

        self.dq_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.dq_input.setAlignment(Qt.AlignRight)
        self.dq_input.setMaximumWidth(80)
        self.dq_input.setValidator(QDoubleValidator())

    def create_layout(self):
        self.grid_layout = QGridLayout()
        self.grid_layout.setContentsMargins(25, 10, 25, 7)
        self.grid_layout.setSpacing(10)
        self.data_lbl_frame.setWidget(self.data_filename_lbl)
        self.bkg_lbl_frame.setWidget(self.bkg_filename_lbl)
        self.grid_layout.addWidget(self.load_data_btn, 0, 0)
        self.grid_layout.addWidget(self.data_lbl_frame, 0, 1)
        self.grid_layout.addWidget(self.load_bkg_btn, 1, 0)
        self.grid_layout.addWidget(self.bkg_lbl_frame, 1, 1)
        self.grid_layout.addWidget(self.plot_raw_lbl, 2, 0)
        self.grid_layout.addWidget(self.plot_raw_check, 2, 1)
        self.grid_layout.addWidget(self.plot_log_label, 3, 0)
        self.grid_layout.addWidget(self.plot_log_check, 3, 1)
        self.grid_layout.addWidget(self.dq_label, 4, 0)
        self.grid_layout.addWidget(self.dq_input, 4, 1)
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

        self.scale_step.setMaximumWidth(80)
        self.scale_step.setValidator(QDoubleValidator(2.225e-308,np.inf,-1))
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

        self.create_widgets()
        self.style_widgets()
        self.create_layout()
        self.create_signals()

    def create_widgets(self):
        self.load_conv_data_btn = QPushButton('2θ Data')
        self.data_filename_lbl = QLabel('None')
        self.lambda_lbl = QLabel('Wavelength (λ): ')
        self.lambda_input = QLineEdit('')
        self.lambda_input.setValidator(QDoubleValidator(2.225e-308,np.inf,-1))
        self.save_conv_data_btn = QPushButton('Convert')

    def style_widgets(self):
        self.lambda_input.setMaximumWidth(120)
        self.data_filename_lbl.setAlignment(Qt.AlignCenter)
        self.lambda_lbl.setAlignment(Qt.AlignRight)
        self.lambda_input.setAlignment(Qt.AlignLeft)

    def create_layout(self):
        self.outer_layout = QVBoxLayout()
        self.outer_layout.setContentsMargins(25, 25, 25, 25)
        self.outer_layout.setSpacing(10)
        self.inner_layout = QHBoxLayout()
        self.inner_layout.setContentsMargins(0, 0, 0, 0)
        self.inner_layout.setSpacing(10)
        self.outer_layout.addWidget(self.load_conv_data_btn)
        self.outer_layout.addWidget(self.data_filename_lbl)
        self.inner_layout.addStretch()
        self.inner_layout.addWidget(self.lambda_lbl)
        self.inner_layout.addWidget(self.lambda_input)
        self.inner_layout.addStretch()
        self.outer_layout.addLayout(self.inner_layout)
        self.outer_layout.addWidget(self.save_conv_data_btn)
        self.setLayout(self.outer_layout)

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
            self.error_msg.setWindowTitle(__appname__ + ' v' + __version__)
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
            __lambda = np.float64(self.lambda_input.text())
        # Return if no wavelength set
        except ValueError:
            _message = ['No Wavelength Set!', 'Please set wavelength value']
            self.error_msg = QMessageBox()
            self.error_msg.setIcon(QMessageBox.Warning)
            self.error_msg.setStandardButtons(QMessageBox.Ok)
            self.error_msg.setText(_message[0])
            self.error_msg.setInformativeText((_message[1]))
            self.error_msg.setWindowTitle(__appname__ + ' v' + __version__)
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
        __q_data = data_utils.convert_two_theta(self.two_theta_data[:,0], __lambda)
        __out_data = np.column_stack((__q_data, self.two_theta_data[:,1]))
        np.savetxt(self.convert_filename, __out_data)
        # Message user success
        self.success_msg = QMessageBox()
        self.success_msg.setIcon(QMessageBox.Information)
        self.success_msg.setStandardButtons(QMessageBox.Ok)
        self.success_msg.setText('Data converted to Q-Space!')
        self.success_msg.setInformativeText(('Lambda: ' + str(__lambda) + '\nConverted data: ' + self.convert_filename))
        self.success_msg.setWindowTitle(__appname__ + ' v' + __version__)
        self.success_msg.adjustSize()
        self.success_msg.show()
        # Clear variables
        del __lambda, __q_data, __out_data, self.convert_filename


class DataCorrectionsGroupBox(QGroupBox):

    def __init__(self, *args):
        super(DataCorrectionsGroupBox, self).__init__(*args)
        self.setTitle('Data Corrections')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        self.setCheckable(False)

        self.create_widgets()
        self.style_widgets()
        self.create_layout()
        self.create_signals()

    def create_widgets(self):
        self.zero_shift_check = QCheckBox('Zero offset')

        self.zero_shift_button_group = QButtonGroup()
        self.zero_first_val_btn = QRadioButton('Zero I(Q=0)')
        self.zero_min_val_btn = QRadioButton('Zero min I(Q)')
        self.custom_shift_btn = QRadioButton('Custom Shift: ')
        self.zero_shift_button_group.addButton(self.zero_first_val_btn)
        self.zero_shift_button_group.addButton(self.zero_min_val_btn)
        self.zero_shift_button_group.addButton(self.custom_shift_btn)
        self.shift_sb = QDoubleSpinBox()
        self.shift_sb_step = QLineEdit('1')

        self.attenuation_slab_check = QCheckBox('Attenuation/Self-shielding (slab geometry)')
        self.wavelength_label = QLabel('Wavelength, λ (Å): ')
        self.wavelength_input = QLineEdit()
        self.attenuation_factor_label = QLabel('Attenuation factor, μ (cm\u207b\u00b9): ')
        self.attenuation_factor_input = QLineEdit()
        self.thickness_label = QLabel('Sample thickness (mm): ')
        self.thickness_input = QLineEdit()
        self.rotation_label = QLabel('Rotation (deg): ')
        self.rotation_input = QDoubleSpinBox()
        self.inspect_btn = QPushButton('Inspect')
        self.attenuation_slab_widgets = [self.wavelength_label, self.wavelength_input,
                                         self.attenuation_factor_label, self.attenuation_factor_input,
                                         self.thickness_label, self.thickness_input,
                                         self.rotation_label, self.rotation_input,
                                         self.inspect_btn]

    def style_widgets(self):
        self.zero_first_val_btn.setEnabled(False)
        self.zero_min_val_btn.setEnabled(False)
        self.custom_shift_btn.setEnabled(False)
        self.zero_first_val_btn.setChecked(True)
        self.attenuation_slab_check.setChecked(False)
        for _widget in self.attenuation_slab_widgets:
            _widget.setEnabled(False)

        self.shift_sb.setValue(0.0)
        self.shift_sb.setSingleStep(1)
        self.shift_sb.setDecimals(3)
        self.shift_sb.setMaximumWidth(200)
        self.shift_sb.setAlignment(Qt.AlignRight)
        self.shift_sb.setMinimum(-np.inf)
        self.shift_sb.setMaximum(np.inf)
        self.shift_sb.setEnabled(False)
        # Adaptive spin box steps
        #self.shift_sb.setStepType(QAbstractSpinBox.AdaptiveDecimalStepType)

        self.shift_sb_step.setMaximumWidth(100)
        self.shift_sb_step.setValidator(QDoubleValidator(2.225e-308,np.inf,-1))
        self.shift_sb_step.setAlignment(Qt.AlignRight)
        self.shift_sb_step.setEnabled(False)

        self.wavelength_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.attenuation_factor_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.thickness_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.rotation_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)

        self.wavelength_input.setValidator(QDoubleValidator(2.225e-308,np.inf,-1))
        self.attenuation_factor_input.setValidator(QDoubleValidator(2.225e-308,np.inf,-1))
        self.thickness_input.setValidator(QDoubleValidator(2.225e-308,np.inf,-1))
        self.thickness_input.setText('1.0')
        self.rotation_input.setValue(0.0)
        self.rotation_input.setDecimals(2)
        self.rotation_input.setMinimum(0.0)
        self.rotation_input.setMaximum(90.0)

        self.wavelength_input.setMaximumWidth(100)
        self.attenuation_factor_input.setMaximumWidth(100)
        self.thickness_input.setMaximumWidth(100)
        self.rotation_input.setMaximumWidth(100)
        self.inspect_btn.setMaximumWidth(150)

        self.wavelength_label.setToolTip('X-ray wavelength in Å (to calculate 2θ)')
        self.wavelength_input.setToolTip('X-ray wavelength in Å (to calculate 2θ)')
        self.attenuation_factor_label.setToolTip('Attenuation factor in inverse cm\nE.g. see NIST Standard Reference Database 66 or 126')
        self.attenuation_factor_input.setToolTip('Attenuation factor in inverse cm\nE.g. see NIST Standard Reference Database 66 or 126')
        self.thickness_label.setToolTip('Sample thickness in mm')
        self.thickness_input.setToolTip('Sample thickness in mm')
        self.rotation_label.setToolTip('Rotation from normal orientation in degrees (slab normal to the beam = 0°)')
        self.rotation_input.setToolTip('Rotation from normal orientation in degrees (slab normal to the beam = 0°)')
        self.inspect_btn.setToolTip('Inspect attenuation factor')

        self.spacer = QWidget()
        self.spacer.setFixedHeight(5)

    def create_layout(self):

        self.outer_layout = QVBoxLayout()
        self.outer_layout.setContentsMargins(25, 10, 25, 7)
        self.outer_layout.setSpacing(5)

        self.zero_grid_layout = QGridLayout()
        self.zero_grid_layout.setContentsMargins(25, 5, 10, 5)
        self.zero_grid_layout.setSpacing(5)
        self.zero_grid_layout.addWidget(self.zero_first_val_btn, 0, 0)
        self.zero_grid_layout.addWidget(self.zero_min_val_btn, 1, 0)
        self.zero_grid_layout.addWidget(self.custom_shift_btn, 2, 0)
        self.sb_layout = QHBoxLayout()
        self.sb_layout.addWidget(self.shift_sb)
        self.sb_layout.addWidget(self.shift_sb_step)
        self.zero_grid_layout.addLayout(self.sb_layout, 2, 1)
        self.zero_grid_layout.addWidget(QWidget(), 2, 2)

        self.slab_grid_layout = QGridLayout()
        self.slab_grid_layout.setContentsMargins(25, 5, 10, 5)
        self.slab_grid_layout.setSpacing(5)
        self.slab_grid_layout.addWidget(self.wavelength_label, 0, 0)
        self.slab_grid_layout.addWidget(self.wavelength_input, 0, 1)
        self.slab_grid_layout.addWidget(self.attenuation_factor_label, 1, 0)
        self.slab_grid_layout.addWidget(self.attenuation_factor_input, 1, 1)
        self.slab_grid_layout.addWidget(self.thickness_label, 2, 0)
        self.slab_grid_layout.addWidget(self.thickness_input, 2, 1)
        self.slab_grid_layout.addWidget(self.rotation_label, 3, 0)
        self.slab_grid_layout.addWidget(self.rotation_input, 3, 1)
        self.slab_grid_layout.addWidget(self.spacer, 4, 0)
        self.slab_grid_layout.addWidget(self.inspect_btn, 5, 0)
        self.slab_grid_layout.addWidget(QWidget(), 5, 2)

        self.outer_layout.addWidget(self.zero_shift_check)
        self.outer_layout.addLayout(self.zero_grid_layout)
        self.outer_layout.addWidget(self.attenuation_slab_check)
        self.outer_layout.addLayout(self.slab_grid_layout)
        self.setLayout(self.outer_layout)

    def create_signals(self):
        self.zero_shift_check.toggled.connect(self.toggle_zero_correction)
        self.shift_sb_step.editingFinished.connect(self.shift_step_changed)
        self.custom_shift_btn.toggled.connect(self.toggle_custom_shift)
        self.attenuation_slab_check.toggled.connect(self.toggle_attenuation_slab_correction)

    def shift_step_changed(self):
        self.shift_sb.setSingleStep(float(str(self.shift_sb_step.text())))

    def toggle_custom_shift(self, btn_state):
        self.shift_sb.setEnabled(btn_state)
        self.shift_sb_step.setEnabled(btn_state)

    def toggle_zero_correction(self, check_state):
        for btn in self.zero_shift_button_group.buttons():
            btn.setEnabled(check_state)
        if self.custom_shift_btn.isChecked():
            self.shift_sb.setEnabled(check_state)
            self.shift_sb_step.setEnabled(check_state)

    def toggle_attenuation_slab_correction(self, check_state):
        for _widget in self.attenuation_slab_widgets:
            _widget.setEnabled(check_state)
