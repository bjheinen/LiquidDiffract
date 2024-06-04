# -*- coding: utf-8 -*-
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2018-2024, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

try:
    import importlib_resources as resources
except ImportError:
    from importlib import resources
import time
import numpy as np
from qtpy.QtCore import Qt, Signal, QEventLoop
from qtpy.QtGui import QIntValidator, QDoubleValidator, QPixmap, \
                        QTextBlockFormat, QTextCursor, QIcon
from qtpy.QtWidgets import QDialog, QFileDialog, QMessageBox, \
                            QStyledItemDelegate, QFrame, QGroupBox, \
                            QScrollArea, QSplitter, QSizePolicy, \
                            QVBoxLayout, QHBoxLayout, QGridLayout, \
                            QLabel, QLineEdit, QCheckBox, QComboBox, \
                            QSpinBox, QDialogButtonBox, QPushButton, \
                            QListWidget, QButtonGroup, QRadioButton, \
                            QTextBrowser, QPlainTextEdit, QProgressBar, \
                            QWidget, QApplication
from pyqtgraph.exporters import ImageExporter
import LiquidDiffract.core.core as core
from LiquidDiffract.gui import plot_widgets
from LiquidDiffract.version import __appname__, __version__


def get_filename(io='open', caption='Load Data File', directory=None):
    if io == 'open':
        _filter = 'All Files (*);;Chi Files (*.chi);;Data Files (*.dat);;xy Files (*.xy)'
        file_name = QFileDialog.getOpenFileName(caption=caption,
                                                directory=directory,
                                                filter=_filter)
        file_name = file_name[0]
    elif io == 'save':
        _filter = 'All Files (*);;Chi Files (*.chi);;Data Files (*.dat);;xy Files (*.xy)'
        file_name = QFileDialog.getSaveFileName(caption=caption,
                                                directory=directory,
                                                filter=_filter)
        file_name = file_name[0]
        # Only check/add extensiton if filename returned
        if file_name:
            if not file_name.lower().endswith(('.dat', '.chi', '.xy')):
                file_name += '.dat'
    elif io == 'save_fig':
        _filter = 'Image Files (*.png *.jpg *.jpeg *.bmp *.tif *.tiff);;All Files (*)'
        file_name = QFileDialog.getSaveFileName(caption=caption,
                                                directory=directory,
                                                filter=_filter)
        file_name = file_name[0]
        # Check/add extension if fname returned
        if file_name:
            if not file_name.lower().endswith(('.png', '.jpg', '.jpeg', '.bmp', '.tif', '.tiff')):
                file_name += '.png'
    else:
        raise ValueError('Bad argument')

    if file_name:
        return file_name
    else:
        return None


class ErrorMessageBox(QMessageBox):
    def __init__(self, _message):
        super(ErrorMessageBox, self).__init__()
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setIcon(QMessageBox.Warning)
        self.setStandardButtons(QMessageBox.Ok)
        self.setText(_message[0])
        self.setInformativeText((_message[1]))
        self.setWindowTitle(__appname__ + ' v' + __version__)
        with resources.as_file(resources.files('LiquidDiffract.resources.icons').joinpath('gs_icon.png')) as path:
            self.setWindowIcon(QIcon(str(path)))
        self.adjustSize()


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


class CustomIntValidator(QIntValidator):
    '''Subclass of QIntValidator to fixup string on editingFinished'''
    def __init__(self, *args):
        super(CustomIntValidator, self).__init__(*args)

    def fixup(self, text):
        try:
            val = int(text)
        # Catch empty strings and do nothing
        except ValueError:
            return
        # Set to min if < min
        if val < self.bottom():
            return str(self.bottom())
        # Set to max if > max
        if val > self.top():
            return str(self.top())


class CustomDoubleValidator(QDoubleValidator):
    '''Subclass of QDoubleValidator to fixup string on editingFinished'''
    def __init__(self, *args):
        super(CustomDoubleValidator, self).__init__(*args)

    def fixup(self, text):
        try:
            val = float(text)
        # Catch empty strings and do nothing
        except ValueError:
            return
        # Set to min if < min
        if val < self.bottom():
            return str(self.bottom())
        # Set to max if > max
        if val > self.top():
            return str(self.top())


class PreferencesDialog(QDialog):

    fft_check_signal = Signal()
    fft_check_result = 0

    def __init__(self, preferences):
        super(PreferencesDialog, self).__init__()
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.title = 'Additional Preferences | ' + __appname__ + ' v' + __version__
        self.setWindowTitle(self.title)

        self.outer_layout = QVBoxLayout()
        self.outer_layout.setContentsMargins(5, 3, 5, 7)
        self.outer_layout.setSpacing(10)

        self.pref_widget = QWidget()
        self.vlayout = QVBoxLayout()
        self.vlayout.setContentsMargins(5, 3, 5, 7)
        self.vlayout.setSpacing(10)

        self.app_settings_gb = AppSettingsGroupBox(preferences)
        self.data_settings_gb = DataSettingsGroupBox(preferences)
        self.ref_proc_settings_gb = IterativeProcedureSettingsGroupBox(preferences)
        self.refine_settings_gb = SolverSettingsGroupBox(preferences)
        self.global_min_settings_gb = GlobalMinSettingsGroupBox(preferences)
        self.gaussian_fit_settings_gb = GaussianFittingSettingsGroupBox(preferences)

        self.vlayout.addWidget(self.app_settings_gb)
        self.vlayout.addWidget(self.data_settings_gb)
        self.vlayout.addWidget(self.ref_proc_settings_gb)
        self.vlayout.addWidget(self.refine_settings_gb)
        self.vlayout.addWidget(self.global_min_settings_gb)
        self.vlayout.addWidget(self.gaussian_fit_settings_gb)
        self.pref_widget.setLayout(self.vlayout)

        self.scroll_area = QScrollArea()
        self.scroll_area.setFrameShape(QFrame.NoFrame)
        self.scroll_area.setWidget(self.pref_widget)
        self.scroll_area.setWidgetResizable(True)

        self.button_box = QDialogButtonBox()
        self.button_box.addButton('&Cancel', QDialogButtonBox.RejectRole)
        self.button_box.addButton('&Apply', QDialogButtonBox.AcceptRole)

        self.hline = QFrame()
        self.hline.setFrameShape(QFrame.HLine)
        self.hline.setFrameShadow(QFrame.Sunken)
        self.hline.setObjectName('hline')

        self.outer_layout.addWidget(self.scroll_area)
        self.outer_layout.addWidget(self.hline)
        self.outer_layout.addWidget(self.button_box)

        self.setLayout(self.outer_layout)

        self.button_box.accepted.connect(self.accept_preferences)
        self.button_box.rejected.connect(self.rejected)

        self.rejected.connect(self.close)

    def accept_preferences(self):

        # Log mode
        if self.app_settings_gb.log_mode_input.currentText() == 'Append':
            _append_log_mode = 1
        else:
            _append_log_mode = 0

        try:
            _data_units = int(self.data_settings_gb.data_units_input.currentIndex())
            self.data_units_check = _data_units
            _window_length = int(self.data_settings_gb.window_length_input.text())
            _poly_order = int(self.data_settings_gb.poly_order_input.text())
            _rescale_AL = int(self.data_settings_gb.rescale_AL_input.isChecked())
            _fft_N = int(self.data_settings_gb.fft_N_input.value())
            self.fft_check = _fft_N
            self.fft_check_signal.emit()
            if self.fft_check_result == 1:
                raise RuntimeWarning()
            _mod_func_mode = int(self.ref_proc_settings_gb.mod_func_mode_input.isChecked())
            _op_method = self.refine_settings_gb.op_method_input.currentText()
            _disp = int(self.refine_settings_gb.disp_check.isChecked())
            _maxiter = int(self.refine_settings_gb.maxiter_input.text())
            _ftol = np.float64(self.refine_settings_gb.ftol_input.text())
            # Get l_bfgs_b specific options
            if _op_method == 'L-BFGS-B':
                _maxfun = int(self.refine_settings_gb.maxfun_input.text())
                _gtol = np.float64(self.refine_settings_gb.gtol_input.text())
            # Get globam minimisation (basin-hopping) options
            _bh_disp = int(self.global_min_settings_gb.disp_check.isChecked())
            _bh_niter = int(self.global_min_settings_gb.niter_basin_input.text())
            _bh_temp = float(self.global_min_settings_gb.temp_basin_input.text())
            _bh_step_size = float(self.global_min_settings_gb.stepsize_basin_input.text())
            _bh_interval = int(self.global_min_settings_gb.interval_basin_input.text())
            # Get Gaussiant fitting options
            _xray_weight_mode = int(self.gaussian_fit_settings_gb.xray_weight_mode_input.isChecked())

        # Handle for missing values
        except ValueError:
            _message = ['Missing Values!', 'Please ensures all values are set properly']
            _error_msg = ErrorMessageBox(_message)
            _error_msg.exec()
            return

        except RuntimeWarning:
            _message = ['\'N\' out of range!', 'Please increase N\nDefault: 12']
            _error_msg = ErrorMessageBox(_message)
            _error_msg.exec()
            return

        if not _window_length % 2:
            _message = ['Error!', 'Please ensure window_length is positive odd integer!']
            _error_msg = ErrorMessageBox(_message)
            _error_msg.exec()
            return

        if _poly_order >= _window_length:
            _message = ['Error!', 'Please ensure poly_order < window_length!']
            _error_msg = ErrorMessageBox(_message)
            _error_msg.exec()
            return

        # Set minimisation options dictionary
        _minimisation_options = {'disp': _disp,
                                 'maxiter': _maxiter,
                                 'ftol': _ftol
                                 }
        # Add additional optiosn for l_bfgs_b method
        if _op_method == 'L-BFGS-B':
            _minimisation_options['maxfun'] = _maxfun
            _minimisation_options['gtol'] = _gtol

        _global_minimisation = self.global_min_settings_gb.isChecked()
        _global_min_options = {'disp': _bh_disp,
                               'niter': _bh_niter,
                               'T': _bh_temp,
                               'stepsize': _bh_step_size,
                               'interval': _bh_interval
                               }

        # Set preferences dictionary to return
        self._preferences = {'append_log_mode': _append_log_mode,
                             'data_units': _data_units,
                             'window_length': _window_length,
                             'poly_order': _poly_order,
                             'rescale_AL': _rescale_AL,
                             'fft_N': _fft_N,
                             'mod_func_mode': _mod_func_mode,
                             'op_method': _op_method,
                             'minimisation_options': _minimisation_options,
                             'global_minimisation': _global_minimisation,
                             'global_min_options': _global_min_options,
                             'xray_weight_mode': _xray_weight_mode
                             }

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
        self.log_mode_input.setMaximumWidth(150)

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
        self.data_units_label = QLabel('Data file Q units: ')
        self.data_units_input = QComboBox()
        self.data_units_input.insertItem(0, '1/Angstroms')
        self.data_units_input.insertItem(1, '1/nano-metres')
        self.smoothing_label = QLabel('Savitsky-golay filter parameters (Data smoothing): ')
        self.window_length_label = QLabel('Window size')
        self.window_length_input = QLineEdit()
        self.poly_order_label = QLabel('Poly order')
        self.poly_order_input = QLineEdit()
        self.rescale_AL_label = QLabel('Plot rescaled Ashcroft-Langreth S(Q) & g(r)')
        self.rescale_AL_input = QCheckBox()
        self.fft_label = QLabel('FFT Options:')
        self.fft_N_label = QLabel('<b><i>N</i></b> | Size of padded array for FFT = 2<sup>N</sup>: ')
        self.fft_N_input = QSpinBox()

    def set_data(self, preferences):
        self.data_units_input.setCurrentIndex(preferences['data_units'])
        self.window_length_input.setText(str(preferences['window_length']))
        self.poly_order_input.setText(str(preferences['poly_order']))
        self.rescale_AL_input.setChecked(preferences['rescale_AL'])
        self.fft_N_input.setValue(preferences['fft_N'])

    def style_widgets(self):
        self.data_units_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.data_units_input.setMaximumWidth(300)

        self.smoothing_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        self.window_length_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.window_length_input.setAlignment(Qt.AlignRight)
        self.window_length_input.setValidator(QIntValidator())
        self.window_length_input.setMaximumWidth(70)

        self.poly_order_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.poly_order_input.setAlignment(Qt.AlignRight)
        self.poly_order_input.setValidator(QIntValidator())
        self.poly_order_input.setMaximumWidth(70)

        self.rescale_AL_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        self.fft_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        
        self.fft_N_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.fft_N_input.setAlignment(Qt.AlignRight)
        self.fft_N_input.setRange(1,20)
        self.fft_N_input.setMaximumWidth(70)

    def create_layout(self):
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(20, 10, 20, 7)
        self.main_layout.setSpacing(25)

        self.grid_layout = QGridLayout()
        self.grid_layout.setSpacing(15)

        self.grid_layout.addWidget(self.data_units_label, 0, 0)
        self.grid_layout.addWidget(self.data_units_input, 0, 1)
        self.grid_layout.addWidget(self.smoothing_label, 1, 0)
        self.grid_layout.addWidget(self.window_length_label, 2, 0)
        self.grid_layout.addWidget(self.window_length_input, 2, 1)
        self.grid_layout.addWidget(self.poly_order_label, 3, 0)
        self.grid_layout.addWidget(self.poly_order_input, 3, 1)
        self.grid_layout.addWidget(self.rescale_AL_label, 4, 0)
        self.grid_layout.addWidget(self.rescale_AL_input, 4, 1)
        self.grid_layout.addWidget(self.fft_label, 5, 0)
        self.grid_layout.addWidget(self.fft_N_label, 6, 0)
        self.grid_layout.addWidget(self.fft_N_input, 6, 1)

        self.main_layout.addLayout(self.grid_layout)
        self.setLayout(self.main_layout)


class IterativeProcedureSettingsGroupBox(QGroupBox):

    def __init__(self, preferences):
        super(IterativeProcedureSettingsGroupBox, self).__init__()
        self.setTitle('Iterative Refinement Options')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')

        self.create_widgets()
        self.set_data(preferences)
        self.style_widgets()
        self.create_layout()

    def create_widgets(self):
        self.mod_func_mode_label = QLabel('Use modification func in iterative refinement?: ')
        self.mod_func_mode_input = QCheckBox()

    def set_data(self, preferences):
        self.mod_func_mode_input.setChecked(preferences['mod_func_mode'])

    def style_widgets(self):
        self.mod_func_mode_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

    def create_layout(self):
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(20, 10, 20, 7)
        self.main_layout.setSpacing(25)

        self.grid_layout = QGridLayout()
        self.grid_layout.setSpacing(15)
        self.grid_layout.addWidget(self.mod_func_mode_label, 0, 0)
        self.grid_layout.addWidget(self.mod_func_mode_input, 0, 1)

        self.main_layout.addLayout(self.grid_layout)
        self.setLayout(self.main_layout)


class SolverSettingsGroupBox(QGroupBox):

    def __init__(self, preferences):
        super(SolverSettingsGroupBox, self).__init__()
        self.setTitle('Solver Options')
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
        # Default is to use limited-memory BFGS code for optimising rho
        # See http://users.iems.northwestern.edu/~nocedal/lbfgsb.html for details
        self.op_method_input.insertItem(0, 'L-BFGS-B')
        self.op_method_input.insertItem(1, 'SLSQP')
        self.op_method_input.insertItem(2, 'COBYLA')

        _app_url = '"https://github.com/bjheinen/LiquidDiffract#density-ρ-refinement"'
        _text = (f'<a href={_app_url}><span style="color: #0c0263;"><span '
                 f'lang="zxx"><u>More information...</u></span></span></a>')

        self.solver_info_link = QLabel(_text)
        self.solver_info_link.setOpenExternalLinks(True)

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

        self.hline = QFrame()
        self.hline.setFrameShape(QFrame.HLine)
        self.hline.setFrameShadow(QFrame.Sunken)
        self.hline.setObjectName('hline')

    def set_data(self, preferences):
        self.op_method_input.setCurrentText(preferences['op_method'])
        _min_options = preferences['minimisation_options']
        self.disp_check.setChecked(_min_options['disp'])
        self.maxiter_input.setText(str(_min_options['maxiter']))
        self.ftol_input.setText(str(_min_options['ftol']))

        # Handle missing maxfun/gtol if op_method!='L-BFGS-B'
        try:
            self.maxfun_input.setText(str(_min_options['maxfun']))
            self.gtol_input.setText(str(_min_options['gtol']))
        except KeyError:
            pass

    def style_widgets(self):

        self.op_method_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.op_method_label.setToolTip('Solver to use when refining density.')
        self.op_method_input.setMaximumWidth(150)

        self.solver_info_link.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        self.disp_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.disp_label.setToolTip('Print solver specific convergence messages?')

        self.maxiter_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.maxiter_input.setAlignment(Qt.AlignRight)
        self.maxiter_input.setValidator(QIntValidator())
        self.maxiter_input.setMaximumWidth(100)

        self.maxfun_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.maxfun_input.setAlignment(Qt.AlignRight)
        self.maxfun_input.setValidator(QIntValidator())
        self.maxfun_input.setMaximumWidth(100)

        self.ftol_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        _tooltip = ('The function convergence limit is the precision goal for '
                    'the function value in the stopping criterion. \nFor the '
                    'L-BFGS-B method the iteration stops when '
                    '(f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= ftol). \n'
                    'For the COBYLA method ftol is a lower bound on the trust '
                    'region and is not precisely guaranteed.'
                    )
        self.ftol_label.setToolTip(_tooltip)
        self.ftol_input.setAlignment(Qt.AlignRight)
        self.ftol_input.setValidator(QDoubleValidator())
        self.ftol_input.setMaximumWidth(100)

        self.gtol_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        _tooltip = ('Gradient limit for stopping criterion. \nThe iteration '
                    'will stop when max{|proj g_i | i = 1, ..., n} <= gtol \n'
                    'where pg_i is the i-th component of '
                    'the projected gradient.'
                    )
        self.gtol_label.setToolTip(_tooltip)
        self.gtol_input.setAlignment(Qt.AlignRight)
        self.gtol_input.setValidator(QDoubleValidator())
        self.gtol_input.setMaximumWidth(100)

        if self.op_method_input.currentText() != 'L-BFGS-B':
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
        self.grid_layout.addWidget(self.solver_info_link, 1, 0)
        self.grid_layout.addWidget(self.hline, 2, 0)
        self.grid_layout.addWidget(self.disp_label, 3, 0)
        self.grid_layout.addWidget(self.disp_check, 3, 1)
        self.grid_layout.addWidget(self.maxiter_label, 4, 0)
        self.grid_layout.addWidget(self.maxiter_input, 4, 1)
        self.grid_layout.addWidget(self.maxfun_label, 5, 0)
        self.grid_layout.addWidget(self.maxfun_input, 5, 1)
        self.grid_layout.addWidget(self.ftol_label, 6, 0)
        self.grid_layout.addWidget(self.ftol_input, 6, 1)
        self.grid_layout.addWidget(self.gtol_label, 7, 0)
        self.grid_layout.addWidget(self.gtol_input, 7, 1)

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

        elif self.op_method_input.currentText() == 'COBYLA':
            self.maxfun_label.setEnabled(False)
            self.maxfun_input.setEnabled(False)
            self.gtol_label.setEnabled(False)
            self.gtol_input.setEnabled(False)
            # Set default values when the method is changed
            self.disp_check.setChecked(False)
            self.maxiter_input.setText('1000')
            self.ftol_input.setText('1e-7')

        else:
            pass


class GlobalMinSettingsGroupBox(QGroupBox):

    def __init__(self, preferences):
        super(GlobalMinSettingsGroupBox, self).__init__()

        self.setTitle('Global Minimisation?')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        self.setCheckable(True)

        self.create_widgets()
        self.set_data(preferences)
        self.style_widgets()
        self.create_layout()

    def create_widgets(self):

        self.disp_label = QLabel('Output convergence information? ')
        self.disp_check = QCheckBox()

        self.niter_basin_label = QLabel('Number of basin-hopping iterations: ')
        self.niter_basin_input = QLineEdit()

        self.temp_basin_label = QLabel('Basin-hopping temperature parameter (T): ')
        self.temp_basin_input = QLineEdit()

        self.stepsize_basin_label = QLabel('Initial step size for basin-hopping algorithm: ')
        self.stepsize_basin_input = QLineEdit()

        self.interval_basin_label = QLabel('Update interval: ')
        self.interval_basin_input = QLineEdit()

    def set_data(self, preferences):

        self.setChecked(preferences['global_minimisation'])
        _global_min_options = preferences['global_min_options']
        self.disp_check.setChecked(_global_min_options['disp'])

        self.niter_basin_input.setText(str(_global_min_options['niter']))
        self.temp_basin_input.setText(str(_global_min_options['T']))
        self.stepsize_basin_input.setText(str(_global_min_options['stepsize']))
        self.interval_basin_input.setText(str(_global_min_options['interval']))

    def style_widgets(self):

        self.disp_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.disp_label.setToolTip('Print basin-hopping status messages?')

        self.niter_basin_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.niter_basin_input.setAlignment(Qt.AlignRight)
        self.niter_basin_input.setValidator(QIntValidator())
        self.niter_basin_input.setMaximumWidth(80)

        self.temp_basin_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.temp_basin_input.setAlignment(Qt.AlignRight)
        self.temp_basin_input.setValidator(QDoubleValidator())
        self.temp_basin_input.setMaximumWidth(80)

        self.stepsize_basin_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.stepsize_basin_input.setAlignment(Qt.AlignRight)
        self.stepsize_basin_input.setValidator(QDoubleValidator())
        self.stepsize_basin_input.setMaximumWidth(80)

        self.interval_basin_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.interval_basin_input.setAlignment(Qt.AlignRight)
        self.interval_basin_input.setValidator(QIntValidator())
        self.interval_basin_input.setMaximumWidth(80)

    def create_layout(self):
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(20, 10, 20, 7)
        self.main_layout.setSpacing(25)

        self.grid_layout = QGridLayout()
        self.grid_layout.setSpacing(15)
        self.grid_layout.addWidget(self.disp_label, 0, 0)
        self.grid_layout.addWidget(self.disp_check, 0, 1)
        self.grid_layout.addWidget(self.niter_basin_label, 1, 0)
        self.grid_layout.addWidget(self.niter_basin_input, 1, 1)
        self.grid_layout.addWidget(self.temp_basin_label, 2, 0)
        self.grid_layout.addWidget(self.temp_basin_input, 2, 1)
        self.grid_layout.addWidget(self.stepsize_basin_label, 3, 0)
        self.grid_layout.addWidget(self.stepsize_basin_input, 3, 1)
        self.grid_layout.addWidget(self.interval_basin_label, 4, 0)
        self.grid_layout.addWidget(self.interval_basin_input, 4, 1)

        self.main_layout.addLayout(self.grid_layout)

        self.setLayout(self.main_layout)


class GaussianFittingSettingsGroupBox(QGroupBox):

    def __init__(self, preferences):
        super(GaussianFittingSettingsGroupBox, self).__init__()
        self.setTitle('Gaussian Fitting Options')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')

        self.create_widgets()
        self.set_data(preferences)
        self.style_widgets()
        self.create_layout()

    def create_widgets(self):
        self.xray_weight_mode_label = QLabel('Use Q-dependent x-ray weighting factors?: ')
        self.xray_weight_mode_input = QCheckBox()

    def set_data(self, preferences):
        self.xray_weight_mode_input.setChecked(preferences['xray_weight_mode'])

    def style_widgets(self):
        self.xray_weight_mode_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        _tooltip = ('Default behaviour is to use the effective atomic number, Kp, \n'
                    'of each species given by the Warren-Krutter-Morningstar approximation at Q=0.\n'
                    'Check this box to calculate Kp as an average of the entire Q-range used.'
                    )
        self.xray_weight_mode_label.setToolTip(_tooltip)

    def create_layout(self):
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(20, 10, 20, 7)
        self.main_layout.setSpacing(25)

        self.grid_layout = QGridLayout()
        self.grid_layout.setSpacing(15)
        self.grid_layout.addWidget(self.xray_weight_mode_label, 0, 0)
        self.grid_layout.addWidget(self.xray_weight_mode_input, 0, 1)

        self.main_layout.addLayout(self.grid_layout)
        self.setLayout(self.main_layout)


class AboutDialog(QDialog):
    def __init__(self):
        super(AboutDialog, self).__init__()
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setWindowTitle
        self.title = __appname__ + ' v' + __version__
        self.setWindowTitle(self.title)

        self.vlayout = QVBoxLayout()
        self.vlayout.setContentsMargins(5, 15, 5, 7)
        self.vlayout.setSpacing(10)

        self.logo_box = QLabel()

        with resources.as_file(resources.files('LiquidDiffract.resources.icons').joinpath('logo.png')) as path:
            self.logo = QPixmap(str(path))
        self.logo_box.setPixmap(self.logo)
        self.logo_box.setAlignment(Qt.AlignCenter)

        self.text_display = QTextBrowser()
        self.text_display.setReadOnly(True)
        self.text_display.setFrameStyle(QFrame.NoFrame)
        self.text_display.setStyleSheet("* { background-color: rgba(0, 0, 0, 0); }")
        self.text_display.setOpenExternalLinks(True)

        _app_url = 'https://github.com/bjheinen/LiquidDiffract'
        _paper_url = 'https://doi.org/10.1007/s00269-022-01186-6'
        _gpl_url = 'https://www.gnu.org/licenses/gpl.html'
        _pstr = '<p class="western" align="center">'
        _description = ('LiquidDiffract is a graphical application for '
                        'X-ray total-scattering analysis of liquids and '
                        'disordered solids.\n'
                        'It implements procedures to obtain information on '
                        'macroscopic bulk properties and local atomic-scale '
                        'structure from total scattering data.\n'
                        'LiquidDiffract provides tools for  background '
                        'subtraction; calculation, normalisation, and '
                        'refinement of the reciprocal-space structure factor '
                        'and real-space correlation functions; '
                        'and extraction of structural information such as '
                        'bond lengths, coordination number and bulk density.'
                        )
        _paper_str = (f'Heinen, B. J., & Drewitt, J. W. (2022). '
                      f'LiquidDiffract: Software for liquid total scattering analysis. '
                      f'<em>Physics and Chemistry of Minerals, 49:9</em>. doi:10.1007/s00269-022-01186-6')

        _app_url_str = (f'{_pstr}<a class="western" href="{_app_url}">'
                        f'<span style="color: #000080;"><span lang="zxx">'
                        f'<u>{_app_url}</u></span></span></a></p>'
                        )
        _cite_str = (f'{_pstr}If you use LiquidDiffract in your work please cite it as:<br><br>'
                     f'<a class="western" href="{_paper_url}">'
                     f'<span style="color: #000080;"><span lang="zxx">'
                     f'<u>{_paper_str}</u></span></span></a></p>')

        _copyright = 'Copyright &copy; 2018-2024 &ndash; Benedict J. Heinen'
        _warranty = ('This program comes with absolutely no '
                     'warranty or guarantee.')
        _license = '<u>GNU General Public Licence, version 3 or later</u>'

        _l_str = (f'See the </span><a class="western" href="{_gpl_url}">'
                  f'<span style="color: #000080;"><span style='
                  f'"font-size: small;"><span lang="zxx">{_license}</span>'
                  f'</span></span></a><span style="font-size: small;"> '
                  f'for details.</span></p>'
                  )

        _text = (f'{_pstr}<strong>{__appname__}</strong></p>'
                 f'{_pstr}v{__version__}</p>'
                 f'{_pstr}&nbsp;</p>'
                 f'{_pstr}<em>{_description}</em></p>'
                 f'{_pstr}&nbsp;</p>'
                 f'{_app_url_str}'
                 f'{_pstr}&nbsp;</p>'
                 f'{_cite_str}'
                 f'{_pstr}&nbsp;</p>'
                 f'{_pstr}<span style="font-size: small;">{_copyright}</span></p>'
                 f'{_pstr}<span style="font-size: small;">{_warranty}<br>'
                 f'{_l_str}'
                 )
        self.text_display.textCursor().insertHtml(_text)

        self.vlayout.addWidget(self.logo_box)
        self.vlayout.addWidget(self.text_display)
        self.setLayout(self.vlayout)


class CheckFileDialog(QDialog):

    def __init__(self, data_file):
        super(CheckFileDialog, self).__init__()
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.title = 'Check Data File! | ' + __appname__ + ' v' + __version__
        self.setWindowTitle(self.title)
        with resources.as_file(resources.files('LiquidDiffract.resources.icons').joinpath('gs_icon.png')) as path:
            self.setWindowIcon(QIcon(str(path)))

        self.data_file = data_file

        self.create_widgets()
        self.create_layout()
        self.create_signals()
        self.init_data()

    def create_widgets(self):

        self.info_message = QLabel('<b>Error trying to read file! Does it contain a header?</b> \
                                   <br>Please select the header block below, use the hash (#) \
                                    symbol to mark header lines, or select a new file.')

        self.fname_lbl = QLabel('Filename: <i>' + self.data_file + '</i>')
        self.fname_lbl.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        self.header_len_lbl = QLabel('Header length (lines)')
        self.header_len = QSpinBox()
        self.header_len_lbl.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.header_len.setAlignment(Qt.AlignRight)
        self.header_len.setRange(0, 100)
        self.header_len.setMaximumWidth(70)

        # Text area for file preview
        self.preview_file_area = QPlainTextEdit()
        self.preview_file_area.setReadOnly(True)

        self.button_box = QDialogButtonBox()
        self.button_box.addButton('&Cancel', QDialogButtonBox.RejectRole)
        self.button_box.addButton('&Load File', QDialogButtonBox.AcceptRole)

        self.default_format = QTextBlockFormat(self.preview_file_area.document().firstBlock().blockFormat())
        self.highlight_format = QTextBlockFormat()
        self.highlight_format.setBackground(Qt.yellow)

    def create_layout(self):

        self.outer_layout = QVBoxLayout()
        self.outer_layout.setContentsMargins(5, 3, 5, 7)
        self.outer_layout.setSpacing(10)

        self.vlayout = QVBoxLayout()
        self.vlayout.setContentsMargins(5, 3, 5, 7)
        self.vlayout.setSpacing(10)

        self.grid_layout = QGridLayout()
        self.grid_layout.setSpacing(15)
        self.grid_layout.addWidget(self.fname_lbl, 0 , 0)
        self.grid_layout.addWidget(self.header_len_lbl, 1, 0)
        self.grid_layout.addWidget(self.header_len, 1, 1)

        self.vlayout.addWidget(self.info_message)
        self.vlayout.addLayout(self.grid_layout)
        self.vlayout.addWidget(self.preview_file_area)

        self.preview_widget = QWidget()
        self.preview_widget.setLayout(self.vlayout)

        self.scroll_area = QScrollArea()
        self.scroll_area.setFrameShape(QFrame.NoFrame)
        self.scroll_area.setWidget(self.preview_widget)
        self.scroll_area.setWidgetResizable(True)

        self.hline = QFrame()
        self.hline.setFrameShape(QFrame.HLine)
        self.hline.setFrameShadow(QFrame.Sunken)
        self.hline.setObjectName('hline')

        self.outer_layout.addWidget(self.scroll_area)
        self.outer_layout.addWidget(self.hline)
        self.outer_layout.addWidget(self.button_box)

        self.setLayout(self.outer_layout)

    def init_data(self):
        with open(self.data_file, "r") as f:
            file_contents = f.read()
            self.preview_file_area.setPlainText(file_contents)
        self.header_len.setValue(4)

    def create_signals(self):
        self.button_box.accepted.connect(self.attempt_file_load)
        self.button_box.rejected.connect(self.rejected)
        self.rejected.connect(self.close)
        self.header_len.valueChanged.connect(self.highlight_header)

    def highlight_header(self):
        # Remove existing highlighting from entire document
        _cursor = QTextCursor(self.preview_file_area.document())
        _cursor.select(QTextCursor.Document)
        _cursor.setBlockFormat(self.default_format)
        if self.header_len.value() > 0:
            # Highlight header lines (n-1)
            _cursor = QTextCursor(self.preview_file_area.document().findBlockByLineNumber(self.header_len.value() - 1))
            _cursor.movePosition(QTextCursor.Start, QTextCursor.KeepAnchor)
            _cursor.setBlockFormat(self.highlight_format)
        else:
            return

    def attempt_file_load(self):
        self._header_len = self.header_len.value()
        try:
            _x, _y = np.loadtxt(self.data_file, unpack=True, skiprows=self._header_len)
            self.done(1)
        except ValueError:
            _message = ['Error loading file!', 'Unable to load file.\nPlease check header.']
            ErrorMessageBox(_message).exec()
            return

    def get_header_len(self):
        return self._header_len


class ComputeMapDialog(QDialog):

    def __init__(self, data, ref_data, prefs, bkg_flag):
        super(ComputeMapDialog, self).__init__()
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.title = 'Compute χ\u00b2 map | ' + __appname__ + ' v' + __version__
        self.setWindowTitle(self.title)

        # data, extra refinement data, and app preferences
        self.data = data
        self.ref_data = ref_data
        self.prefs = prefs
        self.bkg_flag = bkg_flag
        # init map_data
        self.map_data = None

        self.layout = QHBoxLayout(self)
        self.layout.setSpacing(0)
        # Make Config Widget
        self.compute_map_config_widget = ComputeMapConfigWidget()
        self.config_scroll_area = QScrollArea()
        self.config_scroll_area.setFrameShape(QFrame.NoFrame)
        self.config_scroll_area.setWidget(self.compute_map_config_widget)
        self.config_scroll_area.setWidgetResizable(True)
        self.config_scroll_area.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.Ignored)

        # Make vertical line separator
        self.vline = QFrame()
        self.vline.setFrameShape(QFrame.VLine)
        self.vline.setFrameShadow(QFrame.Sunken)
        self.vline.setObjectName("vline")
        # Make Plot Widget
        self.map_plot_widget = plot_widgets.ChiSquaredMapWidget()
        self.plot_scroll_area = QScrollArea()
        self.plot_scroll_area.setWidget(self.map_plot_widget)
        self.plot_scroll_area.setWidgetResizable(True)
        self.plot_scroll_area.setFrameShape(QFrame.NoFrame)
        self.hsplitter = QSplitter(Qt.Horizontal)
        self.hsplitter.addWidget(self.config_scroll_area)
        self.hsplitter.addWidget(self.plot_scroll_area)
        self.hsplitter.setStretchFactor(0, 2)
        self.hsplitter.setStretchFactor(1, 5)
        self.layout.addWidget(self.hsplitter)

        self.setLayout(self.layout)

        # Disable background scaling option if no background present
        if not self.bkg_flag:
            _x_b_item = self.compute_map_config_widget.map_options_gb.x_param_input_list.item(1)
            _y_b_item = self.compute_map_config_widget.map_options_gb.y_param_input_list.item(1)
            # (item.flags() & ~Qt.ItemIsSelectable & ~Qt.ItemIsEnabled)
            # Qt.ItemIsSelectable == 1, Qt.ItemIsEnabled == 32
            _x_b_item.setFlags(_x_b_item.flags() & ~Qt.ItemIsSelectable & ~Qt.ItemIsEnabled)
            _y_b_item.setFlags(_y_b_item.flags() & ~Qt.ItemIsSelectable & ~Qt.ItemIsEnabled)

        self.create_signals()

    def create_signals(self):
        self.compute_map_config_widget.map_options_gb.estimate_time_btn.clicked.connect(self.set_time_estimate_text)
        self.compute_map_config_widget.map_options_gb.compute_map_btn.clicked.connect(self.compute_chi_squared_map)
        self.compute_map_config_widget.map_result_gb.map_options_changed.connect(self.plot_data)
        self.compute_map_config_widget.map_result_gb.export_map_btn.clicked.connect(self.export_figure)
        self.compute_map_config_widget.map_result_gb.save_data_btn.clicked.connect(self.save_data)
        self.map_plot_widget.img.mouse_pos_changed.connect(self.compute_map_config_widget.map_pos_gb.set_pos_label)
        self.map_plot_widget.img.mouse_pos_cleared.connect(self.compute_map_config_widget.map_pos_gb.clear_pos_label)

    def make_map_args(self):
        # Get map parameter names/types
        x_param, y_param = self.compute_map_config_widget.map_options_gb.get_param_names()
        # Get x/y parameter arrays to compute with validation
        try:
            x_array, y_array, x_min, y_min, x_step, y_step = self.compute_map_config_widget.map_options_gb.get_param_array()
        except RuntimeError:
            raise RuntimeError()
        # Make arg dict for core.refinement_objfun
        # If varying bkg_scaling additional data is required
        if 'bkg_scaling' in [x_param, y_param]:
            arg_dict = {'Q_data': self.data['cor_x'],
                        'I_data_uncorrected': self.data['uncorrected_y'],
                        'I_bkg': self.data['bkg_y'],
                        'data_correction': self.data['data_correction'],
                        'q_min': self.data['qmin'],
                        'q_max': self.data['qmax'],
                        'smooth_flag': self.ref_data['smooth_flag'],
                        'smooth_window_length': self.prefs['window_length'],
                        'smooth_poly_order': self.prefs['poly_order'],
                        'composition': self.ref_data['composition'],
                        'r_min': self.ref_data['r_min'],
                        'n_iter': self.ref_data['n_iter'],
                        'sq_method': self.ref_data['sq_method'],
                        'mod_func': self.ref_data['mod_func'],
                        'window_start': self.ref_data['window_start'],
                        'fft_N': self.prefs['fft_N'],
                        'opt_rho': 1, 'opt_bkg': 1}

        else:
            arg_dict = {'Q_data': self.data['cor_x_cut'],
                        'I_data': self.data['cor_y_cut'],
                        'composition': self.ref_data['composition'],
                        'r_min': self.ref_data['r_min'],
                        'n_iter': self.ref_data['n_iter'],
                        'sq_method': self.ref_data['sq_method'],
                        'mod_func': self.ref_data['mod_func'],
                        'window_start': self.ref_data['window_start'],
                        'fft_N': self.prefs['fft_N'],
                        'opt_rho': 1, 'opt_bkg': 0}

        # Make iterable array of compute parameters
        # Get cartesian product of x and y
        param_grid = np.array(np.meshgrid(x_array, y_array)).T.reshape(-1,2)
        n_pairs = len(param_grid)

        # Set parameter names
        self.x_param = x_param
        self.y_param = y_param
        # Set x, y arrays
        self.x_array = x_array
        self.y_array = y_array

        # Set map shape info
        self.transform_scale = (x_step, y_step)
        self.transform_translate = ((x_min - (x_step/2.0))/x_step, (y_min - (y_step/2.0))/y_step)
        self.map_shape = (len(x_array), len(y_array))

        # Get arrays for 4 possible parameters
        if 'rho' == x_param:
            rho_arr = param_grid[:,0]
        elif 'rho' == y_param:
            rho_arr = param_grid[:,1]
        else:
            rho_arr = np.repeat(self.ref_data['rho'], n_pairs)

        if 'bkg_scaling' == x_param:
            bkg_scaling_arr = param_grid[:,0]
        elif 'bkg_scaling' == y_param:
            bkg_scaling_arr = param_grid[:,1]
        else:
            bkg_scaling_arr = np.repeat(np.nan, n_pairs)

        if 'r_min' == x_param:
            r_min_arr = param_grid[:,0]
        elif 'r_min' == y_param:
            r_min_arr = param_grid[:,1]
        else:
            r_min_arr = np.repeat(arg_dict['r_min'], n_pairs)

        if 'n_iter' == x_param:
            n_iter_arr = param_grid[:,0]
        elif 'n_iter' == y_param:
            n_iter_arr = param_grid[:,1]
        else:
            n_iter_arr = np.repeat(arg_dict['n_iter'], n_pairs)

        # Make parameter array with cols: rho, bkg_scaling, r_min, n_iter
        f_params = np.column_stack((rho_arr, bkg_scaling_arr, r_min_arr, n_iter_arr))

        return f_params, arg_dict

    def estimate_map_time(self, f_params=None, arg_dict=None):
        """Estimate time for map based on time for average value (used for n_iter)"""
        if f_params is None or arg_dict is None:
            try:
                f_params, arg_dict = self.make_map_args()
            except RuntimeError:
                raise RuntimeError()
        mean_param = np.mean(f_params, axis=0)
        start_time = time.perf_counter()
        # Add some time for list setup etc.
        _est_list = []
        # Single function call
        _est_list.append(core.chi_squared_map_helper(mean_param, arg_dict))
        _est_list = np.asarray(_est_list)
        end_time = time.perf_counter()
        # Time in seconds for calculating average param
        estimated_time = end_time - start_time
        # Estimated time in seconds for total map
        estimated_time *= len(f_params)
        # Prefer overestimate of time
        estimated_time *= 1.2
        return estimated_time

    def set_time_estimate_text(self):
        # For high n_iter getting the estimate may be slow
        # Disable widgets while running
        self.compute_map_config_widget.map_options_gb.enable_widgets(False)
        # Reset progress bar
        self.compute_map_config_widget.progress_gb.reset_progress()
        # Clear Map Inspect labels
        self.compute_map_config_widget.map_pos_gb.clear_all()
        # Clear plot data
        self.clear_map_data()
        # Get time estimate
        try:
            estimated_time = self.estimate_map_time()
        except RuntimeError:
            # Re-enable widgets and return
            self.compute_map_config_widget.map_options_gb.enable_widgets(True)
            return
        # Format time string
        if estimated_time < 60:
            time_str = str(int(estimated_time)) + ' s'
        else:
            time_str = time.strftime('%H:%M:%S', time.gmtime(estimated_time))
        self.compute_map_config_widget.map_options_gb.estimate_time_label.setText(time_str)
        # Re-enable widgets
        self.compute_map_config_widget.map_options_gb.enable_widgets(True)

    def compute_chi_squared_map(self):
        '''Computes a map of χ^2 values for x, y values of two parameters.'''
        # Get parameter array and function arguments
        try:
            f_params, arg_dict = self.make_map_args()
        except RuntimeError:
            return

        # Get total time estimate
        estimated_time = self.estimate_map_time(f_params=f_params, arg_dict=arg_dict)
        # Get start time for map computation
        start_time = time.perf_counter()
        self.compute_map_config_widget.progress_gb.init_progress(len(f_params), estimated_time, start_time)

        chi_sq_map = []
        # Loop over f_params, calculate chi^2 for each value and update progress
        for iteration, parameter in enumerate(f_params):
            chi_sq = core.chi_squared_map_helper(parameter, arg_dict)
            chi_sq_map.append(chi_sq)
            self.compute_map_config_widget.progress_gb.update_progress(iteration+1)
        chi_sq_map = np.asarray(chi_sq_map)
        # Store 1d array
        self.raw_map_data = chi_sq_map
        # Reshape map
        self.map_data = chi_sq_map.reshape(self.map_shape).T

        # Plot data
        self.plot_data()

    def plot_data(self):
        # Return if no data
        if self.map_data is None:
            return
        # Apply data normalisations
        self.normalise_data()
        # Pack parameter names
        _param_names = (self.x_param, self.y_param)
        # Get contour level options
        norm_flag, log_flag, contour_flag, n_levels, _ = self.compute_map_config_widget.map_result_gb.get_map_options()
        # Plot map_data
        self.map_plot_widget.set_data(self.norm_map_data,
                                      self.transform_scale, self.transform_translate,
                                      _param_names,
                                      contour_flag=contour_flag, n_levels=n_levels)
        # Update Inspect Map labels
        self.compute_map_config_widget.map_pos_gb.set_param_label(*_param_names, norm_flag, log_flag)

    def normalise_data(self):
        # Get data normalisation options
        norm_flag, log_flag, *_ = self.compute_map_config_widget.map_result_gb.get_map_options()
        # Normalise data
        if norm_flag == 0:
            norm_map_data = self.map_data
        elif norm_flag == 1:
            norm_map_data = self.map_data / np.min(self.map_data)
        elif norm_flag == 2:
            norm_map_data = (self.map_data.T / np.min(self.map_data, axis=1)).T
        elif norm_flag == 3:
            norm_map_data = self.map_data / np.min(self.map_data, axis=0)
        else:
            raise NotImplementedError()
        # Apply data transforms (log etc.)
        if log_flag == 0:
            norm_map_data = norm_map_data
        elif log_flag == 1:
            norm_map_data = np.log(norm_map_data)
        elif log_flag == 2:
            norm_map_data = np.log(norm_map_data)
            log_zero_mask = np.where(norm_map_data == 0)
            norm_map_data[log_zero_mask] = 1.e-14
            norm_map_data = np.log(norm_map_data)
        else:
            raise NotImplementedError()
        # Set normalised data for plotting
        self.norm_map_data = norm_map_data

    def clear_map_data(self):
        # Clear plots
        self.map_plot_widget.clear_data()
        # Clear map data
        self.raw_map_data = None
        self.map_data = None
        self.norm_map_data = None

    def save_data(self):
        # Get option to save normalised map or raw data
        *_, save_map_flag = self.compute_map_config_widget.map_result_gb.get_map_options()

        if self.map_data is None:
            return
        # Get save filename
        #__default_file_name?
        _file_name = get_filename(io='save', caption='Save χ\u00b2 Map Data')
        if not _file_name:
            return
        # Make file header
        _header = (f'Computed Chi^2 Map\n'
                   f'{__appname__} v{__version__}\n'
                   f'Mapped parameters:\n'
                   f'\tx : {self.x_param}\n'
                   f'\ty : {self.y_param}\n'
                   f'Map size : {self.map_shape}\n'
                   f'Composition {{element: (Z, charge, n)}} : {self.ref_data["composition"]}\n'
                   f'S(Q) formalism : {self.ref_data["sq_method"]}\n'
                   )
        if 'rho' not in [self.x_param, self.y_param]:
            _header += f'rho : {self.ref_data["rho"]}\n'
        if 'r_min' not in [self.x_param, self.y_param]:
            _header += f'r_min : {self.ref_data["r_min"]}\n'
        if 'n_iter' not in [self.x_param, self.y_param]:
            _header += f'n_iter : {self.ref_data["n_iter"]}\n'

        if save_map_flag:
            # Check data operations in self.norm_map_data
            _norm_flag, _log_flag, *_ = self.compute_map_config_widget.map_result_gb.get_map_options()
            if _norm_flag == 0:
                _chi_str = 'Chi^2'
            elif _norm_flag == 1:
                _chi_str = 'Chi^2 (normalised)'
            elif _norm_flag == 2:
                _chi_str = 'Chi^2 (normalised per row)'
            elif _norm_flag == 3:
                _chi_str = 'Chi^2 (normalised per column)'
            else:
                raise NotImplementedError()
            if _log_flag == 0:
                pass
            elif _log_flag == 1:
                _chi_str = 'log(' + _chi_str + ')'
            elif _log_flag == 2:
                _chi_str = 'log(log(' + _chi_str + '))'
            else:
                raise NotImplementedError()
            # Set data header
            _header += f'Data: xx, yy, {_chi_str}'
            # Set data
            _xx, _yy = np.meshgrid(self.x_array, self.y_array)
            _data_out = np.array([_xx, _yy, self.norm_map_data])
            # Save data, one array at a time
            _header + f'xx :\n'
            # Make separate headers for arrays
            _header = [_header, f'yy :', f'{_chi_str} :']
            with open(_file_name, 'ab') as _data_file:
                for _n in range(3):
                    np.savetxt(_data_file, _data_out[_n], header=_header[_n], comments='#')

        else:
            # Set data header
            _header += f' x | y | Chi^2'
            # Set data
            _data_out = np.append(np.array(np.meshgrid(self.x_array, self.y_array)).T.reshape(-1,2),
                                  self.raw_map_data.reshape(-1, 1), axis=1)
            # Save data
            np.savetxt(_file_name, _data_out, header=_header, comments='#')

    def export_figure(self):
        if self.map_data is None:
            return
        # Get filename
        _file_name = get_filename(io='save_fig', caption='Save χ\u00b2 Map Figure')
        if not _file_name:
            return
        # Make pyqtgraph exporter
        _exporter = ImageExporter(self.map_plot_widget.map_plot)
        _check_save = _exporter.export(_file_name)
        # Handle for error saving (likely missing Qt filetype handler)
        if not _check_save:
            _file_name += '.png'
            _check_save = _exporter.export(_file_name)
            # If png doesn't work raise an error
            if not _check_save:
                _message = ['Error exporting figure!', 'Please check your Qt version for supported image formats\nSee QtGui.QImageWriter']
                _error_msg = ErrorMessageBox(_message)
                _error_msg.exec()
                return
        print('χ\u00b2 map figure saved: ', _file_name)


class ComputeMapConfigWidget(QWidget):

    def __init__(self):
        super(ComputeMapConfigWidget, self).__init__()

        self.vlayout = QVBoxLayout()
        self.vlayout.setContentsMargins(0, 0, 5, 0)
        self.vlayout.setSpacing(10)
        self.map_options_gb = ComputeMapOptionsGroupBox()
        self.progress_gb = ProgressBarGroupBox()
        self.map_result_gb = MapResultGroupBox()
        self.map_pos_gb = MapPositionGroupBox()
        self.vlayout.addWidget(self.map_options_gb, 1)
        self.vlayout.addWidget(self.progress_gb, 1)
        self.vlayout.addWidget(self.map_result_gb, 1)
        self.vlayout.addWidget(self.map_pos_gb, 1)
        self.vlayout.addWidget(QWidget(), 3)
        self.setLayout(self.vlayout)


class ComputeMapOptionsGroupBox(QGroupBox):

    def __init__(self, *args):
        super(ComputeMapOptionsGroupBox, self).__init__(*args)
        self.setTitle('Chi-Squared Map Options')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')

        self.create_widgets()
        self.style_widgets()
        self.create_layout()
        self.create_signals()
        self.set_param_array_size_label()

    def create_widgets(self):
        # Four options for variables to compute
        _param_items = ['Density (\u03c1)',
                        'Background Scaling (b)',
                        'r\u2098\u1d62\u2099',
                        'No. iterations']

        # Combo boxes for selecting x and y
        self.x_param_label = QLabel('x: ')
        self.x_param_input = QComboBox()
        self.x_param_input_list = QListWidget()
        self.x_param_input_list.hide()
        self.x_param_input.setModel(self.x_param_input_list.model())
        self.x_param_input_list.addItems(_param_items)
        self.y_param_label = QLabel('y: ')
        self.y_param_input = QComboBox()
        self.y_param_input_list = QListWidget()
        self.y_param_input_list.hide()
        self.y_param_input.setModel(self.y_param_input_list.model())
        self.y_param_input_list.addItems(_param_items)

        self.min_label = QLabel('Min')
        self.max_label = QLabel('Max')
        self.step_label = QLabel('Step')

        self.x_min_input = QLineEdit()
        self.x_max_input = QLineEdit()
        self.x_step_input = QLineEdit()
        self.x_inputs = [self.x_min_input,
                         self.x_max_input,
                         self.x_step_input]
        self.y_min_input = QLineEdit()
        self.y_max_input = QLineEdit()
        self.y_step_input = QLineEdit()
        self.y_inputs = [self.y_min_input,
                         self.y_max_input,
                         self.y_step_input]

        self.inputs = self.x_inputs + self.y_inputs

        self.estimate_time_btn = QPushButton('Estimate Time')
        self.estimate_time_label = QLabel('??:??:??')
        self.map_size_label = QLabel()

        self.compute_map_btn = QPushButton('Compute Map')

    def style_widgets(self):
        self.x_param_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.y_param_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

        self.min_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
        self.max_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
        self.step_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)

        self.estimate_time_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
        self.map_size_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)

        for _widget in self.inputs:
            _widget.setMaximumWidth(80)
            _widget.setValidator(CustomDoubleValidator(0, np.inf, -1))

        self.estimate_time_btn.setStyleSheet("padding-left: 22px; padding-right: 22px;"
                                             "padding-top: 5px; padding-bottom: 5px;");
        self.compute_map_btn.setStyleSheet("padding-left: 22px; padding-right: 22px;"
                                           "padding-top: 5px; padding-bottom: 5px;");

        self.spacer = QWidget()
        self.spacer.setFixedHeight(15)

    def create_layout(self):
        self.grid_layout = QGridLayout()
        self.grid_layout.setContentsMargins(25, 10, 25, 7)
        self.grid_layout.setSpacing(5)

        self.grid_layout.addWidget(self.min_label, 0, 2)
        self.grid_layout.addWidget(self.max_label, 0, 3)
        self.grid_layout.addWidget(self.step_label, 0, 4)
        self.grid_layout.addWidget(self.x_param_label, 1, 0)
        self.grid_layout.addWidget(self.x_param_input, 1, 1)
        self.grid_layout.addWidget(self.x_min_input, 1, 2)
        self.grid_layout.addWidget(self.x_max_input, 1, 3)
        self.grid_layout.addWidget(self.x_step_input, 1, 4)
        self.grid_layout.addWidget(self.y_param_label, 2, 0)
        self.grid_layout.addWidget(self.y_param_input, 2, 1)
        self.grid_layout.addWidget(self.y_min_input, 2, 2)
        self.grid_layout.addWidget(self.y_max_input, 2, 3)
        self.grid_layout.addWidget(self.y_step_input, 2, 4)
        self.grid_layout.addWidget(self.spacer, 3, 1)
        self.grid_layout.addWidget(self.estimate_time_btn, 4, 1, Qt.AlignCenter)
        self.grid_layout.addWidget(self.estimate_time_label, 4, 2, 1, 3)
        self.grid_layout.addWidget(self.compute_map_btn, 5, 1, Qt.AlignCenter)
        self.grid_layout.addWidget(self.map_size_label, 5, 2, 1, 3)

        self.setLayout(self.grid_layout)

    def create_signals(self):
        self.x_param_input.currentIndexChanged.connect(lambda n: self.set_input_validator(n, 0))
        self.y_param_input.currentIndexChanged.connect(lambda n: self.set_input_validator(n, 1))
        for _widget in self.inputs:
            _widget.textChanged.connect(self.set_param_array_size_label)
        for _widget in self.x_inputs:
            _widget.editingFinished.connect(lambda: self.set_validator_range(0))
        for _widget in self.y_inputs:
            _widget.editingFinished.connect(lambda: self.set_validator_range(1))

    def get_param_array_size(self):
        try:
            _x_min = float(self.x_min_input.text())
            _x_max = float(self.x_max_input.text())
            _x_step = float(self.x_step_input.text())
            _x_size = int((_x_max - _x_min)/_x_step + 1)
        except (ValueError, ZeroDivisionError):
            _x_size = None
        try:
            _y_min = float(self.y_min_input.text())
            _y_max = float(self.y_max_input.text())
            _y_step = float(self.y_step_input.text())
            _y_size = int((_y_max - _y_min)/_y_step + 1)
        except (ValueError, ZeroDivisionError):
            _y_size = None
        try:
            _n_points = _x_size * _y_size
        except TypeError:
            _n_points = None
        return _x_size, _y_size, _n_points

    def get_param_array(self):
        try:
            _x_min = float(self.x_min_input.text())
            _x_max = float(self.x_max_input.text())
            _x_step = float(self.x_step_input.text())
            _y_min = float(self.y_min_input.text())
            _y_max = float(self.y_max_input.text())
            _y_step = float(self.y_step_input.text())
        except ValueError:
            raise RuntimeError()
        # Get arange min-max inclusive of max if possible
        _x_array = np.arange(_x_min, _x_max+_x_step, _x_step)
        _x_array = _x_array[np.round(_x_array, 7) <= _x_max]
        _y_array = np.arange(_y_min, _y_max+_y_step, _y_step)
        _y_array = _y_array[np.round(_y_array, 7) <= _y_max]
        return _x_array, _y_array, _x_min, _y_min, _x_step, _y_step

    def set_input_validator(self, _index, _var):
        # Update validator type (int/double)
        if not _var:
            _inputs = self.x_inputs
        else:
            _inputs = self.y_inputs
        if _index == 3:
            for _widget in _inputs:
                _widget.setValidator(CustomIntValidator(1, 2147483647))
                try:
                    _widget.setText(str(int(float(_widget.text()))))
                except ValueError:
                    _widget.setText(None)
            try:
                if int(_inputs[2].text()) == 0:
                    _inputs[2].setText(str(1))
            except ValueError:
                pass
        else:
            for _widget in _inputs:
                _widget.setValidator(CustomDoubleValidator(0, np.inf, -1))

    def set_validator_range(self, _var):
        if not _var:
            _inputs = self.x_inputs
            _param_type = self.x_param_input.currentIndex()
        else:
            _inputs = self.y_inputs
            _param_type = self.y_param_input.currentIndex()
        # Get min/max values if present
        try:
            _min = float(_inputs[0].text())
        except ValueError:
            if _param_type == 3:
                _min = 1
            else:
                _min = 0
        try:
            _max = float(_inputs[1].text())
        except ValueError:
            if _param_type == 3:
                _max = 2147483647
            else:
                _max = np.inf
        if _param_type == 3:
            _min = int(_min)
            _max = int(_max)
        # Update min/max values on validator
        _inputs[0].validator().setTop(_max)
        _inputs[1].validator().setBottom(_min)

    def set_param_array_size_label(self):
        _x_size, _y_size, _n_points = self.get_param_array_size()
        # Convert to strings
        _x_size = str(_x_size or '?')
        _y_size = str(_y_size or '?')
        _n_points = str(_n_points or '?')
        # Make string
        _size_str = _x_size + ' x ' + _y_size + ' (' + _n_points + ')'
        self.map_size_label.setText(_size_str)

    def get_param_names(self):
        _x_param = self.x_param_input.currentIndex()
        _y_param = self.y_param_input.currentIndex()
        _param_label = ['rho', 'bkg_scaling', 'r_min', 'n_iter']
        _x_param, _y_param = _param_label[_x_param], _param_label[_y_param]
        return _x_param, _y_param

    def enable_widgets(self, state):
        # Enable/disable widgets
        self.estimate_time_btn.setEnabled(state)
        self.estimate_time_label.setEnabled(state)
        self.compute_map_btn.setEnabled(state)
        self.map_size_label.setEnabled(state)
        # Force Qt to process events
        QApplication.processEvents(QEventLoop.ExcludeUserInputEvents)


class ProgressBarGroupBox(QGroupBox):

    def __init__(self, *args):
        super(ProgressBarGroupBox, self).__init__(*args)
        self.setAlignment(Qt.AlignLeft)
        self.create_widgets()
        self.style_widgets()
        self.create_layout()

    def create_widgets(self):
        self.progress_bar = QProgressBar()
        self.eta_label = QLabel()

    def style_widgets(self):
        self.eta_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
        self.spacer = QWidget()
        self.spacer.setFixedHeight(15)

    def create_layout(self):
        self.v_layout = QVBoxLayout()
        self.v_layout.addWidget(self.progress_bar)
        self.v_layout.addWidget(self.spacer)
        self.v_layout.addWidget(self.eta_label)
        self.setLayout(self.v_layout)

    def update_progress(self, n):
        elapsed_time = time.perf_counter() - self.start_time
        eta = self.estimated_time - elapsed_time
        # Re-calculate estimate if elpased time > estimate
        if eta < 0:
            eta = (self.n_points - n)/self.n_points * n
        if n == self.n_points:
            time_str = 'Map Complete!'
        else:
            time_str = self.make_time_string(eta)
        self.eta_label.setText(time_str)
        self.progress_bar.setValue(n)

    def init_progress(self, n_points, estimated_time, start_time):
        self.reset_progress()
        self.n_points = n_points
        self.estimated_time = estimated_time
        self.start_time = start_time
        self.progress_bar.setRange(0, self.n_points)
        self.progress_bar.setValue(0)
        _time_str = self.make_time_string(self.estimated_time)
        self.eta_label.setText(_time_str)

    def reset_progress(self):
        self.progress_bar.reset()
        self.eta_label.setText(None)

    def make_time_string(self, eta):
        # eta in secnods
        if eta <= 60:
            time_str = str(int(eta)) + ' s'
        else:
            time_str = time.strftime('%H:%M:%S', time.gmtime(eta))
        return time_str


class MapResultGroupBox(QGroupBox):

    # Create pyqtSignal for updating options
    map_options_changed = Signal()

    def __init__(self, *args):
        super(MapResultGroupBox, self).__init__(*args)
        self.setTitle('Map Options')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')

        self.create_widgets()
        self.style_widgets()
        self.create_layout()
        self.create_signals()

    def create_widgets(self):
        # Options for data normalisation
        self.norm_button_group = QButtonGroup()
        self.raw_data_btn = QRadioButton('Raw data')
        self.norm_btn = QRadioButton('Normalise data')
        self.row_norm_btn = QRadioButton('Normalise rows')
        self.col_norm_btn = QRadioButton('Normalise columns')
        self.norm_button_group.addButton(self.raw_data_btn, 0)
        self.norm_button_group.addButton(self.norm_btn, 1)
        self.norm_button_group.addButton(self.row_norm_btn, 2)
        self.norm_button_group.addButton(self.col_norm_btn, 3)

        # Options for data transformations
        self.trans_button_group = QButtonGroup()
        self.chi_sq_btn = QRadioButton('χ\u00b2')
        self.log_btn = QRadioButton('log(χ\u00b2)')
        self.dbl_log_btn = QRadioButton('log(log(χ\u00b2))')
        self.trans_button_group.addButton(self.chi_sq_btn, 0)
        self.trans_button_group.addButton(self.log_btn, 1)
        self.trans_button_group.addButton(self.dbl_log_btn, 2)

        # Options for plotting contours
        self.plot_contours = QCheckBox('Plot contours?')
        self.levels_label = QLabel('Contour levels: ')
        self.levels = QSpinBox()

        # Data output
        self.export_map_btn = QPushButton('Save Figure')
        self.save_data_btn = QPushButton('Save Data')
        self.save_data_options_button_group = QButtonGroup()
        self.save_raw_data = QRadioButton('Save raw data')
        self.save_map_data = QRadioButton('Save map data')
        self.save_data_options_button_group.addButton(self.save_raw_data, 0)
        self.save_data_options_button_group.addButton(self.save_map_data, 1)


    def style_widgets(self):
        self.levels_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.levels.setMaximumWidth(80)
        self.levels.setAlignment(Qt.AlignLeft)

        self.spacer = QWidget()
        self.spacer.setFixedHeight(15)
        self.small_spacer = QWidget()
        self.small_spacer.setFixedHeight(4)

        self.levels.setRange(1, 99)

        # Setup default options
        self.raw_data_btn.setChecked(True)
        self.chi_sq_btn.setChecked(True)
        self.plot_contours.setChecked(True)
        self.levels.setValue(10)
        self.save_map_data.setChecked(True)


    def create_layout(self):
        self.level_options_layout = QHBoxLayout()
        self.level_options_layout.addWidget(self.levels_label)
        self.level_options_layout.addWidget(self.levels)
        self.level_options_layout.setContentsMargins(0, 0, 0, 0)
        self.level_options_layout.setSpacing(0)

        self.grid_layout = QGridLayout()
        self.grid_layout.setContentsMargins(25, 10, 25, 7)
        self.grid_layout.setSpacing(5)

        self.grid_layout.addWidget(self.raw_data_btn, 0, 0)
        self.grid_layout.addWidget(self.norm_btn, 1, 0)
        self.grid_layout.addWidget(self.row_norm_btn, 2, 0)
        self.grid_layout.addWidget(self.col_norm_btn, 3, 0)
        self.grid_layout.addWidget(self.chi_sq_btn, 0, 1)
        self.grid_layout.addWidget(self.log_btn, 1, 1)
        self.grid_layout.addWidget(self.dbl_log_btn, 2, 1)
        self.grid_layout.addWidget(self.spacer, 4, 0)
        self.grid_layout.addWidget(self.plot_contours, 5, 0)
        self.grid_layout.addLayout(self.level_options_layout, 5, 1)
        self.grid_layout.addWidget(self.spacer, 6, 0)
        self.grid_layout.addWidget(self.export_map_btn, 7, 0, 1, 2)
        self.grid_layout.addWidget(self.small_spacer, 8, 0)
        self.grid_layout.addWidget(self.save_data_btn, 9, 0, 1, 2)
        self.grid_layout.addWidget(self.small_spacer, 10, 0)
        self.grid_layout.addWidget(self.save_raw_data, 11, 0)
        self.grid_layout.addWidget(self.save_map_data, 11, 1)

        self.setLayout(self.grid_layout)

    def create_signals(self):
        self.norm_button_group.buttonClicked.connect(self.map_options_changed)
        self.trans_button_group.buttonClicked.connect(self.map_options_changed)
        self.plot_contours.stateChanged.connect(self.map_options_changed)
        self.levels.valueChanged.connect(self.map_options_changed)

    def get_map_options(self):
        # Get data_norm flag
        norm_flag = self.norm_button_group.checkedId()
        # Get data log flag
        log_flag = self.trans_button_group.checkedId()
        # Get contour flag
        contour_flag = self.plot_contours.isChecked()
        n_levels = self.levels.value()
        # Get save raw/map flag
        save_map_flag = self.save_data_options_button_group.checkedId()
        return norm_flag, log_flag, contour_flag, n_levels, save_map_flag


class MapPositionGroupBox(QGroupBox):

    def __init__(self, *args):
        super(MapPositionGroupBox, self).__init__(*args)
        self.setTitle('Inspect Map')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')

        self.param_labels = {'rho': 'Density (\u03c1\u2080): ',
                             'bkg_scaling': 'Background Scaling (b): ',
                             'r_min': 'r\u2098\u1d62\u2099: ',
                             'n_iter': 'No. iterations: ',
                             0: 'χ\u00b2',
                             1: 'χ\u00b2<sub>norm</sub>',
                             2: 'χ\u00b2<sub>norm-row</sub>',
                             3: 'χ\u00b2<sub>norm-col</sub>'
                             }

        self.create_widgets()
        self.style_widgets()
        self.create_layout()

    def create_widgets(self):
        self.x_param_label = QLabel()
        self.y_param_label = QLabel()
        self.x_pos_label = QLabel()
        self.y_pos_label = QLabel()
        self.chi_param_label = QLabel()
        self.chi_val_label = QLabel()

    def style_widgets(self):
        self.x_param_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.x_pos_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.y_param_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.y_pos_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.chi_param_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.chi_val_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)

    def create_layout(self):
        self.grid_layout = QGridLayout()
        self.grid_layout.addWidget(self.x_param_label, 0, 0)
        self.grid_layout.addWidget(self.x_pos_label, 0, 1)
        self.grid_layout.addWidget(self.y_param_label, 1, 0)
        self.grid_layout.addWidget(self.y_pos_label, 1, 1)
        self.grid_layout.addWidget(self.chi_param_label, 2, 0)
        self.grid_layout.addWidget(self.chi_val_label, 2, 1)
        self.setLayout(self.grid_layout)

    def set_pos_label(self, _x, _y, _val):
        self.x_pos_label.setText(f'{_x:.5f}')
        self.y_pos_label.setText(f'{_y:.5f}')
        self.chi_val_label.setText(f'{_val:.5f}')

    def clear_pos_label(self):
        self.x_pos_label.setText(None)
        self.y_pos_label.setText(None)
        self.chi_val_label.setText(None)

    def set_param_label(self, _x_p, _y_p, _norm_flag, _log_flag):
        self.x_param_label.setText(self.param_labels[_x_p])
        self.y_param_label.setText(self.param_labels[_y_p])
        _chi_norm_label = self.param_labels[_norm_flag]
        _chi_log_labels = {0: f'{_chi_norm_label}: ',
                           1: f'log({_chi_norm_label}): ',
                           2: f'log(log({_chi_norm_label})): '
                           }
        self.chi_param_label.setText(_chi_log_labels[_log_flag])

    def clear_param_label(self):
        self.x_param_label.setText(None)
        self.y_param_label.setText(None)
        self.chi_param_label.setText(None)

    def clear_all(self):
        self.clear_pos_label()
        self.clear_param_label()


class AttenuationInspectDialog(QDialog):

    def __init__(self, data):
        super(AttenuationInspectDialog, self).__init__()
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.title = 'Self-shielding Attenuation Factor | ' + __appname__ + ' v' + __version__
        self.setWindowTitle(self.title)
        with resources.as_file(resources.files('LiquidDiffract.resources.icons').joinpath('gs_icon.png')) as path:
            self.setWindowIcon(QIcon(str(path)))
        # Unpack data
        self.q_data, self.two_theta_data, self.attenuation_data = data
        self.layout = QHBoxLayout(self)
        self.layout.setSpacing(0)
        # Make Config Widget
        self.config_widget = AttenuationPlotConfigWidget()
        self.config_scroll_area = QScrollArea()
        self.config_scroll_area.setFrameShape(QFrame.NoFrame)
        self.config_scroll_area.setWidget(self.config_widget)
        self.config_scroll_area.setWidgetResizable(True)
        self.config_scroll_area.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.Ignored)
        # Make vertical line separator
        self.vline = QFrame()
        self.vline.setFrameShape(QFrame.VLine)
        self.vline.setFrameShadow(QFrame.Sunken)
        self.vline.setObjectName("vline")
        # Make Plot Widget
        self.plot_widget = plot_widgets.AttenuationCorrectionPlotWidget()
        self.plot_scroll_area = QScrollArea()
        self.plot_scroll_area.setWidget(self.plot_widget)
        self.plot_scroll_area.setWidgetResizable(True)
        self.plot_scroll_area.setFrameShape(QFrame.NoFrame)
        self.hsplitter = QSplitter(Qt.Horizontal)
        self.hsplitter.addWidget(self.config_scroll_area)
        self.hsplitter.addWidget(self.plot_scroll_area)
        self.hsplitter.setStretchFactor(0, 2)
        self.hsplitter.setStretchFactor(1, 5)
        self.layout.addWidget(self.hsplitter)
        self.setLayout(self.layout)

        self.create_signals()
        self.plot_data()

    def create_signals(self):
        self.config_widget.attenuation_config_gb.x_unit_button_group.buttonToggled.connect(self.plot_data)
        self.config_widget.attenuation_config_gb.save_data_btn.clicked.connect(self.save_data)

    def plot_data(self):
        # Get option (Q-space or two-theta)
        _plot_two_theta = self.config_widget.attenuation_config_gb.plot_two_theta.isChecked()
        if _plot_two_theta:
            _x = self.two_theta_data
        else:
            _x = self.q_data
        # Plot data
        self.plot_widget.plot_data(_x, self.attenuation_data, _plot_two_theta)

    def save_data(self):
        # Get save filename
        _file_name = get_filename(io='save', caption='Save Self-shielding attenuation')
        if not _file_name:
            return
        # Make file header
        _header = (f'Self-shielding attenuation factor\n'
                   f'{__appname__} v{__version__}\n'
                   f'Q | 2-theta | A_s,s\n'
                   )
        _data_out = np.column_stack((self.q_data, self.two_theta_data, self.attenuation_data))
        np.savetxt(_file_name, _data_out, header=_header, comments='#')


class AttenuationPlotConfigWidget(QWidget):

    def __init__(self):
        super(AttenuationPlotConfigWidget, self).__init__()
        self.vlayout = QVBoxLayout()
        self.vlayout.setContentsMargins(0, 0, 5, 0)
        self.vlayout.setSpacing(10)
        self.attenuation_config_gb = AttenuationPlotConfigGroupBox()
        self.vlayout.addWidget(self.attenuation_config_gb, 1)
        self.vlayout.addWidget(QWidget(), 3)
        self.setLayout(self.vlayout)


class AttenuationPlotConfigGroupBox(QGroupBox):

    def __init__(self, *args):
        super(AttenuationPlotConfigGroupBox, self).__init__(*args)
        self.setAlignment(Qt.AlignLeft)
        self.create_widgets()
        self.style_widgets()
        self.create_layout()

    def create_widgets(self):
        self.x_unit_label = QLabel('Units: ')
        self.plot_q = QRadioButton('Q')
        self.plot_two_theta = QRadioButton('2θ')
        self.x_unit_button_group = QButtonGroup()
        self.x_unit_button_group.addButton(self.plot_q, 0)
        self.x_unit_button_group.addButton(self.plot_two_theta, 1)
        self.save_data_btn = QPushButton('Save Data')

    def style_widgets(self):
        self.x_unit_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.plot_q.setStyleSheet('font: italic')
        self.plot_two_theta.setStyleSheet('font: italic')
        self.plot_q.setChecked(True)
        self.spacer = QWidget()
        self.spacer.setFixedHeight(10)

    def create_layout(self):
        self.grid_layout = QGridLayout()
        self.grid_layout.setContentsMargins(25, 10, 25, 7)
        self.grid_layout.setSpacing(7)
        self.grid_layout.addWidget(self.x_unit_label, 0, 0)
        self.grid_layout.addWidget(self.plot_q, 0, 1)
        self.grid_layout.addWidget(self.plot_two_theta, 1, 1)
        self.grid_layout.addWidget(self.spacer, 2, 0)
        self.grid_layout.addWidget(self.save_data_btn, 3, 0, 1, 2)
        self.setLayout(self.grid_layout)