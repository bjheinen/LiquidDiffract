# -*- coding: utf-8 -*-
"""
GUI frontend for LiquidDiffract
<https://github.com/bjheinen/LiquidDiffract>
"""
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2018-2021, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

import os.path
import webbrowser
# importlib.resources only available in python>=3.7
try:
    from importlib import resources as importlib_resources
except ImportError:
    import importlib_resources
import numpy as np
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QTabWidget, QAction
from PyQt5.QtGui import QIcon

from LiquidDiffract.gui import bkg_ui
from LiquidDiffract.gui import optim_ui
from LiquidDiffract.gui import results_ui
from LiquidDiffract.gui import structure_ui
from LiquidDiffract.gui import utility
from LiquidDiffract.version import __appname__, __version__


class App(QMainWindow):

    def __init__(self, screen_size):
        super().__init__()
        # Set options here before running initUI to initialise
        # Get primary (current) screen dimensions
        self.screen_size = screen_size
        self.width = self.screen_size.width()
        self.height = self.screen_size.height()
        self.initUI()

    def initUI(self):
        self.title = __appname__ + ' v' + __version__
        self.setWindowTitle(self.title)
        self.icon_module = 'LiquidDiffract.resources.icons'
        with importlib_resources.path(self.icon_module, 'gs_icon.png') as path:
            self.setWindowIcon(QIcon(str(path)))
        self.setGeometry(0, 0, self.width, self.height)

        self.menu_bar = self.menuBar()
        self.tools_menu = self.menu_bar.addMenu('&Tools')
        self.help_menu = self.menu_bar.addMenu('&Help')

        with importlib_resources.path(self.icon_module, 'config.png') as path:
            self.preferences_action = QAction(QIcon(str(path)),
                                              'Additional Preferences...', self)

        with importlib_resources.path(self.icon_module, 'browser.png') as path:
            self.documentation_action = QAction(QIcon(str(path)),
                                                self.title + ' Documentation', self)

        with importlib_resources.path(self.icon_module, 'info.png') as path:
            self.about_action = QAction(QIcon(str(path)), 'About', self)

        self.tools_menu.addAction(self.preferences_action)
        self.help_menu.addAction(self.documentation_action)
        self.help_menu.addAction(self.about_action)
        self.preferences_action.triggered.connect(self.call_preferences_dialog)
        self.about_action.triggered.connect(self.call_about_dialog)
        self.documentation_action.triggered.connect(self.open_docs)

        self.table_widget = MainContainer(self)
        self.setCentralWidget(self.table_widget)
        self.set_default_preferences()
        self.showMaximized()

    def set_default_preferences(self):
        self.preferences = {'append_log_mode': 1,
                            'data_units': 0,
                            'window_length': 5,
                            'poly_order': 3,
                            'rescale_AL': 1,
                            'fft_N': 12,
                            'mod_func_mode': 0,
                            'op_method': 'L-BFGS-B',
                            'minimisation_options':
                                {'disp': 0,
                                 'maxiter': 15000,
                                 'maxfun': 15000,
                                 'ftol': 2.22e-8,
                                 'gtol': 1e-10},
                            'global_minimisation': 0,
                            'global_min_options':
                                {'disp': 0,
                                 'niter': 100,
                                 'T': 1.0,
                                 'stepsize': 0.01,
                                 'interval': 50},
                            'xray_weight_mode': 0
                            }
        self.set_preferences()

    def call_preferences_dialog(self):
        self.preferences_dialog = utility.PreferencesDialog(self.preferences)
        self.preferences_dialog.fft_check_signal.connect(self.check_fft_N)
        # set window icon
        with importlib_resources.path(self.icon_module, 'gs_icon.png') as path:
            self.preferences_dialog.setWindowIcon(QIcon(str(path)))

        if self.preferences_dialog.exec_() == utility.PreferencesDialog.Accepted:
            self.preferences = self.preferences_dialog.get_preferences()
            # Set preferences in OptimUI
            self.set_preferences()
            # Call method plot_data to update data treatment
            self.table_widget.optim_ui.plot_data()

    def set_preferences(self):
        # Set data and minimisation preferences
        self.table_widget.optim_ui.mod_func_mode = self.preferences['mod_func_mode']
        self.table_widget.optim_ui.op_method = self.preferences['op_method']
        self.table_widget.optim_ui.minimisation_options = self.preferences['minimisation_options']
        self.table_widget.optim_ui.global_minimisation = self.preferences['global_minimisation']
        self.table_widget.optim_ui.global_min_options = self.preferences['global_min_options']
        self.table_widget.optim_ui.append_log_mode = self.preferences['append_log_mode']
        self.table_widget.structure_ui.append_log_mode = self.preferences['append_log_mode']
        self.table_widget.optim_ui.window_length = self.preferences['window_length']
        self.table_widget.optim_ui.poly_order = self.preferences['poly_order']
        # Set N for FFT
        self.table_widget.bkg_ui.fft_N = self.preferences['fft_N']
        self.table_widget.optim_ui.fft_N = self.preferences['fft_N']
        self.table_widget.results_ui.fft_N = self.preferences['fft_N']
        # Set data units for file loading
        self.table_widget.bkg_ui.data_units = self.preferences['data_units']
        # Set preferences for rescaling Ashcroft-Langreth S(Q)/g(r) plots
        self.table_widget.results_ui.rescale_AL = self.preferences['rescale_AL']
        # Set preferences for x-ray weight mode (0: calc Kp at Q=0, 1: calc Q-dependent Kp)
        self.table_widget.structure_ui.xray_weight_mode = self.preferences['xray_weight_mode']

    def call_about_dialog(self):
        self.about_dialog = utility.AboutDialog()
        self.about_dialog.exec_()

    def open_docs(self):
        webbrowser.open_new('https://github.com/bjheinen/LiquidDiffract')

    def check_fft_N(self):
        try:
            _qmax = self.table_widget.bkg_ui.data['data_raw_x'][-1]
        except IndexError:
            try:
                _qmax = self.table_widget.bkg_ui.data['bkg_raw_x'][-1]
            except IndexError:
                self.preferences_dialog.fft_check_result = 0
                return
        _fft_check = self.preferences_dialog.fft_check
        _dx = (self.table_widget.bkg_ui.data['data_raw_x'][1]
               - self.table_widget.bkg_ui.data['data_raw_x'][0])
        if _qmax > ((2**_fft_check / 2) * _dx):
            self.preferences_dialog.fft_check_result = 1


class MainContainer(QWidget):

    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.setStyleSheet('QTabWidget::pane { border: 0; } \
                            QTabBar {font-size: 11pt; text-align: center} \
                            QTabWidget::pane {border-top: 4px solid #444; }')
        self.layout = QVBoxLayout(self)
        self.tabs = QTabWidget()
        self.bkg_ui = bkg_ui.BkgUI(self)
        self.optim_ui = optim_ui.OptimUI(self)
        self.results_ui = results_ui.ResultsUI(self)
        self.structure_ui = structure_ui.StructureUI(self)
        self.tabs.addTab(self.bkg_ui, 'Background Subtraction')
        self.tabs.addTab(self.optim_ui, 'Refine Structure Factor')
        self.tabs.addTab(self.results_ui, 'Calculate PDF')
        self.tabs.addTab(self.structure_ui, 'Structural Information')
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

        self.create_signal_tab_links()

    def create_signal_tab_links(self):
        self.bkg_ui.plots_changed.connect(self.bkg_plots_changed_slot)
        self.optim_ui.results_changed.connect(self.results_changed_slot)
        self.optim_ui.results_cleared.connect(self.results_cleared_slot)
        self.bkg_ui.file_name_changed.connect(self.update_filename)
        self.optim_ui.optim_config_widget.optim_results_gb.bkg_scale_copy_btn.clicked.connect(self.copy_bkg_scale_result)

    def bkg_plots_changed_slot(self):
        # Clear data from Optim UI (as S(Q) etc. need to be recalculated)
        self.optim_ui.data = {'uncorrected_y': np.asarray([]), 'bkg_y': np.asarray([]), 'bkg_scale': None,
                              'cor_x': np.asarray([]), 'cor_y': np.asarray([]),
                              'cor_x_cut': np.asarray([]), 'cor_y_cut': np.asarray([]),
                              'sq_x': np.asarray([]), 'sq_y': np.asarray([]),
                              'dr_x': np.asarray([]), 'dr_y': np.asarray([]),
                              'rescaled_cor_y_cut': np.asarray([]),
                              'mod_func': 'None', 'qmin': None, 'qmax': None
                              }
        # Pass data to Optim UI
        self.optim_ui.data['cor_x'] = self.bkg_ui.data['cor_x']
        self.optim_ui.data['cor_y'] = self.bkg_ui.data['cor_y']
        if self.bkg_ui.bkg_config_widget.bkg_subtract_gb.isChecked() and self.bkg_ui.data['bkg_y'].size:
            self.optim_ui.data['uncorrected_y'] = self.bkg_ui.data['data_y']
            self.optim_ui.data['bkg_y'] = self.bkg_ui.data['bkg_y']
            self.optim_ui.data['bkg_scale'] = self.bkg_ui.bkg_config_widget.bkg_subtract_gb.scale_sb.value()
        self.optim_ui.plot_data()

    def results_changed_slot(self):
        # Clear old s(q), g(r), rdf(r)
        self.results_ui.clear_data()
        # Set rho
        if (self.optim_ui.optim_config_widget.optim_options_gb.opt_check.isChecked()) & (self.optim_ui.optim_config_widget.optim_options_gb.isChecked()):
            self.results_ui.data['rho'] = self.optim_ui.data['refined_rho']
        else:
            self.results_ui.data['rho'] = np.float(self.optim_ui.optim_config_widget.composition_gb.density_input.text())
        if self.optim_ui.optim_config_widget.optim_options_gb.isChecked():
            self.results_ui.data['int_func'] = self.optim_ui.data['impr_int_func']
            self.results_ui.data['sq_x'] = self.optim_ui.data['impr_iq_x']
        else:
            self.results_ui.data['int_func'] = self.optim_ui.data['int_func']
            self.results_ui.data['sq_x'] = self.optim_ui.data['iq_x']
        self.results_ui.data['composition'] = self.optim_ui.optim_config_widget.composition_gb.get_composition_dict()
        self.results_ui.data['mod_func'] = self.optim_ui.data['mod_func']
        self.results_ui.data['window_start'] = self.optim_ui.data['window_start']
        self.results_ui.data['sq_method'] = self.optim_ui.data['sq_method']
        self.results_ui.plot_data()

        # Clear data from structure ui
        self.structure_ui.clear_data()

        self.structure_ui.data['rdf_x'] = self.results_ui.data['rdf_x']
        self.structure_ui.data['rdf_y'] = self.results_ui.data['rdf_y']
        self.structure_ui.data['tr_x'] = self.results_ui.data['tr_x']
        self.structure_ui.data['tr_y'] = self.results_ui.data['tr_y']
        self.structure_ui.data['sq_x'] = self.results_ui.data['sq_x']
        self.structure_ui.data['rho'] = self.results_ui.data['rho']
        self.structure_ui.data['composition'] = self.results_ui.data['composition']

        self.structure_ui.update_obj_fun()
        self.structure_ui.set_atoms()
        self.structure_ui.set_weights()
        self.structure_ui.structure_plot_widget.toggle_fit_limits(True,
                                                                  obj_x=self.structure_ui.data['rdf_x'],
                                                                  sq_x=self.structure_ui.data['sq_x'])
        self.structure_ui.plot_data()

    def results_cleared_slot(self):
        self.results_ui.clear_data()
        self.results_ui.results_plot_widget.update_plots(self.results_ui.data)
        
        self.structure_ui.clear_data()
        self.structure_ui.structure_plot_widget.clear_plots(_clear_all=True)

    def update_filename(self):
        _base_name, _ext = os.path.splitext(self.bkg_ui.data_file)
        self.results_ui.base_filename = _base_name
        self.optim_ui.base_filename = _base_name
        self.structure_ui.base_filename = _base_name
        self.optim_ui.filename_ext = _ext
        self.structure_ui.filename_ext = _ext

    def copy_bkg_scale_result(self):
        _bkg_scale_result = self.optim_ui.optim_config_widget.optim_results_gb.bkg_scale_output.text()
        self.bkg_ui.bkg_config_widget.bkg_subtract_gb.scale_sb.setValue(np.float64(_bkg_scale_result))