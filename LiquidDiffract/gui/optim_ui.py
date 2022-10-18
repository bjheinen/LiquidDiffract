# -*- coding: utf-8 -*-
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2018-2021, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

import os.path
import datetime
# importlib.resources only available in python>=3.7
try:
    from importlib import resources as importlib_resources
except ImportError:
    import importlib_resources
import numpy as np
from scipy.optimize import minimize, basinhopping
from PyQt5.QtCore import Qt, pyqtSignal, QEventLoop
from PyQt5.QtGui import QDoubleValidator, QIntValidator
from PyQt5.QtWidgets import QWidget, QFrame, QGridLayout, QVBoxLayout, \
                            QHBoxLayout, QGroupBox, QPushButton, QLineEdit, \
                            QComboBox, QTableWidget, QTableWidgetItem, \
                            QLabel, QCheckBox, QButtonGroup, QRadioButton, \
                            QScrollArea, QSplitter, QSizePolicy,  \
                            QAbstractScrollArea, QHeaderView, QStyle, \
                            QApplication
from LiquidDiffract.gui import plot_widgets
from LiquidDiffract.gui import utility
from LiquidDiffract.core import data_utils
import LiquidDiffract.core.core as core
from LiquidDiffract.version import __appname__, __version__


class OptimUI(QWidget):
    # Create custom signal to link Optim/Results UI
    results_changed = pyqtSignal()
    results_cleared = pyqtSignal()

    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QHBoxLayout(self)
        self.layout.setSpacing(0)

        # Make Config Widget
        self.optim_config_widget = OptimConfigWidget()

        # Make vertical line separator
        self.vline = QFrame()
        self.vline.setFrameShape(QFrame.VLine)
        self.vline.setFrameShadow(QFrame.Sunken)
        self.vline.setObjectName("vline")

        self.optim_plot_widget = plot_widgets.OptimPlotWidget()
        self.plot_scroll_area = QScrollArea()
        self.plot_scroll_area.setWidget(self.optim_plot_widget)
        self.plot_scroll_area.setWidgetResizable(True)
        self.plot_scroll_area.setFrameShape(QFrame.NoFrame)

        self.config_scroll_area = QScrollArea()
        self.config_scroll_area.setFrameShape(QFrame.NoFrame)
        self.config_scroll_area.setWidget(self.optim_config_widget)
        self.config_scroll_area.setWidgetResizable(True)

        self.hsplitter = QSplitter(Qt.Horizontal)
        self.hsplitter.addWidget(self.config_scroll_area)
        self.hsplitter.addWidget(self.plot_scroll_area)
        self.hsplitter.setStretchFactor(0, 2)
        self.hsplitter.setStretchFactor(1, 5)

        self.layout.addWidget(self.hsplitter)
        self.setLayout(self.layout)

        self.data = {'cor_x': np.asarray([]), 'cor_y': np.asarray([]),
                     'cor_x_cut': np.asarray([]), 'cor_y_cut': np.asarray([]),
                     'sq_x': np.asarray([]), 'sq_y': np.asarray([]),
                     'dr_x': np.asarray([]), 'dr_y': np.asarray([]),
                     'int_func': np.asarray([]), 'impr_int_func': np.asarray([]),
                     'impr_dr_x': np.asarray([]), 'impr_dr_y': np.asarray([]),
                     'impr_iq_x': np.asarray([]), 'mod_func': 'None',
                     'rescaled_cor_y_cut': np.asarray([]),
                     'qmin': None, 'qmax': None}

        self.create_signals()

    def create_signals(self):
        self.optim_config_widget.data_options_gb.qmax_check.stateChanged.connect(self.plot_data)
        self.optim_config_widget.data_options_gb.qmax_input.textChanged.connect(self.plot_data)
        self.optim_config_widget.data_options_gb.qmin_check.stateChanged.connect(self.plot_data)
        self.optim_config_widget.data_options_gb.qmin_input.textChanged.connect(self.plot_data)
        self.optim_config_widget.composition_gb.mass_density.textChanged.connect(self.plot_data)
        self.optim_config_widget.data_options_gb.al_btn.toggled.connect(self.plot_data)
        self.optim_config_widget.data_options_gb.mod_func_input.currentIndexChanged.connect(self.plot_data)
        self.optim_config_widget.data_options_gb.calc_sq_btn.clicked.connect(self.on_click_calc_sq)
        self.optim_config_widget.data_options_gb.smooth_data_check.toggled.connect(self.plot_data)
        self.optim_config_widget.optim_options_gb.opt_button.clicked.connect(self.on_click_refine)

    def plot_data(self):
        # Clear the results tab as well
        self.results_cleared.emit()
        # Plots the data, no through update when this is changed
        # so the other data (S(Q) & D(r)) are cleared first
        _ea = np.asarray([])
        self.data['iq_x'] = _ea
        self.data['impr_iq_x'] = _ea
        self.data['sq_y'] = _ea
        self.data['dr_x'] = _ea
        self.data['dr_y'] = _ea
        self.data['int_func'] = _ea
        self.data['impr_int_func'] = _ea
        self.data['impr_dr_x'] = _ea
        self.data['impr_dr_y'] = _ea
        self.data['cor_x_cut'] = _ea
        self.data['cor_y_cut'] = _ea
        self.data['rescaled_cor_y_cut'] = _ea
        self.data['qmin'] = None
        self.data['qmax'] = None

        qmax_cut = (
            self.optim_config_widget.data_options_gb.qmax_check.isChecked()
            and self.optim_config_widget.data_options_gb.qmax_input.text()
            and self.optim_config_widget.data_options_gb.qmax_input.hasAcceptableInput()
        )
        qmin_cut = (
            self.optim_config_widget.data_options_gb.qmin_check.isChecked()
            and self.optim_config_widget.data_options_gb.qmin_input.text()
            and self.optim_config_widget.data_options_gb.qmin_input.hasAcceptableInput()
        )

        # Cut q at qmax first (if selected) and define cor_x_cut
        if qmax_cut:
            # Get q_max to cut at
            self.data['qmax'] = np.float(self.optim_config_widget.data_options_gb.qmax_input.text())
            # cut q data at qmax
            _cut = np.where(self.data['cor_x'] < self.data['qmax'])
            self.data['cor_x_cut'] = self.data['cor_x'][_cut]
            self.data['cor_y_cut'] = self.data['cor_y'][_cut]

        # Define cor_x_cut = cor_x for simpler plotting logic
        else:
            self.data['cor_x_cut'] = self.data['cor_x']
            self.data['cor_y_cut'] = self.data['cor_y']

        # Cut at qmin if selected
        if qmin_cut:
            self.data['qmin'] = np.float(self.optim_config_widget.data_options_gb.qmin_input.text())
            # Take first intensity value after q_min
            # Catch empty array error caused by no data:
            try:
                _fill_val = self.data['cor_y'][np.argmax(self.data['cor_x_cut'] > self.data['qmin'])]
                _cut = self.data['cor_y'][np.where(self.data['cor_x_cut'] > self.data['qmin'])]
                _padding = np.asarray([_fill_val] * (len(self.data['cor_x_cut']) - len(_cut)))
                self.data['cor_y_cut'] = np.concatenate((_padding, _cut))
            except ValueError:
                pass

        # this method is only used for 'raw' S(Q) actual data so mod_func = 1
        self.data['modification'] = 1
        # Smooth data if option checked. Error handling for empty array
        if self.optim_config_widget.data_options_gb.smooth_data_check.isChecked() and self.data['cor_y_cut'].size:
            self.smooth_data()
        else:
            self.optim_plot_widget.update_plots(self.data)

    def on_click_calc_sq(self):
        # Run only if data present
        if not self.data['cor_x_cut'].size:
            return
        # Run only if composition set
        _composition = self.optim_config_widget.composition_gb.get_composition_dict()
        if not _composition:
            return

        # Delete previous refined data
        self.data['impr_iq_x'] = np.asarray([])
        self.data['impr_dr_x'] = np.asarray([])
        self.data['impr_dr_y'] = np.asarray([])
        self.data['impr_int_func'] = np.asarray([])
        # Get modification function to use
        self.data['mod_func'] = self.optim_config_widget.data_options_gb.mod_func_input.currentText()
        if self.data['mod_func'] == 'Cosine-window':
            try:
                self.data['window_start'] = np.float(self.optim_config_widget.data_options_gb.window_start_input.text())
            except ValueError:
                self.window_func_error()
                return
        else:
            self.data['window_start'] = None
        # Get S(Q) method
        if self.optim_config_widget.data_options_gb.al_btn.isChecked():
            _method = 'ashcroft-langreth'
        elif self.optim_config_widget.data_options_gb.fz_btn.isChecked():
            _method = 'faber-ziman'
        self.data['sq_method'] = _method
        # Get rho 0 - Force intermediate values passed by QValidator to 0
        try:
            _rho_0 = np.float(self.optim_config_widget.composition_gb.density_input.text())
        except ValueError:
            _rho_0 = 0.0
        self.data['iq_x'] = self.data['cor_x_cut']
        self.data['sq_y'] = core.calc_structure_factor(self.data['cor_x_cut'],
                                                       self.data['cor_y_cut'],
                                                       _composition, _rho_0,
                                                       method=_method)
        _S_inf = core.calc_S_inf(_composition, self.data['cor_x_cut'], method=_method)
        self.data['int_func'] = self.data['sq_y'] - _S_inf
        self.data['dr_x'], self.data['dr_y'] = core.calc_correlation_func(self.data['iq_x'], self.data['int_func'], _rho_0,
                                                             N=self.fft_N, mod_func=self.data['mod_func'], window_start=self.data['window_start'])
        self.data['modification'] = core.get_mod_func(self.data['iq_x'], self.data['mod_func'], self.data['window_start'])
        self.optim_plot_widget.update_plots(self.data)
        self.results_cleared.emit()
        # If optimisation gb is not checked pass through Krogh-Moe/Norman S(Q)
        # to results
        if not self.optim_config_widget.optim_options_gb.isChecked():
            self.results_changed.emit()

    def on_click_refine(self):
        # Delete previous chi_sq & refined density
        try:
            del self.data['chi_sq']
            del self.data['refined_rho']
        except KeyError:
            pass
        # Don't run if sq/int_func not calculated yet
        if not self.data['int_func'].size:
            return
        # Check validity of input fields before continuing
        if not (
            self.optim_config_widget.composition_gb.density_input.hasAcceptableInput()
            and self.optim_config_widget.optim_options_gb.rmin_input.hasAcceptableInput()
            and self.optim_config_widget.optim_options_gb.niter_input.hasAcceptableInput()
        ):
            print('Error: Please ensure values are set for density, r_min, and n_iterations')
            return
        # Get modification function to use again
        self.data['mod_func'] = self.optim_config_widget.data_options_gb.mod_func_input.currentText()
        if self.data['mod_func'] == 'Cosine-window':
            try:
                self.data['window_start'] = np.float(self.optim_config_widget.data_options_gb.window_start_input.text())
            except ValueError:
                self.window_func_error()
                return
        _composition = self.optim_config_widget.composition_gb.get_composition_dict()
        # Don't run if no composition set
        if not _composition:
            return
        # Get S(Q) method
        if self.optim_config_widget.data_options_gb.al_btn.isChecked():
            _method = 'ashcroft-langreth'
        elif self.optim_config_widget.data_options_gb.fz_btn.isChecked():
            _method = 'faber-ziman'
        # Get density
        _rho_0 = np.float(self.optim_config_widget.composition_gb.density_input.text())
        # Get r_min
        _r_min = np.float(self.optim_config_widget.optim_options_gb.rmin_input.text())
        # Get no. iterations for Eggert refinement
        _n_iter = np.int(self.optim_config_widget.optim_options_gb.niter_input.text())
        if _n_iter < 2:
            print('Warning: n_iter >= 2 is recommended for convergence!')
        if (
            self.optim_config_widget.optim_options_gb.opt_check.isChecked() or
            self.optim_config_widget.optim_options_gb.opt_check_bkg.isChecked()
            ):
            # Get bounds but don't continue if none set
            try:
                if self.optim_config_widget.optim_options_gb.opt_check_bkg.isChecked():
                    # Check bkg scaling factor, uncorrected data, and bkg are present
                    if not self.data['uncorrected_y'].size or not self.data['bkg_y'].size or not self.data['bkg_scale']:
                        self.bkg_error()
                        return
                    # Get bounds on background scaling factor - return if none set
                    _lb_bkg = np.float(self.optim_config_widget.optim_options_gb.lb_input_bkg.text())
                    _ub_bkg = np.float(self.optim_config_widget.optim_options_gb.ub_input_bkg.text())
                    _bkg_scale_0 = self.data['bkg_scale']
                    opt_bkg = 1
                else:
                    opt_bkg = 0
                if self.optim_config_widget.optim_options_gb.opt_check.isChecked():
                    _lb_rho = np.float(self.optim_config_widget.optim_options_gb.lb_input.text())
                    _ub_rho = np.float(self.optim_config_widget.optim_options_gb.ub_input.text())
                    opt_rho = 1
                else:
                    opt_rho = 0
            except ValueError:
                self.limits_error()
                return
            if _n_iter > 10:
                 print('Warning: n_iter <= 10 is recommended for convergence!')
            # Only use mod_func in refinement if option is set (default is no)
            if self.mod_func_mode:
                _mod_func = self.data['mod_func']
            else:
                _mod_func = None
            # Set arguments for objective function
            # Additional arguments required if refining the background
            if opt_bkg == True:
                _smooth_flag = self.optim_config_widget.data_options_gb.smooth_data_check.isChecked()
                _args = (self.data['cor_x'], self.data['uncorrected_y'], self.data['bkg_y'],
                         self.data['qmin'], self.data['qmax'],
                         _smooth_flag, self.window_length, self.poly_order,
                         _composition, _r_min, _n_iter,
                         _method, _mod_func, self.data['window_start'], self.fft_N, opt_rho, opt_bkg)
                # Setup independent variable to pass to minimisation routine
                if opt_rho == True:
                    _minvar = np.array([_rho_0, _bkg_scale_0])
                else:
                    _minvar = _bkg_scale_0
                    # rho must be passed as argument if *only* refining the background
                    _args = (_rho_0, *_args)
            else:
                _args = (self.data['cor_x_cut'], self.data['cor_y_cut'], _composition, _r_min, _n_iter,
                         _method, _mod_func, self.data['window_start'], self.fft_N, opt_rho, opt_bkg)
                _minvar = _rho_0
            # Set up kw arguments for solver
            _solver_kwargs = {'args': _args,
                              'options': dict(self.minimisation_options),
                              'method': self.op_method}
            # Set up bounds for minimisation
            if self.op_method == 'COBYLA':
                # Construct constraints for COBYLA method. The COBYLA solver
                # does not use the 'bounds' kwarg, and requires lb & ub passed
                # as a constraint function
                if opt_rho == True:
                    _cons = [{'type': 'ineq', 'fun': lambda x: x - _lb_rho},
                             {'type': 'ineq', 'fun': lambda x: _ub_rho - x}
                             ]
                    if opt_bkg == True:
                        _cons.extend([
                                    {'type': 'ineq', 'fun': lambda x: x - _lb_bkg},
                                    {'type': 'ineq', 'fun': lambda x: _ub_bkg - x}
                                    ])
                else:
                    _cons = [{'type': 'ineq', 'fun': lambda x: x - _lb_bkg},
                             {'type': 'ineq', 'fun': lambda x: _ub_bkg - x}
                             ]
                _solver_kwargs['constraints'] = _cons
                # Handle different name of solver option ftol (tol)
                _solver_kwargs['options']['tol'] = _solver_kwargs['options'].pop('ftol')
            else:
                if opt_rho == True:
                    _solver_kwargs['bounds'] = ((_lb_rho, _ub_rho),)
                    if opt_bkg == True:
                        _solver_kwargs['bounds'] = (*_solver_kwargs['bounds'], (_lb_bkg, _ub_bkg))
                else:
                    _solver_kwargs['bounds'] = ((_lb_bkg, _ub_bkg),)

            print('\n*************************\n')
            print('Finding optimal density/background...')

            if self.global_minimisation == 1:
                # Grey out config panel while refining
                self.enable_config_panel(False)
                print('\nRunning basin-hopping algorithm to find global minimum')
                _global_min_kwargs = dict(self.global_min_options)
                _global_min_kwargs['minimizer_kwargs'] = _solver_kwargs
                _opt_result = basinhopping(core.refinement_objfun,
                                           _minvar,
                                           **_global_min_kwargs)
                # Re-enable config panel
                self.enable_config_panel(True)
            else:
                _opt_result = minimize(core.refinement_objfun,
                                       _minvar,
                                       **_solver_kwargs)
            # Extract optimised variables from _opt_result
            if opt_rho == True:
                if opt_bkg == True:
                    self.data['refined_rho'] = _opt_result.x[0]
                    _refined_bkg_scale = _opt_result.x[1]
                    print('Refined background scaling = ', _refined_bkg_scale)
                else:
                    # Handle for anomalous solver behaviour
                    # e.g. (COBYLA opt_result.x is scalar)
                    try:
                        self.data['refined_rho'] = _opt_result.x[0]
                    except IndexError:
                        self.data['refined_rho'] = _opt_result.x
                # Set optimised rho result
                _rho_temp = self.data['refined_rho']
                print('Refined density = ', _rho_temp)
                self.optim_config_widget.optim_results_gb.density_output.setText('{:.4e}'.format(self.data['refined_rho']))
                _mass_density = core.conv_density(_rho_temp, _composition)
                self.optim_config_widget.optim_results_gb.mass_density.setText('{0:.3f}'.format(_mass_density))
            else:
                # Set _rho_temp to _rho_0 if not refining rho
                _rho_temp = _rho_0
                # Handle for anomalous solver behaviour
                # e.g. (COBYLA opt_result.x is scalar)
                try:
                    _refined_bkg_scale = _opt_result.x[0]
                except IndexError:
                    _refined_bkg_scale = _opt_result.x
                print('Refined background scaling = ', _refined_bkg_scale)
            # Set optimised bkg scaling output
            if opt_bkg == True:
                self.optim_config_widget.optim_results_gb.bkg_scale_output.setText('{:.4e}'.format(_refined_bkg_scale))

            # Print Chi^2 to screen
            print('Chi^2 = ', _opt_result.fun, '\n')

            # Recalculate interference function with new background scaling and/or rho
            if opt_bkg == True:
                _I_data_corrected = self.data['uncorrected_y'] - (self.data['bkg_y'] * _refined_bkg_scale)
                # Cut data at qmax/qmin
                if self.data['qmax']:
                    _cut_max = np.where(self.data['cor_x'] < self.data['qmax'])
                    _cut_q_data = self.data['cor_x'][_cut_max]
                    _I_data_corrected = _I_data_corrected[_cut_max]
                else:
                    _cut_q_data = self.data['cor_x']
                if self.data['qmin']:
                    # Cut data to qmin - setting I(q<qmin)=I(qmin)
                    _fill_val = _I_data_corrected[np.argmax(_cut_q_data > self.data['qmin'])]
                    _cut_min = _I_data_corrected[np.where(_cut_q_data > self.data['qmin'])]
                    _padding = np.asarray([_fill_val] * (len(_cut_q_data) - len(_cut_min)))
                    _I_data_corrected = np.concatenate((_padding, _cut_min))
                # smooth data if smooth_flag == 1
                if _smooth_flag == True:
                    _I_data_corrected = data_utils.smooth_data(_I_data_corrected,
                                            window_length=self.window_length,
                                            poly_order=self.poly_order)
                _structure_factor = core.calc_structure_factor(_cut_q_data, _I_data_corrected,
                                                               _composition, _rho_temp, method=_method)
                _int_func_temp = _structure_factor - core.calc_S_inf(_composition, _cut_q_data, method=_method)
                self.data['rescaled_cor_y_cut'] = _I_data_corrected
            else:
                _structure_factor = core.calc_structure_factor(self.data['cor_x_cut'], self.data['cor_y_cut'],
                                                               _composition, _rho_temp, method=_method)
                _int_func_temp = _structure_factor - core.calc_S_inf(_composition, self.data['cor_x_cut'], method=_method)
        else:
            _rho_temp = _rho_0
            _int_func_temp = self.data['int_func']
        # Only use mod_func in refinement if option is set (default is no)
        if self.mod_func_mode:
            _mod_func = self.data['mod_func']
        else:
            _mod_func = None
        _args = (self.data['cor_x_cut'], _int_func_temp,
                 _composition, _rho_temp, _r_min, _n_iter, _method,
                 _mod_func, self.data['window_start'], self.fft_N)
        self.data['impr_int_func'], self.data['chi_sq'] = core.calc_impr_interference_func(*_args)
        self.optim_config_widget.optim_results_gb.chi_sq_output.setText('{:.4e}'.format(self.data['chi_sq']))
        # Calculated improved F_r
        self.data['impr_dr_x'], self.data['impr_dr_y'] = core.calc_correlation_func(self.data['iq_x'], self.data['impr_int_func'], _rho_temp, N=self.fft_N,
                                                                       mod_func=self.data['mod_func'], window_start=self.data['window_start'])
        self.data['impr_iq_x'] = self.data['iq_x']
        # Set modification function to None so it is not plotted this time
        _mod_func = self.data['mod_func']
        self.data['mod_func'] = 'None'
        # Plot data
        self.optim_plot_widget.update_plots(self.data)
        self.data['mod_func'] = _mod_func
        self.data['sq_method'] = _method
        # Reset rescaled I(Q) data to prevent plot persistence
        self.data['rescaled_cor_y_cut'] = np.asarray([])

        # Save refinement parameters to file
        if self.optim_config_widget.data_options_gb.smooth_data_check.isChecked():
            _smooth_bool = 'Y'
        else:
            _smooth_bool = 'N'
        if self.mod_func_mode:
            _mod_func_mode_bool = 'Y'
        else:
            _mod_func_mode_bool = 'N'
        if self.optim_config_widget.optim_options_gb.opt_check.isChecked():
            _refine_density_bool = 'Y'
            _refine_density_bounds = (
                f'Lower bound (rho) : {_lb_rho}\n'
                f'Upper bound (rho) : {_ub_rho}\n'
            )
            _refine_density_res = (
                f'Refined density : {self.data["refined_rho"]} (at/A^3)\n'
                f'{" "*18}{_mass_density} (g/cm3)\n'
            )
        else:
            _refine_density_bool = 'N'
            _refine_density_bounds = ''
            _refine_density_res = ''
        if self.optim_config_widget.optim_options_gb.opt_check_bkg.isChecked():
            _refine_bkg_bool = 'Y'
            _refine_bkg_bounds = (
                f'Lower bound (b) : {_lb_bkg}\n'
                f'Upper bound (b) : {_ub_bkg}\n'
            )
            _refine_bkg_res = (
                f'Refined background scaling : {_refined_bkg_scale}\n'
            )
        else:
            _refine_bkg_bool = 'N'
            _refine_bkg_bounds = ''
            _refine_bkg_res = ''

        if (
            self.optim_config_widget.optim_options_gb.opt_check.isChecked() or
            self.optim_config_widget.optim_options_gb.opt_check_bkg.isChecked()
            ):
            _refine_log = (
                f'Solver : {self.op_method}\n' +
                _refine_density_bounds +
                _refine_bkg_bounds +
                f'{"*"*25}\n'
                f'{_opt_result}\n'
                f'{"*"*25}\n\n' +
                _refine_bkg_res +
                _refine_density_res +
                f'Chi^2 : {self.data["chi_sq"]}\n'
            )
        else:
            _refine_log = (
                f'Chi^2 : {self.data["chi_sq"]}\n'
            )

        log_string = (
            f'refinement_log\n'
            f'{__appname__} v{__version__}\n\n'
            f'{"-"*25}\n'
            f'Data File : {self.base_filename}{self.filename_ext}\n'
            f'Composition [Element: (Z, Charge, n)]: {_composition}\n'
            f'Q_min : {self.optim_config_widget.data_options_gb.qmin_input.text()}\n'
            f'Q_max : {self.optim_config_widget.data_options_gb.qmax_input.text()}\n'
            f'Data smoothing? : {_smooth_bool}\n'
            f'Modification function : {self.data["mod_func"]}\n'
            f'Modification func used in iterative refinement? : {_mod_func_mode_bool}\n'
            f'Cosine window start : {self.data["window_start"]}\n'
            f'S(Q) formulation : {_method}\n'
            f'Density : {_rho_0}\n'
            f'r_min : {_r_min}\n'
            f'Number iterations (Eggert) : {_n_iter}\n'
            f'{"*"*25}\n'
            f'Density refined? : {_refine_density_bool}\n'
            f'Background refined? : {_refine_bkg_bool}\n' +
            _refine_log
        )

        # Append to log file or overwrite?
        if self.append_log_mode:
            _log_file = open(os.path.join(os.path.dirname(self.base_filename), 'refinement.log'), 'ab')
            # Add timestamp to top of log_string
            log_string = (f'{"#"*30}\n'
                          f'{datetime.datetime.now()}\n'
                          f'{"#"*30}\n' +
                          log_string
                          )
            np.savetxt(_log_file, [log_string], fmt='%s')
            _log_file.close()
        else:
            _log_file = self.base_filename + '_refinement.log'
            np.savetxt(_log_file, [log_string], fmt='%s')

        # Emit results changed signal to main_widget
        self.results_changed.emit()

    def smooth_data(self):
        # Update with smoothed data
        self.data['cor_y_cut'] = data_utils.smooth_data(self.data['cor_y_cut'],
                                                        window_length=self.window_length,
                                                        poly_order=self.poly_order)
        self.optim_plot_widget.update_plots(self.data)

    def enable_config_panel(self, state):
        self.optim_config_widget.setEnabled(state)
        QApplication.processEvents(QEventLoop.ExcludeUserInputEvents)

    def window_func_error(self):
        print('Please set limit for Cosine-window function')
        _message = ['Error computing S(Q)!', 'Please set limit for Cosine-window function']
        _error_msg = utility.ErrorMessageBox(_message)
        _error_msg.exec_()

    def limits_error(self):
        print('Warning: Must set bounds to refine density!')
        _message = ['Refinement error!', 'Must set bounds to refine density!']
        _error_msg = utility.ErrorMessageBox(_message)
        _error_msg.exec_()

    def bkg_error(self):
        print('Warning: Must load background file to refine background scaling factor')
        _message = ['Missing background!', 'Must load background to refine background scale factor']
        _error_msg = utility.ErrorMessageBox(_message)
        _error_msg.exec_()


class OptimConfigWidget(QWidget):

    def __init__(self):
        super(OptimConfigWidget, self).__init__()

        self.vlayout = QVBoxLayout()
        self.vlayout.setContentsMargins(0, 0, 5, 0)
        self.vlayout.setSpacing(10)

        self.composition_gb = CompositionGroupBox()
        self.data_options_gb = DataOptionsGroupBox()
        self.optim_options_gb = OptimOptionsGroupBox()
        self.optim_results_gb = OptimResultsGroupBox()

        self.vlayout.addWidget(self.composition_gb)
        self.vlayout.addWidget(self.data_options_gb)
        self.vlayout.addWidget(self.optim_options_gb)
        self.vlayout.addWidget(self.optim_results_gb)

        self.setLayout(self.vlayout)

        self.create_signals()

    def create_signals(self):
        self.optim_options_gb.toggled.connect(self._toggle_results_gb)
        self.optim_options_gb.opt_check.stateChanged.connect(self._toggle_density_refine)
        self.optim_options_gb.opt_check_bkg.stateChanged.connect(self._toggle_bkg_scale_refine)
        self.optim_results_gb.bkg_scale_copy_btn.clicked.connect(self._copy_bkg_scale_output)
        self.optim_results_gb.density_copy_btn.clicked.connect(self._copy_density_output)

    def _toggle_results_gb(self, on):
        self.optim_results_gb.setEnabled(on)

    def _toggle_density_refine(self, state):
        self.optim_results_gb.density_output.setEnabled(state)
        self.optim_results_gb.density_output_label.setEnabled(state)
        self.optim_results_gb.mass_density.setEnabled(state)
        self.optim_results_gb.mass_density_label.setEnabled(state)
        self.optim_results_gb.density_copy_btn.setEnabled(state)

    def _toggle_bkg_scale_refine(self, state):
        self.optim_results_gb.bkg_scale_output.setEnabled(state)
        self.optim_results_gb.bkg_scale_output_label.setEnabled(state)
        self.optim_results_gb.bkg_scale_copy_btn.setEnabled(state)
    
    def _copy_bkg_scale_output(self):
        # Copying of data to bkg_ui is handled by
        # gui.main_widget.table_widget.copy_bkg_scale_result
        # This function notifies the user that the bkg_ui tab has been updated
        print('Background subtraction tab updated with refined scaling factor!')
        _message = ['Background subtraction tab updated!', 'The background subtration tab has been updated with the refined scaling factor']
        _error_msg = utility.ErrorMessageBox(_message)
        _error_msg.exec_()

    def _copy_density_output(self):
        self.composition_gb.density_input.setText(self.optim_results_gb.density_output.text())


class CompositionGroupBox(QGroupBox):

    with importlib_resources.open_binary('LiquidDiffract.resources', 'pt_data.npy') as fp:
        _element_dict = np.load(fp, allow_pickle=True).item()

    def __init__(self, *args):
        super(CompositionGroupBox, self).__init__(*args)
        self.setTitle('Composition')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')

        self.create_widgets()
        self.style_widgets()
        self.create_layout()
        self.create_signals()

    def create_widgets(self):

        self.add_element_btn = QPushButton("Add")
        self.delete_element_btn = QPushButton("Delete")

        self.density_lbl = QLabel("Density:")
        self.density_input = QLineEdit("1.0")

        self.mass_density_label = QLabel('g/cm<sup>3</sup>')
        self.mass_density = QLineEdit('')

        self.composition_table = QTableWidget()

    def style_widgets(self):

        self.density_lbl.setAlignment(Qt.AlignVCenter | Qt.AlignRight)

        self.density_input.setAlignment(Qt.AlignRight)
        self.density_input.setValidator(QDoubleValidator())
        self.density_input.setMaximumWidth(100)

        self.mass_density.setAlignment(Qt.AlignRight)
        self.mass_density.setMaximumWidth(100)

        self.mass_density.setReadOnly(True)
        self.mass_density.isReadOnly()

        self.composition_table.setColumnCount(4)
        self.composition_table.horizontalHeader().setVisible(True)
        self.composition_table.verticalHeader().setVisible(False)
        self.composition_table.setContentsMargins(0, 0, 0, 0)
        self.composition_table.horizontalHeader().setStretchLastSection(False)
        self.composition_table.setHorizontalHeaderLabels(['Element', 'Z', 'Charge', 'n'])
        self.composition_table.horizontalHeaderItem(0).setToolTip('Atomic element ')
        self.composition_table.horizontalHeaderItem(1).setToolTip('Atomic number ')
        self.composition_table.horizontalHeaderItem(2).setToolTip('Charge ')
        self.composition_table.horizontalHeaderItem(3).setToolTip('Proportion of compound ')

        # Set the alignment to the headers
        self.composition_table.horizontalHeaderItem(0).setTextAlignment(Qt.AlignLeft)
        self.composition_table.horizontalHeaderItem(1).setTextAlignment(Qt.AlignHCenter)
        self.composition_table.horizontalHeaderItem(2).setTextAlignment(Qt.AlignHCenter)
        self.composition_table.horizontalHeaderItem(3).setTextAlignment(Qt.AlignHCenter)

        self.composition_table.setColumnWidth(0, 71)
        self.composition_table.setColumnWidth(1, 66)
        self.composition_table.setColumnWidth(2, 66)
        self.composition_table.setColumnWidth(3, 66)

        self.composition_table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.composition_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.composition_table.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Minimum)
        self.composition_table.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)

        self.composition_table.setItemDelegate(utility.ValidatedItemDelegate(self))

    def create_layout(self):
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(10, 10, 10, 7)
        self.main_layout.setSpacing(5)

        self.button_layout = QHBoxLayout()
        self.button_layout.setSpacing(15)
        self.button_layout.addWidget(self.add_element_btn)
        self.button_layout.addWidget(self.delete_element_btn)

        self.density_layout = QGridLayout()
        self.density_layout.addWidget(self.density_lbl, 0, 0)
        self.density_layout.addWidget(self.density_input, 0, 1)
        self.density_layout.addWidget(QLabel('atoms/A<sup>3</sup>'), 0, 2)
        self.density_layout.addWidget(self.mass_density, 1, 1)
        self.density_layout.addWidget(self.mass_density_label, 1, 2)

        self.main_layout.addLayout(self.button_layout)
        self.main_layout.addWidget(self.composition_table)
        self.main_layout.addLayout(self.density_layout)

        self.setLayout(self.main_layout)

    def create_signals(self):
        self.delete_element_btn.clicked.connect(self.delete_row)
        self.add_element_btn.clicked.connect(self.add_row)
        self.density_input.textChanged.connect(self.update_mass_density)
        self.composition_table.cellChanged.connect(self.update_mass_density)

    def add_row(self):
        _row_position = self.composition_table.rowCount()
        self.composition_table.insertRow(_row_position)
        _element_editor = QComboBox()
        _element_editor.setStyleSheet('QComboBox {border: 0px ;} ')
        for index, element in enumerate(CompositionGroupBox._element_dict):
            _element_editor.insertItem(index, element)
        # Could use a default style? e.g.:
        # _element_editor.setStyle(QStyleFactory.create('Cleanlooks'))
        # Block composition_table signals while setting the cells
        self.composition_table.blockSignals(True)
        self.composition_table.setCellWidget(_row_position, 0, _element_editor)
        _element_editor.setCurrentIndex(30)
        # No need to set default here
        # self.composition_table.setItem(_row_position , 0, QTableWidgetItem('Ga'))
        self.composition_table.setItem(_row_position, 1, QTableWidgetItem('31'))
        self.composition_table.setItem(_row_position, 2, QTableWidgetItem('0'))
        self.composition_table.setItem(_row_position, 3, QTableWidgetItem('1'))
        # Re-enable signals
        self.composition_table.blockSignals(False)

        self.composition_table.item(_row_position, 1).setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
        self.composition_table.item(_row_position, 1).setTextAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.composition_table.item(_row_position, 1).setTextAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
        self.composition_table.item(_row_position, 2).setTextAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
        self.composition_table.item(_row_position, 3).setTextAlignment(Qt.AlignHCenter | Qt.AlignVCenter)

        # Create Signal
        _element_editor.currentIndexChanged.connect(self.update_cb_val)

        self.update_mass_density()

    def delete_row(self):
        # Selects last row if none selected by user
        _selected_row = self.composition_table.currentRow()
        if _selected_row == -1:
            _selected_row = self.composition_table.rowCount() - 1
        else:
            pass
        self.composition_table.removeRow(_selected_row)
        self.update_mass_density()

    def update_cb_val(self):
        _cb_widget = self.sender()
        _current_row = self.composition_table.indexAt(_cb_widget.pos()).row()
        _new_element = str(_cb_widget.currentText())
        _new_Z_val = str(CompositionGroupBox._element_dict[_new_element])
        self.composition_table.item(_current_row, 1).setText(_new_Z_val)
        self.update_mass_density()

    def get_composition_dict(self):
        '''Return composition dictionary'''
        # Form of dictionary is col 0 is key, val is tuple(col1,col2,col3)
        _composition_dict = {}
        _key_list = list(CompositionGroupBox._element_dict.keys())
        _val_list = list(CompositionGroupBox._element_dict.values())
        for _row_index in range(self.composition_table.rowCount()):
            _Z = np.int(self.composition_table.item(_row_index, 1).text())
            _charge = np.int(self.composition_table.item(_row_index, 2).text())
            _n = np.int(self.composition_table.item(_row_index, 3).text())
            _key = str(_key_list[_val_list.index(_Z)])
            _dict_entry = {_key: [_Z, _charge, _n]}
            _composition_dict.update(_dict_entry)

        # This code snippet convert _n from integer number of atoms in the
        # formula unit to a fraction of the total - e.g. SiO2 from Si=1 >> 1/3
        # _n_total = sum([_composition_dict[_el][2] for _el in _composition_dict])
        # for _el in _composition_dict:
        #     _composition_dict[_el][2] /= _n_total
        #     _composition_dict[_el] = tuple(_composition_dict[_el])
        # print(_composition_dict)
        return _composition_dict

    def update_mass_density(self):
        _composition = self.get_composition_dict()
        # Handle errors caused by QDoubleValidator allowing intermediate
        # values to pass e.g. '.', ' ' etc.
        try:
            _atomic_density = np.float(self.density_input.text())
        except ValueError:
            _atomic_density = 0.0
        _mass_density = core.conv_density(_atomic_density, _composition)
        self.mass_density.setText('{0:.3f}'.format(_mass_density))


class DataOptionsGroupBox(QGroupBox):

    def __init__(self, *args):
        super(DataOptionsGroupBox, self).__init__(*args)
        self.setTitle('Data Options')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')

        self.create_widgets()
        self.style_widgets()
        self.create_layout()
        self.create_signals()

    def create_widgets(self):

        self.qmax_label = QLabel('Q-Max cutoff: ')
        self.qmax_input = QLineEdit()
        self.qmax_check = QCheckBox()

        self.qmin_label = QLabel('Q-Min cutoff: ')
        self.qmin_input = QLineEdit()
        self.qmin_check = QCheckBox()

        self.smooth_label = QLabel('Smooth Data? ')
        self.smooth_data_check = QCheckBox()

        self.mod_func_lbl = QLabel('Use modification function?')
        self.mod_func_input = QComboBox()
        self.mod_func_input.insertItem(0, 'None')
        self.mod_func_input.insertItem(1, 'Lorch')
        self.mod_func_input.insertItem(2, 'Cosine-window')
        self.mod_func_input.setCurrentIndex(0)

        self.window_start_input = QLineEdit()

        self.method_button_group = QButtonGroup()
        self.al_btn = QRadioButton('Ashcroft-Langreth')
        self.fz_btn = QRadioButton('Faber-Ziman')
        self.method_button_group.addButton(self.fz_btn)
        self.method_button_group.addButton(self.al_btn)

        self.method_lbl = QLabel('S(Q) formulation: ')

        self.calc_sq_btn = QPushButton('Calc S(Q)')

    def style_widgets(self):

        self.qmax_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.qmax_input.setAlignment(Qt.AlignRight)
        self.qmax_input.setValidator(QDoubleValidator())
        self.qmax_input.setMaximumWidth(70)
        self.qmax_input.setEnabled(False)
        self.qmax_check.setChecked(False)
        self.qmax_label.setEnabled(False)

        self.qmin_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.qmin_input.setAlignment(Qt.AlignRight)
        self.qmin_input.setValidator(QDoubleValidator())
        self.qmin_input.setMaximumWidth(70)
        self.qmin_input.setEnabled(False)
        self.qmin_check.setChecked(False)
        self.qmin_label.setEnabled(False)

        self.smooth_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.smooth_data_check.setChecked(False)

        self.window_start_input.setValidator(QDoubleValidator())
        self.window_start_input.setEnabled(False)
        self.fz_btn.setChecked(True)

    def create_layout(self):
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(20, 10, 20, 7)
        self.main_layout.setSpacing(5)

        self.grid_layout = QGridLayout()
        self.grid_layout.setSpacing(15)
        self.grid_layout.setColumnStretch(0, 4)
        self.grid_layout.setColumnStretch(1, 4)
        self.grid_layout.setColumnStretch(2, 2)

        self.grid_layout.addWidget(self.qmin_label, 0, 0)
        self.grid_layout.addWidget(self.qmin_input, 0, 1)
        self.grid_layout.addWidget(self.qmin_check, 0, 2)

        self.grid_layout.addWidget(self.qmax_label, 1, 0)
        self.grid_layout.addWidget(self.qmax_input, 1, 1)
        self.grid_layout.addWidget(self.qmax_check, 1, 2)

        self.grid_layout.addWidget(self.smooth_label, 2, 0)
        self.grid_layout.addWidget(self.smooth_data_check, 2, 2)

        self.grid_layout.addWidget(self.mod_func_lbl, 3, 0, 1, 2)
        self.grid_layout.addWidget(self.mod_func_input, 4, 0, 1, 2)
        self.grid_layout.addWidget(self.window_start_input, 4, 2)
        self.grid_layout.addWidget(self.method_lbl, 5, 0)

        self.hbtn_layout = QHBoxLayout()
        self.hbtn_layout.addWidget(self.fz_btn)
        self.hbtn_layout.addWidget(self.al_btn)

        self.main_layout.addLayout(self.grid_layout)
        self.main_layout.addLayout(self.hbtn_layout)
        self.main_layout.addWidget(self.calc_sq_btn)

        self.setLayout(self.main_layout)

    def create_signals(self):
        self.qmax_check.stateChanged.connect(self.qmax_state_changed)
        self.qmin_check.stateChanged.connect(self.qmin_state_changed)
        self.mod_func_input.currentIndexChanged.connect(self.mod_func_changed)

    def qmax_state_changed(self):
        if self.qmax_check.isChecked():
            self.qmax_input.setEnabled(True)
            self.qmax_label.setEnabled(True)
        else:
            self.qmax_input.setEnabled(False)
            self.qmax_label.setEnabled(False)

    def qmin_state_changed(self):
        if self.qmin_check.isChecked():
            self.qmin_input.setEnabled(True)
            self.qmin_label.setEnabled(True)
        else:
            self.qmin_input.setEnabled(False)
            self.qmin_label.setEnabled(False)

    def mod_func_changed(self):
        if self.mod_func_input.currentText() == 'Cosine-window':
            self.window_start_input.setEnabled(True)
        else:
            self.window_start_input.setEnabled(False)


class OptimOptionsGroupBox(QGroupBox):

    def __init__(self, *args):
        super(OptimOptionsGroupBox, self).__init__(*args)
        self.setTitle('Optimisation Options')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        self.setCheckable(True)
        self.setChecked(True)

        self.create_widgets()
        self.style_widgets()
        self.create_layout()
        self.create_signals()

    def create_widgets(self):
        self.rmin_label = QLabel('R-Min cutoff: ')
        self.rmin_input = QLineEdit('2.3')
        self.niter_label = QLabel('No. iterations: ')
        self.niter_input = QLineEdit('5')
        self.opt_check = QCheckBox('Refine density? ')
        self.lb_label = QLabel('Lower bound: ')
        self.lb_input = QLineEdit()
        self.ub_label = QLabel('Upper bound: ')
        self.ub_input = QLineEdit()
        self.opt_check_bkg = QCheckBox('Refine background scale factor?')
        self.lb_label_bkg = QLabel('Lower bound: ')
        self.lb_input_bkg = QLineEdit()
        self.ub_label_bkg = QLabel('Upper bound: ')
        self.ub_input_bkg = QLineEdit()
        self.opt_button = QPushButton('Refine S(Q)')

    def style_widgets(self):

        self.rmin_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.niter_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.lb_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.ub_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.lb_label_bkg.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.ub_label_bkg.setAlignment(Qt.AlignVCenter | Qt.AlignRight)

        self.rmin_input.setAlignment(Qt.AlignRight)
        self.rmin_input.setValidator(QDoubleValidator())
        self.rmin_input.setMaximumWidth(70)
        self.rmin_input.setToolTip('Intramolecular distance cut-off')
        self.rmin_label.setToolTip('Intramolecular distance cut-off')

        self.niter_input.setAlignment(Qt.AlignRight)
        self.niter_input.setValidator(QIntValidator())
        self.niter_input.setMaximumWidth(70)

        self.lb_input.setAlignment(Qt.AlignRight)
        self.lb_input.setValidator(QDoubleValidator())
        self.lb_input.setMaximumWidth(70)

        self.ub_input.setAlignment(Qt.AlignRight)
        self.ub_input.setValidator(QDoubleValidator())
        self.ub_input.setMaximumWidth(70)

        self.lb_input_bkg.setAlignment(Qt.AlignRight)
        self.lb_input_bkg.setValidator(QDoubleValidator())
        self.lb_input_bkg.setMaximumWidth(70)

        self.ub_input_bkg.setAlignment(Qt.AlignRight)
        self.ub_input_bkg.setValidator(QDoubleValidator())
        self.ub_input_bkg.setMaximumWidth(70)

        self.opt_check.setChecked(False)
        self.opt_check_bkg.setChecked(False)

        self.lb_label.setEnabled(False)
        self.lb_input.setEnabled(False)
        self.ub_label.setEnabled(False)
        self.ub_input.setEnabled(False)
        self.lb_label_bkg.setEnabled(False)
        self.lb_input_bkg.setEnabled(False)
        self.ub_label_bkg.setEnabled(False)
        self.ub_input_bkg.setEnabled(False)

    def create_layout(self):

        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(20, 10, 20, 7)
        self.main_layout.setSpacing(25)

        self.top_grid = QGridLayout()
        self.top_grid.setSpacing(15)
        # Can manually set column stretch:
        # self.grid_layout.setColumnStretch(0, 4)
        # self.grid_layout.setColumnStretch(1, 4)
        # self.grid_layout.setColumnStretch(2, 2)
        self.top_grid.addWidget(self.rmin_label, 0, 0)
        self.top_grid.addWidget(self.rmin_input, 0, 1)
        self.top_grid.addWidget(self.niter_label, 1, 0)
        self.top_grid.addWidget(self.niter_input, 1, 1)

        self.bottom_grid = QGridLayout()
        self.bottom_grid.setSpacing(15)
        self.bottom_grid.addWidget(self.opt_check, 0, 0)
        self.bottom_grid.addWidget(self.lb_label, 1, 0)
        self.bottom_grid.addWidget(self.lb_input, 1, 1)
        self.bottom_grid.addWidget(self.ub_label, 2, 0)
        self.bottom_grid.addWidget(self.ub_input, 2, 1)

        self.bottom_grid.addWidget(self.opt_check_bkg, 3, 0)
        self.bottom_grid.addWidget(self.lb_label_bkg, 4, 0)
        self.bottom_grid.addWidget(self.lb_input_bkg, 4, 1)
        self.bottom_grid.addWidget(self.ub_label_bkg, 5, 0)
        self.bottom_grid.addWidget(self.ub_input_bkg, 5, 1)

        self.main_layout.addLayout(self.top_grid)
        self.main_layout.addLayout(self.bottom_grid)
        self.main_layout.addWidget(self.opt_button)

        self.setLayout(self.main_layout)

    def create_signals(self):
        self.opt_check.stateChanged.connect(self.opt_state_changed)
        self.opt_check_bkg.stateChanged.connect(self.opt_bkg_state_changed)

    def opt_state_changed(self):
        if self.opt_check.isChecked():
            self.lb_input.setEnabled(True)
            self.lb_label.setEnabled(True)
            self.ub_input.setEnabled(True)
            self.ub_label.setEnabled(True)
        else:
            self.lb_input.setEnabled(False)
            self.lb_label.setEnabled(False)
            self.ub_input.setEnabled(False)
            self.ub_label.setEnabled(False)

    def opt_bkg_state_changed(self):
        if self.opt_check_bkg.isChecked():
            self.lb_input_bkg.setEnabled(True)
            self.lb_label_bkg.setEnabled(True)
            self.ub_input_bkg.setEnabled(True)
            self.ub_label_bkg.setEnabled(True)
        else:
            self.lb_input_bkg.setEnabled(False)
            self.lb_label_bkg.setEnabled(False)
            self.ub_input_bkg.setEnabled(False)
            self.ub_label_bkg.setEnabled(False)

class OptimResultsGroupBox(QGroupBox):

    def __init__(self, *args):
        super(OptimResultsGroupBox, self).__init__(*args)
        self.setTitle('Results')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        self.setEnabled(True)

        self.create_widgets()
        self.style_widgets()
        self.create_layout()

    def create_widgets(self):

        self.chi_sq_label = QLabel('Final Chi-squared: ')
        self.chi_sq_output = QLineEdit()
        self.bkg_scale_output_label = QLabel('Refined background scale factor: ')
        self.bkg_scale_output = QLineEdit()
        self.bkg_scale_copy_btn = QPushButton()
        self.density_output_label = QLabel('Refined density (at/A<sup>3</sup>): ')
        self.density_output = QLineEdit()        
        self.density_copy_btn = QPushButton()
        self.mass_density_label = QLabel('(g/cm<sup>3</sup>): ')
        self.mass_density = QLineEdit()

    def style_widgets(self):

        self.chi_sq_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.bkg_scale_output_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.density_output_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.mass_density_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.bkg_scale_copy_btn.setIcon(self.style().standardIcon(getattr(QStyle, 'SP_BrowserReload')))
        self.density_copy_btn.setIcon(self.style().standardIcon(getattr(QStyle, 'SP_BrowserReload')))

        self.chi_sq_output.isReadOnly()
        self.bkg_scale_output.isReadOnly()
        self.density_output.isReadOnly()
        self.mass_density.isReadOnly()
        self.chi_sq_output.setMaximumWidth(82)
        self.bkg_scale_output.setMaximumWidth(82)
        self.density_output.setMaximumWidth(82)
        self.mass_density.setMaximumWidth(82)
        self.bkg_scale_copy_btn.setMaximumWidth(30)
        self.density_copy_btn.setMaximumWidth(30)

        self.bkg_scale_output.setEnabled(False)
        self.bkg_scale_output_label.setEnabled(False)
        self.density_output.setEnabled(False)
        self.density_output_label.setEnabled(False)
        self.mass_density.setEnabled(False)
        self.mass_density_label.setEnabled(False)
        self.bkg_scale_copy_btn.setEnabled(False)
        self.density_copy_btn.setEnabled(False)

        self.bkg_scale_output.setReadOnly(True)
        self.density_output.setReadOnly(True)
        self.chi_sq_output.setReadOnly(True)
        self.mass_density.setReadOnly(True)

    def create_layout(self):

        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(20, 10, 20, 7)
        self.main_layout.setSpacing(25)

        self.grid_layout = QGridLayout()
        self.grid_layout.setSpacing(15)

        self.grid_layout.addWidget(self.chi_sq_label, 0, 1)
        self.grid_layout.addWidget(self.chi_sq_output, 0, 2)
        self.grid_layout.addWidget(self.bkg_scale_copy_btn, 1, 0)
        self.grid_layout.addWidget(self.bkg_scale_output_label, 1, 1)
        self.grid_layout.addWidget(self.bkg_scale_output, 1, 2)
        self.grid_layout.addWidget(self.density_copy_btn, 2, 0)
        self.grid_layout.addWidget(self.density_output_label, 2, 1)
        self.grid_layout.addWidget(self.density_output, 2, 2)
        self.grid_layout.addWidget(self.mass_density_label, 3, 1)
        self.grid_layout.addWidget(self.mass_density, 3, 2)

        self.main_layout.addLayout(self.grid_layout)

        self.setLayout(self.main_layout)
