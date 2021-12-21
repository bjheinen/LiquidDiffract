# -*- coding: utf-8 -*-
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2021, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

import os.path
import datetime
import numpy as np
import lmfit
from PyQt5.QtCore import Qt, pyqtSignal, QSignalMapper, QObject
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QWidget, QFrame, QGridLayout, QVBoxLayout, \
                            QHBoxLayout, QGroupBox, QPushButton, QRadioButton, \
                            QLineEdit, QDoubleSpinBox, QLabel, QScrollArea, \
                            QCheckBox, QSplitter, QComboBox, QStyle
# LiquidDiffract imports
from LiquidDiffract.gui import plot_widgets
from LiquidDiffract.core import data_utils
from LiquidDiffract.core import peak_fit
import LiquidDiffract.core.core as core
from LiquidDiffract.version import __appname__, __version__


class StructureUI(QWidget):

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

        # Initialise data
        self.active_func = 'rdf'
        self.clear_data()

        self.create_signals()

    def create_signals(self):

        self.structure_config_widget.plot_view_gb.rdf_btn.toggled.connect(self.toggle_plot_view)

        self.structure_config_widget.monatomic_gb.r0_input.valueChanged.connect(self.update_plot_data)
        self.structure_config_widget.monatomic_gb.rpmax_input.valueChanged.connect(self.update_plot_data)
        self.structure_config_widget.monatomic_gb.rmax_input.valueChanged.connect(self.update_plot_data)
        self.structure_config_widget.monatomic_gb.rmin_input.valueChanged.connect(self.update_plot_data)

        self.structure_config_widget.monatomic_gb.r0_input.valueChanged.connect(self.calc_integrals)
        self.structure_config_widget.monatomic_gb.rpmax_input.valueChanged.connect(self.calc_integrals)
        self.structure_config_widget.monatomic_gb.rmax_input.valueChanged.connect(self.calc_integrals)
        self.structure_config_widget.monatomic_gb.rmin_input.valueChanged.connect(self.calc_integrals)

        self.int_limit_signalMapper = QSignalMapper()
        self.int_limit_signalMapper.mapped[QObject].connect(self.update_int_limits)

        self.structure_config_widget.monatomic_gb.find_limits_btn.clicked.connect(self.auto_integration_limits)
        self.structure_config_widget.monatomic_gb.calc_N_btn.clicked.connect(self.calc_integrals)

        self.structure_config_widget.monatomic_gb.toggled.connect(self.toggle_monatomic_gb)
        self.structure_config_widget.polyatomic_gb.toggled.connect(self.toggle_polyatomic_gb)

        self.structure_plot_widget.fit_limits.sigRegionChanged.connect(self.update_fit_limits_sb)
        self.structure_config_widget.polyatomic_gb.min_limit_input.valueChanged.connect(self.update_fit_limits_plot)
        self.structure_config_widget.polyatomic_gb.max_limit_input.valueChanged.connect(self.update_fit_limits_plot)

        self.structure_config_widget.polyatomic_gb.fit_peaks_btn.clicked.connect(self.on_click_fit_peaks)
        self.structure_config_widget.polyatomic_gb.guess_peak.connect(self.guess_peak_params)
        self.structure_config_widget.polyatomic_gb.gauss_params_changed_relay.connect(self.make_peak_plots)

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
                     'r0': 0.0, 'rmax': 0.0, 'rpmax': 0.0, 'rmin': 0.0,
                     'fit_r': np.asarray([]),
                     'obj_fun': np.asarray([]),
                     'gauss_peaks': [],
                     'gauss_model': np.asarray([]),
                     'gauss_residuals': np.asarray([]),
                     'gauss_residuals_full': np.asarray([]),
                     'rho': None, 'composition': None, 'sq_x': np.asarray([])}
        self.weights = {}
        self.c_dict = {}

    def set_atoms(self):

        _atom_list = self.data['composition'].keys()
        if self.structure_config_widget.polyatomic_gb.atom_list == _atom_list:
            pass

        else:
            self.structure_config_widget.polyatomic_gb.atom_list = _atom_list
            for _peak in self.structure_config_widget.polyatomic_gb.peak_dict.values():
                _peak.populate_atom_list(_atom_list)

    def plot_data(self, _reset_view=1):
        self.update_plot_data()
        self.structure_plot_widget.update_plot_windows(self.data, _reset_view)

    def update_plot_data(self):
        self.data['r0'] = self.structure_config_widget.monatomic_gb.r0_input.value()
        self.data['rpmax'] = self.structure_config_widget.monatomic_gb.rpmax_input.value()
        self.data['rmax'] = self.structure_config_widget.monatomic_gb.rmax_input.value()
        self.data['rmin'] = self.structure_config_widget.monatomic_gb.rmin_input.value()

        if not self.data['rdf_y'].size:
            return

        self.structure_plot_widget.update_plots(self.data)
        self.create_int_limit_signals()

    def update_obj_fun(self):
        if self.active_func == 'rdf':
            self.data['obj_fun'] = self.data['rdf_y']
        elif self.active_func == 'tr':
            self.data['obj_fun'] = self.data['tr_y']

    def toggle_plot_view(self):
        # Change active view
        if self.structure_config_widget.plot_view_gb.rdf_btn.isChecked():
            self.active_func = 'rdf'
            # Set label on curve fitting plot
            self.structure_plot_widget.fit_plot.setLabel('left', text='RDF(r)')
        elif self.structure_config_widget.plot_view_gb.tr_btn.isChecked():
            self.active_func = 'tr'
            self.structure_plot_widget.fit_plot.setLabel('left', text='T(r)')
        self.update_obj_fun()

        # Toggle plots if monatomic gb selected
        if self.structure_config_widget.monatomic_gb.isChecked():
            if self.structure_config_widget.plot_view_gb.rdf_btn.isChecked():
                self.structure_plot_widget.fit_layout_widget.setVisible(False)
                self.structure_plot_widget.tr_int_layout_widget.setVisible(False)
                self.structure_plot_widget.rdf_int_layout_widget.setVisible(True)
            elif self.structure_config_widget.plot_view_gb.tr_btn.isChecked():
                self.structure_plot_widget.fit_layout_widget.setVisible(False)
                self.structure_plot_widget.rdf_int_layout_widget.setVisible(False)
                self.structure_plot_widget.tr_int_layout_widget.setVisible(True)

        elif self.structure_config_widget.polyatomic_gb.isChecked():
            self.structure_plot_widget.fit_layout_widget.setVisible(True)
            self.structure_plot_widget.tr_int_layout_widget.setVisible(False)
            self.structure_plot_widget.rdf_int_layout_widget.setVisible(False)
            self.make_peak_plots()


    def toggle_monatomic_gb(self):
        if self.structure_config_widget.monatomic_gb.isChecked():
            # Toggled to enable monatomic (integration) gb
            # Disable gaussian fitting gb
            self.structure_config_widget.polyatomic_gb.setChecked(False)
            self.toggle_plot_view()

    def toggle_polyatomic_gb(self):
        if self.structure_config_widget.polyatomic_gb.isChecked():
            # Toggled to enable polyatomic (gaussian fitting) gb
            # Disable monatomic (integration) gb
            self.structure_config_widget.monatomic_gb.setChecked(False)
            self.toggle_plot_view()

    def update_fit_limits_sb(self):
        _min, _max = self.structure_plot_widget.fit_limits.getRegion()
        self.structure_config_widget.polyatomic_gb.min_limit_input.setValue(_min)
        self.structure_config_widget.polyatomic_gb.max_limit_input.setValue(_max)

    def update_fit_limits_plot(self):
        _min = self.structure_config_widget.polyatomic_gb.min_limit_input.value()
        _max = self.structure_config_widget.polyatomic_gb.max_limit_input.value()
        self.structure_plot_widget.fit_limits.setRegion((_min, _max))

    def auto_integration_limits(self):
        if not self.data['rdf_y'].size:
            return
        _r_0, _rp_max, _r_max, _r_min = data_utils.find_integration_limits(self.data['rdf_x'], self.data['rdf_y'])
        self.structure_config_widget.monatomic_gb.r0_input.setValue(_r_0)
        self.structure_config_widget.monatomic_gb.rpmax_input.setValue(_rp_max)
        self.structure_config_widget.monatomic_gb.rmax_input.setValue(_r_max)
        self.structure_config_widget.monatomic_gb.rmin_input.setValue(_r_min)

    def calc_integrals(self):
        if not self.data['rdf_y'].size:
            return
        _Na, _Nb, _Nc = core.integrate_coordination_sphere(self.data['rdf_x'], self.data['rdf_y'],
                                                           r_0=self.data['r0'], rp_max=self.data['rpmax'],
                                                           r_max=self.data['rmax'], r_min=self.data['rmin'],
                                                           method=0)
        self.structure_config_widget.monatomic_gb.Na_output.setText('{0:.3f}'.format(_Na))
        self.structure_config_widget.monatomic_gb.Nb_output.setText('{0:.3f}'.format(_Nb))
        self.structure_config_widget.monatomic_gb.Nc_output.setText('{0:.3f}'.format(_Nc))
        self.data['N_a'] = _Na
        self.data['N_b'] = _Nb
        self.data['N_c'] = _Nc

    def set_weights(self):
        if self.xray_weight_mode == 1:
            self.weights, self.c_dict = core.calculate_weights(self.data['composition'], self.data['sq_x'])
        elif self.xray_weight_mode == 0:
            self.weights, self.c_dict = core.calculate_weights(self.data['composition'], [0])
        else:
            raise ValueError

    # TODO - Proper guess peak implement plain gauss!
    def guess_peak_params(self, _idx):
        _peak_dict = self.structure_config_widget.polyatomic_gb.peak_dict
        _peak_widget = _peak_dict[_idx]
        if self.data['rdf_y'].size:
            # When data is present, guess the peak parameters
            # Not yet implemented - currently just set r to middle of
            # fitting range
            _lb, _ub = self.structure_plot_widget.fit_limits.getRegion()
            _lb = round(_lb, 1)
            _ub = round(_ub, 1)
            _pos = round((_lb + _ub)/2, 1)
            _peak_widget.r_input.setText(str(_pos))
            _peak_widget.r_lb.setText(str(_lb))
            _peak_widget.r_ub.setText(str(_ub))

        else:
            # If not data present set default parameters here
            # Currently default parameters are set for N/s in groupbox
            # Here we just set r to 3.0 as a default
            _peak_widget.r_input.setText('3.0')
            _peak_widget.r_lb.setText('1.0')
            _peak_widget.r_ub.setText('4.0')

    def make_peak_params(self):
        # Clear peak_fit_dict each time
        self.peak_fit_dict = {}
        self.fit_result = None

        _peak_dict = self.structure_config_widget.polyatomic_gb.peak_dict
        self.peak_dict = _peak_dict
        _peak_idx_list = list(_peak_dict.keys())
        _peak_idx_remove = []
        for _peak_idx in _peak_idx_list:
            _peak_idx_str = str(_peak_idx)
            _peak_widget = _peak_dict[_peak_idx]
            _alpha = _peak_widget.alpha_input.currentText()
            _beta = _peak_widget.beta_input.currentText()

            try:
                _N = np.float(_peak_widget.N_input.text())
            except ValueError:
                _N = 0.0
            _N_lb = np.float(_peak_widget.N_lb.text())
            _N_ub = np.float(_peak_widget.N_ub.text())
            _N_refine = _peak_widget.N_refine.isChecked()

            try:
                _r = np.float(_peak_widget.r_input.text())
            except ValueError:
                _r = 0.0
            _r_lb = np.float(_peak_widget.r_lb.text())
            _r_ub = np.float(_peak_widget.r_ub.text())
            _r_refine = _peak_widget.r_refine.isChecked()

            try:
                _s = np.float(_peak_widget.s_input.text())
            except ValueError:
                _s = 0.0
            _s_lb = np.float(_peak_widget.s_lb.text())
            _s_ub = np.float(_peak_widget.s_ub.text())
            _s_refine = _peak_widget.s_refine.isChecked()

            if _peak_widget.skew_toggle.isChecked():
                try:
                    _xi = np.float(_peak_widget.xi_input.text())
                except ValueError:
                    _xi = 0.0
                _xi_lb = np.float(_peak_widget.xi_lb.text())
                _xi_ub = np.float(_peak_widget.xi_ub.text())
                _xi_refine = _peak_widget.xi_refine.isChecked()

            else:
                _xi = 0.0
                _xi_lb = 0.0
                _xi_ub = 1.0
                _xi_refine = False

            # Make peak name
            _name = '#' + _peak_idx_str + ' ' + _alpha + '-' + _beta

            try:
                # Make kwarg dict
                _W = self.weights[(_alpha,_beta)]
                _cb = self.c_dict[_beta]
                _kwargs = {'W'+_peak_idx_str: _W, 'c_b'+_peak_idx_str: _cb}
                # Make parameter object
                _params = lmfit.Parameters()
                # Form for lmfit parameters
                # (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
                _params.add_many(
                    ('N'+_peak_idx_str, _N, _N_refine, _N_lb, _N_ub, None, None),
                    ('r'+_peak_idx_str, _r, _r_refine, _r_lb, _r_ub, None, None),
                    ('s'+_peak_idx_str, _s, _s_refine, _s_lb, _s_ub, None, None),
                    ('xi'+_peak_idx_str, _xi, _xi_refine, _xi_lb, _xi_ub, None, None))
                # Store name, params, kwargs in peak_fit_dict
                self.peak_fit_dict['name'+_peak_idx_str] = _name
                self.peak_fit_dict['kwargs'+_peak_idx_str] = _kwargs
                self.peak_fit_dict['params'+_peak_idx_str] = _params

            except KeyError:
                # If alpha/beta not set remove peak from the list
                _peak_idx_remove.append(_peak_idx)

        # Remove peaks with alpha-beta not set
        for _peak_to_remove in _peak_idx_remove:
            _peak_idx_list.remove(_peak_to_remove)
        # Store peak_idx_list in peak_fit_dict
        self.peak_fit_dict['peak_idx_list'] = _peak_idx_list

        # Combine param objects
        _obj_params = lmfit.Parameters()
        for _peak_idx in _peak_idx_list:
            _peak_idx_str = str(_peak_idx)
            _peak_params = self.peak_fit_dict['params'+_peak_idx_str]
            for _param in list(_peak_params):
                _obj_params.add(_peak_params[_param])
        # Store combined parameters for objective function
        self.peak_fit_dict['obj_params'] = _obj_params

        # Combine kwarg dicts for objective function
        _obj_kwargs = {}
        for _peak_idx in _peak_idx_list:
            _peak_idx_str = str(_peak_idx)
            _obj_kwargs.update(self.peak_fit_dict['kwargs'+_peak_idx_str])
            # Stor peak idx (as peak_idx_list) in each kwarg dict
            self.peak_fit_dict['kwargs'+_peak_idx_str]['peak_idx_list'] = [_peak_idx]
            self.peak_fit_dict['kwargs'+_peak_idx_str]['active_func'] = self.active_func

        # Store combined kwargs in peak_fit_dict
        self.peak_fit_dict['obj_kwargs'] = _obj_kwargs
        # Store peak_idx_list in objective kwargs also
        self.peak_fit_dict['obj_kwargs']['peak_idx_list'] = _peak_idx_list
        # And the active obj func
        self.peak_fit_dict['obj_kwargs']['active_func'] = self.active_func

        # Cut data at the fit range
        _fit_min = self.structure_config_widget.polyatomic_gb.min_limit_input.value()
        _fit_max = self.structure_config_widget.polyatomic_gb.max_limit_input.value()
        _fit_range = np.where((self.data['rdf_x'] >= _fit_min) & (self.data['rdf_x'] <= _fit_max))
        self.data['fit_r'] = self.data['rdf_x'][_fit_range]
        self.data['fit_y'] = self.data['obj_fun'][_fit_range]

    def on_click_fit_peaks(self):

        self.make_peak_params()
        # Exit if no peaks present
        if not self.peak_fit_dict['peak_idx_list']:
            return

        # Double check there are no NaNs in data-set - can occur if x-lim is 0
        if np.isnan(self.data['fit_y']).any():
            self.data['fit_y'] = data_utils.interp_nan(self.data['fit_y'])

        self.fit_result = lmfit.minimize(peak_fit.gauss_obj_func,
                                         self.peak_fit_dict['obj_params'],
                                         args=(self.data['fit_r'], self.data['fit_y']),
                                         kws=self.peak_fit_dict['obj_kwargs'])

        # Save peak fit results to file
        # Make log of peak IDs/atom-pairs
        _peak_log = (f'Peak IDs\n{"-"*8}\n')
        for _peak_idx in self.peak_fit_dict['peak_idx_list']:
            _peak_id_str = f'{_peak_idx} :  {self.peak_fit_dict["name"+str(_peak_idx)]}\n'
            _peak_log += _peak_id_str
        # Make log string
        log_string = (
            f'refinement_log (fit peaks)\n'
            f'{__appname__} v{__version__}\n\n'
            f'{"-"*25}\n'
            f'Data File : {self.base_filename}{self.filename_ext}\n'
            f'Composition [Element: (Z, Charge, n)]: {self.data["composition"]}\n'
            f'{_peak_log}'
            f'{"*"*25}\n'
            f'Using lmfit v{lmfit.__version__}\n'
            f'Peak fit success?: {self.fit_result.success} | {self.fit_result.message}\n'
            f'{lmfit.fit_report(self.fit_result)}\n'
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

        # Set refined peak parameters
        self.set_peak_params()

    def set_peak_params(self):

        # Get parameter object
        _params = self.fit_result.params
        # can loop through peak list
        _peak_idx_list = self.peak_fit_dict['peak_idx_list']
        for _peak_idx in _peak_idx_list:
            _peak_idx_str = str(_peak_idx)
            _peak_widget = self.peak_dict[_peak_idx]
            # Reset alpha and beta
            _name = self.peak_fit_dict['name'+_peak_idx_str]
            _alpha, _beta = _name.split(' ')[-1].split('-')
            _peak_widget.alpha_input.setCurrentText(_alpha)
            _peak_widget.beta_input.setCurrentText(_beta)
            # Set N
            _N_param = _params['N'+_peak_idx_str]
            _peak_widget.N_input.setText('{0:.2f}'.format(_N_param.value))
            _peak_widget.N_lb.setText('{0:.2f}'.format(_N_param.min))
            _peak_widget.N_ub.setText('{0:.2f}'.format(_N_param.max))
            _peak_widget.N_refine.setChecked(_N_param.vary)

            # Set r
            _r_param = _params['r'+_peak_idx_str]
            _peak_widget.r_input.setText('{0:.2f}'.format(_r_param.value))
            _peak_widget.r_lb.setText('{0:.2f}'.format(_r_param.min))
            _peak_widget.r_ub.setText('{0:.2f}'.format(_r_param.max))
            _peak_widget.r_refine.setChecked(_r_param.vary)

            # Set sigma
            _s_param = _params['s'+_peak_idx_str]
            _peak_widget.s_input.setText('{0:.2f}'.format(_s_param.value))
            _peak_widget.s_lb.setText('{0:.2f}'.format(_s_param.min))
            _peak_widget.s_ub.setText('{0:.2f}'.format(_s_param.max))
            _peak_widget.s_refine.setChecked(_s_param.vary)

            # Set xi (skewness) only if in use
            if _peak_widget.skew_toggle.isChecked():
                _xi_param = _params['xi'+_peak_idx_str]
                _peak_widget.xi_input.setText('{0:.2f}'.format(_xi_param.value))
                _peak_widget.xi_lb.setText('{0:.2f}'.format(_xi_param.min))
                _peak_widget.xi_ub.setText('{0:.2f}'.format(_xi_param.max))
                _peak_widget.xi_refine.setChecked(_xi_param.vary)

        # Clear the peak fit dict
        self.peak_fit_dict = {}
        # Clear fitting results
        self.fit_result = None

    def make_peak_plots(self):
        # Make peak_fit_dict
        self.make_peak_params()
        # Do not plot if no peaks present
        if not self.peak_fit_dict['peak_idx_list']:
            self.structure_plot_widget.clear_gauss_curves()
            self.data['gauss_model'] = np.asarray([])
            self.data['gauss_peaks'] = []
            self.plot_data(_reset_view=0)
            return
        # Evaluate objective function over whole r range
        self.data['gauss_model'] = peak_fit.gauss_obj_func(self.peak_fit_dict['obj_params'],
                                                           self.data['rdf_x'], data=None,
                                                           **self.peak_fit_dict['obj_kwargs'])

        # Calculate residuals from objective function over min < r < max
        self.data['gauss_residuals'] = peak_fit.gauss_obj_func(self.peak_fit_dict['obj_params'],
                                                               self.data['fit_r'], data=self.data['fit_y'],
                                                               **self.peak_fit_dict['obj_kwargs'])
        # Calculate residuals from objective function over total r range
        self.data['gauss_residuals_full'] = peak_fit.gauss_obj_func(
                                                self.peak_fit_dict['obj_params'],
                                                self.data['rdf_x'], data=self.data['obj_fun'],
                                                **self.peak_fit_dict['obj_kwargs'])

        # Evaluate each individual peak over r range
        # List to store each individual peak
        self.data['gauss_peaks'] = []
        self.data['gauss_peaks_names'] = []
        for _peak_idx in self.peak_fit_dict['peak_idx_list']:
            _peak_idx_str = str(_peak_idx)
            _peak_model = peak_fit.gauss_obj_func(self.peak_fit_dict['params'+_peak_idx_str],
                                                  self.data['rdf_x'], data=None,
                                                  **self.peak_fit_dict['kwargs'+_peak_idx_str])
            self.data['gauss_peaks'].append(_peak_model)
            self.data['gauss_peaks_names'].append(self.peak_fit_dict['name'+_peak_idx_str])

        # Call update_plots
        self.structure_plot_widget.update_plots(self.data)
        self.structure_plot_widget.set_res_window(self.data)


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
        self.plot_view_label = QLabel('Select function: ')
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

        self.hbtn_layout = QGridLayout()
        self.hbtn_layout.setContentsMargins(0, 0, 0, 0)
        self.hbtn_layout.addWidget(self.rdf_btn, 0, 0)
        self.hbtn_layout.addWidget(self.tr_btn, 0, 1)

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
        self.setChecked(True)

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

    # Create custom signals emmitted when peaks added/changed
    guess_peak = pyqtSignal(int)
    # Relay signal from GaussianPeakGroupBox
    gauss_params_changed_relay = pyqtSignal()

    def __init__(self, *args):
        super(PolyatomicGroupBox, self).__init__(*args)
        self.setTitle('Fit Gaussians (Polyatomic)')
        self.setAlignment(Qt.AlignLeft)
        self.setStyleSheet('GroupBox::title{subcontrol-origin: margin; subcontrol-position: top left;}')
        self.setCheckable(True)
        self.setChecked(False)

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

        self.min_limit_input = QDoubleSpinBox()
        self.max_limit_input = QDoubleSpinBox()

        self.fit_peaks_btn = QPushButton('Fit')

        self.add_peak_btn = QPushButton('Add Peak')

    def style_widgets(self):

        self.fit_limit_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        self.min_limit_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.max_limit_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)

        self.min_limit_input.setMaximumWidth(82)
        self.max_limit_input.setMaximumWidth(82)

        self.min_limit_input.setSingleStep(0.1)
        self.min_limit_input.setDecimals(2)
        self.min_limit_input.setAlignment(Qt.AlignRight)

        self.max_limit_input.setSingleStep(0.1)
        self.max_limit_input.setDecimals(2)
        self.max_limit_input.setAlignment(Qt.AlignRight)

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
        self.outer_layout.addWidget(self.fit_peaks_btn)
        self.outer_layout.addWidget(self.add_peak_btn)

        self.setLayout(self.outer_layout)

    def init_peaks(self):
        self.peak_count = 1
        self.peak_dict = {}

    def create_signals(self):

        self.add_peak_btn.clicked.connect(self.add_peak)

        self.del_peak_signalMapper = QSignalMapper()
        self.del_peak_signalMapper.mapped[int].connect(self.del_peak)

        self.rename_peak_signalMapper = QSignalMapper()
        self.rename_peak_signalMapper.mapped[int].connect(self.rename_peak)

        self.min_limit_input.valueChanged.connect(self.gauss_params_changed_relay)
        self.max_limit_input.valueChanged.connect(self.gauss_params_changed_relay)

    def add_peak(self):

        if len(self.peak_dict) == 0:
            self.init_peaks()

        _new_peak = GaussianPeakGroupBox()
        _new_peak._unique_idx = self.peak_count
        _new_peak.title.setText(_new_peak.title.text() + ' #' + str(_new_peak._unique_idx) + ' ')
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

        # Re-emit parameter changed signal
        self.peak_dict[_new_peak._unique_idx].gauss_params_changed.connect(self.gauss_params_changed_relay)

        # Emit signal to guess peak parameters of new peak
        self.guess_peak.emit(_new_peak._unique_idx)

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
        # Update the plots
        self.gauss_params_changed_relay.emit()

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

    # Create signal emmitted on parameter change
    gauss_params_changed = pyqtSignal()

    def __init__(self, *args):
        super(GaussianPeakGroupBox, self).__init__(*args)

        self.create_widgets()
        self.style_widgets()
        self.create_layout()
        self.create_signals()

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
        self.lb_label = QLabel('min')
        self.ub_label = QLabel('max')

        self.N_label = QLabel('N: ')
        self.N_input = QLineEdit('2.0')
        self.N_lb = QLineEdit('0')
        self.N_ub = QLineEdit('10.0')
        self.N_refine = QCheckBox()

        self.r_label = QLabel('r: ')
        self.r_input = QLineEdit('2.0')
        self.r_lb = QLineEdit('1.0')
        self.r_ub = QLineEdit('3.0')
        self.r_refine = QCheckBox()

        self.s_label = QLabel('σ: ')
        self.s_input = QLineEdit('0.1')
        self.s_lb = QLineEdit('0')
        self.s_ub = QLineEdit('1.0')
        self.s_refine = QCheckBox()

        self.skew_toggle = QCheckBox('Skew?: ')
        self.xi_label = QLabel('ξ: ')
        self.xi_input = QLineEdit('0.01')
        self.xi_lb = QLineEdit('0')
        self.xi_ub = QLineEdit('2.0')
        self.xi_refine = QCheckBox()

    def style_widgets(self):

        self.title.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        self.del_peak_btn.setIcon(self.style().standardIcon(getattr(QStyle, 'SP_TitleBarCloseButton')))
        self.del_peak_btn.setMaximumWidth(25)
        self.del_peak_btn.setFlat(True)

        self.alpha_label.setAlignment(Qt.AlignVCenter | Qt.AlignCenter)
        self.beta_label.setAlignment(Qt.AlignVCenter | Qt.AlignCenter)

        self.refine_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        self.lb_label.setAlignment(Qt.AlignVCenter | Qt.AlignCenter)
        self.ub_label.setAlignment(Qt.AlignVCenter | Qt.AlignCenter)

        self.N_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.r_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.s_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.xi_label.setAlignment(Qt.AlignVCenter | Qt.AlignRight)

        self.skew_toggle.setLayoutDirection(Qt.RightToLeft)
        #self.skew_toggle.text.setAlignment(Qt.AlignVCenter)

        self.N_input.setMaximumWidth(100)
        self.r_input.setMaximumWidth(100)
        self.s_input.setMaximumWidth(100)
        self.xi_input.setMaximumWidth(100)

        self.N_lb.setMaximumWidth(100)
        self.N_ub.setMaximumWidth(100)
        self.r_lb.setMaximumWidth(100)
        self.r_ub.setMaximumWidth(100)
        self.s_lb.setMaximumWidth(100)
        self.s_ub.setMaximumWidth(100)
        self.xi_lb.setMaximumWidth(100)
        self.xi_ub.setMaximumWidth(100)

        # Set validators
        self.N_input.setValidator(QDoubleValidator())
        self.r_input.setValidator(QDoubleValidator())
        self.s_input.setValidator(QDoubleValidator())
        self.xi_input.setValidator(QDoubleValidator())
        self.N_lb.setValidator(QDoubleValidator())
        self.N_ub.setValidator(QDoubleValidator())
        self.r_lb.setValidator(QDoubleValidator())
        self.r_ub.setValidator(QDoubleValidator())
        self.s_lb.setValidator(QDoubleValidator())
        self.s_ub.setValidator(QDoubleValidator())
        self.xi_lb.setValidator(QDoubleValidator())
        self.xi_ub.setValidator(QDoubleValidator())

        self.N_refine.setChecked(True)
        self.r_refine.setChecked(True)
        self.s_refine.setChecked(True)
        self.xi_refine.setChecked(True)

        self.skew_toggle.setChecked(False)
        self.toggle_snd()

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

        self.params_grid_layout.addWidget(self.lb_label, 0, 2)
        self.params_grid_layout.addWidget(self.ub_label, 0, 3)

        self.params_grid_layout.addWidget(self.refine_label, 0, 4)

        self.params_grid_layout.addWidget(self.N_label, 1, 0)
        self.params_grid_layout.addWidget(self.N_input, 1, 1)
        self.params_grid_layout.addWidget(self.N_lb, 1, 2)
        self.params_grid_layout.addWidget(self.N_ub, 1, 3)
        self.params_grid_layout.addWidget(self.N_refine, 1, 4)

        self.params_grid_layout.addWidget(self.r_label, 2, 0)
        self.params_grid_layout.addWidget(self.r_input, 2, 1)
        self.params_grid_layout.addWidget(self.r_lb, 2, 2)
        self.params_grid_layout.addWidget(self.r_ub, 2, 3)
        self.params_grid_layout.addWidget(self.r_refine, 2, 4)

        self.params_grid_layout.addWidget(self.s_label, 3, 0)
        self.params_grid_layout.addWidget(self.s_input, 3, 1)
        self.params_grid_layout.addWidget(self.s_lb, 3, 2)
        self.params_grid_layout.addWidget(self.s_ub, 3, 3)
        self.params_grid_layout.addWidget(self.s_refine, 3, 4)

        self.params_grid_layout.addWidget(self.skew_toggle, 4, 0, 1, 2)

        self.params_grid_layout.addWidget(self.xi_label, 5, 0)
        self.params_grid_layout.addWidget(self.xi_input, 5, 1)
        self.params_grid_layout.addWidget(self.xi_lb, 5, 2)
        self.params_grid_layout.addWidget(self.xi_ub, 5, 3)
        self.params_grid_layout.addWidget(self.xi_refine, 5, 4)

        self.vlayout.addLayout(self.title_layout)
        self.vlayout.addLayout(self.atoms_grid_layout)
        self.vlayout.addLayout(self.params_grid_layout)

        self.setLayout(self.vlayout)

    def create_signals(self):
        self.alpha_input.currentIndexChanged.connect(self.emit_param_change)
        self.beta_input.currentIndexChanged.connect(self.emit_param_change)
        self.N_input.textChanged.connect(self.emit_param_change)
        self.r_input.textChanged.connect(self.emit_param_change)
        self.s_input.textChanged.connect(self.emit_param_change)
        self.xi_input.textChanged.connect(self.emit_param_change)
        self.skew_toggle.stateChanged.connect(self.toggle_snd)

    def emit_param_change(self):
        self.gauss_params_changed.emit()

    def populate_atom_list(self, _atom_list):
        self.alpha_input.clear()
        self.beta_input.clear()
        self.alpha_input.addItems(_atom_list)
        self.beta_input.addItems(_atom_list)
        self.alpha_input.setCurrentIndex(-1)
        self.beta_input.setCurrentIndex(-1)

    def toggle_snd(self):
        if self.skew_toggle.isChecked():
            self.xi_label.setEnabled(True)
            self.xi_input.setEnabled(True)
            self.xi_lb.setEnabled(True)
            self.xi_ub.setEnabled(True)
            self.xi_refine.setEnabled(True)
        else:
            self.xi_label.setEnabled(False)
            self.xi_input.setEnabled(False)
            self.xi_lb.setEnabled(False)
            self.xi_ub.setEnabled(False)
            self.xi_refine.setEnabled(False)
        self.emit_param_change()
