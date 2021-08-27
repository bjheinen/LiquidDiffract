# -*- coding: utf-8 -*-
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2018-2021, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

import numpy as np
from PyQt5.QtWidgets import QWidget, QFrame, QVBoxLayout, QHBoxLayout, \
                            QPushButton, QScrollArea
from LiquidDiffract.gui import plot_widgets
from LiquidDiffract.gui import utility
import LiquidDiffract.core.core as core
from LiquidDiffract.version import __appname__, __version__


class ResultsUI(QWidget):

    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(0)

        # Make plot widget
        self.results_plot_widget = plot_widgets.ResultsPlotWidget()
        self.plot_scroll_area = QScrollArea()
        self.plot_scroll_area.setWidget(self.results_plot_widget)
        self.plot_scroll_area.setWidgetResizable(True)
        self.plot_scroll_area.setFrameShape(QFrame.NoFrame)

        # Make horizontal line separator
        self.hline = QFrame()
        self.hline.setFrameShape(QFrame.HLine)
        self.hline.setFrameShadow(QFrame.Sunken)
        self.hline.setObjectName("hline")

        # Make save buttons
        self.save_widget = SaveButtonsWidget()

        self.layout.addWidget(self.plot_scroll_area)
        self.layout.addWidget(self.hline)
        self.layout.addWidget(self.save_widget)

        self.setLayout(self.layout)

        self.clear_data()

        self.create_signals()

    def create_signals(self):
        self.save_widget.save_sq_btn.clicked.connect(self.save_sq)
        self.save_widget.save_gr_btn.clicked.connect(self.save_gr)
        self.save_widget.save_rdf_btn.clicked.connect(self.save_rdf)

    def clear_data(self):

        self.data = {'sq_x': np.asarray([]), 'sq_y': np.asarray([]),
                     'int_func': np.asarray([]),
                     'mod_func': None, 'sq_y_mod': np.asarray([]),
                     'gr_x': np.asarray([]), 'gr_y': np.asarray([]),
                     'rdf_x': np.asarray([]), 'rdf_y': np.asarray([]),
                     'rho': None, 'composition': None}

    def plot_data(self):
        self.data['mod_int_func'] = self.data['int_func'] * core.get_mod_func(self.data['sq_x'], self.data['mod_func'], self.data['window_start'])
        # Re-calculate S(Q) from interference func
        self.data['S_inf'] = core.calc_S_inf(self.data['composition'], self.data['sq_x'], method=self.data['sq_method'])
        self.data['sq_y'] = self.data['int_func'] + self.data['S_inf']
        self.data['sq_y_mod'] = self.data['mod_int_func'] + self.data['S_inf']
        # Calculate g(r) & rdf(r) from s(q)
        # Calculate g(r) - pair-distribution function
        self.data['gr_x'], self.data['gr_y'] = core.calc_correlation_func(self.data['sq_x'], self.data['int_func'], self.data['rho'], N=self.fft_N,
                                                             mod_func=self.data['mod_func'], window_start=self.data['window_start'],
                                                             function='pair_dist_func')
        # Calculate RDF - radial distribution function
        self.data['rdf_x'], self.data['rdf_y'] = core.calc_correlation_func(self.data['sq_x'], self.data['int_func'], self.data['rho'], N=self.fft_N,
                                                               mod_func=self.data['mod_func'], window_start=self.data['window_start'],
                                                               function='radial_dist_func')
        # Optional rescale of Ashcroft-Langreth S(Q)/g(r)
        if self.rescale_AL == 1:
            # Renormalise Ashcroft-Langreth S(Q) to 1 at high Q
            self.data['sq_y'] = core.normalise_achroft_langreth_func(self.data['S_inf'], sq=self.data['sq_y'])
            self.data['sq_y_mod'] = core.normalise_achroft_langreth_func(self.data['S_inf'], sq=self.data['sq_y_mod'])
            # Remormalise Ashcroft-Langreth g(r) to between 0 and 1
            self.data['gr_y'] = core.normalise_achroft_langreth_func(self.data['S_inf'], gr=self.data['gr_y'])
        # Renormalise Ashcroft-Langreth RDF(r) to standard (Faber-Ziman) scaling
        self.data['rdf_y'] = core.normalise_achroft_langreth_func(self.data['S_inf'],
                                                                  rdf=self.data['rdf_y'],
                                                                  r=self.data['rdf_x'], rho_0=self.data['rho'])
        # Create data required for structure_ui
        self.data['tr_x'] = self.data['rdf_x']
        with np.errstate(divide='ignore', invalid='ignore'):
            self.data['tr_y'] = self.data['rdf_y'] / self.data['tr_x']

        self.results_plot_widget.update_plots(self.data)

    def save_sq(self):
        if not self.data['sq_y'].size:
            return
        # Base filename is the name of the input file (set in main_widget)
        __default_file_name = self.base_filename + '_SQ.dat'
        __file_name = utility.get_filename(io='save', directory=__default_file_name, caption='Save S(Q)')
        if not __file_name:
            return
        __header = (f'Refined S(Q)\n'
                    f'{__appname__} v{__version__}\n'
                    f'See refinement_log for more info\n'
                    f'Composition {{element: (Z, charge, n)}}: {self.data["composition"]}\n'
                    f'Rho = {self.data["rho"]}\n'
                    f'S(Q) formalism : {self.data["sq_method"]}\n'
                    f'S(Q) = i(Q) + S_inf\n'
                    f'S_inf = {np.str(self.data["S_inf"])}\n'
                    )
        if self.data['mod_func'] == 'None':
            __header += 'Q|S(Q)'
            __data = np.column_stack((self.data['sq_x'], self.data['sq_y']))
        else:
            __header += f'M(Q) : {self.data["mod_func"]}\n'
            if self.data['mod_func'] == 'Cosine-window':
                __header += f'Window-start : {self.data["window_start"]}\n'
            __header += f'Q|S(Q)|i(Q)*M(Q)'
            __data = np.column_stack((self.data['sq_x'], self.data['sq_y'], self.data['mod_int_func']))
        np.savetxt(__file_name, __data, header=__header, comments='#')

    def save_gr(self):
        if not self.data['gr_y'].size:
            return
        __default_file_name = self.base_filename + '_gr.dat'
        __file_name = utility.get_filename(io='save', directory=__default_file_name, caption='Save g(r)')
        if not __file_name:
            return
        __header = (f'Refined g(r)\n'
                    f'{__appname__} v{__version__}\n'
                    f'See refinement_log for more info\n'
                    f'Composition {{element: (Z, charge, n)}}: {self.data["composition"]}\n'
                    f'Rho = {self.data["rho"]}\n'
                    f'r|g(r)'
                    )
        __data = np.column_stack((self.data['gr_x'], self.data['gr_y']))
        np.savetxt(__file_name, __data, header=__header, comments='#')

    def save_rdf(self):
        if not self.data['rdf_y'].size:
            return
        __default_file_name = self.base_filename + '_rdf.dat'
        __file_name = utility.get_filename(io='save', directory=__default_file_name, caption='Save RDF')
        if not __file_name:
            return
        __header = (f'Refined rdf\n'
                    f'{__appname__} v{__version__}\n'
                    f'See refinement_log for more info\n'
                    f'Composition {{element: (Z, charge, n)}}: {self.data["composition"]}\n'
                    f'Rho = {self.data["rho"]}\n'
                    f'r|rdf'
                    )
        __data = np.column_stack((self.data['rdf_x'], self.data['rdf_y']))
        np.savetxt(__file_name, __data, header=__header, comments='#')


class SaveButtonsWidget(QWidget):
    def __init__(self):
        super(SaveButtonsWidget, self).__init__()

        self.horizontal_layout = QHBoxLayout()
        self.horizontal_layout.setContentsMargins(500, 9, 500, 0)
        self.horizontal_layout.setSpacing(250)

        self.save_sq_btn = QPushButton('Save S(Q)')
        self.save_sq_btn.setFlat(True)
        self.save_gr_btn = QPushButton('Save g(r)')
        self.save_gr_btn.setFlat(True)
        self.save_rdf_btn = QPushButton('Save RDF')
        self.save_rdf_btn.setFlat(True)

        self.horizontal_layout.addWidget(self.save_sq_btn)
        self.horizontal_layout.addWidget(self.save_gr_btn)
        self.horizontal_layout.addWidget(self.save_rdf_btn)
        self.setLayout(self.horizontal_layout)
