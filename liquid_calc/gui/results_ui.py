# -*- coding: utf-8 -*-
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018, Benedict J Heinen"
__email__ = "benedict.heinen@gmail.com"


from PyQt5.QtWidgets import QWidget, QFrame, QVBoxLayout, QHBoxLayout, \
                            QPushButton

import numpy as np
from . import plot_widgets
from . import utility
import core.core as core

class ResultsUI(QWidget):
        
    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(0)

        # Make plot widget
        self.results_plot_widget = plot_widgets.ResultsPlotWidget()

        # Make horizontal line separator
        self.hline = QFrame()
        self.hline.setFrameShape(QFrame.HLine)
        self.hline.setFrameShadow(QFrame.Sunken)
        self.hline.setObjectName("hline")
        
        # Make save buttons
        self.save_widget = SaveButtonsWidget()
        
        self.layout.addWidget(self.results_plot_widget)        
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
                     'gr_x': np.asarray([]), 'gr_y': np.asarray([]),
                     'rdf_x':  np.asarray([]), 'rdf_y':  np.asarray([]),
                     'rho': None, 'composition': None}
    
    def plot_data(self):
        # Re-calculate S(Q) from interference func
        self.data['sq_y'] = self.data['int_func'] + core.calc_S_inf(self.data['composition'], self.data['sq_x'])
        # Calculate g(r) & rdf(r) from s(q)
        # Calculate g(r) - pair-distribution function
        self.data['gr_x'], self.data['gr_y'] = core.calc_F_r(self.data['sq_x'], self.data['int_func'], self.data['rho'], function='pair_dist_func')
        # Calculate RDF - radial distribution function
        self.data['rdf_x'], self.data['rdf_y'] = core.calc_F_r(self.data['sq_x'], self.data['int_func'], self.data['rho'], function='radial_dist_func')
        self.results_plot_widget.update_plots(self.data)
        
        
        
        
    def save_sq(self):
        if not self.data['sq_y'].size:
            return
        __file_name = utility.get_filename(io='save', caption='Save S(Q)')
        if not __file_name:
            return
        __data = np.column_stack((self.data['sq_x'], self.data['sq_y']))
        np.savetxt(__file_name, __data, header='Refined S(Q)\nLiquidCalc v0.1', comments='#')

        
    def save_gr(self):
        if not self.data['gr_y'].size:
            return
        __file_name = utility.get_filename(io='save', caption='Save g(r)')
        if not __file_name:
            return
        __data = np.column_stack((self.data['gr_x'], self.data['gr_y']))
        np.savetxt(__file_name, __data, header='Refined g(r)\nLiquidCalc v0.1', comments='#')

        
    def save_rdf(self):
        if not self.data['rdf_y'].size:
            return
        __file_name = utility.get_filename(io='save', caption='Save RDF')
        if not __file_name:
            return 
        __data = np.column_stack((self.data['rdf_x'], self.data['rdf_y']))
        np.savetxt(__file_name, __data, header='Refined RDF\nLiquidCalc v0.1', comments='#')

        
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
