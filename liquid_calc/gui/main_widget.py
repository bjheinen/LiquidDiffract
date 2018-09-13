# -*- coding: utf-8 -*-
"""
GUI frontend for python implementation of Eggert method.
"""
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018, Benedict J Heinen"
__email__ = "benedict.heinen@gmail.com"


from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QTabWidget
from PyQt5.QtGui import QIcon
import numpy as np
import os

from . import bkg_ui
from . import optim_ui
from . import results_ui


class App(QMainWindow):
    
    def __init__(self, screen_size):
        super().__init__()
        # Set options here before running initUI to initialise
        self.title = 'LiquidCalc v0.1'
        # Get primary (current) screen dimensions
        self.screen_size = screen_size
        self.left = 0
        self.top = 0
        self.width = self.screen_size.width()
        self.height = self.screen_size.height()
        
        self.setWindowIcon(QIcon(os.path.join(os.path.abspath(os.getcwd()), 'data', 'icons', 'gs_icon.png')))

        self.initUI()
        
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left,self.top,self.width,self.height)
        
        self.table_widget = MainContainer(self)
        self.setCentralWidget(self.table_widget)
        self.showMaximized()


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
        #self.tabs.resize(300,200)
        self.tabs.addTab(self.bkg_ui,'Background Subtraction')
        self.tabs.addTab(self.optim_ui,'Refine Structure Factor')
        self.tabs.addTab(self.results_ui,'Calculate PDF')
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
        
        self.create_signal_tab_links()

    def create_signal_tab_links(self):
        self.bkg_ui.plots_changed.connect(self.bkg_plots_changed_slot)
        self.optim_ui.results_changed.connect(self.results_changed_slot)
        
    def bkg_plots_changed_slot(self):
        # Clear data from Optim UI (as S(Q) etc. need to be recalculated)
        self.optim_ui.data = {'cor_x': np.asarray([]), 'cor_y': np.asarray([]),
                              'cor_x_cut': np.asarray([]), 'cor_y_cut': np.asarray([]),
                              'sq_x':  np.asarray([]), 'sq_y':  np.asarray([]),
                              'fr_x':  np.asarray([]), 'fr_y':  np.asarray([]),
                              'mod_func': 'None'
                              }
        # Pass data to Optim UI
        self.optim_ui.data['cor_x'] = self.bkg_ui.data['cor_x']
        self.optim_ui.data['cor_y'] = self.bkg_ui.data['cor_y']
        self.optim_ui.plot_data()

    def results_changed_slot(self):
        # Clear old s(q), g(r), rdf(r)
        self.results_ui.clear_data()
        # Set rho
        if self.optim_ui.optim_config_widget.optim_options_gb.opt_check.isChecked():
            self.results_ui.data['rho'] = self.optim_ui.data['refined_rho']
        else:
            self.results_ui.data['rho'] = np.float(self.optim_ui.optim_config_widget.composition_gb.density_input.text())
        self.results_ui.data['int_func'] = self.optim_ui.data['impr_int_func']
        self.results_ui.data['sq_x'] = self.optim_ui.data['impr_iq_x']
        self.results_ui.data['composition'] = self.optim_ui.optim_config_widget.composition_gb.get_composition_dict()
        self.results_ui.data['mod_func'] = self.optim_ui.data['mod_func']
        self.results_ui.data['window_start'] = self.optim_ui.data['window_start']
        self.results_ui.plot_data()
