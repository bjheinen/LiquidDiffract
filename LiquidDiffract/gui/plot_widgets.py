# -*- coding: utf-8 -*-
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2018-2021, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import QWidget, QVBoxLayout
import pyqtgraph as pg
import numpy as np
from LiquidDiffract.core import data_utils

pg_options = {'leftButtonPan': False, 'background': 0.9, 'foreground': 0.15,
              'antialias': True}
pg.setConfigOptions(**pg_options)


class BkgPlotWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super(BkgPlotWidget, self).__init__(*args, **kwargs)

        self.layout = QVBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(8)
        self.create_plots()
        self.style_plots()

        self.create_signals()

        self.setLayout(self.layout)

    def create_plots(self):
        self.pg_layout_widget = pg.GraphicsLayoutWidget()
        self.pg_layout = pg.GraphicsLayout()
        self.pg_layout.setContentsMargins(0, 0, 0, 0)
        self.pg_layout_widget.setContentsMargins(0, 0, 0, 0)

        self.data_plot = CustomPlotItem()
        self.bkg_corrected_plot = CustomPlotItem()

        self.data_plot.plot(x=[], y=[], pen={'color': 0.1, 'width': 1.2})
        self.bkg_corrected_plot.plot(x=[], y=[], pen={'color': 0.1, 'width': 1.2})

        self.pg_layout.addItem(self.data_plot, row=1, col=0)
        self.pg_layout.addItem(self.bkg_corrected_plot, row=2, col=0)

        self.pg_layout_widget.addItem(self.pg_layout)

        self.layout.addWidget(self.pg_layout_widget)

    def style_plots(self):
        self.data_plot.setLabel('bottom', text='Q (1/A)')
        self.data_plot.setLabel('left', text='Intensity (a.u.)')

        self.bkg_corrected_plot.setLabel('bottom', text='Q (1/A)')
        self.bkg_corrected_plot.setLabel('left', text='Intensity (a.u.)')

        self.pos_label = pg.LabelItem(justify='right')
        self.pg_layout.addItem(self.pos_label, col=0, row=0)

    def update_plots(self, _data, _plot_raw):
        try:
            self.p1_a.clear()
            self.p2_a.clear()
            self.p3.clear()
        except AttributeError:
            pass
        try:
            self.p1_b.clear()
            self.p2_b.clear()
        except AttributeError:
            pass

        self.p1_a = self.data_plot.plot(x=_data['data_x'], y=_data['data_y'], pen={'color': 0.1, 'width': 1.2})
        self.p2_a = self.data_plot.plot(x=_data['bkg_x'], y=_data['bkg_y_sc'], pen={'color': '#342256', 'width': 1.2, 'style': Qt.DashLine})
        self.p3 = self.bkg_corrected_plot.plot(x=_data['cor_x'], y=_data['cor_y'], pen={'color': 0.1, 'width': 1.2})

        if _plot_raw:
            self.p1_b = self.data_plot.plot(x=_data['data_raw_x'], y=_data['data_raw_y'], pen=None, symbolPen={'color': 0.1}, symbolBrush=0.1, symbol='x', symbolSize=7)
            self.p2_b = self.data_plot.plot(x=_data['bkg_raw_x'], y=_data['bkg_raw_y_sc'], pen=None, symbolPen={'color': '#342256'}, symbolBrush='#342256', symbol='x', symbolSize=7)

        self.data_plot.vb.autoRange()
        self.bkg_corrected_plot.vb.autoRange()

    def create_signals(self):
        self.mouse_proxy = pg.SignalProxy(self.pg_layout.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved)

    def mouse_moved(self, __evt):
        # using signal proxy turns original arguments into a tuple
        __pos = __evt[0]
        if self.data_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.data_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)

            self.data_plot.vline.setPos(__mousePoint.x())
            self.data_plot.hline.setPos(__mousePoint.y())

            self.data_plot.vline.setPen((0, 135, 153), width=0.75)
            self.data_plot.hline.setPen((0, 135, 153), width=0.75)

            self.bkg_corrected_plot.vline.setPen(None)
            self.bkg_corrected_plot.hline.setPen(None)

        elif self.bkg_corrected_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.bkg_corrected_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)

            self.bkg_corrected_plot.vline.setPos(__mousePoint.x())
            self.bkg_corrected_plot.hline.setPos(__mousePoint.y())

            self.bkg_corrected_plot.vline.setPen((0, 135, 153), width=0.75)
            self.bkg_corrected_plot.hline.setPen((0, 135, 153), width=0.75)

            self.data_plot.vline.setPen(None)
            self.data_plot.hline.setPen(None)

    def set_mouse_pos_label(self, pos):
        _pos_str = (f'<span style="font-size: 11pt; color:#008799">x='
                    f'{pos.x():.2f}, y={pos.y():.2f}</span'
                    )
        self.pos_label.setText(_pos_str)


class OptimPlotWidget(QWidget):

    def __init__(self, *args, **kwargs):
        super(OptimPlotWidget, self).__init__(*args, **kwargs)

        self.layout = QVBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(8)
        self.create_plots()
        self.style_plots()

        self.create_signals()

        self.setLayout(self.layout)

    def create_plots(self):
        self.pg_layout_widget = pg.GraphicsLayoutWidget()
        self.pg_layout = pg.GraphicsLayout()
        self.pg_layout.setContentsMargins(0, 0, 0, 0)
        self.pg_layout_widget.setContentsMargins(0, 0, 0, 0)

        self.data_plot = CustomPlotItem()
        self.iq_plot = CustomPlotItem()
        self.dr_plot = CustomPlotItem()

        self.data_plot.plot(x=[], y=[])
        self.iq_plot.plot(x=[], y=[])
        self.dr_plot.plot(x=[], y=[], pen={})

        self.pg_layout.addItem(self.data_plot, row=1, col=0)
        self.pg_layout.addItem(self.iq_plot, row=2, col=0)
        self.pg_layout.addItem(self.dr_plot, row=3, col=0)

        self.pg_layout_widget.addItem(self.pg_layout)

        self.layout.addWidget(self.pg_layout_widget)

    def style_plots(self):
        self.data_plot.setLabel('bottom', text='Q (1/A)')
        self.data_plot.setLabel('left', text='Intensity (a.u.)')

        self.iq_plot.setLabel('bottom', text='Q (1/A)')
        self.iq_plot.setLabel('left', text='i(Q)')

        self.dr_plot.setLabel('bottom', text='r (A)')
        self.dr_plot.setLabel('left', text='D(r)')

        self.pos_label = pg.LabelItem(justify='right')
        self.pg_layout.addItem(self.pos_label, col=0, row=0)

    def update_plots(self, _data):
        try:
            self.p1_a.clear()
            self.p2_a.clear()
            self.p3_a.clear()
        except AttributeError:
            pass
        try:
            self.p1_b.clear()
        except AttributeError:
            pass
        try:
            self.p2_b.clear()
            self.p3_b.clear()
        except AttributeError:
            pass
        try:
            self.p2_c.clear()
        except AttributeError:
            pass

        # Some versions of pyqtgraph cannot produce plot if nan values present
        # First value in some arrays is nan (e.g. int func)
        # For interference function this should be == 0 - S_inf
        # Fix nan values by interpolation
        if np.isnan(_data['int_func']).any():
             _data['int_func'] = data_utils.interp_nan(_data['int_func'])
        if np.isnan(_data['impr_int_func']).any():
             _data['impr_int_func'] = data_utils.interp_nan(_data['impr_int_func'])
        if np.isnan(_data['dr_x']).any():
             _data['dr_x'] = data_utils.interp_nan(_data['dr_x'])
        if np.isnan(_data['impr_dr_x']).any():
             _data['impr_dr_x'] = data_utils.interp_nan(_data['impr_dr_x'])

        # For the D(r) the r step = pi/q_max. Because the data is padded this
        # q_max is larger than the original q_max. i.e. q_max = dq * 2**N/2
        # Then: r_max = 2**N/2 * dr = pi/(dq * 2**N/2) * 2**N/2 = pi/dq
        # Displaying the full length data is not necessary as most of the 
        # useful information is contained at fairly low r, and the length of r
        # is controlled only by sampling frequency in Q-space. At high r-values
        # the D(r) is dominated by ripples from the truncated integral (0-qmax)
        # The max value of r with 'real' q resolution is 1/dq.
        _window = 0
        try:
            _dq = _data['iq_x'][1] - _data['iq_x'][0]
            try:
                _window = np.argmax(_data['dr_x'] >= (1/_dq))
            except ValueError:
                 pass
        except IndexError:
            pass

        self.p1_a = self.data_plot.plot(x=_data['cor_x_cut'], y=_data['cor_y_cut'], pen={'color': 0.1, 'width': 1.2})
        self.p1_b = self.data_plot.plot(x=_data['cor_x_cut'], y=_data['rescaled_cor_y_cut'], pen={'color': '#342256', 'width': 1.2, 'style': Qt.DashLine})
        self.p2_a = self.iq_plot.plot(x=_data['iq_x'], y=_data['int_func'], pen={'color': 0.1, 'width': 1.2})
        self.p2_b = self.iq_plot.plot(x=_data['impr_iq_x'], y=_data['impr_int_func'], pen={'color': '#342256', 'width': 1.2, 'style': Qt.DashLine})
        self.p3_a = self.dr_plot.plot(x=_data['dr_x'][:_window], y=_data['dr_y'][:_window], pen={'color': 0.1, 'width': 1.2})
        self.p3_b = self.dr_plot.plot(x=_data['impr_dr_x'][:_window], y=_data['impr_dr_y'][:_window], pen={'color': '#342256', 'width': 1.2, 'style': Qt.DashLine})

        if _data['mod_func'] != 'None':
            # Apply modification function and handle any nans that arise 
            _modified_int_func = _data['modification'] * _data['int_func']
            if np.isnan(_modified_int_func).any():
                 _modified_int_func = data_utils.interp_nan(_modified_int_func)
            self.p2_c = self.iq_plot.plot(x=_data['iq_x'], y=_modified_int_func, pen={'color': '#342256', 'width': 0.8, 'style': Qt.DashLine})

        self.data_plot.vb.autoRange()
        self.iq_plot.vb.autoRange()
        self.dr_plot.vb.autoRange()

    def create_signals(self):
        self.mouse_proxy = pg.SignalProxy(self.pg_layout.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved)

    def mouse_moved(self, __evt):
        # using signal proxy turns original arguments into a tuple
        __pos = __evt[0]
        if self.data_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.data_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)

            self.data_plot.vline.setPos(__mousePoint.x())
            self.data_plot.hline.setPos(__mousePoint.y())

            self.data_plot.vline.setPen((0, 135, 153), width=0.75)
            self.data_plot.hline.setPen((0, 135, 153), width=0.75)

            self.iq_plot.vline.setPen(None)
            self.iq_plot.hline.setPen(None)
            self.dr_plot.vline.setPen(None)
            self.dr_plot.hline.setPen(None)

        elif self.iq_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.iq_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)

            self.iq_plot.vline.setPos(__mousePoint.x())
            self.iq_plot.hline.setPos(__mousePoint.y())

            self.iq_plot.vline.setPen((0, 135, 153), width=0.75)
            self.iq_plot.hline.setPen((0, 135, 153), width=0.75)

            self.data_plot.vline.setPen(None)
            self.data_plot.hline.setPen(None)
            self.dr_plot.vline.setPen(None)
            self.dr_plot.hline.setPen(None)

        elif self.dr_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.dr_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)

            self.dr_plot.vline.setPos(__mousePoint.x())
            self.dr_plot.hline.setPos(__mousePoint.y())

            self.dr_plot.vline.setPen((0, 135, 153), width=0.75)
            self.dr_plot.hline.setPen((0, 135, 153), width=0.75)

            self.data_plot.vline.setPen(None)
            self.data_plot.hline.setPen(None)
            self.iq_plot.hline.setPen(None)
            self.iq_plot.vline.setPen(None)

    def set_mouse_pos_label(self, pos):
        _pos_str = (f'<span style="font-size: 11pt; color:#008799">x='
                    f'{pos.x():.2f}, y={pos.y():.2f}</span'
                    )
        self.pos_label.setText(_pos_str)


class ResultsPlotWidget(QWidget):

    def __init__(self, *args, **kwargs):
        super(ResultsPlotWidget, self).__init__(*args, **kwargs)

        self.layout = QVBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(8)
        self.create_plots()
        self.style_plots()

        self.create_signals()

        self.setLayout(self.layout)

    def create_plots(self):
        self.pg_layout_widget = pg.GraphicsLayoutWidget()
        self.pg_layout = pg.GraphicsLayout()
        self.pg_layout.setContentsMargins(0, 0, 0, 0)
        self.pg_layout_widget.setContentsMargins(0, 0, 0, 0)

        self.sq_plot = CustomPlotItem()
        self.gr_plot = WindowedPlotItem()
        self.rdf_plot = WindowedPlotItem()

        self.sq_plot.plot(x=[], y=[])
        self.gr_plot.plot(x=[], y=[])
        self.rdf_plot.plot(x=[], y=[])
        self.rdf_plot.setXLink(self.gr_plot)

        self.pg_layout.addItem(self.sq_plot, row=1, col=0)
        self.pg_layout.addItem(self.gr_plot, row=2, col=0)
        self.pg_layout.addItem(self.rdf_plot, row=3, col=0)

        self.pg_layout_widget.addItem(self.pg_layout)

        self.layout.addWidget(self.pg_layout_widget)

    def style_plots(self):
        self.sq_plot.setLabel('bottom', text='Q (1/A)')
        self.sq_plot.setLabel('left', text='S(Q)')

        self.gr_plot.setLabel('bottom', text='r (A)')
        self.gr_plot.setLabel('left', text='g(r)')

        self.rdf_plot.setLabel('bottom', text='r (A)')
        self.rdf_plot.setLabel('left', text='RDF(r)')

        self.pos_label = pg.LabelItem(justify='right')
        self.pg_layout.addItem(self.pos_label, col=0, row=0)

    def update_plots(self, _data):
        try:
            self.p1.clear()
            self.p2.clear()
            self.p3.clear()
        except AttributeError:
            pass
        try:
            self.p1_mod.clear()
        except AttributeError:
            pass

        # Some versions of pyqtgraph cannot produce plot if nan values present
        # Fix nan values by interpolation
        if np.isnan(_data['sq_y']).any():
             _data['sq_y'] = data_utils.interp_nan(_data['sq_y'])
        if np.isnan(_data['gr_y']).any():
             _data['gr_y'] = data_utils.interp_nan(_data['gr_y'])
        if np.isnan(_data['rdf_y']).any():
             _data['rdf_y'] = data_utils.interp_nan(_data['rdf_y'])
        if np.isnan(_data['sq_y_mod']).any():
            _data['sq_y_mod'] = data_utils.interp_nan(_data['sq_y_mod'])
        _window = 0
        # Determine data window for Q-space resolution, dq
        try:
            _dq = _data['sq_x'][1] - _data['sq_x'][0]
            try:
                _window = np.argmax(_data['gr_x'] >= (1/_dq))
            except ValueError:
                 pass
        except IndexError:
            pass

        self.p1 = self.sq_plot.plot(x=_data['sq_x'], y=_data['sq_y'], pen={'color': 0.1, 'width': 1.2})
        if _data['mod_func'] != None:
            self.p1_mod = self.sq_plot.plot(x=_data['sq_x'], y=_data['sq_y_mod'], pen={'color': '#342256', 'width': 0.8, 'style': Qt.DashLine})
        self.p2 = self.gr_plot.plot(x=_data['gr_x'][:_window], y=_data['gr_y'][:_window], pen={'color': 0.1, 'width': 1.2})
        self.p3 = self.rdf_plot.plot(x=_data['rdf_x'][:_window], y=_data['rdf_y'][:_window], pen={'color': 0.1, 'width': 1.2})

        # Limit the inital view to important information
        try:
            self.x_max = _data['sq_x'][-1]
        except IndexError:
            return

        _gr_cut = np.nan_to_num(_data['gr_y'][np.where(_data['gr_x'] < self.x_max)])
        self.y_min_gr = np.min(_gr_cut)
        self.y_max_gr = np.max(_gr_cut)

        _rdf_cut = np.nan_to_num(_data['rdf_y'][np.where(_data['rdf_x'] < self.x_max)])
        self.y_min_rdf = np.min(_rdf_cut)
        self.y_max_rdf = np.max(_rdf_cut)

        self.set_gr_window()
        self.set_rdf_window()

        self.sq_plot.vb.autoRange()

    def set_gr_window(self):
        try:
            self.gr_plot.vb.setRange(xRange=(0, self.x_max),
                                     yRange=(self.y_min_gr, self.y_max_gr))
        except:
            return

    def set_rdf_window(self):
        try:
            self.rdf_plot.vb.setRange(xRange=(0, self.x_max),
                                      yRange=(self.y_min_rdf, self.y_max_rdf))
        except:
            return

    def create_signals(self):
        self.mouse_proxy = pg.SignalProxy(self.pg_layout.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved)
        self.gr_plot.reset_window.connect(self.set_gr_window)
        self.rdf_plot.reset_window.connect(self.set_rdf_window)

    def mouse_moved(self, __evt):
        # Using signal proxy turns original args into tuple
        __pos = __evt[0]
        if self.sq_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.sq_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)

            self.sq_plot.vline.setPos(__mousePoint.x())
            self.sq_plot.hline.setPos(__mousePoint.y())

            self.sq_plot.vline.setPen((0, 135, 153), width=0.75)
            self.sq_plot.hline.setPen((0, 135, 153), width=0.75)

            self.gr_plot.vline.setPen(None)
            self.gr_plot.hline.setPen(None)
            self.rdf_plot.vline.setPen(None)
            self.rdf_plot.hline.setPen(None)

        elif self.gr_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.gr_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)

            self.gr_plot.vline.setPos(__mousePoint.x())
            self.gr_plot.hline.setPos(__mousePoint.y())

            self.gr_plot.vline.setPen((0, 135, 153), width=0.75)
            self.gr_plot.hline.setPen((0, 135, 153), width=0.75)

            self.sq_plot.vline.setPen(None)
            self.sq_plot.hline.setPen(None)
            self.rdf_plot.vline.setPen(None)
            self.rdf_plot.hline.setPen(None)

        elif self.rdf_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.rdf_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)

            self.rdf_plot.vline.setPos(__mousePoint.x())
            self.rdf_plot.hline.setPos(__mousePoint.y())

            self.rdf_plot.vline.setPen((0, 135, 153), width=0.75)
            self.rdf_plot.hline.setPen((0, 135, 153), width=0.75)

            self.sq_plot.vline.setPen(None)
            self.sq_plot.hline.setPen(None)
            self.gr_plot.hline.setPen(None)
            self.gr_plot.vline.setPen(None)

    def set_mouse_pos_label(self, pos):
        _pos_str = (f'<span style="font-size: 11pt; color:#008799">x='
                    f'{pos.x():.2f}, y={pos.y():.2f}</span'
                    )
        self.pos_label.setText(_pos_str)


class StructurePlotWidget(QWidget):

    def __init__(self, *args, **kwargs):
        super(StructurePlotWidget, self).__init__(*args, **kwargs)

        self.layout = QVBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(8)
        self.create_plots()
        self.style_plots()

        self.create_signals()

        self.setLayout(self.layout)

    def create_plots(self):

        # Seperate widgets for integrating RDF, int Tr, and curve fitting
        # Widget for integrating across RDF
        self.rdf_int_layout_widget = pg.GraphicsLayoutWidget()
        self.rdf_int_layout_widget.setContentsMargins(0, 0, 0, 0)
        self.rdf_int_layout = pg.GraphicsLayout()
        self.rdf_int_layout.setContentsMargins(0, 0, 0, 0)
        # Widget for integrating T(r)
        # Separate widget used as easier than handling all the children
        self.tr_int_layout_widget = pg.GraphicsLayoutWidget()
        self.tr_int_layout_widget.setContentsMargins(0, 0, 0, 0)
        self.tr_int_layout = pg.GraphicsLayout()
        self.tr_int_layout.setContentsMargins(0, 0, 0, 0)

        # Curve fitting widget
        # One widget for RDF and Tr and only data and label change - all the
        # user added curves stay the same on switch
        self.fit_layout_widget = pg.GraphicsLayoutWidget()
        self.fit_layout_widget.setContentsMargins(0, 0, 0, 0)
        self.fit_layout = pg.GraphicsLayout()
        self.fit_layout.setContentsMargins(0, 0, 0, 0)
        self.fit_layout.setSpacing(0)

        # RDF(r) plot
        self.rdf_plot = WindowedPlotItem()
        self.rdf_plot.plot(x=[], y=[])
        self.rdf_int_layout.addItem(self.rdf_plot, row=1, col=0)
        self.rdf_int_layout_widget.addItem(self.rdf_int_layout)

        # T(r) plot
        self.tr_plot = WindowedPlotItem()
        self.tr_plot.plot(x=[], y=[])
        self.tr_int_layout.addItem(self.tr_plot, row=1, col=0)
        self.tr_int_layout_widget.addItem(self.tr_int_layout)

        # Curve (Gaussian) fitting plot
        self.fit_plot = WindowedPlotItem()
        self.fit_plot.plot(x=[], y=[])

        # Create fit range selection tool
        self.fit_limits = pg.LinearRegionItem(values=(0, 10), orientation='vertical', span=(0, 1), movable=True, swapMode='block', brush=pg.mkBrush(None), hoverBrush=pg.mkBrush(None))
        self.deselect_lower = pg.LinearRegionItem(values=(0, 0), orientation='vertical', span=(0, 1), movable=False, swapMode='block', pen=pg.mkPen(None), brush=pg.mkBrush((229,229,229,205)), hoverPen=pg.mkPen(None), hoverBrush=None)
        self.deselect_upper = pg.LinearRegionItem(values=(10, 10), orientation='vertical', span=(0, 1), movable=False, swapMode='block', pen=pg.mkPen(None), brush=pg.mkBrush((229,229,229,205)), hoverPen=pg.mkPen(None), hoverBrush=None)
        self.fit_plot.addItem(self.fit_limits)
        self.fit_plot.addItem(self.deselect_lower)
        self.fit_plot.addItem(self.deselect_upper)

        # Fitting residuals plot
        self.res_plot = WindowedPlotItem()
        self.res_plot.plot(x=[], y=[])
        self.res_plot.setXLink(self.fit_plot)
        self.fit_layout.addItem(self.fit_plot, row=1, col=0)
        self.fit_layout.addItem(self.res_plot, row=2, col=0)
        self.fit_layout.layout.setRowStretchFactor(1, 3)
        self.fit_layout.layout.setRowStretchFactor(2, 2)
        self.fit_layout_widget.addItem(self.fit_layout)

        # Add widgets to the layout
        self.layout.addWidget(self.rdf_int_layout_widget)
        self.layout.addWidget(self.tr_int_layout_widget)
        self.layout.addWidget(self.fit_layout_widget)

        # Set RDF(r) plot (integration) as default view
        self.rdf_int_layout_widget.setVisible(True)
        self.tr_int_layout_widget.setVisible(False)
        self.fit_layout_widget.setVisible(False)

    def style_plots(self):

        self.rdf_plot.setLabel('bottom', text='r (A)')
        self.rdf_plot.setLabel('left', text='RDF(r)')

        self.tr_plot.setLabel('bottom', text='r (A)')
        self.tr_plot.setLabel('left', text='T(r)')

        # No label on fit_plot - just on residuals
        self.fit_plot.setLabel('left', text='RDF(r)')
        self.res_plot.setLabel('bottom', text='r (A)')
        self.res_plot.setLabel('left', text='residuals')
        # Fix axis width so plots are aligned
        self.fit_plot.getAxis('bottom').setStyle(showValues=False)
        self.fit_plot.getAxis('left').setWidth(55)
        self.res_plot.getAxis('left').setWidth(55)

        self.rdf_xaxis = pg.InfiniteLine(pos=0, angle=0, movable=False, pen={'color': 'k', 'width': 0.75})
        self.rdf_plot.addItem(self.rdf_xaxis)
        self.tr_xaxis = pg.InfiniteLine(pos=0, angle=0, movable=False, pen={'color': 'k', 'width': 0.75})
        self.tr_plot.addItem(self.tr_xaxis)
        self.fit_xaxis = pg.InfiniteLine(pos=0, angle=0, movable=False, pen={'color': 'k', 'width': 0.75})
        self.fit_plot.addItem(self.fit_xaxis)
        self.res_xaxis = pg.InfiniteLine(pos=0, angle=0, movable=False, pen={'color': 'k', 'width': 0.75})
        self.res_plot.addItem(self.res_xaxis)
        self.fit_xaxis.setZValue(3)

        # Make fit_limits region invisible to start
        self.toggle_fit_limits(False)

        # Add position label to all plot views
        self.pos_label_rdf = pg.LabelItem(justify='right')
        self.rdf_int_layout.addItem(self.pos_label_rdf, col=0, row=0)

        self.pos_label_tr = pg.LabelItem(justify='right')
        self.tr_int_layout.addItem(self.pos_label_tr, col=0, row=0)

        self.pos_label_fit = pg.LabelItem(justify='right')
        self.fit_layout.addItem(self.pos_label_fit, col=0, row=0)

        self.rdf_plot.vline.setPen((0, 135, 153), width=0.75)
        self.rdf_plot.hline.setPen((0, 135, 153), width=0.75)

        self.tr_plot.vline.setPen((0, 135, 153), width=0.75)
        self.tr_plot.hline.setPen((0, 135, 153), width=0.75)

        self.fit_plot.vline.setPen((0, 135, 153), width=0.75)
        self.fit_plot.hline.setPen((0, 135, 153), width=0.75)
        self.res_plot.vline.setPen((0, 135, 153), width=0.75)
        self.res_plot.hline.setPen((0, 135, 153), width=0.75)

        # Set Z-value of selection region to plot above curves
        self.deselect_lower.setZValue(1)
        self.deselect_upper.setZValue(1)
        self.fit_limits.setZValue(2)
        self.fit_plot.vline.setZValue(6)
        self.fit_plot.hline.setZValue(6)
        self.fit_plot.getAxis('bottom').setZValue(7)

        # Colors for gaussian curves
        self.gauss_colors = [(26,121,199), (22,171,69), (246,174,45),
                             (242,100,25), (99,18,153), (5,117,59),
                             (196,35,72), (90,73,57)]

    def toggle_fit_limits(self, toggle_bool, obj_x=None, sq_x=None):
        self.fit_limits.setVisible(toggle_bool)
        self.deselect_lower.setVisible(toggle_bool)
        self.deselect_upper.setVisible(toggle_bool)
        if toggle_bool:
            # Set fit_limits default region and bounds to new data
            self.deselect_lower.setRegion([0,0])
            self.deselect_upper.setRegion([sq_x[-1],obj_x[-1]])
            self.fit_limits.setRegion([obj_x[0],sq_x[-1]])
            self.fit_limits.setBounds([0, obj_x[-1]])

    def clear_gauss_curves(self):
        try:
            self.p_gauss_model.clear()
            self.p_res.clear()
            self.p_res_full.clear()
            # Clear peak curves
            for _p_peak in self.p_peaks:
                _p_peak.clear()
                del _p_peak
            self.p_peaks = []
            # Delete peak labels
            for _peak_label in self.p_peaks_labels:
                _peak_label.deleteLater()
                del _peak_label
            self.p_peaks_labels = []
            # Delete arrows (pg.ArrowItem has no deleteLater method)
            for _peak_arrow in self.p_peaks_arrows:
                self.fit_plot.removeItem(_peak_arrow)
                del _peak_arrow
            self.p_peaks_arrows = []

        except AttributeError:
            pass

    def clear_plots(self, _clear_all=False):
        try:
            self.p_rdf.clear()
            self.p_tr.clear()
            self.p_fit.clear()
        except AttributeError:
            pass

        self.clear_gauss_curves()

        try:
            self.p_Na.clear()
        except AttributeError:
            pass

        try:
            self.p_Nb.clear()
        except AttributeError:
            pass

        try:
            self.p_Nc.clear()
        except AttributeError:
            pass

        if _clear_all:
            self.toggle_fit_limits(False)
            try:
                self.rdf_plot.removeItem(self.r0_line_rdf)
                self.tr_plot.removeItem(self.r0_line_tr)
                self.rdf_plot.removeItem(self.rpmax_line_rdf)
                self.tr_plot.removeItem(self.rpmax_line_tr)
                self.rdf_plot.removeItem(self.rmax_line_rdf)
                self.tr_plot.removeItem(self.rmax_line_tr)
                self.rdf_plot.removeItem(self.rmin_line_rdf)
                self.tr_plot.removeItem(self.rmin_line_rdf)

                self.r0_line_rdf.deleteLater()
                self.r0_line_tr.deleteLater()
                self.rpmax_line_rdf.deleteLater()
                self.rpmax_line_tr.deleteLater()
                self.rmax_line_rdf.deleteLater()
                self.rmax_line_tr.deleteLater()
                self.rmin_line_rdf.deleteLater()
                self.rmin_line_rdf.deleteLater()

                del self.r0_line_rdf
                del self.r0_line_tr
                del self.rpmax_line_rdf
                del self.rpmax_line_tr
                del self.rmax_line_rdf
                del self.rmax_line_tr
                del self.rmin_line_rdf
                del self.rmin_line_rdf
            except AttributeError:
                pass


    def update_plots(self, _data):

        # First clear the plots
        self.clear_plots()

        # Some versions of pyqtgraph cannot produce plot if nan values present
        # Fix nan values by interpolation
        if np.isnan(_data['rdf_y']).any():
             _data['rdf_y'] = data_utils.interp_nan(_data['rdf_y'])
        if np.isnan(_data['tr_y']).any():
             _data['tr_y'] = data_utils.interp_nan(_data['tr_y'])

        # Interpolate data for smoother plots
        # Ignoring for now
        #_data['rdf_x'], _data['rdf_y'] = data_utils.rebin_data(_data['rdf_x'], _data['rdf_y'], dx=0.01)
        #_data['tr_x'], _data['tr_y'] = data_utils.rebin_data(_data['tr_x'], _data['tr_y'], dx=0.01)

        # Create areas of N integrals
        _Na_area_idx = np.where((_data['tr_x'] > _data['r0']) & (_data['tr_x'] < _data['rpmax']))        
        _Nb_area_idx = np.where((_data['rdf_x'] > _data['r0']) & (_data['rdf_x'] < _data['rmax']))
        _Nc_area_idx = np.where((_data['rdf_x'] > _data['r0']) & (_data['rdf_x'] < _data['rmin']))

        _r0_tr_pt = data_utils.interp_data(_data['tr_x'], _data['tr_y'], _data['r0'])
        _r0_rdf_pt = data_utils.interp_data(_data['rdf_x'], _data['rdf_y'], _data['r0'])
        _rpmax_tr_pt = data_utils.interp_data(_data['tr_x'], _data['tr_y'], _data['rpmax'])        
        _rmax_rdf_pt = data_utils.interp_data(_data['rdf_x'], _data['rdf_y'], _data['rmax'])
        _rmin_rdf_pt = data_utils.interp_data(_data['rdf_x'], _data['rdf_y'], _data['rmin'])

        _Na_area_x = np.concatenate(([_data['r0']], _data['tr_x'][_Na_area_idx], [_data['rpmax']]))
        _Na_area_y = np.concatenate((_r0_tr_pt, _data['tr_y'][_Na_area_idx], _rpmax_tr_pt))
        _Nb_area_x = np.concatenate(([_data['r0']], _data['rdf_x'][_Nb_area_idx], [_data['rmax']]))
        _Nb_area_y = np.concatenate((_r0_rdf_pt, _data['rdf_y'][_Nb_area_idx], _rmax_rdf_pt))
        _Nc_area_x = np.concatenate(([_data['r0']], _data['rdf_x'][_Nc_area_idx], [_data['rmin']]))
        _Nc_area_y = np.concatenate((_r0_rdf_pt, _data['rdf_y'][_Nc_area_idx], _rmin_rdf_pt))

        # Plot fill areas for integrals behind function/data curve
        # Nc is plotted first/behind as it extends beyond Nb
        # Only plot if r0 < r
        if _data['r0'] < _data['rmin']:
            self.p_Nc = self.rdf_plot.plot(x=_Nc_area_x, y=_Nc_area_y, brush=(26,121,199, 245), pen=None, fillLevel=0)
        if _data['r0'] < _data['rpmax']:
            self.p_Na = self.tr_plot.plot(x=_Na_area_x, y=_Na_area_y, brush=(199,26,74, 245), pen=None, fillLevel=0)
        if _data['r0'] < _data['rmax']:
            self.p_Nb = self.rdf_plot.plot(x=_Nb_area_x, y=_Nb_area_y, brush=(199,26,74, 245), pen=None, fillLevel=0)

        # Plot data/functions RDF(r), T(r)
        self.p_rdf = self.rdf_plot.plot(x=_data['rdf_x'], y=_data['rdf_y'], pen={'color': 0.1, 'width': 1.2})
        self.p_tr = self.tr_plot.plot(x=_data['tr_x'], y=_data['tr_y'], pen={'color': 0.1, 'width': 1.2})
        self.p_fit = self.fit_plot.plot(x=_data['rdf_x'], y=_data['obj_fun'], pen={'color': 0.1, 'width': 1.2})

        # Limit the inital view to important information
        try:
            self.x_max = _data['sq_x'][-1]
        except IndexError:
            return

        # Plot integral limits as movable InfLines 
        try:
            # Set line positions
            self.r0_line_rdf.setPos(_data['r0'])
            self.r0_line_tr.setPos(_data['r0'])

            self.rpmax_line_rdf.setPos(_data['rpmax'])
            self.rpmax_line_tr.setPos(_data['rpmax'])

            self.rmax_line_rdf.setPos(_data['rmax'])
            self.rmax_line_tr.setPos(_data['rmax'])

            self.rmin_line_rdf.setPos(_data['rmin'])
            self.rmin_line_tr.setPos(_data['rmin'])

        except AttributeError:
            # Initial creation of lines
            # Create draggable lines for r0, rmax, rpmax and rmin
            self.r0_line_rdf = pg.InfiniteLine(pos=_data['r0'], movable=True, label='r0', bounds=[0, self.x_max], 
                                               labelOpts={'position': 0.75, 'movable': True, 'color': (0,10,40)}, name='r0')
            self.r0_line_tr = pg.InfiniteLine(pos=_data['r0'], movable=True, label='r0', bounds=[0, self.x_max], 
                                              labelOpts={'position': 0.75, 'movable': True, 'color': (0,10,40)}, name='r0')

            self.rpmax_line_rdf = pg.InfiniteLine(pos=_data['rpmax'], movable=True, label='r\'max', bounds=[0, self.x_max], 
                                                  labelOpts={'position': 0.8, 'movable': True, 'color': (0,10,40)}, name='rpmax')
            self.rpmax_line_tr = pg.InfiniteLine(pos=_data['rpmax'], movable=True, label='r\'max', bounds=[0, self.x_max], 
                                                 labelOpts={'position': 0.8, 'movable': True, 'color': (0,10,40)}, name='rpmax')

            self.rmax_line_rdf = pg.InfiniteLine(pos=_data['rmax'], movable=True, label='rmax', bounds=[0, self.x_max], 
                                                 labelOpts={'position': 0.85, 'movable': True, 'color': (0,10,40)}, name='rmax')
            self.rmax_line_tr = pg.InfiniteLine(pos=_data['rmax'], movable=True, label='rmax', bounds=[0, self.x_max], 
                                                labelOpts={'position': 0.85, 'movable': True, 'color': (0,10,40)}, name='rmax')

            self.rmin_line_rdf = pg.InfiniteLine(pos=_data['rmin'], movable=True, label='rmin', bounds=[0, self.x_max], 
                                                 labelOpts={'position': 0.75, 'movable': True, 'color': (0,10,40)}, name='rmin')
            self.rmin_line_tr = pg.InfiniteLine(pos=_data['rmin'], movable=True, label='rmin', bounds=[0, self.x_max], 
                                                labelOpts={'position': 0.75, 'movable': True, 'color': (0,10,40)}, name='rmin')

            self.rdf_plot.addItem(self.r0_line_rdf)
            self.tr_plot.addItem(self.r0_line_tr)

            self.rdf_plot.addItem(self.rpmax_line_rdf)
            self.tr_plot.addItem(self.rpmax_line_tr)

            self.rdf_plot.addItem(self.rmax_line_rdf)
            self.tr_plot.addItem(self.rmax_line_tr)

            self.rdf_plot.addItem(self.rmin_line_rdf)
            self.tr_plot.addItem(self.rmin_line_tr)        

        if _data['gauss_model'].size:
            # Plot individual gaussian peaks
            self.p_peaks = []
            self.p_peaks_labels = []
            self.p_peaks_arrows = []
            for i, _peak_model in enumerate(_data['gauss_peaks']):
                # Get color - use i%8 to cycle list of colors
                _color = self.gauss_colors[i%len(self.gauss_colors)]
                _p_peak = self.fit_plot.plot(x=_data['rdf_x'], y=_peak_model, pen={'color': _color, 'width': 1.2})
                _p_peak.setZValue(4)
                self.p_peaks.append(_p_peak)
                # Label peak at location of curve maximum
                _max_idx = np.argmax(_peak_model)
                _peak_loc = (_data['rdf_x'][_max_idx], _peak_model[_max_idx])
                _peak_arrow = pg.ArrowItem(pos=_peak_loc, angle=260, pen=_color, brush=_color, tailLen=30, tailWidth=2, headLen=10, headWidth=4, pxMode=True)
                _peak_label = pg.TextItem(text=_data['gauss_peaks_names'][i], color=_color, anchor=(1.0,2.47))
                _peak_label.setPos(*_peak_loc)
                self.fit_plot.addItem(_peak_label)
                self.fit_plot.addItem(_peak_arrow)
                self.p_peaks_labels.append(_peak_label)
                self.p_peaks_arrows.append(_peak_arrow)

            # Plot gaussian model (sum of individual peaks)
            # change color of this
            self.p_gauss_model = self.fit_plot.plot(x=_data['rdf_x'], y=_data['gauss_model'], pen={'color': (199,26,74), 'width': 1.8, 'style':Qt.DashLine})
            self.p_gauss_model.setZValue(5)
            # Plot residuals over whole r range in grey
            self.p_res_full = self.res_plot.plot(x=_data['rdf_x'], y=_data['gauss_residuals_full'], pen={'color': 0.75, 'width': 1.2})
            # Plot residuals over fit range in colour
            self.p_res = self.res_plot.plot(x=_data['fit_r'], y=_data['gauss_residuals'], pen={'color': (26,121,199), 'width': 1.5})

        # Force repaint
        self.rdf_int_layout.update()
        self.tr_int_layout.update()
        self.fit_layout.update()

    def update_plot_windows(self, _data, _reset_view):
        _rdf_cut = np.nan_to_num(_data['rdf_y'][np.where(_data['rdf_x'] < self.x_max)])
        self.y_min_rdf = np.min(_rdf_cut)
        self.y_max_rdf = np.max(_rdf_cut)

        _tr_cut = np.nan_to_num(_data['tr_y'][np.where(_data['tr_x'] < self.x_max)])
        self.y_min_tr = np.min(_tr_cut)
        self.y_max_tr = np.max(_tr_cut)

        _obj_fun_cut = np.nan_to_num(_data['obj_fun'][np.where(_data['rdf_x'] < self.x_max)])
        self.y_min_fit = np.min(_obj_fun_cut)
        self.y_max_fit = np.max(_obj_fun_cut)

        if _reset_view:
            self.set_rdf_window()
            self.set_tr_window()
            self.set_fit_window()
            self.set_res_window(_data)

    def set_rdf_window(self):
        try:
            self.rdf_plot.vb.setRange(xRange=(0, self.x_max),
                                      yRange=(self.y_min_rdf, self.y_max_rdf))
        except:
            return

    def set_tr_window(self):
        try:
            self.tr_plot.vb.setRange(xRange=(0, self.x_max),
                                     yRange=(self.y_min_tr, self.y_max_tr))
        except:
            return

    def set_fit_window(self):
        try:
            self.fit_plot.vb.setRange(xRange=(0, self.x_max),
                                     yRange=(self.y_min_fit, self.y_max_fit))
        except:
            return

    def set_res_window(self, _data):
        try:
            self.y_min_res = np.min(_data['gauss_residuals'])
            self.y_max_res = np.max(_data['gauss_residuals'])
            self.res_plot.vb.setRange(yRange=(self.y_min_res, self.y_max_res))
        except:
            return

    def create_signals(self):
        self.rdf_mouse_proxy = pg.SignalProxy(self.rdf_int_layout.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved_rdf)
        self.tr_mouse_proxy = pg.SignalProxy(self.tr_int_layout.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved_tr)
        self.gauss_mouse_proxy = pg.SignalProxy(self.fit_layout.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved_fit)

        self.rdf_plot.reset_window.connect(self.set_rdf_window)
        self.tr_plot.reset_window.connect(self.set_tr_window)
        self.fit_plot.reset_window.connect(self.set_fit_window)
        self.res_plot.reset_window.connect(self.set_fit_window)

        self.fit_limits.sigRegionChanged.connect(self.update_linear_region)

    def update_linear_region(self):
        _lb, _ub = self.fit_limits.getRegion()
        _min, _ = self.deselect_lower.getRegion()
        _, _max = self.deselect_upper.getRegion()
        self.deselect_lower.setRegion((_min, _lb))
        self.deselect_upper.setRegion((_ub, _max))

    def mouse_moved_rdf(self, __evt):
        # Use different slots as mouse always in scene bounding rectangle
        # Using signal proxy turns original args into tuple
        __pos = __evt[0]
        __mousePoint = self.rdf_plot.vb.mapSceneToView(__pos)
        self.set_mouse_pos_label(__mousePoint)
        self.rdf_plot.vline.setPos(__mousePoint.x())
        self.rdf_plot.hline.setPos(__mousePoint.y())

    def mouse_moved_tr(self, __evt):
        __pos = __evt[0]
        __mousePoint = self.tr_plot.vb.mapSceneToView(__pos)
        self.set_mouse_pos_label(__mousePoint)
        self.tr_plot.vline.setPos(__mousePoint.x())
        self.tr_plot.hline.setPos(__mousePoint.y())

    def mouse_moved_fit(self, __evt):
        __pos = __evt[0]
        if self.fit_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.fit_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)
            self.fit_plot.vline.setPos(__mousePoint.x())
            self.fit_plot.hline.setPos(__mousePoint.y())
            self.fit_plot.vline.setPen((0, 135, 153), width=0.75)
            self.fit_plot.hline.setPen((0, 135, 153), width=0.75)
            self.res_plot.vline.setPen(None)
            self.res_plot.hline.setPen(None)

        elif self.res_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.res_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)
            self.res_plot.vline.setPos(__mousePoint.x())
            self.res_plot.hline.setPos(__mousePoint.y())
            self.res_plot.vline.setPen((0, 135, 153), width=0.75)
            self.res_plot.hline.setPen((0, 135, 153), width=0.75)
            self.fit_plot.vline.setPen(None)
            self.fit_plot.hline.setPen(None)

    def set_mouse_pos_label(self, pos):
        _pos_str = (f'<span style="font-size: 11pt; color:#008799">x='
                    f'{pos.x():.2f}, y={pos.y():.2f}</span'
                    )
        self.pos_label_rdf.setText(_pos_str)
        self.pos_label_tr.setText(_pos_str)
        self.pos_label_fit.setText(_pos_str)


class CustomPlotItem(pg.PlotItem):

    def __init__(self, *args, **kwargs):
        super(CustomPlotItem, self).__init__(*args, **kwargs)

        self.enableAutoRange()
        self.setMenuEnabled(False)
        self.vline = pg.InfiniteLine(angle=90, movable=False)
        self.hline = pg.InfiniteLine(angle=0, movable=False)

        self.vline.setPen(width=0)
        self.hline.setPen(width=0)

        self.addItem(self.vline, ignoreBounds=True)
        self.addItem(self.hline, ignoreBounds=True)

        self.vb.mouseDoubleClickEvent = self.mouse_double_click_event

    def mouse_double_click_event(self, __evt):
        if __evt.button() == Qt.RightButton:
            self.vb.autoRange()


class WindowedPlotItem(CustomPlotItem):

    reset_window = pyqtSignal()

    def autoBtnClicked(self):
        self.reset_window.emit()
