# -*- coding: utf-8 -*-
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2018, Benedict J. Heinen"
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
        self.fr_plot = CustomPlotItem()

        self.data_plot.plot(x=[], y=[])
        self.iq_plot.plot(x=[], y=[])
        self.fr_plot.plot(x=[], y=[], pen={})

        self.pg_layout.addItem(self.data_plot, row=1, col=0)
        self.pg_layout.addItem(self.iq_plot, row=2, col=0)
        self.pg_layout.addItem(self.fr_plot, row=3, col=0)

        self.pg_layout_widget.addItem(self.pg_layout)

        self.layout.addWidget(self.pg_layout_widget)

    def style_plots(self):
        self.data_plot.setLabel('bottom', text='Q (1/A)')
        self.data_plot.setLabel('left', text='Intensity (a.u.)')

        self.iq_plot.setLabel('bottom', text='Q (1/A)')
        self.iq_plot.setLabel('left', text='i(Q)')

        self.fr_plot.setLabel('bottom', text='r (A)')
        self.fr_plot.setLabel('left', text='F(r)')

        self.pos_label = pg.LabelItem(justify='right')
        self.pg_layout.addItem(self.pos_label, col=0, row=0)

    def update_plots(self, _data):
        try:
            self.p1.clear()
            self.p2_a.clear()
            self.p3_a.clear()
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
        if np.isnan(_data['fr_x']).any():
             _data['fr_x'] = data_utils.interp_nan(_data['fr_x'])
        if np.isnan(_data['impr_fr_x']).any():
             _data['impr_fr_x'] = data_utils.interp_nan(_data['impr_fr_x'])

        # For the F(r) the r step = pi/q_max. Because the data is padded this
        # q_max is larger than the original q_max. i.e. q_max = dq * 2**N/2
        # Then: r_max = 2**N/2 * dr = pi/(dq * 2**N/2) * 2**N/2 = pi/dq
        # Displaying the full length data is not necessary as most of the 
        # useful information is contained at fairly low r, and the length of r
        # is controlled only by sampling frequency in Q-space. At high r-values
        # the F(r) is dominated by ripples from the truncated integral (0-qmax)
        # The max value of r with 'real' q resolution is 1/dq.
        _window = 0
        try:
            _dq = _data['iq_x'][1] - _data['iq_x'][0]
            try:
                _window = np.argmax(_data['fr_x'] >= (1/_dq))
            except ValueError:
                 pass
        except IndexError:
            pass

        self.p1 = self.data_plot.plot(x=_data['cor_x_cut'], y=_data['cor_y_cut'], pen={'color': 0.1, 'width': 1.2})
        self.p2_a = self.iq_plot.plot(x=_data['iq_x'], y=_data['int_func'], pen={'color': 0.1, 'width': 1.2})
        self.p2_b = self.iq_plot.plot(x=_data['impr_iq_x'], y=_data['impr_int_func'], pen={'color': '#342256', 'width': 1.2, 'style': Qt.DashLine})
        self.p3_a = self.fr_plot.plot(x=_data['fr_x'][:_window], y=_data['fr_y'][:_window], pen={'color': 0.1, 'width': 1.2})
        self.p3_b = self.fr_plot.plot(x=_data['impr_fr_x'][:_window], y=_data['impr_fr_y'][:_window], pen={'color': '#342256', 'width': 1.2, 'style': Qt.DashLine})

        if _data['mod_func'] != 'None':
            # Apply modification function and handle any nans that arise 
            _modified_int_func = _data['modification'] * _data['int_func']
            if np.isnan(_modified_int_func).any():
                 _modified_int_func = data_utils.interp_nan(_modified_int_func)
            self.p2_c = self.iq_plot.plot(x=_data['iq_x'], y=_modified_int_func, pen={'color': '#342256', 'width': 0.8, 'style': Qt.DashLine})

        self.data_plot.vb.autoRange()
        self.iq_plot.vb.autoRange()
        self.fr_plot.vb.autoRange()

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
            self.fr_plot.vline.setPen(None)
            self.fr_plot.hline.setPen(None)

        elif self.iq_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.iq_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)

            self.iq_plot.vline.setPos(__mousePoint.x())
            self.iq_plot.hline.setPos(__mousePoint.y())

            self.iq_plot.vline.setPen((0, 135, 153), width=0.75)
            self.iq_plot.hline.setPen((0, 135, 153), width=0.75)

            self.data_plot.vline.setPen(None)
            self.data_plot.hline.setPen(None)
            self.fr_plot.vline.setPen(None)
            self.fr_plot.hline.setPen(None)

        elif self.fr_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.fr_plot.vb.mapSceneToView(__pos)
            self.set_mouse_pos_label(__mousePoint)

            self.fr_plot.vline.setPos(__mousePoint.x())
            self.fr_plot.hline.setPos(__mousePoint.y())

            self.fr_plot.vline.setPen((0, 135, 153), width=0.75)
            self.fr_plot.hline.setPen((0, 135, 153), width=0.75)

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

        # Some versions of pyqtgraph cannot produce plot if nan values present
        # Fix nan values by interpolation
        if np.isnan(_data['sq_y']).any():
             _data['sq_y'] = data_utils.interp_nan(_data['sq_y'])
        if np.isnan(_data['gr_y']).any():
             _data['gr_y'] = data_utils.interp_nan(_data['gr_y'])
        if np.isnan(_data['rdf_y']).any():
             _data['rdf_y'] = data_utils.interp_nan(_data['rdf_y'])
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

        # Seperate widgets for rdf, tr

        self.pg_layout_widget_rdf = pg.GraphicsLayoutWidget()
        self.pg_layout_widget_rdf.setContentsMargins(0, 0, 0, 0)
        self.pg_layout_rdf = pg.GraphicsLayout()
        self.pg_layout_rdf.setContentsMargins(0, 0, 0, 0)

        self.pg_layout_widget_tr = pg.GraphicsLayoutWidget()
        self.pg_layout_widget_tr.setContentsMargins(0, 0, 0, 0)
        self.pg_layout_tr = pg.GraphicsLayout()
        self.pg_layout_tr.setContentsMargins(0, 0, 0, 0)        

        # RDF(r) plot
        self.rdf_plot = WindowedPlotItem()
        self.rdf_plot.plot(x=[], y=[])
        self.pg_layout_rdf.addItem(self.rdf_plot, row=1, col=0)
        self.pg_layout_widget_rdf.addItem(self.pg_layout_rdf)

        # T(r) plot
        self.tr_plot = WindowedPlotItem()
        self.tr_plot.plot(x=[], y=[])
        self.pg_layout_tr.addItem(self.tr_plot, row=1, col=0)
        self.pg_layout_widget_tr.addItem(self.pg_layout_tr)

        # Set RDF(r) plot as default view
        self.layout.addWidget(self.pg_layout_widget_rdf)
        self.layout.addWidget(self.pg_layout_widget_tr)

        self.pg_layout_widget_rdf.setVisible(True)
        self.pg_layout_widget_tr.setVisible(False)

    def style_plots(self):

        self.rdf_plot.setLabel('bottom', text='r (A)')
        self.rdf_plot.setLabel('left', text='RDF(r)')

        self.tr_plot.setLabel('bottom', text='r (A)')
        self.tr_plot.setLabel('left', text='T(r)')

        self.rdf_xaxis = pg.InfiniteLine(pos=0, angle=0, movable=False, pen={'color': 'k', 'width': 0.75})
        self.rdf_plot.addItem(self.rdf_xaxis)
        self.tr_xaxis = pg.InfiniteLine(pos=0, angle=0, movable=False, pen={'color': 'k', 'width': 0.75})
        self.tr_plot.addItem(self.tr_xaxis)

        # Add position label to all plot views
        self.pos_label_rdf = pg.LabelItem(justify='right')
        self.pg_layout_rdf.addItem(self.pos_label_rdf, col=0, row=0)

        self.pos_label_tr = pg.LabelItem(justify='right')
        self.pg_layout_tr.addItem(self.pos_label_tr, col=0, row=0)

        self.rdf_plot.vline.setPen((0, 135, 153), width=0.75)
        self.rdf_plot.hline.setPen((0, 135, 153), width=0.75)

        self.tr_plot.vline.setPen((0, 135, 153), width=0.75)
        self.tr_plot.hline.setPen((0, 135, 153), width=0.75)


    def clear_plots(self, _clear_all=False):
        try:
            self.p_rdf.clear()
            self.p_tr.clear()            
        except AttributeError:
            pass      

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

        self.clear_plots()

        # Some versions of pyqtgraph cannot produce plot if nan values present
        # Fix nan values by interpolation
        if np.isnan(_data['rdf_y']).any():
             _data['rdf_y'] = data_utils.interp_nan(_data['rdf_y'])
        if np.isnan(_data['tr_y']).any():
             _data['tr_y'] = data_utils.interp_nan(_data['tr_y'])
        if np.isnan(_data['fr_y']).any():
             _data['fr_y'] = data_utils.interp_nan(_data['fr_y'])

        # Interpolate data for smoother plots
        #_data['rdf_x'], _data['rdf_y'] = data_utils.rebin_data(_data['rdf_x'], _data['rdf_y'], dx=0.01)
        #_data['tr_x'], _data['tr_y'] = data_utils.rebin_data(_data['tr_x'], _data['tr_y'], dx=0.01)

        self.p_rdf = self.rdf_plot.plot(x=_data['rdf_x'], y=_data['rdf_y'], pen={'color': 0.1, 'width': 1.2})
        self.p_tr = self.tr_plot.plot(x=_data['tr_x'], y=_data['tr_y'], pen={'color': 0.1, 'width': 1.2})
        #self.p_fr = self.rdf_plot.plot(x=_data['rdf_x'][:_window], y=_data['rdf_y'][:_window], pen={'color': 0.1, 'width': 1.2})

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

        # Plot fill areas for integrals
        # Nc is plotted first/behind as it extends beyond Nb
        # Only plot if r0 < r
        if _data['r0'] < _data['rmin']:
            self.p_Nc = self.rdf_plot.plot(x=_Nc_area_x, y=_Nc_area_y, brush=(26,121,199), pen=None, fillLevel=0)
        if _data['r0'] < _data['rpmax']:
            self.p_Na = self.tr_plot.plot(x=_Na_area_x, y=_Na_area_y, brush=(199,26,74), pen=None, fillLevel=0)
        if _data['r0'] < _data['rmax']:
            self.p_Nb = self.rdf_plot.plot(x=_Nb_area_x, y=_Nb_area_y, brush=(199,26,74), pen=None, fillLevel=0)


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


    def update_plot_windows(self, _data):
        _rdf_cut = np.nan_to_num(_data['rdf_y'][np.where(_data['rdf_x'] < self.x_max)])
        self.y_min_rdf = np.min(_rdf_cut)
        self.y_max_rdf = np.max(_rdf_cut)

        _tr_cut = np.nan_to_num(_data['tr_y'][np.where(_data['tr_x'] < self.x_max)])
        self.y_min_tr = np.min(_tr_cut)
        self.y_max_tr = np.max(_tr_cut)

        self.set_rdf_window()
        self.set_tr_window()       


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


    def create_signals(self):
        self.rdf_mouse_proxy = pg.SignalProxy(self.pg_layout_rdf.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved_rdf)
        self.tr_mouse_proxy = pg.SignalProxy(self.pg_layout_tr.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved_tr)

        self.rdf_plot.reset_window.connect(self.set_rdf_window)

    def mouse_moved_rdf(self, __evt):
        # Use different slots as mouse always in scene bounding rectangle
        # Using signal proxy turns original args into tuple

        __pos = __evt[0]
        #if self.rdf_plot.sceneBoundingRect().contains(__pos):
        __mousePoint = self.rdf_plot.vb.mapSceneToView(__pos)
        self.set_mouse_pos_label(__mousePoint)

        self.rdf_plot.vline.setPos(__mousePoint.x())
        self.rdf_plot.hline.setPos(__mousePoint.y())


    def mouse_moved_tr(self, __evt):
        __pos = __evt[0]
        #if self.tr_plot.sceneBoundingRect().contains(__pos):
        __mousePoint = self.tr_plot.vb.mapSceneToView(__pos)
        self.set_mouse_pos_label(__mousePoint)

        self.tr_plot.vline.setPos(__mousePoint.x())
        self.tr_plot.hline.setPos(__mousePoint.y())


    def set_mouse_pos_label(self, pos):
        _pos_str = (f'<span style="font-size: 11pt; color:#008799">x='
                    f'{pos.x():.2f}, y={pos.y():.2f}</span'
                    )
        self.pos_label_rdf.setText(_pos_str)
        self.pos_label_tr.setText(_pos_str)


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
            # self.vb.enableAutoRange()
            # self._auto_range = True
            # self.vb.sigRangeChangedManually.emit(self.vb.state['mouseEnabled'])


class WindowedPlotItem(CustomPlotItem):

    reset_window = pyqtSignal()

    def autoBtnClicked(self):
        self.reset_window.emit()
