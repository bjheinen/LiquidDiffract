# -*- coding: utf-8 -*-
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018, Benedict J Heinen"
__email__ = "benedict.heinen@gmail.com"

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import QWidget, QVBoxLayout
import pyqtgraph as pg
import numpy as np
                
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
        #self.create_items()

        #self.mouse_position_widget = MousePositionWidget()
        #self.layout.addWidget(self.mouse_position_widget)

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
        #self.data_plot.setMenuEnabled(False)
        #self.data_plot.enableAutoRange()
        
        self.bkg_corrected_plot.setLabel('bottom', text='Q (1/A)')
        self.bkg_corrected_plot.setLabel('left', text='Intensity (a.u.)')
        #self.bkg_corrected_plot.setMenuEnabled(False)
        #self.bkg_corrected_plot.enableAutoRange()

        self.pos_label = pg.LabelItem(justify='right')
        self.pg_layout.addItem(self.pos_label,col=0, row=0)
    
    def update_plots(self, _data):
        try:
            self.p1.clear()
            self.p2.clear()
            self.p3.clear()
        except AttributeError:
            pass
        self.p1 = self.data_plot.plot(x=_data['data_x'], y=_data['data_y'], pen={'color': 0.1, 'width': 1.2})
        self.p2 = self.data_plot.plot(x=_data['bkg_x'], y=_data['bkg_y_sc'], pen={'color': '#342256', 'width': 1.2, 'style': Qt.DashLine})
        self.p3 = self.bkg_corrected_plot.plot(x=_data['cor_x'], y=_data['cor_y'], pen={'color': 0.1, 'width': 1.2})
        self.data_plot.vb.autoRange()
        self.bkg_corrected_plot.vb.autoRange()
    
    def create_signals(self):
        self.mouse_proxy = pg.SignalProxy(self.pg_layout.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved)

        
    def mouse_moved(self, __evt):
        __pos = __evt[0]  ## using signal proxy turns original arguments into a tuple
        if self.data_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.data_plot.vb.mapSceneToView(__pos)
            self.pos_label.setText("<span style='font-size: 11pt; color:#008799'>x=%0.1f, y=%0.1f</span>" % (__mousePoint.x(), __mousePoint.y()))

            self.data_plot.vline.setPos(__mousePoint.x())
            self.data_plot.hline.setPos(__mousePoint.y())
                        
            self.data_plot.vline.setPen((0, 135, 153), width=0.75)
            self.data_plot.hline.setPen((0, 135, 153), width=0.75)

            self.bkg_corrected_plot.vline.setPen(None)
            self.bkg_corrected_plot.hline.setPen(None)

        elif self.bkg_corrected_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.bkg_corrected_plot.vb.mapSceneToView(__pos)
            self.pos_label.setText("<span style='font-size: 11pt; color:#008799'>x=%0.1f, y=%0.1f</span>" % (__mousePoint.x(), __mousePoint.y()))

            self.bkg_corrected_plot.vline.setPos(__mousePoint.x())
            self.bkg_corrected_plot.hline.setPos(__mousePoint.y())       

            self.bkg_corrected_plot.vline.setPen((0, 135, 153), width=0.75)
            self.bkg_corrected_plot.hline.setPen((0, 135, 153), width=0.75)

            self.data_plot.vline.setPen(None)
            self.data_plot.hline.setPen(None)


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
        self.pg_layout.addItem(self.pos_label,col=0, row=0)
    
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

        _window = len(_data['cor_x_cut'])

        self.p1 = self.data_plot.plot(x=_data['cor_x_cut'], y=_data['cor_y_cut'], pen={'color': 0.1, 'width': 1.2})
        self.p2_a = self.iq_plot.plot(x=_data['iq_x'], y=_data['int_func'], pen={'color': 0.1, 'width': 1.2})
        self.p2_b = self.iq_plot.plot(x=_data['impr_iq_x'], y=_data['impr_int_func'], pen={'color': '#342256', 'width': 1.2, 'style': Qt.DashLine})
        self.p3_a = self.fr_plot.plot(x=_data['fr_x'][:_window], y=_data['fr_y'][:_window], pen={'color': 0.1, 'width': 1.2})
        self.p3_b = self.fr_plot.plot(x=_data['impr_fr_x'][:_window], y=_data['impr_fr_y'][:_window], pen={'color': '#342256', 'width': 1.2, 'style': Qt.DashLine})
    
        self.data_plot.vb.autoRange()
        self.iq_plot.vb.autoRange()
        self.fr_plot.vb.autoRange()
    
    
    def create_signals(self):
        self.mouse_proxy = pg.SignalProxy(self.pg_layout.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved)

        
    def mouse_moved(self, __evt):
        __pos = __evt[0]  ## using signal proxy turns original arguments into a tuple
        if self.data_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.data_plot.vb.mapSceneToView(__pos)
            self.pos_label.setText("<span style='font-size: 11pt; color:#008799'>x=%0.1f, y=%0.1f</span>" % (__mousePoint.x(), __mousePoint.y()))

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
            self.pos_label.setText("<span style='font-size: 11pt; color:#008799'>x=%0.1f, y=%0.1f</span>" % (__mousePoint.x(), __mousePoint.y()))

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
            self.pos_label.setText("<span style='font-size: 11pt; color:#008799'>x=%0.1f, y=%0.1f</span>" % (__mousePoint.x(), __mousePoint.y()))

            self.fr_plot.vline.setPos(__mousePoint.x())
            self.fr_plot.hline.setPos(__mousePoint.y())       

            self.fr_plot.vline.setPen((0, 135, 153), width=0.75)
            self.fr_plot.hline.setPen((0, 135, 153), width=0.75)

            self.data_plot.vline.setPen(None)
            self.data_plot.hline.setPen(None)
            self.iq_plot.hline.setPen(None)
            self.iq_plot.vline.setPen(None)

#############################################################################
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
        self.pg_layout.addItem(self.pos_label,col=0, row=0)
    
    def update_plots(self, _data):
        try:
            self.p1.clear()
            self.p2.clear()
            self.p3.clear()
        except AttributeError:
            pass

        # Determine data window for length of sq (pre fft)
        _window = len(_data['sq_x'])

        self.p1 = self.sq_plot.plot(x=_data['sq_x'], y=_data['sq_y'], pen={'color': 0.1, 'width': 1.2})
        self.p2 = self.gr_plot.plot(x=_data['gr_x'][:_window], y=_data['gr_y'][:_window], pen={'color': 0.1, 'width': 1.2})
        self.p3 = self.rdf_plot.plot(x=_data['rdf_x'][:_window], y=_data['rdf_y'][:_window], pen={'color': 0.1, 'width': 1.2})
            
        
        self.x_max_gr = 12
        _gr_cut = np.nan_to_num(_data['gr_y'][np.where(_data['gr_x']<self.x_max_gr)])
        self.y_min_gr = np.min(_gr_cut)
        self.y_max_gr = np.max(_gr_cut)
        
        self.x_max_rdf = 8
        _rdf_cut =  np.nan_to_num(_data['rdf_y'][np.where(_data['rdf_x']<self.x_max_rdf)])
        self.y_min_rdf = np.min(_rdf_cut)
        self.y_max_rdf = np.max(_rdf_cut)
        
        self.set_gr_window()
        self.set_rdf_window()
        
        self.sq_plot.vb.autoRange()

    def set_gr_window(self):
        try:
            self.gr_plot.vb.setRange(xRange=(0, self.x_max_gr), 
                                     yRange=(self.y_min_gr, self.y_max_gr))
        except:
            return
    
    def set_rdf_window(self):
        try:
            self.rdf_plot.vb.setRange(xRange=(0, self.x_max_rdf), 
                                      yRange=(self.y_min_rdf, self.y_max_rdf))
        except:
            return

    def create_signals(self):
        self.mouse_proxy = pg.SignalProxy(self.pg_layout.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved)
        self.gr_plot.reset_window.connect(self.set_gr_window)
        self.rdf_plot.reset_window.connect(self.set_rdf_window)

        
    def mouse_moved(self, __evt):
        __pos = __evt[0]  ## using signal proxy turns original arguments into a tuple
        if self.sq_plot.sceneBoundingRect().contains(__pos):
            __mousePoint = self.sq_plot.vb.mapSceneToView(__pos)
            self.pos_label.setText("<span style='font-size: 11pt; color:#008799'>x=%0.1f, y=%0.1f</span>" % (__mousePoint.x(), __mousePoint.y()))

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
            self.pos_label.setText("<span style='font-size: 11pt; color:#008799'>x=%0.1f, y=%0.1f</span>" % (__mousePoint.x(), __mousePoint.y()))

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
            self.pos_label.setText("<span style='font-size: 11pt; color:#008799'>x=%0.1f, y=%0.1f</span>" % (__mousePoint.x(), __mousePoint.y()))

            self.rdf_plot.vline.setPos(__mousePoint.x())
            self.rdf_plot.hline.setPos(__mousePoint.y())       

            self.rdf_plot.vline.setPen((0, 135, 153), width=0.75)
            self.rdf_plot.hline.setPen((0, 135, 153), width=0.75)

            self.sq_plot.vline.setPen(None)
            self.sq_plot.hline.setPen(None)
            self.gr_plot.hline.setPen(None)
            self.gr_plot.vline.setPen(None)





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
            #self.vb.enableAutoRange()
            #self._auto_range = True
            #self.vb.sigRangeChangedManually.emit(self.vb.state['mouseEnabled'])
            
            
class WindowedPlotItem(CustomPlotItem):

    reset_window = pyqtSignal()
    
    def autoBtnClicked(self):
        self.reset_window.emit()