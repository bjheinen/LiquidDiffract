# -*- coding: utf-8 -*-
"""
Accessory functions for data manipulation
"""
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018-2019, Benedict J Heinen"
__email__ = "benedict.heinen@gmail.com"

import numpy as np
import scipy.interpolate
from scipy.signal import savgol_filter


def rebin_data(x, y, dx=0.02):
    '''
    Rebins x,y data via interpolation to a set step size (dx).
    Also extrapolates to zero.
    '''
    x_rebin = np.arange(0, x[-1], dx)
    f_interp = scipy.interpolate.interp1d(x, y, kind='cubic', fill_value='extrapolate')
    y_rebin = f_interp(x_rebin)
    return x_rebin, y_rebin
    
def convert_two_theta(two_theta, wavelength):
    '''
    Convert 2theta data to Q space. The X-ray wavelength must be known
    '''
    q_data = (4 * np.pi / wavelength) * np.sin(np.radians(two_theta)/2)
    return q_data

def zero_norm(int_func, _S_inf):
    '''
    Shift data to fit first value to zero - S_inf
    '''
    shift = int_func[0] 
    int_func = int_func - shift
    return int_func

def bkg_scaling_residual(bkg_scaling, *args):
    '''
    Returns the mean deviation between data and a scaled background
    '''
    (data, bkg) = args
    return np.mean(np.abs(data - (bkg*bkg_scaling)))


def smooth_data(data, method='savitzky-golay', window_length=31, poly_order=3):
    '''
    Applies a savitzky-golay filter to smooth noise from x, y data.
    
    Arguments:
        
    'data' must be a 1d array of y data only.
    
    'window_length' is the length of the filter window (i.e. the number of coefficients). 
    This must be a positive odd integer.
    
    'poly_order' is the order of the polynomial used to fit the samples. 
    This must be less than window_length.
    '''
    if method == 'savitzky-golay':
        return savgol_filter(data, window_length, poly_order)
