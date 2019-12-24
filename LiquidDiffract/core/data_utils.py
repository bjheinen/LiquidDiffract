# -*- coding: utf-8 -*-
"""
Utility functions for common data operations
"""
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018-2019, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

from functools import lru_cache, wraps
import numpy as np
import scipy.interpolate
from scipy.signal import savgol_filter


def rebin_data(x, y, dx=0.02):
    '''
    Re-bins x,y data (using interpolation) to achieve a set step size of x (dx)
    Extrapolates x < x_min to zero
    '''
    x_rebin = np.arange(0, x[-1], dx)
    f_interp = scipy.interpolate.interp1d(x, y,
                                          kind='cubic',
                                          fill_value='extrapolate')
    y_rebin = f_interp(x_rebin)
    return x_rebin, y_rebin


def convert_two_theta(two_theta, wavelength):
    '''Convert 2theta values to Q-space values at a known wavelength'''
    q_data = (4 * np.pi / wavelength) * np.sin(np.radians(two_theta)/2)
    return q_data


def zero_norm(y, shift=None):
    '''Shift data on y axis. If no 'shift' given, first value is set to zero'''
    if shift == None:
         shift = y[0]
    return y - shift


def interp_nan(y):
     '''
     Uses interpolation to fix nan values in an array
     This is useful for nan values at start of interference function array
     which cannot be set to zero with np.nan_to_num
     '''
     # Array of True/False values for nan location
     nans = np.isnan(y)
     # Return nonzero (non-nan) values
     f = lambda z: z.nonzero()[0]
     y[nans] = np.interp(f(nans), f(~nans), y[~nans])
     return y


def bkg_scaling_residual(bkg_scaling, *args):
    '''
    Returns the mean deviation between data and a scaled background

    Uses *args syntax so the function can be passed to a solver

    It is expected that *args == (data, bkg)

    Assumes data/bkg x-axis values are the same

    Arguments:

    bkg_scaling -- Background scaling factor

    data -- y values only

    bkg -- unscaled background data (y values only)

    '''
    (data, bkg) = args
    return np.mean(np.abs(data - (bkg*bkg_scaling)))


def smooth_data(data, method='savitzky-golay', window_length=31, poly_order=3):
    '''
    Applies a savitzky-golay filter to smooth noise from x, y data

    Arguments:

    data -- Must be a 1d array of y data only.

    window_length -- The length of the filter window (number of coefficients)
                     Must be a positive odd integer

    poly_order -- The order of the polynomial used to fit the samples
                  This must be an integer less than window_length.
    '''
    if method == 'savitzky-golay':
        return savgol_filter(data, window_length, poly_order)


def data_cache(*args, **kwargs):
    '''
    Provides a decorator to cache the result of functions of the form
    f(hashable_type, 1d_Array). The numpy array is converted to a tuple so
    that functools.lru_cache can be used.

    Here it is used to cache the computation of element/Q specific form-
    factors and compton scattering.
    '''
    def decorator(function):
        @wraps(function)
        def wrapper(hashable_arg, np_array):
            hashable_array = tuple(np_array)
            return cached_wrapper(hashable_arg, hashable_array)

        @lru_cache(*args, **kwargs)
        def cached_wrapper(hashable_arg, hashable_array):
            array = np.array(hashable_array)
            return function(hashable_arg, array)

        wrapper.cache_info = cached_wrapper.cache_info
        wrapper.cache_clear = cached_wrapper.cache_clear

        return wrapper
    return decorator