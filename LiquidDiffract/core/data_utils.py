# -*- coding: utf-8 -*-
"""
Utility functions for common data operations in LiquidDiffract
"""
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018-2021, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

from functools import lru_cache, wraps
import numpy as np
import scipy.interpolate
import scipy.optimize
import scipy.signal

def rebin_data(x, y, dx=0.02, x_lim=None):
    '''
    Re-bins x,y data (using interpolation) to achieve a set step size of x (dx)
    Extrapolates x < x_min to zero
    '''
    if x_lim != None:
        x_rebin = np.arange(x_lim[0], x_lim[1], dx)
    else:
        x_rebin = np.arange(0, x[-1], dx)
    f_interp = scipy.interpolate.interp1d(x, y,
                                          kind='cubic',
                                          fill_value='extrapolate')
    y_rebin = f_interp(x_rebin)
    return x_rebin, y_rebin

def interp_data(x, y, x2):
    '''
    Helper function to return interpolated y-value at given x2 for f(x)=y
    '''
    f_interp = scipy.interpolate.interp1d(x, y, kind='cubic',
                                          fill_value='extrapolate')
    return f_interp(x2).flatten()


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
        return scipy.signal.savgol_filter(data, window_length, poly_order)


def find_integration_limits(r, rdf, rho=None, peak_search_limit=10.0, search_method='first'):
    '''
    Computes the integration limits used in the calculation of the 
    first coordination number for a monatomic system.

    The function returns several values for use in different methods of
    calculating the coordination number. 
    See 'core.integrate_coordination_sphere' for more information on
    the methods implemented in LiquidDiffract.

    The location of the peak corresponding to the first coordination sphere 
    (the 1st peak in RDF(r)) is first approximated from the discrete data 
    using a peak finding algorithm. A continuous function is then generated
    using cubic interpolation and the position of r_0 refined using a root
    finding algorithm and the positions of r_max, rp_max, and r_min are 
    refined through calls to scipy.minimize

    Args:
        r   - r values (numpy array)
        rdf - Corresponding values of the radial distribution function, RDF(r) (numpy array)

        rho - optional kwarg to provide number density (rho). This is used to
              calculate g(r) if the peak search in rdf(r) fails

        peak_search_limit - optional kwarg to set limit on intial peak search
                            this is useful as the most prominent (topographic)
                            peak is taken and the function can be confused by
                            choosing a peak at high r. If the main peak is the
                            nth peak, there must be at least n+1 peaks found

    Returns:
        r_0 - Leading edge of the 1st peak in RDF(r), taken as the nearest root
        rp_max - Peak centre in T(r) (equivalent to centre in r*g(r))
        r_max - Peak centre in RDF(r) (equivalent to centre in r^2*g(r))
        r_min - Position of 1st minimum after the 1st peak in RDF(r)
    '''
    # Calculate T(r)
    with np.errstate(divide='ignore', invalid='ignore'):
        Tr = np.nan_to_num(rdf/r)

    # Create continuous functions of RDF(r) and T(r) via cubic interpolation
    rdf_interp = scipy.interpolate.interp1d(r, rdf, kind='cubic', fill_value='extrapolate')
    Tr_interp = scipy.interpolate.interp1d(r, Tr, kind='cubic', fill_value='extrapolate')

    # Take the last root in RDF(r) as the leading edge of the 1st peak, r_0
    # Uses Brent method to find a zero of RDF(r) on the sign changing interval [a , b].
    # First find last sign-changing interval in discrete data
    brent_a = np.argwhere(np.sign(rdf) == -1)[-1]
    # Find root using scipy.optimize.brentq
    r_0 = scipy.optimize.brentq(rdf_interp, r[brent_a], r[brent_a+1])

    # Generate list of indices of approximate peak positions in the discrete 
    # data. Here using scipy.signal but using np.diff(np.sign(np.diff(rdf)) to
    # avoid the extra import may be less expensive
    peak_list, _ = scipy.signal.find_peaks(rdf[r<peak_search_limit])

    # Select 1st peak after r0 as 1st coordination sphere
    peak_idx_first = np.argmax(peak_list>brent_a)

    # Select most prominent peak as 1st coordination sphere and get positions
    # at the base either side
    peak_prom, left_bases, right_bases = scipy.signal.peak_prominences(rdf, peak_list)
    peak_idx_prom = np.argmax(peak_prom)
    
    if search_method == 'first':
        peak_idx = peak_idx_first
    elif search_method == 'prominent':
        peak_idx = peak_idx_prom
    else:
        raise AttributeError('\'method\' must be \'first\' or \'prominent\'')
    
    peak_bases = (r[left_bases[peak_idx]], r[right_bases[peak_idx]])

    # Get approximate position of next peak to provide upper bound on r_min
    try:
        next_peak = r[peak_list[peak_idx+1]]
    # Handle IndexError if chosen peak is last in peak_list
    except IndexError:
        # Re-do peak search in g(r) if density (rho) available
        if rho != None:
            with np.errstate(divide='ignore', invalid='ignore'):
                gr = rdf/(4*np.pi*rho*r**2)
            gr_peak_list, _ = scipy.signal.find_peak(gr[r<peak_search_limit])
            # Take the first peak after r0
            next_peak = r[gr_peak_list[np.argmax(gr_peak_list>brent_a)+1]]
        # If rho not available use
        else:
            next_peak = peak_search_limit

    # Refine rp_max (peak centre in r g(r))
    rp_max = scipy.optimize.minimize_scalar(lambda x: -Tr_interp(x), method='bounded', bounds=peak_bases).x
    # Refine r_max (peak centre in r^2 g(r))
    r_max = scipy.optimize.minimize_scalar(lambda x: -rdf_interp(x), method='bounded', bounds=peak_bases).x

    # Refine r_min (position of 1st minimum after 1st peak)
    # Note: r_min should be global minimum but optimisation may return local min
    r_min = scipy.optimize.minimize_scalar(rdf_interp, method='bounded', bounds=(r_max, next_peak)).x

    return r_0, rp_max, r_max, r_min


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