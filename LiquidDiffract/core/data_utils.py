# -*- coding: utf-8 -*-
"""
Utility functions for common data operations in LiquidDiffract
"""
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018-2024, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

from functools import lru_cache, wraps
import numpy as np
import scipy.interpolate
import scipy.optimize
import scipy.signal

def rebin_data(x, y, dx=0.02, x_lim=None, extrapolate_mode='fill'):
    '''
    Re-bins x,y data (using interpolation) to achieve a set step size of x (dx)
    Extrapolates x < x_min to zero
    '''
    if x_lim != None:
        x_rebin = np.arange(x_lim[0], x_lim[1], dx)
    else:
        x_rebin = np.arange(0, x[-1], dx)
    # Check for any nans in y and remove via interpolation
    if np.isnan(np.sum(y)):
        y = interp_nan(y)
    if extrapolate_mode == 'fill':
        f_interp = scipy.interpolate.interp1d(x, y,
                                              kind='cubic',
                                              fill_value=y[0],
                                              bounds_error=False)
    elif extrapolate_mode == 'extrapolate':
        f_interp = scipy.interpolate.interp1d(x, y,
                                              kind='cubic',
                                              fill_value='extrapolate')
    else:
        raise NotImplementedError('Extrapolate mode not valid')
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


def convert_q_space(Q, wavelength):
    '''Convert Q-space values to 2 theta (in radians)'''
    return 2 * np.arcsin(wavelength * Q / (4*np.pi))


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
    else:
        raise NotImplementedError('Smooth method not implemented, please check keyword argument')


def clean_base_indices(lefts, rights):
    '''
    Clean peak base limits
    See https://github.com/scipy/scipy/issues/19232
    '''
    _lefts = np.copy(lefts)
    _rights = np.copy(rights)
    for i in range(len(lefts)-1):
        if lefts[i] == lefts[i+1]:
            _lefts[i+1] = rights[i]
        if rights[i] == rights[i+1]:
            _rights[i] = lefts[i+1]
    return _lefts, _rights


def find_integration_limits(r, rdf, rho=None, peak_search_limit=10.0, search_method=None):
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

    # Get index of first peak after r0 (last RDF = 0)
    peak_idx_first = np.argmax(peak_list>brent_a)

    # Get peak prominences and positions at the bases either side
    peak_prom, left_bases, right_bases = scipy.signal.peak_prominences(rdf, peak_list)
    # Check if left bases merge peaks and fix if necessary
    if len(np.unique(left_bases)) < len(left_bases):
        left_bases, right_bases = clean_base_indices(left_bases, right_bases)
    # Get index of most prominent peak
    peak_idx_prom = np.argmax(peak_prom)
    # Could check relative prominence here? e.g. np.diff(peak_prom)/peak_prom[0:-1]
    # Currently, simply assume most prominent peak is the right target. It should
    # be the first peak with significant prominence.

    # If the most prominent peak is not the first after last RDF = 0 it
    # is usually because small oscillation at RDF ~ 0 need to be discounted
    if peak_idx_prom > peak_idx_first:
        # Take most prominent peak as position to use
        peak_idx = peak_idx_prom
        # By definition peak_idx here must be >= 1
        # Re-define r_0 as minimum to left of peak to ignore prior oscillations
        # Find minimum in Tr as sharper
        r_0 = scipy.optimize.minimize_scalar(Tr_interp,
                                            method='bounded',
                                            bounds=(r[peak_list[peak_idx-1]],
                                                    r[peak_list[peak_idx]])).x

    # If the most prominent peak occurs before the 'first' peak after the last
    # RDF = 0 it is usually because the right base of the prominent peak dips
    # below RDF = 0. The left base may also be at RDF < 0.
    elif peak_idx_prom < peak_idx_first:
        # Take the most prominent peak as position to use
        peak_idx = peak_idx_prom
        # Re-define r_0
        # Find minimum preceding peak
        # Set upper bound as peak position
        # (-1 r-step to account for exact peak position)
        r_0_ub = r[peak_list[peak_idx]-1]
        if peak_idx > 0:
            # If preceding peaks have been found take the max of the next one as
            # lower bound to find r_0 (+1 r-step to account for exact position)
            r_0_lb = r[peak_list[peak_idx-1]+1]
        else:
            # If no previous peaks found take the last change in sign of gradient
            # (-1 r-step again)
            r_0_lb = r[np.argwhere(np.sign(np.diff(Tr[np.where(r <= r_0_ub)])) == -1)[-2]]
        # Find minimum
        preceding_minimum = scipy.optimize.minimize_scalar(Tr_interp,
                                                           method='bounded',
                                                           bounds=(r_0_lb, r_0_ub)).x
        # Check if minimum at RDF < 0
        if rdf_interp(preceding_minimum) < 0:
            # Redefine r_0 as last RDF = 0 crosing before peak
            # First find last sign-changing interval in discrete data from
            # preceding minimum to peak position
            check_interval = np.ravel(np.where((r >= preceding_minimum) & (r <= r_0_ub)))
            # Pad check_interval if len < 3
            if len(check_interval) < 3:
                check_interval = np.pad(check_interval, 1, 'edge')
                check_interval[0] -= 1
                check_interval[-1] += 1
            brent_a = check_interval[np.argwhere(np.sign(rdf[check_interval]) == -1)[-1]]
            # Find root using scipy.optimize.brentq
            r_0 = scipy.optimize.brentq(rdf_interp, r[brent_a], r[brent_a+1])
        else:
            r_0 = preceding_minimum

    else:
        # If peak idx match use located r_0
        peak_idx = peak_idx_first
    
    # Get left and right base of peak
    peak_bases = (r[left_bases[peak_idx]], r[right_bases[peak_idx]])
    # Refine rp_max (peak centre in r g(r))
    rp_max = scipy.optimize.minimize_scalar(lambda x: -Tr_interp(x), method='bounded', bounds=peak_bases).x
    # Refine r_max (peak centre in r^2 g(r))
    r_max = scipy.optimize.minimize_scalar(lambda x: -rdf_interp(x), method='bounded', bounds=peak_bases).x

    # Get approximate position of next peak to provide upper bound on r_min
    try:
        next_peak = r[peak_list[peak_idx+1]]
    # Handle IndexError if chosen peak is last in peak_list
    except IndexError:
        # Re-do peak search in g(r) if density (rho) available
        # This seems useful in rare cases, but usually the search_limit needs to be increased
        if rho != None:
            with np.errstate(divide='ignore', invalid='ignore'):
                gr = rdf/(4*np.pi*rho*r**2)
            gr_peak_list, _ = scipy.signal.find_peaks(gr[r<peak_search_limit])
            # Take the next peak after peak_list[peak_idx]
            try:
                next_peak = r[np.where((r[gr_peak_list] > r_max) & (r[gr_peak_list] > rp_max))][0]
            except IndexError:
                # Using peak_search_limit is non-ideal
                # Consider warning here to suggest increasing peak_search_limit
                next_peak = peak_search_limit
        # If rho not available use
        else:
            # Using peak_search_limit is non-ideal
            # Consider warning here to suggest increasing peak_search_limit
            # Programmatically increasing psl may lead to recursion issues
            next_peak = peak_search_limit

    # Refine r_min (position of 1st minimum after 1st peak)
    # Note: r_min should be global minimum but optimisation may return local min
    r_min = scipy.optimize.minimize_scalar(rdf_interp, method='bounded', bounds=(r_max, next_peak)).x

    # Check r_min is not at RDF < 0 in case where peak precedes final crossing of RDF = 0
    if (peak_idx_prom < peak_idx_first) and (rdf_interp(r_min) < 0):
        # Re-define r_min to preceding crossing of RDF = 0
        # Find first sign changing interval between peak and r_min
        check_interval = np.where((r >= rp_max) & (r <= r_min))
        brent_b = np.ravel(check_interval)[np.argwhere(np.sign(rdf[check_interval]) == -1)[0]]
        # Find root using scipy.optimize.brentq
        r_min = scipy.optimize.brentq(rdf_interp, r[brent_b-1], r[brent_b])

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
