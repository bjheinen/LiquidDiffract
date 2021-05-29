# -*- coding: utf-8 -*-
"""
Objective and utility functions for peak fitting in LiquidDiffract
"""
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2021, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

import numpy as np
#import lmfit

def gauss_func(r, N_ab, r_ab, sigma_ab, W_ab, c_b, sqrt_2pi):
    '''
    Gaussian type function from Cormier, 2019 (pg. 1053) 
    
    Args:
        r - r-values to evaluate the function over (np array)
        N_ab* - Avg. coordination number for pair alpha-beta
        r_ab* - Avg. bond length for pair alpha-beta
        sigma_ab* - Bond length distribution
        W_ab - X-ray pair correlation weighting for pair alpha-beta
        c_b - Atomic fraction of beta in composition
        sqrt_2pi - sqrt(2*pi) for calculation
        
        *fitted parameter
    
    Returns:
        gauss_r - gauss_func(r)
    '''
    scaling_factor = ((W_ab * N_ab) / (c_b * r_ab * sigma_ab * sqrt_2pi))
    gauss_r = scaling_factor * np.exp(-(r - r_ab)**2 / (2 * sigma_ab**2))
    return gauss_r


def gauss_obj_func(params, r, data=None, **kwargs):
    '''
    lmfit objective function to evaluate an arbitrary number of gaussian-type
    peaks (gauss_func(r)) against data
    
    Args:
        
        params - lmfit.Parameters() object containing values, bounds, and
                 refinement options for parameters [N, r, s] for each peak
        
        r - r values to evaluate over
        
        data - f(r) values to fit model to (optional) 
               if None: returns model, else: returns model - data
        
        kwargs - dict containing W_ab, and c_b for each peak, along with
                 pre-computed value of sqrt(2pi) and peak_list
                 
        kwargs['peak_list'] must be a list of unique peak ids. These are used
        to determine corresponding params/kwargs for each peak. The naming 
        convention for these is 'key'+str(peak_id). i.e. for peak id 0, params
        are N0, r0, s0, etc. kwargs are W0, c_b0, etc.
    
    Returns:
        
        model (if data == None) - sum of gaussian peaks evaluated over r
        
        residual (if data != None) - difference between model and data

    '''
    # Unpack lmfit parameters into dictionary
    param_values = params.valuesdict()
    # Get list of uniqe peak indices
    peak_idx_list = kwargs['peak_idx_list']
    # Unpack sqrt2pi from kwargs
    sqrt_2pi = kwargs['sqrt_2pi']
    # Initiate list to store each individual peak evaluated over r
    model = []
    # Loop over peak_idx_list and evaluate gauss_func for each peak
    for peak_idx in peak_idx_list:
        # Get id string
        peak_idx_str = str(peak_idx)
        # Unpack args for clarity
        # (r,) N_ab, r_ab, sigma_ab, W_ab, c_b (, sqrt_2pi)
        N_ab = param_values['N'+peak_idx_str]
        r_ab = param_values['r'+peak_idx_str]
        s_ab = param_values['s'+peak_idx_str]
        W_ab = kwargs['W'+peak_idx_str]
        c_b = kwargs['c_b'+peak_idx_str]
        # r_ab and s_ab may be zero if lower-bounds are 0
        # Catch these cases here as otherwise resulting NaNs cause 
        # lmfit to fail
        if r_ab == 0:
            r_ab += 1.e-6
        if s_ab == 0:
            s_ab += 1.e-6
        peak_model = gauss_func(r, N_ab, r_ab, s_ab, W_ab, c_b, sqrt_2pi)
        if np.isnan(peak_model).any():
            print('nans here (internal) in one of the gauss peaks!!')
        model.append(peak_model)
    # Convet to np.array
    model = np.asarray(model)
    # TODO - Convolve gaussian peak models with step function
    # Sum peaks
    model = np.sum(model, axis=0)
    
    if data is None:
        return model
    else:
        return model - data
