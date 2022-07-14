#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This script computes the function X^2(n, rho).
It does a simple brute force calculation of X^2 across the range of (n, rho),
where n is the number of iterations in the Eggert procedure and rho is the
density.

'''
__author__ = 'Benedict J. Heinen'
__copyright__ = 'Copyright 2019, Benedict J. Heinen'
__license__ = 'GNU GPL v3'
__email__ = 'benedict.heinen@gmail.com'


import numpy as np
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import LiquidDiffract.core.core as liquid
from LiquidDiffract.core.data_utils import rebin_data
import time

def update_progress_bar(iteration_number, total_iterations, elapsed_time, estimated_total_time):
    '''
    Prints a progress bar, with % completed and estimated time remaining, to the terminal.
    '''
    percentage_done = 100 * (iteration_number / total_iterations)
    if iteration_number == total_iterations: 
        percentage_done = 100
    percent_str = ('{0:.1f}').format(percentage_done)
    # Char length of progress bar = 50
    bar_length = 50
    fill_char = '█'
    fill_length = bar_length * iteration_number // total_iterations
    progress_bar_str = fill_char * fill_length + '-' * (bar_length - fill_length)
    eta = estimated_total_time - elapsed_time
    if eta < 60:
        eta_str = str(int(eta)) + 's' + ' '*12
    else:
        eta_str = time.strftime('%H:%M:%S', time.gmtime(eta))
    if iteration_number == total_iterations:
        eta_str = '0s'
    if total_iterations == 1:
        return
    print('\rProgress: |%s| %s%% Complete  |  ETA: %s' % (progress_bar_str, percent_str, eta_str), end='\r')


def brute_compute(brute_array, func_args, 
                  brute_iterations, start_time, estimated_total_time,
                  slice_number = 0):
    '''
    Brute force optimisation, computing CHI^2 at each value in 'brute_array'.
    Warning!: Uses recursion!
        This allows the function to be used for brute force calculations
        that include a 3rd variable - e.g. background scaling, or r_min.
        But watch out for python's recursion limit!
    '''
    if np.shape(brute_array)[0] == 2:
        outer = brute_array[0]
        inner = brute_array[1]
        res = []
        for slice_number, n_iter in enumerate(outer):
            func_args['n_iter'] = n_iter
            brute_slice = brute_compute(inner, func_args, brute_iterations,
                                        start_time, estimated_total_time,
                                        slice_number=slice_number)
            res.append(brute_slice)
        return np.asarray(res)
    else:
        pass
    func_args = tuple(func_args.values())
    res = []
    for iteration, rho in enumerate(brute_array):
        chi_sq = liquid.refinement_objfun(rho, *func_args)
        elapsed_time = time.time() - start_time
        iteration_num = (iteration + 1) + (len(brute_array)*slice_number)
        update_progress_bar(iteration_num, brute_iterations, elapsed_time, estimated_total_time)
        res.append(chi_sq)
    res = np.asarray(res)
    return res

# Data Options
data_file = 'example_data.dat'
q_max = 11.8
comp = {'Ga':(31,0,1)}
dq = 0.02
sq_form = 'faber-ziman'
mod_func = None
window_start = None
r_min = 2.3
fft_N = 12

# Set the array of values to compute X^2
brute_rho_array = np.arange(0.01, 0.08, 0.001)
brute_eggert_iterations_array = np.arange(1, 50, 1)

# Load data
x, y = np.loadtxt(data_file, unpack=True)
# rebin with set dq (steps of q values)
x, y = rebin_data(x, y, dx=dq)
# Cut data at Q_max
x, y = x[np.where(x < q_max)], y[np.where(x < q_max)]

func_args = {'Q_data': x, 'I_data': y, 'comp': comp,
             'r_min': r_min, 'n_iter': int(np.mean(brute_eggert_iterations_array)),
             'method': sq_form, 'mod_func': mod_func,
             'window_start': window_start, 'fft_N': fft_N,
             'opt_rho': 1, 'opt_bkg': 0}

brute_array = (brute_eggert_iterations_array, brute_rho_array)
brute_iterations = len(brute_rho_array) * len(brute_eggert_iterations_array)
# Here we calculate the time to make 1 iteration
t_start = time.time()
brute_compute([1], func_args, 1, 1, 1, 1)
t_end = time.time()
iteration_time = t_end - t_start
# Estimate total time to complete calculations in seconds
estimated_total_time = (iteration_time * brute_iterations)
start_time = time.time()
# Calculate Chi^2 at all values of rho
brute_res = brute_compute(brute_array, func_args,
                          brute_iterations, start_time, estimated_total_time)

brute_res = np.nan_to_num(brute_res)
np.save('brute_res.npy', brute_res, allow_pickle=True)
np.save('brute_rho.npy', brute_rho_array, allow_pickle=True)
np.save('brute_niter.npy', brute_eggert_iterations_array, allow_pickle=True)

fig1 = plt.figure()
for i in range(len(brute_res)):
    plt.plot(brute_rho_array, brute_res[i], label=('n_iter = ' + str(brute_eggert_iterations_array[i])))
plt.legend()
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\chi^2$')
plt.show()
#plt.savefig('fitting region spread.svg')

map_x, map_y = np.meshgrid(brute_rho_array, brute_eggert_iterations_array)
brute_res_normalised = []
for niter_slice in brute_res:
     brute_res_normalised.append(niter_slice/np.min(niter_slice))
brute_res_normalised = np.asarray(brute_res_normalised)

log_map = np.log(brute_res_normalised)
mask = np.where(log_map == 0)
dbl_log_map = log_map
dbl_log_map[mask] = 1e-12
dbl_log_map = np.log(dbl_log_map)

fig2 = plt.figure()
opt_density = []
for index, niter in enumerate(brute_eggert_iterations_array):
    opt_density.append(brute_rho_array[np.where(dbl_log_map[index] == np.min(dbl_log_map[index]))][0])
opt_density = np.asarray(opt_density)
plt.plot(brute_eggert_iterations_array, 1/opt_density)
plt.xlabel('Number of iterations in refinement procedure')
plt.ylabel(r'$\rho^{-1}$')
plt.show()
#plt.savefig('refined density vs n iterations.png', dpi=1200)


fig3 = plt.figure()
# Plot the minimum points
cp = plt.contourf(map_x, map_y, dbl_log_map, 50)
plt.title(r'Minimum $\chi^2$ values')
plt.ylabel('No. iterations in refinement procedure')
plt.xlabel(r'Refined $\rho_0$ $(atoms/\AA^3)$')
plt.show()
#plt.savefig('rho_vs_iterations', dpi=2400)

# Calculate curvature at each Chi^2 minimum
second_derivatives = []
for res in brute_res:
    rho_min = brute_rho_array[np.where(res == np.min(res))]
    spline = UnivariateSpline(brute_rho_array, res/1e6)
    second_derivatives.append(spline(rho_min))
second_derivatives = np.asarray(second_derivatives)

# Plot d^2chi^2 / drho^2 vs n_iter
fig4, ax1 = plt.subplots()
ax1.plot(brute_eggert_iterations_array, second_derivatives)
plt.title('Second derivative at $\chi^2$ minima vs No. of iterations')
plt.xlabel('Number of iterations in refinement procedure')
plt.ylabel(r'$\frac{d^2\chi^2}{d\rho^2}$ $at$ $\rho{}={}\rho_{min}$')
ax2 = fig4.add_axes([0.35, 0.45, 0.3, 0.2])
ax2.plot(brute_eggert_iterations_array[3:17], second_derivatives[3:17])
plt.xticks(brute_eggert_iterations_array[3:17][::2])
plt.show()
#plt.savefig('second_derivatives', dpi=2400)