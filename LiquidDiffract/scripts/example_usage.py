#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script shows an example usage of the LiquidDiffract library. Whilst LiquidDiffract
was primarily designed as a graphical application, the functions in its core module can
be easily called from a python script. This is useful for batch processing files (where
an ongoing graphical view of the steps is not neccessary), or computing probability 
density functions from pre-existing structure factors (fourier transformation). It can 
also be used to call functions that are intermediary calculation steps in the graphical
application - e.g. retreiving atomic form factors or calculating compton scattering 
intensities for an element or multi-component composition.

The two primary functions in normal usage are:

    LiquidDiffract.core.core.calc_structure_factor
    LiquidDiffract.core.core.calc_impr_interference_func

calc_structure_factor computes the structure factor, S(Q), for given sample
composition, density, and experimentally measured sample scattering data.

calc_impr_interference_func optimises the interference function, S(Q) - S_inf,
based on the iterative method developed by Eggert et al., 2002, Phy Rev B, 65(17).

See <https://github.com/bjheinen/LiquidDiffract> for more details.
"""
__author__ = 'Benedict J. Heinen'
__copyright__ = 'Copyright 2018-2022, Benedict J. Heinen'
__email__ = 'benedict.heinen@gmail.com'
__license__='GPLv3'

# Import required python modules
import numpy as np
# Import optional python modules used in this example
import matplotlib.pyplot as plt

from scipy.optimize import minimize
# import the LiquidDiffract core module
import LiquidDiffract.core.core as liquid
import LiquidDiffract.core.data_utils as data_utils

# First set composition
# The composition should be a dictionary dictionary with short element names as keys,
# and values as tuples of the form (Z, charge, number_of_atoms)
composition = {'Ga': (31,0,1)}
# Set initial density - in atoms per cubic Angstrom
rho = 0.055

# Load data
# LiquidDiffract expects x data to be Q values in Angstroms (not nm!)
# The file 'example_data.dat' has already been background corrected
q_raw, I_raw = np.loadtxt('example_data.dat', unpack=True, skiprows=0)

# Next rebin data and trim if necessary
# Re-bin via interpolation so dq (steps of q) is consistent to allow FT
# This can be done using the rebin_data function from core.data_utils
dq=0.02
q_data, I_data = data_utils.rebin_data(q_raw, I_raw, dx=dq, x_lim=[q_raw[0], 11.8])

# Apply any other data treatment
# e.g. a savitsky-golay filter to smooth the data
I_data = data_utils.smooth_data(I_data, method='savitzky-golay', 
                                window_length=31, poly_order=3)

# First calculate the interference function i(Q)
# i(Q) = S(Q) - S_inf (S_inf = 1 for monatomic samples, or in the Faber-Ziman formalism)
structure_factor = liquid.calc_structure_factor(q_data, I_data, composition, rho, method='faber-ziman')
interference_func = structure_factor - liquid.calc_S_inf(composition, q_data, method='faber-ziman')

# store original interference function
interference_func_0 = interference_func

# A low r region cutoff must be chosen for the refinement method of Eggert et al. (2002)
# A density value in atoms(or molecules) / A^3 must also be set. This will be
# used as the starting estimate if the density is refined.
r_min = 2.3
rho_0 = rho

# The function 'calc_impr_interference_func' calculates an improved estimate of
# the interference function via the iterative procedure described in Eggert et al., 2002.
#
# calc_impr_interference_func requires the following arguments:
# q_data - q values
# I_data / i(Q) - the treated intensity data or interference function (depends on opt_flag used)
# composition - composition dictionary
# r_min - low r cut-off
# iter_limit - iteration limit for Eggert procedure
# method - method of calculating S(Q) ('ashcroft-langreth' and 'faber-ziman' are currently supported)
# mod_func - modification function to use ('None', 'Cosine-window' or 'Lorch')
# window_start - window_start of cosine function (if in use)
# fft_N sets the size of the array during the fourier transform to 2**N
# opt_flag - set to 0 if not running from solver
#
# The last positional argument (opt_flag) signals that you only want the chi_squared value to be returned.
# This is useful if using a solver to refine the density by minimising the chi_squared value

# e.g.
# Return the improved i(Q) & chi^2 at density = rho_0
iter_limit = 25
method = 'faber-ziman'
mod_func = 'Cosine-window'
window_start = 7
fft_N = 12
args = (q_data, interference_func_0, composition, rho_0,
        r_min, iter_limit, method, mod_func, window_start, fft_N)
# Store the refined interference function at rho_0
interference_func_1, chi_sq_1 = liquid.calc_impr_interference_func(*args)

# Next we pass the function calc_impr_interference_func to
# a solver to estimate the density. We use a separate objective function
# core.redinement_objfun to make this simpler.
# The argument opt_rho is set to 1 to indicate we are refining the density
opt_rho = 1
# The argument opt_bkg is set to 0 to indicate we are not refining
# the background sclaing factor.
opt_bkg = 0
# When refinind the density an iter_limit <= 10 is recommended
iter_limit_refine = 7
args = (q_data, I_data, composition, r_min,
        iter_limit_refine, method,
        mod_func, window_start, fft_N,
        opt_rho, opt_bkg)
# Set-up bounds and other options according to the documentation of solver/minimisation routine
bounds = ((0.045, 0.065),)
op_method = 'L-BFGS-B'
optimisation_options = {'disp': 0,
                        'maxiter': 15000,
                        'maxfun': 15000,
                        'ftol': 2.22e-12,
                        'gtol': 1e-12
                        }
opt_result = minimize(liquid.refinement_objfun, rho_0,
                      bounds=bounds, args=args,
                      options=optimisation_options,
                      method=op_method)
# The solver finds the value of rho that gives the smallest chi^2
rho_refined = opt_result.x[0]

# The interference function can then be re-calculated using the new value for rho,
# and the F(r) refined as above.
interference_func = (liquid.calc_structure_factor(q_data, I_data, composition, rho_refined) -
                     liquid.calc_S_inf(composition, q_data))
args = (q_data, interference_func_0, composition, rho_refined,
        r_min, iter_limit, method, mod_func, window_start, fft_N)
interference_func_2, chi_sq_2 = liquid.calc_impr_interference_func(*args)

# Calculate the corresponding pair distribution functions g(r)
r, g_r_0 = liquid.calc_correlation_func(q_data, interference_func_0, rho_0, dx=dq,
                                        mod_func=mod_func, window_start=window_start,
                                        function='pair_dist_func')
# r values will be the same
_, g_r_1 = liquid.calc_correlation_func(q_data, interference_func_1, rho_0, dx=dq,
                                        mod_func=mod_func, window_start=window_start,
                                        function='pair_dist_func')
_, g_r_2 = liquid.calc_correlation_func(q_data, interference_func_2, rho_refined, dx=dq,
                                        mod_func=mod_func, window_start=window_start,
                                        function='pair_dist_func')

# Plot the data
fig1 = plt.figure()
plt.xlabel(r'Q ($\AA^{-1}$)')
plt.ylabel(r'i(Q)')
# Initial interference function calculation
plt.plot(q_data, interference_func_0, color='g', label=r'Initial i(Q) | $ρ = {rho:.2f}$'.format(rho=rho_0))
# Optimised at rho_0
plt.plot(q_data, interference_func_1, color='r', label=r'Optimised i(Q) | $ρ = {rho:.2f}$ & $χ^{{2}} = {chisq:.2f}$'.format(rho=rho_0, chisq=chi_sq_1))
# Optimised, using refined density estimate
plt.plot(q_data, interference_func_2, color='b', label=r'Optimised i(Q) | $ρ = {rho:.3f}$ & $χ^{{2}} = {chisq:.2f}$'.format(rho=rho_refined, chisq=chi_sq_2))
# Add a legend
plt.legend(loc='best')
# Show the plots
plt.show()

fig2 = plt.figure()
plt.xlabel(r'r ($\AA$)')
plt.ylabel(r'g(r)')
plt.ylim((np.nanmin(g_r_2), np.nanmax(g_r_2)))
window = len(q_data)
# Initial g(r)
plt.plot(r[:window], g_r_0[:window], color='g', label=r'Initial g(r) | $ρ = {rho:.2f}$'.format(rho=rho_0))
# Optimised g(r) at rho_0
plt.plot(r[:window], g_r_1[:window], color='r', label=r'Optimised g(r) | $ρ = {rho:.2f}$ & $χ^{{2}} = {chisq:.2f}$'.format(rho=rho_0, chisq=chi_sq_1))
# Optimised g(r), using refined density estimate
plt.plot(r[:window], g_r_2[:window], color='b', label=r'Optimised g(r) | $ρ = {rho:.3f}$ & $χ^{{2}} = {chisq:.2f}$'.format(rho=rho_refined, chisq=chi_sq_2))
# Add a legend
plt.legend(loc='best')
# Show the plots
plt.show()
