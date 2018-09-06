#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python implementation of Eggert method for estimating liquid structure factor
by optimising rho/alpha.
"""
__author__ = "Benedict J Heinen"
__copyright__ = "Copyright 2018, Benedict J Heinen"
__email__ = "benedict.heinen@gmail.com"

import os
import numpy as np
#from scipy.optimize import minimize
from scipy.integrate import simps
import scipy.interpolate
import scipy.fftpack

def calc_Z_sum(composition):
    '''
    Calculates the sum of the atomic number of all elements in the composition
    Args: 
        composition - dictionary of atomic elements with values as tuples of
        the form (Z,charge,n)
    Returns:
        z_tot - total atomic number
    '''
    Z_tot = 0
    for element in composition:
        Z, _, n = composition[element] 
        Z_tot += Z * n
    return Z_tot


def calc_atomic_ff(element,Q):
    '''
    Calculates the atomic form factors fp(Q) using the analytic tabulation
    of ...???....
    
    Args:
        element - tuple of form (Z,charge)
        Q - Data in Q space and units Angstroms^-1
    
    Returns: 
        form_factor - NumPy array of atomic form factors at each Q value
    '''
    atomic_number, charge, *_ = element
    __data_path = os.path.join(os.path.abspath(os.getcwd()), 'data')
    ff_data = np.loadtxt(os.path.join(__data_path, 'Mcad-formfactors.dat'))
    ff_data = ff_data[np.where((ff_data[:,0]==charge) & (ff_data[:,10]==atomic_number))][0,1:10]        
    s = (Q / (4 * np.pi))**2
    form_factor = (ff_data[0] * np.exp(-ff_data[1] * s) + 
                   ff_data[2] * np.exp(-ff_data[3] * s) + 
                   ff_data[4] * np.exp(-ff_data[5] * s) + 
                   ff_data[6] * np.exp(-ff_data[7] * s) + ff_data[8]
                   )
    return form_factor


def calc_effective_ff(composition,Q):
    '''
    Calculates the effective electronic form factor as defined in Eq. 10 of
    Eggert et al., 2002.
    
    i.e.
    
    f_e(Q) = SUM f_p(Q) / Z_tot
    where f_p(Q) are the atomic form factors and Z_tot is the total atomic 
    number of the compositional unit
    
    Args:
        composition - composition dictionary where values are tuples of the 
        form (Z, charge, n)
        Q - q space values to evaluate over
    
    Returns:
        effective_ff - effective electronic form factor at Q values
        atomic_ff - atomic form factors at Q values for each element
    '''    
    Z_tot = calc_Z_sum(composition)
    atomic_ff = []       
    for element in composition:
        atomic_ff.append(calc_atomic_ff(composition[element],Q))
    atomic_ff = np.asarray(atomic_ff)
    effective_ff = np.sum(atomic_ff,0) / Z_tot
    return effective_ff, atomic_ff
        

def calc_compton_scattering(element, Q):
    '''
    Interpolates compton scattering from tabulated data to estimate scattering 
    intensity at each Q value
    
    Args:
        element - dictionary entry from 'comp' where keys are element symbols
        Q - Q values to estimate compton scattering intensity at
    '''    
    __data_path = os.path.join(os.path.abspath(os.getcwd()), 'data')    
    filename = os.path.join(__data_path, 'hubbel-compton', (element + '.cmp'))
    cs_Q, _, cs_comp = np.loadtxt(filename, unpack=True, skiprows=1)
    cs_Q *= (np.pi * 4)
    # Interpolate compton scattering to required Q spacings
    finterp = scipy.interpolate.interp1d(cs_Q, cs_comp, kind='cubic', fill_value='extrapolate')
    return finterp(Q)


def calc_J(composition, Q):
    '''
    Calculates J value as defined in Eq 35 of Eggert et al., 2002
    i.e.
        J(Q) = SUM {Compton(Q)} / Z_tot^2 * f_e(Q)^2
        where f_e(Q) is the effective electronic form factor
    
    Args:
        composition - dictionary where keys are element symbols and values
        are tuples of the form (Z, charge, n)
        
        Q - Q values to evaluate
    '''

    effective_ff, _ = calc_effective_ff(composition,Q)
    Z_tot = calc_Z_sum(composition)
    compton_intensity_total=0
    for element in composition:
        compton_intensity_total += calc_compton_scattering(element,Q)
    J = compton_intensity_total / (Z_tot**2 * effective_ff**2)
    return J

 
def calc_K_p(composition, Q):
    '''
    Calculates effective average effective atomic numbers of each atomic
    element in the composition as defined in Eqs 14 & 11 of Eggert et al., 2002
    
    i.e.
    
    K_p(Q) for each element = f_p(Q)/f_e(Q)
    K_p for each element = <K_p(Q)>_Q
    
    Args:
        composition - dictionary of elements in composition where values
        are tuples of the form (Z, charge, n)
        
        Q - Q range to evaluate over
    '''
    effective_ff, atomic_ff = calc_effective_ff(composition,Q)
    K_p = atomic_ff / effective_ff
    K_p = np.mean(K_p,1)
    return K_p  
    

def calc_S_inf(composition, Q):
    '''
    Caclulates the value S_inf as defined in Eq 19 of Eggert et al., 2002
    i.e. S_inf = SUM K^2_p / Z_tot^2
    S_inf defaults to 1 if the composition is monatomic
    '''
    if len(composition) == 1:
        S_inf = 1.0
    else:
        K_p = calc_K_p(composition,Q)
        Z_tot = calc_Z_sum(composition)
        S_inf = np.sum(K_p**2) / Z_tot**2
    return S_inf


def calc_alpha(Q_cor, I_cor, Z_tot, rho, J, S_inf, effective_ff):
    '''
    Calculates 'alpha' - the normalisation factor defined in Eq 34 of 
    Eggert et al., 2002
    
    Args:
        Q_cor - Rebinned Q space data cut at Q_max
        I_cor - Background corrected intensity data
        Z_tot - Atomic number sum of elements in composition
        rho - actual or initial estimate of atomic density in atoms/Angstrom^3
        J - J value as defined in Eq 35 of Eggert et al., 2002
        S_inf - S_inf value as defined in Eq 19 of Eggert et al., 2002
        effective_ff - effective electronic form factor
        
    Returns:
        alpha
    '''
    # define the two integrals in the equation
    int_1 = simps((J + S_inf) * Q_cor**2, Q_cor)
    int_2 = simps((I_cor / effective_ff**2) * Q_cor**2, Q_cor)
    # calculate alpha
    alpha = Z_tot**2 * (((-2 * np.pi**2 * rho) + int_1) / int_2)
    return alpha


def calc_coherent_scattering(Q_cor, I_cor, composition, alpha):
	'''
	Calculate coherent scattering intensity as per Eq 27 of Eggert et al., 2002
	
	Coherent_scattering = N * [alpha * Sample_scattering(Q) - SUM Incoherent_scattering]

	Args:
		Q_cor
		I_cor
		composition
		alpha
	'''
	# First calculate incoherent (compton) scattering	
	incoherent_scattering = []
	for element in composition:
		incoherent_scattering.append(calc_compton_scattering(element,Q_cor))
	incoherent_scattering = np.sum(np.asarray(incoherent_scattering),0)

	# Next calculate coherent scattering
	N = len(composition)
	coherent_scattering = N * ((alpha * I_cor) - incoherent_scattering)
	return coherent_scattering


def calc_structure_factor(Q_cor, I_cor, composition, rho):
    '''
	Calculate the molecular structure factor using Eq 18 of Eggert et al., 2002
    
    Args:
        Q_cor - Rebinned Q space data
        I_cor - Bkgd corrected intensity data
        composition - dictionary of atomic elements, with values as tuples in 
                      the form (Z,charge,n)
        rho - actual or initial estimate of atomic density in atoms/Angstrom^3
	'''
    N = len(composition)
    Z_tot = calc_Z_sum(composition)
    J = calc_J(composition, Q_cor)
    S_inf = calc_S_inf(composition,Q_cor)
    effective_ff, _ = calc_effective_ff(composition,Q_cor)
    alpha = calc_alpha(Q_cor, I_cor, Z_tot, rho, J, S_inf, effective_ff)
    coherent_scattering = calc_coherent_scattering(Q_cor, I_cor, composition, alpha)
    
    structure_factor = coherent_scattering / (N * Z_tot**2 * effective_ff**2)
    return structure_factor
    

def calc_F_r(x, y, rho, dx='check', N=12, lorch=False, function='density_func'):
    '''
    Calculates the fourier transform of the interference function 
    i(Q) = S(Q) - S_inf to obtain the pair-distribution function [g(r)], the
    density function [F(r)] or the radial distribution function [RDF(r)]. This 
    conversion to real-space is useful for PDF-analysis.
    
    Example:
    
    F(q) --> G(r) transformation
    FFT here is sine-only and so assumes F(q) is an odd function
    F(q) is therefore padded with its odd image and zero-padded in the middle 
    for efficiency and for interpolation in r-space (padded len =  2^n).

    i.e. F(q) = [0 1 3 2] --> F(q)_pad = [0 1 3 2 0 -2 -3 -1]

    Values of G(r) are the imaginary components of the fourier transorm of
    F(q). The real components are all zero. Only the first half of the 
    transformed array is used as the second is its odd image.
    
    dr correponds to the minimum frequency available for the entire Q-range.
    i.e.:  
        dr = 2 * np.pi / (len(F(q)_padded) * dq)

    Therefore for the reverse transformation
        dq = np.pi / (len(r) * dr)

    The pair-distribution function is proportional to the probability of 
    finding an atom at a distance r from an average atom taken as the origin.
    The radial distribution function gives atomic coordination numbers when 
    integratd across peaks. The density function is used in the analysis here
    to normalise i(q).

    Pair-distribution function g(r):
        
        g(r) - 1 = 1/(2 * pi^2 * r * rho_0) * INT[ q[i(q)] sin(qr) dq]
    
    Density function, D(r) = F(r):
        
        F(r) = 4 * pi * r * rho_0 * [g(r) - 1]
        F(r) = 2/pi * INT[ q[i(q)] sin(qr) dq]
        
    Radial distribution function, RDF(r):
        
        RDF(r) = 4 * pi * r^2 * rho_0 * g(r)
        RDF(r) = 2r*INT[ q[i(q)] sin(qr) dq] - 4*pi*r^2*rho_0
    

    Args:
        x   - Q values (numpy array)
        y   - Corresponding values of the interference function (numpy array)
        rho - atomic density (float)
        dq  - steps of q values (which need to be constant)
              this defaults to 0.02, but can also be calculated in function if 
              set to 'check'
        N   - Sets size of array for fft where array_size = 2^N, defaults to 12 
        lorch - Flag to use Lorch function (bool)
        function - density_func, pair_dist_func, radial_dist_func
        
    Returns:
        r   - Array of r values
        F_r/g_r/RDF_r - Corresponding array of values for calculated function 

    '''    
    # Define dx (qstep) if not passed to function
    if dx == 'check':
        dx = x[1] - x[0]
    else:
        pass        
    # Calculate qi(q) input array
    z = x * y
    # Make sure using float64
    z = z.astype(np.float64)
    # Pad with enough zeros to make array of size 2**N
    padding_zeros = np.int(2**N - (len(z)*2 - 1))    
    # Pad array with odd image of function
    z = np.concatenate((z, np.zeros(padding_zeros), -np.flip(z[1::],0)))
    # Fourier transform here uses scipy.fftpack over numpy.fft
    # ifft used as converting from frequency (intensity) domain
    # This means no scaling is needed and the result is not negative
    fft_z = scipy.fftpack.ifft(z)
    # Only the imaginary component of the transformed array is used
    # Only the first half of the array is used as the rest is the odd component
    fft_z = np.imag(fft_z)[:len(fft_z)//2]
    # Calculate steps of r and r array
    dr = np.pi*2/(len(z)*dx)
    r = np.arange(len(fft_z))*dr
    if function == 'density_func':
        # Scale by factor 2/pi 
        F_r = fft_z * 2/np.pi
        return r, F_r
    elif function == 'pair_dist_func':
        sf = 1 / (2 * r * rho * np.pi**2)
        g_r = fft_z * sf
        g_r += 1
        return r, g_r
    elif function == 'radial_dist_func':
        RDF_r = (2 * r * fft_z / np.pi) + (4 * np.pi * rho * r**2)
        return r, RDF_r        
    else:
        raise ValueError('arg \'function\' must be valid option')


def calc_F_r_iteration_term(delta_F_r, N=12):
    '''
    Calculates the term responsible for oscillations at small r values.
    Equation 45 in Eggert et al., 2002
    
    i.e. Delta_alpha * Q * S_inf = INT 0>r_min [Delta_F(r) sin(Qr) dr]
    
    Delta_F(r) is the difference between F(r) and its expected behaviour
    
    
    
    '''
    padding_zeros = np.int(2**N - (len(delta_F_r)*2 - 1))
    z = np.concatenate((delta_F_r, 
                        np.zeros(padding_zeros), 
                        -np.flip(delta_F_r[1::],0)))
    fft_z = -np.imag(scipy.fftpack.fft(z))
    
    return fft_z


def calc_model_F_intra_r(Q, r, composition, rho, d_pq):
    '''
    d_pq is list of distances?
    not sure how to do this for polyatmoic cases
    > stick to monatomic for now
        
    divide by 4pi^2 ??? why!?
    
    '''
    # Monatomic case:
    if len(composition) == 1:
        d = d_pq
        r_minus_d = r - d
        r_plus_d = r + d
        K_p = calc_K_p(composition,Q)[0]
        Z_tot = calc_Z_sum(composition)
        Q_max = np.max(Q)
        F_intra_r = K_p**2/(np.pi*d*Z_tot**2) *\
                    ((np.sin(r_minus_d*Q_max)/r_minus_d) - 
                     (np.sin(r_plus_d*Q_max)/r_plus_d))
        model_F_intra_r = F_intra_r - (4*np.pi*r*rho)
        model_F_intra_r /= (4*np.pi**2)
        return model_F_intra_r
    
    else:
        print('Cannot calculate F_intra(r) for polaytomic species')


def calc_chi_squared(r_intra, delta_F_r, method='simps'):
    '''
    Calculates a chi squared figure of merit as in equation 50 of Eggert. This 
    is a measure of how small Delta_F_i(r) is, to be used as a tolerance factor
    to check convergence of i_i+1(q) and in the optimisation routine of rho/s.
    
    Chi^2 = INT 0-r_min [Delta_F(r)]^2 dr
    
    The value chi^2 is defined as the area under the curve Delta_F(r) limited
    to r_min (squared to remove negatives).
    
    Args:
        delta_F_r - numpy array of delta_F_r from r=0 to r=r_min
    
    Returns:
        chi_squared - Value of chi^2
    '''
    # Square to remove negatives
    f = delta_F_r ** 2
    if method == 'simps':
        # Integrate using simpson's rule
        chi_squared = simps(f, x=r_intra)
    return chi_squared


def stop_iteration(chi_squared, count, iter_limit=5):
    if count < iter_limit:
        #print('count', count)
        #print(chi_squared)
        return False
    else:
        #print('stop the loop!')
        return True


def calc_impr_interference_func(rho, *args):
    '''
    Calculates an improved estimate of the interference function via the
    iterative procedure described in Eggert et al., 2002. This is done by
    scaling the data and forcing its behaviour in the low r region to some
    modelled behaviour.
    
    The form of this function also enables it to be passed to an optimisation
    routine to estimate rho by minimising chi^2.
    
    Args:
        rho - density
        
        args - tuple of static parameters and arguments
            If opt_flag == False, expected form is:
            (q_dat, interference_func, composition, r_min, d_pq, opt_flag)
            If opt_flag == True, expected form is:
            (q_dat, I_dat, composition, r_min, d_pq, iter_limit, opt_flag)
        
        q_dat - numpy array of Q space data
        I_dat - Background corrected/scaled intensity data for the sample
        interference_func - initial interference function (if using fixed rho)
        composition - dict of elements, values as tuples of form (Z,charge,n)
        r_min - intramolecular distance cut-off
        d_pq - distance between neighbouring species used to calculate expected
               intramolecular behaviour at low r.
        opt_flag - function returns chi_squared to minimisation routine if true
                   function returns improved interference_func with 
                   chi squared value if false
    '''
    # Unpack static parameters & arguments
    *_, opt_flag = args
    if opt_flag == False:
        (q_dat, interference_func, composition, r_min, d_pq, iter_limit, opt_flag) = args
    elif opt_flag == True:
        (q_dat, I_dat, composition, r_min, d_pq, iter_limit, opt_flag) = args
        # Calculate initial interference func for rho value
        structure_factor = calc_structure_factor(q_dat,I_dat, composition, rho)
        interference_func = structure_factor - calc_S_inf(composition, q_dat)
    else:
        raise ValueError('Arg - opt_flag - must be boolean')

    # Calculate initial F(r)
    r, F_r = calc_F_r(q_dat, interference_func, rho)
    # Calculate expected behaviour of F(r) for intramolecular distances (r<r_min)
    model_F_intra_r = calc_model_F_intra_r(q_dat, r, composition, rho, d_pq)
    
    # Calculate static terms of iterative proceduce
    with np.errstate(divide='ignore', invalid='ignore'):
        t1 = 1 / q_dat
    t2_divisor = calc_S_inf(composition,q_dat) + calc_J(composition,q_dat)

    # Count no. iterations for verbosity or control
    # Function set up to call stop_looping func in case want to modify
    # to use chi-squared value as tolerance value instead
    count = 1
    done_looping = False
    while not done_looping:
        try:
            done_looping = stop_iteration(chi_squared, count, iter_limit=iter_limit)
            interference_func = interference_func_impr
            r, F_r = calc_F_r(q_dat, interference_func, rho)
        except NameError:
            pass
                        
        delta_F_r = (F_r - model_F_intra_r)[np.where(r<r_min)]
        r_intra = r[np.where(r<r_min)]
        chi_squared = calc_chi_squared(r_intra, delta_F_r)
        
        t2 = (interference_func/t2_divisor) + 1
        t3 = calc_F_r_iteration_term(delta_F_r)[:len(interference_func)]
        with np.errstate(invalid='ignore'):
            interference_func_impr = interference_func - (t1 * t2 * t3)
        count+=1

    if opt_flag == True:
        return chi_squared
    elif opt_flag == False:
        return interference_func_impr, chi_squared
    else:
        raise ValueError('Argument - opt_flag - must be boolean')
        


'''
Example use without GUI


# Composition with element name as key, entry as tuple with Z, charge, fraction
composition = {'Ga': (31,0,1)}
# Set density
rho = 0.05

# Load your data
q_raw, I_raw = np.loadtxt('Ga_backcor_origin.dat', unpack=True, skiprows=0)
# Re-bin via interpolation so dq(steps of q) = 0.02
# Cutoff below Q_cutoff
q_cutoff = 9.0
dq = 0.02
q_dat = np.arange(0, q_raw[-1], dq)
q_dat = q_dat[q_dat<q_cutoff]
finterp = scipy.interpolate.interp1d(q_raw, I_raw, kind='cubic', fill_value='extrapolate')
I_dat = finterp(q_dat)

# Smooth the data > using savitsky-golay filter
if you wish

#############################################################

# First calculate the interference function i(Q)
# i(Q) = S(Q) - S_inf
structure_factor = calc_structure_factor(q_dat,I_dat, composition, rho)
interference_func = structure_factor - calc_S_inf(composition, q_dat) 
# store original interference function
interference_func_0 = interference_func

r_min = 2.3
d_pq = 2.9
rho_0 = 0.05
args = (q_dat, interference_func_0, composition, r_min, d_pq, 0)
interference_func_1 = calc_impr_interference_func(rho_0, *args)density_input


args = (q_dat, I_dat, composition, r_min, d_pq, 1)
bounds = ((0.02, 0.08),)
op_method = 'L-BFGS-B'
optimisation_options = {'disp': 1,
                        'maxiter': 15000,
                        'maxfun': 15000,
                        'ftol': 2.22e-8,
                        'gtol': 1e-10
                        }
opt_result = minimize(calc_impr_interference_func, rho_0,
                      bounds=bounds, args=args,
                      options=optimisation_options,
                      method=op_method)

rho_impr = opt_result.x[0]



sf = calc_structure_factor(q_dat,I_dat, composition, rho_impr)
int_f_2 = sf - calc_S_inf(composition, q_dat) 
args = (q_dat, int_f_2, composition, r_min, d_pq, 0)
interference_func_2 = calc_impr_interference_func(rho_impr, *args)

# inital guess
plt.plot(q_dat, interference_func_0, color='g')
# optimised for rho = 0.05
plt.plot(q_dat, interference_func_1, color='r')
# for optimised rho
plt.plot(q_dat, interference_func_2, color='b')

'''
