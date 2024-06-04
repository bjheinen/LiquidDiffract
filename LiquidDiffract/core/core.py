#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LiquidDiffract core functions module

Python implementation of iterative method for liquid structure factor
normalisation based on Eggert et al., 2002.
Optimises the normalisation factor, alpha, and allows refinement of sample density, rho.

Provides functions for:

  * Calculation of liquid structure factor, S(Q)
  * Fourier transformation of S(Q)/g(r)
  * Calculation of pair and radial distribution functions, g(r) and RDF(r)
  * Calculation of atomic form factors, effective electronic form factor, and scattering intensities
  * Utility functions for density and molecular mass conversions/automatic calculation

See <https://github.com/bjheinen/LiquidDiffract#using-liquiddiffract-core-library> for more details

"""
__author__ = 'Benedict J. Heinen'
__copyright__ = 'Copyright 2018-2024, Benedict J. Heinen'
__email__ = 'benedict.heinen@gmail.com'
__license__ = 'GNU GPL v3'

from functools import lru_cache
import itertools
import numpy as np
from scipy.integrate import simpson, quad
import scipy.interpolate
import scipy.fftpack
try:
    import importlib_resources as resources
except ImportError:
    from importlib import resources

# Get version number from version.py
from LiquidDiffract.version import __appname__, __version__
from LiquidDiffract.core.data_utils import data_cache, smooth_data, convert_q_space


def calc_mol_mass(composition):
    '''
    Calculates the molecular mass for a given composition and returns
    the average molecular mass for one atom

    composition is dictionary in the form: (Z, charge, n)
    where n is number of atoms in formula unit
    '''
    with resources.files('LiquidDiffract.resources').joinpath('mass_data.npy').open('rb') as fp:
        mass_dict = np.load(fp, allow_pickle=True).item()
    mol_mass = np.sum([mass_dict[element]*(composition[element][2]) for element in composition])
    return mol_mass


def conv_density(rho, composition):
    '''
    Converts atomic number density - atoms/angstrom^3
    to mass density in g/cm^3
    '''
    # Calculate molecular mass g/mole
    mol_mass = calc_mol_mass(composition)
    # Calculate number of atoms and scale mol_mass
    N = sum([element[2] for element in composition.values()])
    # Calculate atoms/cm3 > then moles/cm3
    with np.errstate(divide='ignore', invalid='ignore'):
        mass_density = (rho * 10 / 6.0221408) * (mol_mass/N)
    return mass_density


def calc_self_shielding(Q, mu, alpha, thickness, wavelength):
    '''
    Calculates the self-shielding attenuation factor for a sample in a simple slab geometry.

    Args:
        Q - Q-space value(s) to evaluate attenuation at
        mu - attenuation coefficient in m^-1
        alpha - angle between slab face and incident beam in degrees (90 is normal)
        thickness - sample/slab thickness in m
        wavelength - x-ray wavelength in Angstroms
    Returns:
        attenuation - Self-shielding attenuation factor A_s,s
    '''
    Q = Q.astype(np.longdouble)
    # Convert alpha to radians
    alpha = np.radians(alpha)
    # Get two theta values for Q-space data
    two_theta = convert_q_space(Q, wavelength)
    # Get direct path length (extended by slab rotation)
    L = thickness / np.cos(np.pi/2 - alpha)
    # Get scattered beam path length (scattering angle dependent)
    x = L * (np.sin(alpha) / np.sin(alpha + two_theta))
    # Get attenuation factor
    # Supress divide by zero warnings
    with np.errstate(divide='ignore', invalid='ignore'):
        _angle_factor = (np.exp(mu*(x-L))-1) / (mu*(x-L))
    # Handle for two_theta = 0 (angle_factor = 1.0)
    _angle_factor = np.nan_to_num(_angle_factor, nan=1.0)
    attenuation = np.exp(-1 * mu * x) * _angle_factor
    attenuation = attenuation.astype(np.float64)
    return attenuation


def calc_Z_sum(composition):
    '''
    Calculates the sum of the atomic numbers of all atoms in a
    compositional unit

    Args:
        composition - dictionary of atomic elements with values as tuples of
        the form (Z,charge,n)
    Returns:
        z_tot - total atomic number
    '''
    # Loop and tuple unpacking used over indices and list comprehension for
    # readability
    # Equivalent to:
    # sum(comp[el][0]*comp[el][2] for el in comp)
    # where comp = composition and el = element
    Z_tot = 0
    for element in composition:
        # Unpack Z and n
        Z, _, n = composition[element]
        Z_tot += (Z * n)
    return Z_tot


@lru_cache(maxsize=1)
def load_ff_data():
    with resources.files('LiquidDiffract.resources').joinpath('form_factors.npy').open('rb') as fp:
        ff_data = np.load(fp, allow_pickle=True)
    return ff_data


@data_cache(maxsize=32)
def calc_atomic_ff(element, Q):
    '''
    Calculates the atomic form factors fp(Q) of an element, p, over the
    Q values, Q. Uses the analytic tabulation published in the International
    Tables of Crystallography.

    This data is also freely available at:
    http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php

    Args:
        element - tuple of form (Z, charge, *_)
        Q - Data in Q space and units Angstroms^-1

    Returns:
        form_factor - NumPy array of atomic form factors at each Q value

    note: calc_atomic_ff ignores n values - atomic fractions or number of
          atoms in formula unit must be treated elsewhere as the function
          calculates form factor for 1 atom in the compositional unit.
    '''
    atomic_number, charge, *_ = element

    # Data is loaded in a separate function so it can be cached
    ff_data = np.copy(load_ff_data())
    ff_data = ff_data[np.where((ff_data[:, 0] == charge) &
                               (ff_data[:, 10] == atomic_number))][0, 1:10]
    # s = Q / 4pi, s here is s^2
    s = (Q / (4 * np.pi))**2
    form_factor = (ff_data[0] * np.exp(-ff_data[1] * s) +
                   ff_data[2] * np.exp(-ff_data[3] * s) +
                   ff_data[4] * np.exp(-ff_data[5] * s) +
                   ff_data[6] * np.exp(-ff_data[7] * s) + ff_data[8]
                   )
    return form_factor


def calc_effective_ff(composition, Q):
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
    # Reduce composition dict to list of values for clarity
    # Prefer explicit conversion to list of tuples over
    # e.g. list(composition.values()) to prevent errors if
    # dictionary values passed as lists instead of tuples
    composition = [tuple(element) for element in composition.values()]
    # Get counts of each species in composition (n values from (Z, charge, n))
    element_counts = np.array([element[2] for element in composition])
    # Get atomic form factors for each species in composition
    atomic_ff = np.array([calc_atomic_ff(atom, Q) for atom in composition])
    # Calculate effective electronic form factor
    # np.dot here is the sum product over the last axis of a and b.
    # np.dot used this way appears faster than e.g.:
    #   np.einsum('i,ij->j', element_counts, atomic_ff)
    effective_ff = np.dot(element_counts, atomic_ff) / Z_tot
    return effective_ff, atomic_ff


def calc_average_scattering(composition, Q):
    '''
    Computes the average scattering functions, <f>2 and <f2> required for the
    faber-ziman formalism of total structure factor S(Q). <f>2 is the square of
    the mean scattering factors, whereas <f2> is the mean square average.

    <f2> = 1/N Σ_p (f_p(Q))^2
    <f>2 = 1/N**2 Σ_p Σ_q f_p(Q) * f_q(Q)

    where:
        N is the number of atoms in the composition (e.g. MgSiO3: N = 5)
        f_p(Q) are the atomic atomic form factors of an element, p, calculated
        over the Q values, Q

    <f2> is equal to <f>2 for the monatomic case

    Note:
        The composition dict can be expanded to a list of tuples repeated for
        the number of atoms of that type with:
            composition_atoms = itertools.chain.from_iterable(
                itertools.repeat(tuple(composition[el]),composition[el][2])
                for el in composition)
        and the atomic ff for each individual atom then with e.g.:
            atomic_ff = []
            for atom in composition_atoms:
                atomic_ff.append(np.copy(calc_atomic_ff(atom, Q)))

        The average scattering functions are then equivalent to:
            1/N * np.einsum('ij,ij->j', atomic_ff, atomic_ff)
            and
            1/N**2 * np.einsum('ij,kj->j', atomic_ff, atomic_ff)
            (alternatively, np.sum(atomic_ff**2, 0) and np.sum(atomic_ff, 0)**2)

        Instead, here, the sum of the weighted form factors and weighted squared
        form factors is computed more efficiently using np.dot

    Args:
        composition - dictionary with keys that are element symbols and values
                      are tuples of the form (Z, charge, n_atoms)
        Q - Q-space values to evaluate over

    Returns:
        average_scattering - list of average scattering functions
            average_scattering[0] --> <f2>
            average_scattering[1] --> <f>2
    '''
    # Reduce composition dict to list of values for clarity
    # Prefer explicit conversion to list of tuples over
    # e.g. *list(composition.values())* to prevent errors if
    # dictionary values passed as lists instead of tuples
    composition = [tuple(element) for element in composition.values()]
    # Make np array of n_atoms for each element and calculate total n_atoms (N)
    element_counts = np.array([element[2] for element in composition])
    N = np.sum(element_counts)
    # Get atomic_ff for each type of atom
    atomic_ff = np.array([calc_atomic_ff(atom, Q) for atom in composition])
    # Create list to hold both average scattering functions
    average_scattering = []
    # 1st function - sum of squared form factors
    average_scattering.append(
            np.dot(element_counts, atomic_ff**2) / N
        )
    # 2nd function - square of summed form factors
    average_scattering.append(
            (np.dot(element_counts, atomic_ff) / N)**2
        )
    return average_scattering


@lru_cache(maxsize=8)
def load_compton_data(element):
    element_data_fname = element + '.npy'
    with resources.files('LiquidDiffract.resources.hubbel_compton').joinpath(element_data_fname).open('rb') as fp:
        cs_Q, _, cs_comp = np.load(fp, allow_pickle=True)
    return cs_Q, cs_comp


@data_cache(maxsize=32)
def calc_compton_scattering(element, Q):
    '''
    Interpolates compton scattering from tabulated data to estimate scattering
    intensity at each Q value.

    The analytical tabulations of compton scattering intensity are taken
    from Hubbel et al., 1975, J. Phys. Chem. Ref. Data 4, 471.

    Args:
        element - dictionary entry from 'comp' where keys are element symbols
        Q - Q values to estimate compton scattering intensity at
    '''
    # Element specific compton data is cached to save i/o calls
    cs_Q, cs_comp = load_compton_data(element)
    cs_Q, cs_comp = np.copy(cs_Q), np.copy(cs_comp)

    # The data is tabulated using s values instead of Q > Q = 4*pi*s
    cs_Q *= (np.pi * 4)
    # Interpolate compton scattering to required Q spacings
    finterp = scipy.interpolate.interp1d(cs_Q, cs_comp, kind='cubic', fill_value='extrapolate')
    return finterp(Q)


def calc_total_compton_scattering(composition, Q):
    '''
    Helper function to calculate total compton (incoherent) scattering for
    all atoms in a formula unit (composition)

    See function core.calc_compton_scattering for details on data used
    '''
    element_counts = np.array([element[2] for element in composition.values()])
    compton_scattering = np.array([calc_compton_scattering(element, Q) for element in composition.keys()])
    compton_scattering = np.dot(element_counts, compton_scattering)
    return compton_scattering


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
    effective_ff, _ = calc_effective_ff(composition, Q)
    Z_tot = calc_Z_sum(composition)
    compton_scattering = calc_total_compton_scattering(composition, Q)
    J = compton_scattering / (Z_tot**2 * effective_ff**2)
    return J


def calc_K_p(composition, Q):
    '''
    Calculates average effective atomic numbers of each atomic element
    in a formula unit over a set Q range as defined in Eqs 11 & 14 of
    Eggert et al., 2002.

    i.e.

    K_p(Q) for each element = f_p(Q)/f_e(Q)
    K_p for each element = <K_p(Q)>_Q

    Where:
        f_e(Q) is the effective electronic form factor of the total composition
        f_p(Q) are the atomic form factors for each atom, p, in the formula unit
        K_p(Q) is the Q-dependent effective atomic number for each atom, p
        K_p is the average effective atomic number


    Args:
        composition - dictionary of elements in composition where values
        are tuples of the form (Z, charge, n)

        Q - Q range to evaluate over

    Returns:
        K_p - numpy array of K_p values for each species in composition
              length is equal to number of species *not* number of atoms
              in a formula unit
    '''
    effective_ff, atomic_ff = calc_effective_ff(composition, Q)
    K_p = atomic_ff / effective_ff
    K_p = np.mean(K_p, 1)
    return K_p


def calculate_weights(composition, Q):
    '''
    Calculate X-ray scattering weightings for atomic pair correlations in
    the composition for a set Q-range. Uses the Warren-Krutter-Morningstar
    approximation for average effective atomic number.
    '''
    # Get average effective atomic number for each species
    K_p = calc_K_p(composition, Q)
    # Get counts of each species
    element_counts = np.array([element[2] for element in composition.values()])
    # Number of atoms in composition
    N = np.sum(element_counts)
    # K_p_norm - normalise to number of atoms in composition
    K_p_norm = np.dot(element_counts/N, K_p)**2
    # Dict of Kp for each species
    K_p_dict = dict(zip(composition.keys(), K_p))
    # List of atoms in composition
    atom_list = list(composition.keys())
    # List of atom (alpha-beta) pairs
    pairs = list(set(itertools.product(atom_list, repeat=2)))
    weights = {}
    c_dict = {}

    for pair in pairs:
        alpha = pair[0]
        beta = pair[1]
        c_alpha = composition[alpha][2]/N
        c_beta = composition[beta][2]/N
        K_alpha = K_p_dict[alpha]
        K_beta = K_p_dict[beta]
        # Kronecker delta - i.e. delta = 1 of a==b and delta = 0 if a!=b
        delta = int(alpha == beta)
        # Calculate weight for alpha-beta pair
        w_alphabeta = (c_alpha*K_alpha * c_beta*K_beta * (2-delta))/K_p_norm

        weights[(alpha,beta)] = w_alphabeta

    for atom in composition:
        c_beta = composition[atom][2]/N
        c_dict[atom] = c_beta

    return weights, c_dict


def calc_S_inf(composition, Q, method='ashcroft-langreth'):
    '''
    Caclulates the value S_inf as defined in Eq 19 of Eggert et al., 2002
    i.e. S_inf = SUM K^2_p / Z_tot^2
    S_inf defaults to 1 if the composition is monatomic

    For the Faber-Ziman structure factor formulation, the interference
    function, i(Q) is defined as S(Q) - 1. calc_S_inf is still called in the
    program to aid readability and avoid too many if statement blocks in the
    main code, though there is a slight performance decrease.
    '''
    if method == 'faber-ziman':
        S_inf = 1.0

    elif method == 'ashcroft-langreth':
        element_counts = np.array([element[2] for element in composition.values()])
        N = np.sum(element_counts)
        if N == 1:
            S_inf = 1.0
        elif N > 1:
            K_p = calc_K_p(composition, Q)
            Z_tot = calc_Z_sum(composition)
            S_inf = np.dot(element_counts, K_p**2) / Z_tot**2
        else:
            raise ValueError('Error in composition - n value < 1')

    else:
        raise ValueError('Please select a valid method for structure factor')
    return S_inf


def calc_alpha(Q_cor, I_cor, rho,
               Z_tot=None, J=None, S_inf=None, effective_ff=None,
               average_scattering=None, compton_scattering=None,
               method='ashcroft-langreth'):
    '''
    Calculates 'alpha' - the normalisation factor

    For the Ashcroft-Langreth formalisation this is defined in
    Eq.34 of Eggert et al., 2002.

    For the Faber-Ziman formalism this is defined in equation A13 of
    Morard, 2013 but note there is a missing bracket - refer to the definition
    in Decremps et al., 2016, PhysRev B 93, 054209 instead.

    Args:
        Q_cor - Rebinned Q space data cut at Q_max
        I_cor - Background corrected intensity data
        rho - actual or initial estimate of atomic density in atoms/Angstrom^3
        method - 'ashcroft-lamgreth' or 'faber-ziman'

        For AL:
        ------
        Z_tot - Atomic number sum of elements in composition
        J - J value as defined in Eq 35 of Eggert et al., 2002
        S_inf - S_inf value as defined in Eq 19 of Eggert et al., 2002
        effective_ff - effective electronic form factor

        For FZ:
        ------
        average_scattering - Average scattering functions <f2> and <f>2
                             see calc_average_scattering for definitions
        compton_scattering - Compton scattering intensity

    Returns:
        alpha
    '''
    if method == 'ashcroft-langreth':
        # define the two integrals in the equation
        int_1 = simpson((J + S_inf) * Q_cor**2, x=Q_cor)
        int_2 = simpson((I_cor / effective_ff**2) * Q_cor**2, x=Q_cor)
        # calculate alpha
        alpha = Z_tot**2 * (((-2 * np.pi**2 * rho) + int_1) / int_2)
        return alpha
    elif method == 'faber-ziman':
        # define the two integrals in the equation
        int_1 = simpson(((compton_scattering + average_scattering[0])/average_scattering[1]) * Q_cor**2, x=Q_cor)
        int_2 = simpson(((Q_cor**2 * I_cor) / average_scattering[1]), x=Q_cor)
        # Calculate alpha
        alpha = ((-2 * np.pi**2 * rho) + int_1)/int_2
        return alpha
    else:
        raise ValueError()


def calc_coherent_scattering(Q_sample, I_sample, composition, alpha,
                             compton_scattering=None, method='ashcroft-langreth'):
    '''
    Calculates coherent scattering intensity for a given composition

    Coherent scattering is defined as:
    Coh_I = N * [alpha * Sample_scattering(Q) - SUM Incoherent_scattering]
    where N is number of atoms, n_atoms

    For the Faber-Ziman formalism there is no need to scale by the number
    of atoms in the compositional unit (N)

    Args:
        Q_sample - The Q-values to evaluate over (np array)
        I_sample - The sample scattering intensity (i.e. bkgd corrected) (np array)
        composition - Dict of sample composition
        alpha - The normalisation factor
        method - S(Q) formalism, 'ashcroft-langreth' or 'faber-ziman'

    Optional:
        compton_scattering - The compton scattering intensity for the sample
    '''
    # If compton scattering data is not provided, retrieve it now
    if compton_scattering is None:
        compton_scattering = calc_total_compton_scattering(composition, Q_sample)
        if method == 'faber-ziman':
             compton_scattering /= sum([element[2] for element in composition.values()])
    else:
        pass
    coherent_scattering = (alpha * I_sample) - compton_scattering
    if method=='ashcroft-langreth':
        # N - number of atoms in compositional unit (molecule)
        n_atoms = sum([element[2] for element in composition.values()])
        # Normalise coherent scattering
        coherent_scattering *= n_atoms
    else:
        pass

    return coherent_scattering


def calc_structure_factor(Q_cor, I_cor, composition, rho, method='ashcroft-langreth', return_alpha=False):
    '''
    Calculates the molecular structure factor in normalised units using the
    method first developed by Krogh-Moe (1956) and Norman (1957).

    Args:
        Q_cor - Rebinned Q space data
        I_cor - Bkgd corrected intensity data
        composition - dictionary of atomic elements, with values as tuples in
                      the form (Z,charge,n)
        rho - actual or initial estimate of atomic density in atoms/Angstrom^3

    Kwargs:
        return_alpha (bool) - optionally return the Krogh-Moe-Norman normalisation factor (α) along with S(Q)
        method - S(Q) formalism (see below)

    The code provides an option to use either the Ashcroft-Langreth or
    Faber-Ziman formulation of the structure factor S(Q)

    S_AL(Q) = I^COHERENT(Q) / N * Z_TOT^2 * f_e^2(Q)

    S_FZ(Q) = I^COHERENT(Q) - (<f^2> - <f>^2 ) / <f>^2

    I^COHERENT (the coherent scattering component of measured scattering
    intensity) has a different formulation for AL & FZ methods.

    I^COH_AL = N [alpha_AL * I^SAMPLE(Q) - SUM_p[I^COMPTON_p(Q)] ]

    I^COH_FZ = alpha_FZ * I^SAMPLE(Q) - SUM_p[I^COMPTON_p(Q)]

    I^COMPTON_p is the compton (incoherent) scattering of each atom
    (see core.calc_compton_scattering)

    alpha is the normalisation factor (see core.calc_alpha)

    N is the total number of atomic species in a formula unit
    (e.g. for SiO2 N=3)

    Z_TOT is the total Z number of a formula unit

    f_e(Q) is the effective electronic form factor

    <f^2> and <f>^2 are intermediate functions which quantify average
    scattering in a similar way to the effective electronic form factor used
    in the Ashcroft-Langreth formulation
    '''
    # N is number of atoms in 1 formula unit
    n_atoms = sum([element[2] for element in composition.values()])
    if method == 'ashcroft-langreth':
        # Z is total atomic number of formula unit
        Z_tot = calc_Z_sum(composition)
        # Calculate intermediate functions for alpha/coherent scattering for
        # computational efficiency
        J = calc_J(composition, Q_cor)
        S_inf = calc_S_inf(composition, Q_cor, method=method)
        # effective electronic form factor also used to calculate alpha
        effective_ff, _ = calc_effective_ff(composition, Q_cor)
        alpha = calc_alpha(Q_cor, I_cor, rho,
                           Z_tot=Z_tot, J=J, S_inf=S_inf,
                           effective_ff=effective_ff, method='ashcroft-langreth')
        coherent_scattering = calc_coherent_scattering(Q_cor, I_cor, composition, alpha, method=method)
        structure_factor = coherent_scattering / (n_atoms * Z_tot**2 * effective_ff**2)

    elif method == 'faber-ziman':
        # Calculate average scattering functions <f^2> & <f>^2
        avg_scattering_f = calc_average_scattering(composition, Q_cor)
        # Calculate compton (incoherent scattering) used in alpha_FZ
        compton_scattering = calc_total_compton_scattering(composition, Q_cor)
        compton_scattering = compton_scattering / n_atoms
        # Calculate alpha
        alpha = calc_alpha(Q_cor, I_cor, rho,
                           average_scattering=avg_scattering_f,
                           compton_scattering=compton_scattering,
                           method='faber-ziman')
        # Calculate coherent scattering
        coherent_scattering = calc_coherent_scattering(Q_cor, I_cor,
                                                       composition, alpha,
                                                       compton_scattering=compton_scattering,
                                                       method=method)
        # Calculate structure factor
        structure_factor = (coherent_scattering -
                            (avg_scattering_f[0] - avg_scattering_f[1])
                           ) / avg_scattering_f[1]

    else:
        raise ValueError('Please select a valid method for structure factor')

    if return_alpha:
        if method == 'ashcroft-langreth':
            # The normalisation factor in AL formalism is scaled
            # by the number of atoms in the compositional unit
            alpha = alpha / n_atoms
        return structure_factor, alpha
    else:
        return structure_factor

def get_mod_func(q, mod_func, window_start):
    '''
    Returns the modification function to scale S(q) before FFT

    Implemented modification functions:
    ----------------------------------

    'lorch' - from Lorch, 1969:

        M(Q) = sin(pi*Q/Q_max) / (pi*Q/Q_max)

    'cosine-window' -  - from Drewitt et al., 2013:

        M(Q) = 1 if Q < Q1
        M(Q) = 0.5[1 + cos(x*pi / N-1)] if Q1 <= Q <= Qmax

        where N is width of the window function
        and x is an integer with values from 0 - (N-1) across the window

        if cosine-window is selected, optional kwarg window-start must
        be passed (Q1)
    '''
    if mod_func == 'Cosine-window':
        pre_window = np.ones(len(q[np.where(q<window_start)]))
        window = np.arange(len(q[np.where(q>=window_start)]))
        window_length = len(window)
        modification = 0.5 * (1 + np.cos(window*np.pi/(window_length-1)))
        modification = np.concatenate((pre_window, modification))

    elif mod_func == 'Lorch':
        with np.errstate(divide='ignore', invalid='ignore'):
            modification = np.sinc(q / np.max(q))
    else:
        modification = 1

    return modification


def calc_correlation_func(x, y, rho, dx='check', N=12, mod_func=None,
             window_start=None, function='density_func'):
    '''
    Calculates the fourier transform of the interference function
    i(Q) = S(Q) - S_inf to obtain the real-space correlations functions:
    1) pair-distribution function [g(r)], 2) the differential correlation
    function (density function) [D(r)] or, 3) the radial distribution function
    [RDF(r)]. The conversion to real-space is useful for structural (PDF) analysis.

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

    Differential correlation function (density function) D(r) = F(r):

        D(r) = 4 * pi * r * rho_0 * [g(r) - 1]
        D(r) = 2/pi * INT[ q[i(q)] sin(qr) dq]

    Radial distribution function, RDF(r):

        RDF(r) = 4 * pi * r^2 * rho_0 * g(r)
        RDF(r) = (2r/pi)*INT[ q[i(q)] sin(qr) dq] + 4*pi*r^2*rho_0

    Args:
        x   - Q values (numpy array)
        y   - Corresponding values of the interference function (numpy array)
        rho - atomic density (float)
        dq  - steps of q values (which need to be constant)
              this defaults to 0.02, but can also be calculated in function if
              set to 'check'
        N   - Sets size of array for fft where array_size = 2^N, defaults to 12
        mod_func - Modification function to scale S(Q) before FFT
        function - density_func, pair_dist_func, radial_dist_func

    Returns:
        r   - Array of r values
        D_r/g_r/RDF_r - Corresponding array of values for calculated function
    '''
    # Define dx (qstep) if not passed to function
    if dx == 'check':
        dx = x[1] - x[0]
    else:
        pass

    # Get modification function
    modification = get_mod_func(x, mod_func, window_start)

    # Calculate qi(q) input array
    z = x * y * modification
    del modification
    # Make sure using float64
    z = z.astype(np.float64, copy=False)
    # Pad with enough zeros to make array of size 2**N
    # Length of array to fourier transform is 2**N
    # The array is mirrored and padding zeros are added to increase
    # the size of the array
    # Therefore len(z) < ((2**N)/2)
    # N=12 allows max 2048 data points (at q-spacing 0f 0.02 range is 40.94)
    if len(z) > (2**N/2):
        # Error checking also takes place in LiquidDiffract.gui
        # so that the graphical app does not crash
        raise ValueError('Length of data array > [(2**N)/2]\n' +
                         '            Please check your data is in ' +
                         'Q-space or increase value of \'N\'')
    else:
        pass
    padding_zeros = int(2**N - (len(z)*2 - 1))
    # Pad array with odd image of function
    z = np.concatenate((z, np.zeros(padding_zeros), -np.flip(z[1::], 0)))
    # Fourier transform here uses scipy.fftpack over numpy.fft
    # ifft used as converting from frequency (intensity) domain
    # This means no array-length scaling is needed and the result is not negative
    fft_z = scipy.fftpack.ifft(z)
    # Only the imaginary component of the transformed array is used
    # Only the first half of the array is used as the rest is the odd component
    fft_z = np.imag(fft_z)[:len(fft_z)//2]
    # Scaling for the inverse transform is QMax = dx * (2**N)/2
    fft_scaling = 2**N / 2 * dx
    fft_z *= fft_scaling
    # Calculate steps of r and r array
    dr = np.pi*2/(len(z)*dx)
    r = np.arange(len(fft_z))*dr
    if function == 'density_func':
        # Scale by factor 2/pi
        D_r = fft_z * 2/np.pi
        return r, D_r
    elif function == 'pair_dist_func':
        with np.errstate(divide='ignore', invalid='ignore'):
            sf = 1 / (2 * r * rho * np.pi**2)
            g_r = fft_z * sf
        g_r += 1
        return r, g_r
    elif function == 'radial_dist_func':
        RDF_r = (2 * r * fft_z / np.pi) + (4 * np.pi * rho * r**2)
        return r, RDF_r
    else:
        raise ValueError('arg \'function\' must be valid option')

def normalise_achroft_langreth_func(S_inf, sq=None, gr=None, rdf=None,
                                    rho_0=None, r=None):
    '''
    Helper function to normalise functions in the Aschroft-Langreth formalism
    for polyatomic samples. For polyatomic samples 0 < S_inf < 1 and so
    S(Q) --> S_inf at infinite Q and g(r) --> (1-S_inf) at r=0. RDF(r) and T(r)
    are similarly scaled. Renormalising these functions by S_inf makes
    comparisons to functions in the Faber-Ziman formalism easier and allows the
    same curve-fitting routines to be used on both.

    S(Q) normalisation:

        S_s(Q) = S_AL(Q) / S_inf

    g(r) normalisation:

        g_s(r) = (g(r) - 1)/S_inf + 1

    RDF(r) normalisation:

        RDF_s(r) = (RDF(r) - 4πr^2ρ_0)/S_inf + 4πr^2ρ_0

    Args:
        S_inf   - S_inf value for normalisation (see core.calc_S_inf)
        sq*     - Values of S(Q) to normalise (numpy array)
        gr*     - Values of g(r) to normalise (numpy array)
        rdf**   - Values of RDF(r) to normalise
        rho_0*  - atomic number density (float)
        r*      - Values of r that correspond to the RDF(r)

        *  - Optional kwargs
        ** - The RDF(r) normalisation requires the additional kwargs rho_0 & r

    Returns 1 of:

        sq_norm  - Renormalised S(Q) (numpy array)
        gr_norm  - Renormalised g(r) (numpy array)
        rdf_norm - Renormalised RDF(r) (numpy array)
    '''
    # S(Q) normalisation
    if sq is not None:
        sq_norm = sq / S_inf
        return sq_norm
    elif gr is not None:
        gr_norm = (gr - 1)/S_inf + 1
        return gr_norm
    elif rdf is not None:
        try:
            constant = 4*np.pi*rho_0*r**2
            rdf_norm = (rdf - constant)/S_inf + constant
            return rdf_norm
        except TypeError:
            raise NameError('RDF(r) normalisation requires keyword arguments rho_0 and r')
    else:
        raise NameError('Must provide function to normalise')


def calc_D_r_iteration_term(delta_D_r, N=12, dq=0.02):
    '''
    Calculates the term responsible for oscillations at small r values.
    Equation 45 in Eggert et al., 2002

    i.e. Delta_alpha * Q * S_inf = INT 0>r_min [Delta_D(r) sin(Qr) dr]

    Delta_D(r) is the difference between D(r) and its expected behaviour

    N defines the size of the array to be fourier transormed (2**N)
    This should be the same as used to calculate D(r) from i(Q) to
    preserve correct scaling factors
    '''
    # Make sure array size is not too large for padding value N
    if len(delta_D_r) > (2**N/2):
        raise ValueError('Length of array > [(2**N)/2]')
    else:
        pass
    padding_zeros = int(2**N - (len(delta_D_r)*2 - 1))
    z = np.concatenate((delta_D_r,
                        np.zeros(padding_zeros),
                        -np.flip(delta_D_r[1::], 0)))
    # For the forward fourier transform we need the negative of the imaginary components
    # in the first half of the array
    fft_z = -np.imag(scipy.fftpack.fft(z))
    # The fourier transform of Qi(Q) was scaled by the artificial Q_max to obtain the corrected
    # magnitudes - i.e. fft_scaling = len(D(r) * dQ = 2**N / 2 * q_step
    # This scaling has to be reversed here. For N=12 and q_step = 0.02 the scaling is 40.96
    fft_scaling = 2**N / 2 * dq
    fft_z = fft_z / fft_scaling
    # For the forward transform here
    return fft_z


@data_cache(maxsize=128)
def calc_model_D_intra_r(rho, r):
    '''
    Returns the model behaviour of D(r) at r < r_min, where r_min is the
    largest distance (0 -- r-min) where no atom can be found.
    In a liquid, this should be the distance of the 1st coordination shell.

    The model behaviour is:

         D_(r<r_min) = D_intra(r) − 4πrρ

         D_intra(r) is the expected contribution of *intra*molecular forces
         i.e. interaction of atoms within the same molecule

         This can be calculated using a fixed model approach
         (see eq. 42 in Eggert et al., 2002)

         There are also further intermolecular contributions at r<r_min
         if the sample has intermolecular coordination between compositional
         units (e.g. polymeric liquids or framework glasses)

     LiquidDiffract assumes D_intra(r) is negligible. This is only true for
     some compositions (e.g. monatomic metallic liquids). This functionality
     may be added in future releases however.
    '''
    model_D_intra_r = -4 * np.pi * r * rho
    return model_D_intra_r


def calc_chi_squared(r_intra, delta_D_r, method='simpson'):
    '''
    Calculates a chi squared figure of merit as in equation 50 of Eggert. This
    is a measure of how small Delta_D_i(r) is, to be used as a tolerance factor
    to check convergence of i_i+1(q) and in the optimisation routine of rho/s.

    Chi^2 = INT 0-r_min [Delta_D(r)]^2 dr

    The value chi^2 is defined as the area under the curve Delta_D(r) limited
    to r_min (squared to remove negatives).

    The Chi^2 value used in LiquidDiffract is scaled by 1e6 to enable easier minimisation.

    Args:
        delta_D_r - numpy array of delta_D_r from r=0 to r=r_min

    Returns:
        chi_squared - Value of chi^2
    '''
    # Square to remove negatives
    f = delta_D_r ** 2
    if method == 'simpson':
        # Integrate using simpson's rule
        chi_squared = simpson(f, x=r_intra)
    else:
        raise NotImplementedError()
    return chi_squared * 1e6


def stop_iteration(stop_condition='count', count=None,
                   iter_limit=5, chi_squared=None):
    '''Helper function for checking stop condition of iteration

    Currently LiquidDiffract uses an iteration limit only so this function
    is very unnecessary/a little expensive. It is left here to make any changes
    easier to implement (i.e. checking for chi^2 convergence)
    '''
    if stop_condition == 'count':
        if count < iter_limit:
            return False
        else:
            return True
    else:
        raise NotImplementedError()


def calc_impr_interference_func(q_data, interference_func,
                                composition, rho,
                                r_min, iter_limit,
                                method, mod_func, window_start, fft_N,
                                return_alpha=False, alpha=None):
    '''
    Calculates an improved estimate of the interference function via the
    iterative procedure described in Eggert et al., 2002. This is done by
    scaling the data and forcing its behaviour in the low r region to some
    modelled behaviour.

    The function also returns a chi^2 figure of merit which can be passed
    to a minimisation routine to estimate rho/b.

    Args:
        q_data - Numpy array of Q space data
        interference_func - Initial interference function
        composition - Dict of elements, values as tuples of form (Z,charge,n)
        rho - Density (average number density in atoms/angstrom^3)
        r_min - Intramolecular distance cut-off
        iter_limit - Iteration limit for S(Q) refinement. iter_limit >= 3 is
                     recommended for convergence. If optimising for minimum
                     chi^2 then 3 <= iter_limit <= 10 is recommended. Please
                     see the LiquidDiffract README for more information
        method - Formalisation of total scattering structure factor to use
                 ('ashcroft-langreth' or 'faber-ziman')
        mod_func - Modification function to scale data before fft
                   None, 'Cosine-window' or 'Lorch'
        window_start - Start of Cosine-window function needed if
                       mod_func == 'Cosine-window'
        fft_N - Sets size of array for fft where array_size = 2^N

    Kwargs:
        return_alpha (bool) - Optionally return the modified Q-dependent normalisation factor
        alpha - Krogh-Moe-Norman normalisation factor (α_0), required if return_alpha==True

    Returns:
        interference_func_impr - Interference function after iter_limit iterations
        chi_sq - Figure of merit (see calc_chi_squared)
    '''
    # Unpack rho - solver passes rho as a single element array
    rho = np.float64(rho)

    dq = q_data[1] - q_data[0]
    # Calculate initial D(r)
    r, D_r = calc_correlation_func(q_data, interference_func, rho,
                      mod_func=mod_func, window_start=window_start,
                      dx=dq, N=fft_N)
    # Calculate expected behaviour of D(r) for intramolecular distances (r<r_min)
    # and scale by S_inf (1 for FZ/monatomic - only necessary for AL-S(Q))
    model_D_intra_r = (np.copy(calc_model_D_intra_r(rho, r)) *
                       calc_S_inf(composition, q_data, method=method))

    # Calculate static terms of iterative proceduce
    with np.errstate(divide='ignore', invalid='ignore'):
        t1 = 1 / q_data

    if method == 'ashcroft-langreth':
        t2_divisor = calc_S_inf(composition, q_data, method=method) + calc_J(composition, q_data)
    elif method == 'faber-ziman':
        t2_divisor = 1.0

    # Count no. iterations for verbosity or control
    # Function set up to call stop_looping func in case want to modify
    # to use chi-squared value as tolerance value instead
    count = 1
    done_looping = False
    while not done_looping:
        try:
            done_looping = stop_iteration(stop_condition='count', count=count, iter_limit=iter_limit)
            interference_func = interference_func_impr
            r, D_r = calc_correlation_func(q_data, interference_func, rho,
                              mod_func=mod_func, window_start=window_start,
                              dx=dq, N=fft_N)
            if return_alpha:
                alpha = alpha_impr
        except NameError:
            pass

        # Calculate the difference between observed and model D(r) along
        # with the chi^2 figure of merit
        delta_D_r = (D_r - model_D_intra_r)[np.where(r < r_min)]
        r_intra = r[np.where(r<r_min)]
        chi_squared = calc_chi_squared(r_intra, delta_D_r)

        # Calculate next iteration of i(Q)
        t3 = calc_D_r_iteration_term(delta_D_r, N=fft_N, dq=dq)[:len(interference_func)]
        t2 = (interference_func/t2_divisor) + 1
        with np.errstate(invalid='ignore'):
            interference_func_impr = interference_func - (t1 * t2 * t3)

        # Calculate next iteration of alpha
        # α(n+1) = α(n)*[1 + Δα(n)]
        #   Δα(n)_FZ = (1/Q) * ∫ΔD(r)sin(Qr)dr
        #   Δα(n)_AL = (1/Q) * [S_inf+J(Q)]^-1 * ∫ΔD(r)sin(Qr)dr
        #   This is: t1 * 1/t2_divisor * t3
        if return_alpha:
            with np.errstate(invalid='ignore'):
                delta_alpha = t1 * (1/t2_divisor) * t3
            alpha_impr = alpha * (1 + delta_alpha)

        count += 1

    if return_alpha:
        return interference_func_impr, chi_squared, alpha_impr
    else:
        return interference_func_impr, chi_squared


def refinement_objfun(minvar, *args):
    '''
    Wrapper objective function to pass to a solver to refine density (rho)
    and/or background scaling factor (b) by minimising a chi^2 figure of merit.

    This function provides set-up and logic for dealing with refinement of
    different variables, including recalculating the background subtraction
    when refining b. The actual calculations for the optimisation objective
    function are in core.calc_impr_interference_func

    Args:

        minvar - independent variable(s) in the optimisation (numpy array with shape (n,))
                 Can be one or both of:
                    rho - average (bulk) number density (atoms/angstrom^3)
                    bkg_scale - background scaling factor

        *args - The expected arguments depend on if the background scaling factor
                is being refined. The last two arguments should be opt_rho and
                opt_bkg, i.e.:

                    *_, opt_rho, opt_bkg = args

                These should be boolean flags indicating if the density and
                background scale factor are being refined respectively.

                If opt_bkg == True (1), the expected args are:

                    (*_, q_data, I_data_uncorrected, I_bkg, data_correction,
                     qmin, qmax,
                     smooth_flag, smooth_window_length, smooth_poly_order,
                     composition, r_min, iter_limit,
                     method, mod_func, window_start, fft_N,
                     opt_rho, opt_bkg) = args

                If opt_rho == False (0), the density should precede q_data

                If opt_bkg == False (0), the expected args are:

                    (q_data, I_data_corrected,
                     composition, r_min, iter_limit,
                     method, mod_func, window_start, fft_N,
                     opt_rho, opt_bkg) = args

            The individual arguments are:

                q_data - Numpy array of Q space data

                I_data_uncorrected - Numpy array of intensity (I(Q)) values
                                     prior to background subtraction
                I_data_corrected - Numpy array of background subtracted I(Q) values
                I_bkg - Numpy array of background intensities
                data_correction - Multiplicative data correction (e.g. attenuation correction)
                qmin - Minimum Q value of useable data (float or None)
                qmax - Maximum Q value to truncate data (float or None)
                smooth_flag - Boolean flag to apply savitzky golay filter
                smooth_window_length - Window length for optional smoothing filter
                smooth_poly_order - Polynomial order for optional smoothing filter
                composition - Dict of elements, values as tuples of form (Z,charge,n)
                r_min - Intramolecular distance cut-off
                iter_limit - Iteration limit for S(Q) refinement. iter_limit >= 3 is
                             recommended for convergence. If optimising for minimum
                             chi^2 then 3 <= iter_limit <= 10 is recommended. Please
                             see the LiquidDiffract README for more information
                method - Formalisation of total scattering structure factor to use
                         ('ashcroft-langreth' or 'faber-ziman')
                mod_func - Modification function to scale data before fft
                           None, 'Cosine-window' or 'Lorch'
                window_start - Start of Cosine-window function needed if
                               mod_func == 'Cosine-window'
                fft_N - Sets size of array for fft where array_size = 2^N

    Returns:
        chi_sq - Figure of merit (see calc_chi_squared)
                 This is the value to be minimised by optimising
                 the independent variable(s) minvar
    '''
    *_, opt_rho, opt_bkg = args
    # Unpack optimisation variables (rho and/or bkg_scale)
    if opt_rho == True and opt_bkg == True:
        rho = np.float64(minvar[0])
        bkg_scale = np.float64(minvar[1])
    # Unpack rho - solver passes rho as a single element array
    elif opt_rho == True and opt_bkg == False:
        rho = np.float64(minvar)
    # Unpack bkg_scale as abovs
    elif opt_rho == False and opt_bkg == True:
        bkg_scale = np.float64(minvar)
    else:
        raise NotImplementedError()

    if opt_rho == False:
        rho, *_ = args
        rho = np.float64(rho)

    if opt_bkg == True:
        # Unpack args - additional args needed when refining bkg_scale
        (*_,
         q_data, I_data_uncorrected, I_bkg, data_correction,
         qmin, qmax, smooth_flag, smooth_window_length, smooth_poly_order,
         composition, r_min, iter_limit,
         method, mod_func, window_start, fft_N, _, _) = args
        I_data_corrected = (I_data_uncorrected - (I_bkg * bkg_scale)) * data_correction
        # Cut at qmax if qmax is != None
        if qmax:
            # Cut data at qmax
            cut_max = np.where(q_data < qmax)
            q_data = q_data[cut_max]
            I_data_corrected = I_data_corrected[cut_max]
        # Similarly cut at qmin if present
        if qmin:
            # Cut data to qmin - setting I(q<qmin)=I(qmin)
            fill_val = I_data_corrected[np.argmax(q_data > qmin)]
            cut_min = I_data_corrected[np.where(q_data > qmin)]
            padding = np.asarray([fill_val] * (len(q_data) - len(cut_min)))
            I_data_corrected = np.concatenate((padding, cut_min))
        # smooth data if smooth_flag == 1
        if smooth_flag == True:
            I_data_corrected = smooth_data(I_data_corrected,
                                      window_length=smooth_window_length,
                                      poly_order=smooth_poly_order)
    else:
        # Unpack args here for opt_bkg == False
        # (if opt_bkg is 0, opt_rho is always 1)
        (q_data, I_data_corrected, composition, r_min, iter_limit,
         method, mod_func, window_start, fft_N, _, _) = args

    structure_factor = calc_structure_factor(q_data, I_data_corrected,
                                             composition, rho,
                                             method=method)
    interference_func = structure_factor - calc_S_inf(composition, q_data, method=method)

    _, chi_squared = calc_impr_interference_func(q_data, interference_func,
                                                 composition, rho, r_min, iter_limit,
                                                 method, mod_func, window_start, fft_N)

    return chi_squared


def chi_squared_map_helper(param, arg_dict):
    '''
    Helper function to compute χ^2 maps from a 2D parameter space with a choice
    from rho, background scaling factor (b), no. iterations, and r-min cutoff).

    Args:

        param - ndarray with shape (4,) in format [rho, b, r_min, n_iter]
                To call refinement_objfun with no background scaling set
                b to None
        arg_dict - Dictionary object of argument to pass to refinement_objfun
                   see refinement_objfun for details

    Returns:

        chi_squared - see refinement_objfun for details
    '''
    # if bkg_scaling is nan - no bkg_scaling
    if np.isnan(param[1]):
        objfun_var = param[0]
    else:
        objfun_var = param[:2]
    arg_dict['r_min'] = param[2]
    arg_dict['n_iter'] = param[3]
    objfun_arg = tuple(arg_dict.values())
    return refinement_objfun(objfun_var, *objfun_arg)


def integrate_coordination_sphere(r, rdf,
                                  r_0=None, rp_max=None,
                                  r_max=None, r_min=None, method=0):
    '''
    Calculates the 1st coordination number, N_1, for a monatomic sample by 
    integrating the area under the first peak in RDF(r)

    The position of the 1st peak in RDF(r) directly gives the average 
    interatomic bond length, r_ab, between an atom 'a' at the origin and the
    next-nearest neighbour. The area under this curve gives the average 
    coordination number of the 1st coordination sphere, N_1. RDF(r) is the 
    radial distribution function, defined as:

        RDF(r) = 4πρr^2g(r), where g(r) is the pair-distribution function

    The coordination number is fairly sensitive to the upper integration limit 
    chosen. There are several methods used in the literature to define the 
    limits based on either physical models or empirical observations. Three 
    methods are implimented in LiquidDiffract:

        1) Integration of area symmetrical in rg(r)

        This method is based on the assumption that the quantity rg(r) is 
        symmetrical for a coordination shell about its average position. The
        coordination number is calculated by integrating across RDF(r) from 
        the leading edge of the 1st peak, r_0, to the peak centre of the
        function r*g(r), r'_max. The integral is doubled to account for the 
        right half of the peak:

            (Eq. 2.5.1 in Waseda, 1980)
            N_a = 2 * INT 4πρr[rg(r)]_sym dr | r=r_0 --> r=r'_max 

        2) Integration of area symmetrical in r^2g(r)

        This method is similar to method 1 but assumes the peak is symmetrical
        in r^2g(r) instead of rg(r). The limit r_max is therefore the peak 
        centre in RDF(r):

            (Eq. 2.5.2 in Waseda, 1980)
            N_b = 2 * INT 4πρ[r^2g(r)]_sym dr | r=r_0 --> r=r_max

        It should be noted that r_max > r'_max

        3) Integration to the firt minimum in RDF(r)

        This is the most commonly used method to estimate N_1. It is based on
        the observation that the 1st peak is generally assymetric in 
        experimental data. The function 4πρr^2g(r) is integrated from r_0 to 
        the first minimum after the peak:

            (Eq. 2.5.3 in Waseda, 1980)
            N_c = INT RDF(r) dr | r=r_0 --> r=r_min 

        N_c is sensitive to the value of r_min used and r_min is not always
        easily refined. Noise and truncation ripples in RDF(r) can make it
        difficult to determine r_min, as can the broadening of RDF(r) 
        associated with an increase in temperature.

    It should be noted that: N_a < N_b but N_c is not always > N_b
    N_a can be considered a lower-bound on the coordination numeber as the 
    peak is never truly symmetrical. In some cases N_c can provide an upper
    limit, but the value is much more sensitive as r_min is never as well
    defined as r_0 and r_max/r'_max

    Args:
        r   - r values (numpy array)
        rdf - Corresponding values of the radial distribution function, RDF(r) (numpy array)
        r_0 - Leading edge of the 1st peak in RDF(r)
        rp_max - r'_max, peak centre in rg(r)
        r_max - peak centre in r^2g(r)
        r_min - 1st minimum to the right of the 1st peak in RDF(r)

    Returns:
        N_a - 1st coordination number computed via method 1.
        N_b - 1st coordination number computed via method 2.
        N_c - 1st coordination number computed via method 3.
    '''

    rdf_interp = scipy.interpolate.interp1d(r, rdf, kind='cubic', fill_value='extrapolate')

    if method == 0:
        N_a = integrate_coordination_sphere(r, rdf, r_0=r_0, rp_max=rp_max, method=1)
        N_b = integrate_coordination_sphere(r, rdf, r_0=r_0, r_max=r_max, method=2)
        N_c = integrate_coordination_sphere(r, rdf, r_0=r_0, r_min=r_min, method=3)
        return N_a, N_b, N_c

    elif method == 1:
        N_a, _ = quad(rdf_interp, r_0, rp_max, limit=200)
        N_a *= 2.0
        return N_a

    elif method == 2:
        N_b, _ = quad(rdf_interp, r_0, r_max, limit=200)
        N_b *= 2.0
        return N_b

    elif method == 3:
        N_c, _ = quad(rdf_interp, r_0, r_min, limit=200)
        return N_c

    else:
        raise NameError('\'method\' must be one of either 1, 2, 3, or 0')
