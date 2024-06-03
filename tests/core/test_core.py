#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for LiquidDiffract.core.core
"""
from __future__ import absolute_import
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2023-2024, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

import unittest
from importlib import resources
import os.path
import warnings
import math
from packaging.version import parse as parse_version
import numpy as np
from scipy import __version__ as sp_ver
from platform import python_version
# import core module to test
import LiquidDiffract.core.core as core
import LiquidDiffract.core.data_utils as data_utils
# import data path and custom assertions from tests/util
from util import data_path
from util import CustomAssertions

class TestCalcMolMass(unittest.TestCase, CustomAssertions):
    def test_calc_mol_mass(self):
        composition_Ga = {'Ga': (31,0,1)}
        mass_Ga = core.calc_mol_mass(composition_Ga)
        composition_GaGaGa = {'Ga': (31,0,3)}
        mass_GaGaGa = core.calc_mol_mass(composition_GaGaGa)
        composition_H2O = {'H': (1,0,2), 'O': (8,0,1)}
        mass_H2O = core.calc_mol_mass(composition_H2O)
        self.assertFloatEqual(mass_Ga, 69.723)
        self.assertFloatEqual(mass_GaGaGa, 209.169)
        self.assertFloatEqual(mass_H2O, 18.0152)

    def test_bad_element_key(self):
        composition = {'BAD_SYMBOL': (1, 0, 2)}
        self.assertRaises(KeyError, core.calc_mol_mass, composition)


class TestConvDensity(unittest.TestCase):
    def test_conv_density(self):
        composition_H = {'H': (1,0,1)}
        rho_H = 1.0
        mass_density_H = core.conv_density(rho_H, composition_H)

        composition_HH = {'H': (1,0,2)}
        rho_HH = 1.0
        mass_density_HH = core.conv_density(rho_HH, composition_HH)

        composition_CaSiO3 = {'Ca': (20,0,1), 'Si': (14,0,1), 'O': (8,0,3)}
        rho_CaSiO3 = 100
        mass_density_CaSiO3 = core.conv_density(rho_CaSiO3, composition_CaSiO3)

        self.assertEqual(mass_density_H, 1.6736573146878266)
        self.assertEqual(mass_density_HH, 1.6736573146878266)
        self.assertEqual(mass_density_CaSiO3, 3857.8207935623163)

    def test_zero_density(self):
        composition = {'Ga': (31,0,1)}
        rho = 0.0
        mass_density = core.conv_density(rho, composition)
        self.assertEqual(mass_density, 0)


class TestAttenuationFactor(unittest.TestCase, CustomAssertions):

    def test_calc_self_shielding(self):
        # two_theta --> [0, 30, 45, 60]
        Q = np.array([0, (math.sqrt(6) - math.sqrt(2))/4, math.sqrt(2 - math.sqrt(2))/2, 0.5])
        expected_attenuation = np.array([1.0/math.e, 0.34083619263158144144, 0.30120381868339322352, (math.e-1)/math.e**2])
        attenuation = core.calc_self_shielding(Q, 1.0, 90.0, 1.0, np.pi*4)
        self.assertFloatArrayEqual(attenuation, expected_attenuation)


class TestCalcZSum(unittest.TestCase):
    def test_calc_Z_sum(self):
        composition_H = {'H': (1,0,1)}
        Z_tot_H = core.calc_Z_sum(composition_H)

        composition_H100 = {'H': (1,0,100)}
        Z_tot_H100 = core.calc_Z_sum(composition_H100)

        composition_H0 = {'H': (1,0,0)}
        Z_tot_H0 = core.calc_Z_sum(composition_H0)

        composition_CaSiO3 = {'Ca': (20,0,1), 'Si': (14,0,1), 'O': (8,0,3)}
        Z_tot_CaSiO3 = core.calc_Z_sum(composition_CaSiO3)

        self.assertEqual(Z_tot_H, 1.0)
        self.assertEqual(Z_tot_H100, 100.0)
        self.assertEqual(Z_tot_H0, 0.0)
        self.assertEqual(Z_tot_CaSiO3, 58)

@unittest.skipIf(parse_version(python_version()) < parse_version('3.9'), 'Tests unsupported for python<3.9 - Suppport will be removed in a future version of LiquidDiffract')
class TestLoadData(unittest.TestCase):
    def setUp(self):
        with resources.files('LiquidDiffract.resources').joinpath('pt_data.npy').open('rb') as fp:
            self.element_dict = np.load(fp, allow_pickle=True).item()

    def test_load_ff_data(self):
        ff_data = core.load_ff_data()
        self.assertEqual(ff_data.shape, (210, 11))
        # Z column ([:,10]) is used to look up elements in ff_data
        ff_lookup = ff_data[:,10]
        for atomic_number in iter(self.element_dict.values()):
            self.assertIn(atomic_number, ff_lookup)

    def test_load_compton_data(self):        
        for element in self.element_dict:
            try:
                self.assertTrue(core.load_compton_data(element))
            except FileNotFoundError:
                self.fail("core.load_compton_data() raised FileNotFoundError!")

    def test_load_mass_data(self):
        with resources.files('LiquidDiffract.resources').joinpath('mass_data.npy').open('rb') as fp:
            mass_dict = np.load(fp, allow_pickle=True).item()
        for element in self.element_dict:
            self.assertIn(element, mass_dict)


class TestCalcAtomicFF(unittest.TestCase, CustomAssertions):
    def test_calc_atomic_ff(self):
        expected_atomic_ff_Ga = np.load(os.path.join(data_path, 'atomic_ff_Ga_0-12.npy'))
        expected_atomic_ff_Si = np.load(os.path.join(data_path, 'atomic_ff_Si_0-12.npy'))
        element_Ga = (31,0,1) # Gallium
        element_Si = (14,0,1) # Silicon
        Q = np.arange(0, 12, 0.02)
        atomic_ff_Ga = core.calc_atomic_ff(element_Ga,Q)
        atomic_ff_Si = core.calc_atomic_ff(element_Si,Q)
        self.assertFloatArrayEqual(atomic_ff_Ga, expected_atomic_ff_Ga)
        self.assertFloatArrayEqual(atomic_ff_Si, expected_atomic_ff_Si)


class TestCalcEffectiveFF(unittest.TestCase, CustomAssertions):
    def test_calc_effective_ff(self):
        Q = np.arange(0, 12, 0.02)
        Q_long = np.arange(0, 24, 0.02)
        composition_Ga = {'Ga': (31,0,1)}
        composition_CaSiO3 = {'Ca': (20,0,1), 'Si': (14,0,1), 'O': (8,0,3)}
        composition_GaSn = {'Ga': [31, 0, 915], 'Sn': [50, 0, 85]}
        expected_aff_Ga = np.load(os.path.join(data_path, 'atomic_ff_Ga_0-12.npy'))
        expected_eff_CaSiO3 = np.load(os.path.join(data_path, 'effective_ff_CaSiO3_0-12.npy'))
        expected_eff_GaSn = np.load(os.path.join(data_path, 'effective_ff_Ga915Sn85_0-24.npy'))

        eff_Ga, aff_Ga = core.calc_effective_ff(composition_Ga, Q)
        self.assertEqual(len(aff_Ga), 1)
        self.assertEqual(len(aff_Ga[0]), 600)
        self.assertFloatArrayEqual(aff_Ga[0], expected_aff_Ga)
        self.assertFloatArrayEqual(aff_Ga/eff_Ga/31.0, 1.0)

        eff_CaSiO3, _ = core.calc_effective_ff(composition_CaSiO3, Q)
        self.assertFloatArrayEqual(eff_CaSiO3, expected_eff_CaSiO3)

        eff_GaSn, _ = core.calc_effective_ff(composition_GaSn, Q_long)
        self.assertFloatArrayEqual(eff_GaSn, expected_eff_GaSn)


class TestCalcAverageScattering(unittest.TestCase, CustomAssertions):
    def test_calc_average_scattering(self):
        Q = np.arange(0, 12, 0.02)
        composition_Ga = {'Ga': (31,0,1)}
        composition_CaSiO3 = {'Ca': (20,0,1), 'Si': (14,0,1), 'O': (8,0,3)}
        composition_GaSn = {'Ga': [31, 0, 915], 'Sn': [50, 0, 85]}
        expected_avg_scattering_CaSiO3 = np.load(os.path.join(data_path, 'average_scattering_functions_CaSiO3_0-12.npy'))
        expected_avg_scattering_GaSn = np.load(os.path.join(data_path, 'average_scattering_functions_Ga915Sn85_0-12.npy'))
        average_scattering_Ga = core.calc_average_scattering(composition_Ga, Q)
        average_scattering_CaSiO3 = core.calc_average_scattering(composition_CaSiO3, Q)
        average_scattering_GaSn = core.calc_average_scattering(composition_GaSn, Q)
        # <f2> == <f>2 for monatomic case
        self.assertFloatArrayEqual(average_scattering_Ga[0], average_scattering_Ga[1])
        self.assertEqual(len(average_scattering_CaSiO3), 2)
        self.assertFloatArrayEqual(average_scattering_CaSiO3[0], expected_avg_scattering_CaSiO3[0])
        self.assertFloatArrayEqual(average_scattering_CaSiO3[1], expected_avg_scattering_CaSiO3[1])
        self.assertFloatArrayEqual(average_scattering_GaSn[0], expected_avg_scattering_GaSn[0])
        self.assertFloatArrayEqual(average_scattering_GaSn[1], expected_avg_scattering_GaSn[1])


class TestCalcTotalComptonScattering(unittest.TestCase, CustomAssertions):
    def test_calc_compton_scattering(self):
        Q = np.arange(0, 12, 0.02)
        composition_Ga = {'Ga': (31,0,1)}
        composition_GaGa = {'Ga': (31,0,2)}
        composition_CaSiO3 = {'Ca': (20,0,1), 'Si': (14,0,1), 'O': (8,0,3)}
        expected_compton_CaSiO3 = np.load(os.path.join(data_path, 'compton_scattering_CaSiO3_0-12.npy'))

        compton_Ga = core.calc_total_compton_scattering(composition_Ga, Q)
        compton_GaGa = core.calc_total_compton_scattering(composition_GaGa, Q)
        compton_CaSiO3 = core.calc_total_compton_scattering(composition_CaSiO3, Q)

        self.assertFloatArrayEqual(compton_Ga, compton_GaGa/2.0)
        self.assertFloatArrayEqual(compton_CaSiO3, expected_compton_CaSiO3)


class TestCalcJ(unittest.TestCase, CustomAssertions):
    def test_calc_j(self):
        composition_CaSiO3 = {'Ca': (20,0,1), 'Si': (14,0,1), 'O': (8,0,3)}
        Q = np.array([0, 1, 50])
        expected_J_CaSiO3 = np.array([0, 0.002139397870661679, 5.261489234552006])
        calc_J_CaSiO3 = core.calc_J(composition_CaSiO3, Q)
        self.assertFloatArrayEqual(calc_J_CaSiO3, expected_J_CaSiO3)


class TestCalcKP(unittest.TestCase, CustomAssertions):
    def test_calc_k_p(self):
        Q = np.arange(0, 12, 0.02)
        composition_CaSiO3 = {'Ca': (20,0,1), 'Si': (14,0,1), 'O': (8,0,3)}
        composition_Ga = {'Ga': (31,0,1)}
        composition_GaGa = {'Ga': (31,0,2)}
        expected_K_CaSiO3 = np.array([22.28649729, 15.07790776,  6.87853165])
        K_Ga = core.calc_K_p(composition_Ga, Q)
        K_GaGa = core.calc_K_p(composition_GaGa, Q)
        K_CaSiO3 = core.calc_K_p(composition_CaSiO3, Q)
        self.assertEqual(len(K_Ga), 1)
        self.assertEqual(len(K_GaGa), 1)
        self.assertEqual(len(K_CaSiO3), 3)
        self.assertEqual(K_Ga, 31.0)
        self.assertFloatArrayEqual(expected_K_CaSiO3, K_CaSiO3)


class TestCalcWeights(unittest.TestCase, CustomAssertions):
    def test_calculate_weights(self):
        Q = np.arange(0, 12, 0.02)
        composition_Ga = {'Ga': (31, 0, 1)}
        composition_H2O = {'H': (1,0,2), 'O':(8,0,1)}
        composition_test = {'Mg': (12,0,1), 'Al': (13,0,1), 'Pd': (46,0,1), 'W': (74,0,1), 'Sn': (50,0,1), 'N': (7,0,5)}
        weights_Ga, c_dict_Ga = core.calculate_weights(composition_Ga, Q)
        weights_H2O, c_dict_H2O = core.calculate_weights(composition_H2O, Q)
        weights_test, c_dict_test = core.calculate_weights(composition_test, Q)
        self.assertEqual(weights_Ga[('Ga', 'Ga')], 1.0)
        self.assertEqual(c_dict_Ga['Ga'], 1.0)
        self.assertFloatEqual(weights_H2O[('O', 'H')], 0.14955153186983297)
        self.assertFloatEqual(weights_H2O[('H', 'O')], 0.14955153186983297)
        self.assertFloatEqual(weights_H2O[('H', 'H')], 0.006626295626438653)
        self.assertFloatEqual(weights_H2O[('O', 'O')], 0.8438221725037282)
        self.assertEqual(len(weights_test), 36)
        self.assertEqual(sum(c_dict_test.values()), 1.0)


class TestCalcSInf(unittest.TestCase, CustomAssertions):
    def test_calc_s_inf(self):
        Q = np.arange(0, 12, 0.02)
        composition_Ga = {'Ga': (31,0,1)}
        composition_H2O = {'H': (1,0,2), 'O': (8,0,1)}
        expected_S_inf_H2O = 0.8471353203169478
        S_inf_Ga = core.calc_S_inf(composition_Ga, Q)
        S_inf_H2O = core.calc_S_inf(composition_H2O, Q)
        self.assertFloatEqual(S_inf_H2O, expected_S_inf_H2O)
        self.assertEqual(S_inf_Ga, 1.0)

    def test_calc_s_inf_fz(self):
        test_1 = core.calc_S_inf({'Ga': (31,0,1)}, np.arange(0,12,0.02))
        test_2 = core.calc_S_inf('ARG', {'ARG2': (0,0,False)}, method='faber-ziman')
        self.assertEqual(test_1, test_2)
        self.assertEqual(test_1, 1.0)


class TestCalcAlpha(unittest.TestCase, CustomAssertions):
    @unittest.skipIf(parse_version(sp_ver) < parse_version('1.11.0'), 'scipy.integrate.simpson behaviour changed in scipy=1.11.0 - Test data assumes scipy>1.11.0 and your version is earlier')
    @unittest.skipIf(parse_version(python_version()) < parse_version('3.9'), 'Tests unsupported for python<3.9 - Suppport will be removed in a future version of LiquidDiffract')
    def test_calc_alpha(self):
        with resources.files('LiquidDiffract.scripts').joinpath('example_data.dat').open('r') as fp:
            q_test, I_test = np.loadtxt(fp, unpack=True, skiprows=0)
        q_test, I_test = data_utils.rebin_data(q_test, I_test, dx=0.02, extrapolate_mode='extrapolate')
        composition_Ga = {'Ga': (31,0,1)}
        J_Ga = core.calc_J(composition_Ga, q_test)
        Z_tot_Ga = core.calc_Z_sum(composition_Ga)
        S_inf_Ga = core.calc_S_inf(composition_Ga, q_test)
        eff_Ga, _ = core.calc_effective_ff(composition_Ga, q_test)
        avg_scattering_Ga = core.calc_average_scattering(composition_Ga, q_test)
        compton_Ga = core.calc_total_compton_scattering(composition_Ga, q_test) / 1.0
        alpha_Ga_al = core.calc_alpha(q_test, I_test, 0.05, Z_tot=Z_tot_Ga, J=J_Ga, S_inf=S_inf_Ga, effective_ff=eff_Ga, method='ashcroft-langreth')
        alpha_Ga_fz = core.calc_alpha(q_test, I_test, 0.05, average_scattering=avg_scattering_Ga, compton_scattering=compton_Ga, method='faber-ziman')
        composition_CaSiO3 = {'Ca': (20,0,1), 'Si': (14,0,1), 'O': (8,0,3)}
        J_CaSiO3 = core.calc_J(composition_CaSiO3, q_test)
        Z_tot_CaSiO3 = core.calc_Z_sum(composition_CaSiO3)
        S_inf_CaSiO3 = core.calc_S_inf(composition_CaSiO3, q_test)
        eff_CaSiO3, _ = core.calc_effective_ff(composition_CaSiO3, q_test)
        avg_scattering_CaSiO3 = core.calc_average_scattering(composition_CaSiO3, q_test)
        compton_CaSiO3 = core.calc_total_compton_scattering(composition_CaSiO3, q_test) / 5.0
        alpha_CaSiO3_al = core.calc_alpha(q_test, I_test, 0.12, Z_tot=Z_tot_CaSiO3, J=J_CaSiO3, S_inf=S_inf_CaSiO3, effective_ff=eff_CaSiO3, method='ashcroft-langreth')
        alpha_CaSiO3_fz = core.calc_alpha(q_test, I_test, 0.12, average_scattering=avg_scattering_CaSiO3, compton_scattering=compton_CaSiO3, method='faber-ziman')

        self.assertFloatEqual(alpha_Ga_al, alpha_Ga_fz)
        self.assertFloatEqual(alpha_Ga_fz, 12.385819924876388)
        self.assertFloatEqual(alpha_CaSiO3_al, 10.54252314476738)
        self.assertFloatEqual(alpha_CaSiO3_fz, 2.20134540142743)


class TestCalcCoherentScattering(unittest.TestCase, CustomAssertions):
    def test_calc_coherent_scattering(self):
        composition_H2O = {'H': (1,0,2), 'O': (8,0,1)}
        fz_test = core.calc_coherent_scattering(None, 2, composition_H2O, 3, compton_scattering=3.5, method='faber-ziman')
        al_test = core.calc_coherent_scattering(None, 2, composition_H2O, 3, compton_scattering=3.5, method='ashcroft-langreth')
        fz_comp_test = core.calc_coherent_scattering([0], 2, composition_H2O, 3, method='faber-ziman')
        al_comp_test = core.calc_coherent_scattering([0], 2, composition_H2O, 3, method='ashcroft-langreth')
        self.assertEqual(fz_test, 2.5)
        self.assertEqual(al_test, 7.5)
        self.assertEqual(fz_comp_test, 6.0)
        self.assertEqual(al_comp_test, 18.0)

class TestNormaliseALFunc(unittest.TestCase, CustomAssertions):
    def test_normalise_al(self):
        sq_norm = core.normalise_achroft_langreth_func(0.5, sq=np.array([1, 2]))
        gr_norm = core.normalise_achroft_langreth_func(0.5, gr=np.array([1, 3]))
        rdf_norm = core.normalise_achroft_langreth_func(0.5, rdf=np.array([1, 3]), rho_0=(1/np.pi), r=np.array([0, 0.5]))
        expected_sq_norm = np.array([2.0, 4.0])
        expected_gr_norm = np.array([1.0, 5.0])
        expected_rdf_norm = np.array([2.0, 5.0])
        self.assertArrayEqual(sq_norm, expected_sq_norm)
        self.assertArrayEqual(gr_norm, expected_gr_norm)
        self.assertArrayEqual(rdf_norm, expected_rdf_norm)


class TestGetModFunc(unittest.TestCase, CustomAssertions):
    def test_get_lorch(self):
        test_Q = np.array([0, 0.02, 0.5, 1, 9.9, 10.0])
        expected_lorch = np.array([1.0,
                                   0.9999934202767204, 0.9958927352435614,
                                   0.983631643083466, 0.010099348633439868, 0])
        lorch = core.get_mod_func(test_Q, 'Lorch', None)
        self.assertFloatArrayEqual(lorch[1:], expected_lorch[1:])

        if not math.isclose(lorch[0], expected_lorch[0], rel_tol=1e-7, abs_tol=1e-16):
            warnings.warn(f'First values in lorch function arrays not equal | {lorch[0]} != {expected_lorch[0]}'
                          f' with relative tolerance of {1e-7:.3g} and absolute tolerance of {1e-16:.3g}')


    def test_get_cosine_window(self):
        test_Q = np.array([0, 0.02, 0.5, 1, 2.2, 4.2,
                           4.5, 6.5, 7.0, 9.9, 10.0])
        expected_window = np.array([1, 1, 1, 1, 1, 1,
                                    1.0,
                                    0.5 + 2**0.5/4,
                                    0.5,
                                    0.5 - 2**0.5/4,
                                    0])
        cosine_window = core.get_mod_func(test_Q, 'Cosine-window', 4.4)
        self.assertFloatArrayEqual(cosine_window, expected_window)

    def test_no_mod(self):
        Q = np.arange(0, 12, 0.02)
        self.assertTrue(core.get_mod_func(Q, None, None) == 1)


class TestCalcModelDIntraR(unittest.TestCase, CustomAssertions):
    def test_calc_model(self):
        # model = −4πrρ
        rho = 0.5
        r = np.array([0, 1.5, 2])
        test_model = core.calc_model_D_intra_r(rho, r)
        expected_model = np.array([0, -9.42477796076938, -12.566370614359172])
        test_model_unity = core.calc_model_D_intra_r(1/np.pi, [-1/4])
        expected_model_unity = np.array([1.0])
        self.assertFloatArrayEqual(test_model, expected_model)
        self.assertArrayEqual(test_model_unity, expected_model_unity)


class TestCalcChiSquared(unittest.TestCase, CustomAssertions):
    def test_calc_chi_squared(self):
        x = np.array([1, 3, 4])
        y = x
        chi_squared = core.calc_chi_squared(x, y)
        self.assertEqual(chi_squared, 21e6)

    def test_bad_method(self):
        x = np.array([1, 3, 4])
        y = x
        self.assertRaises(NotImplementedError, core.calc_chi_squared, x, y, method='quadrature')


class TestStopIteration(unittest.TestCase):
    def test_stop_iteration(self):
        self.assertTrue(core.stop_iteration(count=5, iter_limit=5))
        self.assertFalse(core.stop_iteration(count=5, iter_limit=6))
    def test_bad_method(self):
        self.assertRaises(NotImplementedError, core.stop_iteration, chi_squared=1.0, stop_condition='chi2_converge')


class TestIntegrateCoordSphere(unittest.TestCase, CustomAssertions):
    def test_integrate_coordination_sphere(self):
        test_r = np.arange(0, 3.5, 0.01)
        test_rdf = np.cos(test_r)
        r_0 = 0.0
        rp_max = np.pi/2.0
        r_max = np.pi/2.0
        r_min = np.pi
        N_a, N_b, N_c = core.integrate_coordination_sphere(test_r, test_rdf, r_0=r_0, rp_max=rp_max, r_max=r_max, r_min=r_min)
        self.assertFloatEqual(N_a, 2)
        self.assertFloatEqual(N_b, 2)
        self.assertFloatEqual(N_c, 0, atol=1e-10)


@unittest.skipIf(parse_version(python_version()) < parse_version('3.9'), 'Tests unsupported for python<3.9 - Suppport will be removed in a future version of LiquidDiffract')
class TestCalcStructureFactor(unittest.TestCase, CustomAssertions):
    def setUp(self):
        with resources.files('LiquidDiffract.scripts').joinpath('example_data.dat').open('r') as fp:
            _q_test, _I_test = np.loadtxt(fp, unpack=True, skiprows=0)
        self.q_test, self.I_test = data_utils.rebin_data(_q_test, _I_test, dx=0.02, extrapolate_mode='extrapolate')
        self.composition_Ga = {'Ga': (31,0,1)}
        self.composition_CaSiO3 = {'Ca': (20,0,1), 'Si': (14,0,1), 'O': (8,0,3)}
        self.S_inf_CaSiO3 = core.calc_S_inf(self.composition_CaSiO3, self.q_test)
        self.expected_SQ_Ga = np.load(os.path.join(data_path, 'SQ_initial_Ga.npy'))
        self.expected_SQ_CaSiO3_AL = np.load(os.path.join(data_path, 'SQ_initial_CaSiO3_AL.npy'))
        self.expected_SQ_CaSiO3_FZ = np.load(os.path.join(data_path, 'SQ_initial_CaSiO3_FZ.npy'))

    @unittest.skipIf(parse_version(sp_ver) < parse_version('1.11.0'), 'scipy.integrate.simpson behaviour changed in scipy=1.11.0 - Test data assumes scipy>1.11.0 and your version is earlier')
    def test_calc_sq_mono(self):
        SQ_Ga_AL = core.calc_structure_factor(self.q_test, self.I_test, self.composition_Ga, 0.048, method='ashcroft-langreth')
        SQ_Ga_FZ = core.calc_structure_factor(self.q_test, self.I_test, self.composition_Ga, 0.048, method='faber-ziman')
        self.assertFloatArrayEqual(SQ_Ga_AL, SQ_Ga_FZ)
        self.assertFloatArrayEqual(SQ_Ga_FZ, self.expected_SQ_Ga)

    @unittest.skipIf(parse_version(sp_ver) < parse_version('1.11.0'), 'scipy.integrate.simpson behaviour changed in scipy=1.11.0 - Test data assumes scipy>1.11.0 and your version is earlier')
    def test_calc_sq_poly(self):
        SQ_CaSiO3_AL = core.calc_structure_factor(self.q_test, self.I_test, self.composition_CaSiO3, 0.048, method='ashcroft-langreth')
        SQ_CaSiO3_AL_rescaled = core.normalise_achroft_langreth_func(SQ_CaSiO3_AL, self.S_inf_CaSiO3)
        SQ_CaSiO3_FZ = core.calc_structure_factor(self.q_test, self.I_test, self.composition_CaSiO3, 0.048, method='faber-ziman')
        self.assertFloatArrayEqual(SQ_CaSiO3_AL, self.expected_SQ_CaSiO3_AL)
        self.assertFloatArrayEqual(SQ_CaSiO3_FZ, self.expected_SQ_CaSiO3_FZ)
        self.assertNotFloatArrayEqual(SQ_CaSiO3_AL_rescaled, SQ_CaSiO3_FZ)

    def test_bad_method(self):
        self.assertRaises(ValueError, core.calc_structure_factor, self.q_test, self.I_test, self.composition_Ga, 0.048, method='BAD')

    def test_prescale_IQ(self):
        SQ_CaSiO3_FZ = core.calc_structure_factor(self.q_test, self.I_test, self.composition_CaSiO3, 0.048, method='faber-ziman')
        SQ_CaSiO3_FZ_a = core.calc_structure_factor(self.q_test, self.I_test*10.0, self.composition_CaSiO3, 0.048, method='faber-ziman')
        SQ_CaSiO3_FZ_b = core.calc_structure_factor(self.q_test, self.I_test*1e9, self.composition_CaSiO3, 0.048, method='faber-ziman')
        SQ_CaSiO3_FZ_c = core.calc_structure_factor(self.q_test, self.I_test*1e-9, self.composition_CaSiO3, 0.048, method='faber-ziman')
        self.assertFloatArrayEqual(SQ_CaSiO3_FZ_a, SQ_CaSiO3_FZ)
        self.assertFloatArrayEqual(SQ_CaSiO3_FZ_b, SQ_CaSiO3_FZ)
        self.assertFloatArrayEqual(SQ_CaSiO3_FZ_c, SQ_CaSiO3_FZ)

        SQ_CaSiO3_AL = core.calc_structure_factor(self.q_test, self.I_test, self.composition_CaSiO3, 0.048, method='ashcroft-langreth')
        SQ_CaSiO3_AL_a = core.calc_structure_factor(self.q_test, self.I_test*10.0, self.composition_CaSiO3, 0.048, method='ashcroft-langreth')
        SQ_CaSiO3_AL_b = core.calc_structure_factor(self.q_test, self.I_test*1e9, self.composition_CaSiO3, 0.048, method='ashcroft-langreth')
        SQ_CaSiO3_AL_c = core.calc_structure_factor(self.q_test, self.I_test*1e-9, self.composition_CaSiO3, 0.048, method='ashcroft-langreth')
        self.assertFloatArrayEqual(SQ_CaSiO3_AL_a, SQ_CaSiO3_AL)
        self.assertFloatArrayEqual(SQ_CaSiO3_AL_b, SQ_CaSiO3_AL)
        self.assertFloatArrayEqual(SQ_CaSiO3_AL_c, SQ_CaSiO3_AL)


@unittest.skipIf(parse_version(python_version()) < parse_version('3.9'), 'Tests unsupported for python<3.9 - Suppport will be removed in a future version of LiquidDiffract')
class TestCalcCorrelationFunc(unittest.TestCase, CustomAssertions):

    def setUp(self):
        # Can combine with fixtures from other tests
        with resources.files('LiquidDiffract.scripts').joinpath('example_data.dat').open('r') as fp:
            _q_test, _I_test = np.loadtxt(fp, unpack=True, skiprows=0)
        self.q_test, self.I_test = data_utils.rebin_data(_q_test, _I_test, dx=0.02,extrapolate_mode='extrapolate')
        self.composition_Ga = {'Ga': (31,0,1)}
        self.rho = 0.048
        self.intf_func_Ga = core.calc_structure_factor(self.q_test, self.I_test, self.composition_Ga, self.rho, method='faber-ziman') - core.calc_S_inf(self.composition_Ga, self.q_test)
        self.expected_r_Ga, self.expected_Dr_Ga = np.load(os.path.join(data_path, 'Dr_initial_Ga.npy'))

    def test_bad_method(self):
        self.assertRaises(ValueError, core.calc_correlation_func, self.q_test, self.intf_func_Ga, self.rho, function='BAD')

    @unittest.skipIf(parse_version(sp_ver) < parse_version('1.11.0'), 'scipy.integrate.simpson behaviour changed in scipy=1.11.0 - Test data assumes scipy>1.11.0 and your version is earlier')
    def test_calc_Dr(self):
        # Test calculation of D(r) and r-values
        r_Ga, Dr_Ga = core.calc_correlation_func(self.q_test, self.intf_func_Ga, self.rho, mod_func='Cosine-window', window_start=8.0)
        self.assertFloatArrayEqual(Dr_Ga, self.expected_Dr_Ga)
        self.assertFloatArrayEqual(r_Ga, self.expected_r_Ga)

    def test_calc_correlation_funcs(self):
        # Test calculation of g(r), D(r), and RDF(r)
        r_Ga_from_Dr, Dr_Ga = core.calc_correlation_func(self.q_test, self.intf_func_Ga, self.rho, function='density_func')
        r_Ga_from_gr, gr_Ga = core.calc_correlation_func(self.q_test, self.intf_func_Ga, self.rho, function='pair_dist_func')
        r_Ga_from_RDF, RDF_Ga = core.calc_correlation_func(self.q_test, self.intf_func_Ga, self.rho, function='radial_dist_func')
        # Check r-values calculated in the same way
        self.assertFloatArrayEqual(r_Ga_from_Dr, r_Ga_from_gr)
        self.assertFloatArrayEqual(r_Ga_from_gr, r_Ga_from_RDF)
        with np.errstate(divide='ignore', invalid='ignore'):
            # Test function returned g(r) vs recalculating from D(r)
            test_gr_Ga = (Dr_Ga / (4*np.pi*r_Ga_from_Dr*self.rho)) + 1
        self.assertFloatArrayEqual(test_gr_Ga, gr_Ga)
        # Test function returned RDF vs recalculating from g(r)
        test_RDF_Ga = 4 * np.pi * r_Ga_from_gr**2 * self.rho * gr_Ga
        # g(r)[0] is nan, due to division by r (r[0] = 0)
        self.assertFloatArrayEqual(test_RDF_Ga[1:], RDF_Ga[1:])
        if not math.isclose(test_RDF_Ga[0], RDF_Ga[0], rel_tol=1e-7, abs_tol=1e-16):
            warnings.warn(f'First values in RDF arrays not equal | {test_RDF_Ga[0]} != {RDF_Ga[0]}'
                          f' with relative tolerance of {1e-7:.3g} and absolute tolerance of {1e-16:.3g}')

    def test_calc_correlation_FT(self):
        # Sinc function should give peak centred at pi
        x = np.arange(0, 50, 0.5)
        y = np.sinc(x)
        r, Dr = core.calc_correlation_func(x, y, 1.0)
        self.assertEqual(r[np.argmax(Dr)].round(2), 3.14)

    def test_calc_D_r_iteration_term(self):
        # Test backward transform case
        r_Ga, Dr_Ga = core.calc_correlation_func(self.q_test, self.intf_func_Ga, self.rho)
        with np.errstate(divide='ignore', invalid='ignore'):
            test_intf_func_Ga = core.calc_D_r_iteration_term(Dr_Ga*np.pi/2, dq=0.02)[:len(self.q_test)]/self.q_test
        # [0] value not recovered as * 0, then later / 0
        self.assertFloatArrayEqual(test_intf_func_Ga[1:], self.intf_func_Ga[1:])
        if not math.isclose(test_intf_func_Ga[0], self.intf_func_Ga[0], rel_tol=1e-7, abs_tol=1e-16):
            warnings.warn(f'First values in int_func arrays not equal | {test_intf_func_Ga[0]} != {self.intf_func_Ga[0]}'
                          f' with relative tolerance of {1e-7:.3g} and absolute tolerance of {1e-16:.3g}')


@unittest.skipIf(parse_version(python_version()) < parse_version('3.9'), 'Tests unsupported for python<3.9 - Suppport will be removed in a future version of LiquidDiffract')
class TestCalcImprIntFunc(unittest.TestCase, CustomAssertions):

    def setUp(self):
        # Need to combine with fixtures from other tests
        with resources.files('LiquidDiffract.scripts').joinpath('example_data.dat').open('r') as fp:
            _q_test, _I_test = np.loadtxt(fp, unpack=True, skiprows=0)
        self.q_test, self.I_test = data_utils.rebin_data(_q_test, _I_test, dx=0.02, extrapolate_mode='extrapolate')
        self.composition_Ga = {'Ga': (31,0,1)}
        self.rho = 0.048
        self.intf_func_Ga = core.calc_structure_factor(self.q_test, self.I_test, self.composition_Ga, self.rho, method='faber-ziman') - core.calc_S_inf(self.composition_Ga, self.q_test)
        self.expected_refined_intf_func_Ga = np.load(os.path.join(data_path, 'iQ_refined_Ga.npy'))

    @unittest.skipIf(parse_version(sp_ver) < parse_version('1.11.0'), 'scipy.integrate.simpson behaviour changed in scipy=1.11.0 - Test data assumes scipy>1.11.0 and your version is earlier')
    def test_calc_impr_interference_func(self):
        r_min = 2.3
        iter_limit = 5
        method = 'faber-ziman'
        mod_func = None
        window_start = None
        fft_N = 12
        impr_intf_func, chi_sq = core.calc_impr_interference_func(self.q_test, self.intf_func_Ga,
                                                          self.composition_Ga, self.rho,
                                                          r_min, iter_limit,
                                                          method, mod_func, window_start, fft_N)
        self.assertFloatArrayEqual(impr_intf_func, self.expected_refined_intf_func_Ga)


if __name__ == "__main__":
    unittest.main()