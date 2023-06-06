#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for LiquidDiffract.core.core
"""
from __future__ import absolute_import
import unittest
from importlib import resources
import os.path
import numpy as np
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
        with self.assertRaises(KeyError):
            core.calc_mol_mass(composition)

    def test_data_coherence(self):
        with resources.open_binary('LiquidDiffract.resources', 'mass_data.npy') as fp:
            mass_dict = np.load(fp, allow_pickle=True).item()
        with resources.open_binary('LiquidDiffract.resources', 'pt_data.npy') as fp:
            element_dict = np.load(fp, allow_pickle=True).item()
        self.assertEqual(len(mass_dict), len(element_dict))
        for element in element_dict:
            self.assertIn(element, mass_dict)


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


class TestLoadData(unittest.TestCase):
    def test_load_ff_data(self):
        self.assertEqual(core.load_ff_data().shape, (210, 11))
    def test_load_compton_data(self):        
        with resources.open_binary('LiquidDiffract.resources', 'pt_data.npy') as fp:
            element_dict = np.load(fp, allow_pickle=True).item()
        for element in element_dict:
            self.assertTrue(core.load_compton_data(element))


class TestCalcAtomicFF(unittest.TestCase, CustomAssertions):
    def test_calc_atomic_ff(self):
        expected_atomic_ff_Ga = np.load(os.path.join(data_path, 'atomic_ff_Ga_0-12.npy'))
        expected_atomic_ff_Si = np.load(os.path.join(data_path, 'atomic_ff_Si_0-12.npy'))
        element_Ga = (31,0,1) # Gallium
        element_Si = (14,0,1) # Silicon
        Q = np.arange(0, 12, 0.02)
        atomic_ff_Ga = core.calc_atomic_ff(element_Ga,Q)
        atomic_ff_Si = core.calc_atomic_ff(element_Si,Q)
        self.assertArrayEqual(atomic_ff_Ga, expected_atomic_ff_Ga)
        self.assertArrayEqual(atomic_ff_Si, expected_atomic_ff_Si)


class TestCalcEffectiveFF(unittest.TestCase, CustomAssertions):
    def test_calc_effective_ff(self):
        Q = np.arange(0, 12, 0.02)
        composition_Ga = {'Ga': (31,0,1)}
        composition_CaSiO3 = {'Ca': (20,0,1), 'Si': (14,0,1), 'O': (8,0,3)}
        composition_GaSn = {'Ga': [31, 0, 915], 'Sn': [50, 0, 85]}
        expected_aff_Ga = np.load(os.path.join(data_path, 'atomic_ff_Ga_0-12.npy'))
        expected_eff_CaSiO3 = np.load(os.path.join(data_path, 'effective_ff_CaSiO3_0-12.npy'))

        eff_Ga, aff_Ga = core.calc_effective_ff(composition_Ga, Q)
        self.assertTrue(len(aff_Ga) == 1)
        self.assertTrue(len(aff_Ga[0]) == 600)
        self.assertFloatArrayEqual(aff_Ga[0], expected_aff_Ga)
        self.assertFloatArrayEqual(aff_Ga/eff_Ga/31.0, 1.0)

        eff_CaSiO3, _ = core.calc_effective_ff(composition_CaSiO3, Q)
        self.assertFloatArrayEqual(eff_CaSiO3, expected_eff_CaSiO3)

        self.assertIsInstance(core.calc_effective_ff(composition_GaSn, Q), tuple)


class TestCalcAverageScattering(unittest.TestCase, CustomAssertions):
    def test_calc_average_scattering(self):
        Q = np.arange(0, 12, 0.02)
        composition_Ga = {'Ga': (31,0,1)}
        composition_CaSiO3 = {'Ca': (20,0,1), 'Si': (14,0,1), 'O': (8,0,3)}
        expected_avg_scattering_CaSiO3 = np.load(os.path.join(data_path, 'average_scattering_functions_CaSiO3_0-12.npy'))
        average_scattering_Ga = core.calc_average_scattering(composition_Ga, Q)
        average_scattering_CaSiO3 = core.calc_average_scattering(composition_CaSiO3, Q)
        # <f2> == <f>2 for monatomic case
        self.assertFloatArrayEqual(average_scattering_Ga[0], average_scattering_Ga[1])
        self.assertEqual(len(average_scattering_CaSiO3), 2)
        self.assertFloatArrayEqual(average_scattering_CaSiO3[0], expected_avg_scattering_CaSiO3[0])
        self.assertFloatArrayEqual(average_scattering_CaSiO3[1], expected_avg_scattering_CaSiO3[1])


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
        expected_K_CaSiO3 = np.array([22.28649729, 15.07790776,  6.87853165,  6.87853165,  6.87853165])
        K_Ga = core.calc_K_p(composition_Ga, Q)
        K_GaGa = core.calc_K_p(composition_GaGa, Q)
        K_CaSiO3 = core.calc_K_p(composition_CaSiO3, Q)
        self.assertTrue(len(K_Ga), 1)
        self.assertTrue(len(K_GaGa), 2)
        self.assertTrue(len(K_CaSiO3), 5)
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
        self.assertTrue(test_1 == test_2 == 1.0)


class TestCalcAlpha(unittest.TestCase, CustomAssertions):
    def test_calc_alpha(self):
        with resources.files('LiquidDiffract.scripts').joinpath('example_data.dat').open('r') as fp:
            q_test, I_test = np.loadtxt(fp, unpack=True, skiprows=0)
        q_test, I_test = data_utils.rebin_data(q_test, I_test, dx=0.02)
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
        self.assertFloatEqual(alpha_Ga_fz, 12.38770462391852)
        self.assertFloatEqual(alpha_CaSiO3_al, 10.544150126760833)
        self.assertFloatEqual(alpha_CaSiO3_fz, 2.201685124321958)


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



if __name__ == "__main__":
    unittest.main()