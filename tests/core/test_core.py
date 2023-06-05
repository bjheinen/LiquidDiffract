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


if __name__ == "__main__":
    unittest.main()