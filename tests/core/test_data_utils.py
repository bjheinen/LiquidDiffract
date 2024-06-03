#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for LiquidDiffract.core.core
"""
from __future__ import absolute_import
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2024, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

import unittest
import os.path
import math
import numpy as np

# import core module to test
import LiquidDiffract.core.data_utils as data_utils

# import data path and custom assertions from tests/util
from util import data_path
from util import CustomAssertions


class TestInterp(unittest.TestCase, CustomAssertions):

    def setUp(self):
        self.x = np.arange(0, 7.0, 0.0001)
        self.y = np.cos(self.x)

    def test_interp_data(self):
        x2 = np.array([0, np.pi/2.0, np.pi, 2*np.pi])
        expected_interp = np.array([1.0, 0.0, -1.0, 1.0])
        test_interp = data_utils.interp_data(self.x, self.y, x2)
        self.assertFloatArrayEqual(test_interp, expected_interp)

    def test_interp_nan(self):
        test_y = np.copy(self.y)
        test_y[0] = np.nan
        self.assertFloatArrayEqual(data_utils.interp_nan(test_y), self.y)


class TestConvertTwoThetaQSpace(unittest.TestCase, CustomAssertions):

    def setUp(self):
        self.wavelength_a = 1.0
        self.wavelength_b = 0.1722
        self.test_two_theta = np.array([0.0, 180.0, 9.128558416134155])
        self.test_q_space = np.array([0.0, 4*np.pi, 1.0])

    def test_convert_two_theta(self):
        self.assertFloatArrayEqual(data_utils.convert_two_theta(self.test_two_theta, self.wavelength_a), self.test_q_space, atol=1.e-14)
        self.assertFloatArrayEqual(data_utils.convert_two_theta(self.test_two_theta, self.wavelength_b), self.test_q_space/self.wavelength_b, atol=1.e-14)

    def test_convert_q_space(self):
        self.assertFloatArrayEqual(np.degrees(data_utils.convert_q_space(self.test_q_space, self.wavelength_a)), self.test_two_theta, atol=1.e-14)
        self.assertFloatArrayEqual(data_utils.convert_q_space(self.test_q_space/self.wavelength_b, self.wavelength_b), np.radians(self.test_two_theta), atol=1.e-14)

class TestSmoothData(unittest.TestCase, CustomAssertions):

    def test_smooth_data(self):
        x = np.arange(0, 100, 0.01)
        self.assertFloatArrayEqual(data_utils.smooth_data(x), x)
        self.assertTrue(type(data_utils.smooth_data(x)), np.ndarray)
        self.assertRaises(NotImplementedError, data_utils.smooth_data, x, method='any')

class TestZeroNorm(unittest.TestCase, CustomAssertions):

    def test_zero_norm(self):
        y = np.array([1.2, 2.0, 5.0, 0.2])
        self.assertFloatArrayEqual(data_utils.zero_norm(y), np.array([0.0, 0.8, 3.8, -1.0]))
        self.assertFloatArrayEqual(data_utils.zero_norm(y, shift=2.0), np.array([-0.8, 0.0, 3.0, -1.8]))


class TestBkgScalingResidual(unittest.TestCase, CustomAssertions):

    def test_bkg_scaling_residual(self):
        x = np.random.uniform(-5.0, 15.0, 5)
        self.assertFloatEqual(data_utils.bkg_scaling_residual(5, *(x, x*0.2)), 0.0, atol=1e-15)
        self.assertFloatEqual(data_utils.bkg_scaling_residual(1, x, x+1), 1.0, atol=1e-15)
        self.assertFloatEqual(data_utils.bkg_scaling_residual(0.2, x, x*5), 0.0, atol=1e-15)


class TestFindIntegrationLimits(unittest.TestCase, CustomAssertions):

    def test_find_integration_limits(self):
        # Load test RDF data
        rdf_data_a, rdf_data_b, rdf_data_c = np.load(os.path.join(data_path, 'rdf_peak_test_data.npy'), allow_pickle=True)
        # State expected limit positions
        expected_limits_a = (2.239933089082285, 2.9162691080482155,
                             2.9444851167454846, 3.4521306936532428)
        expected_limits_b = (2.1482738727778155, 2.906024136949478,
                             2.9433519729572653, 3.5421030348918663)
        expected_limits_c = (1.4413920965789064, 1.6158820019375424,
                             1.621174657388145, 1.7887945926039175)

        # Test regular case where peak_idx_prom == peak_idx_first
        limits_a = data_utils.find_integration_limits(*rdf_data_a.T)
        # Test kwargs, test higher peak search etc.
        limits_a_kw = data_utils.find_integration_limits(*rdf_data_a.T, rho=0.05, peak_search_limit=10)
        limits_a_psl = data_utils.find_integration_limits(*rdf_data_a.T, rho=0.05, peak_search_limit=50)

        # Test case where no next peak, for rho and no rho
        limits_a_np = data_utils.find_integration_limits(*rdf_data_a.T, peak_search_limit=3.5)
        limits_a_np_rho = data_utils.find_integration_limits(*rdf_data_a.T, rho=0.05, peak_search_limit=3.5)

        # Test case where peak_idx_prom > peak_idx_first (oscillations at RDF ~ 0 need to be discounted)
        limits_b = data_utils.find_integration_limits(*rdf_data_b.T)

        # Test case where peak_idx_prom < peak_idx_first (right base at RDF < 0)
        limits_c = data_utils.find_integration_limits(*rdf_data_c.T)

        self.assertFloatArrayEqual(limits_a, expected_limits_a)
        self.assertFloatArrayEqual(limits_a_kw, expected_limits_a)
        self.assertFloatArrayEqual(limits_a_kw, limits_a_psl)

        self.assertFloatArrayEqual(limits_a_np[:-1], expected_limits_a[:-1])
        self.assertFloatArrayEqual(limits_a_np_rho[:-1], expected_limits_a[:-1])
        self.assertFloatArrayEqual(limits_a_np, expected_limits_a, atol=1.e-6)
        self.assertFloatArrayEqual(limits_a_np_rho, expected_limits_a, atol=1.e-6)

        self.assertFloatArrayEqual(limits_b, expected_limits_b)
        self.assertFloatArrayEqual(limits_c, expected_limits_c)


class TestRebinData(unittest.TestCase, CustomAssertions):

    def test_rebin_data(self):
        x_small_dx = np.arange(0, 12.0005, 0.001)
        y_small_dx = np.cos(x_small_dx)
        x_large_dx = np.arange(0, 12.01, 0.02)
        y_large_dx = np.cos(x_large_dx)
        rebin_x_small_dx, rebin_y_small_dx = data_utils.rebin_data(x_small_dx, y_small_dx, dx=0.02)
        rebin_x_large_dx, rebin_y_large_dx = data_utils.rebin_data(x_large_dx, y_large_dx, dx=0.02)
        # Check dx respected
        self.assertFloatArrayEqual(np.diff(rebin_x_small_dx), 0.02)
        # Check rebinning binned array returns *same*
        self.assertFloatArrayEqual(rebin_y_large_dx, y_large_dx[:-1])
        # Check interpolation in rebinning accurate
        self.assertFloatArrayEqual(rebin_y_small_dx, y_large_dx[:-1])
        # Check auto limits are 0 and x[-1]
        self.assertEqual(rebin_x_small_dx[0], 0.0)
        self.assertEqual(rebin_x_small_dx[-1], 11.98)
        self.assertEqual(rebin_x_large_dx[-1], 11.98)

    def test_rebin_nan(self):
        # Test nan values in y have little effect
        x = np.arange(0, 12.0, 0.001)
        y = np.cos(x)
        # Inject single np.nan
        y[500] = np.nan
        expected_rebin_y = np.cos(np.arange(0, 12, 0.02))
        _, rebin_y = data_utils.rebin_data(x, y, dx=0.02)
        # Check interpolation in rebinning accurate
        # y[500] --> rebin_y[25] (0.02/0.001 = 20)
        # Check arrays equal skipping interpolated value
        self.assertFloatArrayEqual(np.delete(rebin_y, 25), np.delete(expected_rebin_y, 25))
        self.assertFloatEqual(rebin_y[25], expected_rebin_y[25], rtol=1e-6)

    def test_rebin_extrapolate_fill(self):
        x = np.arange(0.34, 12.0, 0.02)
        y = np.cos(x)
        rebin_x, rebin_y = data_utils.rebin_data(x, y, dx=0.02)
        # Check rebin has extrapolated/filled down to 0
        self.assertEqual(rebin_x[0], 0.0)
        # Check y[0] used as fill value at x < x[0]
        fill_values = set(rebin_y[np.where(rebin_x < x[0])])
        self.assertEqual(len(fill_values), 1)
        self.assertEqual(fill_values.pop(), y[0])
