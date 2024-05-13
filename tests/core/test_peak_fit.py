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
from scipy.integrate import simpson

#import core modules to test
import LiquidDiffract.core.peak_fit as peak_fit

# import data path and custom assertions from tests/util
from util import data_path
from util import CustomAssertions


class TestCachedSqrt(unittest.TestCase):
    def test_cached_sqrt(self):
        self.assertEqual(peak_fit.cached_sqrt('2'), 1.4142135623730951)
        self.assertEqual(peak_fit.cached_sqrt('2pi'), 2.5066282746310002)
        self.assertEqual(peak_fit.cached_sqrt('256'), 16.0)


class TestSkewGauss(unittest.TestCase, CustomAssertions):
    def test_skew_gauss(self):
        # skew_gauss(r, N_ab, r_ab, sigma_ab, xi_ab, W_ab, c_b)
        test_r = np.arange(0, 20, 0.002)
        # At r = r_ab the second and third terms cancel
        self.assertEqual(peak_fit.skew_gauss(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1/math.sqrt(2*math.pi))
        self.assertEqual(peak_fit.skew_gauss(1.0, 1.0, 1.0, 1.0, 1.0, math.sqrt(math.pi * 2), 1.0), 1.0)
        self.assertEqual(peak_fit.skew_gauss(1.0, 1.0, 1.0, 1.0, 1.95e6, math.sqrt(math.pi * 2), 1.0), 1.0)
        self.assertEqual(peak_fit.skew_gauss(1.0, 4.23, 1.0, 1.0, 1.0, math.sqrt(math.pi * 2), 1.0), 4.23)
        self.assertEqual(peak_fit.skew_gauss(test_r, 4.23, 1.0, 1.0, 1.0, math.sqrt(math.pi * 2), 1.0)[np.where(test_r == 1.0)], 4.23)
        test_r = np.arange(0, 3.5, 0.0102)
        expected_gauss_Si_O = np.load(os.path.join(data_path, 'gaussian_Si_O.npy'))
        gauss_Si_O = peak_fit.skew_gauss(test_r, 4.0, 1.624, 0.11, 0.0, 0.19974189041782536, 0.6)
        small_skew_gauss_Si_O = peak_fit.skew_gauss(test_r, 4.0, 1.624, 0.11, 0.001, 0.19974189041782536, 0.6)
        big_skew_gauss_Si_O = peak_fit.skew_gauss(test_r, 4.0, 1.624, 0.11, 0.001, 0.19974189041782536, 0.6)
        self.assertArrayEqual(gauss_Si_O, expected_gauss_Si_O)
        self.assertNotFloatArrayEqual(small_skew_gauss_Si_O, expected_gauss_Si_O)
        # Skew should not change area
        self.assertFloatEqual(simpson(gauss_Si_O, x=test_r), simpson(small_skew_gauss_Si_O, x=test_r))
        self.assertFloatEqual(simpson(small_skew_gauss_Si_O, x=test_r), simpson(big_skew_gauss_Si_O, x=test_r))
        # N_ab = c_b/w_ab * Area
        self.assertFloatEqual(0.6 * simpson(small_skew_gauss_Si_O, x=test_r)/0.19974189041782536, 4.0)

