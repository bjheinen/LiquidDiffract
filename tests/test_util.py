# -*- coding: utf-8 -*-
"""
Unit tests for custom test classes in tests/util
"""
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2023, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

import unittest
import os.path
import numpy as np
from util import CustomAssertions
from util import data_path

class TestCustomAssertions(unittest.TestCase, CustomAssertions):
    def test_assert_float_equal(self):
        self.assertFloatEqual(1000.0, 1001, rtol=1e-3, atol=1e-16)
        self.assertFloatEqual(1000.0, 1000.0+1e-5)
        self.assertFloatEqual(1e-5, 1e-5+1e-12)
        self.assertFloatEqual(0, 1e-17)

    def test_assert_array_equal(self):
        a = np.array([1,1e-12,2,1e-5,1e6,0])
        b = np.array([1,1e-12,2,1e-5,1e6,0])
        self.assertArrayEqual(a, b)

    def test_assert_float_array_equal(self):
        a = np.array([1,1e-12,2,1e-5,1e6,0])
        b = a + 1e-20
        self.assertFloatArrayEqual(a,b)

    def test_assert_not_float_array_equal(self):
        a = np.array([1,1e-12,2,1e-5,1e6,0])
        b = a + 1e-3
        self.assertNotFloatArrayEqual(a, b)

class TestDataPath(unittest.TestCase):
    def test_data_path(self):
        fname = os.path.join(data_path, 'atomic_ff_Ga_0-12.npy')
        assert os.path.isfile(fname)

if __name__ == "__main__":
    unittest.main()