# -*- coding: utf-8 -*-
"""
Custom test classes for unit tests
"""
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2023, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

import math
import os.path
import numpy as np

# Get test data path
unittest_path = os.path.dirname(__file__)
data_path = os.path.join(unittest_path, './data')
__unittest = True

class CustomAssertions:
    # Compare floats
    def assertFloatEqual(self, a, b, rtol=1e-7, atol=1e-16, err_msg=''):       
        std_err_msg = (f'{a} != {b} with relative tolerance of {rtol:.3g}'
                       f' and absolute tolerance of {atol:.3g}')
        # equivalent to abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
        assert math.isclose(a, b, rel_tol=rtol, abs_tol=atol), f'{std_err_msg} : {err_msg}'
    # Compare numpy arrays with np.testing
    def assertArrayEqual(self, a, b, err_msg='', verbose=True, strict=True):
        np.testing.assert_array_equal(a, b, err_msg=err_msg, verbose=verbose, strict=strict)
    def assertFloatArrayEqual(self, a, b, rtol=1e-7, atol=1e-16, err_msg='', equal_nan=True, verbose=True):
        # np.allclose calculates abs(a - b) <= (atol + rtol * abs(b)), which is not great
        np.testing.assert_allclose(a, b, rtol=rtol, atol=atol, err_msg=err_msg, equal_nan=equal_nan, verbose=verbose)
    def assertNotFloatArrayEqual(self, a, b, rtol=1e-7, atol=1e-16, err_msg='', equal_nan=True, verbose=True):
        np.testing.assert_raises(AssertionError,
                                 np.testing.assert_allclose,
                                 a, b, rtol=rtol, atol=atol, err_msg=err_msg, equal_nan=equal_nan, verbose=verbose)