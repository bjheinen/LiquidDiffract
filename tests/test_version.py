# -*- coding: utf-8 -*-
__author__ = "Benedict J. Heinen"
__copyright__ = "Copyright 2023, Benedict J. Heinen"
__email__ = "benedict.heinen@gmail.com"

import unittest
import warnings
from packaging.version import Version

class TestVersionFile(unittest.TestCase):
    def test_appname(self):
        from LiquidDiffract.version import __appname__
        self.assertEqual(__appname__, 'LiquidDiffract')
    def test_version(self):
        from LiquidDiffract.version import __version__
        # Version numbers on dev branch like '1.2.7-dev'
        # Check if version number still -dev and warn if true
        self.assertGreater(len(__version__), 4)
        try:
            self.assertTrue(Version(__version__))
        except InvalidVersion:
            self.fail("Version number invalid!")
        if 'dev' in __version__ or not __version__[-1].isdigit():
            warnings.warn(f'Version number not in release format| v: {__version__}')

if __name__ == "__main__":
    unittest.main()