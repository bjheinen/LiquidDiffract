#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest

class TestVersionFile(unittest.TestCase):
    def test_appname(self):
        from LiquidDiffract.version import __appname__
        self.assertEqual(__appname__, 'LiquidDiffract')
    def test_version(self):
        from LiquidDiffract.version import __version__
        # Version numbers on dev branch like '1.2.7-dev'
        # Check if version number still -dev
        self.assertNotIn('dev', __version__)
        self.assertTrue(__version__[-1].isdigit())

if __name__ == "__main__":
    unittest.main()