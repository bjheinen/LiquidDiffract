#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
LiquidDiffract application loader script
<https://github.com/bjheinen/LiquidDiffract>
'''
__author__ = 'Benedict J. Heinen'
__copyright__ = 'Copyright 2018-2024, Benedict J. Heinen'
__license__ = 'Gnu GPL v3'
__email__ = 'benedict.heinen@gmail.com'

import sys
import os
from importlib.util import find_spec
# Check for QT_API environment variable, default to PyQt6 if available
if os.environ.get('QT_API') is None:
    if find_spec('PyQt6'):
        os.environ['QT_API'] = 'PyQt6'
from qtpy.QtWidgets import QApplication
# Get Qt API info from qtpy
from qtpy import QT_VERSION, API_NAME
import LiquidDiffract.gui.main_widget
# Get version number from version.py
from LiquidDiffract.version import __version__, __appname__

def main():
    '''Launch the GUI'''
    app = QApplication(sys.argv)
    screen_size = app.primaryScreen().size()
    print(f'{__appname__} v{__version__}')
    print(f'Using Qt v{QT_VERSION} with {API_NAME}')
    print(f'{__copyright__}\n')
    _ = LiquidDiffract.gui.main_widget.App(screen_size)
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
