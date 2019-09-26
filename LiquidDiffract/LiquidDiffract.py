#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
LiquidDiffract application loader script
<https://github.com/bjheinen/LiquidDiffract>
'''
__author__ = 'Benedict J Heinen'
__copyright__ = 'Copyright 2018-2019, Benedict J Heinen'
__license__ = 'Gnu GPL v3'
__email__ = 'benedict.heinen@gmail.com'
# Get version number from version.py
from LiquidDiffract.version import __version__

import sys
from PyQt5.QtWidgets import QApplication
import LiquidDiffract.gui.main_widget


def main():
    app = QApplication(sys.argv)
    screen_size = app.primaryScreen().size()
    ex = LiquidDiffract.gui.main_widget.App(screen_size)    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()