#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'Benedict J Heinen'
__copyright__ = 'Copyright 2018, Benedict J Heinen'
__license__ = 'Gnu GPL v3'
__email__ = 'benedict.heinen@gmail.com'


import sys
from PyQt5.QtWidgets import QApplication
import gui.main_widget


def main():
    app = QApplication(sys.argv)
    screen_size = app.primaryScreen().size()
    ex = gui.main_widget.App(screen_size)    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
