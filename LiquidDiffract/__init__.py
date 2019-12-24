'''
LiquidDiffract
==============

LiquidDiffract is a graphical data analysis application for processing
experimental total X-ray diffraction data of liquids and disordered solids.

LiquidDiffract can perform background subtraction,fourier transformation,
numerical optimisation/normalisation of the quantitative structure factor,
and estimation of sample density.

The sub-package LiquidDiffract.core provides useful functions for common
numerical operations on total diffraction data.

Core functions available from LiquidDiffract.core.core

    >>> import LiquidDiffract.core.core
    >>> help(LiquidDiffract.core.core)

Additional utilities from LiquidDiffract.core.data_utils

    >>> import LiquidDiffract.core.data_utils
    >>> help(LiquidDiffract.core.data_utils)


Please see the documentation on the LiquidDiffract homepage for more details.
<https://www.github.com/bjheinen/LiquidDiffract>

Copyright (c) 2019 Benedict J. Heinen

Licensed under the GNU General Public License (GPL), version 3 or later.

'''
from LiquidDiffract.version import __version__
