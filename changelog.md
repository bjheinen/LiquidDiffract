# Changelog
---------

## v1.1.12 (May 17, 2023)

### Bug Fixes

- Fixed a breaking bug due to the deprecation of np.float/np.int/np.str etc. as of numpy v1.20.0 in favour of built-in types.
- Fixed some possible buggy behaviour in user inputs

### Minor Changes

- An option has been added to plot log(_I(Q)_) in the background subtraction tab to make small differences between the measured and background _I(Q)_ more obvious.
- Added a prompt to select a file header when a file load error occurs.


## v1.1.11 (March 2, 2023)

### Patch release to fix errors in tabulated Compton scattering data

* The data in LiquidDiffract/LiquidDiffract/resources/hubbel_compton/ was incorrect for several elements due to a filename issue. Some filenames dropped a character and overwrote data for elements with single-character symbols. The affected elements were  C, F, N, P, and S, and data for Ac, Am, Cm, Fm, Np, Pm, Sc, Tc, and Tm were missing. Using Hg resulted in a (crashing) error due to a mislabeling as hg. Data has been corrected for all elements and now matches the human-readable tabulations in LiquidDiffract/LiquidDiffract/resources/human_readable/hubbel_compton/


## v1.1.10 (October 18, 2022)

### Patch release to fix minor (but crashing) bug in v1.1.9 when using Python v>=3.10

* LiquidDiffract/gui/optim_ui.py line 705: setResizeMode --> setSectionResizeMode
  The QHeaderView attribute setResizeMode is deprecated in PyQt5 and causes an AttributeError in some version combinations of PyQt5 and Python. setSectionResizeMode replaces it and is compatible with all versions.


## v1.1.9 (July 14, 2022)

### Major Changes

- **Background refinement**
  The background scaling factor (_b_) can now be refined alongside the density in the _Refinement Tab_
  Background refinement is optional and can be done independently from or simultaneously with density refinement.
  The background scaling factor set in the _Background Subtraction Tab_ is used as the initial estimate by the solver.
  The figure of merit is redefined to take into account the background scaling - _&Chi;_^2(&rho;;b)

- **LiquidDIffraction publication**
  There is now a paper describing LiquidDiffract:
  [Heinen, B. J., & Drewitt, J. W. (2022). LiquidDiffract: Software for liquid total scattering analysis. Physics and Chemistry of Minerals, 49:9. doi:10.1007/s00269-022-01186-6](https://link.springer.com/content/pdf/10.1007/s00269-022-01186-6.pdf)
  Please cite this paper if you use LiquidDiffraction in your work!

### Minor Changes

- LiquidDiffract.core.core now uses a separate objective function for density and/or background refinement - core.refinement_objfun. 
As a result, the arguments taken by core.calc_impr_interference_function have been redefined. The docs and example scripts have been updated to reflect this.
- Fixed some documentation issues and removed some dead code


## v1.1.8 (December 21, 2021)

### Major Changes

- Gaussian fitting toolbox
The calculation of the x-ray weighting factors for partial pair correlations used by the curve-fitting toolbox has been altered slightly. The default behaviour is now to calculate the WKM approximation for the effective atomic number of each species, _K_p,_ at _Q=0_. The previous behaviour of calculating the average _K_p_ across the whole _Q_-range of the data is still available by selecting the option in the _Additional Preferences_ dialog. However, this option may lead to unphysical results for some sample if the _Q_-range is >~10, as the x-ray weights are affected by the _Q_-range.

### Minor Changes

- Fixed a bug in the gaussian fitting toolbox where the plot view reset when the first peak was added


## v1.1.7 (September 16, 2021)

### Major Changes

- The _Structural Information_ Tab
  This new tab provides functionality to extract information on average bond lengths and coordination number from the calculated _RDF(r)_ and/or _T(r)_ functions.

  There are two toolboxes provided:
  - An _Integration Toolbox_: Integrate over the first peak in the _RDF(r)_ or _T(r)_ to estimate the bond-length and average coordination number of the first coordination shell for monatomic samples. Three different methods of calculating the coordination number are provided. The limits can be set manually, but there is also an auto-find feature to automatically refine the integration limits using optimisation/root-finding routines.
  - A _Curve-fitting Toolbox_: Fit an arbitrary number of gaussian-type peaks to the _RDF(r)_ or _T(r)_ to estimate bond lengths and coordination numbers of correlated atomic pairs in polyatomic samples. X-ray weighting factors are applied automatically and the bond length and coordination number are actual fitting parameters. Any combination of individual peaks or parameters can be fixed or refined during the fit. There is also an optional skewness parameter to fit a skew-normal distribution that can be toggled for any peak. 

- Fixed several bugs in the calculation of the S(Q) and g(r) for polyatomic samples

- LiquidDiffract is now distributed on PyPI

### Minor Changes

- Added an option to calculate PDF functions directly without refining the _g(r)_ / _S(Q)_
- Added an option to load data in nm<sup>-1</sup>
- Added an option to plot re-scaled Ashcroft-Langreth functions
- Added an option to use a modification function in the iterative refinement procedure
- Added a button to copy the refined value of &rho; to the inital value input
- Faber-Ziman formalism is now the default
- The composition table is now scroll-able and easier to use
- Fixed several minor bugs in the GUI
- Fixed several Qt garbage collection errors when running LiquidDiffract on Windows
- Simplified the iterative refinement procedure in the code
- Changed references to the differential correlation function from _F(r)_ to _D(r)_ to be more consistent with the literature
- Updated the documentation


## Release v1.0.0 (December 24, 2019)