
<p align="center"><img src="https://raw.githubusercontent.com/bjheinen/LiquidDiffract/master/LiquidDiffract/resources/icons/logo.png"></p>

<p align="center">
A GUI program to treat experimental X-ray diffraction data of liquid structures. 
LiquidDiffract can perform background subtraction, fourier transformation, data (quantitative structure factor) optimisation/normalisation, and estimation of density.
</p>


### Maintainer

Benedict J. Heinen (benedict.heinen@gmail.com)


## Table of Contents

* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)
  * [Loading the GUI](#loading-the-gui)
  * [Basic Usage](#basic-usage)
  * [Background Subtraction Tab](#background-subtraction-tab)
  * [Refinement Tab](#refinement-tab)
    * [Composition Toolbox](#composition-toolbox)
    * [Data Options](#data-options)
    * [Iterative Structure Factor Refinement](#iterative-structure-factor-refinement)
    * [Density (&rho;) Refinement](#density-%CF%81-refinement)
      * [Global optimisation capability](#global-optimisation-capability)
    * [Terminal & Log-file Output](#terminal--log-file-output)
  * [PDF Calculation (Output) Tab](#pdf-calculation-output-tab)
  * [Using LiquidDiffract Core Library](#using-liquiddiffract-core-library)
* [References](#references)
* [License](#license)



## Requirements

* [Python >= 3.5](https://www.python.org)
* [SciPy >= 1.2.1](https://www.scipy.org)
* [NumPy >= 1.16.2](https://numpy.org)
* [PyQt5 >= 5.12](https://riverbankcomputing.com/software/pyqt/intro)
* [pyqtgraph >= 0.10.0](http://www.pyqtgraph.org)
* [importlib_resources](https://pypi.org/project/importlib_resources/) is required if using Python<3.7

LiquidDiffract should run with earlier versions of these python packages but is untested. If you do not have a python installation, we recommend [Anaconda](https://www.anaconda.com/distribution/).
LiquidDiffract is system-independent and has been tested on Linux, Mac, and Windows.

## Installation

The simplest way to install LiquidDiffract is with pip and git:

```$ pip install git+https://github.com/bjheinen/LiquidDiffract.git```

which should also deal with dependencies automatically.

If you have issues clone this repository or download the source and install locally.

The source is directly available [here](https://github.com/bjheinen/LiquidDiffract/archive/master.zip)

 We still reccomend pip for installing from a local directory.
***Invoke setup.py directly at your own risk!***

<details><summary><b>Testing the software or modifying source code</b></summary>

You can make a development (editable) install by using pip's *-e* flag:

```$ pip install -e git+https://github.com/bjheinen/LiquidDiffract#egg=LiquidDiffract```

From a local directory:

```
$ cd /path/to/local/directory/
$ curl -Lo LiquidDiffract.zip https://github.com/bjheinen/LiquidDiffract/archive/master.zip
$ unzip LiquidDiffract.zip 
$ mv LiquidDiffract-*/ LiquidDiffract/
$ pip install -e /path/to/local/directory/LiquidDiffract/
```

This will not actually install anything, but create a special .egg-link file in the deployment directory, that links to LiquidDiffract's source code.

It is useful if you want to mess around with, or make changes to the source code without having to re-install. It is also good for testing the software, as it makes uninstalling simpler. Use Python's site-packages directory as the deployment directory if you want your editable install of LiquidDiffract available on your sys.path for other programs using your Python installation. To do this from the github page use the *-t* flag ```-t /path/to/directory```

</details>

## Usage

### Loading the GUI

To run the software use the command

```$ LiquidDiffract```

or run *'LiquidDiffract.py'* from the installation directory ```$ python LiquidDiffract.py```

### Basic Usage

LiquidDiffract is designed to enable a streamlined and visual workflow to treat experimental X-Ray diffraction data of liquid and amorphous solids (glasses).

There are three main tabs which provide a selection of toolboxes for data operations at different stages of the workflow.

* Background scaling and subtraction
* Data operations, structure factor calculation, interference function optimisation, and density refinement
* PDF calculation and data output

Data are automatically plotted and tabs updated as operations are made. The graphical plots display coordinates in the upper-right corner. Click and drag or use the scroll-wheel to zoom in on a region. Double right-click to reset the view.


<p align="center"><img src="https://raw.githubusercontent.com/bjheinen/LiquidDiffract/master/LiquidDiffract/resources/icons/LiquidDiffract v0.1_refinement_tab.png"></p>


### Background Subtraction Tab

This tab allows data and (optionally) background files to be loaded in. The auto-scale feature speeds up workflow but is not recommended for a close fit. 

All data must be in Q-space. A toolbox is provided to convert raw experimental data of 2&theta; values to Q-space.

### Refinement Tab

This is the tab where most of the data processing takes place.


#### Composition Toolbox

The sample composition must be set before the structure factor, S(Q), can be calculated. The sample density (in atoms per cubic angstroms) should also be set here. If the density is to be refined, this is used as the initial value passed to the solver.

#### Data Options

At high Q-values experimental data often becomes increasingly noisy because the relative contribution of coherent scattering to the experimental signal decreases. This increasingly noisy signal can lead to dramatic and anomalous oscillations near the first peak in F(r). It can therefore be beneficial to truncate the data at high-Q. However truncation can cause spurious peaks in F(r) as a result of the fourier transform. The positions of these peaks are a function of Q-max, so the use of different values should be investigated. Discussion on the effect of Q-max can be found in ref. [1] and [5].

Spurious peaks or oscillations (e.g. from a slightly mislocated beam stop) in the low-Q region can be removed by setting a Q-min cut-off.

The option to smooth the data will apply a Savitzky-Golay filter to the data. Options to change the length of the filter window and the order of the polynomial used to fit samples can be found in the *Additional Preferences* dialog, accessible from the *Tools* menu.

LiquidDiffract provides the option to use either Ashcroft-Langreth [6] or Faber-Ziman [7] formalisms of the structure factor. See Keen et al., 2001 [8] for discussions of differences in structure factor formalisms, or the LiquidDiffract source code for specific implementations. 

The total molecular structure factor, S(Q), is converted to normalised units using the method of Krogh-Moe [9] and Norman [10].

Applying a modification function to S(Q) before FFT can help suppress truncation ripples and can correct a gradient in the high-Q region. Two functions are currently implemented:

**Lorch [13]:**

><img src='http://latex.codecogs.com/svg.latex?M\left(Q\right)=\frac{\sin\left(\tfrac{\pi%20Q}{Q_{max}}\right)}{\left(\tfrac{\pi%20Q}{Q_{max}}\right)}'/>


**Cosine Window (Drewitt et al.) [14]:**

><img src='http://latex.codecogs.com/svg.latex?M\left(Q\right)%20=\begin{cases}1&if\;\,Q%3CQ_{window\_start}\\0.5\left[1+\cos\left(\frac{x\pi}{N-1}\right)\right]&if\;\,Q_{window\_start}%3CQ\leq%20Q_{max}\\0&if\;\,Q%3EQ_{max}\end{cases}'/>
>
>where N is width of the window function
>and x is an integer with values from 0 : (N-1) across the window

After calculating the structure factor the final tab can be used to output S(Q), the pair distribution function g(r), and the radial distribution function RDF(r), as is.


#### Iterative Structure Factor Refinement

The numerical iterative procedure used by LiquidDiffract to minimize the error in the determination of g(r) follows the one proposed by Eggert et al., 2002 [1-4]. This procedure is based on the assumption that a minimum distance, r-min, can be defined, which represents the largest distance (0 -- r-min) where no atom can be found. In a liquid, this should be the distance of the 1st coordination shell. Because no atom can be present in this region, no oscillation should be observed in the g(r) function. As a result, the function *F(r < r-min)* = -4&pi;r&rho; However oscillations are commonly observed in this region, due to systematic errors such as the effect of an experimentally limited Q range (Q-max < &inf;) on the determination of the normalisation factor, &alpha;. The iterative procedure calculates the difference between real and model data in the low-r region and scales S(Q) accordingly to reduce this.

To refine S(Q) the value of r-min should be set carefully, as it has a strong influence. The position of r-min should correspond to the base of the first coordinence sphere in the g(r).

The number of iterations in the procedure can also be set; a minimum of 3 is normally required for convergence.

A *&Chi;<sup>2</sup>* figure of merit, defined as the area under the curve *&Delta;F(r) for r<r-min*, is used to rate the refinement.

#### Density (&rho;) Refinement

The sample density can be determined by finding the value of &rho; that provides the best convergence of the iterative procedure described above. This is done by minimising the resultant value of *&Chi;<sup>2</sup>*. LiquidDiffract supports several different solvers to do this. The solver in use, along with specific options like convergence criteria and number of iterations, can be selected from the *Additional Preferences* dialog. The solvers currently supported are:

* L-BFGS-B [16-17]

* SLSQP [18]

* COBYLA [19-21]

All solvers require upper and lower bounds on the density to be set.

For more information on these optimisation algorithms please see the [SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html), or consult the references listed.

##### Global optimisation capability

The optimisation algorithms used to estimate density work best when the initial guess, &rho;<sub>0</sub> is close to the true value. If the density of the material is very poorly known it is recommended to use the global optimisation option provided. It is not recommended if a good initial guess can be made, due to slower, and sometimes poor, convergence.

The global optimisation procedure used is the basin-hopping algorithm provided by the SciPy library [22-23]. Basin-hopping is a stochastic algorithm which attempts to find the global minimum of a function by applying an iterative process with each cycle composed of the following features:

1. random perturbation of the coordinates
2. local minimization
3. accept or reject the new coordinates based on the minimized function value

The acceptance test used here is the Metropolis criterion of standard Monte Carlo algorithms [23].

For more information see the [SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html), or consult the references listed.


#### Terminal & Log-file Output

An log is automatically generated for any refinement made. This log includes information on the data file, sample composition, data and refinement options used, solver output/convergence info (if refining density), and the final *&Chi;<sup>2</sup>* and *&rho;*.

The log for each refinement is automatically written to file. The default behaviour is to store each log in a file named 'refinement.log' within the current data directory. Each log is preceded by a time-stamp. The log-mode can be changed to *Overwrite* in the *Additional Preferences* dialog. This creates a new log file for each data file loaded, which will be overwritten if already present. The file-names generated are of the form 'DATAFILENAME_refinement.log' and are similarly created in the source directory of the loaded data file.

LiquidDiffract outputs some density refinement information to the terminal. There is an option provided to print convergence progress messages to the terminal also. This can be useful for real-time monitoring (e.g. when using the slow global optimisation with a known number of iterations).

### PDF Calculation (Output) Tab

The final tab displays the optimised S(Q), g(r), and RDF(r). The buttons at the bottom of the window allow each one to be saved to a text file. If a modification function has been used in the data treatment then information on this will also be saved, along with the raw S(Q).


### Using LiquidDiffract Core Library

The module LiquidDiffract.core.core is the core library which provides all the specific functions for processing (liquid) diffraction data. This can be imported and used in your own Python code. There are also some data helper functions (e.g. to rebin data) found in LiquidDiffract.core.data_utils

A brief example of using LiquidDiffract in a custom data processing script is given below:

<details><summary><b>Example usage</b></summary>

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Import required python modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.optimize import minimize

# import the LiquidDiffract core module
import LiquidDiffract.core.core as liquid

# First set composition
# Composition dictionary should have element names as keys, with values as tuples of the form (Z, charge, fraction)
composition = {'Ga': (31,0,1)}
# Set initial density
rho = 0.05

# Load your background subtracted data
q_raw, I_raw = np.loadtxt('example_data.dat', unpack=True, skiprows=0)

# Next rebin data and trim if necessary

# Re-bin via interpolation so dq (steps of q) is consistent
# Suggested dq = 0.02
dq = 0.02
# Cut-off below Q_cutoff
q_cutoff = 9.0
q_data = np.arange(0, q_raw[-1], dq)
q_data = q_data[q_data<q_cutoff]
finterp = interp1d(q_raw, I_raw, kind='cubic', fill_value='extrapolate')
I_data = finterp(q_data)

# Apply any other data treatment
# e.g. a savitsky-golay filter to smooth the data
window_length = 31
polyorder = 3
I_data = savgol_filter(I_data, window_length, polyorder)

# First calculate the interference function i(Q)
# i(Q) = S(Q) - S_inf
structure_factor = liquid.calc_structure_factor(q_data, I_data, composition, rho)
interference_func = structure_factor - liquid.calc_S_inf(composition, q_data)

# store original interference function
interference_func_0 = interference_func

# A low r region cutoff must be chosen for the refinement method of Eggert et al. (2002)
# A density value in atoms(or molecules) / A^3 must also be set. This will be
# used as the starting estimate if the density is refined.
r_min = 2.3
rho_0 = 0.06

# The function 'calc_impr_interference_func' calculates an improved estimate of
# the interference function via the iterative procedure described in Eggert et al., 2002.
#
# calc_impr_interference_func requires the following arguments:
# q_data - q values
# I_data / i(Q) - the treated intensity data or interference function (depends on opt_flag used)
# composition - composition dictionary
# r_min - low r cut-off
# iter_limit - iteration limit for Eggert procedure
# method - method of calculating S(Q) ('ashcroft-langreth' and 'faber-ziman' are currently supported)
# mod_func - modification function to use ('None', 'Cosine-window' or 'Lorch')
# window_start - window_start of cosine function (if in use)
# opt_flag - set to 0 if not running from solver
#
# The last positional argument (opt_flag) signals that you only want the chi_squared value to be returned.
# This is useful if using a solver to refine the density by minimising the chi_squared value

# e.g.
# Return the improved i(Q) & chi^2 at density = rho_0
iter_limit = 20
method = 'ashcroft-langreth'
mod_func = 'Cosine-window'
window_start = 9
args = (q_data, interference_func_0, composition, r_min, iter_limit, method, mod_func, window_start, 0)
# Store the refined interference function at rho_0
interference_func_1, chi_sq_1 = liquid.calc_impr_interference_func(rho_0, *args)

# Next we use opt_flag = 1 and pass the function calc_impr_interference_func to
# a solver to estimate the density
args = (q_data, I_data, composition, r_min, iter_limit, method, mod_func, window_start, 1)
# Set-up bounds and other options according to the documentation of solver/minimisation routine
bounds = ((0.03, 0.08),)
op_method = 'L-BFGS-B'
optimisation_options = {'disp': 1,
                        'maxiter': 15000,
                        'maxfun': 15000,
                        'ftol': 2.22e-12,
                        'gtol': 1e-12
                        }
opt_result = minimize(liquid.calc_impr_interference_func, rho_0,
                      bounds=bounds, args=args,
                      options=optimisation_options,
                      method=op_method)
# The solver finds the value of rho that gives the smallest chi^2
rho_refined = opt_result.x[0]

# The interference function can then be re-calculated using the new value for rho,
# and the F(r) refined as above.
interference_func = (liquid.calc_structure_factor(q_data,I_data, composition, rho_refined) - 
                     liquid.calc_S_inf(composition, q_data))
args = (q_data, interference_func, composition, r_min, iter_limit, method, mod_func, window_start, 0)
interference_func_2, chi_sq_2 = liquid.calc_impr_interference_func(rho_refined, *args)

# Plot the data
# Initial interference function calculation
plt.plot(q_data, interference_func_0, color='g', label=r'Initial i(Q) | $ρ = {rho:.2f}$'.format(rho=rho_0))
# Optimised at rho_0
plt.plot(q_data, interference_func_1, color='r', label=r'Optimised i(Q) | $ρ = {rho:.2f}$ & $χ^{{2}} = {chisq:.2f}$'.format(rho=rho_0, chisq=chi_sq_1))
# Optimised, using refined density estimate
plt.plot(q_data, interference_func_2, color='b', label=r'Optimised i(Q) | $ρ = {rho:.3f}$ & $χ^{{2}} = {chisq:.2f}$'.format(rho=rho_refined, chisq=chi_sq_2))
# Add a legend
plt.legend(loc='best')
# Show the plots
plt.show()
```

</details>


## References


[1] Eggert, J. H., Weck, G., Loubeyre, P., & Mezouar, M. (2002). Quantitative structure factor and density measurements of high-pressure fluids in diamond anvil cells by x-ray diffraction: Argon and water. Physical Review B, 65(17), 174105.

[2] Morard, G., Siebert, J., Andrault, D., Guignot, N., Garbarino, G., Guyot, F., & Antonangeli, D. (2013). The Earth's core composition from high pressure density measurements of liquid iron alloys. Earth and Planetary Science Letters, 373, 169-178.

[3] Decremps, F., Morard, G., Garbarino, G., & Casula, M. (2016). Polyamorphism of a Ce-based bulk metallic glass by high-pressure and high-temperature density measurements. Physical Review B, 93(5), 054209.

[4] Kaplow, R., Strong, S. L., & Averbach, B. L. (1965). Radial density functions for liquid mercury and lead. Physical Review, 138(5A), A1336.

[5] Shen, G., Rivers, M. L., Sutton, S. R., Sata, N., Prakapenka, V. B., Oxley, J., & Suslick, K. S. (2004). The structure of amorphous iron at high pressures to 67 GPa measured in a diamond anvil cell. Physics of the Earth and Planetary Interiors, 143, 481-495.

[6] Ashcroft, N. W., & Langreth, D. C. (1967). Structure of binary liquid mixtures. I. Physical Review, 156(3), 685.

[7] Faber, T. E., & Ziman, J. M. (1965). A theory of the electrical properties of liquid metals: III. The resistivity of binary alloys. Philosophical Magazine, 11(109), 153-173.

[8] Keen, D. A. (2001). A comparison of various commonly used correlation functions for describing total scattering. Journal of Applied Crystallography, 34(2), 172-177.

[9] Krogh-Moe, J. (1956). A method for converting experimental X-ray intensities to an absolute scale. Acta Crystallographica, 9(11), 951-953.

[10] Norman, N. (1957). The Fourier transform method for normalizing intensities. Acta Crystallographica, 10(5), 370-373.

[11] Hubbell, J. H., Veigele, W. J., Briggs, E. A., Brown, R. T., Cromer, D. T., & Howerton, D. R. (1975). Atomic form factors, incoherent scattering functions, and photon scattering cross sections. Journal of physical and chemical reference data, 4(3), 471-538.

[12a] Brown, P. J., Fox, A. G., Maslen, E. N., O'Keefe, M. A., & Willis, B. T. M. (2006). 'Intensity of diffracted intensities' in *International Tables for Crystallography*, Vol. C., ch. 6.1, pp. 554-595 

[12b] This form factor data is also available freely at <http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php>

[13] Lorch, E. (1969). Neutron diffraction by germania, silica and radiation-damaged silica glasses. Journal of Physics C: Solid State Physics, 2(2), 229.

[14] Drewitt, J. W., Sanloup, C., Bytchkov, A., Brassamin, S., & Hennet, L. (2013). Structure of (FexCa1−xO)y(SiO 2)1−y liquids and glasses from high-energy x-ray diffraction: Implications for the structure of natural basaltic magmas. Physical Review B, 87(22), 224201.

[15] Morard, G., Garbarino, G., Antonangeli, D., Andrault, D., Guignot, N., Siebert, J., Roberge, M., Boulard, E., Lincot, A., Denoeud, A. & Petitgirard, S. (2014). Density measurements and structural properties of liquid and amorphous metals under high pressure. High Pressure Research, 34(1), 9-21.

[16] Byrd, R H and P Lu and J. Nocedal. 1995. A Limited Memory Algorithm for Bound Constrained Optimization. SIAM Journal on Scientific and Statistical Computing 16 (5): 1190-1208.

[17] Zhu, C and R H Byrd and J Nocedal. 1997. L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization. ACM Transactions on Mathematical Software 23 (4): 550-560.

[18] Kraft, D. A software package for sequential quadratic programming. 1988. Tech. Rep. DFVLR-FB 88-28, DLR German Aerospace Center – Institute for Flight Mechanics, Koln, Germany.

[19] Powell, M J D. A direct search optimization method that models the objective and constraint functions by linear interpolation. 1994. Advances in Optimization and Numerical Analysis, eds. S. Gomez and J-P Hennart, Kluwer Academic (Dordrecht), 51-67.

[20] Powell M J D. Direct search algorithms for optimization calculations. 1998. Acta Numerica 7: 287-336.

[21] Powell M J D. A view of algorithms for optimization without derivatives. 2007.Cambridge University Technical Report DAMTP 2007/NA03

[22] Wales, D J, and Doye J P K, Global Optimization by Basin-Hopping and the Lowest Energy Structures of Lennard-Jones Clusters Containing up to 110 Atoms. Journal of Physical Chemistry A, 1997, 101, 5111.

[23] Li, Z. and Scheraga, H. A., Monte Carlo-minimization approach to the multiple-minima problem in protein folding, Proc. Natl. Acad. Sci. USA, 1987, 84, 6611.



## License

Licensed under the GNU General Public License (GPL), version 3 or later.

See the [license file](../master/LICENSE) for more information.

This program comes with absolutely no warranty or guarantee.

Copyright © 2018-2019 – Benedict J Heinen


