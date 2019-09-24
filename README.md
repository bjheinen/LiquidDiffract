<p align="center"><img src="https://raw.githubusercontent.com/bjheinen/LiquidDiffract/master/liquid_diffract/data/icons/logo.png"></p>

<p align="center">
A GUI program to treat experimental X-ray diffraction data of liquid structures. 
LiquidDiffract can perform background subtraction, fourier transformation, data (quantitative structure factor) optimisation/normalisation, and estimation of density.
</p>


### Maintainer

Benedict J. Heinen (benedict.heinen@gmail.com)


## Table of Contents

* [Requirements](#requirements)
* [Licence](#licence)
* [Usage](#usage)
  * [Loading the GUI](#loading-the-gui)
  * [Basic Usage](#basic-usage)
  * [Background Subtraction Tab](#background-subtraction-tab)
  * [Refinement Tab](#refinement-tab)
    * [Composition toolbox](#composition-toolbox)
    * [Data Options](#data-options)
    * [Iterative structure factor refinement](#iterative-structure-factor-refinement)
    * [Density (&rho;) refinement](#density-rho-refinement)
      * [Global optimisation capability](#global-optimisation-capability)
    * [Terminal & Log-file output](#terminal-log-file-output)
  * [PDF Calculation (Outut) Tab](#pdf-calculation-output-tab)
* [References](#references)

## Requirements

* Python >= 3.7       | https://www.python.org
* Scipy >= 1.2.1      | https://www.scipy.org
* Numpy >= 1.16.2     | https://numpy.org
* PyQt5 >= 5.12       | https://riverbankcomputing.com/software/pyqt/intro
* pyqtgraph >= 0.10.0 | http://www.pyqtgraph.org

LiquidDiffract may run with earlier versions of these python modules but is untested.
The software is system-independent and does not require installation. 
It has been tested on Linux, Mac, and Windows.

## Licence

LiquidDiffract is issued under the GNU General Public Licence, v3.
This program comes with absolutely no warranty or guarantee.
Copyright © 2019 – Benedict J Heinen
See the [licence file](../blob/master/LICENSE) for more information.

## Usage

### Loading the GUI

### Basic Usage

LiquidDiffract is designed to enable a streamlined and visual workflow to treat experimental X-Ray diffraction data of liquid and amorphous solids (glasses).

There are three main tabs which provide a selection of toolboxes for data operations at different stages of the workflow.

* Background scaling and subtraction
* Data operations, structure factor calculation, interference function optimisation, and density refinement
* PDF calculation and data output

Data are automatically plotted and tabs updated as operations are made. The graphical plots display coordinates in the upper-right corner. Click and drag or use the scroll-wheel to zoom in on a region. Double right-click to reset the view.


<p align="center"><img src="https://raw.githubusercontent.com/bjheinen/LiquidDiffract/master/liquid_diffract/data/icons/LiquidDiffract v0.1_refinement_tab.png"></p>


### Background Subtraction Tab

This tab allows data and (optionally) background files to be loaded in. The auto-scale feature sppeds up workflow but is not reccomended for a close fit. 

All data must be in Q-space. A toolbox is provided to convert raw experimental data of 2&theta; values to Q-space.

### Refinement Tab

This is the tab where most of the data processing takes place.


#### Composition toolbox

The sample composition must be set before the structure factor, S(Q), can be calculated. The sample density (in atoms per cubic angstrom) should also be set here. If the density is to be refined, this is used as the initial value passed to the solver.

#### Data Options

At high Q-values experimental data often becomes increasingly noisy because the relative contribution of coherent scattering to the experimental signal decreases. This increasingly noisy signal can lead to dramatic and anomalous oscillations near the first peak in F(r). It can therefore be benificial to truncate the data at high-Q. However truncation can cause spurious peaks in F(r) as a result of the fourier transform. The positions of these peaks are a function of Q-max, so the use of different values should be investigated. Discussion on the effect of Q-max can be found in ref. [1] and [5].

Spurious peaks or oscillations (e.g. from a slightly mislocated beam stop) in the low-Q region can be removed by setting a Q-min cuttoff.

The option to smooth the data will apply a Savitzky-Golay filter to the data. Options to change the length of the filter window and the order of the polynomial used to fit samples can be found in the *Additional Preferences* dialog, accessible from the *Tools* menu.

LiquidDiffract provides the option to use either Ashcroft-Langreth [6] or Faber-Ziman [7] formalisms of the structure factor. See Keen et al., 2001 [8] for discussions of differences in structure factor formalisms, or the LiquidDiffract source code for specific implementations. 

Applying a modification function to S(Q) before FFT can help suppress truncation ripples amd can correct a gradient in the high-Q region. Two functions are currently implements:

##### Lorch [11]

M(Q) = sin(pi*Q/Q_max) / (pi*Q/Q_max)

##### Cosine Window (Drewitt et al.) [12]

	M(Q) = 1 if Q < Q1
     	M(Q) = 0.5[1 + cos(x*pi / N-1)] if Q1 <= Q <= Qmax
            
      	where N is width of the window function
       	and x is an integer with values from 0 - (N-1) across the window

After calculating the structure factor the final tab can be used to output S(Q), the pair distribution function g(r), and the radial distiribution function RDF(r), as is.


#### Iterative structure factor refinement

The numerical iterative procedure used by LiquidDiffract to minimize the error in the determination of g(r) follows the one proposed by Eggert et al., 2002 [1-4]. This procedure is based on the assumption that a minimum distance, r-min, can be defined, which represents the largest distance (0 -- r-min) where no atom can be found. In a liquid, this should be the distance of the 1st coordination shell. Because no atom can be present in this region, no oscillation should be observed in the g(r) function. As a result, the function F(r < r-min) = -4&pi;r&rho; However oscillations are commonly observed in this region, due to systematic errors such as the effect of an experimentally limited Q range (Q-max < &inf;) on the determination of the normalisation factor, &alpha;. The iterative procedure calculates the difference between real and model data in the low-r region and scales S(Q) accordingly to reduce this.

To refine S(Q) the value of r-min should be set carefully, as it has a strong influence. The position of r-min should correspond to the base of the first coordinence sphere in the g(r).

The number of iterations in the procedure can also be set; a minimum of 3 is normally required for convergence.

A &Chi;^2 figure of merit, defined as the area under the curve &Delta;F(r) for r<r-min, is used to rate the refinement.

#### Density (&rho;) refinement

The sample density can be determined by finding the value of &rho; that provides the best convergence of the iterative procedure described above. This is done by minimising the resultant value of &Chi;^2. LiquidDiffract supports several different solvers to do this. The solver in use, along with specific options like convergence criteria and number of iterations, can be selected from the *Additional Preferences* dialog. The solvers currently supported are:

L-BFGS-B [14-15]

SLSQP [16]

COBYLA [17-19]

All solvers require upper and lower bounds on the density to be set.

For more information on these optimisation algorithms please see the [SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html), or consult the references listed.

##### Global optimisation capability

The optimisation algorithms used to estimate density work best when the initial guess, &rho;_0 is close to the true value. If the density of the material is very poorly known it is reccomended to use the global optimisation option provided. It is not reccomended if a good initial guess can be made, due to slower, and sometimes poor, convergence.

The global optimisation procedure used is the basin-hopping algorithm provided by the SciPy library [20-21]. Basin-hopping is a stochastic algorithm which attempts to find the global minimum of a function by applying an iterative process with each cycle composed of the following features:

1. random perturbation of the coordinates
2. local minimization
3. accept or reject the new coordinates based on the minimized function value

The acceptance test used here is the Metropolis criterion of standard Monte Carlo algorithms [21].

For more information see the [SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html), or consult the references listed.


#### Terminal & Log-file output

An log is automatically generated for any refinement made. This log includes information on the data file, sample composition, data and refinement options used, solver ouput/convergence info (if refining density), and the final &Chi;^2 and &rho;.

The log for each refinement is automatically written to file. The default behaviour is to store each log in a file named 'refinement.log' within the current data directory. Each log is preceded by a timestamp. The log-mode can be changed to *Overwrite* in the *Additional Preferences* dialog. This creates a new log file for each data file loaded, which will be overwritten if already present. The filenames generated are of the form 'DATAFILENAME_refinement.log' and are similarly created in the source directory of the loaded data file.

LiquidDiffract outputs some density refinement information to the terminal. There is an option provided to print convergence progress messages to the terminal also. This can be useful for real-time monitoring (e.g. when using the slow global optimisation with a known number of iterations).

### PDF Calculation (Outut) Tab

The final tab displays the optimised S(Q), g(r), and RDF(r). The buttons at the bottom of the window allow each one to be saved to a text file. If a modification function has been used in the data treatment then information on this will also be saved, along with the raw S(Q).





## References


[1] Eggert, J. H., Weck, G., Loubeyre, P., & Mezouar, M. (2002). Quantitative structure factor and density measurements of high-pressure fluids in diamond anvil cells by x-ray diffraction: Argon and water. Physical Review B, 65(17), 174105.

[2] Morard, G., Siebert, J., Andrault, D., Guignot, N., Garbarino, G., Guyot, F., & Antonangeli, D. (2013). The Earth's core composition from high pressure density measurements of liquid iron alloys. Earth and Planetary Science Letters, 373, 169-178.

[3] Decremps, F., Morard, G., Garbarino, G., & Casula, M. (2016). Polyamorphism of a Ce-based bulk metallic glass by high-pressure and high-temperature density measurements. Physical Review B, 93(5), 054209.

[4] Kaplow, R., Strong, S. L., & Averbach, B. L. (1965). Radial density functions for liquid mercury and lead. Physical Review, 138(5A), A1336.

[5] Shen, G., Rivers, M. L., Sutton, S. R., Sata, N., Prakapenka, V. B., Oxley, J., & Suslick, K. S. (2004). The structure of amorphous iron at high pressures to 67 GPa measured in a diamond anvil cell. Physics of the Earth and Planetary Interiors, 143, 481-495.

[6] Ashcroft, N. W., & Langreth, D. C. (1967). Structure of binary liquid mixtures. I. Physical Review, 156(3), 685.

[7] Faber, T. E., & Ziman, J. M. (1965). A theory of the electrical properties of liquid metals: III. The resistivity of binary alloys. Philosophical Magazine, 11(109), 153-173.

[8] Keen, D. A. (2001). A comparison of various commonly used correlation functions for describing total scattering. Journal of Applied Crystallography, 34(2), 172-177.

[9] hubbel compton data

[10] form factors data

[11] Lorch, E. (1969). Neutron diffraction by germania, silica and radiation-damaged silica glasses. Journal of Physics C: Solid State Physics, 2(2), 229.

[12] Drewitt, J. W., Sanloup, C., Bytchkov, A., Brassamin, S., & Hennet, L. (2013). Structure of (FexCa1−xO)y(SiO 2)1−y liquids and glasses from high-energy x-ray diffraction: Implications for the structure of natural basaltic magmas. Physical Review B, 87(22), 224201.

[13] Morard, G., Garbarino, G., Antonangeli, D., Andrault, D., Guignot, N., Siebert, J., Roberge, M., Boulard, E., Lincot, A., Denoeud, A. & Petitgirard, S. (2014). Density measurements and structural properties of liquid and amorphous metals under high pressure. High Pressure Research, 34(1), 9-21.

[14] Byrd, R H and P Lu and J. Nocedal. 1995. A Limited Memory Algorithm for Bound Constrained Optimization. SIAM Journal on Scientific and Statistical Computing 16 (5): 1190-1208.

[15] Zhu, C and R H Byrd and J Nocedal. 1997. L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization. ACM Transactions on Mathematical Software 23 (4): 550-560.

[16] Kraft, D. A software package for sequential quadratic programming. 1988. Tech. Rep. DFVLR-FB 88-28, DLR German Aerospace Center – Institute for Flight Mechanics, Koln, Germany.

[17] Powell, M J D. A direct search optimization method that models the objective and constraint functions by linear interpolation. 1994. Advances in Optimization and Numerical Analysis, eds. S. Gomez and J-P Hennart, Kluwer Academic (Dordrecht), 51-67.

[18] Powell M J D. Direct search algorithms for optimization calculations. 1998. Acta Numerica 7: 287-336.

[19] Powell M J D. A view of algorithms for optimization without derivatives. 2007.Cambridge University Technical Report DAMTP 2007/NA03

[20] Wales, D J, and Doye J P K, Global Optimization by Basin-Hopping and the Lowest Energy Structures of Lennard-Jones Clusters Containing up to 110 Atoms. Journal of Physical Chemistry A, 1997, 101, 5111.

[21] Li, Z. and Scheraga, H. A., Monte Carlo-minimization approach to the multiple-minima problem in protein folding, Proc. Natl. Acad. Sci. USA, 1987, 84, 6611.



[Features](#features) | [Installation](#installation) | [Usage](#usage) | [Examples](#examples) | [Command-line options](#options) | [Configuration](#configuration)




