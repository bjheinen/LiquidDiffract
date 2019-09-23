<p align="center"><img src="https://github.com/bjheinen/LiquidDiffract/blob/master/liquid_diffract/data/icons/logo.png"></p>

A GUI program to treat experimental X-ray diffraction data of liquid structures. LiquidDiffract can perform background subtraction, fourier transformation, data (quantitative structure factor) optimisation/normalisation, and estimation of density.


### Maintainer

Benedict J. Heinen (benedict.heinen@gmail.com)


## Table of Contents

* [Requirements](#requirements)
* [Usage](#usage)
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


## Usage

### Loading the GUI



load data, scale (including autoscale) and subtract background data.

Also a widget to convert 2theta data to Q-space



Adding composition
setting intial density


options to apply a Q_Max, minimum q cutoff

smooth the data applies a savitzky golay filter

use a modification function (need refs for this one too)

use one of two structure factor formalisms


optimisation options toolbox

set r_min cut off and the number of iterations for the eggert procedure


refine density using non-linear solver

(and apply bounds)

these are shown in the results toolbox


last tab displays and allows to save S(Q), g(r), and RDF(r)









Click and drag or use the scroll-wheel to zoom in on a region. Double right-click to reset the view.








## References

Morard, 2013 but note there is a missing bracket - refer to the definition
    in Decremps et al., 2016
    
    
drewitt et al., 2013



[1] Eggert, J. H., Weck, G., Loubeyre, P., & Mezouar, M. (2002). Quantitative structure factor and density measurements of high-pressure fluids in diamond anvil cells by x-ray diffraction: Argon and water. Physical Review B, 65(17), 174105.

[2] Morard, G., Siebert, J., Andrault, D., Guignot, N., Garbarino, G., Guyot, F., & Antonangeli, D. (2013). The Earth's core composition from high pressure density measurements of liquid iron alloys. Earth and Planetary Science Letters, 373, 169-178.

[3] Decremps, F., Morard, G., Garbarino, G., & Casula, M. (2016). Polyamorphism of a Ce-based bulk metallic glass by high-pressure and high-temperature density measurements. Physical Review B, 93(5), 054209.

references for data used


[4] Kaplow, R., Strong, S. L., & Averbach, B. L. (1965). Radial density functions for liquid mercury and lead. Physical Review, 138(5A), A1336.



[5] Shen, G., Rivers, M. L., Sutton, S. R., Sata, N., Prakapenka, V. B., Oxley, J., & Suslick, K. S. (2004). The structure of amorphous iron at high pressures to 67 GPa measured in a diamond anvil cell. Physics of the Earth and Planetary Interiors, 143, 481-495.

[5] shows variations arising from qmax and rmin changes


[6] Drewitt, J. W., Sanloup, C., Bytchkov, A., Brassamin, S., & Hennet, L. (2013). Structure of (FexCa1−xO)y(SiO 2)1−y liquids and glasses from high-energy x-ray diffraction: Implications for the structure of natural basaltic magmas. Physical Review B, 87(22), 224201.


[7] Keen, D. A. (2001). A comparison of various commonly used correlation functions for describing total scattering. Journal of Applied Crystallography, 34(2), 172-177.


See Keen et al., 2001 [7] for discussions of differences in structure factor formalisms.



blems and even
allows us to experimentally determine r0. Following the
general outline of a method for minimizing errors in the
determination of g(r) pioneered by Kaplow et al.,
23 we force
the behavior of F(r) at small r ~below the first intermolecular peak! to match the expected behavior for a given sample.



Faber, T. E., & Ziman, J. M. (1965). A theory of the electrical properties of liquid metals: III. The resistivity of binary alloys. Philosophical Magazine, 11(109), 153-173.


Ashcroft, N. W., & Langreth, D. C. (1967). Structure of binary liquid mixtures. I. Physical Review, 156(3), 685.




other reading

Morard, G., Garbarino, G., Antonangeli, D., Andrault, D., Guignot, N., Siebert, J., ... & Petitgirard, S. (2014). Density measurements and structural properties of liquid and amorphous metals under high pressure. High Pressure Research, 34(1), 9-21.










