# SARA-Calib-RI
Joint DDE calibration and imaging via non-convex optimisation applied to real RI observations.
To launch Joint calibration and imaging, an initial estimate of the model image need to be provided in *./image_init/*. Consider running imaging only beforehand. 

-Data and observation specs need to be saved in *.mat* files and saved in *./data*.
A python script from pyxis library is used to build these *.mat* files from  a given MS (https://github.com/ska-sa/pyxis)

-The source code relies on an external library for NUFFT (Fessler et al. 2003).

-A toy example can be launched by running **main_example.m**. 
Considered data are a subset of  3c391 data from VLA tutorial https://casaguides.nrao.edu/index.php/VLA_Continuum_Tutorial_3C391-CASA4.6


**Associated paper:**
> A. Dabbech, A. Repetti, R. Perley, O. M. Smirnov, & Y. Wiaux, Cygnus A jointly calibrated and imaged via non-convex optimisation from JVLA data</a>, <i>Monthly Notices of the Royal Astronomical Society </i>, Jan. 2021 (Submitted).
# <a href="https://arxiv.org/abs/1701.03689">
