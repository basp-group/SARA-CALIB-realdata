# SARA-Calib-RI
Joint DDE calibration and imaging via non-convex optimisation applied to real RI observations.
To launch Joint calibration and imaging, an initial estimate of the model image need to be provided in *./image_init/*. Consider running imaging only beforehand. 

1.Data and observation specs need to be saved in *.mat* files and saved in *./data*.
A python script *./pyxis/py-ckat.py* from pyxis library is used to build these *.mat* files from  a given MS. 
*Pyxis* need to be installed from .

3. A toy example can be launched by running **main_example.m**. 
Considered data are a subset of  3c391 data from VLA tutorial https://casaguides.nrao.edu/index.php/VLA_Continuum_Tutorial_3C391-CASA4.6

External lib requirements
1. Non-Uniform FFT (Fessler et al. 2003): The needed files are included in the source code (*./lib/lib_external/*). The full toolbox can be found at http://web.eecs.umich.edu/~fessler/irt/fessler.tgz.
2. Pyxis library can be used to generate the *.mat* files of data and observations specs. The library can be found at https://github.com/ska-sa/pyxis. The script used to generate the *.mat* files are in *./pyxis/* .
**Associated paper:**
> A. Dabbech, A. Repetti, R. Perley, O. M. Smirnov, & Y. Wiaux, Cygnus A jointly calibrated and imaged via non-convex optimisation from JVLA data</a>, <i>Monthly Notices of the Royal Astronomical Society </i>, Jan. 2021 (Submitted).
# <a href="https://arxiv.org/abs/1701.03689">
