## example to launch from command line : pyxis MS='My_MS'   SRC='3c391'   TAG='ch1' get_data
## Rodrigues compatibale simulation pipeline 
## sphesihle makhathini [sphemakh@gmail.com]

import Pyxis
import ms
import mqt 
import im
import imager
import lsm
import im.argo
import im.lwimager
from Pyxis.ModSupport import *
from astropy.io import fits
from astropy.io import ascii
from casacore.tables import *
import scipy.io as sio
import pyrap.tables
import pyfits
import Tigger
import numpy
import os
import math
import json
import time


PI = math.pi
FWHM = math.sqrt( math.log(256) )
LightSpeed = 299792458


def get_data(msfile='$MS',src_name='$SRC',out_data_set_tag='$TAG'):

    dataSetFolder = II('%s'%out_data_set_tag)
    msname = II('%s'%msfile)
    srcFolder = II('%s'%src_name)
    x.sh("mkdir -p  ./%s"%srcFolder)
    x.sh("mkdir -p  ./%s/%s"%(srcFolder,dataSetFolder))

    tab = ms.ms(msname,write=False);
    info("columns are",*(tab.colnames()));

    uvw = tab.getcol("UVW")


    
    spwtab = ms.ms(subtable="SPECTRAL_WINDOW")
    freq0 = spwtab.getcol("CHAN_FREQ")[ms.SPWID,0]#get the first freq
    wavelength = LightSpeed/freq0
    freqStepVect = spwtab.getcol("CHAN_WIDTH")
    freqsVect=spwtab.getcol("CHAN_FREQ")
    info(freqsVect)
    

    #Specs
    nchan = len(freqsVect[0,:]);
    bandwidth = spwtab.getcol("CHAN_WIDTH")[ms.SPWID,0]
    spwtab.close()

    field = tab.getcol("FIELD_ID");
    ant1=tab.getcol("ANTENNA1");
    ant2=tab.getcol("ANTENNA2");
    time    = tab.getcol("TIME");
    scantab = tab.getcol("SCAN_NUMBER")
    integrationTimeVect = tab.getcol("EXPOSURE");
    dt = tab.getcol("EXPOSURE",0,1)[0]
    dtf = (tab.getcol("TIME",tab.nrows()-1,1)-tab.getcol("TIME",0,1))[0]
    info(">>>>   Freq %.2f MHz (lambda=%.2fm), bandwidth %.2g kHz, %.2fs integrations, %.2fh synthesis"%(freq0*1e-6,wavelength,bandwidth*1e-3,dt,dtf/3600))

    
    sio.savemat("./%s/%s/msSpecs.mat"%(srcFolder,dataSetFolder),{'field':field,'ant1':ant1.astype(float),'ant2':ant2.astype(float),'freqsVect':freqsVect,'time':time,'uvw':uvw,'scans':scantab})

 
    #Data
    weightfull = tab.getcol("WEIGHT_SPECTRUM");
    flag_row = tab.getcol("FLAG_ROW")
    flagAll= tab.getcol("FLAG")
    data_full = tab.getcol("DATA");
    
    for ifreq in range(0, nchan):

        data_selected   = data_full[:,ifreq,:];
        weight_selected = weightfull[:,ifreq,:];
        weightI= (weight_selected[:,0]+weight_selected[:,3])
        dataI_w = (weight_selected[:,0]*data_selected[:,0]+weight_selected[:,3]*data_selected[:,3])/weightI ;
 
 
        flag_bis  = numpy.logical_or(flagAll[:,ifreq,0],flagAll[:,ifreq,3]);
        flag = numpy.logical_or(flag_bis,flag_row)
        flag = flag.astype(float)
        
        sio.savemat("./%s/%s/data.mat"%(srcFolder,dataSetFolder),{'y_I':dataI_w, 'weights':weightI.astype(float), 'flag':flag});

       
        info("Read data .. writing file..Freq %s"%(ifreq))
    
    tab.close()


   


