import numpy as np
import os
import sys
from astropy.io import fits
#from scc_funs import secchi_rectify, fill_from_defhdr, rebinIDL, scc_getbkgimg, scc_sebip
#from wcs_funs import get_Suncent, fitshead2wcs
import datetime
from scipy.interpolate import griddata


def euvi_prep(im, hdr, pointingOff=False, ):
    if hdr['detector'] != 'EUVI':
        sys.exit('Calibration for EUVI DETECTOR only')
        
    # Pointing correction on 
    if not pointingOff:
        print (a)
    
    print (sd)