import numpy as np
import sunpy.map
import sys, os
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from sunpy.time import parse_time

efacDir = '/Users/kaycd1/ssw/soho/lasco/idl/expfac/data/'
vigFile = '/Users/kaycd1/ssw/soho/lasco/idl/data/calib/c2vig_final.fts'

#from scc_funs import scc_make_array, scc_zelensky_array, rebinIDL
#from cor_prep import cor_prep
#from hi_prep import hi_prep
#from sunpy.coordinates import get_horizons_coord
#from sunspyce import get_sunspyce_hpc_point, get_sunspyce_roll, get_sunspyce_coord, get_sunspyce_lonlat, get_sunspyce_p0_angle, get_sunspyce_carr_rot
#from sunpy.coordinates import spice
#import scipy.io
#from wcs_funs import fitshead2wcs, wcs_get_coord, idlsav2wcs

def get_exp_factor(hdr):
    tel = hdr['detector'].lower()
    mjd = hdr['mid_date']
    jd = mjd + 2400000.5
    t = Time(jd, format='jd')
    dt = t.to_datetime()
    yymm = dt.strftime('%y%m%d') # idl has days so we do too despite the name    
    fn = tel + '_expfactor_'+yymm+'.dat'
    
    # Going rogue from here bc IDL is using common blocks but found the source files
    myDir = efacDir + yymm[:4] + '/'
    efac = 1 # bc if we don't find it, it's nearly one anyway
    bias = None
    if 'offset' in hdr:
        bias = hdr['offset']
    if os.path.exists(myDir+fn):
        data = np.genfromtxt(myDir+fn, dtype=str)
        if hdr['filename'] in data[:,0]:
            idx = np.where(data[:,0] == hdr['filename'])[0]
            efac = float(data[idx[0], 1])
        if type(bias) != type(None):
            bias = float(data[idx[0], 2])            
    return efac, bias

def c2_calfactor(hdr, nosum=False):
    filt = hdr['filter'].upper()
    polar = hdr['POLAR'].upper()
    mjd = hdr['mid_date']
    
    if filt == 'ORANGE':
        cal_factor=4.60403e-07*mjd+0.0374116
        polref=cal_factor/.25256
        if polar == 'CLEAR':
            cal_factor = cal_factor
        else:
            cal_factor = polref
    else:
        sys.exit('Other c2_calfactor cases not ported')
        
    if not nosum:
        if (hdr['sumcol'] > 0): cal_factor = cal_factor / hdr['sumcol']
        if (hdr['sumrow'] > 0): cal_factor = cal_factor / hdr['sumrow']
        if (hdr['lebxsum'] > 0): cal_factor = cal_factor / hdr['lebxsum']
        if (hdr['lebysum'] > 0): cal_factor = cal_factor / hdr['lebysum']
    
    cal_factor = cal_factor * 1e-10    
    return cal_factor
    

def c2_calibrate(imIn, hdr):
    im = np.copy(imIn)
    if hdr['detector'] != 'C2':
        sys.exit('Error in c2_calibrate, passed non C2 files')
    
    # Get the exp_factor
    expfac, bias = get_exp_factor(hdr)
    hdr['EXPTIME'] = hdr['EXPTIME'] * expfac
    hdr['offset']  = bias
    calfac = c2_calfactor(hdr)

    # open the vignetting file
    with fits.open(vigFile) as hdulist:
        vig  = hdulist[0].data
        hdrV = hdulist[0].header
    # Checked and dont need to mask hi/lo
    
    if (hdr['r1col'] != 20) or (hdr['r1row'] != 1) or (hdr['r2col'] != 1043) or (hdr['r2row'] != 1024):
        x1 = hdr['r1col'] - 20
        x2 = hdr['r2col'] - 20
        y1 = hdr['r1row'] - 1
        y2 = hdr['r2row'] - 1
        vig =  vig[y1:y2+1,x1:x2+1]
        print ('Hitting uncheck portion of vignetting in c2_calibrate, should double check')
        
    if (hdr['sumcol'] > 1) or (hdr['sumrow'] > 1):
        # lines 170 -178 in IDL
        sys.exit('Need to rebin vignetting and not done yet, byeeee')
        
    # Ignoring some header history updating
    if hdr['polar'] in ['PB', 'TI', 'UP', 'JY', 'JZ', 'Qs', 'Us', 'Qt', 'Jr', 'Jt']:
        im = im / hdr.exptime
        im = im * calfac
        im = im * vig
    else:
        im = (im - bias) * calfac / hdr['exptime']
        im = im * vig
    
    
    return im, hdr
 
if __name__ == '__main__':
         
    fileA = '/Users/kaycd1/wombat/fits/C2_22800178.fts'
    with fits.open(fileA) as hdulist:
        imA  = hdulist[0].data
        hdrA = hdulist[0].header
    imA, hdrA = c2_calibrate(imA, hdrA)

    fileB ='/Users/kaycd1/wombat/fits/C2_22800186.fts'
    with fits.open(fileB) as hdulist:
        imB  = hdulist[0].data
        hdrB = hdulist[0].header
    imB, hdrB = c2_calibrate(imB, hdrB)
    
    diff = imB - imA

