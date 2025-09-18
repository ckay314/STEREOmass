import numpy as np
import sunpy.map
import sys, os
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from sunpy.time import parse_time

efacDir = '/Users/kaycd1/ssw/soho/lasco/idl/expfac/data/'
c2vigFile = '/Users/kaycd1/ssw/soho/lasco/idl/data/calib/c2vig_final.fts'
c3previgFile = '/Users/kaycd1/ssw/soho/lasco/idl/data/calib/c3vig_preint_final.fts'
c3postvigFile = '/Users/kaycd1/ssw/soho/lasco/idl/data/calib/c3vig_postint_final.fts'
c3maskFile = '/Users/kaycd1/ssw/soho/lasco/idl/data/calib/c3_cl_mask_lvl1.fts' 
c3rampFile = '/Users/kaycd1/ssw/soho/lasco/idl/data/calib/C3ramp.fts'
c3bkgFile  = '/Users/kaycd1/ssw/soho/lasco/idl/data/calib/3m_clcl_all.fts'

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
        if polar in ['+60DEG', '0DEG', '-60DEG', 'nd']:
            cal_factor = polref
    elif filt == 'BLUE':
        cal_factor = 0.1033
        polref=cal_factor / 0.25256
        if polar in ['+60DEG', '0DEG', '-60DEG', 'nd']:
            cal_factor = polref
    elif filt == 'DEEPRD':
        cal_factor = 0.1033
        polref=cal_factor / 0.25256
        if polar in ['+60DEG', '0DEG', '-60DEG', 'nd']:
            cal_factor = polref
        else:
            cal_factor = 0.
    elif filt in ['HALPHA', 'LENS']:
        cal_factor = 0.1055 # IDL labels this as wrong but doesn't provide a right
        polref=cal_factor / 0.25256
        if polar in ['+60DEG', '0DEG', '-60DEG', 'nd']:
            cal_factor = polref
       
    else:
        sys.exit('Unrecognized filter in c3_calfactor')
        
    if not nosum:
        if (hdr['sumcol'] > 0): cal_factor = cal_factor / hdr['sumcol']
        if (hdr['sumrow'] > 0): cal_factor = cal_factor / hdr['sumrow']
        if (hdr['lebxsum'] > 1): cal_factor = cal_factor / hdr['lebxsum']
        if (hdr['lebysum'] > 1): cal_factor = cal_factor / hdr['lebysum']
    
    cal_factor = cal_factor * 1e-10    
    return cal_factor
    
def c3_calfactor(hdr, nosum=False):
    filt = hdr['filter'].upper()
    polar = hdr['POLAR'].upper()
    mjd = hdr['mid_date']
    
    if filt == 'ORANGE':
        cal_factor = 0.0297
        polref = cal_factor/.25256	
        if polar == '+60DEG':	
            cal_factor=polref
        elif polar =='0DEG':	
            cal_factor=polref * 0.9648
        elif polar == '-60DEG':	
            cal_factor=polref * 1.0798
    elif filt == 'BLUE':
        cal_factor = 0.0975	
        polref = cal_factor / 0.25256
        if polar == '+60DEG':	
            cal_factor=polref
        elif polar =='0DEG':	
            cal_factor=polref * 0.9734
        elif polar == '-60DEG':	
            cal_factor=polref * 1.0613
    elif filt == 'CLEAR':
        cal_factor=7.43e-8 * (mjd - 50000) + 5.96e-3
        polref = cal_factor / 0.25256
        if polar == '+60DEG':	
            cal_factor=polref
        elif polar =='0DEG':	
            cal_factor=polref * 0.9832
        elif polar == '-60DEG':	
            cal_factor=polref * 1.0235
        elif polar == 'H_ALPHA':
            cal_factor = 1.541
    elif filt == 'DEEPRD':
        cal_factor = 0.0259	
        polref = cal_factor / 0.25256
        if polar == '+60DEG':	
            cal_factor=polref
        elif polar =='0DEG':	
            cal_factor=polref * 0.9983
        elif polar == '-60DEG':	
            cal_factor=polref * 1.0300
    elif filt == 'IR':
        cal_factor = 0.0887	
        polref = cal_factor / 0.25256
        if polar == '+60DEG':	
            cal_factor=polref
        elif polar =='0DEG':	
            cal_factor=polref * 0.9833
        elif polar == '-60DEG':	
            cal_factor=polref * 1.0288
    else:
        sys.exit('Unrecognized filter in c3_calfactor')
        
    
    if not nosum:
        if (hdr['sumcol'] > 0): cal_factor = cal_factor / hdr['sumcol']
        if (hdr['sumrow'] > 0): cal_factor = cal_factor / hdr['sumrow']
        if (hdr['lebxsum'] > 1): cal_factor = cal_factor / hdr['lebxsum']
        if (hdr['lebysum'] > 1): cal_factor = cal_factor / hdr['lebysum']
        
    cal_factor = cal_factor*1.e-10
                     
    return cal_factor
    
def c2_calibrate(imIn, hdr, noCalFac=False):
    im = np.copy(imIn)
    if hdr['detector'] != 'C2':
        sys.exit('Error in c2_calibrate, passed non C2 files')
    
    # Get the exp_factor
    expfac, bias = get_exp_factor(hdr)
    hdr['EXPTIME'] = hdr['EXPTIME'] * expfac
    hdr['offset']  = bias
    if not noCalFac:
        calfac = c2_calfactor(hdr)
    else:
        calfac = 1

    # open the vignetting file
    with fits.open(c2vigFile) as hdulist:
        vig  = hdulist[0].data
        hdrV = hdulist[0].header
    # Checked the file used and dont need to mask hi/lo
     
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
        im = im / hdr['exptime']
        im = im * calfac
        im = im * vig
    else:
        im = (im - bias) * calfac / hdr['exptime']
        im = im * vig
    
    return im, hdr
    
def c3_calibrate(imIn, hdr, noCalFac=False, noMask=False):
    im = np.copy(imIn)
    if hdr['detector'] != 'C3':
        sys.exit('Error in c3_calibrate, passed non C3 files')
        
    # Get the exp_factor
    expfac, bias = get_exp_factor(hdr)
    expfac, bias = get_exp_factor(hdr)
    hdr['EXPTIME'] = hdr['EXPTIME'] * expfac
    hdr['offset']  = bias
    
    if not noCalFac:
        calfac = c3_calfactor(hdr)
    else:
        calfac = 1
        
    # Vignetting    
    mjd = hdr['mid_date']
    if mjd < 51000:
        vigFile = c3previgFile
    else:
        vigFile = c3postvigFile
    with fits.open(vigFile) as hdulist:
        vig  = hdulist[0].data
        hdrV = hdulist[0].header
    
    # Mask file    
    with fits.open(c3maskFile) as hdulist:
        mask  = hdulist[0].data
        hdrM = hdulist[0].header
    
    # Ramp file    
    with fits.open(c3rampFile) as hdulist:
        ramp = hdulist[0].data
        hdrR = hdulist[0].header
     
    # Bkg file for fuzziness?    
    with fits.open(c3bkgFile) as hdulist:
        bkg = hdulist[0].data
        hdrb = hdulist[0].header
    
    bkg = 0.8 * bkg / hdrb['exptime']

    if (hdr['r1col'] != 20) or (hdr['r1row'] != 1) or (hdr['r2col'] != 1043) or (hdr['r2row'] != 1024):
        x1 = hdr['r1col'] - 20
        x2 = hdr['r2col'] - 20
        y1 = hdr['r1row'] - 1
        y2 = hdr['r2row'] - 1
        vig =  vig[y1:y2+1,x1:x2+1]
        mask =  mask[y1:y2+1,x1:x2+1]
        ramp =  ramp[y1:y2+1,x1:x2+1]
        bkg =  bkg[y1:y2+1,x1:x2+1]
        
        print ('Hitting uncheck portion of vignetting in c2_calibrate, should double check')
    
    if (hdr['sumcol'] > 1) or (hdr['sumrow'] > 1):
        # lines 265 -289 in IDL
        sys.exit('Need to rebin vignetting and not done yet, byeeee')
    
    # check if monthly image    
    if hdr['fileorig'] == 0:
        print ('Untested monthly image portion of c3_calibrate, should doublecheck')
        im = im / hdr['exptime']
        im = im * calfac * vig - ramp
        if not noMask:
            im = im * mask
        return im, hdr
    
    # Rm the ramp for the colored filters. You know
    if hdr['FILTER'].upper() != 'CLEAR':
        ramp = 0
    
    # check if a polarization brightness image
    if hdr['polar'] in ['PB', 'TI', 'UP', 'JY', 'JZ', 'Qs', 'Us', 'Qt', 'Jr', 'Jt']:
        im = im / hdr['exptime']
        im = im * calfac * vig
        if not noMask:
            im = im * mask
        return im, hdr
    else:
        # Ignoring fuzzy image things
        im = (im - bias) / hdr['exptime']
        im = im * vig * calfac - ramp
        if not noMask:
            im = im * mask
        return im, hdr
    
        
    
    
def c2_prep(filesIn):
    ims, hdrs = [], []
    for aFile in filesIn:
        with fits.open(aFile) as hdulist:
            im  = hdulist[0].data
            hdr = hdulist[0].header
        im, hdr = c2_calibrate(im, hdr)
        ims.append(im)
        hdrs.append(hdr)
    
    return ims, hdrs
    
def c3_prep(filesIn):
     ims, hdrs = [], []
     for aFile in filesIn:
         with fits.open(aFile) as hdulist:
             im  = hdulist[0].data
             hdr = hdulist[0].header
         im, hdr = c3_calibrate(im, hdr)
         ims.append(im)
         hdrs.append(hdr)
    
     return ims, hdrs
    

 
if __name__ == '__main__':
         
    fileA = '/Users/kaycd1/wombat/fits/C2_22800178.fts'
    fileB ='/Users/kaycd1/wombat/fits/C2_22800186.fts'
    
    #fileA ='/Users/kaycd1/wombat/fits/C3_32048310.fts'
    #fileB ='/Users/kaycd1/wombat/fits/C3_32048311.fts'
    
    with fits.open(fileA) as hdulist:
        imA  = hdulist[0].data
        hdrA = hdulist[0].header
    imA, hdrA = c2_calibrate(imA, hdrA)

    with fits.open(fileB) as hdulist:
        imB  = hdulist[0].data
        hdrB = hdulist[0].header
    imB, hdrB = c2_calibrate(imB, hdrB)
    
    diff = imB - imA

