import numpy as np
import sunpy.map
import sys
from scc_funs import scc_make_array, scc_zelensky_array, rebinIDL
from cor_prep import cor_prep
from hi_prep import hi_prep
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.time import parse_time
from sunpy.coordinates import get_horizons_coord
from sunspyce import get_sunspyce_hpc_point, get_sunspyce_roll
from sunpy.coordinates import spice
import scipy.io
from wcs_funs import fitshead2wcs, wcs_get_coord



# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)
np.seterr(divide='ignore')

spice.install_frame('IAU_SUN')

# python doesn't know about IDL specific env_vars bc runs in diff env
# so hardcode this for now (NEED TO FIX EVENTUALLY!!!!)
global wcalpath
wcalpath = '/Users/kaycd1/ssw/psp/wispr/calibration/'

global dtor 
dtor = np.pi / 180.


def wispr_readfits(fileIn, LASCO=False, isCal=False):
    # Assume we are passed a good file, open it up and skip
    # down to clean up stuff at 128
    with fits.open(fileIn) as hdulist:
        im  = hdulist[0].data
        hdr = hdulist[0].header
    
    
    # Fix type in early def_wispr_hdr.pro
    hdr['CTYPE1A'] = 'RA---ZPN'
    hdr['CTYPE2A'] = 'DEC--ZPN'
    
    # Fix DATE-* keywords from Bug 460 WISPR FITS DATE-OBS incorrect
    if not isCal: # Our cal files are missing some keywords but don't do this anyway
        if hdr['VERS_CAL'][:2] == '0x':
            sys.exit('Hitting unported code for bug 460. TBD')
        
        # Fill in DATE_OBS in headers
        if not LASCO:
            isOld = False
            if hdr['DATE-OBS'] < '2006-02-14T00:00:00.000':
                isOld = True
        
            if hdr['INSTRUME'] != 'WISPR':
                sys.exit('Hitting unported code for not a WISPR header')
            
    return im, hdr
        
def wispr_bias_offset(im, hdr):
    if hdr['rectify'] == True:
        rectrota = hdr['rectrota']
        # not sure the extent of possible options for this
        # and IDL rotates by (4-rectrota) if it is 1 or 3 
        # test case has 6
        # leave other options uncoded for now
        if rectrota == 1: # -> IDL rotate (im,4-1)
            tempIm = np.rot90(im,k=1)
        elif rectrota == 3: # -> IDL rotate (im,4-3)
            tempIm = np.rot90(im,k=-1)
        elif rectrota == 6:
            tempIm = np.rot90(im[:,::-1], k=-1)
        
        rows = hdr['naxis2']
    else:
        imrot = im
        
    mask = np.array([4,8]) / hdr['nbin1'] - 1
    mask = mask.astype(int)
    subIm = tempIm[mask[0]:mask[1]+1,:]                    
    offset = np.median(subIm[np.isfinite(subIm)])
    
    return offset
    
def wispr_get_calfac(hdr):
    if (hdr['detector'] in [1, 2]) and  (hdr['gainmode'] in ['HIGH', 'LOW']):
        if hdr['detector'] == 1:
            if hdr['gainmode'] == 'LOW': calfac = 2.49e-13
            else:  calfac = 4.09e-14
        else:
            if hdr['gainmode'] == 'LOW': calfac = 3.43e-13
            else:  calfac = 7.28e-14
        if hdr['gaincmd'] == 12: calfac = calfac * 1.27
        hdr['history'] = 'Calibration factor ' + str(calfac) +' applied'                    
    else:
        sys.exit('Invalid WISPR header, issue in calfac')
    return calfac
    
def wispr_get_calimg(hdr):
    if hdr['detector'] not in [1, 2]:
        sys.exit('Invalid WISPR header, issue in calimg')
    if hdr['detector'] == 1:
        calFile = 'WISPR_FredVignettingFitted_inner_20180809_01.fits'
    else:
        calFile = 'WISPR_FredVignettingFitted_outer_20190813_01.fits'    
    hdr['history'] = 'Using '+ calFile
        
    calimg, calhdr = wispr_readfits(wcalpath+calFile, isCal=True)
    
    # Check for sub field
    if (hdr['pxend1'] - hdr['pxbeg1'] != 2047) or (hdr['pxend2'] - hdr['pxbeg2'] != 1919):
        print ('Hitting untested portion, should double check cal img in wispr_get_calimg')
        calimg =  calimg[hdr['pxbeg2']-1:hdr['pxend2'], hdr['pxbeg1']-1:hdr['pxend1']]
    
    # chec if input image has been rectified and roate vignetting to mach
    if hdr['rectify']:
        calimg = np.transpose(calimg)
    
    # If input is binned, bin cali image
    if hdr['nbin'] > 1:
        calimg = rebinIDL(calimg, np.array([hdr['naxis2'], hdr['naxis1']]))

    calimg = calimg ** 2
    w = np.where(calimg == 0)
    if len(w[0]) > 0:
        calimg[w] = np.nan
    return calimg, calhdr

    
def wispr_correction(im, hdr, calfacOff=False, calimgOff=False, exptimeOff=False, truncOff=False):
    hdr['history'] = 'Applied wispr_correction (Python port)'
    
    calcnt = 0 
    
    # Correct for truncation
    if not truncOff and (hdr['IPTRUNC'] != 0):
        im = im * 2**hdr['IPTRUNC']
        hdr['history'] = 'Image multiplied by ' + str(2**hdr['IPTRUNC']) + ' due to truncation'
        
    # Correct for DN/s for binned images
    divfactor = hdr['nbin']
    if divfactor > 1:
        im = im / divfactor
        hdr['history'] = 'Image divided by ' + str(divfactor) + ' for binning'
        hdr['dsatval'] = hdr['dsatval'] / divfactor
            
    # Normalize for exposure time
    if not exptimeOff:
        im = im / hdr['xposure']
        hdr['dsatval'] = hdr['dsatval'] / hdr['xposure']
        hdr['bunit'] = 'DN/s'
        calcnt = calcnt + 1
        hdr['history'] = 'Exposure normalized to 1 second'
    
    # Calfac
    calfac = 1.
    if not calfacOff:
        hdr['bunit'] = 'MSB'
        calfac = wispr_get_calfac(hdr)
        hdr['dsatval'] = hdr['dsatval'] * calfac
        calcnt = calcnt + 1
   
    # cal img
    calimg = 1.
    if not calimgOff:
       calimg, calhdr = wispr_get_calimg(hdr) 
       calcnt = calcnt+1
       
    im = im / calimg * calfac
    
    # calcnt
    if calcnt == 3:
        if hdr['level'] != 'L2':
            hdr['level'] = 'L2'
   
    return im, hdr
    
def get_wispr_pointing(shdr, doSpice=True, doCoords=False):
    shdr['VERS_CAL'] = '2020915'
    detect = shdr['DETECTOR']
    instr = ['SPP_WISPR_INNER', 'SPP_WISPR_OUTER']
    
    # From WISPR_FITSHeaderDefinition_20170811_01.pdf
    #8 WCS Compliant FITS Keyword Definition
    #detector depndent information
    
    # Detector 1
    if detect == 1:
        shdr['CRPIX1'] =  (991.547 - (1920 - shdr['PXEND2'])) / np.sqrt(shdr['nbin'])
        shdr['CRPIX2'] =  (1015.11 - (2048 - shdr['PXEND1'])) / np.sqrt(shdr['nbin'])
        
        if ((shdr['pxend2'] - shdr['pxbeg2'] + 1) / shdr['nbin1'] - shdr['naxis1']) == 32:
            shdr['CRPIX1'] = shdr['CRPIX1'] - 32
        
        shdr['CDELT1']  =	 .0211525 * shdr['nbin1'] # deg
        shdr['CDELT2']  =	 .0211525 * shdr['nbin2'] # deg
        shdr['CDELT1A'] =	-.0211525 * shdr['nbin1'] # deg
        shdr['CDELT2A'] =	 .0211525 * shdr['nbin2'] # deg
        
        shdr['PV2_0'] = -1.9879e-8
        shdr['PV2_1'] = 1.00501
        shdr['PV2_2'] = .0729583
        shdr['PV2_3'] = .275292
        shdr['PV2_4'] = -.701881
        shdr['PV2_5'] = 1.97518

        instrument_roll = 0.417683

        skew = np.array([[.99606, .999191], [.998218,.995773]])

        xcor = 0.54789
        ycor = -.3501
    
    # Detector 2    
    elif detect == 2:
        shdr['CRPIX1'] = (984.733 - (1920 - shdr['PXEND2'])) / shdr['nbin1']
        shdr['CRPIX2'] = (1026.37 - (2048 - shdr['PXEND1'])) / shdr['nbin2']
        
        if (shdr['pxend2'] - shdr['pxbeg2'] + 1) / shdr['nbin1'] - shdr['naxis1'] == 32:
             shdr['CRPIX1'] = shdr['CRPIX1'] - 32
             
        shdr['CDELT1']  =	0.0282376 * shdr['nbin1']  # deg
        shdr['CDELT2']  =	0.0282376 * shdr['nbin2']  # deg
        shdr['CDELT1A'] =  -0.0282376 * shdr['nbin1']  # deg
        shdr['CDELT2A'] =	0.0282376 * shdr['nbin2']  # deg
        
        shdr['PV2_0'] = .000168385
        shdr['PV2_1'] = .983801
        shdr['PV2_2'] = .0737626
        shdr['PV2_3'] = -.374471
        shdr['PV2_4'] = 0.585763
        shdr['PV2_5'] = -.410706

        instrument_roll = 0.40121182
        skew = [[1,1],[1,1]]

        xcor = 0
        ycor = 0
        
    else:
        sys.exit('Invalid detector given to get_wispr_pointing')
        
    shdr['CUNIT1']  =	'deg'
    shdr['CUNIT2']  =	'deg'
    shdr['CUNIT1A'] =	'deg'
    shdr['CUNIT2A'] =	'deg'
    shdr['PV2_0A'] = shdr['PV2_0']
    shdr['PV2_1A'] = shdr['PV2_1']
    shdr['PV2_2A'] = shdr['PV2_2']
    shdr['PV2_3A'] = shdr['PV2_3']
    shdr['PV2_4A'] = shdr['PV2_4']
    shdr['PV2_5A'] = shdr['PV2_5']
    
    point = np.zeros(3)   
    if doSpice:
        point = get_sunspyce_hpc_point(shdr['DATE-AVG'], 'psp', instrument=instr[detect-1], doDeg=True) 
    point[2] = point[2] + instrument_roll
    roll1 = point[2]
    shdr['CRVAL1'] = point[0] + xcor
    shdr['CRVAL2'] = point[1] + ycor
        
    pc = np.array([[np.cos(point[2]*dtor),  -np.sin(point[2]*dtor)], [np.sin(point[2]*dtor), np.cos(point[2]*dtor)]]) * skew
    shdr['PC1_1'] = pc[0,0]
    shdr['PC1_2'] = pc[0,1] # swapped index to match IDL 
    shdr['PC2_1'] = pc[1,0] # swapped index to match IDL 
    shdr['PC2_2'] = pc[1,1] 

    shdr['nbin1'] = np.sqrt(shdr['nbin'])
    shdr['nbin2'] = np.sqrt(shdr['nbin'])
        
    # Skipping get att file for now, IDL doesn't seem to pull anything
    #if ('ATT_FILE' in shdr) & doSpice:
        
    shdr['PV1_1']    = 0.0 	    # deg
    shdr['PV1_2']    = 90.0	    # deg
    shdr['PV1_3']    = 180.0	# deg
    shdr['PV1_1A']   = 0.0 	    # deg
    shdr['PV1_2A']   = 90.0	    # deg
    shdr['PV1_3A']   = 180.0	# deg
    shdr['LATPOLE']  = 0.0
    shdr['LONPOLE']  = 180.0
    shdr['LATPOLEA'] = 0.0
    shdr['LONPOLEA'] = 180.0
    roll = 0.
    ra = 0.
    dec = 90.
        
    if doSpice:
        roll, dec, ra = get_sunspyce_roll(shdr['DATE-AVG'], 'psp', instrument=instr[detect-1], system='GEI')   
    roll = roll + instrument_roll
    if ra < 0: ra = ra +360   
    shdr['CRVAL1A'] = ra+xcor
    shdr['CRVAL2A'] = dec+ycor
    shdr['CRPIX1A'] = shdr['CRPIX1']
    shdr['CRPIX2A'] = shdr['CRPIX2']
    pc= np.array([[np.cos(roll*dtor), np.sin(roll*dtor)], [-np.sin(roll*dtor), np.cos(roll*dtor)]]) * skew
    shdr['PC1_1A'] = pc[0,0]
    shdr['PC1_2A'] = pc[1,0] # swapped index to match IDL 
    shdr['PC2_1A'] = pc[0,1] # swapped index to match IDL 
    
    shdr['PC2_2A'] = pc[1,1]

    # On to detector stuff
    if shdr['detector'] == 1:
        # Load up save files from IDL, may want to replace eventually
        if np.abs(roll) <= 150.:
            idl_save = wcalpath + 'rollcomp.sav'
        elif roll > 150:
            idl_save = wcalpath + 'rollcompp180.sav'
        elif roll < -150:
            idl_save = wcalpath + 'rollcompm180.sav'
        idlRoll = scipy.io.readsav(idl_save, python_dict=True)
        
        idlwcso = scipy.io.readsav(wcalpath+'wcso.sav', python_dict=True)
        wcso =idlwcso['wcso']
        wcso.naxis = wcso.naxis / np.sqrt(shdr['nbin'])
        wcso.crpix = wcso.crpix / np.sqrt(shdr['nbin'])
        wcso.cdelt = wcso['cdelt'] * np.sqrt(shdr['nbin'])
        # Things saved as arrays are [array, dtype]
        wcso.pc[0] = [[np.cos(roll1 * dtor), np.sin(roll1 * dtor)], [-np.sin(roll1 * dtor), np.cos(roll1 * dtor)]]
        wcso.crval[0] = [shdr['crval1'] - xcor, shdr['crval2'] - ycor]
        
        # poly(val, coeffs) IDL -> polyval(coeffs[::-1], val)
        shdr['CRVAL1'] = shdr['CRVAL1'] + np.polyval(idlRoll['p1'][::-1], roll)[0]
        shdr['CRVAL2'] = shdr['CRVAL2'] + np.polyval(idlRoll['p2'][::-1], roll)[0]
        roll1 = roll1 + np.polyval(idlRoll['p3'][::-1], roll)[0]
        pc = np.array([[np.cos(roll1*dtor), -np.sin(roll1*dtor)], [np.sin(roll1*dtor), np.cos(roll1*dtor)]]) * skew
        shdr['PC1_1']=pc[0,0]
        shdr['PC1_2']=pc[0,1] # swapped index to match IDL 
        shdr['PC2_1']=pc[1,0]# swapped index to match IDL 
        shdr['PC2_2']=pc[1,1]
        
        wcs = fitshead2wcs(shdr)
        pt1 = wcs_get_coord(wcs, np.array([wcs['naxis'][0], wcs['naxis'][1]]).astype(int)/2)
        print(sd)
        #pt2 = wcs_get_coord(wcso, [wcs.naxis[0], wcso.naxis[1]]/2)
        #diff = pt1 - pt2
        #print (diff)
        
    coords = None
    return coords
    

def wispr_prep(filesIn, outSize=None, silent=False, biasOff=False, biasOffsetOff=False, lin_correct=False, straylightOff=False, pointingOff=False):
    # Port of basic functionality of IDL version
    
    # Want filesIn as a list, even if single
    if isinstance(filesIn, str):
        filesIn = [filesIn]
        
    # Ignore notupdated, write_flag keywords for now (113-6)
            
    # dont need to pre make arrays bc not idl
    
    num = len(filesIn)
    
    # Ignoreing save path
    
    images_out = []
    headers_out = []
    for i in range(num):
        if not silent:
            print ('Processing image ', i, ' out of ', num)
        
        # Read in the fits aka wispreadfits 
        im, hdr = wispr_readfits(filesIn[i])
        
        # Bias subtraction routine
        if ~biasOff & hdr['ipbias'] ==0:
            sys.exit(' wisrp_ bias needs to be ported')
            # add in bias comment to history
            
        # Remove bias offset
        if ~biasOffsetOff & hdr['ipbias'] !=0:
            offset = wispr_bias_offset(im, hdr)
            im = im - offset
            hdr['history'] = 'Subtracted ' + str(offset) +' for image for bias offest'

        # lin_correct not coded bc not called
        
        # -----------------------
        #    Begin Calibration
        # -----------------------
        im, hdr = wispr_correction(im, hdr)
        
        if (hdr['detector'] == 2) & (not straylightOff):
            print ('Havent ported straylight calib yet')
            
        if int(hdr['level'][-1:]) > 1:
            notupdated = False
            
        full = False    
        if hdr['nbin'] == 1:
            full = True
        if (hdr['detector'] == 1 ) & ((hdr['naxis1'] != 1920/hdr['nbin1']) or (hdr['naxis2'] != 2048/hdr['nbin2'])):
            sys.exit('Need to port putin/zelensky array for wispr 1')
        elif (hdr['detector'] == 2 ) & (hdr['ipmask'] ==0) & ((hdr['naxis1'] != 1920/hdr['nbin1']) or (hdr['naxis2'] != 2048/hdr['nbin2'])):
            sys.exit('Need to port putin/zelensky array for wispr 2')
            
        # Pointing
        if not pointingOff:
            coords = get_wispr_pointing(hdr, doCoords=True)
   
       
            
    return 6, 5