import numpy as np
import os
import sys
from astropy.io import fits
import datetime
from scipy.interpolate import griddata
from scc_funs import scc_sebip, scc_hi_diffuse
from cor_prep import get_calfac, get_calimg

#from scc_funs import secchi_rectify, fill_from_defhdr, rebinIDL, scc_getbkgimg
#from wcs_funs import get_Suncent, fitshead2wcs

def sc_inverse(n,diag, below, above):
    wt_above = above / diag
    wt_below = below / diag
    
    wt_above1 = wt_above - 1
    wt_below1 = wt_below - 1
    
    ints = np.array(range(n-2))+1
    power_above = np.zeros(n-1)
    power_below = np.zeros(n-1)
    
    power_above[0] = 1
    power_below[0] = 1
    
    # Can do w/o for loop in python
    power_above[1:] = wt_above1 ** ints
    power_below[1:] = wt_below1 ** ints
    
    v, u = np.zeros(n), np.zeros(n)
    v[1:] = wt_below * (power_below * power_above[::-1])
    u[1:] = wt_above * (power_above * power_below[::-1])

    d = -u[1] / wt_above - (np.sum(v)-v[n-1])
    f = 1. / (diag * (d + wt_above*np.sum(v)))
    
    u[0], v[0] = d, d
    u = u*f
    v = v[::-1] * f
    p = np.empty([n,n])
    # set up p the same as IDL, transverse at tend
    # IDL indexing is very non-pythonic, think just give it starting index?
    p[:,0] = u[0] 
    for row in np.array(range(n-3))+1:
        p[:row,row] = v[n-row-1:n-1]
        p[row:,row] = u[0:n-row]
    p[:,-1] = v

    # now in python indexing which is oppo of IDL
    p = np.transpose(p)

    return p

def hi_remove_saturation(im, hdr, saturation_limit=None, nsaturated=None):
    # ignoring header check
    if saturation_limit == None: saturation_limit = 14000
    if saturation_limit < 0:
        return im, hdr
    if nsaturated == None: nsaturated = 5
    
    n_im = hdr['imgseq'] + 1
    ssum = hdr['summed']
    dsatval = saturation_limit * n_im*(2.**(ssum-1))**2
    
    ii = np.where(im > dsatval)[0]
    nii = len(ii)
    
    if nii > 0:
        print('Hitting coded portion of hi_remove_saturation')
        print (Quit)
        mask = np.copy(im) * 0
        mask[ii] = 1
    else:
        return im, hdr

def hi_cosmics(im, hdr):
    # Assuming ok header
    if ('s4h' not in hdr['filename']):
        cosmics = hdr['cosmics']
    elif (hdr['n_images'] < 1) & (hdr['imgseq'] < 1):
        cosmics = hdr['cosmics']
    else:
        count = hdr['imgseq'] + 1
    
        inverted = False
        if hdr['rectify']:
            if (hdr['date-obs'] > '2015-07-01T00:00:00') & (hdr['date-obs'] < '2023-08-12T00:00:00'):
                inverted = hdr['OBSRVTRY'] == 'STEREO_A'
            else:
                inverted = hdr['OBSRVTRY'] == 'STEREO_B'
                
        # Inverted Case
        if inverted:
            print ('invert case in hi_cosmics not tested, should double check')
            cosmic_counter = im[0,count]
            if cosmic_counter == count:
                cosmics = im[0,:count][::-1]
            else:
                print ('hit un ported part in hi_cosmics, need to do')
                print (Quit)
        
        # Non inverted case
        else:
            naxis1 = hdr['naxis1']
            naxis2 = hdr['naxis2']
            # gotta switch axes for python
            cosmic_counter = im[naxis2-1, naxis1-count-1]
            if cosmic_counter == count:
                cosmics = np.copy(im[naxis2-1,naxis1-count:naxis1])
                # fill image from row below
                im[naxis2-1,naxis1-count-1:naxis1]=im[naxis2-2, naxis1-count-1:naxis1]
            else:
                print ('hit un ported part in hi_cosmics, need to do')
                print (Quit)
    return cosmics

def hi_desmear(im,hdr):
    # Check header for valid values
    if hdr['CLEARTIM'] < 0:
        sys.exit('CLEARTIM invalid in header. Cannot desmear')
    if hdr['RO_DELAY'] < 0:
        sys.exit('RO_DELAY invalid in header. Cannot desmear')
    if hdr['LINE_CLR'] < 0:
        sys.exit('LINE_CLR invalid in header. Cannot desmear')
    if hdr['LINE_RO'] < 0:
        sys.exit('LINE_RO invalid in header. Cannot desmear')
        
    date_obs = hdr['date-obs']    
    post_conj = False
    if (date_obs > '2015-07-01T00:00:00') & (date_obs < '2023-08-12T00:00:00'):
        post_conj = True
    
    # Extract image array if underscan present.
    if (hdr['dstart1'] < 1) or (hdr['naxis1'] == hdr['naxis2']):
        image = im
    else:
        print ('Unchecked line in hi_desmear')
        image = im[hdr['dstart2']-1:hdr['dstop2'],hdr['dstart1']-1:hdr['dstop1']]
    
    # Compute the inverted correction matrix
    clearest=0.70
    exp_eff = hdr['EXPTIME'] + hdr['n_images'] * (clearest-hdr['CLEARTIM'] + hdr['RO_DELAY'])
    
    # Weight correction by number of images and some other words
    dataWeight = hdr['n_images'] * ((2**(hdr['ipsum']-1)))

    inverted = False
    if hdr['rectify']:
        if hdr['OBSRVTRY'] == 'STEREO_B': inverted = True
        if (hdr['OBSRVTRY'] == 'STEREO_A') & (post_conj): inverted = True
    
    if inverted:
        print ('Unchecked line in hi_desmear')
        fixup = sc_inverse(hdr['naxis2'], exp_eff, dataWeight*hdr['line_clr'], dataWeight*hdr['line_ro'])
    else:
        fixup = sc_inverse(hdr['naxis2'], exp_eff, dataWeight*hdr['line_ro'], dataWeight*hdr['line_clr'])
    
    # matrix mult
    image = np.matmul(fixup, image)
 
    # patch the repaired image back in if needed
    if (hdr['dstart1'] < 1) or (hdr['naxis1'] == hdr['naxis2']):
        im = image
    else:
        im[hdr['dstart2']-1:hdr['dstop2'],hdr['dstart1']-1:hdr['dstop1']] = image
        
    return im, hdr

def hi_correction(im, hdr, sebip_off=False, bias_off=False, exptime_off=False, desmear_off=False, calfac_off=False, calimg_off=False):
    # Assuming valid header structure
    
    # Correct for SEB IP 
    if not sebip_off:
        im, hdr, sebipFlag  = scc_sebip(im, hdr)

    # Bias subtraction
    if bias_off:
        biasmean = 0.
    else:
        biasmean = float(hdr['biasmean'])
        # Check if done onboard
        ip19str =hdr['IP_00_19']
        if ('103' in ip19str) or ('37' in ip19str) or ('38' in ip19str):
            biasmean = 0
        if hdr['ipsum'] > 1:
            biasmean = biasmean * (2** (hdr['ipsum']-1))**2
        if biasmean != 0.:
            hdr['history'] = 'Bias subtracted '+ str(biasmean)
            hdr['OFFSETCR'] = biasmean
            im = im - biasmean
            
    # Extract and correct for cosmics
    cosmics = hi_cosmics(im, hdr)    
    im, hdr = hi_remove_saturation(im, hdr)
    
    # Exposure time
    if not exptime_off:
        if desmear_off:
            print('Need to code desmear off in hi_correction in hi_prep')
            print (Quit)
        else:
            im, hdr = hi_desmear(im, hdr)
            if hdr['nmissing'] > 0:
                print ('Need to code hi_fill_missing')
                print (Quit)
            hdr['bunit'] = 'DN/s'
    
    # capture ipsum
    ipkeep = hdr['ipsum']
    
    # Apply calibration factor, cannot be done before hi_desmear
    if calfac_off:
        calfac = 1.
    else:
        calfac = get_calfac(hdr)
    
    diffuse = 1.0
    if (calfac != 1):
        if not calimg_off:
            hdr['history'] = 'Applied calibration factor'
            diffuse = scc_hi_diffuse(hdr, ipsum=ipkeep)
            hdr['history'] = 'Applied diffuse source correction'
    else:
        calfac_off = True
        
    # Correction for flat field and vignetting
    if calimg_off:
        calimg = 1.0
    else:
        calimg = get_calimg(hdr)
        print(calimg[1069,420])
    
    
    sys.exit('Reached here')
    return im, hdr

def hi_prep(im, hdr):
    # Assuming proper header structure again
    det = hdr['DETECTOR']
    if det not in ['HI1', 'HI2']:
        sys.exit('hi_prep only for HI detector')
    
    # fixing 'bugzilla 332'    
    if (hdr['naxis1'] > 1024) & (hdr['imgseq'] !=0) & (hdr['n_images'] == 1):
        hdr['imgseq'] = 0
    
    # Calibration/correction
    # cosmic seems undefined on first pass then saved after for future calls?
    im, hdr = hi_correction(im, hdr)
    
    sys.exit('Reached here')
    return im, hdr