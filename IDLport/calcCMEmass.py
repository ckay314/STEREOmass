import numpy as np
import sunpy.map
import sys
import astropy.units as u
#from astropy.wcs import WCS
from secchi_prep import secchi_prep
from wispr_prep import wispr_prep
from wcs_funs import fitshead2wcs, get_Suncent, wcs_get_coord

import matplotlib.pyplot as plt

# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)

dtor = np.pi / 180.

def elTheory(Rin,theta,limb=0.63,center=False, returnAll=False):
    # units are in solar radii and degrees (based on IDL)
    radeg = 180. / np.pi
    # Port of SSWIDL code
    # 0.63 is default for limb darkening
    u = limb
    
    #const = sigma * pi / 2
    #where sigma is the scattering cross section (7.95e-26 per steradian)
    const = 1.24878e-25
    if not center:
        const = const/(1-u/3) # convert to mean solar brightness
    
    # make sure theta is less than 90 deg
    #if theta >= 90:
    #    print ('PoS angle greater than 90, exiting elTheory without running')
    #    return None,  None,  None,  None,  None
    
    R = Rin/np.cos(theta/radeg)
    sinchi2 = (Rin/R)**2	# angle between sun center, electron, observer
    s = 1./R
    #if s >= 1:
    #    print ('Error in computing s, exiting elTheory without running')
    #    return None,  None,  None,  None,  None
    s2 = s**2
    c2 = (1.-s2)
    #if c2 < 0:
    #    print ('Error in computing c2, exiting elTheory without running')
    #    return None,  None,  None,  None,  None	
    c = np.sqrt(c2)			# cos(omega)
    g = c2*(np.log((1.+s)/c))/s
    
    #  Compute Van de Hulst Coefficients
    #  Expressions are given in Billings (1968) after Minnaert (1930)
    ael = c*s2
    cel = (4.-c*(3.+c2))/3.
    bel = -(1.-3.*s2-g*(1.+3.*s2))/8.
    del0 = (5.+s2-g*(5.-s2))/8.
    
    #  Compute electron brightness
    #  pB is polarized brightness (Bt-Br)
    Bt = const*( cel + u*(del0-cel) )
    pB = const* sinchi2 *( (ael + u*(bel-ael) ) )
    Br = Bt-pB
    B = Bt+Br
    Pol = pB/B
    
    if returnAll:
        return R,B,Bt,Br,Pol
    else:
        return R,B


def calcCMEmass(img, hdr, box=None, onlyNe=False, doPB=False):
    # Port of scc_calc_cme_mass from IDL
    # Inputs
    # img = 2d difference image (containing the CME). units should be mean solar brightness
    # hdr = header strcture
    # box = array containing region coordinates (TBD on specifics)
    #     If box = None then just computes mass over full image
    
    # IDL seems to be set up to potentially continue running this with the 
    # distance factors already saved in a common block. This is not the default
    # and for now process doing full on calc each time
    
    # |---------------------------------------|
    # |------ Set up the distance array ------|
    # |---------------------------------------|
    # First part a little different order than IDL but we have code to 
    # get sunc from wcs already
    wcs = fitshead2wcs(hdr) # not system = A here it seems
    
    #sunc =  get_Suncent(wcs)
    #coord = [sunc[0], sunc[1], hdr['crota']*dtor, hdr['rsun']/hdr['cdelt1']]
    
    # Get the distance factor over the full grid so we can use that in eltheory
    dist = wcs_get_coord(wcs) #[2,naxis,naxis] with axis having usual swap from idl
    
    if hdr['cunit1'] == 'deg':
        dist = dist * 3600.
    dist = np.sqrt(dist[0,:,:]**2 + dist[1,:,:]**2) / hdr['rsun']
    
    # |---------------------------------------|
    # |----------- Apply el Theory -----------|
    # |---------------------------------------|
    # Assume no pos keyword or cmelonlat for now (1-6-126)
    pos_angle = 0.
    if not doPB:
        R,B = elTheory(dist,0)
    else:
        R,B,Bt,Br,Pol = elTheory(dist,0, returnAll=True)
    B[np.where(B == 0)] = 1
    
    
    
    # |---------------------------------------|
    # |--------- Various Conversions ---------|
    # |---------------------------------------|
    # A good portion of IDL seems to be commented out so ignoring that
    solar_radius = hdr['rsun']
    cm_per_arcsec = 6.96e10 / solar_radius
    if hdr['cunit1'] == 'deg':
        cm2_per_pixel= (cm_per_arcsec * hdr['cdelt1']*3600.)**2
    else:
        cm2_per_pixel = (cm_per_arcsec * hdr['cdelt1'])**2
    
    # Electron density or mass? 
    if onlyNe:
        conv = cm2_per_pixel
    else:
        conv = cm2_per_pixel * 1.974e-24 # why is this ne2mass separate function in IDL...
    
    if doPB:    
        mass = img / (Bt - Br)
    else:
        mass = img / B
    
    mass = conv * mass     
    
    # Not doing ROI here
    
    # add a tag into header
    hdr['history'] = 'Converted to mass units using calcCMEmass.py'
    
    return mass, hdr
    

# COR2A test files
#fileA = '/Users/kaycd1/wombat/fits/testing/COR2A_20241028_002330_d4c2A.fts'
#fileB = '/Users/kaycd1/wombat/fits/testing/COR2A_20241028_125330_d4c2A.fts'

# COR2B 
fileA = '/Users/kaycd1/wombat/fits/testing/COR2B_20120712_162400_d4c2B.fts'
fileB = '/Users/kaycd1/wombat/fits/testing/COR2B_20120712_185400_d4c2B.fts'

#COR1A
#fileA = '/Users/kaycd1/wombat/fits/testing/COR1A_20121231_001000_n4c1A.fts'
#fileB = '/Users/kaycd1/wombat/fits/testing/COR1A_20121231_011000_n4c1A.fts'

# COR1B
#fileA = '/Users/kaycd1/wombat/fits/testing/COR1B_20111114_001000_n4c1B.fts'
#fileB = '/Users/kaycd1/wombat/fits/testing/COR1B_20111114_011000_n4c1B.fts'

#HI1A
#fileA = '/Users/kaycd1/wombat/fits/testing/HI1A_20100501_000901_s4h1A.fts'
#fileB = '/Users/kaycd1/wombat/fits/testing/HI1A_20100501_164901_s4h1A.fts'

#HI1B
#fileA = '/Users/kaycd1/wombat/fits/testing/HI1B_20070801_012900_s4h1B.fts'
#fileB = '/Users/kaycd1/wombat/fits/testing/HI1B_20070801_164900_s4h1B.fts'

#HI2A
#fileA = '/Users/kaycd1/wombat/fits/testing/HI2A_20090301_000921_s4h2A.fts'
#fileB = '/Users/kaycd1/wombat/fits/testing/HI2A_20090301_180921_s4h2A.fts'

#HI2B
#fileA = '/Users/kaycd1/wombat/fits/testing/HI2B_20130201_000921_s4h2B.fts'
#fileB = '/Users/kaycd1/wombat/fits/testing/HI2B_20130201_100921_s4h2B.fts'


# Python secchi_prep appears to match IDL within 0.005% 
ims, hdrs = secchi_prep([fileA, fileB], outSize=([1024,1024]))
# PSP WISPR
#fileA = '/Users/kaycd1/wombat/fits/testing/psp_L2_wispr_20250610T000025_V0_1221.fits'
#fileB = '/Users/kaycd1/wombat/fits/testing/psp_L2_wispr_20250610T203026_V0_1221.fits'
#ims, hdrs = wispr_prep([fileA, fileB])

diff = ims[1] - ims[0]


mass, hdr = calcCMEmass(diff, hdrs[1])
print(mass[111,444])

fig = plt.figure()
plt.imshow(np.log(np.abs(mass)) * np.sign(mass))
plt.show()