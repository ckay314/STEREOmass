import numpy as np
import sunpy.map
import sys
import astropy.units as u
from astropy.wcs import WCS


# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)
def calcCMEmass(img, hdr, box=None):
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
    wcs = WCS(hdr)
    
    # IDL uses the crpix a values and they seem to not exist in this
    # version of wcs so pass directly
    crpix = [hdr['crpix1a'], hdr['crpix2a']]
    # IDL seems to use non A versions for this
    cunits = [hdr['CUNIT1'], hdr['CUNIT2']]
    for key in hdr.keys():
        if 'PC' in key:
            print (key, hdr[key])
    dist = wcs_get_coord(wcs, crpix, cunits)
     
    
    # |---------------------------------------|
    # |----------- Apply el Theory -----------|
    # |---------------------------------------|
    
    
    # |---------------------------------------|
    # |--------- Various Conversions ---------|
    # |---------------------------------------|
    
    return
    
def wcs_get_coord(wcs, crpix, cunits):
    # Skipping the distortion part that doesn't get called
    
    # Calculate the 'indicies for each dimension, relative to the ref pixel'  
    naxis = wcs.naxis
    # Assuming we aren't getting the 'PIXEL-LIST' option
    
    # startiing at IDL line 195
    naxis1 = wcs.array_shape[0]
    naxis2 = wcs.array_shape[1]
    num_elements = naxis1 * naxis2
    index = np.arange(num_elements).astype(int)
    coord = np.empty([naxis, num_elements])
    #nn = 1
    #for i in range(2):
    #    coord[i,:] = (index / nn) % wcs.array_shape[i]
    #    nn = nn * wcs.array_shape[i]
    coord[0,:] = index  % naxis1
    coord[1,:] = (index / naxis1 % naxis2)
    
    
    # Skipping apply_dist and jumping to 216
    coord[0,:] = coord[0,:] - (crpix[0] -1)
    coord[1,:] = coord[1,:] - (crpix[1] -1)
    # coord matches at this point, chnages at 263
    # -> # wcs.pc is an operation, not a comment in IDL!
    print (coord[:,0])
    coord = np.matmul(wcs.wcs.pc, coord)
    print (wcs.wcs.pc)
    print (coord[0,0])
    print (sd)
    # Skipping associate and pixel list stuff, down to 356, expect case of TAN
    # Contents of wcs_proj_tan
    halfpi = np.pi / 2.
    cx = np.pi / 180.
    if cunits[0] == 'arcsec':
        cx = cx / 3600.
    # add other cases?
    
    cy = np.pi / 180.
    if cunits[1] == 'arcsec':
        cy = cy / 3600.
    
    # Assuming not helioprojective radial so going to 93
    xmin = np.min(coord[0,:])
    print (coord[0,:])
    

def reclip(aMap, OGmap):
    myCent = [aMap.reference_pixel.x.to_value(), aMap.reference_pixel.y.to_value()]
    OGdim = OGmap.dimensions
    OGx = OGdim.x.to_value()
    OGy = OGdim.y.to_value()
    hwx = OGx / 2 
    hwy = OGy / 2
    myData = aMap.data
    ix1, ix2 = int(myCent[0] - hwx), int(myCent[0] + hwx)
    iy1, iy2 = int(myCent[1] - hwy), int(myCent[1] + hwy)
    reclipData = myData[iy1:iy2, ix1:ix2]
    # Sunpy submap routine only changes crpix/naxis so following that
    aMap.meta['crpix1'] = (aMap.meta['crpix1'] - ix1) 
    aMap.meta['crpix2'] = (aMap.meta['crpix2'] - iy1) 
    aMap.meta['naxis1']  = OGmap.meta['naxis1']
    aMap.meta['naxis2']  = OGmap.meta['naxis2']
    # IF have issues try updating map (bottom_left_coord, dimensions, reference_pixel, top_right_coord)
    # or meta (crval, crota/pc_matrix, crvalA, crpixA, xcen, ycen, pcA )
    reclipMap = sunpy.map.Map(reclipData, aMap.meta)
    return reclipMap
 
def getDiff(aPair):
    my_map1 = sunpy.map.Map(aPair[0])
    my_map2 = sunpy.map.Map(aPair[1])
    
    print (my_map1.fits_header['PC1_1A'])
    print (my_map2.fits_header['PC1_1A'])
    
    if my_map1.dimensions != my_map2.dimensions:
        sys.exit('Dimension mismatch between image and base for '+ aPair[0] + ' and ' + aPair[1])
    
    flData = my_map1.data.astype(np.float32)
        
    my_map1F = sunpy.map.Map(flData, my_map1.meta)
    if 'crota' in my_map1.meta:
        crota = my_map1.meta['crota']
    elif 'crota1' in my_map1.meta:
        crota = my_map1.meta['crota1']
    else:
        crota = 0
     
    if np.abs(crota) > 0.001: 
        # Rotation occurs about crpix (ref pix) which is neither the Sun nor image center
        # After rot sun is centered at crpix-crval/cdelt in fits units (index from 1)
        # and exact match to various sunpy methods of obtaining (accounting for index from 0)
        # If recenter then crpix is in image center (w/indexing diff)
        my_map1FR = my_map1F.rotate(angle=crota*u.deg, missing=0, recenter=True) 
        
        # Not passing an angle to rotate is equiv to rot by angle =crota*u.deg

        # Check the dimensions bc rotation changes it
        if (my_map1.dimensions[0] != my_map1FR.dimensions[0]) or (my_map1.dimensions[1] != my_map1FR.dimensions[1]):
            my_map1FR = reclip(my_map1FR, my_map1)
    else:
        my_map1FR = my_map1F
    

    flData2 = my_map2.data.astype(np.float32)

    if 'crota' in my_map2.meta:
        crota2 = my_map2.meta['crota']
    elif 'crota1' in my_map2.meta:
        crota2 = my_map2.meta['crota1']
    else:
        crota2 = 0
    my_map2F = sunpy.map.Map(flData2, my_map2.meta)
    if np.abs(crota2) > 0.001: 
        my_map2FR = my_map2F.rotate(angle=crota2*u.deg, missing=0, recenter=True)
        if (my_map2.dimensions[0] != my_map2FR.dimensions[0]) or (my_map2.dimensions[1] != my_map2FR.dimensions[1]):
           my_map2FR = reclip(my_map2FR, my_map2)
    else:
        my_map2FR = my_map2F
    
    # check if HI and histogram equalize if so
    if 'HI' in my_map2.meta['detector']:
        
        temp = exposure.equalize_hist(my_map2FR.data - my_map1FR.data)
        rd = sunpy.map.Map(temp, my_map2FR.meta)
       
       
        #fig = plt.figure
        #plt.imshow(temp)
        #plt.show()
        
    else:
        rd = sunpy.map.Map(my_map2FR - my_map1FR.quantity)
    

    return rd

def secchi_prep(fileIn):
    #Port of the most base COR2 version of secchi prep
    print ('secchi preped')  

fileA = '/Users/kaycd1/wombat/fits/20241028_002330_d4c2A.fts'
fileB = '/Users/kaycd1/wombat/fits/20241028_125330_d4c2A.fts'



#diff = getDiff([fileA, fileB])  
#calcCMEmass(diff, diff.fits_header)