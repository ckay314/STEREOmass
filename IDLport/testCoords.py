from wcs_funs import fitshead2wcs, wcs_get_coord, wcs_inv_proj_azp, wcs_inv_proj_tan, wcs_get_pixel
from secchi_prep import secchi_prep
import numpy as np

fileA = '/Users/kaycd1/wombat/fits/testing/HI2A_20090301_000921_s4h2A.fts'
fileA = '/Users/kaycd1/wombat/fits/testing/COR2B_20120712_162400_d4c2B.fts'

'''im, hdr = secchi_prep(fileA) 

myWCS = fitshead2wcs(hdr[0])
coord = wcs_get_coord(myWCS)
aCoord = np.array([coord[0,400:402,230],coord[1,400:402,230]])
print(aCoord)
print ("")
#coord = wcs_inv_proj_azp(myWCS, aCoord)
#print (aCoord)

coord = wcs_inv_proj_tan(myWCS, aCoord)
print (aCoord)'''


# Convert from heliocentric radial to heliocentric cartesian




fileA = '/Users/kaycd1/wombat/fits/20120712_172400_d4c2A.fts'
im, hdr = secchi_prep(fileA) 
print (hdr)
#print(hdr)
'''myWCS = fitshead2wcs(hdr[0])
coord = np.array([0,0])
coord = wcs_get_pixel(myWCS, coord)
print (coord)
coord = np.array([995.83,0])
coord = wcs_get_pixel(myWCS, coord)
print (coord)
coord = np.array([0,995.83])
coord = wcs_get_pixel(myWCS, coord)'''

myWCS = fitshead2wcs(hdr[0])
pt = [-10873.3709868, -260.11013321]
coord = wcs_get_pixel(myWCS, pt, doQuick = False)

print (coord)
