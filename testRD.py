import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from astropy.io import fits


baseDir = '/Users/kaycd1/IDLWorkspace/temp/'

imA0F = '20120712_162400_d4c2A.fts'
imB0F = '20120712_162400_d4c2B.fts'

imAF = '20120712_175400_d4c2A.fts'
imBF = '20120712_175400_d4c2B.fts'


hdul = fits.open(baseDir+imA0F) 
imA0 = hdul[0].data
hdul.close()

hdul = fits.open(baseDir+imB0F) 
imB0 = hdul[0].data
hdul.close()


hdul = fits.open(baseDir+imAF) 
imA = hdul[0].data
hdul.close()

hdul = fits.open(baseDir+imBF) 
imB = hdul[0].data
hdul.close()


dif = imB - imB0

fig = plt.figure()
plt.imshow(dif, origin='lower', cmap='magma')
plt.show()