import numpy as np
from sunpy.coordinates import HeliocentricEarthEcliptic, get_body_heliographic_stonyhurst, get_horizons_coord
from sunpy.time import parse_time


# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)


data = np.genfromtxt('haveMassImage.txt', dtype=str)
data2 = np.genfromtxt('coreGCScases.txt', dtype=str)



for i in range(len(data[:,0])):
    eventTime = data[i,0]
    user      = data[i,1]
    fitTime   = data[i,2]
    corAid    = data[i,4]
    CMElon    = data2[i,5]
        
    
    timeA = parse_time(fitTime)
    staLoc = get_horizons_coord('STEREO-A', timeA)
    stbLoc = get_horizons_coord('STEREO-B', timeA)
    lonA = staLoc.lon.degree
    lonB = stbLoc.lon.degree
    
    sepA = np.abs((float(CMElon) % 360) - (lonA % 360))
    if sepA > 180: sepA = 360 - sepA
    possepA = np.abs(90 - sepA)
    
    sepB = np.abs((float(CMElon) % 360) - (lonB % 360))
    if sepB > 180: sepB = 360 - sepB
    possepB = np.abs(90 - sepB)
   
    
    
    print (eventTime, user.rjust(10), fitTime, CMElon.rjust(7), '{:.1f}'.format(lonA).rjust(7), '{:.1f}'.format(sepA).rjust(7),  '{:.1f}'.format(possepA).rjust(7), '{:.1f}'.format(lonB).rjust(7), '{:.1f}'.format(sepB).rjust(7), '{:.1f}'.format(possepB).rjust(7))