import sunpy.map
import astropy.units as u
from aiapy.calibrate import register, update_pointing
from aiapy.calibrate.util import get_pointing_table
                
def aia_prep(filesIn):
    # No reason to port idl, we have aiapy and wont be doing any WL
    # mass calc so assuming a close enough match
    
    
    
    # |------------------------------------------------------|
    # |------------- Loop to process each image -------------|
    # |------------------------------------------------------|
    num = len(filesIn)
    maps_out = []

    # Assume were working from files and not something loaded from read_sdo
    for i in range(num):
        # aia files are compressed/different from secchi so doensn't work with
        # straight up fits read, but using sunpy map equiv to read_sdo
        aia_map = sunpy.map.Map(filesIn[i])
        
        # Make range wide enough to get closest 3-hour pointing
        pointing_table = get_pointing_table("JSOC", time_range=(aia_map.date - 12 * u.h, aia_map.date + 12 * u.h))
        aia_map_updated_pointing = update_pointing(aia_map, pointing_table=pointing_table)
        
        aia_map_registered = register(aia_map_updated_pointing)
        
        maps_out.append(aia_map_registered)
    return maps_out 
    
fname = ['/Users/kaycd1/wombat/obsFiles/AIA/aia_lev1_304a_2023_03_04t14_00_05_15z_image_lev1.fits']
aiaMap = aia_prep(fname)
print (aiaMap[0].data[1650, 1100])    