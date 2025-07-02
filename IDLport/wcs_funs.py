import numpy as np
from astropy import wcs

c = np.pi / 180.
cunit2rad = {'arcmin': c / 60.,   'arcsec': c / 3600.,  'mas': c / 3600.e3,  'rad':  1.}


def fitshead2wcs(hdr,system=''):
    # port of IDL version bc astropy is different somehow
    # -> mostly seems that this has ability to pull the A tags
    # for a different 'system'
    
    # skipping hdr check, assuming is fine/tbd
    
    tags = list(hdr.keys())
    
    # Skipping Nan/nan check (250-257)
    # Assuming column is not passed (258-325)
    column = ''
    orig_column = ''
    
    # assume system is A for now bc that's whats passed
    #system = 'A'
    
    # Skipping a lot of the safety checks
    n_axis = hdr['naxis']
    
    compliant = True
    # Determine kind of WCS to extract
    crota_present = False
    pc_present    = False
    cd_present    = False
    for i in range(n_axis):
        if 'CROTA'+str(i+1) in tags: crota_present = True
        for j in range(n_axis):
            if 'PC'+str(i+1)+'_'+str(j+1) in tags: pc_present = True
            # dunno what cd is so skipping for now
    
    if pc_present:
        variation = 'PC'
        if crota_present:
            compliant = False
    elif crota_present:
        variation = 'CROTA'
    else:
        print ('Issue in fitshead2wcs, might be uncoded CD version')
        print (Quit)
        
    
    # Extract CTYPE keywords 
    # Assume that the ctype is found so 1 is x and 2 is y
    if 'CTYPE1' not in tags:
        print ('Issue in fitshead2wcs, missing ctype in header and havent coded alt version')
        print (Quit)
    crp1 = hdr['CRPIX1'] 
    crp2 = hdr['CRPIX2'] 
    crval1 = hdr['CRVAL1'+system]
    crval2 = hdr['CRVAL2'+system]
    cunit1 = hdr['CUNIT1'+system]
    cunit2 = hdr['CUNIT2'+system]
    
    # skipping cname and non compliant strings
    
    if variation == 'PC':
        cdelt1 = hdr['CDELT1'+system]
        cdelt2 = hdr['CDELT2'+system]
        pc = np.zeros([n_axis,n_axis])
        for i in range(n_axis):
            for j in range(n_axis):
                pc[i,j] = hdr['PC'+str(i+1)+'_'+str(j+1)+system]                
    elif variation == 'CROTA':
        print ('Issue in fitshead2wcs, havent coded crota parts')
        print (Quit)
                
    # Determine type of coord sys
    if hdr['CTYPE1'+system][:4] == 'RA--':
        coord_type = 'Celestial-Equatorial'
        projection = hdr['CTYPE1'+system][5:]
    elif  hdr['CTYPE1'+system][:4] == 'HPLN':
        coord_type = 'Helioprojective-Cartesian'
        projection = hdr['CTYPE1'+system][5:]
    else:
        print ('Issue in fitshead2wcs, havent coded non RA coord types')
        print (Quit)
        
    # skipping stuff don't think is needed
    wcsname = coord_type
    
    # Ignoring proj names/vals for now bc unneeded and confusing?
    # the test case has lonpole and pv2_1 at 180 and 0
    
    
    # Make the output dictionary
    wcs = {}
    wcs['coord_type'] = coord_type
    wcs['wcsname'] = wcsname
    wcs['naxis'] = [hdr['naxis1'], hdr['naxis2']]
    wcs['variation'] = variation
    wcs['compliant'] = True
    wcs['projection'] = projection
    wcs['ix'] = 0
    wcs['iy'] = 1
    wcs['crpix'] = [crp1, crp2]
    wcs['crval'] = [crval1, crval2]
    wcs['ctype'] = [hdr['ctype1'+system], hdr['ctype2'+system]]
    wcs['cunit'] = [cunit1, cunit2]
    wcs['cdelt'] = [cdelt1, cdelt2]
    wcs['pc'] = pc
    # wcs['proj_names']
    # wcs['proj_values']
    wcs['roll_angle'] = hdr['SC_ROLL'+system]
    # wc['simple']
    # wcs['time']
    # wcs['position']

    return wcs
    
    
        
        
    
    

def get_Suncent(my_wcs):
    c2rx = cunit2rad[my_wcs['cunit'][0]]
    c2ry = cunit2rad[my_wcs['cunit'][1]]
    coord = [0,0]
    cx = (coord[0] - my_wcs['crval'][0]) / my_wcs['cdelt'][0]
    cy = (coord[1] - my_wcs['crval'][1]) / my_wcs['cdelt'][1]
    pc = my_wcs['pc']
    scx = cx * pc[0,0] + cy * pc[1,0] + my_wcs['crpix'][0] - 1
    scy = cx * pc[0,1] + cy * pc[1,1] + my_wcs['crpix'][1] - 1
    
    return [scx, scy]