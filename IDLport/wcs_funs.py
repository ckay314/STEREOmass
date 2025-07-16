import numpy as np
#from astropy import wcs

c = np.pi / 180.
cunit2rad = {'arcmin': c / 60.,   'arcsec': c / 3600.,  'mas': c / 3600.e3,  'rad':  1., 'deg':c}


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
            # add lowercase to work with sunpy map metadata format
            elif 'pc'+str(i+1)+'_'+str(j+1) in tags: 
                pc_present = True
                hdr['PC'+str(i+1)+'_'+str(j+1)] = hdr['pc'+str(i+1)+'_'+str(j+1)]
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
    if 'CTYPE1' in tags:
        crp1 = hdr['CRPIX1'] 
        crp2 = hdr['CRPIX2'] 
        crval1 = hdr['CRVAL1'+system]
        crval2 = hdr['CRVAL2'+system]
        cunit1 = hdr['CUNIT1'+system]
        cunit2 = hdr['CUNIT2'+system]
    elif 'ctype1' in tags:
        crp1 = hdr['crpix1'] 
        crp2 = hdr['crpix2'] 
        crval1 = hdr['crval1'+system]
        crval2 = hdr['crval2'+system]
        cunit1 = hdr['cunit1'+system]
        cunit2 = hdr['cunit2'+system]
    else:
        print ('Issue in fitshead2wcs, missing ctype in header and havent coded alt version')
        print (Quit)
    
    
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
    
    # Proj names/vals 
    proj_names = []
    proj_values = []
    if 'LONPOLE' in hdr:
        proj_names.append('LONPOLE')
        proj_values.append(hdr['LONPOLE'])
    if 'LATPOLE' in hdr:
        proj_names.append('LATPOLE')
        proj_values.append(hdr['LATPOLE'])
    for key in tags:
        if 'PV' in key:
            if key[-1] != 'A': # don't seem to want the A tag version, not certain here tho
                proj_names.append(key)
                proj_values.append(hdr[key])
    
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
    wcs['proj_names'] = np.array(proj_names)
    wcs['proj_values'] = np.array(proj_values)
    wcs['roll_angle'] = hdr['SC_ROLL'+system]
    # wc['simple']
    # wcs['time']
    # wcs['position']

    return wcs
                       
def get_Suncent(my_wcs):
    # this is simple version of wcs_get_coord with coord [0,0]
    # shouldn't use for non TAN projection
    c2rx = cunit2rad[my_wcs['cunit'][0]]
    c2ry = cunit2rad[my_wcs['cunit'][1]]
    coord = [0,0]
    cx = (coord[0] - my_wcs['crval'][0]) / my_wcs['cdelt'][0]
    cy = (coord[1] - my_wcs['crval'][1]) / my_wcs['cdelt'][1]
    pc = my_wcs['pc']
    scx = cx * pc[0,0] + cy * pc[1,0] + my_wcs['crpix'][0] - 1
    scy = cx * pc[0,1] + cy * pc[1,1] + my_wcs['crpix'][1] - 1
    
    return [scx, scy]
    
def wcs_proj_tan(my_wcs, coord, doQuick=False, force_proj=False):
    # Check shape of input array
    singlePt = False
    if len(coord.shape) == 1:
        singlePt = True
        coord = coord.reshape([2,1])
        
    dtor = np.pi / 180.
    halfpi = np.pi / 2
    cx = cunit2rad[my_wcs['cunit'][0]]
    cy = cunit2rad[my_wcs['cunit'][1]]

    # Check if within 3 deg of sun, if so switch to quick proj
    # Setting force_proj to True will overwrite this
    if not force_proj:
        if my_wcs['coord_type'] == 'Helioprojective-Radial':
            ymm = np.array([np.min(coord[1,:]), np.max(coord[1,:])])
            yrange = (ymm + my_wcs['crval'][1]) * cy + halfpi
            if np.max(np.abs(yrange)) <= 3 * dtor: 
                doQuick = True
        else:
            xmm = np.array([np.min(coord[0,:]), np.max(coord[0,:])])
            xxrange = (xmm + my_wcs['crval'][0]) * cx
            ymm = np.array([np.min(coord[1,:]), np.max(coord[1,:])])
            yrange = (ymm + my_wcs['crval'][1]) * cy 
            if (np.max([np.max(np.abs(xxrange)), np.max(np.abs(yrange))])) <= 3 * dtor:
                doQuick = True

    # Quick version
    if doQuick and not force_proj:
        if my_wcs['coord_type'] == 'Helioprojective-Radial':
            x = coord[0,:] * cx
            y = (coord[1,:] + my_wcs['crval'][1]) * cy + halfpi
            coord[0,:] = my_wcs['crval'][0] + np.arctan2(x,y) / cx
            coord[1,:] = (np.sqrt(x**2 + y**2) - halfpi) / cy
            return coord
        else:
            coord[0,:] = coord[0,:] #+ my_wcs['crval'][0]         
            coord[1,:] = coord[1,:] #+ my_wcs['crval'][1]
            return coord
        
    # Full version
    # assume standard phi0, theta0 for now
    phi0, theta0 = 0. * dtor, 90. * dtor
    # Get the celestial longitude and latitude of the fiducial point.
    alpha0 = my_wcs['crval'][0] * cx
    delta0 = my_wcs['crval'][1] * cy
    
    # Get the native long of the celestial pole
    # assuming default values for now
    if delta0 > theta0:
        phip = 0
    else:
        phip = 180 * dtor
        
    # Calculate native spherical coords
    phi = np.arctan2(cx*coord[0,:], -cy*coord[1,:])
    theta = np.sqrt((cx*coord[0,:])**2 + (cy*coord[1,:])**2)
    # Correct the theta
    w0 = np.where(theta == 0)[0]
    w1 = np.where(theta != 0)[0]
    theta[w0] = halfpi
    theta[w1] = np.arctan(1 / theta[w1])

    # Calculate the celestial spherical coordinates
    if delta0 >= halfpi:
        alpha = alpha0 + phi - phip - np.pi
        delta = theta
    elif delta0 <= -halfpi:
        alpha = alpha0 - phi + phip
        delta = -theta
    else:
        dphi = phi - phip
        cos_dphi = np.cos(dphi)
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        alpha = alpha0 + np.arctan2(-cos_theta*np.sin(dphi), sin_theta*np.cos(delta0) - cos_theta*np.sin(delta0)*cos_dphi)
        delta = np.arcsin(sin_theta*np.sin(delta0) + cos_theta * np.cos(delta0) * cos_dphi)
        
    # Convert back to og units
    coord[0,:] = alpha / cx
    coord[1,:] = delta / cy
    
    if singlePt:
        coord = coord.flatten()
    
    return coord

def wcs_inv_proj_tan(my_wcs, coord, doQuick=False, force_proj=False):
    # Check shape of input array
    singlePt = False
    if len(coord.shape) == 1:
        singlePt = True
        coord = coord.reshape([2,1])
        
    dtor = np.pi / 180.
    halfpi = np.pi / 2
    cx = cunit2rad[my_wcs['cunit'][0]]
    cy = cunit2rad[my_wcs['cunit'][1]]
    
    # Quick version
    if not force_proj:
        if my_wcs['coord_type'] == 'Helioprojective-Radial':
            ymm = np.array([np.min(coord[1,:]), np.max(coord[1,:])])
            yrange = (ymm + my_wcs['crval'][1]) * cy + halfpi
            if np.max(np.abs(yrange)) <= 3 * dtor: 
                doQuick = True
        else:
            xmm = np.array([np.min(coord[0,:]), np.max(coord[0,:])])
            xxrange = (xmm + my_wcs['crval'][0]) * cx
            ymm = np.array([np.min(coord[1,:]), np.max(coord[1,:])])
            yrange = (ymm + my_wcs['crval'][1]) * cy 
            if (np.max([np.max(np.abs(xxrange)), np.max(np.abs(yrange))])) <= 3 * dtor:
                doQuick = True
    
    # Quick version
    if doQuick and not force_proj:
        if my_wcs['coord_type'] == 'Helioprojective-Radial':
            r = coord[1,:] * cy + halfpi
            theta = (coord[0,:] - wcs['crval'][0]) * cx
            coord[0,:] = r * np.sin(theta) / cx
            coord[1,:] = r * np.cos(theta) / cy            
            return coord
        else:
            coord[0,:] = coord[0,:] - my_wcs['crval'][0]         
            coord[1,:] = coord[1,:] - my_wcs['crval'][1]
            return coord
            
    phi0 = 0.
    theta0 = 90.
    if 'proj_names' in my_wcs:
        for item in my_wcs['proj_names']:
            if item == 'PV1_1':
                idx = np.where(my_wcs['proj_names'] == item)[0]
                phi0 = my_wcs['proj_values'][idx[0]]
            if item == 'PV1_2':
                idx = np.where(my_wcs['proj_names'] == item)[0]
                phi0 = my_wcs['proj_values'][idx[0]]
    if (phi0 != 0) or (theta0 != 90):
        print ('Non-standard PVi_1 and/or PVi_2 values -- ignored') # so why did we bother checking?
    
    # Convert to rads
    phi0, theta0 = phi0 * dtor, theta0 * dtor

    # Get the celestial longitude and latitude of the fiducial point.
    alpha0, delta0 = my_wcs['crval'][0] * cx , my_wcs['crval'][1] * cy
    
    phip = 180
    if delta0 > theta0: phip = 0
    if 'proj_names' in my_wcs:
        if 'LONPOLE' in my_wcs['proj_names']:
            phip = my_wcs['proj_values'][np.where(my_wcs['proj_names'] == 'LONPOLE')][0]
        if 'PV1_3' in my_wcs['proj_names']:
            phip = my_wcs['proj_values'][np.where(my_wcs['proj_names'] == 'PV1_3')][0]
    if (phip != 180) & (delta0 != halfpi):
        print('Non standard LONPOLE value')
    
    phip = dtor * phip
    
    # Convert from celestial to native spherical coordinates.
    alpha = cx * coord[0,:]
    delta = cy * coord[1,:]
    dalpha = alpha - alpha0
    cos_dalpha = np.cos(dalpha)
    sin_delta = np.sin(delta)
    cos_delta = np.cos(delta)
    phi = phip + np.arctan2(-cos_delta * np.sin(dalpha), sin_delta * np.cos(delta0) - cos_delta * np.sin(delta0) * cos_dalpha)
    theta = np.arcsin(sin_delta * np.sin(delta0) + cos_delta * np.cos(delta0) * cos_dalpha)
    # Calculate the relative coords
    r_theta = np.copy(theta)
    w_good = np.where(r_theta > 0)
    r_theta[w_good] = 1. / np.tan(theta[w_good])
    x = r_theta * np.sin(phi)
    y = -r_theta * np.cos(phi)
    
    # Convert back to og units
    coord[0,:] = x / cx    
    coord[1,:] = y / cy
    
    if singlePt:
        coord = coord.flatten()
    
    
    return coord
    
def wcs_proj_azp(my_wcs, coord):
    dtor = np.pi / 180.
    halfpi = np.pi / 2
    cx = cunit2rad[my_wcs['cunit'][0]]
    cy = cunit2rad[my_wcs['cunit'][1]]
    phi0 = 0.
    theta0 = 90.
    if 'proj_names' in my_wcs:
        for item in my_wcs['proj_names']:
            if item == 'PV1_1':
                idx = np.where(my_wcs['proj_names'] == item)[0]
                phi0 = my_wcs['proj_values'][idx[0]]
            if item == 'PV1_2':
                idx = np.where(my_wcs['proj_names'] == item)[0]
                phi0 = my_wcs['proj_values'][idx[0]]
    if (phi0 != 0) or (theta0 != 90):
        print ('Non-standard PVi_1 and/or PVi_2 values -- ignored') # so why did we bother checking?
    
    # Convert to rads
    phi0, theta0 = phi0 * dtor, theta0 * dtor
    
    mu = 0.
    gamma = 0.
    if 'proj_names' in my_wcs:
        for item in my_wcs['proj_names']:
            if item == 'PV2_1':
                idx = np.where(my_wcs['proj_names'] == item)[0]
                mu = my_wcs['proj_values'][idx[0]]
            if item == 'PV2_2':
                idx = np.where(my_wcs['proj_names'] == item)[0]
                gamma = my_wcs['proj_values'][idx[0]]
    
    # Convert gamma to radians
    gamma = gamma  * dtor
    
    # Get the celestial longitude and latitude of the fiducial point.
    alpha0, delta0 = my_wcs['crval'][0] * cx , my_wcs['crval'][1] * cy
    
    phip = 180
    if delta0 > theta0: phip = 0
    if 'proj_names' in my_wcs:
        if 'LONPOLE' in my_wcs['proj_names']:
            phip = my_wcs['proj_values'][np.where(my_wcs['proj_names'] == 'LONPOLE')][0]
        if 'PV1_3' in my_wcs['proj_names']:
            phip = my_wcs['proj_values'][np.where(my_wcs['proj_names'] == 'PV1_3')][0]
    if (phip != 180) & (delta0 != halfpi):
        print('Non standard LONPOLE value')
    
    phip = dtor * phip
    
    # Calculate the native spherical coords
    phi = np.arctan2(cx*coord[0,:], -cy*coord[1,:])
    r_theta = np.sqrt((cx*coord[0,:])**2 + (cy*coord[1,:]*np.cos(gamma))**2)
    if gamma == 0:
        rho = r_theta / (mu + 1)
    else:
        rho = r_theta / (mu + 1 + cy*coord[1,:]*np.sin(gamma))
    psi = np.atan(1/rho)
    omega = rho * mu / np.sqrt(rho**2 + 1)
    # check for values outside +/- 1
    badIdx = np.where(np.abs(omega) > 1)
    # just replace with signed 1, not sure exactly what IDL does
    omega[badIdx] = np.sign(omega[badIdx])
    omega = np.arcsin(omega)
    theta = psi - omega   
    badIdx = np.where(theta > np.pi)
    theta[badIdx] = theta[badIdx] - 2*np.pi
    theta2 = psi + omega + np.pi
    badIdx = np.where(theta2 > np.pi)
    theta2[badIdx] = theta2[badIdx] - 2*np.pi
    badIdx = np.where((np.abs(theta2-halfpi)< np.abs(theta-halfpi)) & (np.abs(theta2 < halfpi)))
    badIdx2 = np.where(np.abs(theta) > halfpi)
    theta[badIdx] = theta2[badIdx]
    theta[badIdx2] = theta2[badIdx2]
    
    
    # Calculate the celestial spherical coordinates
    if delta0 > halfpi:
        alpha = alpha0 + phi - phip - np.pi
        delta = theta
    elif delta0 < -halfpi:
        alpha = alpha0 - phi + phip
        delta = -theta
    else:
        dphi = phi - phip
        cos_dphi = np.cos(dphi)
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        alpha = alpha0 + np.arctan(-cos_theta * np.sin(dphi) / (sin_theta * np.cos(delta0) - cos_theta * np.sin(delta0) * cos_dphi))
        delta = np.arcsin(sin_theta * np.sin(delta0) + cos_theta * np.cos(delta0) * cos_dphi)
    
    # Convert back into og units    
    coord[0,:] = alpha / cx
    coord[1,:] = delta / cy
    return coord

def wcs_inv_proj_azp(my_wcs, coord):
    # Check shape of input array
    singlePt = False
    if len(coord.shape) == 1:
        singlePt = True
        coord = coord.reshape([2,1])
    
    
    dtor = np.pi / 180.
    halfpi = np.pi / 2
    cx = cunit2rad[my_wcs['cunit'][0]]
    cy = cunit2rad[my_wcs['cunit'][1]]
    
    phi0 = 0.
    theta0 = 90.
    if 'proj_names' in my_wcs:
        for item in my_wcs['proj_names']:
            if item == 'PV1_1':
                idx = np.where(my_wcs['proj_names'] == item)[0]
                phi0 = my_wcs['proj_values'][idx[0]]
            if item == 'PV1_2':
                idx = np.where(my_wcs['proj_names'] == item)[0]
                phi0 = my_wcs['proj_values'][idx[0]]
    if (phi0 != 0) or (theta0 != 90):
        print ('Non-standard PVi_1 and/or PVi_2 values -- ignored') # so why did we bother checking?
    
    # Convert to rads
    phi0, theta0 = phi0 * dtor, theta0 * dtor
    
    mu = 0.
    gamma = 0.
    if 'proj_names' in my_wcs:
        for item in my_wcs['proj_names']:
            if item == 'PV2_1':
                idx = np.where(my_wcs['proj_names'] == item)[0]
                mu = my_wcs['proj_values'][idx[0]]
            if item == 'PV2_2':
                idx = np.where(my_wcs['proj_names'] == item)[0]
                gamma = my_wcs['proj_values'][idx[0]]
    
    # Convert gamma to radians
    gamma = gamma  * dtor
    
    # Get the celestial longitude and latitude of the fiducial point.
    alpha0, delta0 = my_wcs['crval'][0] * cx , my_wcs['crval'][1] * cy
    
    phip = 180
    if delta0 > theta0: phip = 0
    if 'proj_names' in my_wcs:
        if 'LONPOLE' in my_wcs['proj_names']:
            phip = my_wcs['proj_values'][np.where(my_wcs['proj_names'] == 'LONPOLE')][0]
        if 'PV1_3' in my_wcs['proj_names']:
            phip = my_wcs['proj_values'][np.where(my_wcs['proj_names'] == 'PV1_3')][0]
    if (phip != 180) & (delta0 != halfpi):
        print('Non standard LONPOLE value')
    
    phip = dtor * phip

    # Convert from celestial to native spherical coords
    alpha = cx * coord[0,:]
    delta = cy * coord[1,:]
    dalpha = alpha - alpha0
    cos_dalpha = np.cos(dalpha)
    sin_delta = np.sin(delta)
    cos_delta = np.cos(delta)
    phi = phip + np.arctan2(-cos_delta * np.sin(dalpha), sin_delta * np.cos(delta0) - cos_delta * np.sin(delta0) * cos_dalpha)
    theta = np.arcsin(sin_delta * np.sin(delta0) + cos_delta * np.cos(delta0) * cos_dalpha)
    
    # Determine the latitude of divergence
    if mu == 0:
        thetax = 0
    elif np.abs(mu) > 1:
        thetax = np.arcsin(-1/mu)
    else:
        thetax = np.arcsin(-mu)
    
    # Calculate the relative coords
    cos_theta = np.cos(theta)
    if gamma == 0:
        denom = mu + np.sin(theta)
    else:
        denom = mu + np.sin(theta) + cos_theta * np.cos(phi) * np.tan(gamma)
    
    w_good = np.where(denom !=0) or np.where(theta <=thetax) # inverted from IDL bc want good
    if len(w_good[0]) > 0:
        theta[w_good[0]] = (mu + 1) * cos_theta[w_good] / denom[w_good]
    x = theta * np.sin(phi)
    
    if gamma == 0:
        y = -theta * np.cos(phi)
    else:
        y = -theta * np.cos(phi) / np.cos(gamma)
        
    # Convert back to og units
    coord[0,:] = x / cx    
    coord[1,:] = y / cy
    
    if singlePt:
        coord = coord.flatten()
    
    return coord

def wcs_get_coord(my_wcs):
    # Assuming an appropriate header
    # ignoring distortion, associate, apply for now (139-152)
    
    # Get the dimensions/indices
    naxis = my_wcs['naxis']
    n_axis = len(naxis) # these var names bother me but following idl
    naxis1 = naxis[0]
    naxis2 = naxis[1]
    # assuming not pixel list (158-195)
    num_elements = np.prod(naxis)
    index = np.arange(num_elements).astype(int)
    coord = np.empty([n_axis, num_elements])
    coord[0,:] = index  % naxis1
    coord[1,:] = (index / naxis1 % naxis2).astype(int)

    # Skipping distortion
    
    # Apply CRPIX values
    crpix = my_wcs['crpix']
    coord[0,:] = coord[0,:] - (crpix[0] -1)
    coord[1,:] = coord[1,:] - (crpix[1] -1)
        
    # Skipping distortion/associate/pixel-list (218-234)
    
    # Calcualte immedate (relative coordinates)
    # Assuming were doing pc
    coord = np.matmul(my_wcs['pc'], coord)
    # Skipping more distortion (264 - 284)
    coord[0,:] = coord[0,:]*my_wcs['cdelt'][0]
    coord[1,:] = coord[1,:]*my_wcs['cdelt'][1]
    
    # Skipping more distortion (288-301)
    
    # Assume we don't just want relative proj
    
    # Projection table - assume we don't need for now but check and bail if so
    ctypes = my_wcs['ctype']
    for item in ctypes:
        if '-TAB' in item:
            print ('Need to do wcs_proj_tab but havent ported yet')
            print (Quit)
    
    # Assume we dont hit any of the weird cases of proj and crval already dealt with (319-350)
    
    # Apply spherical proj
    proj = my_wcs['projection']
    if proj == 'TAN':
        coord = wcs_proj_tan(my_wcs, coord)
    elif proj == 'AZP':
        coord = wcs_proj_azp(my_wcs, coord)
        
    else:
        print('Other projections not yet ported')
        print(Quit)
       
    # Skipping projextion, pos_long, nowrap since not hit in simple version
    
    # Reformat
    coord = coord.reshape([2, naxis1, naxis2])
    return coord    
    
    
def wcs_get_pixel(my_wcs, coord,  doQuick=False, force_proj=False, noPC=False):
    # believe that these input coords should be heliocartesian (or like wcs)
    # and match the cunit in terms of arcsec...
    if isinstance(coord, list):
        coord = np.array(coord)
    # Check shape of input array
    singlePt = False
    if len(coord.shape) == 1:
        singlePt = True
        coord = coord.reshape([2,1])
    
    if my_wcs['projection'] == 'TAN':
        outpix = wcs_inv_proj_tan(my_wcs, coord, doQuick=doQuick)
    elif my_wcs['projection'] == 'AZP':
        outpix = wcs_inv_proj_azp(my_wcs, coord)
    else:
        print ('Other projections not ported')
        print (Quit)
        
    # Skipping subtract ref values for non-spherical and de-app tablular
    
    coord[0,:] = outpix[0,:] / my_wcs['cdelt'][0]
    coord[1,:] = outpix[1,:] / my_wcs['cdelt'][1]
    #temp = np.copy(coord)
    if not noPC:
        pc = my_wcs['pc']
        coord = np.matmul(np.transpose(pc), coord)
    
    # Add in reference pixel
    coord[0,:] = coord[0,:] + my_wcs['crpix'][0] -1
    coord[1,:] = coord[1,:] + my_wcs['crpix'][1] -1
    
    if singlePt:
        coord = coord.flatten()
    return coord
        
    
    
    
    
    