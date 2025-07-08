from astropy.io import fits
import numpy as np
import datetime
import platform
from astropy import wcs
import sys, os

# Hardcoded secchi backgrounds here  for now
global secchi_bkg
secchi_bkg = '/Volumes/SRP/vourla1/sswdb/stereo/secchi/backgrounds/'
global mjd_epoch
mjd_epoch = datetime.datetime(1858, 11, 17, 0, 0, 0)

def def_secchi_hdr():
    hdr = {}
    hdr['EXTEND'] = 'F'
    hdr['BITPIX'] = 0
    hdr['NAXIS'] = 0
    hdr['NAXIS1'] = 0
    hdr['NAXIS2'] = 0
    hdr['DATE_OBS'] = ''
    hdr['TIME_OBS'] = ''
    hdr['FILEORIG'] = ''
    hdr['SEB_PROG'] = ''
    hdr['SYNC'] = ''
    hdr['SPWX'] = 'F'
    hdr['EXPCMD'] = -1
    hdr['EXPTIME'] = -1
    hdr['DSTART1'] = 0
    hdr['DSTOP1'] = 0
    hdr['DSTART2'] = 0
    hdr['DSTOP2'] = 0
    hdr['P1COL'] = 0
    hdr['P2COL'] = 0
    hdr['P1ROW'] = 0
    hdr['P2ROW'] = 0
    hdr['R1COL'] = 0
    hdr['R2COL'] = 0
    hdr['R1ROW'] = 0
    hdr['R2ROW'] = 0
    hdr['RECTIFY'] = 'F'
    hdr['RECTROTA'] = 0
    hdr['RECTROTA'] = ''
    hdr['LEDPULSE'] = 0
    hdr['OFFSET'] = 9999
    hdr['BIASMEAN'] = 0.
    hdr['BIASSDEV'] = -1.
    hdr['GAINCMD'] = -1
    hdr['GAINMODE'] = '' 
    hdr['SUMMED'] = 0.
    hdr['SUMROW'] = 1
    hdr['SUMCOL'] = 1
    hdr['CEB_T'] = 999
    hdr['TEMP_CCD'] = 9999.
    hdr['POLAR'] = -1
    hdr['ENCODERP'] = -1
    hdr['WAVELNTH'] = 0
    hdr['FILTER'] = ''
    hdr['ENCODERQ'] = -1
    hdr['FPS_ON'] = ''
    hdr['OBS_PROG'] = 'schedule'
    hdr['DOORSTAT'] = -1
    hdr['SHUTTDIR'] = '' 
    hdr['READ_TBL'] = -1
    hdr['CLR_TBL'] = -1
    hdr['READFILE'] = ''
    hdr['DATE_CLR'] = ''
    hdr['DATE_RO'] = ''
    hdr['READTIME'] = -1.
    hdr['CLEARTIM'] = 0.
    hdr['IP_TIME'] = -1
    hdr['COMPRSSN'] = 0
    hdr['COMPFACT'] = 0
    hdr['NMISSING'] = -1.
    hdr['MISSLIST'] = ''
    hdr['SETUPTBL'] = ''
    hdr['EXPOSTBL'] = ''
    hdr['MASK_TBL'] = ''
    hdr['IP_TBL'] = ''
    hdr['COMMENT'] = ''
    hdr['HISTORY'] = ''
    hdr['DIV2CORR'] = 'F'
    hdr['DISTCORR'] = 'F'
    hdr['TEMPAFT1'] = 9999.
    hdr['TEMPAFT2'] = 9999.
    hdr['TEMPMID1'] = 9999.
    hdr['TEMPMID2'] = 9999.
    hdr['TEMPFWD1'] = 9999.
    hdr['TEMPFWD2'] = 9999.
    hdr['TEMPTHRM'] = 9999.
    hdr['TEMP_CEB'] = 9999.
    hdr['ORIGIN'] = ''
    hdr['DETECTOR'] = ''
    hdr['IMGCTR'] = 0
    hdr['TIMGCTR'] = 0
    hdr['OBJECT'] = ''
    hdr['FILENAME'] = ''
    hdr['DATE'] = ''
    hdr['INSTRUME'] = 'SECCHI'
    hdr['OBSRVTRY'] = ''
    hdr['TELESCOP'] = 'STEREO'
    hdr['WAVEFILE'] = ''
    hdr['CCDSUM'] = 0.
    hdr['IPSUM'] = 0.
    hdr['DATE_CMD'] = '' 
    hdr['DATE_AVG'] = ''
    hdr['DATE_END'] = ''
    hdr['OBT_TIME'] = 0.
    hdr['APID'] = 0
    hdr['OBS_ID'] = 0
    hdr['OBSSETID'] = 0
    hdr['IP_PROG0'] = 0
    hdr['IP_PROG1'] = 0
    hdr['IP_PROG2'] = 0
    hdr['IP_PROG3'] = 0
    hdr['IP_PROG4'] = 0
    hdr['IP_PROG5'] = 0
    hdr['IP_PROG6'] = 0
    hdr['IP_PROG7'] = 0
    hdr['IP_PROG8'] = 0
    hdr['IP_PROG8'] = 0
    hdr['IP_00_19'] = ''
    hdr['OBSERVER'] = ''
    hdr['BUNIT'] = ''
    hdr['BLANK'] = 0
    hdr['FPS_CMD'] = ''
    hdr['VERSION'] = ''
    hdr['CEB_STAT'] = -1
    hdr['CAM_STAT'] = -1
    hdr['READPORT'] = ''
    hdr['CMDOFFSE'] = 0.
    hdr['RO_DELAY'] = -1.
    hdr['LINE_CLR'] = -1.
    hdr['RAVG'] = -999.
    hdr['BSCALE'] = 1.0
    hdr['BZERO'] = 0.
    hdr['SCSTATUS'] = -1
    hdr['SCANT_ON'] = ''
    hdr['SCFP_ON'] = ''
    hdr['CADENCE'] = 0
    hdr['CRITEVT'] = ''
    hdr['EVENT'] = 'F'
    hdr['EVCOUNT'] = ''
    hdr['EVROW'] = 0
    hdr['EVCOL'] = 0
    hdr['COSMICS'] = 0
    hdr['N_IMAGES'] = 0
    hdr['VCHANNEL'] = 0
    hdr['OFFSETCR'] = 0.
    hdr['DOWNLINK'] = ''
    hdr['DATAMIN'] = -1.0
    hdr['DATAZER'] = -1.0
    hdr['DATASAT'] = -1.0
    hdr['DSATVAL'] = -1.0
    hdr['DATAAVG'] = -1.0
    hdr['DATASIG'] = -1.0
    hdr['DATAP01'] = -1.0
    hdr['DATAP10'] = -1.0
    hdr['DATAP25'] = -1.0
    hdr['DATAP50'] = -1.0
    hdr['DATAP75'] = -1.0
    hdr['DATAP90'] = -1.0
    hdr['DATAP95'] = -1.0
    hdr['DATAP98'] = -1.0
    hdr['DATAP99'] = -1.0
    hdr['CALFAC'] = 0.
    hdr['CRPIX1'] = 0
    hdr['CRPIX2'] = 0
    hdr['CRPIX1A'] = 0
    hdr['CRPIX2A'] = 0
    hdr['RSUN'] =  0.
    hdr['CTYPE1'] = 'HPLN-TAN'
    hdr['CTYPE2'] = 'HPLN-TAN'
    hdr['CRVAL1'] = 0.
    hdr['CRVAL2'] = 0.
    hdr['CROTA'] = 0.
    hdr['PC1_1'] = 1.
    hdr['PC1_2'] = 0.
    hdr['PC2_1'] = 0.
    hdr['PC2_2'] = 1.
    hdr['CUNIT1'] = ''
    hdr['CUNIT2'] = ''
    hdr['CDELT1'] = 0.
    hdr['CDELT2'] = 0.
    hdr['PV2_1'] = 0.
    hdr['PV2_1A'] = 0.
    hdr['SC_ROLL'] = 9999.
    hdr['SC_PITCH'] = 9999.
    hdr['SC_YAW'] = 9999.
    hdr['SC_PITA'] = 9999.
    hdr['SC_YAWA'] = 9999.
    hdr['INS_R0'] = 0.
    hdr['INS_Y0'] = 0.
    hdr['INS_X0'] = 0.
    hdr['CTYPE1A'] = 'RA---TAN'
    hdr['CTYPE2A'] = 'RA---TAN'
    hdr['CUNIT1A'] = 'deg'
    hdr['CUNIT1A'] = 'deg'
    hdr['CRVAL1A'] = 0.
    hdr['CRVAL2A'] = 0.
    hdr['PC1_1A'] = 1.
    hdr['PC1_2A'] = 0.
    hdr['PC2_1A'] = 0.
    hdr['PC2_2A'] = 1.
    hdr['CDELT1A'] = 0.
    hdr['CDELT2A'] = 0.
    hdr['CRLN_OBS'] = 0.
    hdr['CRLT_OBS'] = 0.
    hdr['XCEN'] = 9999.
    hdr['YCEN'] = 9999.
    hdr['EPHEMFIL'] = ''
    hdr['ATT_FILE'] = ''
    hdr['DSUN_OBS'] = 0.
    hdr['HCIX_OBS'] = 0.
    hdr['HCIY_OBS'] = 0.
    hdr['HCIZ_OBS'] = 0.
    hdr['HAEX_OBS'] = 0.
    hdr['HAEY_OBS'] = 0.
    hdr['HAEZ_OBS'] = 0.
    hdr['HEEX_OBS'] = 0.
    hdr['HEEY_OBS'] = 0.
    hdr['HEEZ_OBS'] = 0.
    hdr['HEQX_OBS'] = 0.
    hdr['HEQY_OBS'] = 0.
    hdr['HEQZ_OBS'] = 0.
    hdr['LONPOLE'] = 180
    hdr['HGLN_OBS'] = 0.
    hdr['HGLT_OBS'] = 0.
    hdr['EAR_TIME'] = 0.
    hdr['SUN_TIME'] = 0.
    #Skipping the EUV only keywords for now 
    
    return hdr

def fill_from_defhdr(hdr):
    mthdr   = def_secchi_hdr()
    allKeys = np.array(mthdr.keys())
    keys    = np.array(hdr.keys())
    for key in mthdr.keys():
        if key not in hdr.keys():
            hdr[key] = mthdr[key]
    return hdr

def sccrorigin(hdr):
    # Full port of sccrorigin
    # Could probably just make a dictionary
    p1col=51 # why are these in IDL?
    p1row=1 
    if hdr['rectify']:
        if hdr['OBSRVTRY'] == 'STEREO_A':
            det = hdr['detector']
            if det == 'EUVI':
                r1col, r1row = 129, 79
            elif det == 'COR1':
                r1col, r1row = 1, 79
            elif det == 'COR2':
                r1col, r1row = 129, 51
            elif det == 'HI1':
                r1col, r1row = 51, 1
            elif det == 'HI2':
                r1col, r1row = 51, 1
        elif hdr['OBSRVTRY'] == 'STEREO_B':
            det = hdr['detector']
            if det == 'EUVI':
                r1col, r1row = 1, 79
            elif det == 'COR1':
                r1col, r1row = 129, 51
            elif det == 'COR2':
                r1col, r1row = 1, 79
            elif det == 'HI1':
                r1col, r1row = 79, 129
            elif det == 'HI2':
                r1col, r1row = 79, 129
        else:
            # asuuming LASCO/EIT
            r1col, r1row = 20, 1
    else:
        r1col, r1row = 51, 1
    return [r1col, r1row]

def scc_make_array(filesIn, outSize=None, trim_off=False):
    # Port of IDL code
    # Starting at line 50
    
    # Cannot find where output_array common is defined but seems 
    # to always have the following values    
    out =  {'outsize':[2048.00, 2048.00], 'offset':[129, 2176, 51, 2098], 'readsize':[2048, 2048], 'binned':1.00000}
    
    num = len(filesIn)
    
    # IDL used readfits to just get the header and passes -1 to void (55)
    # void = sccreadfits(filenames,mhdr,/nodata)
    mhdr = []
    for i in range(num):
        with fits.open(filesIn[i]) as hdulist:
               mhdr.append(hdulist[0].header)
    
    if trim_off:
        # These seem to be the contents of the out common block...
        outsize  = [2176,2176]
        readsize = [2176,2176]
        offset   = [1,2176,1,2176]
        summed   = 1.
    
    else:
        # call sccrorigin
        # gets the 'rectified lower left (origin) value of full im area
        start = sccrorigin(mhdr[0])
        offset = [start[0],start[0]+2047,start[1],start[1]+2047]
            
        # Exclude over/underscan
        for i in range(num):
            r = mhdr[i]['R2COL']-mhdr[i]['R1COL']+mhdr[i]['R2ROW']-mhdr[i]['R1ROW']
            w = r < 2047*2
            
            # Getting max extent if there are sub_fields, doesn't seem to be triggered
            # so ignore for now but set flag for later (74-84)
            if w:
                print ('Have sub fields in ssc_make_array, need to add code')
                print (quit) # force a break
            
        readsize = np.array([offset[1]-offset[0]+1,offset[3]-offset[2]+1])
            
    if outSize:
        summed = np.max(readsize) / outSize[0]
        outSize = readsize / summed
    else:
        summeds = [mhdr[i]['SUMMED'] for i in range(num)]
        if len(np.unique(summeds)) != 1:
            print ('Have different summed values in scc_make_array, need to add code')
            print (quit)
        summed = 2 ** (summeds[0] -1)
        outSize = readsize / summed
            
            
    # Skipping polarized for now (102 - 108)
    
    out = {'outsize':outSize,'offset':offset,'readsize':readsize,'binned':summed}
    outout = np.max(outSize)
    
    # Create output arrays
    imgs = np.empty([int(outSize[0]), int(outSize[1]), num], dtype=float)
    # makes an empty header
    headers = []
    for i in range(num):
        aHdr = def_secchi_hdr()
        headers.append(aHdr)
        
    return imgs,  headers, int(outout), outSize.astype(int), out
        
def scc_zelensky_array(im, hdr, outsize, out):
    # Port of scc_putin_array
    
    # Assume we have been given proper header and out is defined from make_array
    
    # Set up some things
    ccdfac = 1.
    ccdsizeh = 2048
    if hdr['INSTRUME'] != 'SECCHI':
        print ('Cannot guarantee scc_zelelnsky_array works for not secchi cases')
        print (Quit)
    
    # Skipping to 141
    output_img = im
    
    if (hdr['naxis1'] != out['outsize'][0]) or (hdr['naxis2'] != out['outsize'][1]):
        print ('Havent addded resize code in scc_zelelnsky_array yet')
        print (Quit)
    
    return output_img
    
def secchi_rectify(a, scch, norotrate=False):
    # check not already rectified
    if scch['rectify']:
        print('We already done did the rectifying. Returning original img in secchi_rectify')
        return a, scch
    
    # Check if post conjunction
    post_conj = False 
    # Have to check if key actually exist bc stripped from some calibration files
    # Haven't tested a case that hits this
    if 'date_obs' in np.array(scch.keys):
        date_obs = scch['date_obs']
        try:
            doDT = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M" )
            cut1 = datetime.datetime(2015,7,1,0,0,0)
            cut2 = datetime.datetime(2023,8,12,0,0,0)
            if (doDT > cut1) & (doDT < cut2):
                post_conj = True
        except:
            print('No date in header to use in rectify. Assuming not post conjunction')
    
    stch = scch
    if ~norotrate:
        scch['rectify'] = True
        obs = scch['obsrvtry']
        if (obs == 'STEREO_A') & ~post_conj:
            det = scch['detector']  
            
            if det == 'COR2':
                # IDL -> b = rotate(a,1) same as np.rot90(k=3)
                b = np.rot90(a,k=3)
                stch['r1row'] = scch['p1col']
                stch['r1row']=scch['p1col']
                stch['r2row']=scch['p2col']
                stch['r1col']=2176-scch['p2row']+1
                stch['r2col']=2176-scch['p1row']+1
                stch['crpix1']=scch['naxis2']-scch['crpix2']+1
                stch['crpix2']=scch['crpix1']
                stch['naxis1']=scch['naxis2'] 
                stch['naxis2']=scch['naxis1']
                stch['sumrow']=scch['sumcol']
                stch['sumcol']=scch['sumrow']
                stch['rectrota']=1
                rotcmt='rotate 90 deg CCW'
                #--indicate imaging area - rotate 1
                #  naxis1
                stch['dstart1']	=(129-stch['r1col']+1)>1
                stch['dstop1']	=stch['dstart1']-1+((stch['r2col']-stch['r1col']+1)<2048)
                #  naxis2
                stch['dstart2']	=(51-stch['r1row']+1)>1
                stch['dstop2']	=stch['dstart2']-1+((stch['r2row']-stch['r1row']+1)<2048)
                
            else:
                print('Havent ported rectify things beyond COR2 so far')
                print(quit)
                
        if (obs == 'STEREO_B'):
            det = scch['detector']  
            
            if det == 'COR2':
                # IDL -> b = rotate(a,2) same as np.rot90(,k=1)
                b = np.rot90(a)
                stch['r1row']=2176-scch['p2col']+1
                stch['r2row']=2176-scch['p1col']+1
                stch['r1col']=scch['p1row']
                stch['r2col']=scch['p2row']
                stch['crpix1']=scch['crpix2']
                stch['crpix2']=scch['naxis1']-scch['crpix1']+1
                stch['naxis1']=scch['naxis2']
                stch['naxis2']=scch['naxis1']
                stch['sumrow']=scch['sumcol']
                stch['sumcol']=scch['sumrow']
                stch['rectrota']=3
                rotcmt='rotate 270 deg CCW'
                #--indicate imaging area - rotate 3
                #  naxis1
                stch['dstart1']	=1
                stch['dstop1']	=stch['dstart1']-1+((stch['r2col']-stch['r1col']+1)<2048)
                #  naxis2
                stch['dstart2']	=(79-stch['r1row']+1)>1
                stch['dstop2']	=stch['dstart2']-1+((stch['r2row']-stch['r1row']+1)<2048)
            
            else:
                print('Havent ported rectify things beyond COR2 so far')
                print(quit)
                
        if (obs == 'STEREO_A') & post_conj:
            det = scch['detector']
            if det == 'COR2':
                # IDL -> b = rotate(a,2) same as np.rot90(,k=1)
                b = np.rot90(a)
                stch['r1row']=2176-scch['p2col']+1
                stch['r2row']=2176-scch['p1col']+1
                stch['r1col']=scch['p1row']
                stch['r2col']=scch['p2row']
                stch['crpix1']=scch['crpix2']
                stch['crpix2']=scch['naxis1']-scch['crpix1']+1
                stch['naxis1']=scch['naxis2'] 
                stch['naxis2']=scch['naxis1']
                stch['sumrow']=scch['sumcol']
                stch['sumcol']=scch['sumrow']
                stch['rectrota']=3
                rotcmt='rotate 270 deg CCW'
                #--indicate imaging area - rotate 3
                #  naxis1
                stch['dstart1']	=1
                stch['dstop1']	=stch['dstart1']-1+((stch['r2col']-stch['r1col']+1)<2048)
                #  naxis2
                stch['dstart2']	=(79-stch['r1row']+1)>1
                stch['dstop2']	=stch['dstart2']-1+((stch['r2row']-stch['r1row']+1)<2048)
    else:
        stch.rectify = 'F'
        b = a 	    	# no change
        stch['r1row']=scch['p1row']
        stch['r2row']=scch['p2row']
        stch['r1col']=scch['p1col']
        stch['r2col']=scch['p2col']
        #stch.naxis1=scch.naxis1 & stch.naxis2=scch.naxis2
        stch['rectrota']=0
        rotcmt='no rotation'
        
    if stch['r1col'] < 1:
        stch['r2col'] = stch['r2col']+np.abs(stch['r1col'])+1
        stch['r1col'] = 1
    if stch['r1row'] < 1:
        stch['r2row'] = stch['r2row']+np.abs(stch['r1row'])+1
        stch['r1row'] = 1            
    
    # Don't have to explicitly add to header in python (skipping 485-509)
    
    # Reset the hdr to new values
    scch = stch
        
    return b, scch
    
def scc_sun_center(hdr):
    # assuming proper header, doing single not array of headers
    
    sunc = {'xcen':0, 'ycen':0.}
    my_wcs = wcs.WCS(hdr)
    
    # assuming scale is 1 bc not setting keywords right now
    scale = 1
    
    scen = wcs_get_pixel(my_wcs, [0,0])
    
    return suncen
    
def rebinIDL(arr, new_shape):
    factors = arr.shape // new_shape
    outarr = arr.reshape(new_shape[0], factors[0], new_shape[0], factors[0]).mean(3).mean(1)
    return outarr
    
def scc_getbkgimg(hdr, doRot=False):
    # Assume no interp for now
    
    # port of get_delim, untested on windows
    myos = platform.system()
    if myos == 'Windows':
        delim = '\\'
    else:
        delim = '/'
        
    # Assume proper header
    
    tel = hdr['DETECTOR'].upper().strip()
    if tel == 'EUVI':
        cam = 'eu'
    else:
        cam = tel[0].lower()+tel[-1]
    isHI = (cam == 'h1') or (cam == 'h2')
    
    doshift = False # is this still used? commented out below
    
    # Get sc name
    tags = np.array([key for key in hdr])
    filename = hdr['FILENAME']
    if'OBSRVTRY' not in tags:
        #MVI header
        sc = filename[-4].lower()
        print ('Not certain on pulling from MVI check below is expected for sc')
        print (sc)
    else:
        sc = hdr['OBSRVTRY'][7].lower()

    # Check cs name
    if sc not in ['a', 'b']:
        sys.exit('Error in scc_getbkgim. Invalid spacecraft name.')
        
    # Check processing level
    isL1 = False
    levchar = filename[16]
    if (levchar ==1) or ((levchar == 0) & isHI):
        isL1 = True
        match = False

    # Not setting up keywords raw_calroll, calroll, daily so skipping (401-413)
    ndays = 30
    rootdir = secchi_bkg + sc + delim + 'monthly_min' + delim
    fchar = 'm'

    # Form the polar search string
    polstr = '_pTBr_'
    # not using the totalb, double_totalB keys for now so skipping (423-427)
    
    # Get the date
    if 'DATE-AVG' in tags:
        dtin = hdr['DATE-AVG']
        daTag = 'DATE-AVG'
    elif 'DATE_AVG' in tags:
        dtin = hdr['DATE_AVG']
        daTag = 'DATE_AVG'
    cal = dtin[:4]+dtin[5:7]+dtin[8:10] # should pull out YYYYMMDD for reasonable strings
    sdir = cal[0:6]
    sfil = cal[2:8]    
    # assume we can do what we need with DT obj
    avgDT = datetime.datetime.strptime(dtin, "%Y-%m-%dT%H:%M:%S.%f" )
    mjdin = (avgDT - mjd_epoch).total_seconds()/(24*3600)
    
    # Check against major sc repointings
    if sc == 'a':
        repoint = ['2006-12-21T13:15', '2007-02-03T13:15','2015-05-19T00:00', '2023-08-12T00:00']
        if tel == 'COR1':
            repoint.extend(['2010-01-27T16:49', '2010-11-19T16:00', '2011-01-12T12:23', '2011-02-11T04:23', '2011-03-08T17:00', '2011-12-05T12:03', '2012-02-19T02:33', '2012-03-16T00:00', '2014-07-08T00:00', '2014-08-23T17:00', '2015-11-16T00:00', '2016-03-15T00:00', '2016-04-13T01:08', '2016-05-25T05:40', '2016-06-01T12:00', '2016-06-08T00:00', '2016-09-02T16:19', '2018-09-27T16:44', '2019-02-05T21:07', '2019-02-17T05:57', '2020-01-14T23:52', '2021-05-07T23:53', '2023-05-15T23:00'])
        else:
            repoint.extend(['2014-08-19T00:00'])
        nomrollmjd = 54160
    if sc == 'b':
        repoint = ['2007-02-03T18:20', '2007-02-21T20:00']
        if tel == 'COR1':
            repoint.extend(['2009-01-30T16:20', '2010-03-24T01:17', '2011-03-11T00:00'])
        if tel == 'COR2':
            repoint.extend(['2010-02-23T08:12', '2011-01-27T03:47','2011-04-25T18:30'])
        nomrollmjd = 54313
    repoint = np.array(sorted(repoint))        
    # add in check for no exist
    i1 = np.where(repoint < dtin)[0]
    i2 = np.where(repoint > dtin)[0]
    if len(i1) > 0:
        minDT = datetime.datetime.strptime(repoint[np.max(i1)], "%Y-%m-%dT%H:%M")
        mjdmin = (minDT - mjd_epoch).total_seconds()/(24*3600)
    else:
        mjdmin = 0
    if len(i2) > 0:        
        maxDT = datetime.datetime.strptime(repoint[np.min(i2)], "%Y-%m-%dT%H:%M")
        mjdmax = (maxDT - mjd_epoch).total_seconds()/(24*3600)
    else:
        mjdmax = 99999

    pntroll = hdr['sc_roll']
    if (mjdin > 57160) & (mjdin < 60168):
        pntroll = hdr['crota']
    if pntroll < 0:
        titleangle = 360 + pntroll
    else:
        titleangle = pntroll

    if (np.abs(pntroll) < 10) or (tel == 'COR1') or (mjdin < nomrollmjd):
            postd = ''
    else:
        print ("Reached uncoded part, need to test with example")
        print (QUit)
    
    if isHI:
        roll = 0
    
    # Look for the correct or closest file. IDL makes this a headache
    # This is equiv to 531-662 but ignoring interp and other stuff
    exactFile = rootdir + sdir + delim + fchar + cam + sc.upper() + polstr + sfil + postd + '.fts'
    if os.path.exists(exactFile):
        bkgFile = exactFile
    else:
        # loop through n days before and after to find closest friend
        dday = 0
        while dday <= ndays:
            dday += 1
            plusDT = avgDT + datetime.timedelta(days=dday)
            pstr = plusDT.strftime('%Y%m%d')  
            psdir, psfil =  pstr[:-2], pstr[2:]
            pPath = rootdir + psdir + delim + fchar + cam + sc.upper() + polstr + psfil + postd + '.fts'
            hasP =  os.path.exists(pPath)
            minusDT = avgDT + datetime.timedelta(days=-dday)
            mstr = minusDT.strftime('%Y%m%d')  
            msdir, msfil =  mstr[:-2], mstr[2:]
            mPath = rootdir + msdir + delim + fchar + cam + sc.upper() + polstr + msfil + postd + '.fts'
            hasM = os.path.exists(mPath)
            if hasP and hasM:
                dday = 9999
                delP = (plusDT - avgDT).total_seconds()
                delM = (avgDT - minusDT).total_seconds
                if delM <= delP:
                    bkgFile = mPath
                else:
                    bkgFile = pPath
            elif hasP:
                dday = 9999
                bkgFile = pPath
            elif hasM:
                dday = 9999
                bkgFile = mPath
                
            else:
                if dday > ndays:
                    sys.exit('Cannot fine appropriate background file')
                    
    # Skipping more for interp, postd, doubles 
    
    # Read in file, skipping 666-711
    if isL1 and not isHI:
        print ('Need to add secchi prep the background, not coded')
        print (Quit) 
    else:
        with fits.open(bkgFile) as hdulist:
            bim  = hdulist[0].data
            bhdr = hdulist[0].header
    # Skipping median keyword
    btags = np.array([key for key in bhdr])
    date_avg = ''
    if 'DATE-AVG' in btags:
        date_avg = bhdr['DATE-AVG']
    if 'date-avg' in btags:
        date_avg = bhdr['date-avg']
    elif 'DATE_AVG' in btags:
        date_avg = bhdr['DATE_AVG']
    elif 'date_avg' in btags:
        date_avg = bhdr['date_avg']
    if date_avg == '':
        if 'DATE_OBS' in btags:
            date_avg = bhdr['DATE_OBS']
        elif 'date_obs' in btags:
            date_avg = bhdr['date_obs']
        elif 'DATE-OBS' in btags:
            date_avg = bhdr['DATE-OBS']
        elif 'date-obs' in btags:
            date_avg = bhdr['date-obs']
        else:
            sys.exit('Issue getting background image date')
            
    # Match image size
    i_reduce = bhdr['naxis1'] / hdr['naxis1']
    j_reduce = bhdr['naxis2'] / hdr['naxis2']
        
    if (i_reduce != 1) or (j_reduce != 1):
        bim = rebinIDL(bIm, [hdr['naxis1'], hdr['naxis2']])
  
    # Copy header info for some reason (prob rebin?)
    if not isL1:
        bhdr['dstop1'] = hdr['dstop1']
        bhdr['dstop2'] = hdr['dstop2']
        bhdr['naxis1'] = hdr['naxis1']
        bhdr['naxis2'] = hdr['naxis2']
        bhdr['summed'] = hdr['summed']
        bhdr['crpix1'] = hdr['crpix1']
        bhdr['crpix2'] = hdr['crpix2']
        bhdr['crpix1a'] = hdr['crpix1a']
        bhdr['crpix2a'] = hdr['crpix2a']
        bhdr['CDELT1'] = hdr['CDELT1']
        bhdr['CDELT2'] = hdr['CDELT2']
        bhdr['CDELT1A'] = hdr['CDELT1A']
        bhdr['CDELT2A'] = hdr['CDELT2A']
        
    # Shift things for HI, ignore for now 760-786
    maxshift = 0
    if doshift:
        print ('havent coded HI shifting in scc_getbkgimg')
        print (Quit)
        
    # No interp skipping 791 - 854
    
    # skipping another shift thing 856-859
    
    # skipping match check 864
    
    # skipping header correction
    
    # Check for rotate correction
    rolldif = hdr['crota'] - bhdr['crota']
    if (np.abs(rolldif) > 1) and doRoll:
        print ('Havent coded background roll')
        print (Quit)
        
    # skipping double exposure correction
    
    # skipping match corrction
    
    # skipping HI L1 correction
    
    return bim, bhdr
                        
                
            
            
            
            
    
    
    
    
    
    