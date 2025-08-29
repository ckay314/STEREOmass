import numpy as np
import sys
from sunpy.time import parse_time
import os
import glob
import spiceypy as spice
from spiceTest import setupPSPkernels
import math

global scDict, pspOrbit
scDict = {'sta':'-234', 'stereoa':'-234', 'stereoahead':'-234', 'stb':'-235', 'stereob':'-235', 'stereobehind':'-235', 'solo':'-144', 'solarorbiter':'-144', 'psp':'-96','parkersolarprobe':'-96'}
pspOrbit = ['/Users/kaycd1/ssw/psp/gen/data/spice/orbit/spp_nom_20180812_20250831_v040_RO7.bsp']
spiceDir = {'psp':'/Users/kaycd1/ssw/psp/gen/data/spice'}


def get_sunspyce_hpc_point(date, spacecraft, instrument=None, doDeg=False, doRad=False):
    # returns yaw (arcsec), pitch (arcsec), and roll angle (deg)
    # If doDeg then all three params returned in deg
    
    # Define the unilts
    roll_units = 180 / np.pi
    xy_units   = roll_units * 3600
    if doDeg: xy_units = roll_units
    if doRad:
        roll_units = 1
        xy_units   = 1
        
    # Determine which spacecraft was requested and make it spicy    
    if spacecraft.lower() in scDict:
        sc = scDict[spacecraft.lower()]
    else:
        sys.exit(spacecraft.lower() + ' not in scDict for sunspyce. Pick from sta, stb, solo, psp.')
    
    sc_stereo = False
    if sc in ['-234', '-235']: sc_stereo = True
    
    cmat = get_sunspyce_cmat(date, spacecraft, system='HPC', instrument=instrument)

    # Skipping error stuff
    # don't need to predefine pointing
    
    halfpi = np.pi / 2.
    twopi  = np.pi * 2.
    
    if sc == '-96':
        cmat0 = np.matmul([[0,0,1],[-1,0,0],[0,-1,0]], cmat)
        roll, pitch, yaw = spice.m2eul(cmat0, 1,3,2)
        yaw = halfpi - yaw
        roll = roll + halfpi
        if np.abs(roll) > np.pi:
            roll = roll - math.copysign(twopi, roll)

    # Ignoring stereo post conjunction
    
    # correct any cases where pitch is greater than 90
    if np.abs(pitch) > halfpi:
        pitch = math.copysign(np.pi, pitch) - pitch
        yaw = yaw - math.copysign(np.pi, yaw)
        roll = roll - math.copysign(np.pi, roll)
        
    # Apply the units
    pointing = np.zeros(3)
    pointing[0] = yaw * xy_units
    pointing[1] = pitch * xy_units
    pointing[2] = -roll * roll_units   
    
    return pointing


def get_sunspyce_cmat(date, spacecraft, system=None, instrument=None, tolerance=None, sixVec=False):
    # Determine which spacecraft was requested and make it spicy    
    if spacecraft.lower() in scDict:
        sc = scDict[spacecraft.lower()]
    else:
        sys.exit(spacecraft.lower() + ' not in scDict for sunspyce. Pick from sta, stb, solo, psp.')
    
    sc_base = ''
    if sc == '-96':
        sc_base = 'SPP_SPACECRAFT'
    elif sc == '-144':
        sc_base = 'SOLO_SRF'
        
    time = parse_time(date).utc
    
    # load sunspice
    #load_sunspyce(sc)    
    # load sunspice att
    #load_sunspyce_att(sc, time)
    
    setupPSPkernels()
    
    # Determine which coord system is specified
    # assuming single value for now 
    if system:
        system = system.upper()
        if system == 'HEQ': system = 'HEEQ'
        elif system in 'CARRINGTON': system = 'CARRINGTON'
    else:
        system == 'RTN'
        
    # Assume not passed frame bc don't give it the option yet...
    frame = None
    if system in ['HGRTN', 'RTN', 'HPC']:
        if sc == '-96': frame = 'PSPHGRTN'
        elif sc == '-144': frame = 'SOLOHGRTN'
        elif sc == '-234': frame = 'STAHGRTN'
        elif sc == '-235': frame = 'STBHGRTN'
        else:
            sys.exit('Cannot pull frame from sc in get_sunspyce_cmat')
    # ignoring the other systems
    
    # Determine the tolerance
    if tolerance:
        tol = tolerance
    else:
        tol = 1000
        
    # Determine if use ITRF93 kernels - skipping bc don't give keyword
    
    # Convert date/time to eph time and then to sc clock time
    et = spice.str2et(date)
    
    nVec = 3
    if sixVec:
        nVec = 6
    # again single time val for now
    cmat = np.zeros([nVec, nVec])
    
    sclkdp = spice.sce2c(int(sc), et)
     
    # Adding frcode that gets hit by roll GEI code 
    if not frame:
        if system == 'GEI':
            frame = 'J2000' 

    cmat, clkout = spice.ckgp(int(sc)*1000, sclkdp, tol, frame)

    # Modify the c-matrix based on the instrument keyword
    rotMat = spice.pxform(sc_base, instrument, et)
    if np.abs(np.linalg.det(rotMat) - 1) > 1e-5:
        sys.exit('Invalid rotation matrix for instrument')
        
    # Solar orbiter thing ignoring for now
    ccmat = np.matmul(rotMat, cmat)
    
    # Assume c matrix was found
    
    # Apply any additional processing
    if system == 'HPC':
        ccmat = np.matmul(ccmat, [[0, 0, 1.], [1., 0, 0], [0, 1., 0]])
    
    # ignoring weird storing stuff
    return ccmat    
        
def get_sunspyce_roll(date, spacecraft, system=None, instrument=None, doRad=False, tolerance=None):
    # Assuming passed correct things
    units = 180. / np.pi
    
    if doRad:
        units = 1.
        
    if spacecraft in scDict:
        sc = scDict[spacecraft]
    else:
        sys.exit('Spacecraft not in spice codes')    
    
    sc_stereo = sc in ['-234', '-235']
    
    if system:
        system = system.upper()
    else:
        system = 'RTN'
    
    cmat = get_sunspyce_cmat(date, spacecraft, system=system, instrument=instrument, tolerance=tolerance)
    
    roll, pitch, yaw = 0., 0., 0.
    twopi = np.pi * 2.
    halfpi = np.pi / 2.
    
    if sc == '-96':
        cmat0 = np.matmul([[0,0,1],[-1,0,0],[0,-1,0]], cmat)
        roll, pitch, yaw = spice.m2eul(cmat0, 1,2,3)
        roll = -roll
        pitch = -pitch
        if np.abs(roll) > np.pi:
            roll = roll - math.copysign(twopi, roll)
    else:
        sys.exit('Only PSP ported in get_sunspyce_roll')
        
    # Skipping post conjuction
    
    # Correct any cases where pitch > 90 deg
    if np.abs(pitch) > halfpi:
        pitch = math.copysign(np.pi, pitch) - pitch
        yaw   = yaw - math.copysign(np.pi, yaw) 
        roll  = roll - math.copysign(np.pi, roll) 
    
    # Apply the units
    roll  = units * roll
    pitch = units * pitch
    yaw   = units * yaw        
        
    return roll, pitch, yaw














def load_sunspyce(sc):
    # assume things are available for now
    
    # load sunspice gen seems to run and return quickly without doing anything here  
    # add sunspice seems make sure if can find the right function file which we dont care about
    
    if sc == '-234':
        sys.exit('Have not ported STA portion in load sunspyce')
    elif sc == '-235':
        sys.exit('Have not ported STB portion in load sunspyce')
    elif sc == '-184':
        sys.exit('Have not ported SolO portion in load sunspyce')
    elif sc == '-96':
        load_sunspyce_psp()
    else:
        sys.exit('Unrecognized sc key in load_sunspyce')
        
def load_sunspyce_psp():
    global psp_spice, psp_spice_gen, psp_spice_sclk, psp_spice_orbit
    orbit = pspOrbit
    n_kernels = len(orbit)
    psp_spice = spiceDir['psp']
    psp_spice_gen = psp_spice + '/gen'
    psp_spice_sclk = psp_spice + '/operations_sclk_kernel'
    
    aFile = psp_spice+'/frame_files.dat'
    if os.path.isfile(aFile):
        frames = np.genfromtxt(aFile, dtype=str)

    # Get the spacecraft clock file
    allFiles = os.listdir(psp_spice_sclk)
    keepers = []
    for aFile in allFiles:
        if 'spp_sclk' in aFile:
            keepers.append(aFile)
    slck = psp_spice_sclk+'/'+np.sort(keepers)[-1]
    spice.furnsh(slck) # well that is an equiv at least
    
    # throw in the leap second file - CK add
    # IDL probably does elsewhere?
    lsk = psp_spice + '/naif0012.tls'
    if os.path.isfile(lsk):
        spice.furnsh(lsk)
    
    
    # Determine the default predictive orbit file
    psp_spice_orbit = psp_spice + '/orbit'
    def_orbit = psp_spice_orbit + '/spp_nom_20180812_20250831_v034_RO1_TCM1.bsp'
    
    ephFile = psp_spice + '/ephemerides.dat'
    if os.path.isfile(ephFile):
        ephData = np.genfromtxt(ephFile, dtype=str)
        testFile = psp_spice_orbit + '/' + ephData
        if os.path.isfile(testFile):
            def_orbit = testFile
            
    # Determine the predictive orbit file to use -> assume we are not passed an 
    # orbit file for now
    orbit = def_orbit
    
    # Load the predictive orbit file
    if not os.path.isfile(orbit):
        sys.exit('Cannot find orbit file ' + orbit)
    else:
        spice.furnsh(orbit)
        
    # Load any short term predictive orbit files
    predict_list = psp_spice + '/ephemeris_predict.dat'
    if os.path.isfile(predict_list):
        predict_path = psp_spice + '/ephemeris_predict'
        predData = np.genfromtxt(predict_list, dtype=str)
        # Assuming predData is single line
        predFile = predict_path + '/' + predData
        if os.path.isfile(predFile):
            spice.furnsh, predFile
            ephem_predict = [predFile]
        else:
            print ('Predict file ' + predFile + ' not found')
            
    # Load any recon orbit files
    recon_ephem = []
    recon_list = psp_spice + '/reconstructed_ephemeris.dat'
    if os.path.isfile(recon_list):
        recon_path = psp_spice + '/reconstructed_ephemeris'
        reconData = np.genfromtxt(recon_list, dtype=str)
        for i in range(len(reconData)):
            reconFile = recon_path + '/' + reconData[i]
            if os.path.isfile(reconFile):
                spice.furnsh, reconFile
                recon_ephem.append(reconFile)
            else:
                print('Recon file ' + reconFile + ' not found')
                
    # Load any long term pred attitude history files
    att_pred = []
    att_pred_list = psp_spice + '/attitude_long_term_predict.dat'
    if os.path.isfile(att_pred_list):
        att_path = psp_spice + '/attitude_long_term_predict'
        attData = np.genfromtxt(att_pred_list, dtype=str)
        for i in range(len(attData)):
            attFile = att_path + '/' + attData[i]
            if os.path.isfile(attFile):
                spice.furnsh, attFile
                att_pred.append(attFile)
                
    # Also load any short term att files
    att_pred_list = psp_spice + '/attitude_short_term_predict.dat'
    if os.path.isfile(att_pred_list):
        att_path = psp_spice + '/attitude_short_term_predict'
        attData = np.genfromtxt(att_pred_list, dtype=str)
        for i in range(len(attData)):
            attFile = att_path + '/' + attData[i]
            if os.path.isfile(attFile):
                spice.furnsh, attFile
                att_pred.append(attFile)

def load_sunspyce_att(sc, date):
    if sc == '-234':
        sys.exit('Have not ported STA portion in load sunspyce')
    elif sc == '-235':
        sys.exit('Have not ported STB portion in load sunspyce')
    elif sc == '-184':
        sys.exit('Have not ported SolO portion in load sunspyce')
    elif sc == '-96':
        load_sunspyce_psp()
        load_sunspyce_att_psp(date)
    else:
        sys.exit('Unrecognized sc key in load_sunspyce')
    
        
def load_sunspyce_att_psp(date):
    dt = date.datetime
    doy = dt.timetuple().tm_yday
    sdate = str(dt.year)+'_'+str(doy).zfill(3)
    # Assuming single val not array for now
    # assuming no att date
    
    psp_spice = spiceDir['psp']
    attDir = psp_spice + '/attitude_history'
    
    # IDL isn't finding any files here so skipping the rest of this
    
    
        
    