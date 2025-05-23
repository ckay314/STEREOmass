# Following example on how to combine coordinate frameworks and an obs
# from https://docs.sunpy.org/en/latest/generated/gallery/showcase/los_simulation_box_intersection.html

import itertools

import matplotlib.pyplot as plt
import numpy as np

import astropy.time
import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
import sunpy.sun.constants
from sunpy.coordinates import Heliocentric, Helioprojective

# Helper function to make a fake map of data
def make_fake_map(ref_coord):
    # make a grid pf nans
    instrument_data = np.nan * np.ones((25, 25))
    # make a fake header using the ref_coords for the observer
    instrument_header = sunpy.map.make_fitswcs_header(instrument_data,ref_coord,
scale=u.Quantity([20,20])*u.arcsec/u.pix)
    return sunpy.map.Map(instrument_data, instrument_header)
    
# Helper function to determine if a coord is within the box (as defined by edges)
# Returns a boolean mask
def is_coord_in_box(box_edges, coord):
    box_edges = SkyCoord(box_edges)
    coord_hcc = coord.transform_to(box_edges.frame)
    in_x = np.logical_and(coord_hcc.x<box_edges.x.max(), coord_hcc.x>box_edges.x.min())
    in_y = np.logical_and(coord_hcc.y<box_edges.y.max(), coord_hcc.y>box_edges.y.min())
    in_z = np.logical_and(coord_hcc.z<box_edges.z.max(), coord_hcc.z>box_edges.z.min())
    return np.all([in_x, in_y, in_z], axis=0)
 


def doIt():
    # Set up a coord in HGS frame that we will use to define a box in a heliocentric cartesian frame
    hcc_orientation = SkyCoord(lon=0*u.deg, lat=20*u.deg,
                               radius=sunpy.sun.constants.radius,
                               frame='heliographic_stonyhurst')
    # Get that heliocentric frame at a specific time                           
    date = astropy.time.Time('2020-01-01')
    frame_hcc = Heliocentric(observer=hcc_orientation, obstime=date)
    
    # Make a fake observer. Could replace this with a real sat
    observer = SkyCoord(lon=45*u.deg, lat=10*u.deg, radius=1*u.AU, frame='heliographic_stonyhurst')
    
    # Get that observers projected frame
    frame_hpc = Helioprojective(observer=observer, obstime=date)
    
    # Make a set of fake data centered at the first point we defined
    instrument_map = make_fake_map(hcc_orientation.transform_to(frame_hpc))
    
    # use build in map functions to get the coords of the center of each pixel 
    map_indices = sunpy.map.all_pixel_indices_from_map(instrument_map).value.astype(int)
    map_indices = map_indices.reshape((2, map_indices.shape[1]*map_indices.shape[2]))
    
    # convert each pixel index to a line of sight
    lines_of_sight = []
    # take a chunk of distances to make it 3D
    distance = np.linspace(0.99, 1.01, 10000)*observer.radius
    for indices in map_indices.T:
        # Go from pixels to Helioprojective coords
        coord = instrument_map.wcs.pixel_to_world(*indices)
        # Get the line of sight for each pixel within the observer frame
        lines_of_sight.append(SkyCoord(Tx=coord.Tx, Ty=coord.Ty, distance=distance, frame=coord.frame))
        
        
    # Make a fake box to overlay    
    # Set up dimensions
    box_dimensions = u.Quantity([100,100,200])*u.Mm
    # Make a list of the corners (center currently at 0,0,0)
    corners = list(itertools.product(box_dimensions[0]/2*[-1,1], box_dimensions[1]/2*[-1,1], box_dimensions[2]/2*[-1,1]))
    # Take combos of corners and find pairs with only one dim varying to make edges
    # -> list of paired points
    edges = []
    # Make all the combos
    for possible_edges in itertools.combinations(corners,2):
        # Take the difference of each pair
        diff_edges = u.Quantity(possible_edges[0])-u.Quantity(possible_edges[1])
        # Keep if only one of the three pts is nonzero
        if np.count_nonzero(diff_edges) == 1:
            edges.append(possible_edges)
        
    # Place the box in space
    # Start by setting at the point we first defined
    box_origin = hcc_orientation.transform_to(frame_hcc)
    # Move it so back end of box at solar surface instead
    box_origin = SkyCoord(x=box_origin.x, y=box_origin.y, z=box_origin.z+box_dimensions[2]/2, frame=box_origin.frame)
    
    # Use the box origin to convert the edges   
    edge_coords = []
    for edge in edges:
        edge_coords.append(SkyCoord(x=box_origin.x+u.Quantity([edge[0][0],edge[1][0]]),
                                    y=box_origin.y+u.Quantity([edge[0][1],edge[1][1]]),
                                    z=box_origin.z+u.Quantity([edge[0][2],edge[1][2]]),
                                    frame=box_origin.frame))
        
    # Plot the field of view with the solar map and the lines of sight
    if False:
        fig = plt.figure()
        # Project from the fake observer
        ax = fig.add_subplot(projection=instrument_map)
        # Something to set up plot? Affects style but not content
        instrument_map.plot(axes=ax, title=False)
        # Draw the surface grid
        instrument_map.draw_grid(axes=ax, color='k')
        
        # Draw one point at each LoS
        #for los in lines_of_sight:
        #    ax.plot_coord(los[0], color='C0', marker='.', ls='', markersize=1,)   
         
        # Draw the edges of the box
        for edge in edge_coords:
            ax.plot_coord(edge, color='k', ls='-', marker='')
        ax.plot_coord(box_origin, color='r', marker='x')
        ax.plot_coord(hcc_orientation, color='b', marker='x')
        
        # Plot the fake image grid
        # Plot the pixel center positions
        for i,los in enumerate(lines_of_sight):
            ax.plot_coord(los[0], color='C0', marker='.', label='LOS' if i==0 else None, markersize=1)
        # Plot the edges of the pixel grid
        xpix_edges = np.array(range(int(instrument_map.dimensions.x.value)+1))-0.5
        ypix_edges = np.array(range(int(instrument_map.dimensions.y.value)+1))-0.5
        ax.vlines(x=xpix_edges,
                  ymin=ypix_edges[0],
                  ymax=ypix_edges[-1],
                  color='k', ls='--', lw=.5, label='pixel grid')
        ax.hlines(y=ypix_edges,
                  xmin=xpix_edges[0],
                  xmax=xpix_edges[-1],
                  color='k', ls='--', lw=.5,)
        
        plt.show()
        
    # Second part to determine whether a LoS intersects the simulation box or not
    # Set up another map from a different perspective for better vieing
    new_obs = SkyCoord(lon=25*u.deg, lat=0*u.deg, radius=1*u.AU, frame='heliographic_stonyhurst')
    earth_map = make_fake_map(hcc_orientation.transform_to(Helioprojective(observer=new_obs, obstime=date)))
    
    fig = plt.figure()
    ax = fig.add_subplot(projection=earth_map)
    earth_map.plot(axes=ax, title=False)
    earth_map.draw_grid(axes=ax, color='k')
    for los in lines_of_sight:
        ax.plot_coord(los, color='C0', marker='', ls='-', alpha=0.2)
    
    for edge in edge_coords:
        ax.plot_coord(edge, color='k', ls='-')
    
    ax.set_xlim(-10, 35)
    ax.set_ylim(-10, 35)

    plt.show()
    
    return
    
doIt()

# Notes for making the GCS work
# 1. Set up observer at stereo location
# 2. rotate wireframe from theory to Stony, make skycoord
# 3. just project it using observer HPC