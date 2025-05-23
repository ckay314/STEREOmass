import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.coordinates import (CartesianRepresentation, Distance, SkyCoord,
                                 SphericalRepresentation, concatenate, Angle)
import astropy.units as u                                 
from numpy.linalg import norm
from scipy.spatial.transform import Rotation
from sunpy.coordinates import frames
import datetime

dtor = np.pi / 180

# Make a fake axis
thetas = np.linspace(-np.pi/2, np.pi/2, 100)
dcent  = 100
rCS   = 50
xs = dcent + rCS * np.cos(thetas)
ys = np.zeros(100)
zs = rCS * np.sin(thetas)


# Rotate into real world
lat = 15
lon = 80
tilt = 69

aTime = datetime.datetime(2012,2,22,2,22,2)
c1 = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, radius=1*u.AU, frame='heliographic_carrington', observer='earth', obstime=aTime)
print (c1)
c2 = c1.transform_to('heliographic_stonyhurst')
print (c2)

'''v = np.array([xs, ys, zs]).transpose()
#v = Rotation.from_euler('z', [-tilt*dtor, -lat*dtor + np.pi/2, lon*dtor]).apply(v)
v = Rotation.from_euler('xyz', [tilt*dtor-np.pi/2, -lat*dtor, lon*dtor]).apply(v)


spher_rep = SphericalRepresentation(Angle(lon, unit=u.deg), Angle(lat, unit=u.deg), Distance((dcent+rCS)*7e10, unit=u.m))
print (spher_rep)
print (spher_rep.to_cartesian())
center = SkyCoord(spher_rep.to_cartesian(),frame=frames.HeliographicStonyhurst, observer='earth')
print (center.observer)

aTime = datetime.datetime(2012,2,22,2,22,2)
stuff = SkyCoord(CartesianRepresentation(v[:,0], v[:,1], v[:,2]) ,frame=frames.HeliographicStonyhurst, observer='earth', obstime=aTime)
print (stuff[3])

observer1 = SkyCoord(lon=45*u.deg, lat=10*u.deg, radius=1*u.AU, frame='heliographic_stonyhurst')
stuff = SkyCoord(CartesianRepresentation(v[:,0], v[:,1], v[:,2]) ,frame=frames.HeliographicStonyhurst, observer=observer1, obstime=aTime)
print (stuff[3])

if False:
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(xs, ys, zs)
    ax.scatter(v[:,0], v[:,1], v[:,2])
    ax.set_aspect('equal')
    ax.view_init(elev=0, azim=-0, roll=0)
    plt.show()'''