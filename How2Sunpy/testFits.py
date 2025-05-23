from astropy.io import fits
import sunpy.map
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, SqrtStretch
from astropy.coordinates import SkyCoord
import astropy.units as u



path1 = '/Users/kaycd1/IDLWorkspace/L1/A/20120712/20120712_162400_14c2A.fts'
path2 = '/Users/kaycd1/IDLWorkspace/L1/A/20120712/20120712_183900_14c2A.fts'

my_map1 = sunpy.map.Map(path1)
my_map2 = sunpy.map.Map(path2)

rd = sunpy.map.Map(my_map2 - my_map1.quantity)

hcc_orientation = SkyCoord(lon=2*u.deg, lat=-13*u.deg, radius=10*sunpy.sun.constants.radius, frame='heliographic_stonyhurst')

hcc_orientation2 = SkyCoord(lon=40*u.deg, lat=-13*u.deg, radius=10*sunpy.sun.constants.radius, frame='heliographic_stonyhurst')

hcc_orientation3 = SkyCoord(lon=-40*u.deg, lat=-13*u.deg, radius=10*sunpy.sun.constants.radius, frame='heliographic_stonyhurst')

box_origin = hcc_orientation.transform_to(my_map2.coordinate_frame)
box_origin2 = hcc_orientation2.transform_to(my_map2.coordinate_frame)
box_origin3 = hcc_orientation3.transform_to(my_map2.coordinate_frame)

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection=rd)
im = rd.plot(axes=ax, vmin=-5e-13, vmax=5e-13)
ax.plot_coord(box_origin, color='r', marker='x')
ax.plot_coord(box_origin2, color='r', marker='x')
ax.plot_coord(box_origin3, color='r', marker='x')

plt.show()