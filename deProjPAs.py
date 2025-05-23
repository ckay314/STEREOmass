import numpy as np
import sunpy.map
from scipy.io import readsav


dtor = np.pi / 180.

fnameA1 = '/Users/kaycd1/pyGCS/fits/20120712_162400_d4c2A.fts'
fnameA2 = '/Users/kaycd1/pyGCS/fits/20120712_185400_d4c2A.fts' 

fnameB1 = '/Users/kaycd1/pyGCS/fits/20120712_162400_d4c2B.fts'
fnameB2 = '/Users/kaycd1/pyGCS/fits/20120712_185400_d4c2B.fts' 

corA = sunpy.map.Map(fnameA2)
corB = sunpy.map.Map(fnameB2)

corAlon = corA.observer_coordinate.lon.deg * dtor
corAlat = corA.observer_coordinate.lat.deg * dtor

corBlon = corB.observer_coordinate.lon.deg * dtor
corBlat = corB.observer_coordinate.lat.deg * dtor

savFileA = '2087.1938_0A.sav'
savFileB = '2087.2042_0B.sav'

sav_dataA = readsav(savFileA)
rA = sav_dataA['cat'][0][8][0][0][4]
cpaA = sav_dataA['cat'][0][4][0][1][4] * dtor


sav_dataB = readsav(savFileB)
rB = sav_dataB['cat'][0][8][0][0][4]
cpaB = sav_dataB['cat'][0][4][0][1][4] * dtor



#corBlon = - corAlon

'''denom = np.sin(corAlon) * np.cos(corBlon) + np.cos(corAlon) * np.sin(corBlon)
x = (rA * np.sin(cpaA) * np.cos(corBlon) - rB * np.sin(cpaB) * np.cos(corAlon)) / denom
y = - (rA * np.sin(cpaA) * np.sin(corBlon) + rB * np.sin(cpaB) * np.sin(corAlon)) / denom
print (rA, rB)
print (cpaA, cpaB)
print (corAlon, corBlon)
print ('')
print (denom, np.sin(cpaA) * np.cos(corBlon) -  np.sin(cpaB) * np.cos(corAlon) )
print (x, y)
print(rA * np.cos(cpaA), rB * np.cos(cpaB))'''

print (corAlon/dtor, corBlon/dtor)
x = (rA * np.sin(cpaA) * np.cos(corBlon) - rB * np.sin(cpaB) * np.cos(corAlon)) 
# 99% certain the - in Laura's eq 4 should be in numerator first term, not on whole thing
# the calc matches the Mierla version if we correct as such
y = -(rA * np.sin(cpaA) * np.sin(corBlon) - rB * np.sin(cpaB) * np.sin(corAlon)) 
lon = np.arctan(y/x)/dtor
z = 0.5*(rB * np.cos(cpaB) + rA * np.cos(cpaA))
lat = np.arctan(z/np.sqrt(x**2 + y**2))/dtor
r3d = np.sqrt(x**2 + y**2 + z**2)
print (r3d, lat, lon)

'''gamma = corAlon - corBlon
alpha = rB * np.sin(cpaB) / np.sin(gamma)
beta = - rA * np.sin(cpaA) / np.sin(gamma)
a = - rA * np.sin(cpaA)
b = rB * np.sin(cpaB)
z = 0.5*(rB * np.cos(cpaB) + rA * np.cos(cpaA))

print (y,(a-b)*np.sin(corBlon))
print (x,(a+b)*np.cos(corBlon))

R2d = np.sqrt(alpha**2 + beta**2 + 2 * alpha * beta * np.cos(gamma))
R3d = np.sqrt(R2d**2 + z**2)
lat = np.arctan(z/ R2d)
lon = np.arctan( np.tan(gamma/2) * (a-b) / (a+b))
print (R3d, R2d, lat/dtor, (0.5*(corAlon+corBlon) + lon)/dtor, lat2, lon2)'''

# looking for late -13.4, lon 1.5
