from scipy.io import readsav
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'
date = '20120712'
corID = '2087.1938_0A'
yr = date[0:4]    
date2 = date[2:]
mn = date[4:6]
sat = 'A'
fullPath = basePath + sat + '/' + yr + '/' + mn + '/' + date2 + '/' + corID + '/' + corID
savFile = fullPath +'.sav'
sav_data = readsav(savFile)
for item in sav_data['cat'][0]:
    print (item)
    print ('')