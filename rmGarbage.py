import numpy as np
import os
import sys
import shutil

mns = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
#mns = [ '05', '06', '07', '08', '09', '10', '11', '12']import numpy as np
import os
import sys
import shutil

mns = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
#mns = [ '05', '06', '07', '08', '09', '10', '11', '12']
#mns = ['03', '04', '05', '06','07', '08', '09', '10', '11', '12']

def rmGarbage(yr):
    basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'
    morePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/' + sat + '/' + yr + '/'
    for mn in mns:
        mnPath = morePath + mn + '/'
        if os.path.isdir(mnPath):
            allFiles = os.listdir(mnPath)
            allFolds = []
            for aFile in allFiles:
                if '.' not in aFile: 
                    allFolds.append(aFile)
            # Alphabetize for my own happiness
            allFolds = np.sort(np.array([aFold for aFold in allFolds]))     
            for aFold in allFolds:
                print (aFold)
                allFiles = os.listdir(mnPath+aFold+'/')
                subDirs = []
                for aFile in allFiles:
                    if ('.DS' not in aFile) & ('._' not in aFile): 
                        subDirs.append(aFile)
                for subDir in subDirs:
                    print ('   ', subDir)
                    allFiles = os.listdir(mnPath+aFold+'/'+subDir+'/')
                    massFile = subDir+'_20'+aFold+'_mass.txt'
                    if os.path.isfile(mnPath+aFold+'/'+subDir+'/'+massFile):
                        os.remove(mnPath+aFold+'/'+subDir+'/'+massFile)
                    if 'fts' in allFiles:
                        myFold = mnPath+aFold+'/'+subDir+'/'+'fts'
                        theseFiles = os.listdir(myFold)
                        for thisFile in theseFiles:
                            if os.path.isfile(myFold+'/'+thisFile):
                                try:
                                    print (thisFile)
                                    os.remove(myFold+'/'+thisFile)
                                except:
                                    print ('error rming ' + thisFile)
                                    pass
                    if 'png' in allFiles:
                        myFold = mnPath+aFold+'/'+subDir+'/'+'png'
                        theseFiles = os.listdir(myFold)
                        for thisFile in theseFiles:
                            if os.path.isfile(myFold+'/'+thisFile):
                                try:
                                    if thisFile[0] != '.':
                                        print(thisFile)        
                                        os.remove(myFold+'/'+thisFile)
                                except:
                                    print ('error rming ' + thisFile)
                                    pass
            
        
'''sat = 'B'    
rmGarbage('2013')
rmGarbage('2014')'''
sat = 'A'    
#rmGarbage('2013')
rmGarbage('2014')

