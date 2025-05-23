import numpy as np
import os
import sys
import shutil

CopyOrMove = 'move'

tempPath = '/Users/kaycd1/STEREO_Mass/MassFits/'
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'


corsetDict = {}
cdFile = np.genfromtxt('corsetDirectory.txt', dtype=None)
for line in cdFile:
    corsetDict[line[0]] = line[1]


allFiles = os.listdir(tempPath)
if '.DS_Store' in allFiles:
    allFiles.remove('.DS_Store')
    
txtFiles = []
pngFiles = []
ftsFiles = []
uniqDates = []
codes     = []
for aFile in allFiles:
    if aFile[-4:] in ['.fts', '.txt', '.png']:
        if aFile[-4:] == '.fts':
            ftsFiles.append(aFile)
        elif aFile[-4:] == '.png':
            pngFiles.append(aFile)
        elif aFile[-4:] == '.txt':
            txtFiles.append(aFile)
        
            
# Assuming that every png/fits has a text friend
code2Date = {}
code2Path = {}
nF = len(txtFiles)
counter = 0
for aFile in txtFiles:
    counter +=1
    print ('On text file ', counter, ' of ', nF)
    lessFile = aFile.replace('_CORSETmass.txt', '')
    splitidx = lessFile.find('_')
    corID = lessFile[:splitidx+3]
    sat   = lessFile[splitidx+2]
    myDate = lessFile[splitidx+4:]
    code2Date[corID] = myDate
    
    # Set up the path files
    yr = myDate[0:4]
    date2 = myDate[2:]
    mn = myDate[4:6]
    if sat == "A": satL = 'a'
    if sat == "B": satL = 'b'
    
    newPath = basePath + sat + '/' + yr + '/' + mn + '/' + date2 + '/' + corID + '/' 
    code2Path[corID] = newPath
    
    # Move the txt file
    if CopyOrMove in ['copy', 'Copy', 'COPY', 'cp']:
        shutil.copyfile(tempPath+aFile, newPath+aFile)
    elif CopyOrMove in ['move', 'Move', 'MOVE', 'mv']:
        shutil.move(tempPath+aFile, newPath+aFile)
    #print (tempPath+aFile,' moves to ', newPath+aFile)
    

# do the png files
nF = len(pngFiles)
counter = 0
for aFile in pngFiles:
    counter +=1
    print ('On png file ', counter, ' of ', nF)
    splitidx = aFile.find('_')
    corID = aFile[:splitidx+3]
    if corID in code2Path.keys():
        myPath = code2Path[corID]
    else:
        sys.exit('Missing .txt file for '+corID)
    print (myPath)

    # Check if the /png exists
    pngPath = myPath + 'png/'
    if not os.path.exists(pngPath):
        os.mkdir(pngPath)
    
    if CopyOrMove in ['copy', 'Copy', 'COPY', 'cp']:
        shutil.copyfile(tempPath+aFile, pngPath+aFile)
    elif CopyOrMove in ['move', 'Move', 'MOVE', 'mv']:
        shutil.move(tempPath+aFile, pngPath+aFile)
    print (tempPath+aFile, pngPath+aFile)
    #print (tempPath+aFile + ' goes to ' + pngPath+aFile)
    
    
# do the fts files
nF = len(ftsFiles)
counter = 0
for aFile in ftsFiles:
    counter +=1
    print ('On fts file ', counter, ' of ', nF)
    splitidx = aFile.find('_')
    corID = aFile[:splitidx+3]
    myPath = corsetDict[corID] + corID + '/'
    '''if corID in code2Path.keys():
        myPath = code2Path[corID]
    else:
        idx = aFile.find('_')
        sat  = aFile[idx+2]
        date = aFile[(idx+4):(idx+12)]
        yr = date[:4]
        mn = date[4:6]
        date2 = date[2:]
        dateB4 = str(int(date) - 1)[2:]
        myPath = basePath + sat + '/' + yr + '/' + mn + '/' + date2 + '/' + corID + '/' 
        if not os.path.exists(myPath):
            myPath = basePath + sat + '/' + yr + '/' + mn + '/' + dateB4 + '/' + corID + '/' '''
        #print (newPath)
        #sys.exit('Missing .txt file for '+corID)

    # Check if the /fts exists
    ftsPath = myPath + 'fts/'
    if os.path.exists(myPath):
        if not os.path.exists(ftsPath):
            os.mkdir(ftsPath)
    
        try:
            if CopyOrMove in ['copy', 'Copy', 'COPY', 'cp']:
                shutil.copyfile(tempPath+aFile, ftsPath+aFile)
            elif CopyOrMove in ['move', 'Move', 'MOVE', 'mv']:
                shutil.move(tempPath+aFile, ftsPath+aFile)
            print (tempPath+aFile, ftsPath+aFile)
        except:
            print('error moving ', aFile)

