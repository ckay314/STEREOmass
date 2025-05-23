import numpy as np

def cleanTheCSV(fileIn, holder):
    data = np.genfromtxt(fileIn, dtype=str, delimiter=',', skip_header=1)
    names = data[:,1]
    for i in range(len(names)):
        aName = names[i]
        if 'NULL' not in aName:
            holder.append(data[i,:11])

# all the CSV files
allCSV = ['corset_cme_list_2007.csv', 'corset_cme_list_2008.csv', 'corset_cme_list_2009.csv', 'corset_cme_list_2010.csv', 'corset_cme_list_2011.csv', 'corset_cme_list_2012.csv', 'corset_cme_list_2013.csv', 'corset_cme_list_2014.csv']

# File to save all the things
f1 = open('allGood.txt', 'w')

for aCSV in allCSV:
    holder = []
    
    # Get the year from the name
    myYr = aCSV.replace('corset_cme_list_', '').replace('.csv', '')
    f2 = open('allCORSET'+myYr+'.txt', 'w')
    
    # Get the ones tagged with a c, which seem to be the
    # ones that worked and have a mask
    cleanTheCSV(aCSV, holder)
    
    # Save these ones to a big file and a year file
    for line in holder:
        outline = ''
        for item in line:
            outline = outline + item + ' '
        outline = outline +'\n'
        f1.write(outline)
        f2.write(outline)
    f2.close()
f1.close()
