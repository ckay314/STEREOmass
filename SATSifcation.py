import numpy as np
import matplotlib.pyplot as plt


def coveragePlots():
    nPts = 3600
    sepLim = 50
    nSatPts = 1801

    allVals = np.zeros([nPts,nSatPts])
    satLons = np.linspace(0,90, nSatPts)

    for i in range(nSatPts):
        satA = satLons[i]
        satB = -satLons[i]

        posAl = (satA - 90) % 360
        posAr = (satA + 90) % 360 
        posBl = (satB - 90) % 360
        posBr = (satB + 90) % 360


        #allAngsA = np.linspace(0,359,nPts)
        #allAngsB = np.linspace(0,359,nPts)
        allAngsA = np.linspace(-180,180,nPts)
        allAngsB = np.linspace(-180,180,nPts)

        diffAl = np.abs(allAngsA - posAl) % 360
        diffAr =  np.abs(allAngsA - posAr) % 360
        diffBl = np.abs(allAngsB - posBl) % 360
        diffBr =  np.abs(allAngsB - posBr) % 360

        diffAl[np.where(diffAl > 180)] = 360 - diffAl[np.where(diffAl > 180)]
        diffAr[np.where(diffAr > 180)] = 360 - diffAr[np.where(diffAr > 180)]
        diffA = np.min([diffAl, diffAr], axis=0)

        diffBl[np.where(diffBl > 180)] = 360 - diffBl[np.where(diffBl > 180)]
        diffBr[np.where(diffBr > 180)] = 360 - diffBr[np.where(diffBr > 180)]
        diffB = np.min([diffBl, diffBr], axis=0)
    

        goodA = np.where(diffA <= sepLim)[0]
        goodB = np.where(diffB <= sepLim)[0]
        goodAB = np.where((diffA <= sepLim) & (diffB <= sepLim))[0]
        #print (goodAB)
        isGood = np.zeros(nPts)
        isGood[goodA] = 1
        isGood[goodB] = 2
        isGood[goodAB] = 3
    
        allVals[:,i] = isGood


    # recollect things to group color
    allOrd = allVals*0
    for i in range(nSatPts):
        myline = allVals[:,i]
        nBad = len(np.where(myline == 0)[0])
        nA   = len(np.where(myline == 1)[0])
        nB   = len(np.where(myline == 2)[0])
        nGood = len(np.where(myline == 3)[0])
        allOrd[:nBad,i] = 0
        allOrd[nBad:(nBad+nA),i] = 1
        allOrd[(nBad+nA):(nBad+nB+nA),i] = 2
        allOrd[(nBad+nB+nA):,i] = 3
        #print (i, nBad, nA, nB, nGood, nBad + nA + nB +nGood)
    allOrd = allOrd[::-1,:]

    plt.rcParams.update({'font.size': 14})
    ccols = ['w', '#FF9999', '#9999FF', '#995BC1']
    fig = plt.figure(constrained_layout=True, figsize=(10, 5))
    gs = fig.add_gridspec(1, 2)
    f1 = fig.add_subplot(gs[0,0])   
    f2 = fig.add_subplot(gs[0,1], sharex=f1)   
    f1.contourf(satLons, np.linspace(-179,180,nPts), allVals, levels=[-1, 0, 1,2,3], colors=ccols)
    f1.contour(satLons, np.linspace(-179,180,nPts), allVals, levels=[0, 1,2,3], colors='k')
    f1.plot([0,90],[0,90], '--', c='r', lw=1)
    f1.plot([0,90],[0,-90], '--', c='b', lw=1)
    f1.plot([0,90],[180,90], ':', c='b', lw=1)
    f1.plot([0,90],[-180,-90], ':', c='r', lw=1)
    f1.plot([0,90], [0,0], '--', c='k', lw=2)
    f1.plot([60,60], [-180,180], '--', c='#696969', lw=2)
    f2.plot([60,60], [0,100], '--', c='#696969', lw=2)
    f2.contourf(satLons, np.linspace(0,100,nPts), allOrd, levels=[-1, 0, 1,2,3], colors=ccols)
    f2.contour(satLons, np.linspace(0,100,nPts), allOrd, levels=[-1, 0, 1,2,3], colors='k')

    f1.set_ylim([-180,180])
    f1.set_xlabel('Sat Lons ($\\pm ^{\\circ}$)')
    f1.set_ylabel('CME Stony Lon ($^{\\circ}$)')
    f2.set_xlabel('Sat Lons ($\\pm ^{\\circ}$)')
    f2.set_ylabel('Percent Coverage')
    
    f2.text(0.5, 0.04, 'A+B', horizontalalignment='center', transform = f2.transAxes, fontsize=14, color='k')
    f2.text(0.5, 0.4, 'B', horizontalalignment='center', transform = f2.transAxes, fontsize=14, color='k')
    f2.text(0.5, 0.7, 'A', horizontalalignment='center', transform = f2.transAxes, fontsize=14, color='k')

    plt.savefig('paperFigs/coverage.png')
    #plt.show()
    
def sampleViews(satA, satB, sepLim, name):
    dtor = np.pi/180.
    satA = satA * dtor
    satB = satB * dtor
    sep  = sepLim * dtor
    fig = plt.figure(constrained_layout=True, figsize=(8, 8))
    gs = fig.add_gridspec(1, 1)
    ax = fig.add_subplot(gs[0,0], projection='polar')
    #ax = fig.add_subplot(111, projection='polar')
    npts = 1000
    thetas = np.linspace(0,2*np.pi, npts)
    ax.plot(thetas, np.ones(npts), 'k', lw=3)
    ax.plot([0,0], [0,1], 'k--', lw=2 )
    ax.scatter(0,0, s=200, c='y', ec='k', zorder=10)
    ax.scatter(satA, 1, s=200, c='r', ec='k', zorder=10)
    ax.scatter(satB, 1, s=200, c='b', ec='k', zorder=10)
    ax.scatter(0, 1, s=200, c='lightblue', ec='k', zorder=10)

    p2 = np.pi / 2
    posAl, posAr = satA - p2, satA + p2
    posBl, posBr = satB - p2, satB + p2
    
    ax.plot([posAl, posAr], [1,1], 'r--', lw=2, zorder=9)
    ax.plot([posBl, posBr], [1,1], 'b--', lw=2, zorder=9)
    
    angAr1, angAr2 = posAr + sep, posAr - sep,
    angAl1, angAl2 = posAl + sep, posAl - sep,
    if (angAr1 < angAr2):
        angAr1 += 2 * np.pi   
    angsAr = np.linspace(angAr2,angAr1, npts)    
    if (angAl1 < angAl2):
        angAl1 += 2 * np.pi   
    angsAl = np.linspace(angAl2,angAl1, npts)    
 
    angBr1, angBr2 = posBr + sep, posBr - sep,
    angBl1, angBl2 = posBl + sep, posBl - sep,
    if (angBr1 < angBr2):
        angBr1 += 2 * np.pi  
    angsBr = np.linspace(angBr2,angBr1, npts)            
    if (angBl1 < angBl2):
        angBl1 += 2 * np.pi  
    angsBl = np.linspace(angBl2,angBl1, npts)    
         
    
    ax.fill_between(angsAl, 0, 1, color='r', alpha=0.4)
    ax.fill_between(angsAr, 0, 1, color='r', alpha=0.4)
    ax.fill_between(angsBl, 0, 1, color='b', alpha=0.4)
    ax.fill_between(angsBr, 0, 1, color='b', alpha=0.4)

    

    '''x1AU = np.cos(thetas)
    y1AU = np.sin(thetas)
    
    f1.plot(x1AU, y1AU, 'k--')
    f1.scatter(0,0, s=200, c='y', ec='k', zorder=10)
    # plot positions as [y,-x] wrt to "normal" to get earth at bottom
    f1.scatter(np.sin(satA), -np.cos(satA),  s=200, c='r', ec='k', zorder=10)
    f1.scatter(np.sin(satB), -np.cos(satB),  s=200, c='b', ec='k', zorder=10)
    f1.scatter(0,-1,  s=200, ec='k', c='lightblue', zorder=10)
    
    # PoS lines
    p2 = np.pi / 2
    f1.plot([np.sin(satA-p2), np.sin(satA+p2)], [-np.cos(satA-p2), -np.cos(satA+p2)], c='r', zorder=9)
    f1.plot([np.sin(satA-p2-sep), np.sin(satA+p2-sep)], [-np.cos(satA-p2-sep), -np.cos(satA+p2-sep)], '--', c='r', zorder=9)
    f1.plot([np.sin(satA-p2+sep), np.sin(satA+p2+sep)], [-np.cos(satA-p2+sep), -np.cos(satA+p2+sep)], '--', c='r', zorder=9)
    f1.plot([np.sin(satB-p2), np.sin(satB+p2)], [-np.cos(satB-p2), -np.cos(satB+p2)], c='b', zorder=9)'''
    
    
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.savefig(name)
    
    
coveragePlots()    
satLoc = 80
sepLim = 50
sampleViews(satLoc, -satLoc, sepLim, 'paperFigs/sv80.png')