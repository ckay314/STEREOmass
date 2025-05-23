import numpy as np
import matplotlib.pyplot as plt
import datetime 
#from sunpy.coordinates import HeliocentricEarthEcliptic, get_body_heliographic_stonyhurst, get_horizons_coord
#from sunpy.time import parse_time
import pickle
import matplotlib.gridspec as gridspec
from makeMegaDataStructure import AllRes
from scipy import stats
import matplotlib.colors as mcolors
import matplotlib.dates as mdates


plotAll = False
myFold = 'paperFigs/'

cCORA = 'r'
cCORB = 'b'
cGCSA = '#882255'
cGCSB = '#88CCEE'


# |---------------------------------------|
# |-------- Repeated plot types ----------|
# |---------------------------------------|

def time2histo4(p1, p2, p3, p4, axLabs, doLog=True, picName='temp.png', pLabs=[None], colors=['r', 'b', 'maroon', 'lightblue'], timeYL=[None], bins=[[None]], plotFits=True):
    # assume we are given four variables (ps) as [time, val]
    # and if labs as [label1, label2]
    # the variables are not in log but will plot as such 
    # unless the log flag set to false
    fig = plt.figure(constrained_layout=True, figsize=(7,10))
    gs = fig.add_gridspec(4, 2)
    f1 = fig.add_subplot(gs[0,:])   
    f3a = fig.add_subplot(gs[2,0])
    f3b = fig.add_subplot(gs[2,1], sharex=f3a, sharey=f3a)
    if axLabs[0] == axLabs[1]: 
        f2 = fig.add_subplot(gs[1,:], sharex=f1, sharey=f1)    
        f4a = fig.add_subplot(gs[3,0], sharex=f3a, sharey=f3a)
        f4b = fig.add_subplot(gs[3,1], sharex=f3a, sharey=f3a)
    else:
        f2 = fig.add_subplot(gs[1,:], sharex=f1)    
        f4a = fig.add_subplot(gs[3,0])
        f4b = fig.add_subplot(gs[3,1], sharey=f4a)
    
    f1.scatter(p1[0], p1[1], c=colors[0], s=20)
    f1.scatter(p2[0], p2[1], c=colors[1], s=20)
    f2.scatter(p3[0], p3[1], c=colors[2], s=20)
    f2.scatter(p4[0], p4[1], c=colors[3], s=20)
    
    if pLabs[0]:
        f1.text(0.01, 0.9, pLabs[0], horizontalalignment='left', transform = f1.transAxes, fontsize=12, color=colors[0])
        f1.text(0.01, 0.82, pLabs[1], horizontalalignment='left', transform = f1.transAxes, fontsize=12, color=colors[1])
        f2.text(0.01, 0.9, pLabs[2], horizontalalignment='left', transform = f2.transAxes, fontsize=12, color=colors[2])
        f2.text(0.01, 0.82, pLabs[3], horizontalalignment='left', transform = f2.transAxes, fontsize=12, color=colors[3])
    
    if timeYL[0]:
        f1.set_ylim(timeYL[0])
        f2.set_ylim(timeYL[1])
    # Add horiz lines at nice values
    if doLog:
        f1.set_yscale('log')
        yrs0 = f1.get_ylim()
        yrs = [int(np.log10(yrs0[0])), int(np.log10(yrs0[1]))]
        diff = yrs[1] - yrs[0]
        toLine = [np.power(10, yrs[0]+1+x) for x in range(diff)]
        xrs = f1.get_xlim()
        for x in toLine:
            f1.plot(xrs, [x, x], 'k--',  alpha=0.5, zorder=0)
        f1.set_xlim(xrs)
        f1.set_ylim(yrs0)
        
        if axLabs[0] != axLabs[1]:
            f2.set_yscale('log')
            yrs0 = f2.get_ylim()
            yrs = [int(np.log10(yrs0[0])), int(np.log10(yrs0[1]))]
            diff = yrs[1] - yrs[0]
            toLine = [np.power(10, yrs[0]+1+x) for x in range(diff)]
            xrs = f2.get_xlim()
            
        for x in toLine:
            f2.plot(xrs, [x, x], 'k--',  alpha=0.5, zorder=0) 
        f2.set_xlim(xrs)
        f2.set_ylim(yrs0)
                       
        
    if (bins[0][0]) or (bins[0][0] ==0):
        bins1 = bins[0]
        bins2 = bins[1]
        f3a.set_xlim(bins1[0], bins1[-1])
        if axLabs[0] == axLabs[1]:
            f4a.set_xlim(bins2[0], bins2[-1])
        
    else:
        ll = int(np.log10(np.percentile(p1[1], 10))*2)/2
        ul = int(np.log10(np.percentile(p1[1], 99))*2 +1)/2
        bins1 = np.arange(ll,ul,0.1, dtype=float)
        if axLabs[0] == axLabs[1]:
            bins2 = bins1
            f3a.set_xlim([ll, ul])
        else:
            ll2 = int(np.log10(np.percentile(p3[1], 10))*2)/2
            ul2 = int(np.log10(np.percentile(p3[1], 90))*2 +1)/2
            bins2 = np.arange(ll2,ul2,0.1, dtype=float)
            f4a.set_xlim([ll2, ul2])
        
    if doLog:
        h1, h2, h3, h4 = np.log10(p1[1]), np.log10(p2[1]), np.log10(p3[1]), np.log10(p4[1])
    else:
        h1, h2, h3, h4 = p1[1], p2[1], p3[1], p4[1]
        
    f3a.hist(h1, bins=bins1, density=True, color=colors[0], ec='k')
    f3b.hist(h2, bins=bins1, density=True, color=colors[1], ec='k')
    f4a.hist(h3, bins=bins2, density=True, color=colors[2], ec='k')
    f4b.hist(h4, bins=bins2, density=True, color=colors[3], ec='k')
    
    axs = [f3a, f3b, f4a, f4b]
    hs  = [h1, h2, h3, h4]
    lls = [bins1[0], bins1[0],bins2[0], bins2[0]]
    uls = [bins1[-1], bins1[-1],bins2[-1], bins2[-1]]
    friendDict = {0:1, 1:0, 2:3, 3:2}
    
    if plotFits:
        for i in range(4):
            myh, myax, myF = hs[i], axs[i], axs[friendDict[i]]
            # add the skew normal distrubutions
            x = np.linspace(lls[i], uls[i], 1000)
            ae3a, loce3a, scalee3a = stats.skewnorm.fit(myh)
            p3a = stats.skewnorm.pdf(x, ae3a, loce3a, scalee3a)
            # mode
            Za = x[np.where(p3a == np.max(p3a))[0]][0]
            # pdf mean
            m0a  = np.sum(x*p3a) / np.sum(p3a)
       
        
            myax.plot(x, p3a, 'k--', linewidth=3, zorder=4)
            myF.plot(x, p3a, ':', c=colors[i], linewidth=3, zorder=3)
            myax.text(0.05, 0.9, ' N: '+str(len(myh)), horizontalalignment='left', transform = myax.transAxes, fontsize=12)
            myax.text(0.05, 0.82, '$\\mu_F$: '+ '{:.2f}'.format(loce3a), horizontalalignment='left', transform = myax.transAxes, fontsize=12)
            myax.text(0.05, 0.74, '$\\alpha_F$: '+ '{:.2f}'.format(ae3a), horizontalalignment='left', transform = myax.transAxes, fontsize=12)
            myax.text(0.05, 0.66, '$\\sigma_F$: '+ '{:.2f}'.format(scalee3a), horizontalalignment='left', transform = myax.transAxes, fontsize=12)
            myax.text(0.99, 0.82, '$\\mu$: '+ '{:.2f}'.format(m0a), horizontalalignment='right', transform = myax.transAxes, fontsize=12)
            myax.text(0.99, 0.74, '$Z$: '+ '{:.2f}'.format(Za), horizontalalignment='right', transform = myax.transAxes, fontsize=12)

    
    f1.set_ylabel(axLabs[0])
    f2.set_ylabel(axLabs[1])
    f1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    
    f3a.set_ylabel('Prob. Density')
    f4a.set_ylabel('Prob. Density')
    pref=''
    if doLog:
        pref = 'log '
    f3a.set_xlabel(pref+axLabs[0])
    f3b.set_xlabel(pref+axLabs[0])
    f4a.set_xlabel(pref+axLabs[1])
    f4b.set_xlabel(pref+axLabs[1])
    
    #f1.set_xlim([datetime.datetime(2007,1,1),datetime.datetime(2015,1,1)])
    
    #plt.show()
    plt.savefig(picName)
    #plt.close()



def scatter2histo2(p1, p2, p3, p4, axLabs, doLog=True, picName='temp.png', pLabs=[None], colors=['r', 'b', 'maroon', 'lightblue'], axLim=[None], bins=[[None]], scatCols=[[None]]):
    fig = plt.figure(constrained_layout=True, figsize=(12, 8.5))
    gs = fig.add_gridspec(3, 2)
    f1 = fig.add_subplot(gs[0:2,0]) 
    f2 = fig.add_subplot(gs[0:2,1], sharex=f1, sharey=f1) 
    f3a = fig.add_subplot(gs[2,0]) 
    f3b = fig.add_subplot(gs[2,1], sharex=f3a, sharey=f3a) 
    
    pn1, pn2, pn3, pn4 = p1, p2, p3, p4 # plot now
    if doLog:
        lp1, lp2, lp3, lp4 = np.log10(p1), np.log10(p2), np.log10(p3), np.log10(p4) 
        pn1, pn2, pn3, pn4 = lp1, lp2, lp3, lp4 

    if scatCols[0][0]:
        normCmap = mcolors.Normalize(vmin=scatCols[2][0], vmax=scatCols[2][1])
        if scatCols[0][1] == None:
            im = f1.scatter(pn1, pn2, c=scatCols[0][0], s=20, cmap='plasma', norm=normCmap)
        else:
            f1.scatter(pn1, pn2, c='lightgray', s=20)
            im = f1.scatter(pn1[scatCols[0][1]], pn2[scatCols[0][1]], c=scatCols[0][0], s=20, cmap='plasma', norm=normCmap)
        if scatCols[1][1] == None:
            im = f2.scatter(pn3, pn4, c=scatCols[1][0], s=20, cmap='plasma', norm=normCmap)
        else:
            f2.scatter(pn3, pn4, c='lightgray', s=20)
            im = f2.scatter(pn3[scatCols[1][1]], pn4[scatCols[1][1]], c=scatCols[1][0], s=20, cmap='plasma', norm=normCmap)
                    
        cb = plt.colorbar(im, shrink=0.75)
        cb.ax.set_title(scatCols[3],fontsize=12)
        cb.ax.tick_params(labelsize=12)
    else:
        f1.scatter(pn1, pn2, c=colors[0], s=20)
        f2.scatter(pn1, pn2, c=colors[1], s=20)
              
    # add lines at +50%
    fakex = np.linspace(axLim[0][0], axLim[0][1],30)
    f1.plot(fakex, fakex,  'k--',  alpha=0.75, zorder=0)
    f2.plot(fakex, fakex,  'k--',  alpha=0.75, zorder=0)
    if doLog:
        fakex15 = np.log10(1.5*np.power(10,fakex))
        fakex2 = np.log10(2*np.power(10,fakex))
        fakex3 = np.log10(3*np.power(10,fakex))
        fakex4 = np.log10(4*np.power(10,fakex))
    else:
        fakex = np.linspace(axLim[0][0], axLim[0][1],30)
        fakex15 = 1.5*fakex
        fakex2 = 2.*fakex
        fakex3 = 3.*fakex
        fakex4 = 4.*fakex
    
    for af in [f1, f2]:
        af.plot(fakex, fakex15, 'k--',  alpha=0.4, zorder=0)
        af.plot(fakex, fakex2, 'k--',  alpha=0.3, zorder=0)
        af.plot(fakex, fakex3, 'k--',  alpha=0.2, zorder=0)
        af.plot(fakex, fakex4, 'k--',  alpha=0.1, zorder=0)
        af.plot(fakex15, fakex, 'k--',  alpha=0.4, zorder=0)
        af.plot(fakex2, fakex, 'k--',  alpha=0.3, zorder=0)
        af.plot(fakex3, fakex, 'k--',  alpha=0.2, zorder=0)
        af.plot(fakex4, fakex, 'k--',  alpha=0.1, zorder=0)
    
    f1.set_aspect('equal')
    f1.set_xlim(axLim[0][0], axLim[0][1])
    f1.set_ylim(axLim[1][0], axLim[1][1])
    f1.set_title(pLabs[0])
    f2.set_title(pLabs[1])
    
    f1.set_ylabel(axLabs[1])
    f1.set_xlabel(axLabs[0])
    f2.set_xlabel(axLabs[0])
    
    difRat1 = (p1-p2) / (0.5*(p1+p2))
    f3a.hist(difRat1, bins=bins[0], density=True, color=colors[0], ec='k')

    difRat2 = (p3-p4) / (0.5*(p3+p4))
    f3b.hist(difRat2, bins=bins[1], density=True, color=colors[1], ec='k')
    
    xA, xB = np.power(10,CORmassesA), np.power(10,CORmassesB)
    
    x = np.linspace(bins[0][0], bins[0][-1], 1000)
    ae3a, loce3a, scalee3a = stats.skewnorm.fit(difRat1)
    p3a = stats.skewnorm.pdf(x, ae3a, loce3a, scalee3a)
    Za = x[np.where(p3a == np.max(p3a))[0]][0]
    m0a  = np.sum(x*p3a) / np.sum(p3a)
    
    ae3b, loce3b, scalee3b = stats.skewnorm.fit(difRat2)
    p3b = stats.skewnorm.pdf(x, ae3b, loce3b, scalee3b)
    Zb = x[np.where(p3b == np.max(p3b))[0]][0]
    m0b  = np.sum(x*p3b) / np.sum(p3b)
    
    f3a.plot(x, p3a, '--', c='k', linewidth=3, zorder=3)
    f3b.plot(x, p3b, '--', c='k', linewidth=3, zorder=3)
    f3b.plot(x, p3a, '--', c=colors[0], linewidth=3, zorder=2)
    f3a.plot(x, p3b, '--', c=colors[1], linewidth=3, zorder=2)
    
    f3a.text(0.05, 0.9, ' N: '+str(len(difRat1)), horizontalalignment='left', transform = f3a.transAxes, fontsize=12)
    f3a.text(0.05, 0.82, '$\\mu_f$: '+ '{:.2f}'.format(loce3a), horizontalalignment='left', transform = f3a.transAxes, fontsize=12)
    f3a.text(0.05, 0.74, '$\\alpha_f$: '+ '{:.2f}'.format(ae3a), horizontalalignment='left', transform = f3a.transAxes, fontsize=12)
    f3a.text(0.05, 0.66, '$\\sigma_f$: '+ '{:.2f}'.format(scalee3a), horizontalalignment='left', transform = f3a.transAxes, fontsize=12)
    f3a.text(0.99, 0.82, '$\\mu$: '+ '{:.2f}'.format(m0a), horizontalalignment='right', transform = f3a.transAxes, fontsize=12)
    f3a.text(0.99, 0.74, 'Z: '+ '{:.2f}'.format(Za), horizontalalignment='right', transform = f3a.transAxes, fontsize=12)
    
    f3b.text(0.05, 0.9, ' N: '+str(len(difRat2)), horizontalalignment='left', transform = f3b.transAxes, fontsize=12)
    f3b.text(0.05, 0.82, '$\\mu_F$: '+ '{:.2f}'.format(loce3b), horizontalalignment='left', transform = f3b.transAxes, fontsize=12)
    f3b.text(0.05, 0.74, '$\\alpha_F$: '+ '{:.2f}'.format(ae3b), horizontalalignment='left', transform = f3b.transAxes, fontsize=12)
    f3b.text(0.05, 0.66, '$\\sigma_F$: '+ '{:.2f}'.format(scalee3b), horizontalalignment='left', transform = f3b.transAxes, fontsize=12)
    f3b.text(0.99, 0.82, '$\\mu$: '+ '{:.2f}'.format(m0b), horizontalalignment='right', transform = f3b.transAxes, fontsize=12)
    f3b.text(0.99, 0.74, 'Z: '+ '{:.2f}'.format(Zb), horizontalalignment='right', transform = f3b.transAxes, fontsize=12)
    
    
    yl = f3a.get_ylim()
    f3a.plot([0,0], yl, 'k--', alpha=0.75, zorder=0)
    f3b.plot([0,0], yl, 'k--', alpha=0.75, zorder=0)
    f3a.set_xlim(-1.5, 1.5)
    f3a.set_ylim(yl)
    f3a.set_xlabel(axLabs[2])
    f3b.set_xlabel(axLabs[2])
    f3a.set_ylabel('Prob. Density')
    
    # put lines corresponding to the others
    # mathed out that factors of 1.5, 2, 3, 4 correspond to 
    # [0.4, 0.67, 1, 1.2] for diff/avg
    yl = f3a.get_ylim()
    for af in [f3a, f3b]:
        af.plot([0.4,0.4], yl, 'k--',  alpha=0.4, zorder=0)
        af.plot([0.67,0.67], yl, 'k--',  alpha=0.3, zorder=0)
        af.plot([1.,1.], yl, 'k--',  alpha=0.2, zorder=0)
        af.plot([1.2,1.2], yl, 'k--',  alpha=0.1, zorder=0)
        af.plot([-0.4,-0.4], yl, 'k--',  alpha=0.4, zorder=0)
        af.plot([-0.67,-0.67], yl, 'k--',  alpha=0.3, zorder=0)
        af.plot([-1.,-1.], yl, 'k--',  alpha=0.2, zorder=0)
        af.plot([-1.2,-1.2], yl, 'k--',  alpha=0.1, zorder=0)
        
    
    f1.text(0.95, 0.21, 'Equal', horizontalalignment='right', transform = f1.transAxes, fontsize=12)
    f1.text(0.95, 0.17, '1.5x', horizontalalignment='right', transform = f1.transAxes, fontsize=12, alpha=0.4)
    f1.text(0.95, 0.13, '2x', horizontalalignment='right', transform = f1.transAxes, fontsize=12, alpha=0.3)
    f1.text(0.95, 0.09, '3x', horizontalalignment='right', transform = f1.transAxes, fontsize=12, alpha=0.2)
    f1.text(0.95, 0.05, '4x', horizontalalignment='right', transform = f1.transAxes, fontsize=12, alpha=0.1)        
    
    #plt.show()
    plt.savefig(picName)
    plt.close()
    
    

# |---------------------------------------|
# |---------- Basic setup stuff ----------|
# |---------------------------------------|
def setItUp():
    # import the data structure
    f =  open('allRes.pkl', 'rb')
    res, id2time, time2time = pickle.load(f)
    f.close() 


    # Make key arrays for each different set (CORSET, GCS, both)
    global CORkeys, GCSkeys, bothKeys, allKeys
    CORkeys  = []
    GCSkeys  = []
    bothkeys = []
    allKeys  = np.sort(np.array([key for key in res.keys()]))

    for key in allKeys:
        myRes = res[key]
        hasCOR, hasGCS, hasBoth = False, False, False
        if myRes.CORSETmassMA or myRes.CORSETmassMB:
            CORkeys.append(key)
        if (type(myRes.GCSmassesA) != type(None)) or (type(myRes.GCSmassesB) != type(None)):  
            GCSkeys.append(key) 
        if (type(myRes.GCSprofMassA) != type(None)) or (type(myRes.GCSprofMassA) != type(None)):  
            bothkeys.append(key) 

    CORkeys = np.sort(np.array(CORkeys))    
    GCSkeys = np.sort(np.array(GCSkeys))    
    bothkeys = np.sort(np.array(bothkeys))   


    # |---------------------------------------|
    # |--- Collect the data for full sets ----|
    # |---------------------------------------|
    global CORmassesA, CORmassesB, CORdatesA, CORdatesB, CORpixA, CORpixB 
    CORmassesA, CORmassesB = [], []
    CORpixA, CORpixB = [], []
    CORdatesA, CORdatesB   = [], []
    for key in CORkeys:
        thisRes = res[key]
        if thisRes.CORSETmassMA:
            if thisRes.CORSETmassMA > 0:
                CORmassesA.append(thisRes.CORSETmassMA)
                CORdatesA.append(thisRes.CMEtimeDT)
                # Might have multiple times at max mass, take earliest
                pixIdx = np.where(thisRes.CORSETmassesA == thisRes.CORSETmassMA)[0]
                myPix =  thisRes.CORSETpixsA[pixIdx[0]]
                CORpixA.append(myPix)
        if thisRes.CORSETmassMB:
            if thisRes.CORSETmassMB > 0:
                CORmassesB.append(thisRes.CORSETmassMB)
                CORdatesB.append(thisRes.CMEtimeDT)
                pixIdx = np.where(thisRes.CORSETmassesB == thisRes.CORSETmassMB)[0]
                myPix =  thisRes.CORSETpixsB[pixIdx[0]]
                CORpixB.append(myPix)
   
    
    global GCSmassesA, GCSmassesB,GCSmassesMA, GCSmassesMB, GCSdatesMA, GCSdatesMB, GCSdatesA, GCSdatesB, GCSpixA, GCSpixB  
    GCSmassesA, GCSmassesB = [], []
    GCSmassesMA, GCSmassesMB = [], []
    GCSpixA, GCSpixB = [], []
    GCSdatesA, GCSdatesB   = [], []
    GCSdatesMA, GCSdatesMB   = [], []
    for key in GCSkeys:
        thisRes = res[key]
        if thisRes.GCSmassesA:
            theseM = np.array(thisRes.GCSmassesA)
            goodMass = np.where(theseM > 0)[0]
            if len(goodMass) > 0:
                myMassA  = np.max(theseM[goodMass])
                GCSmassesMA.append(myMassA)
                GCSdatesMA.append(thisRes.CMEtimeDT)
                for idx in goodMass:
                    myMassA  = theseM[idx]
                    GCSdatesA.append(thisRes.CMEtimeDT)
                    GCSmassesA.append(myMassA)
                    GCSpixA.append(thisRes.GCSpixsA[idx])
                
        if thisRes.GCSmassesB:
            theseM = np.array(thisRes.GCSmassesB)
            goodMass = np.where(theseM > 0)[0]
            if len(goodMass) > 0:
                myMassB  = np.max(theseM[goodMass])
                GCSmassesMB.append(myMassB)
                GCSdatesMB.append(thisRes.CMEtimeDT)
                for idx in goodMass:
                    myMassB  = theseM[idx]
                    GCSdatesB.append(thisRes.CMEtimeDT)
                    GCSmassesB.append(myMassB)
                    GCSpixB.append(thisRes.GCSpixsB[idx])
    
    
    # |---------------------------------------|
    # |-- Collect the data for paired sets ---|
    # |---------------------------------------|
    global pCORmassesA, pCORmassesB, pCORpixA, pCORpixB, pCORvelsA, pCORvelsB, pCORdates 
    pCORmassesA, pCORmassesB = [], []
    pCORpixA, pCORpixB = [], []
    pCORvelsA, pCORvelsB = [], []
    pCORdates = []

    for key in CORkeys:
        thisRes = res[key]
        if (type(thisRes.CORSETmassMA) != type(None)) & (type(thisRes.CORSETmassMB) != type(None)):
            if (thisRes.CORSETmassMA > 0) & (thisRes.CORSETmassMB > 0):
                pCORmassesA.append(thisRes.CORSETmassMA)
                # Might have multiple times at max mass, take earliest
                pixIdx = np.where(thisRes.CORSETmassesA == thisRes.CORSETmassMA)[0]
                myPix =  thisRes.CORSETpixsA[pixIdx[0]]
                pCORpixA.append(myPix)
                pCORvelsA.append(thisRes.CORSETvelMA)
                pCORdates.append(thisRes.CMEtimeDT)
                
                pCORmassesB.append(thisRes.CORSETmassMB)
                pixIdx = np.where(thisRes.CORSETmassesB == thisRes.CORSETmassMB)[0]
                myPix =  thisRes.CORSETpixsB[pixIdx[0]]
                pCORpixB.append(myPix)
                pCORvelsB.append(thisRes.CORSETvelMB)
                
    global pGCSmassesA, pGCSmassesB, pGCSpixA, pGCSpixB, pGCSvels, pGCSidx, pGCSdates
    pGCSmassesA, pGCSmassesB = [], []
    pGCSpixA, pGCSpixB = [], []
    pGCSvels, pGCSidx = [], []
    pGCSdates = []
    counter = 0
    for key in GCSkeys:
        thisRes = res[key]
        myKeys = np.array([key for key in thisRes.GCSvals.keys()])
        for i in range(len(thisRes.GCSmassesA)):
            mA, mB = thisRes.GCSmassesA[i], thisRes.GCSmassesB[i]
            myVel = thisRes.GCSvals[myKeys[i]][-1]
            if (mA > 0)  & (mB > 0):
                pGCSmassesA.append(mA)
                pGCSmassesB.append(mB)
                pGCSpixA.append(thisRes.GCSpixsA[i])
                pGCSpixB.append(thisRes.GCSpixsB[i])
                pGCSdates.append(thisRes.CMEtimeDT)
                if myVel != 'None':
                    pGCSvels.append(float(myVel))
                    pGCSidx.append(counter)
                counter += 1
     
    global cgGCSmassesA, cgGCSmassesB, cgCORmassesA, cgCORmassesB, cgvelsA, cgvelsB            
    cgGCSmassesA, cgGCSmassesB = [], []
    cgCORmassesA, cgCORmassesB = [], []
    cgvelsA, cgvelsB = [], []
    for key in GCSkeys:
        thisRes = res[key]
        if thisRes.GCSmassesA:
            theseM = np.array(thisRes.GCSmassesA)
            goodMass = np.where(theseM > 0)[0]
            if len(goodMass) > 0:
                if thisRes.CORSETmassMA:
                    CORmassA = float(thisRes.CORSETmassMA)
                    if CORmassA > 0:
                        for myMassA in theseM[goodMass]:
                            cgGCSmassesA.append(myMassA)
                            cgCORmassesA.append(CORmassA)
                            cgvelsA.append(thisRes.CORSETvelMA)
        if thisRes.GCSmassesB:
            theseM = np.array(thisRes.GCSmassesB)
            goodMass = np.where(theseM > 0)[0]
            if len(goodMass) > 0:
                if thisRes.CORSETmassMB:
                    CORmassB = float(thisRes.CORSETmassMB)
                    if CORmassB > 0:
                        for myMassB in theseM[goodMass]:
                            cgGCSmassesB.append(np.min([myMassB,20]))
                            cgCORmassesB.append(CORmassB)
                            cgvelsB.append(thisRes.CORSETvelMB)
    
                
    # |---------------------------------------|
    # |-- Collect the data for deproj sets ---|
    # |---------------------------------------|
    global dpTimes, dpMasses, dpHeightsA, dpHeightsB, dpLons
    global dpSepA, dpSepB, dpPMA, dpPHA, dpPMB, dpPHB, dpGCSlons, dpLatsA, dpLatsB 
    dpTimes, dpMasses, dpHeightsA, dpHeightsB, dpLons, dpLatsA, dpLatsB = [], [], [], [], [], [], []
    dpSepA, dpSepB, dpPMA, dpPHA, dpPMB, dpPHB, dpLons  = [], [], [], [], [], [], []
    for key in CORkeys:
        thisRes = res[key]
        if thisRes.deprojTimes[0]:
            dts = []
            for j in range(len(thisRes.deprojTimes)):
                dts.append(datetime.datetime.strptime(thisRes.deprojTimes[j], "%Y%m%d_%H%M%S" ))
            dpTimes.append(thisRes.deprojTimes)
            dpMasses.append(thisRes.deprojMasses)
            dpHeightsA.append(thisRes.deprojHeightsA)
            dpHeightsB.append(thisRes.deprojHeightsB)
            dpLons.append(thisRes.deprojLons)
            dpSepA.append(thisRes.deprojSepA)
            dpSepB.append(thisRes.deprojSepB)
            dpPMA.append(thisRes.projMassesA)
            dpPMB.append(thisRes.projMassesB)
            dpPHA.append(thisRes.projHeightsA)
            dpPHB.append(thisRes.projHeightsB)
            dpLatsA.append(thisRes.deprojLatA)
            dpLatsB.append(thisRes.deprojLatB)
            
            if len(thisRes.GCSvals) > 0:
                theseGLons = []
                for name in thisRes.GCSvals.keys():
                    theseGLons.append(float(thisRes.GCSvals[name][1]))
                dpLons.append(theseGLons)
            else:
                dpLons.append([None])
            
    # |---------------------------------------|
    # |--- Collect the set w/deproj & GCS ----|
    # |---------------------------------------|
    global reconLonC, reconLatCA, reconLatCB, reconLonG, reconLatG, recont,reconErrA, reconErrB
    reconLonC, reconLatCA, reconLatCB, reconLonG, reconLatG, recont = [], [], [], [], [], []
    reconErrA, reconErrB = [], []
    
    for key in GCSkeys:
        thisRes = res[key]
        if thisRes.deprojTimes[0]:
            # CORSET values
            maxIdxC = np.where(thisRes.deprojMasses == np.max(thisRes.deprojMasses))[0]
            # Longitude things
            maxLonC = thisRes.deprojLons[maxIdxC[0]]
            # Move it so the lon at max mass is at 0 then move all to +-180
            # and take std. Easiest way to account for circular range
            diffLonsC = thisRes.deprojLons - maxLonC
            diffLonsC[np.where(diffLonsC > 180)] -= 360
            diffLonsC[np.where(diffLonsC < -180)] += 360
            stdLonC = np.std(diffLonsC)
            meanDelLonC = np.mean(diffLonsC)
            meanLonC = maxLonC + meanDelLonC
            sortLons = np.sort(diffLonsC) + maxLonC
            lonStuff = [maxLonC, meanLonC, sortLons[0], sortLons[-1], stdLonC]
                        
            # Latitude things - don't need to worry about circles
            if len(thisRes.deprojLatA) > 1:
                sortLatsA = np.sort(thisRes.deprojLatA)
            else:
                sortLatsA = [thisRes.deprojLatA]
            latStuffA = [thisRes.deprojLatA[maxIdxC[0]], np.mean(thisRes.deprojLatA), sortLatsA[0], sortLatsA[-1], np.std(sortLatsA)]
            # Latitude things - don't need to worry about circles
            if len(thisRes.deprojLatB) > 1:
                sortLatsB = np.sort(thisRes.deprojLatB)
            else:
                sortLatsB = [thisRes.deprojLatB]
            latStuffB = [thisRes.deprojLatB[maxIdxC[0]], np.mean(thisRes.deprojLatB), sortLatsB[0], sortLatsB[-1], np.std(sortLatsB)]
            
            
            # GCS values
            lonsG = []
            latsG = []
            for gKey in thisRes.GCSvals.keys():
                lonsG.append((float(thisRes.GCSvals[gKey][1])+360)%360)
                latsG.append(float(thisRes.GCSvals[gKey][0]))
            # Need to check if -180 to 180 better than 0 to 360    
            lonsAlt = np.copy(lonsG)
            lonsAlt[np.where(lonsAlt > 180)] -= 360
            if np.std(lonsAlt) < np.std(lonsG):
                lonsG = lonsAlt
            if len(lonsG) > 1:
                lonsG = np.array(lonsG)
                medLonG = np.median(lonsG)
                shiftLons = lonsG - medLonG
                shiftLons[np.where(shiftLons > 180)] -= 360
                shiftLons[np.where(shiftLons < -180)] += 360
                meanLonG = np.mean(shiftLons) + medLonG
                stdLonG  = np.std(shiftLons)
                sortLons = np.sort(shiftLons) + medLonG
                lonStuffG = [np.mean(lonsG), sortLons[0], sortLons[-1], np.std(sortLons)]
                
                latsG = np.array(latsG)
                sortLats = np.sort(latsG)
                latStuffG = [np.mean(latsG), sortLats[0], sortLats[-1], np.std(sortLats)]
                          
            else:
                lonStuffG = [lonsG[0], lonsG[0], lonsG[0], 0]
                latStuffG = [latsG[0], latsG[0], latsG[0], 0]
            
            recont.append(key)
            reconLonC.append(lonStuff)
            reconLatCA.append(latStuffA)
            reconLatCB.append(latStuffB)
            reconLonG.append(lonStuffG)
            reconLatG.append(latStuffG)
            reconErrA.append(thisRes.projCPAstdA)
            reconErrB.append(thisRes.projCPAstdB)
            
            
            
            
  
setItUp()        
                
# |---------------------------------------|
# |------ Timeline/Histogram figure ------|
# |---------------------------------------|
if plotAll:  
    p1 = [np.array(CORdatesA), np.array(CORmassesA)*1e15] 
    p2 = [np.array(CORdatesB), np.array(CORmassesB)*1e15] 
    p3 = [np.array(GCSdatesMA), np.array(GCSmassesMA)*1e15] 
    p4 = [np.array(GCSdatesMB), np.array(GCSmassesMB)*1e15] 
    axLabs = ['Mass (g)', 'Mass (g)']
    pLabs  = ['CORSET STA', 'CORSET STB', 'GCS STA', 'GCS STB']
    colors = [cCORA, cCORB, cGCSA, cGCSB]
    b1 = np.arange(13.5,16.7, 0.2, dtype=float)
    b1 = np.append(b1, 2*b1[-1] - b1[-2])
    b2 = np.arange(13.5,16.7,0.2, dtype=float)
    b2 = np.append(b2, 2*b2[-1] - b2[-2])
    bins = [b1, b2]
    
    time2histo4(p1, p2, p3, p4, axLabs, picName=myFold+'timelineHisto.png', pLabs=pLabs, timeYL=[[1e13,8e16],[1e13,8e16]], colors=colors, bins=bins )


# |---------------------------------------|
# |------ COR Pix Time/Hist figure -------|
# |---------------------------------------|
if plotAll:
    p1 = [np.array(CORdatesA), np.array(CORpixA)] 
    p2 = [np.array(CORdatesB), np.array(CORpixB)] 
    p3 = [np.array(CORdatesA), np.array(CORmassesA)*1e15 / np.array(CORpixA)] 
    p4 = [np.array(CORdatesB), np.array(CORmassesB)*1e15 / np.array(CORpixB)] 
    axLabs = ['Pixels', 'Mass /Pix (g)']
    pLabs  = ['CORSET STA', 'CORSET STB', '', '']
    colors = [cCORA, cCORB, cCORA, cCORB]
    b1 = np.arange(3,6.1, 0.2, dtype=float)
    b1 = np.append(b1, 2*b1[-1] - b1[-2])
    b2 = np.arange(9.8,11.4,0.1, dtype=float)
    b2 = np.append(b2, 2*b2[-1] - b2[-2])
    bins = [b1, b2]
    time2histo4(p1, p2, p3, p4, axLabs, picName=myFold+'pixHistoCOR.png', pLabs=pLabs, timeYL=[[7e2,9e5], [2e9,7e11]], colors=colors, bins=bins )
      

# |---------------------------------------|
# |------ GCS Pix Time/Hist figure -------|
# |---------------------------------------|
if plotAll:
    p1 = [np.array(GCSdatesA), np.array(GCSpixA)] 
    p2 = [np.array(GCSdatesB), np.array(GCSpixB)] 
    p3 = [np.array(GCSdatesA), np.array(GCSmassesA)*1e15 / np.array(GCSpixA)] 
    p4 = [np.array(GCSdatesB), np.array(GCSmassesB)*1e15 / np.array(GCSpixB)] 
    axLabs = ['Pixels', 'Mass /Pix (g)']
    pLabs  = ['GCS STA', 'GCS STB', '', '']
    colors = [cGCSA, cGCSB, cGCSA, cGCSB]
    b1 = np.arange(3,6.1, 0.2, dtype=float)
    b1 = np.append(b1, 2*b1[-1] - b1[-2])
    b2 = np.arange(9.8,11.4,0.1, dtype=float)
    b2 = np.append(b2, 2*b2[-1] - b2[-2])
    bins = [b1, b2]
    time2histo4(p1, p2, p3, p4, axLabs, picName=myFold+'pixHistoGCS.png', pLabs=pLabs, timeYL=[[7e2,9e5], [2e9,7e11]], colors=colors, bins=bins )


# |---------------------------------------|
# |------ COR Paired Mass/Pix figure -----|
# |---------------------------------------|
if plotAll:
    p3 = [np.array(pCORdates), np.array(pCORpixA)] 
    p4 = [np.array(pCORdates), np.array(pCORpixB)] 
    p1 = [np.array(pCORdates), np.array(pCORmassesA)*1e15] 
    p2 = [np.array(pCORdates), np.array(pCORmassesB)*1e15] 
    axLabs = ['Mass (g)', 'Pixels']
    pLabs  = ['CORSET STA', 'CORSET STB', '', '']
    colors = [cCORA, cCORB, cCORA, cCORB]
    b1 = np.arange(13.5,16.7, 0.2, dtype=float)
    b1 = np.append(b1, 2*b1[-1] - b1[-2])
    b2 = np.arange(3,6.1, 0.2, dtype=float)
    b2 = np.append(b2, 2*b2[-1] - b2[-2])
    bins = [b1, b2]
    time2histo4(p1, p2, p3, p4, axLabs, picName=myFold+'timelineHistoPCOR.png', pLabs=pLabs, timeYL=[[1e13,8e16],[7e2,9e5]], colors=colors, bins=bins )



# |---------------------------------------|
# |-------------- A2B Comp ---------------|
# |---------------------------------------|
if plotAll:
    p1 = np.array(pCORmassesA)*1e15
    p2 = np.array(pCORmassesB)*1e15
    p3 = np.array(pGCSmassesA)*1e15 
    p4 = np.array(pGCSmassesB)*1e15 
    c1 = [pCORvelsA, None]
    c2 = [pGCSvels, pGCSidx]
    axLabs = ['log Mass STA (g)', 'log Mass STB (g)', '(M$_A$ - M$_B$) / <M>']
    pLabs  = ['CORSET', 'GCS']
    colors = [cCORB, cGCSB]
    b1 = np.arange(-1.5,1.5, 0.2, dtype=float)
    b1 = np.append(b1, 2*b1[-1] - b1[-2])
    b2 = np.arange(-1.5,1.5, 0.2, dtype=float)
    b2 = np.append(b2, 2*b2[-1] - b2[-2])
    bins = [b1, b2]
    
    scatter2histo2(p1, p2, p3, p4, axLabs, picName=myFold+'A2Bcomp.png', pLabs=pLabs, axLim=[[13.5,16.7],[13.5,16.7]], colors=colors, scatCols=[c1,c2,[400,1000],'v\n(km/s)\n'], bins=bins )
 
    
# |---------------------------------------|
# |----------- COR 2 GCS Comp ------------|
# |---------------------------------------|
if plotAll:
    p1 = np.array(cgCORmassesA)*1e15
    p2 = np.array(cgGCSmassesA)*1e15
    p3 = np.array(cgCORmassesB)*1e15 
    p4 = np.array(cgGCSmassesB)*1e15 
    c1 = [cgvelsA, None]
    c2 = [cgvelsB, None]
    axLabs = ['CORSET log Mass (g)', 'GCS log Mass (g)', '(M$_{COR}$ - M$_{GCS}$) / <M>']
    pLabs  = ['STEREO A', 'STEREOB']
    colors = [cGCSA, cGCSB]
    b1 = np.arange(-1.5,1.5, 0.2, dtype=float)
    b1 = np.append(b1, 2*b1[-1] - b1[-2])
    b2 = np.arange(-1.5,1.5, 0.2, dtype=float)
    b2 = np.append(b2, 2*b2[-1] - b2[-2])
    bins = [b1, b2]
    
    scatter2histo2(p1, p2, p3, p4, axLabs, picName=myFold+'COR2GCScomp.png', pLabs=pLabs, axLim=[[14.5,16.5],[14.5,16.5]], colors=colors, scatCols=[c1,c2,[400,1000],'COR v\n(km/s)\n'], bins=bins )

# |---------------------------------------|
# |------------ A2B Comp Pix -------------|
# |---------------------------------------|
if plotAll:
    p1 = np.array(pCORpixA)
    p2 = np.array(pCORpixB)
    p3 = np.array(pGCSpixA)
    p4 = np.array(pGCSpixB)
    c1 = [pCORvelsA, None]
    c2 = [pGCSvels, pGCSidx]
    axLabs = ['log Pix STA', 'log Pix STB', '(Pix$_A$ - Pix$_B$) / <Pix>']
    pLabs  = ['CORSET', 'GCS']
    colors = [cCORB, cGCSB]
    b1 = np.arange(-1.5,1.5, 0.2, dtype=float)
    b1 = np.append(b1, 2*b1[-1] - b1[-2])
    b2 = np.arange(-1.5,1.5, 0.2, dtype=float)
    b2 = np.append(b2, 2*b2[-1] - b2[-2])
    bins = [b1, b2]
    
    scatter2histo2(p1, p2, p3, p4, axLabs, picName=myFold+'A2BpixComp.png', pLabs=pLabs, axLim=[[3,6],[3,6]], colors=colors, scatCols=[c1,c2,[400,1000],'v\n(km/s)\n'], bins=bins )

# |---------------------------------------|
# |--------- A2B Comp Mass/Pix -----------|
# |---------------------------------------|
if plotAll:    
    p1 = np.array(pCORmassesA)*1e15/np.array(pCORpixA)
    p2 = np.array(pCORmassesB)*1e15/np.array(pCORpixB)
    p3 = np.array(pGCSmassesA)*1e15/np.array(pGCSpixA)
    p4 = np.array(pGCSmassesB)*1e15/np.array(pGCSpixB)
    c1 = [pCORvelsA, None]
    c2 = [pGCSvels, pGCSidx]
    axLabs = ['log Mass/Pix STA (g)', 'log Mass/Pix STB', '(M/P$_A$ - M/P$_B$) / <M/P>']
    pLabs  = ['CORSET', 'GCS']
    colors = [cCORB, cGCSB]
    b1 = np.arange(-1.5,1.5, 0.2, dtype=float)
    b1 = np.append(b1, 2*b1[-1] - b1[-2])
    b2 = np.arange(-1.5,1.5, 0.2, dtype=float)
    b2 = np.append(b2, 2*b2[-1] - b2[-2])
    bins = [b1, b2]
    
    scatter2histo2(p1, p2, p3, p4, axLabs, picName=myFold+'A2BmppComp.png', pLabs=pLabs, axLim=[[9.5,11.5],[9.5,11.5]], colors=colors, scatCols=[c1,c2,[400,1000],'v\n(km/s)\n'], bins=bins )


# |---------------------------------------|
# |---------- deProj Mass Comp -----------|
# |---------------------------------------|
if plotAll:
    trueM = np.array([np.max(x) for x in dpMasses])
    projMA = np.array([np.max(x) for x in dpPMA])
    projMB = np.array([np.max(x) for x in dpPMB])
    sepA = []
    sepB = []
    hRatioA = []
    hRatioB = []
    times = []
    for i in range(len(trueM)):
        theseM = dpMasses[i]
        idx = np.where(theseM == trueM[i])[0]
        sepA.append(dpSepA[i][idx[0]])
        sepB.append(dpSepB[i][idx[0]])
        hRatioA.append(np.abs(dpPHA[i][idx[0]] / dpHeightsA[i][idx[0]] ))
        hRatioB.append(np.abs(dpPHB[i][idx[0]] / dpHeightsB[i][idx[0]]  ))
        times.append(datetime.datetime.strptime(dpTimes[i][idx[0]], "%Y%m%d_%H%M%S" ))

    MratioA = projMA / trueM
    MratioB = projMB / trueM
    p1 = [times, MratioA] 
    p2 = [times, MratioB] 
    p3 = [times, hRatioA] 
    p4 = [times, hRatioB] 
    axLabs = ['Proj / Deproj Mass', 'Proj / Deproj Height']
    pLabs  = ['COR STA', 'COR STB', '', '']
    colors = [cCORA, cCORB, cCORA, cCORB]
    b1 = np.arange(0, 1., 0.1, dtype=float)
    b1 = np.append(b1, 2*b1[-1] - b1[-2])
    bins = [b1, b1]
    time2histo4(p1, p2, p3, p4, axLabs, picName=myFold+'deProjMasscomp.png', pLabs=pLabs, timeYL=[[0,1.3], [0,1.3]], colors=colors, bins=bins, doLog=False, plotFits=False)
        


# |---------------------------------------|
# |------- deProj Position Comp ----------|
# |---------------------------------------|
if True:
    
    fig = plt.figure(constrained_layout=True, figsize=(10,5))
    gs = fig.add_gridspec(2, 4)
    f1 = fig.add_subplot(gs[0,0])  
    f2 = fig.add_subplot(gs[0,1])
    f3 = fig.add_subplot(gs[0,2], sharex=f2, sharey=f2)
    f4 = fig.add_subplot(gs[0,3], sharex=f2, sharey=f2)
    f1b =  fig.add_subplot(gs[1,0])  
    f2b =  fig.add_subplot(gs[1,1])  
    f3b =  fig.add_subplot(gs[1,2], sharex=f2b)  
    f4b =  fig.add_subplot(gs[1,3], sharex=f2b)  
    
    
    lonsC = np.array([x[0] for x in reconLonC])
    latsA = np.array([x[0] for x in reconLatCA])
    latsB = np.array([x[0] for x in reconLatCB])
    lonsG = np.array([x[0] for x in reconLonG])
    latsG = np.array([x[0] for x in reconLatG])
    
    lonDifs = np.abs(lonsG - lonsC)
    for i in range(len(lonDifs)):
        if lonDifs[i] > 180:
            if lonsC[i] < 180:
                lonsG[i] -= 360
            else:
                lonsG[i] += 360
            lonDifs[i] = np.abs(lonsC[i] - lonsG[i])
    lonDifs = lonsC - lonsG
    f1b.hist(lonDifs)
    
    f1.scatter(lonsC, lonsG)
    f2.scatter(latsA, latsB)
    f2b.hist(latsA - latsB)
    f3.scatter(latsA, latsG)
    f3b.hist(latsA - latsG)
    f4.scatter(latsB, latsG)
    f4b.hist(latsB - latsG)
    
    fs = [f1, f2, f3, f4]
    doneLon = False
    for f in fs:
        f.set_aspect('equal')
        xl = f.get_xlim()
        yl = f.get_ylim()
        f.plot([-200,600], [-200,600], 'k--', zorder=0, lw=2)
        
        delt = 30
        if doneLon:
            delt = 10
        doneLon = True
        for i in range(3):
            j = i+1
            f.plot([-200,600], [-200+delt*j,600+delt*j], 'k--', zorder=0, alpha=(0.5-0.1*j))
            f.plot([-200,600], [-200-delt*j,600-delt*j], 'k--', zorder=0, alpha=(0.5-0.1*j))
        f.set_xlim(xl)
        f.set_ylim(yl)
        
    f1.set_xlim([-50,400])
    plt.show()
    