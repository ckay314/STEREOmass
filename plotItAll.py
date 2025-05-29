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
    global res
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
            dpHeightsA.append(thisRes.deprojHeightsA.astype(float))
            dpHeightsB.append(thisRes.deprojHeightsB.astype(float))
            dpLons.append(thisRes.deprojLons)
            dpSepA.append(thisRes.deprojSepA)
            dpSepB.append(thisRes.deprojSepB)
            dpPMA.append(thisRes.projMassesA)
            dpPMB.append(thisRes.projMassesB)
            dpPHA.append(thisRes.projHeightsA.astype(float))
            dpPHB.append(thisRes.projHeightsB.astype(float))
            dpLatsA.append(thisRes.deprojLatA)
            dpLatsB.append(thisRes.deprojLatB)
            
            if len(thisRes.GCSvals) > 0:
                theseGLons = []
                for name in thisRes.GCSvals.keys():
                    theseGLons.append(float(thisRes.GCSvals[name][1]))
                dpLons.append(theseGLons)
            else:
                dpLons.append([None])

    global CdpTimes, CdpMassesA, CdpMassesB, CdpHeightsA, CdpHeightsB, CdpLons
    global CdpSepA, CdpSepB, CdpPMA, CdpPHA, CdpPMB, CdpPHB, CdpGCSlons, CdpLatsA, CdpLatsB 
    CdpTimes, CdpMassesA, CdpMassesB, CdpHeightsA, CdpHeightsB, CdpLons, CdpLatsA, CdpLatsB = [], [], [], [], [], [], [], []
    CdpSepA, CdpSepB, CdpPMA, CdpPHA, CdpPMB, CdpPHB, CdpLons  = [], [], [], [], [], [], []
    for key in CORkeys:
        thisRes = res[key]
        if thisRes.CdeprojTimes[0]:
            Cdts = []
            for j in range(len(thisRes.CdeprojTimes)):
                Cdts.append(datetime.datetime.strptime(thisRes.CdeprojTimes[j], "%Y%m%d_%H%M%S" ))
            CdpTimes.append(thisRes.CdeprojTimes)
            CdpMassesA.append(thisRes.CdeprojMassesA)
            CdpMassesB.append(thisRes.CdeprojMassesB)
            CdpHeightsA.append(thisRes.CdeprojHeightsA.astype(float))
            CdpHeightsB.append(thisRes.CdeprojHeightsB.astype(float))
            CdpLons.append(thisRes.CdeprojLons)
            CdpSepA.append(thisRes.CdeprojSepA)
            CdpSepB.append(thisRes.CdeprojSepB)
            CdpPMA.append(thisRes.CprojMassesA)
            CdpPMB.append(thisRes.CprojMassesB)
            CdpPHA.append(thisRes.CprojHeightsA.astype(float))
            CdpPHB.append(thisRes.CprojHeightsB.astype(float))
            CdpLatsA.append(thisRes.CdeprojLatA)
            CdpLatsB.append(thisRes.CdeprojLatB)
            
            if len(thisRes.GCSvals) > 0:
                CtheseGLons = []
                for name in thisRes.GCSvals.keys():
                    CtheseGLons.append(float(thisRes.GCSvals[name][1]))
                CdpLons.append(CtheseGLons)
            else:
                CdpLons.append([None])

            
    # |---------------------------------------|
    # |--- Collect the set w/deproj & GCS ----|
    # |---------------------------------------|
    global compTimes, compLons, compLats, compMasses, compHasGCS
    compTimes, compLons, compLats, compMasses, compHasGCS = [], {}, {}, {}, []
    counter = 0
    setupNames = ['COR', 'GCS',  'DP1', 'DP2' ]
    holders = [compLons, compLats, compMasses]
    for holder in holders:
        for name in setupNames:
            holder[name+'A'] = []
            holder[name+'B'] = []
    
    for key in CORkeys:
        thisRes = res[key]
        # Check if had matched A/B and could deproj
        if thisRes.deprojTimes[0]:
            # add to the time array
            compTimes.append(thisRes.CMEtimeDT)
            
            # Get basic corset values (mass + lat from CPA)
            #compMasses['CORA'].append(thisRes.CORSETmassMA)
            compMasses['CORA'].append(np.max(thisRes.CORSETmassesA))
            compMasses['CORB'].append(thisRes.CORSETmassMB)
            # assume PoS so CPA 
            latA = 90 - thisRes.CORSETcpaMA
            latB = 90 - thisRes.CORSETcpaMB
            # move everyone into -180 to 180
            if latA > 180: latA -= 360
            elif latA < -180: latA += 360
            if latB > 180: latB -= 360
            elif latB < -180: latB += 360
            # move everyone into -90
            if latA > 90: latA = 180 - latA
            elif latA < -90: latA = -180 - latA
            if latB > 90: latB = 180 - latB
            elif latB < -90: latB = -180 - latB
            # put in the holders
            compLats['CORA'].append(latA)
            compLats['CORB'].append(latB)
            # No lons here
            compLons['CORA'].append(None)
            
            
            # Get DP1 (using mass only) values
            # Masses are almost identical so treat A/B same
            maxIdx = np.where(thisRes.deprojMasses == np.max(thisRes.deprojMasses))[0]
            compMasses['DP1A'].append(thisRes.deprojMasses[maxIdx[0]])
            compMasses['DP1B'].append(thisRes.deprojMasses[maxIdx[0]])
            # Lon is exact same as output of DP
            maxLon = thisRes.deprojLons[maxIdx[0]]
            shiftLonDiff = thisRes.deprojLons - maxLon
            shiftLonDiff[np.where(shiftLonDiff > 180)] -= 360
            shiftLonDiff[np.where(shiftLonDiff < -180)] += 360
            newMedLon = np.median(shiftLonDiff) + maxLon
            compLons['DP1A'].append(newMedLon)
            compLats['DP1A'].append(np.median(thisRes.deprojLatA))
            compLats['DP1B'].append(np.median(thisRes.deprojLatB))
            
            
            
            # Get DP2 (mass+lat+height) values
            maxIdxA = np.where(thisRes.CdeprojMassesA == np.max(thisRes.CdeprojMassesA))[0]
            maxIdxB = np.where(thisRes.CdeprojMassesB == np.max(thisRes.CdeprojMassesB))[0]
            if maxIdxA[0] != maxIdxB[0]:
                maxIdx = np.max([maxIdxA[0], maxIdxB[0]])
            else:
                maxIdx = maxIdxA[0]
            compMasses['DP2A'].append(thisRes.CdeprojMassesA[maxIdxA[0]])
            compMasses['DP2B'].append(thisRes.CdeprojMassesB[maxIdxB[0]])
            maxLon = thisRes.CdeprojLons[maxIdx]
            shiftLonDiff = thisRes.CdeprojLons - maxLon
            shiftLonDiff[np.where(shiftLonDiff > 180)] -= 360
            shiftLonDiff[np.where(shiftLonDiff < -180)] += 360
            newMedLon = np.median(shiftLonDiff) + maxLon
            compLons['DP2A'].append(newMedLon)
            
            compLats['DP2A'].append(np.median(thisRes.CdeprojLatA))
            compLats['DP2B'].append(np.median(thisRes.CdeprojLatB))
                
                
            
            # check if has GCS values
            if key in GCSkeys:
                compHasGCS.append(counter)
                # GCS values
                lonsG = []
                latsG = []
                
                for gKey in thisRes.GCSvals.keys():
                    lonsG.append((float(thisRes.GCSvals[gKey][1])+360)%360)
                    latsG.append(float(thisRes.GCSvals[gKey][0]))
                massAG = np.array(thisRes.GCSmassesA)
                massAG = massAG[np.where(massAG > 0)] 
                massBG = np.array(thisRes.GCSmassesB)
                massBG = massBG[np.where(massBG > 0)] 
                    
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
                else:
                    meanLonG = lonsG[0]
                compMasses['GCSA'].append(np.mean(massAG))
                compMasses['GCSB'].append(np.mean(massBG))
                compLons['GCSA'].append(meanLonG)
                compLats['GCSA'].append(np.mean(latsG))
                compLats['GCSB'].append(np.mean(latsG))
            else:
                compMasses['GCSA'].append(None)
                compMasses['GCSB'].append(None)
                compLons['GCSA'].append(None)
                compLats['GCSA'].append(None)
                compLats['GCSB'].append(None)
            #print (compMasses['CORA'][-1], compMasses['CORB'][-1], compMasses['GCSA'][-1], compMasses['GCSB'][-1], compMasses['DP1A'][-1], compMasses['DP1B'][-1], compMasses['DP2A'][-1], compMasses['DP2B'][-1])  
            '''if (compMasses['DP1A'][-1] < compMasses['CORA'][-1]):
                print (key, compMasses['DP1A'][-1], compMasses['CORA'][-1])
                print ('   ', thisRes.depSatLonsA[0], thisRes.depSatLonsB[0])
                for i in range(len(thisRes.projHeightsA)):
                    print ('       ', thisRes.projHeightsA[i], thisRes.projMassesA[i], thisRes.projCPAsA[i],'    ',thisRes.projHeightsB[i], thisRes.projMassesB[i], thisRes.projCPAsB[i])'''
            '''if np.abs(compLats['DP1A'][-1]) == 0:
                print (key, thisRes.projHeightsA[i], thisRes.projMassesA[i], thisRes.projCPAsA[i],'    ',thisRes.projHeightsB[i], thisRes.projMassesB[i], thisRes.projCPAsB[i])
                print ('   ', thisRes.depSatLonsA[0], thisRes.depSatLonsB[0])'''
            counter += 1
    
    for key in compLats:
        compLats[key] = np.array(compLats[key])
        compLons[key] = np.array(compLons[key])
        compMasses[key] = np.array(compMasses[key])
    compTimes = np.array(compTimes)
    compHasGCS = np.array(compHasGCS)
            
            
            
  
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
        #if hRatioB[-1] <0.1:
        #    print(times[-1], hRatioA[-1], hRatioB[-1], sepA[-1], sepB[-1])

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
# |------- deProjCombo Mass Comp ---------|
# |---------------------------------------|
if plotAll:
    trueMA = np.array([np.max(x) for x in CdpMassesA])
    trueMB = np.array([np.max(x) for x in CdpMassesB])
    projMA = np.array([np.max(x) for x in CdpPMA])
    projMB = np.array([np.max(x) for x in CdpPMB])
    sepA = []
    sepB = []
    hRatioA = []
    hRatioB = []
    times = []
    for i in range(len(trueM)):
        theseMA = CdpMassesA[i]
        idxA = np.where(theseMA == trueMA[i])[0]
        theseMB = CdpMassesB[i]
        idxB = np.where(theseMB == trueMB[i])[0]
        sepA.append(CdpSepA[i][idxA[0]])
        sepB.append(CdpSepB[i][idxB[0]])
        hRatioA.append(np.abs(CdpPHA[i][idxA[0]] / CdpHeightsA[i][idxA[0]] ))
        hRatioB.append(np.abs(CdpPHB[i][idxB[0]] / CdpHeightsB[i][idxB[0]]  ))
        times.append(datetime.datetime.strptime(CdpTimes[i][idxA[0]], "%Y%m%d_%H%M%S" ))

    MratioA = projMA / trueMA
    MratioB = projMB / trueMB
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
    time2histo4(p1, p2, p3, p4, axLabs, picName=myFold+'deProjMassCombocomp.png', pLabs=pLabs, timeYL=[[0,1.3], [0,1.3]], colors=colors, bins=bins, doLog=False, plotFits=False)



# |---------------------------------------|
# |--------- MLL Comp Figs (2) -----------|
# |---------------------------------------|
if plotAll:
    # globals =  compTimes, compLons, compLats, compMasses, compHasGCS
    mkr = '.'
    cAB = '#882255'
    cAA = '#332288'
    cMix = '#AA4499'
    c1sat = '#66CCEE'
    
    # Position Plot
    fig = plt.figure(figsize=(10,10))
    gs = fig.add_gridspec(6, 6, hspace=0, wspace=0)
    ff1a = fig.add_subplot(gs[0,0]) 
    ff1b = fig.add_subplot(gs[0,1], sharex=ff1a, sharey=ff1a) 
    ff2a = fig.add_subplot(gs[1,0], sharex=ff1a, sharey=ff1a) 
    
    
    f1a = fig.add_subplot(gs[0,-1]) 
    f2a = fig.add_subplot(gs[1,-1], sharex=f1a, sharey=f1a) 
    f2b = fig.add_subplot(gs[1,-2], sharex=f1a, sharey=f1a) 
    f3a = fig.add_subplot(gs[2,-1], sharex=f1a, sharey=f1a) 
    f3b = fig.add_subplot(gs[2,-2], sharex=f1a, sharey=f1a) 
    f3c = fig.add_subplot(gs[2,-3], sharex=f1a, sharey=f1a) 
    f4a = fig.add_subplot(gs[3,-1], sharex=f1a, sharey=f1a) 
    f4b = fig.add_subplot(gs[3,-2], sharex=f1a, sharey=f1a) 
    f4c = fig.add_subplot(gs[3,-3], sharex=f1a, sharey=f1a) 
    f4d = fig.add_subplot(gs[3,-4], sharex=f1a, sharey=f1a) 
    f5a = fig.add_subplot(gs[4,-1], sharex=f1a, sharey=f1a) 
    f5b = fig.add_subplot(gs[4,-2], sharex=f1a, sharey=f1a) 
    f5c = fig.add_subplot(gs[4,-3], sharex=f1a, sharey=f1a) 
    f5d = fig.add_subplot(gs[4,-4], sharex=f1a, sharey=f1a) 
    f5e = fig.add_subplot(gs[4,-5], sharex=f1a, sharey=f1a) 
    f6a = fig.add_subplot(gs[5,-1], sharex=f1a, sharey=f1a) 
    f6b = fig.add_subplot(gs[5,-2], sharex=f1a, sharey=f1a) 
    f6c = fig.add_subplot(gs[5,-3], sharex=f1a, sharey=f1a) 
    f6d = fig.add_subplot(gs[5,-4], sharex=f1a, sharey=f1a) 
    f6e = fig.add_subplot(gs[5,-5], sharex=f1a, sharey=f1a) 
    f6f = fig.add_subplot(gs[5,-6], sharex=f1a, sharey=f1a) 

    allfs = [f1a, f2a, f2b, f3a, f3b, f3c, f4a, f4b, f4c, f4d, f5a, f5b, f5c, f5d, f5e, f6a, f6b, f6c, f6d, f6e, f6f]
    leftfs = [f1a, f2a, f3a, f4a, f5a]
    botfs  = [f6b, f6c, f6d, f6e, f6f]
    midfs = [f2b, f3b, f3c, f4b, f4c, f4d, f5b, f5c, f5d, f5e]
    edges = [f1a, f2b, f3c, f4d, f5e, f6f]
    
    allffs = [ff1a, ff1b, ff2a]
    
    for f in leftfs:
        f.tick_params(labelbottom=False)
        f.yaxis.set_label_position("right")
        f.yaxis.tick_right()
        f.set_ylabel('Lat ($^{\\circ}$)')
    for f in botfs:
        f.tick_params(labelleft=False)
        f.set_xlabel('Lat ($^{\\circ}$)')
        f.yaxis.tick_right()
    for f in midfs:
        f.tick_params(labelbottom=False, labelleft=False)
        f.yaxis.tick_right()
        
    f6a.yaxis.tick_right()
    f6a.set_xlabel('Lat ($^{\\circ}$)')
    f6a.yaxis.set_label_position("right")
    f6a.set_ylabel('Lat ($^{\\circ}$)')
    
    ff1a.set_ylabel('Lon ($^{\\circ}$)')
    ff2a.set_ylabel('Lon ($^{\\circ}$)')
    ff1a.xaxis.set_label_position("top")
    ff1a.xaxis.tick_top()
    ff1a.tick_params(labelbottom=False)
    ff1a.set_xlabel('Lon ($^{\\circ}$)')
    ff1b.xaxis.set_label_position("top")
    ff1b.xaxis.tick_top()
    ff1b.tick_params(labelbottom=False)
    ff1b.tick_params(labelleft=False)
    ff1b.set_xlabel('Lon ($^{\\circ}$)')
    ff2a.tick_params(labelbottom=False)
    ff2a.xaxis.tick_top()
    
    pairs = [[compLons['DP1A'][compHasGCS], compLons['GCSA'][compHasGCS]], [compLons['DP1A'], compLons['DP2A']], [compLons['DP2A'][compHasGCS], compLons['GCSA'][compHasGCS]]]
    i = 0
    
    ll, ul = -50, 410
    ff1a.set_xlim(ll,ul)
    ff1a.set_ylim(ll,ul)
    ff1a.set_aspect('equal')
    #ff1a.set_xticks([0,20,40])

    for ff in allffs:
        x, y = pairs[i][0], pairs[i][1]
        lonDifs = x - y
        for j in range(len(lonDifs)):
            if np.abs(lonDifs[j]) > 180:
                if x[j] < 180:
                    y[j] -= 360
                else:
                    y[j] += 360
        ff.plot([ll, ul], [ll, ul], 'k--', zorder=0, alpha=0.75)
        ff.plot([ll, ul], [ll+20, ul+20], 'k--', zorder=0, alpha=0.4)
        ff.plot([ll, ul], [ll-20, ul-20], 'k--', zorder=0, alpha=0.4)
        ff.scatter(x,y, marker=mkr, c=c1sat)
        i +=1
    
    ff1b.text(1.5, 0.5, 'Deproj1', transform=ff1b.transAxes, fontsize=12, ha='center', va='center')
    ff2a.text(1.5, 0.5, 'Deproj2', transform=ff2a.transAxes, fontsize=12, ha='center', va='center')
    ff2a.text(0.5, -0.2, 'GCS', transform=ff2a.transAxes, fontsize=12, ha='center', va='center')
    
    
    ll, ul = -70, 70
    f1a.set_xlim(ll,ul)
    f1a.set_ylim(ll,ul)
    f1a.set_aspect('equal')
    for f in allfs:
        f.plot([ll, ul], [ll, ul], 'k--', zorder=0, alpha=0.75)
    
    toptits = ['CORSET B', 'GCS', 'Deproj1 A', 'Deproj1 B' 'Deproj2 A', 'deproj2 B']    
    righttits = ['GCS', 'Deproj1 A', 'Deproj1 B', 'Deproj2 A', 'Deproj2 B', 'CORSET\nA']    
    i = 0
    for f in edges:
        if i != 5:
            f.text(-0.5, 0.5, righttits[i], transform=f.transAxes, fontsize=12, ha='center', va='center')
        else:
            f.text(-0.3, 0.5, righttits[i], transform=f.transAxes, fontsize=12, ha='center', va='center')
        i += 1
    f1a.set_title('CORSET B', fontsize=12)


    f1a.scatter(compLats['CORB'], compLats['GCSA'], marker=mkr, c=c1sat)

    f2a.scatter(compLats['CORB'], compLats['DP1A'], marker=mkr, c=cMix)
    f2b.scatter(compLats['GCSA'], compLats['DP1A'], marker=mkr, c=c1sat)

    f3a.scatter(compLats['CORB'], compLats['DP1B'], marker=mkr, c=cAA)
    f3b.scatter(compLats['GCSA'], compLats['DP1B'], marker=mkr, c=c1sat)
    f3c.scatter(compLats['DP1A'], compLats['DP1B'], marker=mkr, c=cAB)

    f4a.scatter(compLats['CORB'], compLats['DP2A'], marker=mkr, c=cMix)
    f4b.scatter(compLats['GCSA'], compLats['DP2A'], marker=mkr, c=c1sat)
    f4c.scatter(compLats['DP1A'], compLats['DP2A'], marker=mkr, c=cAA)
    f4d.scatter(compLats['DP1B'], compLats['DP2A'], marker=mkr, c=cMix)

    f5a.scatter(compLats['CORB'], compLats['DP2B'], marker=mkr, c=cAA)
    f5b.scatter(compLats['GCSA'], compLats['DP2B'], marker=mkr, c=c1sat)
    f5c.scatter(compLats['DP1A'], compLats['DP2B'], marker=mkr, c=cMix)
    f5d.scatter(compLats['DP1B'], compLats['DP2B'], marker=mkr, c=cAA)
    f5e.scatter(compLats['DP2A'], compLats['DP2B'], marker=mkr, c=cAB)

    f6a.scatter(compLats['CORB'], compLats['CORA'], marker=mkr, c=cAB)
    f6b.scatter(compLats['GCSA'], compLats['CORA'], marker=mkr, c=c1sat)
    f6c.scatter(compLats['DP1A'], compLats['CORA'], marker=mkr, c=cAA)
    f6d.scatter(compLats['DP1B'], compLats['CORA'], marker=mkr, c=cAB)
    f6e.scatter(compLats['DP2A'], compLats['CORA'], marker=mkr, c=cMix)
    f6f.scatter(compLats['DP2B'], compLats['CORA'], marker=mkr, c=cMix)


    f6f.text(0.1, 0.55, 'Same Method, Diff Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cAB, )
    f6f.text(0.1, 0.53, 'Diff Method, Same Sat', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cAA)
    f6f.text(0.1, 0.51, 'Diff Method, Diff Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cMix)
    f6f.text(0.1, 0.49, 'One Value For Both Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=c1sat)
    plt.savefig(myFold+'MegaPosComp.png')
    
    # Mass Plot
    fig = plt.figure(figsize=(10,10))
    gs = fig.add_gridspec(6, 6, hspace=0, wspace=0)
    f1a = fig.add_subplot(gs[0,0]) 
    f2a = fig.add_subplot(gs[1,0], sharex=f1a, sharey=f1a) 
    f2b = fig.add_subplot(gs[1,1], sharex=f1a, sharey=f1a) 
    f3a = fig.add_subplot(gs[2,0], sharex=f1a, sharey=f1a) 
    f3b = fig.add_subplot(gs[2,1], sharex=f1a, sharey=f1a) 
    f3c = fig.add_subplot(gs[2,2], sharex=f1a, sharey=f1a) 
    f4a = fig.add_subplot(gs[3,0], sharex=f1a, sharey=f1a) 
    f4b = fig.add_subplot(gs[3,1], sharex=f1a, sharey=f1a) 
    f4c = fig.add_subplot(gs[3,2], sharex=f1a, sharey=f1a) 
    f4d = fig.add_subplot(gs[3,3], sharex=f1a, sharey=f1a) 
    f5a = fig.add_subplot(gs[4,0], sharex=f1a, sharey=f1a) 
    f5b = fig.add_subplot(gs[4,1], sharex=f1a, sharey=f1a) 
    f5c = fig.add_subplot(gs[4,2], sharex=f1a, sharey=f1a) 
    f5d = fig.add_subplot(gs[4,3], sharex=f1a, sharey=f1a) 
    f5e = fig.add_subplot(gs[4,4], sharex=f1a, sharey=f1a) 
    f6a = fig.add_subplot(gs[5,0], sharex=f1a, sharey=f1a) 
    f6b = fig.add_subplot(gs[5,1], sharex=f1a, sharey=f1a) 
    f6c = fig.add_subplot(gs[5,2], sharex=f1a, sharey=f1a) 
    f6d = fig.add_subplot(gs[5,3], sharex=f1a, sharey=f1a) 
    f6e = fig.add_subplot(gs[5,4], sharex=f1a, sharey=f1a) 
    f6f = fig.add_subplot(gs[5,5], sharex=f1a, sharey=f1a) 

    allfs = [f1a, f2a, f2b, f3a, f3b, f3c, f4a, f4b, f4c, f4d, f5a, f5b, f5c, f5d, f5e, f6a, f6b, f6c, f6d, f6e, f6f]
    leftfs = [f1a, f2a, f3a, f4a, f5a]
    botfs  = [f6b, f6c, f6d, f6e, f6f]
    midfs = [f2b, f3b, f3c, f4b, f4c, f4d, f5b, f5c, f5d, f5e]
    edges = [f1a, f2b, f3c, f4d, f5e, f6f]

    ll, ul = 0, 50
    f1a.set_xlim(ll,ul)
    f1a.set_ylim(ll,ul)
    f1a.set_aspect('equal')
    f1a.set_xticks([0,20,40])
    
    for f in leftfs:
        f.tick_params(labelbottom=False)
        f.set_ylabel('M (10$^{15}$ g)')
    for f in botfs:
        f.tick_params(labelleft=False)
        f.set_xlabel('M (10$^{15}$ g)')
    for f in midfs:
        f.tick_params(labelbottom=False, labelleft=False)
    f6a.set_xlabel('M (10$^{15}$ g)')
    f6a.set_ylabel('M (10$^{15}$ g)')
    
    for f in allfs:
        f.plot([ll, ul], [ll, ul], 'k--', zorder=0, alpha=0.5)
        
    toptits = ['CORSET B', 'GCS A', 'GCS B', 'Deproj1', 'Deproj2 A', 'deproj2 B']    
    righttits = ['GCS A', 'GCSB', 'Deproj1', 'Deproj2 A', 'Deproj2 B', 'CORSET\nA']    
    i = 0
    for f in edges:
        if i != 5:
            f.text(1.5, 0.5, righttits[i], transform=f.transAxes, fontsize=12, ha='center', va='center')
        else:
            f.text(1.3, 0.5, righttits[i], transform=f.transAxes, fontsize=12, ha='center', va='center')
        i += 1
    f1a.set_title('CORSET B', fontsize=12)
    
    mkr = '.'
    cAB = '#882255'
    cAA = '#332288'
    cMix = '#AA4499'
    c1sat = '#66CCEE'
    f1a.scatter(compMasses['CORB'], compMasses['GCSA'], marker=mkr, c=cMix)

    f2a.scatter(compMasses['CORB'], compMasses['GCSB'], marker=mkr, c=cAA)
    f2b.scatter(compMasses['GCSA'], compMasses['GCSB'], marker=mkr, c=cAB)

    f3a.scatter(compMasses['CORB'], compMasses['DP1A'], marker=mkr, c= c1sat)
    f3b.scatter(compMasses['GCSA'], compMasses['DP1A'], marker=mkr, c= c1sat)
    f3c.scatter(compMasses['GCSB'], compMasses['DP1A'], marker=mkr, c= c1sat)

    f4a.scatter(compMasses['CORB'], compMasses['DP2A'], marker=mkr, c=cMix)
    f4b.scatter(compMasses['GCSA'], compMasses['DP2A'], marker=mkr, c=cAA)
    f4c.scatter(compMasses['GCSB'], compMasses['DP2A'], marker=mkr, c=cMix)
    f4d.scatter(compMasses['DP1A'], compMasses['DP2A'], marker=mkr, c= c1sat)

    f5a.scatter(compMasses['CORB'], compMasses['DP2B'], marker=mkr, c=cAA)
    f5b.scatter(compMasses['GCSA'], compMasses['DP2B'], marker=mkr, c=cMix)
    f5c.scatter(compMasses['GCSB'], compMasses['DP2B'], marker=mkr, c=cAA)
    f5d.scatter(compMasses['DP1A'], compMasses['DP2B'], marker=mkr, c= c1sat)
    f5e.scatter(compMasses['DP2A'], compMasses['DP2B'], marker=mkr, c=cAB)

    f6a.scatter(compMasses['CORB'], compMasses['CORA'], marker=mkr, c=cAB)
    f6b.scatter(compMasses['GCSA'], compMasses['CORA'], marker=mkr, c=cAA)
    f6c.scatter(compMasses['GCSB'], compMasses['CORA'], marker=mkr, c=cMix)
    f6d.scatter(compMasses['DP1A'], compMasses['CORA'], marker=mkr, c= c1sat)
    f6e.scatter(compMasses['DP2A'], compMasses['CORA'], marker=mkr, c=cAA)
    f6f.scatter(compMasses['DP2B'], compMasses['CORA'], marker=mkr, c=cMix)
    
    f6f.text(0.6, 0.75, 'Same Method, Diff Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cAB)
    f6f.text(0.6, 0.73, 'Diff Method, Same Sat', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cAA)
    f6f.text(0.6, 0.71, 'Diff Method, Diff Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cMix)
    f6f.text(0.6, 0.69, 'Single Value For Both Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=c1sat)

    plt.savefig(myFold+'MegaMassComp.png')
    
    

# |---------------------------------------|
# |------------ Mass Growth --------------|
# |---------------------------------------|
if True:
    
    for key in CORkeys:
        print (key)
        thisRes = res[key]
        
        # Check if had matched A/B and could deproj
        if thisRes.deprojTimes[0]: 
            fig = plt.figure()
            masses = thisRes.deprojMasses
            massesC = thisRes.CdeprojMassesA
            #times  = thisRes.deprojTimes
            heights = np.abs(thisRes.deprojHeightsA.astype(float))
            heightsC = np.abs(thisRes.CdeprojHeightsA.astype(float))
            goodIdx = np.where((heights >= 1) & (heights < 21.5) & (masses > 0) )[0]
            goodIdxC = np.where((heightsC >= 1) & (heightsC < 21.5) & (massesC > 0) )[0]
            if len(goodIdx) > 0:
                subHeight = heights[goodIdx]
                subMass   = masses[goodIdx]
                sclMass = subMass / np.max(subMass)
            
                #plt.plot(heights[goodIdx], masses[goodIdx])
                plt.plot(subHeight, sclMass)
                print (subHeight)
                plt.show()
            plt.close()
                
        #print (sd)
                
    plt.show()
            
            