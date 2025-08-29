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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from sklearn.metrics import r2_score


plotAll = False
myFold = 'paperFigs/'

cCORA = 'r'
cCORB = 'b'
cGCSA = '#882255'
cGCSB = '#88CCEE'
cLBA  = '#ff66cc'
cLBB  = '#1f3d7a'

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

# Generic scatter/histo 2x2 plot
def scatter2histo2(p1, p2, p3, p4, axLabs, doLog=True, picName='temp.png', pLabs=[None], colors=['r', 'b', 'maroon', 'lightblue'], axLim=[None], bins=[[None]], scatCols=[[None]]):
    fig = plt.figure(constrained_layout=True, figsize=(12, 8.5))
    gs = fig.add_gridspec(3, 2)
    f1 = fig.add_subplot(gs[0:2,0]) 
    f2 = fig.add_subplot(gs[0:2,1], sharex=f1, sharey=f1) 
    f3a = fig.add_subplot(gs[2:,0]) 
    f3b = fig.add_subplot(gs[2:,1], sharex=f3a, sharey=f3a) 
    
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
        
        #cbaxes = inset_axes(f2, width="5%", height="30%", loc=4)           
        cbaxes = f2.inset_axes([0.91, 0.03, 0.05, 0.3])  
        cb = plt.colorbar(im, shrink=0.5, cax=cbaxes)
        cbaxes.yaxis.set_ticks_position('left')
        cbaxes.yaxis.set_label_position('left')
        #cb.ax.set_ylabel(scatCols[3],fontsize=12, labelpad=10, rotation=0)
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


# |--------------------------------|
# |------ Individual plots --------|
# |--------------------------------|

# Profile of mass vs distance for each case
def aMassEvol(nowKey, fold=''):
    myRes = res[nowKey]
    fig = plt.figure(constrained_layout=True, figsize=(5,5))
    gs = fig.add_gridspec(1, 1)
    f1 = fig.add_subplot(gs[0,0])   
    
    if nowKey in GCSkeys:
        if myRes.GCSprofHtsA is not None:
            for j in range(len(myRes.GCSprofHtsA)-1):
                i = j+1
                f1.plot(myRes.GCSprofHtsA[i].astype(float), myRes.GCSprofMassA[i].astype(float), '--', c=cGCSA, lw=1.5)
                f1.scatter(myRes.GCSprofHtsA[i].astype(float), myRes.GCSprofMassA[i].astype(float), c=cGCSA)
            f1.text(0.99, 0.2, 'GCS A', c=cGCSA, transform=f1.transAxes, ha='right')
        if myRes.GCSprofHtsB is not None:    
            for j in range(len(myRes.GCSprofHtsB)-1):
                i = j+1
                f1.plot(myRes.GCSprofHtsB[i].astype(float), myRes.GCSprofMassB[i].astype(float), '--', c=cGCSB, lw=1.5)
                f1.scatter(myRes.GCSprofHtsB[i].astype(float), myRes.GCSprofMassB[i].astype(float), c=cGCSB)
        f1.text(0.99, 0.15, 'GCS B', c=cGCSB, transform=f1.transAxes, ha='right')

    if myRes.CORSETheightsA is not None:
        f1.plot(myRes.CORSETheightsA, myRes.CORSETmassesA, 'r--', lw=2, zorder=10)
        f1.scatter(myRes.CORSETheightsA, myRes.CORSETmassesA, c='r', zorder=10)
        f1.text(0.99, 0.10, 'CORSET A', c=cCORA, transform=f1.transAxes, ha='right')
    if myRes.CORSETheightsB is not None:
        f1.plot(myRes.CORSETheightsB, myRes.CORSETmassesB, 'b--', lw=2, zorder=10)
        f1.scatter(myRes.CORSETheightsB, myRes.CORSETmassesB, c='b', zorder=10)
        f1.text(0.99, 0.05, 'CORSET B', c=cCORB, transform=f1.transAxes, ha='right')
    f1.text(0.01, 0.95, nowKey, c='k', transform=f1.transAxes, ha='left')
    
    f1.set_xlabel('Height (Rs)')
    f1.set_ylabel('Mass (10$^{15}$ g)')
    #plt.show()
    plt.savefig(fold+'massProfile_'+nowKey+'.png')
    plt.close()

def allMassEvols():
    for key in CORkeys:
        print (key)
        aMassEvol(key, fold='massProfFigs/')

# Timeline showing CORSET, GCS, and CDAW masses
def make3timeline():
    # Pull in CORSET Masses
    CORmassesA, CORmassesB = [], []
    CORdatesA, CORdatesB   = [], []
    for key in CORkeys:
        thisRes = res[key]
        if thisRes.CORSETmassMA:
            if thisRes.CORSETmassMA > 0:
                CORmassesA.append(thisRes.CORSETmassMA)
                CORdatesA.append(thisRes.CMEtimeDT)
        if thisRes.CORSETmassMB:
            if thisRes.CORSETmassMB > 0:
                CORmassesB.append(thisRes.CORSETmassMB)
                CORdatesB.append(thisRes.CMEtimeDT)
    # Pull in GCS Masses 
    GCSmassesA, GCSmassesB = [], []
    GCSmassesMA, GCSmassesMB = [], []
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
                    
    # get min/max date of CORSET
    minDT = np.min([np.min(CORdatesA), np.min(CORdatesB)])
    maxDT = np.max([np.max(CORdatesA), np.max(CORdatesB)])
                    
    # Pull in LASCO masses
    data = np.genfromtxt('fromAngelos/clean.list', dtype=None, skip_header=1)
    massCO = []
    widCO  = []
    dateCO = []
    for i in range(len(data)):      
        thisDate = data[i][0] + 'T' + data[i][1]
        #1996/01/26 12:40:19
        thisDT = datetime.datetime.strptime(thisDate, "%Y/%m/%dT%H:%M:%S" )
        if (thisDT >= minDT) & (thisDT <= maxDT):
            massCO.append(data[i][13])
            widCO.append(data[i][3])
            dateCO.append(thisDT)
    widCO = np.array(widCO)
    massCO = np.array(massCO)
    dateCO = np.array(dateCO)
    subMassCO = massCO[np.where(widCO > 80)]
    subDateCO = dateCO[np.where(widCO > 80)]
    
    m1A = np.array(CORmassesA)*1e15
    m1B = np.array(CORmassesB)*1e15
    m2A = np.array(GCSmassesA)*1e15
    m2B = np.array(GCSmassesB)*1e15    
    
    # Make the figure
    fig = plt.figure(constrained_layout=True, figsize=(10,7))
    gs = fig.add_gridspec(3, 1)
    f1 = fig.add_subplot(gs[0,0])   
    f2 = fig.add_subplot(gs[1,0], sharex=f1, sharey=f1)
    f3 = fig.add_subplot(gs[2,0], sharex=f1, sharey=f1)
    
    f1.scatter(CORdatesA, m1A, c=cCORA)
    f1.scatter(CORdatesB, m1B, c=cCORB)
    f2.scatter(GCSdatesA, m2A, c=cGCSA)
    f2.scatter(GCSdatesB, m2B, c=cGCSB)
    f3.scatter(dateCO, massCO, c='#696969')
    f3.scatter(subDateCO, subMassCO, c='darkgray')
    
    f1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    f1.set_yscale('log')
    f2.set_yscale('log')
    f3.set_yscale('log')
    
    f1.set_ylim([1e12,6e16])
    
    for f in [f1, f2, f3]:
        yrs0 = f1.get_ylim()
        yrs = [int(np.log10(yrs0[0])), int(np.log10(yrs0[1]))]
        diff = yrs[1] - yrs[0]
        toLine = [np.power(10, yrs[0]+1+x) for x in range(diff)]
        xrs = f.get_xlim()
        for x in toLine:
            f.plot(xrs, [x, x], 'k--',  alpha=0.5, zorder=0)
        f.set_xlim(xrs)
        f.set_ylim(yrs0)
        f.set_ylabel('Mass (g)')
    f1.set_ylim([1e12,6e16])
    
    f1.text(0.01, 0.9, 'CORSET A', horizontalalignment='left', transform = f1.transAxes, fontsize=12, color=cCORA)
    f1.text(0.12, 0.9, 'CORSET B', horizontalalignment='left', transform = f1.transAxes, fontsize=12, color=cCORB)
    f2.text(0.01, 0.9, 'GCS A', horizontalalignment='left', transform = f2.transAxes, fontsize=12, color=cGCSA)
    f2.text(0.08, 0.9, 'GCS B', horizontalalignment='left', transform = f2.transAxes, fontsize=12, color=cGCSB)
    f3.text(0.01, 0.9, 'CDAW', horizontalalignment='left', transform = f3.transAxes, fontsize=12, color='#696969')
    f3.text(0.08, 0.9, 'CDAW*', horizontalalignment='left', transform = f3.transAxes, fontsize=12, color='darkgray')
    
    plt.savefig('paperFigs/triTimeline.png')

# Timeline showing CORSET, GCS, and CDAW masses (with line plot of means )
def make4timeline():
    # Pull in CORSET Masses
    CORmassesA, CORmassesB = [], []
    CORdatesA, CORdatesB   = [], []
    for key in CORkeys:
        thisRes = res[key]
        if thisRes.CORSETmassMA:
            if thisRes.CORSETmassMA > 0:
                CORmassesA.append(thisRes.CORSETmassMA)
                CORdatesA.append(thisRes.CMEtimeDT)
        if thisRes.CORSETmassMB:
            if thisRes.CORSETmassMB > 0:
                CORmassesB.append(thisRes.CORSETmassMB)
                CORdatesB.append(thisRes.CMEtimeDT)
    # Pull in GCS Masses 
    GCSmassesA, GCSmassesB = [], []
    GCSmassesMA, GCSmassesMB = [], []
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
                    
    # get min/max date of CORSET
    minDT = np.min([np.min(CORdatesA), np.min(CORdatesB)])
    maxDT = np.max([np.max(CORdatesA), np.max(CORdatesB)])
                    
    # Pull in LASCO masses
    data = np.genfromtxt('fromAngelos/clean.list', dtype=None, skip_header=1)
    massCO = []
    widCO  = []
    dateCO = []
    for i in range(len(data)):      
        thisDate = data[i][0] + 'T' + data[i][1]
        #1996/01/26 12:40:19
        thisDT = datetime.datetime.strptime(thisDate, "%Y/%m/%dT%H:%M:%S" )
        if (thisDT >= minDT) & (thisDT <= maxDT):
            massCO.append(data[i][13])
            widCO.append(data[i][3])
            dateCO.append(thisDT)
    widCO = np.array(widCO)
    massCO = np.array(massCO)
    dateCO = np.array(dateCO)
    subMassCO = massCO[np.where(widCO > 80)]
    subDateCO = dateCO[np.where(widCO > 80)]
    
    m1A = np.array(CORmassesA)*1e15
    m1B = np.array(CORmassesB)*1e15
    m2A = np.array(GCSmassesA)*1e15
    m2B = np.array(GCSmassesB)*1e15    
    
    # Make the figure
    fig = plt.figure(constrained_layout=True, figsize=(10,10))
    gs = fig.add_gridspec(4, 1)
    f1 = fig.add_subplot(gs[0,0])   
    f2 = fig.add_subplot(gs[1,0], sharex=f1, sharey=f1)
    f3 = fig.add_subplot(gs[2,0], sharex=f1, sharey=f1)
    f4 = fig.add_subplot(gs[3,0], sharex=f1)
    
    f1.scatter(CORdatesA, m1A, c=cCORA)
    f1.scatter(CORdatesB, m1B, c=cCORB)
    f2.scatter(GCSdatesA, m2A, c=cGCSA)
    f2.scatter(GCSdatesB, m2B, c=cGCSB)
    f3.scatter(dateCO, massCO, c='#696969')
    f3.scatter(subDateCO, subMassCO, c='darkgray')
    
    f1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    f1.set_yscale('log')
    f2.set_yscale('log')
    f3.set_yscale('log')
    f4.set_yscale('log')
    
    
    for f in [f1, f2, f3, f4]:
        yrs0 = f1.get_ylim()
        yrs = [int(np.log10(yrs0[0])), int(np.log10(yrs0[1]))]
        diff = yrs[1] - yrs[0]
        toLine = [np.power(10, yrs[0]+1+x) for x in range(diff)]
        xrs = f.get_xlim()
        for x in toLine:
            f.plot(xrs, [x, x], 'k--',  alpha=0.5, zorder=0)
        f.set_xlim(xrs)
        f.set_ylim(yrs0)
        f.set_ylabel('Mass (g)')
    f1.set_ylim([1e12,6e16])
    
    f1.text(0.01, 0.9, 'CORSET A', horizontalalignment='left', transform = f1.transAxes, fontsize=12, color=cCORA)
    f1.text(0.12, 0.9, 'CORSET B', horizontalalignment='left', transform = f1.transAxes, fontsize=12, color=cCORB)
    f2.text(0.01, 0.9, 'GCS A', horizontalalignment='left', transform = f2.transAxes, fontsize=12, color=cGCSA)
    f2.text(0.08, 0.9, 'GCS B', horizontalalignment='left', transform = f2.transAxes, fontsize=12, color=cGCSB)
    f3.text(0.01, 0.9, 'CDAW', horizontalalignment='left', transform = f3.transAxes, fontsize=12, color='#696969')
    f3.text(0.08, 0.9, 'CDAW*', horizontalalignment='left', transform = f3.transAxes, fontsize=12, color='darkgray')
    
    # Add in the means
    meanMs = [[] for i in range(6)]
    stdMs  = [[] for i in range(6)]
    stdMus  = [[] for i in range(6)]
    stdMls  = [[] for i in range(6)]
    myTs   = [[] for i in range(6)]
    allDates = [CORdatesA, CORdatesB, GCSdatesA, GCSdatesB, dateCO, subDateCO]
    allMasses = [m1A, m1B, m2A, m2B, massCO, subMassCO]
    
    for i in range(6):
        allDates[i] = np.array(allDates[i])
        allMasses[i] = np.array(allMasses[i])
    
    minDT = datetime.datetime(2007,1,1)
    midDates = []
    for i in range(16):
        maxDT = minDT + datetime.timedelta(days=(365.25/2))
        #midDates.append(minDT + datetime.timedelta(days=(365.25/4)))
        for j in range(6):
            theseIdx = np.where((allDates[j] >= minDT) & (allDates[j] <= maxDT))[0]
            if len(theseIdx) > 0:
                meanMs[j].append(np.median(allMasses[j][theseIdx]))
                stdMs[j].append(np.std(allMasses[j][theseIdx]))
                stdMus[j].append(np.percentile(allMasses[j][theseIdx], 68) - np.median(allMasses[j][theseIdx]))
                stdMls[j].append( np.median(allMasses[j][theseIdx]) - np.percentile(allMasses[j][theseIdx], 32))
                myTs[j].append(minDT + datetime.timedelta(days=(365.25/4)))
            #else:
            #    meanMs[j].append(None)
            #    stdMs[j].append(0)
        minDT = maxDT
    
    cols = [cCORA, cCORB, cGCSA, cGCSB, '#696969', 'darkgray']
    for i in range(6)[::-1]:
        f4.plot(myTs[i], meanMs[i], c=cols[i], lw=3)   
        #print (np.array(stdMs[i])/1e15)
        f4.errorbar(myTs[i], meanMs[i], yerr=[stdMls[i], stdMus[i]], c=cols[i], capsize=5) 
        f4.scatter(myTs[i], meanMs[i], c=cols[i])    
    f4.set_ylim([3e13,3e16])
    f1.set_ylim([1e12,1e17])
    
    
    
    plt.savefig('paperFigs/triTimelinePlus.png')

# Histograms of even more different masses
def megaHisto():
    # Pull in CORSET Masses
    CORmassesA, CORmassesB = [], []
    CORdatesA, CORdatesB   = [], []
    for key in CORkeys:
        thisRes = res[key]
        if thisRes.CORSETmassMA:
            if thisRes.CORSETmassMA > 0:
                CORmassesA.append(thisRes.CORSETmassMA)
                CORdatesA.append(thisRes.CMEtimeDT)
        if thisRes.CORSETmassMB:
            if thisRes.CORSETmassMB > 0:
                CORmassesB.append(thisRes.CORSETmassMB)
                CORdatesB.append(thisRes.CMEtimeDT)
    # Pull in GCS Masses 
    GCSmassesA, GCSmassesB = [], []
    GCSmassesMA, GCSmassesMB = [], []
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
                    
    LBmassesA, LBmassesB = [], []
    for key in CORkeys:
        thisRes = res[key]
        if thisRes.LBdeprojMA:
            if thisRes.LBdeprojMA > 0:
                LBmassesA.append(thisRes.LBdeprojMA)
        if thisRes.LBdeprojMB:
            if thisRes.LBdeprojMB > 0:
                LBmassesB.append(thisRes.LBdeprojMB)
        
        


    fig = plt.figure(constrained_layout=True, figsize=(10,10))
    gs = fig.add_gridspec(3, 4, hspace=0, wspace=0)
    f1a = fig.add_subplot(gs[0,0]) 
    f1b = fig.add_subplot(gs[0,1], sharex=f1a, sharey=f1a) 
    f1c = fig.add_subplot(gs[0,2], sharex=f1a, sharey=f1a) 
    f1d = fig.add_subplot(gs[0,3], sharex=f1a, sharey=f1a) 
    f2a = fig.add_subplot(gs[1,0], sharex=f1a, sharey=f1a) 
    f2b = fig.add_subplot(gs[1,1], sharex=f1a, sharey=f1a) 
    f2c = fig.add_subplot(gs[1,2], sharex=f1a, sharey=f1a) 
    f2d = fig.add_subplot(gs[1,3], sharex=f1a, sharey=f1a) 
    f3a = fig.add_subplot(gs[2,1], sharex=f1a, sharey=f1a) 
    f3b = fig.add_subplot(gs[2,2], sharex=f1a, sharey=f1a) 
    f3c = fig.add_subplot(gs[2,3], sharex=f1a, sharey=f1a) 
    f3 = fig.add_subplot(gs[2,0], sharex=f1a, sharey=f1a) 
    
    allfs = [f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c]
    rts = [f1c, f2c, f3c]
    lts = [f1a, f2a]
    notBots = [f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d]
    
    massCORA = np.log10(np.array(CORmassesA)*1e15) 
    massCORB = np.log10(np.array(CORmassesB)*1e15)
    massGCSA = np.log10(np.array(GCSmassesMA)*1e15)
    massGCSB = np.log10(np.array(GCSmassesMB)*1e15)
    LBmassesA = np.log10(np.array(LBmassesA)*1e15)
    LBmassesB = np.log10(np.array(LBmassesB)*1e15)
    
    massDP1, massDP2A, massDP2B = [], [], []
    for key in CORkeys:
            thisRes = res[key]
            # Check if had matched A/B and could deproj
            if thisRes.deprojTimes[0]:
                maxIdx = np.where(thisRes.deprojMasses == np.max(thisRes.deprojMasses))[0]
                massDP1.append(thisRes.deprojMasses[maxIdx[0]])
                
                maxIdxA = np.where(thisRes.CdeprojMassesA == np.max(thisRes.CdeprojMassesA))[0]
                maxIdxB = np.where(thisRes.CdeprojMassesB == np.max(thisRes.CdeprojMassesB))[0]
                massDP2A.append(thisRes.CdeprojMassesA[maxIdxA[0]])
                massDP2B.append(thisRes.CdeprojMassesB[maxIdxB[0]])
    massDP1 = np.array(massDP1)
    massDP2A = np.array(massDP2A)
    massDP2B = np.array(massDP2B)
    massDP1 = np.log10(np.array(massDP1[np.where(massDP1 > 0)])*1e15) 
    massDP2A = np.log10(np.array(massDP2A[np.where(massDP2A > 0)])*1e15) 
    massDP2B = np.log10(np.array(massDP2B[np.where(massDP2B > 0)])*1e15)
    
    # read in lasco
    data = np.genfromtxt('fromAngelos/clean.list', dtype=None, skip_header=1)
    massCO = []
    widCO  = []
    for i in range(len(data)):
        massCO.append(data[i][13])
        widCO.append(data[i][3])
    widCO = np.array(widCO)
    massCO = np.log10(np.array(massCO))
    subMassCO = massCO[np.where(widCO > 80.)]
    
    bins1 = np.linspace(13,17,21)
    f1a.set_xlim([13,17])
    f1a.set_ylim([0,1.3])
    f1a.hist(massCO, bins=bins1, density=True, color='#696969', ec='k')
    f1b.hist(subMassCO, bins=bins1, density=True, color='darkgray', ec='k')
    f1c.hist(massCORA, bins=bins1, density=True, color=cCORA, ec='k')
    f1d.hist(massCORB, bins=bins1, density=True, color=cCORB, ec='k')
    f2a.hist(massGCSA, bins=bins1, density=True, color=cGCSA, ec='k')
    f2b.hist(massGCSB, bins=bins1, density=True, color=cGCSB, ec='k')
    f2c.hist(LBmassesA, bins=bins1, density=True, color=cLBA, ec='k')
    f2d.hist(LBmassesB, bins=bins1, density=True, color=cLBB, ec='k')

    f3a.hist(massDP1, bins=bins1, density=True, color='#6f10be', ec='k')
    f3b.hist(massDP2A, bins=bins1, density=True, color='#AA4499', ec='k')
    f3c.hist(massDP2B, bins=bins1, density=True, color='#367cc3', ec='k')
    
    f1a.text(0.05, 0.9, 'CDAW', transform=f1a.transAxes, fontsize=12, ha='left', va='center', c='#696969')
    f1a.text(0.05, 0.83, 'N: '+str(len(massCO)), transform=f1a.transAxes, fontsize=12, ha='left', va='center', c='#696969')
    f1b.text(0.05, 0.9, 'CDAW$^*$', transform=f1b.transAxes, fontsize=12, ha='left', va='center', c='darkgray')
    f1b.text(0.05, 0.83, 'N: '+str(len(subMassCO)), transform=f1b.transAxes, fontsize=12, ha='left', va='center', c='darkgray')
    f1c.text(0.05, 0.9, 'CORSET A', transform=f1c.transAxes, fontsize=12, ha='left', va='center', c=cCORA)
    f1c.text(0.05, 0.83, 'N: '+str(len(massCORA)), transform=f1c.transAxes, fontsize=12, ha='left', va='center', c=cCORA)
    f1d.text(0.05, 0.9, 'CORSET B', transform=f1d.transAxes, fontsize=12, ha='left', va='center', c=cCORB)
    f1d.text(0.05, 0.83, 'N: '+str(len(massCORB)), transform=f1d.transAxes, fontsize=12, ha='left', va='center', c=cCORB)
    
    f2a.text(0.05, 0.9, 'GCS A', transform=f2a.transAxes, fontsize=12, ha='left', va='center', c=cGCSA)
    f2a.text(0.05, 0.83, 'N: '+str(len(massGCSA)), transform=f2a.transAxes, fontsize=12, ha='left', va='center', c=cGCSA)
    f2b.text(0.05, 0.9, 'GCS B', transform=f2b.transAxes, fontsize=12, ha='left', va='center', c=cGCSB)
    f2b.text(0.05, 0.83, 'N: '+str(len(massGCSB)), transform=f2b.transAxes, fontsize=12, ha='left', va='center', c=cGCSB)
    f2c.text(0.05, 0.9, 'DP1A', transform=f2c.transAxes, fontsize=12, ha='left', va='center', c=cLBA)
    f2c.text(0.05, 0.83, 'N: '+str(len(LBmassesA)), transform=f2c.transAxes, fontsize=12, ha='left', va='center', c=cLBA)
    f2d.text(0.05, 0.9, 'DP1B', transform=f2d.transAxes, fontsize=12, ha='left', va='center', c=cLBB)
    f2d.text(0.05, 0.83, 'N: '+str(len(LBmassesB)), transform=f2d.transAxes, fontsize=12, ha='left', va='center', c=cLBB)
    
    f3a.text(0.05, 0.9, 'DP2', transform=f3a.transAxes, fontsize=12, ha='left', va='center', c='#6f10be')
    f3a.text(0.05, 0.83, 'N: '+str(len(massDP1)), transform=f3a.transAxes, fontsize=12, ha='left', va='center', c='#6f10be')
    f3b.text(0.05, 0.9, 'DP3A', transform=f3b.transAxes, fontsize=12, ha='left', va='center', c='#AA4499')
    f3b.text(0.05, 0.83, 'N: '+str(len(massDP2A)), transform=f3b.transAxes, fontsize=12, ha='left', va='center', c='#AA4499')
    f3c.text(0.05, 0.9, 'DP3B', transform=f3c.transAxes, fontsize=12, ha='left', va='center', c='#367cc3')
    f3c.text(0.05, 0.83, 'N: '+str(len(massDP2B)), transform=f3c.transAxes, fontsize=12, ha='left', va='center', c='#367cc3')
                   
    for f in allfs:
        xl = f.get_xlim()
        yl = f.get_ylim()
        for y in [0.25, 0.5, 0.75, 1, 1.25]:
            f.plot(xl, [y, y], 'k--', alpha=0.25, zorder=0)
        for x in [14, 15, 16]:
            f.plot([x,x], yl, 'k--', alpha=0.25, zorder=0)
        if f not in lts:
            f.tick_params(labelleft=False)   
        if f in notBots:
            f.tick_params(labelbottom=False)   
        if f in lts:
            f.set_ylabel('Prob. Density')
        if f in [f3a, f3b, f3c]:
            f.set_xlabel('Log Mass (g)')
        
        f3.yaxis.set_label_position("right")
        f3.yaxis.tick_right()
        f3.tick_params(axis="y",direction="in", pad=-22)
        f3.set_ylabel('Prob. Density', labelpad=-40)
        f3.spines['top'].set_visible(False)
        f3.spines['right'].set_visible(False)
        f3.spines['left'].set_visible(False)
        f3.spines['bottom'].set_visible(False)
        f3.get_xaxis().set_visible(False)
        f3.tick_params(axis='y', length=0)
    
        f.set_xlim(xl)
        f.set_ylim(yl)
    plt.savefig(myFold+'MegaMassHisto.png')

# Histograms of all the different masses
def megaHistoOld():
    # Pull in CORSET Masses
    CORmassesA, CORmassesB = [], []
    CORdatesA, CORdatesB   = [], []
    for key in CORkeys:
        thisRes = res[key]
        if thisRes.CORSETmassMA:
            if thisRes.CORSETmassMA > 0:
                CORmassesA.append(thisRes.CORSETmassMA)
                CORdatesA.append(thisRes.CMEtimeDT)
        if thisRes.CORSETmassMB:
            if thisRes.CORSETmassMB > 0:
                CORmassesB.append(thisRes.CORSETmassMB)
                CORdatesB.append(thisRes.CMEtimeDT)
    # Pull in GCS Masses 
    GCSmassesA, GCSmassesB = [], []
    GCSmassesMA, GCSmassesMB = [], []
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


    fig = plt.figure(constrained_layout=True, figsize=(10,10))
    gs = fig.add_gridspec(3, 3, hspace=0, wspace=0)
    f1a = fig.add_subplot(gs[0,0]) 
    f1b = fig.add_subplot(gs[0,1], sharex=f1a, sharey=f1a) 
    f1c = fig.add_subplot(gs[0,2], sharex=f1a, sharey=f1a) 
    f2a = fig.add_subplot(gs[1,0], sharex=f1a, sharey=f1a) 
    f2b = fig.add_subplot(gs[1,1], sharex=f1a, sharey=f1a) 
    f2c = fig.add_subplot(gs[1,2], sharex=f1a, sharey=f1a) 
    f3a = fig.add_subplot(gs[2,0], sharex=f1a, sharey=f1a) 
    f3b = fig.add_subplot(gs[2,1], sharex=f1a, sharey=f1a) 
    f3c = fig.add_subplot(gs[2,2], sharex=f1a, sharey=f1a) 
    
    allfs = [f1a, f1b, f1c, f2a, f2b, f2c, f3a, f3b, f3c]
    rts = [f1c, f2c, f3c]
    lts = [f1a, f2a, f3a]
    notBots = [f1a, f1b, f1c, f2a, f2b, f2c]
    
    massCORA = np.log10(np.array(CORmassesA)*1e15) 
    massCORB = np.log10(np.array(CORmassesB)*1e15)
    massGCSA = np.log10(np.array(GCSmassesMA)*1e15)
    massGCSB = np.log10(np.array(GCSmassesMB)*1e15)
    
    massDP1, massDP2A, massDP2B = [], [], []
    for key in CORkeys:
            thisRes = res[key]
            # Check if had matched A/B and could deproj
            if thisRes.deprojTimes[0]:
                maxIdx = np.where(thisRes.deprojMasses == np.max(thisRes.deprojMasses))[0]
                massDP1.append(thisRes.deprojMasses[maxIdx[0]])
                
                maxIdxA = np.where(thisRes.CdeprojMassesA == np.max(thisRes.CdeprojMassesA))[0]
                maxIdxB = np.where(thisRes.CdeprojMassesB == np.max(thisRes.CdeprojMassesB))[0]
                massDP2A.append(thisRes.CdeprojMassesA[maxIdxA[0]])
                massDP2B.append(thisRes.CdeprojMassesB[maxIdxB[0]])
    massDP1 = np.array(massDP1)
    massDP2A = np.array(massDP2A)
    massDP2B = np.array(massDP2B)
    massDP1 = np.log10(np.array(massDP1[np.where(massDP1 > 0)])*1e15) 
    massDP2A = np.log10(np.array(massDP2A[np.where(massDP2A > 0)])*1e15) 
    massDP2B = np.log10(np.array(massDP2B[np.where(massDP2B > 0)])*1e15)
    
    # read in lasco
    data = np.genfromtxt('fromAngelos/clean.list', dtype=None, skip_header=1)
    massCO = []
    widCO  = []
    for i in range(len(data)):
        massCO.append(data[i][13])
        widCO.append(data[i][3])
    widCO = np.array(widCO)
    massCO = np.log10(np.array(massCO))
    subMassCO = massCO[np.where(widCO > 80.)]
    
    bins1 = np.linspace(13,17,21)
    f1a.set_xlim([13,17])
    f1a.hist(massCO, bins=bins1, density=True, color='#696969', ec='k')
    f1b.hist(massCORA, bins=bins1, density=True, color=cCORA, ec='k')
    f1c.hist(massCORB, bins=bins1, density=True, color=cCORB, ec='k')
    f2a.hist(subMassCO, bins=bins1, density=True, color='darkgray', ec='k')
    f2b.hist(massGCSA, bins=bins1, density=True, color=cGCSA, ec='k')
    f2c.hist(massGCSB, bins=bins1, density=True, color=cGCSB, ec='k')
    f3a.hist(massDP1, bins=bins1, density=True, color='#6f10be', ec='k')
    f3b.hist(massDP2A, bins=bins1, density=True, color='#AA4499', ec='k')
    f3c.hist(massDP2B, bins=bins1, density=True, color='#367cc3', ec='k')
    
    f1a.text(0.05, 0.9, 'CDAW', transform=f1a.transAxes, fontsize=12, ha='left', va='center', c='#696969')
    f1a.text(0.05, 0.81, 'N: '+str(len(massCO)), transform=f1a.transAxes, fontsize=12, ha='left', va='center', c='#696969')
    f1b.text(0.05, 0.9, 'CORSET A', transform=f1b.transAxes, fontsize=12, ha='left', va='center', c=cCORA)
    f1b.text(0.05, 0.81, 'N: '+str(len(massCORA)), transform=f1b.transAxes, fontsize=12, ha='left', va='center', c=cCORA)
    f1c.text(0.05, 0.9, 'CORSET B', transform=f1c.transAxes, fontsize=12, ha='left', va='center', c=cCORB)
    f1c.text(0.05, 0.81, 'N: '+str(len(massCORB)), transform=f1c.transAxes, fontsize=12, ha='left', va='center', c=cCORB)
    
    f2a.text(0.05, 0.9, 'CDAW >80$^{\\circ}$', transform=f2a.transAxes, fontsize=12, ha='left', va='center', c='darkgray')
    f2a.text(0.05, 0.81, 'N: '+str(len(subMassCO)), transform=f2a.transAxes, fontsize=12, ha='left', va='center', c='darkgray')
    f2b.text(0.05, 0.9, 'GCS A', transform=f2b.transAxes, fontsize=12, ha='left', va='center', c=cGCSA)
    f2b.text(0.05, 0.81, 'N: '+str(len(massGCSA)), transform=f2b.transAxes, fontsize=12, ha='left', va='center', c=cGCSA)
    f2c.text(0.05, 0.9, 'GCS B', transform=f2c.transAxes, fontsize=12, ha='left', va='center', c=cGCSB)
    f2c.text(0.05, 0.81, 'N: '+str(len(massGCSB)), transform=f2c.transAxes, fontsize=12, ha='left', va='center', c=cGCSB)
    
    f3a.text(0.05, 0.9, 'Deproj1', transform=f3a.transAxes, fontsize=12, ha='left', va='center', c='#6f10be')
    f3a.text(0.05, 0.81, 'N: '+str(len(massDP1)), transform=f3a.transAxes, fontsize=12, ha='left', va='center', c='#6f10be')
    f3b.text(0.05, 0.9, 'Deproj2 A', transform=f3b.transAxes, fontsize=12, ha='left', va='center', c='#AA4499')
    f3b.text(0.05, 0.81, 'N: '+str(len(massDP2A)), transform=f3b.transAxes, fontsize=12, ha='left', va='center', c='#AA4499')
    f3c.text(0.05, 0.9, 'Deproj2 B', transform=f3c.transAxes, fontsize=12, ha='left', va='center', c='#367cc3')
    f3c.text(0.05, 0.81, 'N: '+str(len(massDP2B)), transform=f3c.transAxes, fontsize=12, ha='left', va='center', c='#367cc3')
                   
    for f in allfs:
        xl = f.get_xlim()
        yl = f.get_ylim()
        for y in [0.5, 1]:
            f.plot(xl, [y, y], 'k--', alpha=0.25, zorder=0)
        for x in [14, 15, 16]:
            f.plot([x,x], yl, 'k--', alpha=0.25, zorder=0)
        if f in rts:
            f.tick_params(labelleft=False)   
        if f in notBots:
            f.tick_params(labelbottom=False)   
        if f in lts:
            f.set_ylabel('Prob. Density')
        if f in [f3a, f3b, f3c]:
            f.set_xlabel('Log Mass (g)')
    
        f.set_xlim(xl)
        f.set_ylim(yl)
    plt.savefig(myFold+'MegaMassHisto.png')

    

# Comparison of STA to STB for COR & GCS
def A2Bcomp():
    pCORmassesA, pCORmassesB = [], []
    pCORvelsA, pCORvelsB = [], []
    pCORdates = []

    for key in CORkeys:
        thisRes = res[key]
        if (type(thisRes.CORSETmassMA) != type(None)) & (type(thisRes.CORSETmassMB) != type(None)):
            if (thisRes.CORSETmassMA > 0) & (thisRes.CORSETmassMB > 0):
                pCORmassesA.append(thisRes.CORSETmassMA)
                pCORvelsA.append(thisRes.CORSETvelMA)
                pCORdates.append(thisRes.CMEtimeDT)
                pCORmassesB.append(thisRes.CORSETmassMB)
                pCORvelsB.append(thisRes.CORSETvelMB)
                
    pGCSmassesA, pGCSmassesB = [], []
    pGCSvels, pGCSidx = [], []
    pGCSdates = []
    counter = 0
    for key in GCSkeys:
        thisRes = res[key]
        myKeys = np.array([key for key in thisRes.GCSvals.keys()])
        if (isinstance(thisRes.GCSmassesA,list)) & (isinstance(thisRes.GCSmassesB, list)):
            for i in range(len(thisRes.GCSmassesA)):
                mA, mB = thisRes.GCSmassesA[i], thisRes.GCSmassesB[i]
                myVel = thisRes.GCSvals[myKeys[i]][-1]
                if (mA > 0)  & (mB > 0):
                    pGCSmassesA.append(mA)
                    pGCSmassesB.append(mB)
                    pGCSdates.append(thisRes.CMEtimeDT)
                    if myVel != 'None':
                        pGCSvels.append(float(myVel))
                        pGCSidx.append(counter)
                    counter += 1
    
    p1 = np.array(pCORmassesA)*1e15
    p2 = np.array(pCORmassesB)*1e15
    p3 = np.array(pGCSmassesA)*1e15 
    p4 = np.array(pGCSmassesB)*1e15 
    print (r2_score(np.log10(p2),np.log10(p1)))
    print (r2_score(np.log10(p4),np.log10(p3)))
    
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
    
    scatter2histo2(p1, p2, p3, p4, axLabs, picName=myFold+'A2Bcomp.png', pLabs=pLabs, axLim=[[13.5,16.7],[13.5,16.7]], colors=colors, scatCols=[c1,c2,[400,1000],'v\n(km/s)'], bins=bins )
    
# Comparision of COR to GCS for A & B    
def C2Gcomp():
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
    
    p1 = np.array(cgCORmassesA)*1e15
    p2 = np.array(cgGCSmassesA)*1e15
    p3 = np.array(cgCORmassesB)*1e15 
    p4 = np.array(cgGCSmassesB)*1e15 
    c1 = [cgvelsA, None]
    c2 = [cgvelsB, None]
    
    axLabs = ['CORSET log Mass (g)', 'GCS log Mass (g)', '(M$_{COR}$ - M$_{GCS}$) / <M>']
    pLabs  = ['STEREO A', 'STEREO B']
    colors = [cGCSA, cGCSB]
    b1 = np.arange(-1.5,1.5, 0.2, dtype=float)
    b1 = np.append(b1, 2*b1[-1] - b1[-2])
    b2 = np.arange(-1.5,1.5, 0.2, dtype=float)
    b2 = np.append(b2, 2*b2[-1] - b2[-2])
    bins = [b1, b2]
    
    scatter2histo2(p1, p2, p3, p4, axLabs, picName=myFold+'COR2GCScomp.png', pLabs=pLabs, axLim=[[14.5,16.5],[14.5,16.5]], colors=colors, scatCols=[c1,c2,[400,1000],'v\n(km/s)'], bins=bins )
    
# Scatter plot based on CME-PoS and A-B separations
def sepPlot():
    satSep = []
    posSepA2, posSepB2 = [], []
    hRatioA2, hRatioB2 = [], []
    mRatioA2, mRatioB2 = [], []
    for key in CORkeys:
        thisRes = res[key]
        if thisRes.CdeprojTimes[0]:
            # deproj 2
            maxIdxA = np.where(thisRes.CdeprojMassesA == np.max(thisRes.CdeprojMassesA))[0]
            maxIdxA = np.max(maxIdxA)
            maxIdxB = np.where(thisRes.CdeprojMassesB == np.max(thisRes.CdeprojMassesB))[0]
            maxIdxB = np.max(maxIdxB)
            maxIdx = np.max([maxIdxA, maxIdxB])
            sep = np.abs(thisRes.depSatLonsB[maxIdx] - thisRes.depSatLonsA[maxIdx])
            if sep > 180:
                sep = 360 - sep
            
            phAs = thisRes.CprojHeightsA.astype(float)
            phBs = thisRes.CprojHeightsB.astype(float)
            dphAs = thisRes.CdeprojHeightsA.astype(float)
            dphBs = thisRes.CdeprojHeightsB.astype(float)
            pMAs  = thisRes.CprojMassesA
            pMBs  = thisRes.CprojMassesB
            dpMAs  = thisRes.CdeprojMassesA
            dpMBs  = thisRes.CdeprojMassesB
            
            # check for zeros
            if (phAs[maxIdx] == 0) or (phBs[maxIdx] == 0):
                newIdxs = np.where((phAs > 0) & (phAs < 35) & (phBs >0) & (phBs < 35))[0]
                if len(newIdxs) > 0:
                    maxIdx = np.max(newIdxs)
                else:
                    maxIdx = -1
            
            if maxIdx != -1:
                satSep.append(sep)
                posSepA2.append(thisRes.CdeprojSepA[maxIdx])
                posSepB2.append(thisRes.CdeprojSepB[maxIdx])
                hRatioA2.append(np.abs(phAs[maxIdx] / dphAs[maxIdx]))
                hRatioB2.append(np.abs(phBs[maxIdx] / dphBs[maxIdx]))
                mRatioA2.append(pMAs[maxIdx] / dpMAs[maxIdx])
                mRatioB2.append(pMBs[maxIdx] / dpMBs[maxIdx])
            
    
    fig = plt.figure(constrained_layout=True, figsize=(7,6))
    gs = fig.add_gridspec(1, 1)
    f1 = fig.add_subplot(gs[0,0])   
    normCmap = mcolors.Normalize(vmin=0, vmax=1)
    
    f1.scatter(satSep, posSepA2, c=mRatioA2, cmap='plasma', norm=normCmap)
    im = f1.scatter(satSep, posSepB2, c=mRatioB2, cmap='plasma', norm=normCmap)
    
    cb = plt.colorbar(im, shrink=0.75)
    cb.ax.set_title('Mass\nRatio',fontsize=12)
    cb.ax.tick_params(labelsize=12)
    
    f1.set_xlabel('A/B Separation ($^{\\circ}$)')
    f1.set_ylabel('CME/PoS Separation ($^{\\circ}$)')
    f1.set_xlim([0,180])
    f1.set_ylim([0,90])
    
    plt.savefig('paperFigs/sepFig.png')

# Geometric visualization of diff sat views
def coveragePlot(critSep):
    dtor = np.pi/180.
    # Set up figure
    fig = plt.figure(constrained_layout=True, figsize=(12,10))
    gs = fig.add_gridspec(7, 6)  
    f1a = fig.add_subplot(gs[0:3,0:3])
    f1b = fig.add_subplot(gs[0:3,3:])
    f2a = fig.add_subplot(gs[3:5,0:2], projection='polar')
    f2b = fig.add_subplot(gs[3:5,2:4], projection='polar')
    f2c = fig.add_subplot(gs[3:5,4:], projection='polar')
    f3a = fig.add_subplot(gs[5:,0:2], projection='polar')
    f3b = fig.add_subplot(gs[5:,2:4], projection='polar')
    f3c = fig.add_subplot(gs[5:,4:], projection='polar')
    
    # Do the top plots
    nPts = 3600
    sepLim = critSep
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
    
    ccols = ['w', '#FF9999', '#9999FF', '#995BC1']
    
    f1a.contourf(satLons, np.linspace(-179,180,nPts), allVals, levels=[-1, 0, 1,2,3], colors=ccols)
    f1a.contour(satLons, np.linspace(-179,180,nPts), allVals, levels=[0, 1,2,3], colors='k')
    f1a.plot([0,90],[0,90], '--', c='r', lw=1)
    f1a.plot([0,90],[0,-90], '--', c='b', lw=1)
    f1a.plot([0,90], [0,0], '--', c='k', lw=2)
    #f1a.plot([60,60], [-180,180], '--', c='#696969', lw=2)
    #f1b.plot([60,60], [0,100], '--', c='#696969', lw=2)
    f1b.contourf(satLons, np.linspace(0,100,nPts), allOrd, levels=[-1, 0, 1,2,3], colors=ccols)
    f1b.contour(satLons, np.linspace(0,100,nPts), allOrd, levels=[-1, 0, 1,2,3], colors='k')

    f1a.set_ylim([-180,180])
    f1a.set_xlabel('Sat Lons ($\\pm ^{\\circ}$)')
    f1a.set_ylabel('CME Stony Lon ($^{\\circ}$)')
    f1b.set_xlabel('Sat Lons ($\\pm ^{\\circ}$)')
    f1b.set_ylabel('Percent Coverage')
    
    f1b.text(0.5, 0.04, 'A+B', horizontalalignment='center', transform = f1b.transAxes, fontsize=14, color='k')
    f1b.text(0.5, 0.4, 'B', horizontalalignment='center', transform = f1b.transAxes, fontsize=14, color='k')
    f1b.text(0.5, 0.7, 'A', horizontalalignment='center', transform = f1b.transAxes, fontsize=14, color='k')
    
    # Do the polar plots
    pfs = [f2a, f2b, f2c, f3a, f3b, f3c]
    sLons = [40, 60, 80, 40, 60, 80]
    for i in range(6):
        ax = pfs[i]
        satA = sLons[i] * dtor
        satB = -sLons[i] * dtor
        sep  = critSep * dtor
        
        npts = 1000
        thetas = np.linspace(0,2*np.pi, npts)
        ax.plot(thetas, np.ones(npts), 'k', lw=3)
        #ax.plot([0,0], [0,1], 'k--', lw=2 )
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
        
        if i > 2:
            angsL1a = np.linspace(np.pi/2-sep,np.pi/2+sep, npts)
            ax.fill_between(angsL1a, 0, 1, color='k', alpha=0.3)
            angsL1b = np.linspace(-np.pi/2-sep,-np.pi/2+sep, npts)
            ax.fill_between(angsL1b, 0, 1, color='k', alpha=0.3)
            ax.plot([-np.pi/2, np.pi/2], [1,1], 'k--', lw=2, zorder=9)
            
        
        ax.set_aspect('equal')
        ax.set_ylim([0,1.1])
        ax.axis('off')
        ax.text(0.5, 1, '$\\pm$'+ str(sLons[i])+'$^{\\circ}$', horizontalalignment='center', transform = ax.transAxes, fontsize=14)
    
    plt.savefig('paperFigs/coverage60.png')
    #plt.show()      

# Triangle plots comparing reconstructed values
def MLLscattersOld():
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
            compMasses['CORB'].append(np.max(thisRes.CORSETmassMB))
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
            counter += 1
    
    for key in compLats:
        compLats[key] = np.array(compLats[key])
        compLons[key] = np.array(compLons[key])
        compMasses[key] = np.array(compMasses[key])
    compTimes = np.array(compTimes)
    compHasGCS = np.array(compHasGCS)
        
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

    dlon = 30
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
        ff.plot([ll, ul], [ll+dlon, ul+dlon], 'k--', zorder=0, alpha=0.4)
        ff.plot([ll, ul], [ll-dlon, ul-dlon], 'k--', zorder=0, alpha=0.4)
        ff.scatter(x,y, marker=mkr, c=c1sat)
        i +=1
    
    ff1b.text(1.5, 0.5, 'Deproj1', transform=ff1b.transAxes, fontsize=12, ha='center', va='center')
    ff2a.text(1.5, 0.5, 'Deproj2', transform=ff2a.transAxes, fontsize=12, ha='center', va='center')
    ff2a.text(0.5, -0.2, 'GCS', transform=ff2a.transAxes, fontsize=12, ha='center', va='center')
    
    
    ll, ul = -70, 70
    f1a.set_xlim(ll,ul)
    f1a.set_ylim(ll,ul)
    f1a.set_aspect('equal')
    dlat = 10
    for f in allfs:
        f.plot([ll, ul], [ll, ul], 'k--', zorder=0, alpha=0.75)
        f.plot([ll, ul], [ll+dlat, ul+dlat], 'k--', zorder=0, alpha=0.4)
        f.plot([ll, ul], [ll-dlat, ul-dlat], 'k--', zorder=0, alpha=0.4)
    
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

# Triangle plots comparing reconstructed values (incl LB recons)
def MLLscatters():
    compTimes, compLons, compLats, compMasses, compHasGCS, compHasLB, compHasLBGCS = [], {}, {}, {}, [], [], []
    counter = 0
    setupNames = ['COR', 'GCS',  'DP2', 'DP3', 'DP1' ]
    holders = [compLons, compLats, compMasses]
    for holder in holders:
        for name in setupNames:
            holder[name+'A'] = []
            holder[name+'B'] = []
    
    compSeps = {'DP3A':[], 'DP3B':[]}
    for key in CORkeys:
        thisRes = res[key]
        # Check if had matched A/B and could deproj
        if thisRes.deprojTimes[0]:
            # add to the time array
            compTimes.append(thisRes.CMEtimeDT)
            
            # Get basic corset values (mass + lat from CPA)
            #compMasses['CORA'].append(thisRes.CORSETmassMA)
            compMasses['CORA'].append(np.max(np.abs(thisRes.CORSETmassesA)))
            compMasses['CORB'].append(np.max(np.abs(thisRes.CORSETmassMB)))
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
            compMasses['DP2A'].append(np.abs(thisRes.deprojMasses[maxIdx[0]]))
            compMasses['DP2B'].append(np.abs(thisRes.deprojMasses[maxIdx[0]]))
            # Lon is exact same as output of DP
            maxLon = thisRes.deprojLons[maxIdx[0]]
            shiftLonDiff = thisRes.deprojLons - maxLon
            shiftLonDiff[np.where(shiftLonDiff > 180)] -= 360
            shiftLonDiff[np.where(shiftLonDiff < -180)] += 360
            newMedLon = np.median(shiftLonDiff) + maxLon
            compLons['DP2A'].append(newMedLon)
            compLats['DP2A'].append(np.median(thisRes.deprojLatA))
            compLats['DP2B'].append(np.median(thisRes.deprojLatB))
            
            
            
            # Get DP2 (mass+lat+height) values
            maxIdxA = np.where(thisRes.CdeprojMassesA == np.max(thisRes.CdeprojMassesA))[0]
            maxIdxB = np.where(thisRes.CdeprojMassesB == np.max(thisRes.CdeprojMassesB))[0]
            if maxIdxA[0] != maxIdxB[0]:
                maxIdx = np.max([maxIdxA[0], maxIdxB[0]])
            else:
                maxIdx = maxIdxA[0]
            compMasses['DP3A'].append(np.abs(thisRes.CdeprojMassesA[maxIdxA[0]]))
            compMasses['DP3B'].append(np.abs(thisRes.CdeprojMassesB[maxIdxB[0]]))
            maxLon = thisRes.CdeprojLons[maxIdx]
            shiftLonDiff = thisRes.CdeprojLons - maxLon
            shiftLonDiff[np.where(shiftLonDiff > 180)] -= 360
            shiftLonDiff[np.where(shiftLonDiff < -180)] += 360
            newMedLon = np.median(shiftLonDiff) + maxLon
            compLons['DP3A'].append(newMedLon)
            
            compLats['DP3A'].append(np.median(thisRes.CdeprojLatA))
            compLats['DP3B'].append(np.median(thisRes.CdeprojLatB))
            
            compSeps['DP3A'].append(thisRes.CdeprojSepA[maxIdx])
            compSeps['DP3B'].append(thisRes.CdeprojSepB[maxIdx])
            
            
            # Get Laura values
            if thisRes.LBdeprojMA:
                compHasLB.append(counter)
                compLats['DP1A'].append(thisRes.LBdeprojLat)
                compLons['DP1A'].append(thisRes.LBdeprojLon)
                compMasses['DP1A'].append(np.abs(thisRes.LBdeprojMA))
                compMasses['DP1B'].append(np.abs(thisRes.LBdeprojMB))
            else:
                compLats['DP1A'].append(None)
                compLons['DP1A'].append(None)
                compMasses['DP1A'].append(None)
                compMasses['DP1B'].append(None)
                
                
            
            # check if has GCS values
            if key in GCSkeys:
                if (isinstance(thisRes.GCSmassesA,list)) & (isinstance(thisRes.GCSmassesB, list)):
                    compHasGCS.append(counter)
                    if counter in compHasLB:
                        compHasLBGCS.append(counter)
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
            else:
                compMasses['GCSA'].append(None)
                compMasses['GCSB'].append(None)
                compLons['GCSA'].append(None)
                compLats['GCSA'].append(None)
                compLats['GCSB'].append(None)
            counter += 1
    
    for key in compMasses:
        compLats[key] = np.array(compLats[key])
        compLons[key] = np.array(compLons[key])
        compMasses[key] = np.array(compMasses[key])
    compTimes = np.array(compTimes)
    compHasGCS = np.array(compHasGCS)
        
    mkr = '.'
    cAB = '#882255'
    cAA = '#332288'
    cMix = '#AA4499'
    c1sat = '#66CCEE'
    aRatio = np.power(10,compMasses['DP3A']) / np.power(10,compMasses['CORA'])
    bRatio = np.power(10,compMasses['DP3B']) / np.power(10,compMasses['CORB'])
    allRatio = np.array([aRatio, bRatio]).reshape(-1)
    allSeps  = np.array([compSeps['DP3A'], compSeps['DP3A']]).reshape(-1)
    goodIdx1 = np.where((allRatio < 99999999) & (allSeps < 80.))[0]
    allRatio = allRatio[goodIdx1]
    allSeps = allSeps[goodIdx1]
    print(np.mean(aRatio), np.median(aRatio))
    print(np.mean(bRatio), np.median(bRatio))
    print(np.median(allRatio), np.percentile(allRatio, 32), np.percentile(allRatio, 68))
    #print(allRatio)
    
    # Position Plot
    fig = plt.figure(figsize=(10,10))
    gs = fig.add_gridspec(7, 7, hspace=0, wspace=0)
    ff1a = fig.add_subplot(gs[0,0]) 
    ff1b = fig.add_subplot(gs[0,1], sharex=ff1a, sharey=ff1a) 
    ff1c = fig.add_subplot(gs[0,2], sharex=ff1a, sharey=ff1a) 
    ff2a = fig.add_subplot(gs[1,0], sharex=ff1a, sharey=ff1a) 
    ff2b = fig.add_subplot(gs[1,1], sharex=ff1a, sharey=ff1a) 
    ff3a = fig.add_subplot(gs[2,0], sharex=ff1a, sharey=ff1a) 
    
    
    f1a = fig.add_subplot(gs[-7,-1]) 
    f2a = fig.add_subplot(gs[-6,-1], sharex=f1a, sharey=f1a) 
    f2b = fig.add_subplot(gs[-6,-2], sharex=f1a, sharey=f1a) 
    f3a = fig.add_subplot(gs[-5,-1], sharex=f1a, sharey=f1a) 
    f3b = fig.add_subplot(gs[-5,-2], sharex=f1a, sharey=f1a) 
    f3c = fig.add_subplot(gs[-5,-3], sharex=f1a, sharey=f1a) 
    f4a = fig.add_subplot(gs[-4,-1], sharex=f1a, sharey=f1a) 
    f4b = fig.add_subplot(gs[-4,-2], sharex=f1a, sharey=f1a) 
    f4c = fig.add_subplot(gs[-4,-3], sharex=f1a, sharey=f1a) 
    f4d = fig.add_subplot(gs[-4,-4], sharex=f1a, sharey=f1a) 
    f5a = fig.add_subplot(gs[-3,-1], sharex=f1a, sharey=f1a) 
    f5b = fig.add_subplot(gs[-3,-2], sharex=f1a, sharey=f1a) 
    f5c = fig.add_subplot(gs[-3,-3], sharex=f1a, sharey=f1a) 
    f5d = fig.add_subplot(gs[-3,-4], sharex=f1a, sharey=f1a) 
    f5e = fig.add_subplot(gs[-3,-5], sharex=f1a, sharey=f1a) 
    f56a = fig.add_subplot(gs[-2,-1], sharex=f1a, sharey=f1a) 
    f56b = fig.add_subplot(gs[-2,-2], sharex=f1a, sharey=f1a) 
    f56c = fig.add_subplot(gs[-2,-3], sharex=f1a, sharey=f1a) 
    f56d = fig.add_subplot(gs[-2,-4], sharex=f1a, sharey=f1a) 
    f56e = fig.add_subplot(gs[-2,-5], sharex=f1a, sharey=f1a) 
    f56ef = fig.add_subplot(gs[-2,-6], sharex=f1a, sharey=f1a) 

    f6a = fig.add_subplot(gs[-1,-1], sharex=f1a, sharey=f1a) 
    f6b = fig.add_subplot(gs[-1,-2], sharex=f1a, sharey=f1a) 
    f6c = fig.add_subplot(gs[-1,-3], sharex=f1a, sharey=f1a) 
    f6d = fig.add_subplot(gs[-1,-4], sharex=f1a, sharey=f1a) 
    f6e = fig.add_subplot(gs[-1,-5], sharex=f1a, sharey=f1a) 
    f6ef = fig.add_subplot(gs[-1,-6], sharex=f1a, sharey=f1a)     
    f6f = fig.add_subplot(gs[-1,-7], sharex=f1a, sharey=f1a) 

    allfs = [f1a, f2a, f2b, f3a, f3b, f3c, f4a, f4b, f4c, f4d, f5a, f5b, f5c, f5d, f5e, f56a, f56b, f56c, f56d, f56e, f56ef, f6a, f6b, f6c, f6d, f6e, f6ef, f6f]
    leftfs = [f1a, f2a, f3a, f4a, f5a, f56a]
    botfs  = [f6b, f6c, f6d, f6e, f6ef, f6f]
    midfs = [f2b, f3b, f3c, f4b, f4c, f4d, f5b, f5c, f5d, f5e, f56b, f56c, f56d, f56e, f56ef]
    edges = [f1a, f2b, f3c, f4d, f5e, f56ef, f6f]
    
    allffs = [ff1a, ff1b, ff1c, ff2a, ff2b, ff3a]
    
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
    ff3a.set_ylabel('Lon ($^{\\circ}$)')
    
    ff1a.xaxis.set_label_position("top")
    ff1a.xaxis.tick_top()
    ff1a.tick_params(labelbottom=False)
    ff1a.set_xlabel('Lon ($^{\\circ}$)')
    
    ff1b.xaxis.set_label_position("top")
    ff1b.xaxis.tick_top()
    ff1b.tick_params(labelbottom=False)
    ff1b.tick_params(labelleft=False)
    ff1b.set_xlabel('Lon ($^{\\circ}$)')

    ff1c.xaxis.set_label_position("top")
    ff1c.xaxis.tick_top()
    ff1c.tick_params(labelbottom=False)
    ff1c.tick_params(labelleft=False)
    ff1c.set_xlabel('Lon ($^{\\circ}$)')

    ff2a.tick_params(labelbottom=False)
    ff2a.xaxis.tick_top()
    
    ff2b.tick_params(labelbottom=False)
    ff2b.tick_params(labelleft=False)

    ff3a.tick_params(labelbottom=False)
    ff3a.xaxis.tick_top()
    
    # this is for the lons only     
    pairs = [[compLons['DP3A'][compHasGCS], compLons['GCSA'][compHasGCS]], [compLons['DP3A'][compHasLB], compLons['DP1A'][compHasLB]], [compLons['DP3A'], compLons['DP2A']],     [compLons['DP2A'][compHasGCS], compLons['GCSA'][compHasGCS]], [compLons['DP2A'][compHasLB], compLons['DP1A'][compHasLB]],    [compLons['DP1A'][compHasLBGCS], compLons['GCSA'][compHasLBGCS]]]
    
    
    i = 0
    
    ll, ul = -50, 410
    ff1a.set_xlim(ll,ul)
    ff1a.set_ylim(ll,ul)
    ff1a.set_aspect('equal')
    #ff1a.set_xticks([0,20,40])

    dlon = 30
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
        ff.plot([ll, ul], [ll+dlon, ul+dlon], 'k--', zorder=0, alpha=0.4)
        ff.plot([ll, ul], [ll-dlon, ul-dlon], 'k--', zorder=0, alpha=0.4)
        ff.scatter(x,y, marker=mkr, c=c1sat)
        i +=1
    
    ff1c.text(1.5, 0.5, 'DP3', transform=ff1c.transAxes, fontsize=12, ha='center', va='center')
    ff2b.text(1.5, 0.5, 'DP2', transform=ff2b.transAxes, fontsize=12, ha='center', va='center')
    ff3a.text(1.5, 0.5, 'DP1', transform=ff3a.transAxes, fontsize=12, ha='center', va='center')
    ff3a.text(0.5, -0.2, 'GCS', transform=ff3a.transAxes, fontsize=12, ha='center', va='center')
    
    
    ll, ul = -70, 70
    f1a.set_xlim(ll,ul)
    f1a.set_ylim(ll,ul)
    f1a.set_aspect('equal')
    dlat = 10
    for f in allfs:
        f.plot([ll, ul], [ll, ul], 'k--', zorder=0, alpha=0.75)
        f.plot([ll, ul], [ll+dlat, ul+dlat], 'k--', zorder=0, alpha=0.4)
        f.plot([ll, ul], [ll-dlat, ul-dlat], 'k--', zorder=0, alpha=0.4)
    
    #toptits = ['CORSET B', 'GCS', 'DP2A', 'DP2B' 'DP3A', 'DP3B']    
    righttits = ['GCS', 'DP1', 'DP2A', 'DP2B', 'DP3A', 'DP3B',  'CORSET\nA']    
    i = 0
    for f in edges:
        if i != 5:
            f.text(-0.5, 0.5, righttits[i], transform=f.transAxes, fontsize=12, ha='center', va='center')
        else:
            f.text(-0.3, 0.5, righttits[i], transform=f.transAxes, fontsize=12, ha='center', va='center')
        i += 1
    f1a.set_title('CORSET B', fontsize=12)


    f1a.scatter(compLats['CORB'], compLats['GCSA'], marker=mkr, c=c1sat)

    f2a.scatter(compLats['CORB'], compLats['DP1A'], marker=mkr, c=c1sat)
    f2b.scatter(compLats['GCSA'], compLats['DP1A'], marker=mkr, c=c1sat)

    f3a.scatter(compLats['CORB'], compLats['DP2A'], marker=mkr, c=cAB)
    f3b.scatter(compLats['GCSA'], compLats['DP2A'], marker=mkr, c=c1sat)
    f3c.scatter(compLats['DP1A'], compLats['DP2A'], marker=mkr, c=c1sat)

    f4a.scatter(compLats['CORB'], compLats['DP2B'], marker=mkr, c=cAA)
    f4b.scatter(compLats['GCSA'], compLats['DP2B'], marker=mkr, c=c1sat)
    f4c.scatter(compLats['DP1A'], compLats['DP2B'], marker=mkr, c=c1sat)
    f4d.scatter(compLats['DP2A'], compLats['DP2B'], marker=mkr, c=cAB)

    f5a.scatter(compLats['CORB'], compLats['DP3A'], marker=mkr, c=cMix)
    f5b.scatter(compLats['GCSA'], compLats['DP3A'], marker=mkr, c=c1sat)
    f5c.scatter(compLats['DP1A'], compLats['DP3A'], marker=mkr, c=c1sat)
    f5d.scatter(compLats['DP2A'], compLats['DP3A'], marker=mkr, c=cAA)
    f5e.scatter(compLats['DP2B'], compLats['DP3A'], marker=mkr, c=cMix)

    f56a.scatter(compLats['CORB'], compLats['DP3B'], marker=mkr, c=cAA)
    f56b.scatter(compLats['GCSA'], compLats['DP3B'], marker=mkr, c=c1sat)
    f56c.scatter(compLats['DP1A'], compLats['DP3B'], marker=mkr, c=c1sat)
    f56d.scatter(compLats['DP2A'], compLats['DP3B'], marker=mkr, c=cMix)
    f56e.scatter(compLats['DP2B'], compLats['DP3B'], marker=mkr, c=cAA)
    f56ef.scatter(compLats['DP3A'], compLats['DP3B'], marker=mkr, c=cAB)


    f6a.scatter(compLats['CORB'], compLats['CORA'], marker=mkr, c=cAB)
    f6b.scatter(compLats['GCSA'], compLats['CORA'], marker=mkr, c=c1sat)
    f6c.scatter(compLats['DP1A'], compLats['CORA'], marker=mkr, c=c1sat)
    f6d.scatter(compLats['DP2A'], compLats['CORA'], marker=mkr, c=cAA)
    f6e.scatter(compLats['DP2B'], compLats['CORA'], marker=mkr, c=cMix)
    f6ef.scatter(compLats['DP3A'], compLats['CORA'], marker=mkr, c=cAA)
    f6f.scatter(compLats['DP3B'], compLats['CORA'], marker=mkr, c=cMix)


    f6f.text(0.03, 0.47, 'Same Method, Diff Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cAB, )
    f6f.text(0.03, 0.45, 'Diff Method, Same Sat', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cAA)
    f6f.text(0.03, 0.43, 'Diff Method, Diff Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cMix)
    f6f.text(0.03, 0.41, 'One Value For Both Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=c1sat)
    plt.savefig(myFold+'MegaPosComp.png')
    
    '''lonKeys = ['GCSA', 'DP1A', 'DP2A', 'DP3A']
    for key1 in lonKeys:
        print ('')
        for key2 in lonKeys:
            if key1 != key2:
                lon1 =  compLons[key1].astype(float)
                lon2 =  compLons[key2].astype(float)
                lon1[np.isnan(lon1)] = None
                lon2[np.isnan(lon2)] = None
                goodIdx = np.where((lon1 != None) & (~np.isnan(lon1)) & (lon2 != None) & (~np.isnan(lon2)))
                x = lon1[goodIdx]
                y = lon2[goodIdx]
                lonDifs = x - y
                for j in range(len(lonDifs)):
                    if np.abs(lonDifs[j]) > 180:
                        if x[j] < 180:
                            y[j] -= 360
                        else:
                            y[j] += 360
                
                err  = np.abs(x-y)
                #print (err)
                print (key1, key2, np.median(err), stats.pearsonr(x,y))'''
                
    latKeys = ['CORA', 'CORB', 'GCSA', 'DP1A', 'DP2A', 'DP2B', 'DP3A', 'DP3B']
    for key1 in latKeys:
        print ('')
        for key2 in latKeys:
            if key1 != key2:
                lat1 =  compLats[key1].astype(float)
                lat2 =  compLats[key2].astype(float)
                goodIdx = np.where((lat1 != None) & (~np.isnan(lat1)) & (lat2 != None) & (~np.isnan(lat2)))
                x = lat1[goodIdx]
                y = lat2[goodIdx]
                err  = np.abs(x-y)
                #print (err)
                #print (key1, key2, np.median(err), stats.pearsonr(x,y))
    
    
    
    # Mass Plot
    fig = plt.figure(figsize=(10,10))
    gs = fig.add_gridspec(8, 8, hspace=0, wspace=0)
    f1a = fig.add_subplot(gs[0,0]) 
    f2a = fig.add_subplot(gs[1,0], sharex=f1a, sharey=f1a) 
    f2b = fig.add_subplot(gs[1,1], sharex=f1a, sharey=f1a) 
    
    f22a = fig.add_subplot(gs[2,0], sharex=f1a, sharey=f1a) 
    f22b = fig.add_subplot(gs[2,1], sharex=f1a, sharey=f1a) 
    f22c = fig.add_subplot(gs[2,2], sharex=f1a, sharey=f1a) 
    
    f23a = fig.add_subplot(gs[3,0], sharex=f1a, sharey=f1a) 
    f23b = fig.add_subplot(gs[3,1], sharex=f1a, sharey=f1a) 
    f23c = fig.add_subplot(gs[3,2], sharex=f1a, sharey=f1a) 
    f23d = fig.add_subplot(gs[3,3], sharex=f1a, sharey=f1a) 
    
    f3a = fig.add_subplot(gs[4,0], sharex=f1a, sharey=f1a) 
    f3b = fig.add_subplot(gs[4,1], sharex=f1a, sharey=f1a) 
    f3c = fig.add_subplot(gs[4,2], sharex=f1a, sharey=f1a) 
    f3cc = fig.add_subplot(gs[4,3], sharex=f1a, sharey=f1a) 
    f3cd = fig.add_subplot(gs[4,4], sharex=f1a, sharey=f1a) 

    f4a = fig.add_subplot(gs[5,0], sharex=f1a, sharey=f1a) 
    f4b = fig.add_subplot(gs[5,1], sharex=f1a, sharey=f1a) 
    f4c = fig.add_subplot(gs[5,2], sharex=f1a, sharey=f1a) 
    f4cc = fig.add_subplot(gs[5,3], sharex=f1a, sharey=f1a) 
    f4cd = fig.add_subplot(gs[5,4], sharex=f1a, sharey=f1a) 
    f4d = fig.add_subplot(gs[5,5], sharex=f1a, sharey=f1a) 

    f5a = fig.add_subplot(gs[6,0], sharex=f1a, sharey=f1a) 
    f5b = fig.add_subplot(gs[6,1], sharex=f1a, sharey=f1a) 
    f5c = fig.add_subplot(gs[6,2], sharex=f1a, sharey=f1a) 
    f5cc = fig.add_subplot(gs[6,3], sharex=f1a, sharey=f1a) 
    f5cd = fig.add_subplot(gs[6,4], sharex=f1a, sharey=f1a) 
    f5d = fig.add_subplot(gs[6,5], sharex=f1a, sharey=f1a) 
    f5e = fig.add_subplot(gs[6,6], sharex=f1a, sharey=f1a) 

    f6a = fig.add_subplot(gs[7,0], sharex=f1a, sharey=f1a) 
    f6b = fig.add_subplot(gs[7,1], sharex=f1a, sharey=f1a) 
    f6c = fig.add_subplot(gs[7,2], sharex=f1a, sharey=f1a) 
    f6cc = fig.add_subplot(gs[7,3], sharex=f1a, sharey=f1a) 
    f6cd = fig.add_subplot(gs[7,4], sharex=f1a, sharey=f1a) 
    f6d = fig.add_subplot(gs[7,5], sharex=f1a, sharey=f1a) 
    f6e = fig.add_subplot(gs[7,6], sharex=f1a, sharey=f1a) 
    f6f = fig.add_subplot(gs[7,7], sharex=f1a, sharey=f1a) 

    allfs = [f1a, f2a, f2b, f22a, f22b, f22c, f23a, f23b, f23c, f23d, f3a, f3b, f3c, f3cc, f3cd, f4a, f4b, f4c, f4cc, f4cd, f4d, f5a, f5b, f5c, f5cc, f5cd, f5d, f5e, f6a, f6b, f6c, f6cc, f6cd, f6d, f6e, f6f]
    leftfs = [f1a, f2a, f22a, f23a, f3a, f4a, f5a]
    botfs  = [f6b, f6c, f6cc, f6cd, f6d, f6e, f6f]
    midfs = [f2b, f22b, f22c, f23b, f23c, f23d, f3b, f3c, f3cc, f3cd, f4b, f4c, f4cc, f4cd, f4d, f5b, f5c, f5cc, f5cd, f5d, f5e]
    edges = [f1a, f2b, f22c, f23d, f3cd, f4d, f5e, f6f]

    ll, ul = 13, 18
    f1a.set_xlim(ll,ul)
    f1a.set_ylim(ll,ul)
    f1a.set_aspect('equal')
    f1a.set_xticks([14,15,16,17])
    f1a.set_yticks([14,15,16,17])
    f1a.set_xticklabels([14,'',16,''])
    
    
    for key in compMasses:
        idx = np.where(compMasses[key] == None)
        compMasses[key][idx] = 0
        compMasses[key] = (compMasses[key]* 1e15).astype(float)
        idx = np.where(compMasses[key] > 0)
        compMasses[key][idx] = np.log10(np.array(compMasses[key][idx]))
    
    
    for f in leftfs:
        f.tick_params(labelbottom=False)
        f.set_ylabel('log M (g)')
    for f in botfs:
        f.tick_params(labelleft=False)
        f.set_xlabel('log M (g)')
    for f in midfs:
        f.tick_params(labelbottom=False, labelleft=False)
    f6a.set_xlabel('log M (g)')
    f6a.set_ylabel('log M (g)')
    
    for f in allfs:
        f.plot([ll, ul], [ll, ul], 'k--', zorder=0, alpha=0.5)
        
    #toptits = ['CORSET B', 'GCS A', 'GCS B', 'Deproj1', 'Deproj2 A', 'deproj2 B']    
    righttits = ['GCS A', 'GCSB', 'DP1A', 'DP1B', 'DP2', 'DP3A', 'DP3B', 'CORSET\nA']    
    i = 0
    for f in edges:
        if i != 7:
            f.text(1.5, 0.5, righttits[i], transform=f.transAxes, fontsize=12, ha='center', va='center')
        else:
            f.text(1.4, 0.5, righttits[i], transform=f.transAxes, fontsize=12, ha='center', va='center')
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
    
 
    f22a.scatter(compMasses['CORB'], compMasses['DP1A'], marker=mkr, c=cMix)
    f22b.scatter(compMasses['GCSA'], compMasses['DP1A'], marker=mkr, c=cAA)
    f22c.scatter(compMasses['GCSB'], compMasses['DP1A'], marker=mkr, c=cMix)

    f23a.scatter(compMasses['CORB'], compMasses['DP1B'], marker=mkr, c=cAA)
    f23b.scatter(compMasses['GCSA'], compMasses['DP1B'], marker=mkr, c=cMix)
    f23c.scatter(compMasses['GCSB'], compMasses['DP1B'], marker=mkr, c=cAA)
    f23d.scatter(compMasses['DP1A'], compMasses['DP1B'], marker=mkr, c=cAB)
    
    

    f3a.scatter(compMasses['CORB'], compMasses['DP2A'], marker=mkr, c= c1sat)
    f3b.scatter(compMasses['GCSA'], compMasses['DP2A'], marker=mkr, c= c1sat)
    f3c.scatter(compMasses['GCSB'], compMasses['DP2A'], marker=mkr, c= c1sat)
    f3cc.scatter(compMasses['DP1A'], compMasses['DP2A'], marker=mkr, c= c1sat)
    f3cd.scatter(compMasses['DP1B'], compMasses['DP2A'], marker=mkr, c= c1sat)

    f4a.scatter(compMasses['CORB'], compMasses['DP3A'], marker=mkr, c=cMix)
    f4b.scatter(compMasses['GCSA'], compMasses['DP3A'], marker=mkr, c=cAA)
    f4c.scatter(compMasses['GCSB'], compMasses['DP3A'], marker=mkr, c=cMix)
    f4cc.scatter(compMasses['DP1A'], compMasses['DP3A'], marker=mkr, c= cAA)
    f4cd.scatter(compMasses['DP1B'], compMasses['DP3A'], marker=mkr, c= cMix)
    f4d.scatter(compMasses['DP2A'], compMasses['DP3A'], marker=mkr, c= c1sat)

    f5a.scatter(compMasses['CORB'], compMasses['DP3B'], marker=mkr, c=cAA)
    f5b.scatter(compMasses['GCSA'], compMasses['DP3B'], marker=mkr, c=cMix)
    f5c.scatter(compMasses['GCSB'], compMasses['DP3B'], marker=mkr, c=cAA)
    f5cc.scatter(compMasses['DP1A'], compMasses['DP3B'], marker=mkr, c= cMix)
    f5cd.scatter(compMasses['DP1B'], compMasses['DP3B'], marker=mkr, c= cAA)
    f5d.scatter(compMasses['DP2A'], compMasses['DP3B'], marker=mkr, c= c1sat)
    f5e.scatter(compMasses['DP3A'], compMasses['DP3B'], marker=mkr, c=cAB)

    f6a.scatter(compMasses['CORB'], compMasses['CORA'], marker=mkr, c=cAB)
    f6b.scatter(compMasses['GCSA'], compMasses['CORA'], marker=mkr, c=cAA)
    f6c.scatter(compMasses['GCSB'], compMasses['CORA'], marker=mkr, c=cMix)
    f6cc.scatter(compMasses['DP1A'], compMasses['CORA'], marker=mkr, c= cAA)
    f6cd.scatter(compMasses['DP1B'], compMasses['CORA'], marker=mkr, c= cMix)
    f6d.scatter(compMasses['DP2A'], compMasses['CORA'], marker=mkr, c= c1sat)
    f6e.scatter(compMasses['DP3A'], compMasses['CORA'], marker=mkr, c=cAA)
    f6f.scatter(compMasses['DP3B'], compMasses['CORA'], marker=mkr, c=cMix)
    
    f6f.text(0.6, 0.75, 'Same Method, Diff Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cAB)
    f6f.text(0.6, 0.73, 'Diff Method, Same Sat', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cAA)
    f6f.text(0.6, 0.71, 'Diff Method, Diff Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=cMix)
    f6f.text(0.6, 0.69, 'Single Value For Both Sats', transform=fig.transFigure, fontsize=12, ha='left', va='center', c=c1sat)

    plt.savefig(myFold+'MegaMassComp.png')
    
    '''for key1 in compMasses:
        print ('')
        for key2 in compMasses:
            M1 =  np.power(10,compMasses[key1])
            M1[M1 == None] = 0
            M2 = np.power(10,compMasses[key2])
            M2[M2 == None] = 0
            avgM = (M1 + M2) * 0.5
            goodIdx = np.where((M1 > 0) & (M2 >0))[0]
            percErr = '{:.2f}'.format(100*np.median(np.abs(M1[goodIdx] - M2[goodIdx]) / avgM[goodIdx]))
            print (key1, key2, percErr, stats.pearsonr(np.log10(M1[goodIdx]), np.log10(M2[goodIdx])))'''
    
# Plot the recon lat/lon onto solar disk
def circlePlot():
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
            compMasses['CORB'].append(np.max(thisRes.CORSETmassMB))
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
                if (isinstance(thisRes.GCSmassesA,list)) & (isinstance(thisRes.GCSmassesB, list)):
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
            else:
                compMasses['GCSA'].append(None)
                compMasses['GCSB'].append(None)
                compLons['GCSA'].append(None)
                compLats['GCSA'].append(None)
                compLats['GCSB'].append(None)
            counter += 1
    
    for key in compLats:
        compLats[key] = np.array(compLats[key])
        compLons[key] = np.array(compLons[key])
        compMasses[key] = np.array(compMasses[key])
    compTimes = np.array(compTimes)
    compHasGCS = np.array(compHasGCS)
        
    dtor = np.pi / 180.
    fig = plt.figure(constrained_layout=True, figsize=(7,11))
    gs = fig.add_gridspec(3, 2)
    f1a = fig.add_subplot(gs[0,0]) 
    f1b = fig.add_subplot(gs[0,1], sharex=f1a, sharey=f1a) 
    f2a = fig.add_subplot(gs[1,0], sharex=f1a, sharey=f1a) 
    f2b = fig.add_subplot(gs[1,1], sharex=f1a, sharey=f1a) 
    f3a = fig.add_subplot(gs[2,0], sharex=f1a, sharey=f1a) 
    f3b = fig.add_subplot(gs[2,1], sharex=f1a, sharey=f1a) 
   
    allfs = [f1a, f1b, f2a, f2b, f3a, f3b]
    lefts = [f1a,  f2a, f3a]
    rights = [f1b, f2b, f3b]
   
    deg10s = np.array([np.pi/2 - 15*dtor*i for i in range(13)])
   
    thetas = np.linspace(0, 2*np.pi, 180)
    cxs = np.cos(thetas)
    czs = np.sin(thetas)
    for f in allfs:
        f.plot(cxs, czs, 'k')
        for i in range(len(deg10s)):
            myy = np.sin(deg10s[i])
            myx = np.cos(deg10s[i])
            f.plot([-myx, myx], [myy, myy], 'k--', alpha=0.1, zorder=0)
            f.plot([myy, myy], [-myx, myx], 'k--', alpha=0.1, zorder=0)
            
        f.axis('off')
        f.set_aspect('equal')
    
    for f in lefts:
        f.text(-0.05, 0.5, 'E', transform=f.transAxes, fontsize=16, ha='left', va='center', c='k')
        f.text(1.18, 0.5, 'W', transform=f.transAxes, fontsize=16, ha='center', va='center', c='k')
    for f in rights:
        f.text(1.05, 0.5, 'E', transform=f.transAxes, fontsize=16, ha='left', va='center', c='k')
    for f in allfs:
        f.text(0.5, 1.05, 'N', transform=f.transAxes, fontsize=16, ha='center', va='center', c='k')
        f.text(0.5, -0.05, 'S', transform=f.transAxes, fontsize=16, ha='center', va='center', c='k')
    f1a.text(1.18, 1.05, 'GCS', transform=f1a.transAxes, fontsize=16, ha='center', va='center', c='#696969')
    f2a.text(1.18, 1.05, 'Deproj1', transform=f2a.transAxes, fontsize=16, ha='center', va='center', c='#6f10be')
    f3a.text(1.18, 1.05, 'Deproj2', transform=f3a.transAxes, fontsize=16, ha='center', va='center', c='#367cc3')
        
    
   
    # GCS
    '''mod = 'GCSA'
    lats = compLats[mod][compHasGCS]
    lons = compLons[mod][compHasGCS]'''
    lats = []
    lons = []
    # check if has GCS values
    for key in GCSkeys:
        thisRes = res[key]
        lonsG = []
        latsG = []
        
        for gKey in thisRes.GCSvals.keys():
            lonsG.append((float(thisRes.GCSvals[gKey][1])+360)%360)
            latsG.append(float(thisRes.GCSvals[gKey][0]))
        if (isinstance(thisRes.GCSmassesA,list)) & (isinstance(thisRes.GCSmassesB, list)):    
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
        lons.append(meanLonG)
        lats.append(np.mean(latsG))
    lats = np.array(lats)
    lons = np.array(lons)

        
        
    xs, ys, zs = [], [], []
    for i in range(len(lats)):
        xs.append(np.cos(lons[i]*dtor)*np.cos(lats[i]*dtor))
        zs.append(np.sin(lats[i]*dtor))
        ys.append(np.sin(lons[i]*dtor)*np.cos(lats[i]*dtor))
    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    front = np.where(ys >= 0)[0]
    back = np.where(ys < 0)[0]
    f1a.scatter(xs[front], zs[front], c='#696969')    
    f1b.scatter(-xs[back], zs[back], c='#696969')
   
   
    # Deproj 1
    mod = 'DP1A'
    lats = compLats[mod]
    lons = compLons[mod]
    xs = np.cos(lons*dtor)*np.cos(lats*dtor)
    zs = np.sin(lats*dtor)
    ys = np.sin(lons*dtor)*np.cos(lats*dtor)
    front = np.where(ys >= 0)[0]
    back = np.where(ys < 0)[0]
    f2a.scatter(xs[front], zs[front], c='#6f10be')    
    f2b.scatter(-xs[back], zs[back], c='#6f10be')
   
    # Deproj 2
    mod = 'DP2A'
    lats = compLats[mod]
    lons = compLons[mod]
    xs = np.cos(lons*dtor)*np.cos(lats*dtor)
    zs = np.sin(lats*dtor)
    ys = np.sin(lons*dtor)*np.cos(lats*dtor)
    front = np.where(ys >= 0)[0]
    back = np.where(ys < 0)[0]
    f3a.scatter(xs[front], zs[front], c='#367cc3')    
    f3b.scatter(-xs[back], zs[back], c='#367cc3')
    #plt.show()
    plt.savefig(myFold+'CirclePlots.png')
    
#aMassEvol('20120712_165400')
#allMassEvols() # loop aMassEvol over all keys 
#make3timeline() 
#make4timeline() 
#megaHisto()
#A2Bcomp()
#C2Gcomp()
#sepPlot()
#coveragePlot(60)
MLLscatters()
#circlePlot()
