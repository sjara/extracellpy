#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Some plotting functions.

plotintervalbars: Plot with thin errorbars given by interval (not errors). 
'''

import matplotlib.pyplot as mpl
import numpy as np

def boxoff(ax,keep='left'):
    ax.spines['top'].set_visible(False)
    if keep=='left':
        ax.spines['right'].set_visible(False)
    else:
        ax.spines['left'].set_visible(False)        
    xtlines = ax.get_xticklines()
    ytlines = ax.get_yticklines()
    for t in xtlines[1::2]+ytlines[1::2]:
        t.set_visible(False)

def set_log_ticks(ax,xTickValues,axis='x'):
    xTickLogValues = np.log10(xTickValues);
    xTickLabels = ['%d'%(1e-3*xt) for xt in xTickValues];
    if axis=='x':
        ax.set_xticks(xTickLogValues)
        ax.set_xticklabels(xTickLabels)
    else:
        ax.set_yticks(xTickLogValues)
        ax.set_yticklabels(xTickLabels)

def set_logy_ticks(ax,xTickValues):
    xTickLogValues = np.log10(xTickValues);
    xTickLabels = ['%d'%(1e-3*xt) for xt in xTickValues];
    ax.set_yticks(xTickLogValues)
    ax.set_yticklabels(xTickLabels)


def plotintervalbars(Xvalues,Yvalues,Yinterval,PlotColor='k'):
    '''Show a plot with big markers and thin errorbars.
    (hpmain,hpcap,hperr) = plotintervalbars(Xvalues,Yvalues,Yinterval,PlotColor)
    
    Xvalues, Yvalues must be arrays of shape (Nvalues,)
    Yinterval must be an array of shape (Nvalues,2)
    PlotColor must be a matplotlib color
    '''
    #Yerr = np.c_[[1,-1]]*Yvalues + np.c_[[-1,1]]*Yinterval
    Yerr = Yvalues[:,np.newaxis]*np.r_[1,-1] + Yinterval*np.r_[-1,1]
    (hpmain,hpcap,hperr)=mpl.errorbar(Xvalues,Yvalues,yerr=Yerr.T)
    mpl.setp(hpmain,marker='o',linewidth=3,markersize=10,markeredgewidth=3)
    mpl.setp(hpmain,mfc='w',mec=PlotColor,color=PlotColor)
    mpl.setp(hperr,color=PlotColor)
    mpl.setp(hpcap,color=PlotColor)

    return (hpmain,hpcap,hperr)

'''
def stairshistogram(VarValue,BinsEdges,Nvalues=None,Xlims=None):
    (nvals, lowedges) = histogram(VarValue,BinsEdges)
    # New version uses bin edges instead of the left bin edges.
    #(nvals, lowedges) = histogram(VarValue,BinsEdges,new=True) # New version
    if(Nvalues==None):
        Nvalues = len(VarValue)
    if(Xlims==None):
        SamplesToPlot= range(len(nvals))
    else:
        SamplesToPlot = (lowedges>=Xlims[0])&(lowedges<=Xlims[-1])
    hp = stairsplot(lowedges[SamplesToPlot],\
                    nvals[SamplesToPlot]/float(len(VarValue)))
    #print (sum(nvals/float(len(VarValue))))
    return hp
'''

def stairsplot(x,y):
    '''Plot stairs. A horizontal line stays at the same value until next Y point'''
    hp = mpl.plot(np.repeat(x,2)[1:],np.repeat(y,2)[0:-1])
    return hp


def rasterplot(spikeTimesFromEventOnset,indexLimitsEachTrial,timeRange,trialsEachCond=[],
               colorEachCond=None,fillWidth=50):
    '''Plot spikes raster plot grouped by condition

    return (pRaster,hcond)

    First trial is plotted at y=0
    '''
    import extrafuncs
    import matplotlib.pyplot as plt

    #fillWidth = 50
    if trialsEachCond==[]:
        nCond=1
        trialsEachCond = [np.arange(indexLimitsEachTrial.shape[1])]
    else:
        nCond = len(trialsEachCond)
    
    if colorEachCond is None:
        colorEachCond = ['0.5']*nCond
        #colorEachCond = ['0.5']*nCond

    nSpikesEachTrial = np.diff(indexLimitsEachTrial,axis=0)[0]
    nSpikesEachTrial = nSpikesEachTrial*(nSpikesEachTrial>0) # Some are negative
    trialIndexEachCond = []
    spikeTimesEachCond = []
    for indcond,trialsThisCond in enumerate(trialsEachCond):
        spikeTimesThisCond = np.empty(0,dtype='float64')
        trialIndexThisCond = np.empty(0,dtype='int')
        for indtrial,thisTrial in enumerate(trialsThisCond):
            indsThisTrial = slice(indexLimitsEachTrial[0,thisTrial],
                                  indexLimitsEachTrial[1,thisTrial])
            spikeTimesThisCond = np.concatenate((spikeTimesThisCond,
                                                 spikeTimesFromEventOnset[indsThisTrial]))
            trialIndexThisCond = np.concatenate((trialIndexThisCond,
                                                 np.repeat(indtrial,nSpikesEachTrial[thisTrial])))
        trialIndexEachCond.append(np.copy(trialIndexThisCond))
        spikeTimesEachCond.append(np.copy(spikeTimesThisCond))

    nTrialsEachCond = [len(x) for x in trialsEachCond]
    #pR = plot(1e3*spikeTimesFromEventOnset,trialIndexForEachSpike,'.k')
    #pR[0].set_markersize(markerSize)
    xpos = timeRange[0]+np.array([0,fillWidth,fillWidth,0])
    lastTrialEachCond = np.cumsum(nTrialsEachCond)
    firstTrialEachCond = np.r_[0,lastTrialEachCond[:-1]-1]

    hcond=[]
    pRaster = []
    for indcond in range(nCond):
        pRasterOne, = plt.plot(spikeTimesEachCond[indcond],
                            trialIndexEachCond[indcond]+firstTrialEachCond[indcond],'.k')
        pRaster.append(pRasterOne)
        plt.hold(True)
        ypos = [firstTrialEachCond[indcond],firstTrialEachCond[indcond],
                lastTrialEachCond[indcond],lastTrialEachCond[indcond]]
        hcond.extend(plt.fill(xpos,ypos,ec='none',fc=colorEachCond[indcond]))
        hcond.extend(plt.fill(xpos+np.diff(timeRange)-fillWidth,ypos,ec='none',fc=colorEachCond[indcond]))
    plt.hold(False)
    return(pRaster,hcond)

def plot_psth(PSTH,smoothWinSize,binsStartTime,colorEachCond=None,linewidth=3):
    import matplotlib.pyplot as plt
    from scipy.signal import hanning
    winShape = np.ones(smoothWinSize) # Square
    #winShape = hanning(smoothWinSize) # Hanning
    winShape = winShape/np.sum(winShape)
    nCond = PSTH.shape[0]
    if colorEachCond is None:
        colorEachCond = ['0.5']*nCond
    pPSTH = []
    for (indc,thisPSTH) in enumerate(PSTH):
        smoothPSTH = np.convolve(thisPSTH,winShape,mode='same')
        sSlice = slice(0,len(smoothPSTH),50) # 50
        ph, = plt.plot(1e3*binsStartTime[sSlice],smoothPSTH[sSlice],'-')
        pPSTH.append(ph)
        pPSTH[-1].set_linewidth(linewidth)
        pPSTH[-1].set_color(colorEachCond[indc])    
        plt.hold(True)
    return pPSTH

def loop_keycontrol(fig,indrange,ind):
    '''Block figure until key-press and return next (or previous) index
    depending on key (left/right arrows, or Q to quit).
    '''
    from matplotlib import blocking_input
    bi = blocking_input.BlockingInput(fig,eventslist=('key_press_event',))
    while True:
        bi(timeout=-1)
        e=bi.events[-1]
        #print e.key
        if e.key in ['q','escape']:
            raise StopIteration
        elif e.key in ['right','.','d']:
            newind = min(ind+1,indrange[-1])
            break
        elif e.key in ['left',',','a']:
            newind = max(ind-1,indrange[0])
            break
    return newind

'''
THE FOLLOWING CODE NEEDS trialIndexForEachSpike AS PARAMETER.
    trialIndexEachCond = []
    spikeTimesEachCond = []
    for indcond,trialsThisCond in enumerate(trialsEachCond):
        spikeIndsThisCond = extrafuncs.ismember(trialIndexForEachSpike,trialsThisCond)
        tt = trialIndexForEachSpike[spikeIndsThisCond]
        trialIndexThisCond = np.cumsum(np.r_[0,np.diff(tt)]>0)
        ###########
        ### FIXME: what happens if a trial has no spikes???
        ###########
        spikeTimesThisCond = spikeTimesFromEventOnset[spikeIndsThisCond]
        trialIndexEachCond.append(trialIndexThisCond)
        spikeTimesEachCond.append(spikeTimesThisCond)

'''


def plot_index_histogram(val1,val2,pValue,nBins=14,fontSize=12):
    MI = (val1-val2)/(val1+val2)
    MIsignif = MI[pValue<0.05]
    if np.any(np.isnan(MI)):
        print 'MI contains NaN'

    ax = mpl.gca()#mpl.axes([0.2,0.2,0.6,0.6])
    histRange = [-1,1]
    mpl.hist(MI,bins=nBins,range=histRange,color='0.75')
    mpl.hold(True)
    mpl.hist(MIsignif,bins=nBins,range=histRange,color='0')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xtlines = ax.get_xticklines()
    ytlines = ax.get_yticklines()
    for t in xtlines[1::2]+ytlines[1::2]:
        t.set_visible(False)
    yLims = mpl.ylim()
    mpl.plot([0,0],[0,yLims[1]],color='k',linestyle=':')
    mpl.hold(False)
    mpl.xlabel('Index: (x1-x2)/(x1+x2)',fontsize=fontSize)
    mpl.ylabel('Number of cells',fontsize=fontSize)
    mpl.setp(ax.get_xticklabels(),fontsize=fontSize-2)
    mpl.setp(ax.get_yticklabels(),fontsize=fontSize-2)
    ax.set_yticks(range(0,yLims[1]+1,10))
    mpl.draw()
    mpl.show()

def plot_scatter(val1,val2,pValue,fontSize=12):
    hp = mpl.plot(val1,val2,'o',mfc='none',mec='0.5')
    signifPoints = pValue<0.05
    mpl.hold(True)
    hp = mpl.plot(val1[signifPoints],val2[signifPoints],'o',mfc='k',mec='k')
    maxRate = np.nanmax([val1,val2])
    hline = mpl.plot([0,maxRate],[0,maxRate],'0.0',linestyle='--')
    mpl.hold(False)
    mpl.axis([0,maxRate,0,maxRate])
    #mpl.axis('equal')
    mpl.xlabel('x1')
    mpl.ylabel('x2')
    mpl.title('N=%d (%0.0f%% p<0.05)'%(len(pValue),100*float(sum(signifPoints))/len(pValue)))
    mpl.draw()
    mpl.show()



if __name__ == "__main__":

    Xvalues = np.arange(10)
    Yvalues = Xvalues**1.2
    #Yerr = np.c_[[-1,1]]*np.ones((2,10))
    Yerr = np.ones((10,2))*[-1,1]
    Yinterval = Yvalues[:,np.newaxis] + Yerr
    PlotColor = 'y'
    (hpmain,hpcap,hperr) = plotintervalbars(Xvalues,Yvalues,Yinterval,PlotColor)
    mpl.show()
    mpl.draw()
    
