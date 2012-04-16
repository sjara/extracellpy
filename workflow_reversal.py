'''
Workflow for analysis of extracellular data from reversal task
'''

import numpy as np
from extracellpy import settings
reload(settings)
from extracellpy import sessionanalysis
reload(sessionanalysis)
from extracellpy import spikesanalysis
reload(spikesanalysis)
from extracellpy import extraplots
import os
import sys
import matplotlib.pyplot as plt

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

RASTER_MARKERSIZE = 2

def save_raster_data(animalName,lockedTo=None):
    '''
    Save spike times aligned to a task event.
    Data directories are defined in .../extracellpy/settings.py
    
    Parameters
    ----------
    animalName: string. Name of animal.For example 'saja000'
    lockedTo  : string. Either 'SoundOn', 'Cout' or 'SideIn'
                If empty, all of them are evaluated.
    '''
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    dataModule = 'allcells_%s'%(animalName)
    allcells = __import__(dataModule)
    reload(allcells)
    outputPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
    # -- Define events to lock data to --
    if isinstance(lockedTo,str):
        lockedTo = [lockedTo]
    elif lockedTo is None:
        lockedTo = ['SoundOn', 'Cout', 'SideIn']
    # -- Save results for each event type --
    for oneLock in lockedTo:
        outputDir = os.path.join(outputPath,'lockedTo%s'%(oneLock))
        if oneLock=='SoundOn':
            timeRange = np.array([-0.3,0.9]) #sec
            allcells.cellDB.save_locked_spikes(outputDir,timeRange=timeRange,lockTo=1)
        elif oneLock=='Cout':
            timeRange = np.array([-0.3,0.9]) #sec
            allcells.cellDB.save_locked_spikes(outputDir,timeRange=timeRange,lockTo=2)
        elif oneLock=='SideIn':
            timeRange = np.array([-0.6,0.6]) #sec
            allcells.cellDB.save_locked_spikes(outputDir,timeRange=timeRange,lockTo=3)

def save_raster_plots(animalName):
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    dataModule = 'allcells_%s'%(animalName)
    allcells = __import__(dataModule)
    reload(allcells)

    groupedBy = '4cond'
    #oneLock = 'SoundOn'
    for oneLock in ['SoundOn', 'Cout', 'SideIn']:
        dataPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
        dataDir = os.path.join(dataPath,'lockedTo%s'%(oneLock))
        outputDir = os.path.join(dataPath,'rasters_lockedTo%s_%s'%(oneLock,groupedBy))

        if not os.path.exists(outputDir):
            print 'Creating output directory: %s'%(outputDir)
            os.makedirs(outputDir)

        prevSession = ''
        for indcell,onecell in enumerate(allcells.cellDB):
            cellStr = str(onecell).replace(' ','_')
            fileName = os.path.join(dataDir,cellStr+'_'+oneLock+'.npz')
            # -- Load data saved by save_raster_data() --
            ephysData = np.load(fileName)
            # -- Load behavior --
            if ephysData['behavSession']==prevSession:
                print 'Behavior already loaded (%s)'%cellStr,
            else:
                print 'Loading %s [%s] - %s'%(ephysData['animalName'],ephysData['behavSession'],cellStr),
                behavData = sessionanalysis.load_behavior_session(ephysData['animalName'],
                                                                  ephysData['behavSession'])
                prevSession = ephysData['behavSession']

            plot_raster_reversal(onecell,ephysData,behavData,oneLock,groupedBy)    

            plt.gcf().set_size_inches((8,6))
            figFormat = 'png' #'png' #'pdf' #'svg'
            figName = 'raster_%s.%s'%(cellStr,figFormat)
            plt.savefig(os.path.join(outputDir,figName),format=figFormat)
            print '... figure saved.'
            #plt.draw(); plt.show(); break


def plot_raster_reversal(onecell,ephysData,behavData,lockedTo,groupedBy):
    '''
    NOTE: This function is not well designed. It requires inputs in very specific/complex formats.
    It just makes the code more readable.
    '''
    windowSizePSTH = 0.0025          # sec
    rasterDeltaT = 0.1e-3            # sec
    smoothWinSize = int(round(0.025/rasterDeltaT))
    ######MARKERSIZE = 2

    selectedTrials = np.ones(behavData['nTrials'],dtype=bool)
    selectedTrials[onecell.trialsToExclude] = False

    spikeTimesFromEventOnset = ephysData['spikeTimesFromEventOnset']
    indexLimitsEachTrial = ephysData['indexLimitsEachTrial']
    timeRange = ephysData['timeRange']
    # -- Define x-label --
    try:
        xLabelStr = ephysData['xLabelStr']
    except:
        xLabelStr = 'Time'
        print "This file does not have 'xLabelStr'"
    timeVec = np.arange(timeRange[0],timeRange[-1]+rasterDeltaT,rasterDeltaT)

    if groupedBy=='4cond':
        (trialsEachCond,condInfo) = sessionanalysis.trials_by_condition(behavData,1,
                                                                        outcome='correct',
                                                                        selected=selectedTrials)
    elif 0:
        (trialsEachCond,condInfo) = sessionanalysis.trials_by_condition(behavData,1,
                                                                        outcome='correctPerBlock',
                                                                        selected=selectedTrials)
    elif groupedBy=='correctness':
        (trialsEachCond,condInfo) = sessionanalysis.trials_by_condition(behavData,1,
                                                                        outcome='correctness',
                                                                        selected=selectedTrials)
    colorEachCond = condInfo['colorEachCond']

    (PSTH,binsStartTime,spikeRasterMat) = \
        spikesanalysis.calculate_psth_per_condition(spikeTimesFromEventOnset,
                                                    indexLimitsEachTrial,
                                                    timeVec,windowSizePSTH,
                                                    trialsEachCond)
    
    # -- Plot spikes raster --
    plt.clf()
    plt.gcf().set_facecolor('w')
    ax1 = plt.axes([0.1,0.48,0.8,0.45])
    boxAx = ax1.get_position()
    (pRaster,hcond) = extraplots.rasterplot(1e3*spikeTimesFromEventOnset,indexLimitsEachTrial,
                                            1e3*timeRange,trialsEachCond,colorEachCond)
    plt.setp(pRaster,ms=RASTER_MARKERSIZE)
    plt.ylim([0,np.sum([len(x) for x in trialsEachCond])-1])
    plt.title(str(onecell),fontsize=12)
    plt.setp(ax1.get_xticklabels(),visible=False)

    # -- Plot reaction time --
    ax3 = plt.axes([boxAx.xmin,0.1+0.35+0.0075,boxAx.width,0.015],sharex=ax1)
    ax3.set_frame_on(False)
    ax3.set_yticks([])
    ax3.set_xticks([])
    # These times are measured from sound offset
    reactionTime = behavData.centerOutTime-behavData.targetOnsetTime-behavData['TargetDuration']
    sideInTime = behavData.sideInTime-behavData.targetOnsetTime-behavData['TargetDuration']
    sideInFromCoutTime = behavData.sideInTime-behavData.centerOutTime
    if lockedTo=='SoundOn':
        sortedTrials = np.concatenate(trialsEachCond);
        sortedRT = reactionTime[sortedTrials] + behavData['TargetDuration'][sortedTrials]
        sortedSideIn = sideInTime[sortedTrials] + behavData['TargetDuration'][sortedTrials]
        plt.hold(True)
        #plt.plot(1e3*sortedRT,range(len(sortedRT)),'o',mfc='none',mec='y',ms=3)
        #plt.plot(1e3*sortedSideIn,range(len(sortedSideIn)),'o',mfc='none',mec='g',ms=3)
        plt.plot(1e3*sortedRT,range(len(sortedRT)),'.y',ms=1)
        plt.plot(1e3*sortedSideIn,range(len(sortedSideIn)),'.g',ms=1)
        plt.hold(False)
    if lockedTo==2: # Locked to Cout
        sortedTrials = np.concatenate(trialsEachCond);
        sortedSideIn = behavData.sideInTime[sortedTrials]-behavData.centerOutTime[sortedTrials]
        plt.hold(True)
        plt.plot(1e3*sortedSideIn,range(len(sortedSideIn)),'o',mfc='none',mec='g',ms=3)
        plt.hold(False)

    #plt.clf()
    #plt.axes([0.12,0.15,0.8,0.7])
    ax2 = plt.axes([boxAx.xmin,0.1,boxAx.width,0.35],sharex=ax1)
    extraplots.plot_psth(PSTH,smoothWinSize,binsStartTime,colorEachCond=colorEachCond,linewidth=2)
    plt.xlim(1e3*timeRange+[0.05,-0.05])
    plt.xlabel('%s (ms)'%xLabelStr)
    plt.ylabel('Firing rate (spk/s)')

    #break



