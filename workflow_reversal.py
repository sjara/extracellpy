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
from extracellpy import frequencytuning
reload(frequencytuning)
from extracellpy import celldatabase
import os
import sys
import shutil
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

RASTER_MARKERSIZE = 2

def behavior_summary(animalName,sessionList,trialslim=[]):
    '''
    trialslim: array to set xlim() of dynamics' plot
    '''
    #dateRange
    nSessions = len(sessionList)
    gs = gridspec.GridSpec(nSessions, 3)
    gs.update(hspace=0.5,wspace=0.4)
    plt.clf()
    for inds,thisSession in enumerate(sessionList):
        try:
            behavData = sessionanalysis.load_behavior_session(animalName,thisSession)
        except:
            print thisSession+' does not exist'
            continue
        print 'Loaded %s %s'%(animalName,thisSession)
        plt.subplot(gs[3*inds])
        behavData.plot_summary(fontsize=10)
        plt.subplot(gs[3*inds+1:3*inds+3])
        behavData.plot_dynamics(winsize=40,fontsize=10)
        if trialslim:
            plt.xlim(trialslim)
        plt.draw()
        plt.show()
    
def save_freqtuning_plots(animalName,copytogether=True):
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    dataModule = 'alltuning_%s'%(animalName)
    allMU = __import__(dataModule)
    reload(allMU)
    
    outputPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
    outputDir = os.path.join(outputPath,'freqtuning_rasters')
    if not os.path.exists(outputDir):
        print 'Creating output directory: %s'%(outputDir)
        os.makedirs(outputDir)
    print 'Saving figures to %s'%(outputDir)
    
    for indmu,onemu in enumerate(allMU.muDB):
        titleString = '%s [%s] T%d'%(onemu.animalName,onemu.ephysSession,onemu.tetrode)
        (spikeTimesFromEventOnset,trialIndexForEachSpike,
            freqEachTrial,timeRange) = frequencytuning.estimate_frequency_tuning(onemu)
        print titleString,
        frequencytuning.plot_frequency_tuning_raster(spikeTimesFromEventOnset,
                                                     trialIndexForEachSpike,
                                                     freqEachTrial,timeRange)
        plt.title(titleString,fontsize=12)
        plt.draw()
        plt.show()

        tetrodeOutputDir = os.path.join(outputDir,'T%02d'%onemu.tetrode)
        if not os.path.exists(tetrodeOutputDir):
            print 'Creating output directory: %s'%(tetrodeOutputDir)
            os.makedirs(tetrodeOutputDir)
        
        plt.gcf().set_size_inches((8,6))
        figFormat = 'png' #'png' #'pdf' #'svg'
        figName = 'freqtuning_%s_T%02d.%s'%(onemu.behavSession,onemu.tetrode,figFormat)
        plt.savefig(os.path.join(tetrodeOutputDir,figName),format=figFormat)
        if copytogether:
            togetherOutputDir = os.path.join(outputDir,'all')
            if not os.path.exists(togetherOutputDir):
                print 'Creating output directory (for all): %s'%(togetherOutputDir)
                os.makedirs(togetherOutputDir)
            shutil.copyfile(os.path.join(tetrodeOutputDir,figName),
                            os.path.join(togetherOutputDir,figName))
            print '... figure saved and copied.'
        else:
            print '... figure saved.'
        

    
def save_raster_data(animalName,lockedTo=None):
    '''
    Save spike times aligned to a task event.
    Data directories are defined in .../extracellpy/settings.py
    It will read the list of cells from a file called allcells_XXXX.py
    (where XXXX is the animal name)

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
            #timeRange = np.array([-0.5,2]) #sec
            allcells.cellDB.save_locked_spikes(outputDir,timeRange=timeRange,lockTo=3)

def save_raster_plots(animalName,groupedBy='4cond'):
    '''
    Parameters
    ----------
    groupedBy: string. Either '4cond','outcome'
    '''
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    dataModule = 'allcells_%s'%(animalName)
    allcells = __import__(dataModule)
    reload(allcells)

    #groupedBy = '4cond'
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
                print 'Behavior already loaded (%s)'%cellStr
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
    elif groupedBy=='outcome':
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


def save_zscores_response(animalName,lockedTo='SoundOn'):
    '''
    Calculate z-scores for firing rates between each time bin and a pre-stimulus bin.
    Save spike times aligned to a task event.
    Data directories are defined in .../extracellpy/settings.py

    Parameters
    ----------
    animalName: string. Name of animal.For example 'saja000'
    lockedTo  : string. Currently only 'SoundOn' is supported.
    '''
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    dataModule = 'allcells_%s'%(animalName)
    allcells = __import__(dataModule)
    reload(allcells)
    dataPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
    dataDir = os.path.join(dataPath,'lockedTo%s'%(lockedTo))
    outputDir = os.path.join(dataPath,'zscores_response_%s'%lockedTo)

    baseRange = np.r_[0,0.025]-0.100
    rangeLength = np.diff(baseRange)
    rangeStart = np.arange(-8,24)*rangeLength # HARDCODED

    if not os.path.exists(outputDir):
        print 'Creating output directory: %s'%(outputDir)
        os.makedirs(outputDir)

    nCells = len(allcells.cellDB)
    prevSession = ''
    for indcell,onecell in enumerate(allcells.cellDB):
        cellStr = str(onecell).replace(' ','_')
        fileName = os.path.join(dataDir,cellStr+'_'+lockedTo+'.npz')
        ephysData = np.load(fileName)
        spikeTimesFromEventOnset = ephysData['spikeTimesFromEventOnset']
        indexLimitsEachTrial = ephysData['indexLimitsEachTrial']

        # -- Load behavior --
        if ephysData['behavSession']==prevSession:
            print 'Behavior already loaded (%s)'%cellStr,
        else:
            print 'Loading %s [%s] - %s'%(ephysData['animalName'],ephysData['behavSession'],cellStr)
            behavData = sessionanalysis.load_behavior_session(ephysData['animalName'],
                                                              ephysData['behavSession'])
            prevSession = ephysData['behavSession']
            (trialsEachCond,condInfo) = sessionanalysis.trials_by_condition(behavData,1,outcome='correct')
        nCond = len(trialsEachCond)
        zStats = np.empty((len(rangeStart),nCond))
        for indCond in range(nCond):
            (zStats[:,indCond],pValues) = spikesanalysis.evaluate_responsiveness(spikeTimesFromEventOnset,
                                                   indexLimitsEachTrial[:,trialsEachCond[indCond]],
                                                   baseRange,rangeStart)
        outputFileName = os.path.join(outputDir,'zscore_resp_'+cellStr+'_'+lockedTo+'.npz')
        np.savez(outputFileName,zStats=zStats, rangeStart=rangeStart,baseRange=baseRange,
                 trialsEachCondLabels=condInfo['trialsEachCondLabels'],cellInfo=onecell,
                 colorEachCond=condInfo['colorEachCond'])
        #zStatsEachCell[:,:,indcell] = zStats


def plot_zscores_response(animalName,lockedTo='SoundOn'):
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    dataModule = 'allcells_%s'%(animalName)
    allcells = __import__(dataModule)
    reload(allcells)
    dataPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
    dataDir = os.path.join(dataPath,'zscores_response_%s'%lockedTo)
    #outputDir = os.path.join(dataPath,'zplots_response_%s'%lockedTo)
    
    plt.clf()
    for indcell,onecell in enumerate(allcells.cellDB):
        #for ind in range(nCells):
        cellStr = str(onecell).replace(' ','_')
        fileName = os.path.join(dataDir,'zscore_'+cellStr+'_'+lockedTo+'.npz')
        zScoreData = np.load(fileName)
        zStats = zScoreData['zStats']
        rangeStart = zScoreData['rangeStart']
        #cellInfo = zScoreData['cellInfo']
        colorEachCond = zScoreData['colorEachCond']
        nCond = zStats.shape[1]

        for indc in range(nCond):
            thisColor = colorEachCond[indc]
            hp = plt.plot(1e3*(rangeStart+np.diff(rangeStart[:2])/2),zStats[:,indc],'o-')
            plt.setp(hp,color=thisColor,mec=thisColor,mfc=thisColor)
            plt.hold(True)
        h1 = plt.axhline(3)
        h2 = plt.axhline(-3)
        plt.setp((h1,h2),linestyle=':',color='0.5')
        plt.hold(False)

        plt.xlabel('Time from XXX (ms)')
        plt.ylabel('z-score')

        plt.draw()
        plt.show()
        plt.waitforbuttonpress()


def load_cells_database(animalsNames):
    '''
    Load and append cell definitions from multiple animals.
    '''
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    cellDB = celldatabase.CellDatabase()
    for animalName in animalsNames:
        dataModule = 'allcells_%s'%(animalName)
        allcells = __import__(dataModule)
        reload(allcells)
        cellDB.extend(allcells.cellDB)
    return cellDB
    

def OLD_save_summary_responsiveness(animalName):
    '''
    THIS WORKS ONLY FOR ONE ANIMALS AT THE TIME
    Evaluate which cells show a response to each of the stimuli (given z-scores).
    Results are saved in sajaXXX_processed/zscores_response_SoundOn/zscoresAll_sajaXXX...
    '''
    zThreshold = 3
    rangeToEvaluate = [0,0.2] # sec from event onset
    lockedTo = 'SoundOn'
    
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    dataModule = 'allcells_%s'%(animalName)
    allcells = __import__(dataModule)
    reload(allcells)
    dataPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
    dataDir = os.path.join(dataPath,'zscores_response_%s'%lockedTo)
    outputFileName = os.path.join(dataDir,'zscore_resp_All_%s_%s.npz'%(animalName,lockedTo))

    # -- Load first cell to get dimensions of data --
    cellStr = str(allcells.cellDB[0]).replace(' ','_')
    zScoreData = np.load(os.path.join(dataDir,'zscore_'+cellStr+'_'+lockedTo+'.npz'))
    nCells = len(allcells.cellDB)
    nCond = zScoreData['zStats'].shape[1]
    nBins = len(zScoreData['rangeStart'])
    zStatsEachCell = np.empty((nBins,nCond,nCells))
    strEachCell = []

    print 'Loading all zscores... ',
    for indcell,onecell in enumerate(allcells.cellDB):
        #for ind in range(nCells):
        cellStr = str(onecell).replace(' ','_')
        fileName = os.path.join(dataDir,'zscore_'+cellStr+'_'+lockedTo+'.npz')
        zScoreData = np.load(fileName)
        zStatsEachCell[:,:,indcell] = zScoreData['zStats']
        strEachCell.append(cellStr)
    print 'done!'
    rangeStart = zScoreData['rangeStart']
    colorEachCond = zScoreData['colorEachCond']

    rangesOfInterest = (rangeStart>=rangeToEvaluate[0])&(rangeStart<rangeToEvaluate[-1])

    responsiveLowFreq = (np.sum(zStatsEachCell[rangesOfInterest,0,:]>zThreshold,axis=0)>0) | \
                        (np.sum(zStatsEachCell[rangesOfInterest,0,:]<-zThreshold,axis=0)>0)
    responsiveMidFreq = (np.sum(zStatsEachCell[rangesOfInterest,1,:]>zThreshold,axis=0)>0) | \
                        (np.sum(zStatsEachCell[rangesOfInterest,1,:]<-zThreshold,axis=0)>0) | \
                        (np.sum(zStatsEachCell[rangesOfInterest,2,:]>zThreshold,axis=0)>0) | \
                        (np.sum(zStatsEachCell[rangesOfInterest,2,:]<-zThreshold,axis=0)>0)
    responsiveHighFreq = (np.sum(zStatsEachCell[rangesOfInterest,3,:]>zThreshold,axis=0)>0) | \
                         (np.sum(zStatsEachCell[rangesOfInterest,3,:]<-zThreshold,axis=0)>0)

    print 'Saving summary to %s'%outputFileName
    np.savez(outputFileName,zStatsEachCell=zStatsEachCell,zThreshold=zThreshold,
             rangesOfInterest=rangesOfInterest,strEachCell=strEachCell,
             rangeStart=zScoreData['rangeStart'],baseRange=zScoreData['baseRange'],
             trialsEachCondLabels=zScoreData['trialsEachCondLabels'],
             responsiveLowFreq=responsiveLowFreq,responsiveMidFreq=responsiveMidFreq,
             responsiveHighFreq=responsiveHighFreq)

def save_summary_responsiveness(animalsNames):
    '''
    Evaluate which cells show a response to each of the stimuli (given z-scores).
    Results are saved in settings.PROCESSED_REVERSAL_PATH
      where %s is replaced by 'all'.
    '''
    zThreshold = 3
    rangeToEvaluate = [0,0.2] # sec from event onset
    lockedTo = 'SoundOn'
    
    cellDB = load_cells_database(animalsNames)

    # -- Load first cell to get dimensions of data --
    onecell = cellDB[0]
    cellStr = str(onecell).replace(' ','_')
    dataPath = settings.PROCESSED_REVERSAL_PATH%(onecell.animalName)
    dataDir = os.path.join(dataPath,'zscores_response_%s'%lockedTo)
    zScoreData = np.load(os.path.join(dataDir,'zscore_resp_'+cellStr+'_'+lockedTo+'.npz'))
    nCells = len(cellDB)
    nCond = zScoreData['zStats'].shape[1]
    nBins = len(zScoreData['rangeStart'])
    zStatsEachCell = np.empty((nBins,nCond,nCells))
    strEachCell = []

    print 'Loading all zscores... ',
    for indcell,onecell in enumerate(cellDB):
        dataPath = settings.PROCESSED_REVERSAL_PATH%(onecell.animalName)
        dataDir = os.path.join(dataPath,'zscores_response_%s'%lockedTo)
        cellStr = str(onecell).replace(' ','_')
        fileName = os.path.join(dataDir,'zscore_resp_'+cellStr+'_'+lockedTo+'.npz')
        zScoreData = np.load(fileName)
        zStatsEachCell[:,:,indcell] = zScoreData['zStats']
        strEachCell.append(cellStr)
    print 'done!'
    rangeStart = zScoreData['rangeStart']
    colorEachCond = zScoreData['colorEachCond']

    rangesOfInterest = (rangeStart>=rangeToEvaluate[0])&(rangeStart<rangeToEvaluate[-1])

    responsiveLowFreq = (np.sum(zStatsEachCell[rangesOfInterest,0,:]>zThreshold,axis=0)>0) | \
                        (np.sum(zStatsEachCell[rangesOfInterest,0,:]<-zThreshold,axis=0)>0)
    responsiveMidFreq = (np.sum(zStatsEachCell[rangesOfInterest,1,:]>zThreshold,axis=0)>0) | \
                        (np.sum(zStatsEachCell[rangesOfInterest,1,:]<-zThreshold,axis=0)>0) | \
                        (np.sum(zStatsEachCell[rangesOfInterest,2,:]>zThreshold,axis=0)>0) | \
                        (np.sum(zStatsEachCell[rangesOfInterest,2,:]<-zThreshold,axis=0)>0)
    responsiveHighFreq = (np.sum(zStatsEachCell[rangesOfInterest,3,:]>zThreshold,axis=0)>0) | \
                         (np.sum(zStatsEachCell[rangesOfInterest,3,:]<-zThreshold,axis=0)>0)


    strAllAnimals = '-'.join(animalsNames)
    outputDir = settings.PROCESSED_REVERSAL_PATH%('all')
    if not os.path.exists(outputDir):
        print 'Creating output directory: %s'%(outputDir)
        os.makedirs(outputDir)
    outputFileName = os.path.join(outputDir,'summary_resp_%s_%s.npz'%(strAllAnimals,lockedTo))
    print 'Saving summary to %s'%outputFileName
    np.savez(outputFileName,zStatsEachCell=zStatsEachCell,zThreshold=zThreshold,
             rangesOfInterest=rangesOfInterest,strEachCell=strEachCell,
             rangeStart=zScoreData['rangeStart'],baseRange=zScoreData['baseRange'],
             trialsEachCondLabels=zScoreData['trialsEachCondLabels'],
             responsiveLowFreq=responsiveLowFreq,responsiveMidFreq=responsiveMidFreq,
             responsiveHighFreq=responsiveHighFreq)

def print_summary_responsiveness(animalName):
    '''
    Print summary of responses.
    '''
    lockedTo = 'SoundOn'
    
    dataPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
    dataDir = os.path.join(dataPath,'zscores_response_%s'%lockedTo)
    outputFileName = os.path.join(dataDir,'zscore_resp_All_%s_%s.npz'%(animalName,lockedTo))

    rdata = np.load(outputFileName)
    print 'Cells that respond to each freq (|z|>3) out of %d:'%len(rdata['responsiveLowFreq'])
    print 'Low freq: %0.1f%%'%(100*np.mean(rdata['responsiveLowFreq']))
    print 'Mid freq: %0.1f%%'%(100*np.mean(rdata['responsiveMidFreq']))
    print 'High freq: %0.1f%%'%(100*np.mean(rdata['responsiveHighFreq']))

    '''
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    dataModule = 'allcells_%s'%(animalName)
    allcells = __import__(dataModule)
    reload(allcells)
    '''

def save_zscores_modulation(animalName,lockedTo='SoundOn'):
    '''
    Calculate z-scores between midFreq-Right and midFreq-Left
    Data directories are defined in .../extracellpy/settings.py

    Parameters
    ----------
    animalName: string. Name of animal.For example 'saja000'
    lockedTo  : string. Currently only 'SoundOn' is supported.
    '''
    
    MIN_TRIALS_PER_BLOCK = 75
    responseRange = [0.010,0.150]

    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    dataModule = 'allcells_%s'%(animalName)
    allcells = __import__(dataModule)
    reload(allcells)
    dataPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
    dataDir = os.path.join(dataPath,'lockedTo%s'%(lockedTo))
    outputDir = os.path.join(dataPath,'zscores_modulation_%s'%lockedTo)
    
    if not os.path.exists(outputDir):
        print 'Creating output directory: %s'%(outputDir)
        os.makedirs(outputDir)

    typeEachSwitchLabels = {'LowBoundToHighBound':0,'HighBoundToLowBound':1}
    respDuration = np.diff(responseRange)

    indSwitch = 0

    prevSession = ''
    for indcell,onecell in enumerate(allcells.cellDB):
        cellStr = str(onecell).replace(' ','_')
        fileName = os.path.join(dataDir,cellStr+'_'+lockedTo+'.npz')
        ephysData = np.load(fileName)

        # -- Load behavior --
        if ephysData['behavSession']==prevSession:
            print 'Behavior already loaded'
        else:
            print 'Loading %s [%s] - %s'%(ephysData['animalName'],ephysData['behavSession'],cellStr)
            behavData = sessionanalysis.load_behavior_session(ephysData['animalName'],
                                                              ephysData['behavSession'])

            trialsMidFreq = (behavData['TargetFreq']==behavData['FreqMid'][-1])
            prevSession = ephysData['behavSession']

        # -- Find modulation for all blocks merged --
        eachCondLabel = ['LowBoundBlock','HighBoundBlock']
        validMidFreq = trialsMidFreq & ~behavData.early
        validEachCond = np.c_[behavData.lowFreqs,behavData.highFreqs]
        validMidFreqEachCond = validEachCond & validMidFreq[:,np.newaxis]
        trialsToCompare = validMidFreqEachCond
        (meanRespEachCond,pValueMod) = spikesanalysis.evaluate_modulation(ephysData['spikeTimesFromEventOnset'],
                                                                          ephysData['indexLimitsEachTrial'],
                                                                          responseRange,list(trialsToCompare.T))

        # Because it is per block, I don't know how many yet
        meanRespEachSwitch = np.empty((0,2)) # in spk/sec  array with [lowBound,highBound]
        pValueEachSwitch =  np.empty(0)
        cellIDeachSwitch =  np.empty(0,dtype=int)
        typeEachSwitch = np.empty(0,dtype=int)

        #########################################################
        # WARNING!!! check that it works for old misaligned data
        #########################################################

        behavData.find_trials_each_block()
        validTrialsEachBlock = behavData.trialsEachBlock & ~behavData.early[:,np.newaxis]
        nBlocks = validTrialsEachBlock.shape[1]
        nTrialsEachBlock = validTrialsEachBlock.sum(axis=0)

        for indb in range(nBlocks):
            if nTrialsEachBlock[indb]>=MIN_TRIALS_PER_BLOCK:
                if (indb+1)<nBlocks and nTrialsEachBlock[indb+1]>=MIN_TRIALS_PER_BLOCK:
                    if(behavData.lowFreqs[behavData.firstTrialEachBlock[indb]]):
                        thisType = typeEachSwitchLabels['LowBoundToHighBound']
                        trialsToCompare = validTrialsEachBlock[:,[indb,indb+1]] & trialsMidFreq[:,np.newaxis]
                    elif (behavData.highFreqs[behavData.firstTrialEachBlock[indb]]):
                        thisType = typeEachSwitchLabels['HighBoundToLowBound']
                        trialsToCompare = validTrialsEachBlock[:,[indb+1,indb]] & trialsMidFreq[:,np.newaxis]
                    else:
                        thisType = -1
                        trialsToCompare = []
                    (meanSpikes,pValue) = spikesanalysis.evaluate_modulation(ephysData['spikeTimesFromEventOnset'],
                                                                             ephysData['indexLimitsEachTrial'],
                                                                             responseRange,list(trialsToCompare.T))
                    meanRespEachSwitch = np.vstack((meanRespEachSwitch,meanSpikes/respDuration))
                    pValueEachSwitch = np.hstack((pValueEachSwitch,pValue))
                    typeEachSwitch = np.hstack((typeEachSwitch,thisType))
                else:
                    break

            else:
                continue
        # --- Save data for this cell ---
        outputFileName = os.path.join(outputDir,'zscore_mod_'+cellStr+'_'+lockedTo+'.npz')
        np.savez(outputFileName,responseRange=responseRange,
                 meanRespEachSwitch=meanRespEachSwitch,
                 pValueEachSwitch=pValueEachSwitch, typeEachSwitch=typeEachSwitch,
                 typeEachSwitchLabels=typeEachSwitchLabels,
                 eachCondLabel=eachCondLabel, meanRespEachCond=meanRespEachCond,
                 pValueMod=pValueMod)


def OLD_save_summary_modulation(animalName,lockedTo='SoundOn'):
    '''
    THIS WORKS FOR A SINGLE ANIMAL
    Create array with z-scores from all cells.
    Data directories are defined in .../extracellpy/settings.py

    Parameters
    ----------
    animalName: string. Name of animal.For example 'saja000'
    lockedTo  : string. Currently only 'SoundOn' is supported.
    '''
 
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    dataModule = 'allcells_%s'%(animalName)
    allcells = __import__(dataModule)
    reload(allcells)
    dataPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
    dataDir = os.path.join(dataPath,'zscores_modulation_%s'%lockedTo)

    # Because it is per block, I don't know how many yet
    meanRespEachSwitch = np.empty((0,2)) # in spk/sec
    pValueEachSwitch =  np.empty(0)
    cellIDeachSwitch =  np.empty(0,dtype=int)
    typeEachSwitch = np.empty(0,dtype=int)
    strEachCell = []

    nCells = len(allcells.cellDB)
    meanRespEachCond = np.empty((nCells,2),dtype=float)
    pValueMod = np.empty(nCells,dtype=float)
    consistentMod =  np.zeros(nCells,dtype=bool)
    
    print 'Loading all modulation zscores... '
    for indcell,onecell in enumerate(allcells.cellDB):
        cellStr = str(onecell).replace(' ','_')
        fileName = os.path.join(dataDir,'zscore_mod_'+cellStr+'_'+lockedTo+'.npz')
        zScoreData = np.load(fileName)
        nSwitchesThisCell = len(zScoreData['pValueEachSwitch'])

        meanRespEachSwitch = np.vstack((meanRespEachSwitch,zScoreData['meanRespEachSwitch']))
        pValueEachSwitch = np.hstack((pValueEachSwitch,zScoreData['pValueEachSwitch']))
        typeEachSwitch = np.hstack((typeEachSwitch,zScoreData['typeEachSwitch']))
        cellIDeachSwitch = np.hstack((cellIDeachSwitch,nSwitchesThisCell*[indcell]))
        strEachCell.append(cellStr)
        meanRespEachCond[indcell] = zScoreData['meanRespEachCond']
        pValueMod[indcell] = zScoreData['pValueMod']

        # -- Evaluate if consistent change across switches --
        modDirEachCell = np.diff(zScoreData['meanRespEachSwitch'],axis=1)>0
        sumModDir = np.sum(modDirEachCell)
        if (sumModDir==0) | (sumModDir==nSwitchesThisCell):
            consistentMod[indcell]=True

    outputFileName = os.path.join(dataDir,'zscore_mod_All_%s_%s.npz'%(animalName,lockedTo))
    print 'Saving summary to %s'%outputFileName
    np.savez(outputFileName,strEachCell=strEachCell,meanRespEachSwitch=meanRespEachSwitch,
             pValueEachSwitch=pValueEachSwitch,typeEachSwitch=typeEachSwitch,
             typeEachSwitchLabels=zScoreData['typeEachSwitchLabels'],
             cellIDeachSwitch=cellIDeachSwitch,
             eachCondLabel=zScoreData['eachCondLabel'],
             meanRespEachCond=meanRespEachCond,
             pValueMod=pValueMod, consistentMod=consistentMod)

def save_summary_modulation(animalsNames,lockedTo='SoundOn'):
    '''
    Create array with z-scores from all cells.
    Data directories are defined in .../extracellpy/settings.py
    Results are saved in settings.PROCESSED_REVERSAL_PATH
      where %s is replaced by 'all'.

    Parameters
    ----------
    animalsNames: string. Name of animal.For example 'saja000'
    lockedTo  : string. Currently only 'SoundOn' is supported.
    '''

    cellDB = load_cells_database(animalsNames)
    
    # Because it is per block, I don't know how many yet
    meanRespEachSwitch = np.empty((0,2)) # in spk/sec
    pValueEachSwitch =  np.empty(0)
    cellIDeachSwitch =  np.empty(0,dtype=int)
    typeEachSwitch = np.empty(0,dtype=int)
    strEachCell = []

    nCells = len(cellDB)
    meanRespEachCond = np.empty((nCells,2),dtype=float)
    pValueMod = np.empty(nCells,dtype=float)
    consistentMod =  np.zeros(nCells,dtype=bool)
    
    print 'Loading all modulation zscores... '
    for indcell,onecell in enumerate(cellDB):
        dataPath = settings.PROCESSED_REVERSAL_PATH%(onecell.animalName)
        dataDir = os.path.join(dataPath,'zscores_modulation_%s'%lockedTo)
        cellStr = str(onecell).replace(' ','_')
        fileName = os.path.join(dataDir,'zscore_mod_'+cellStr+'_'+lockedTo+'.npz')
        zScoreData = np.load(fileName)
        nSwitchesThisCell = len(zScoreData['pValueEachSwitch'])

        meanRespEachSwitch = np.vstack((meanRespEachSwitch,zScoreData['meanRespEachSwitch']))
        pValueEachSwitch = np.hstack((pValueEachSwitch,zScoreData['pValueEachSwitch']))
        typeEachSwitch = np.hstack((typeEachSwitch,zScoreData['typeEachSwitch']))
        cellIDeachSwitch = np.hstack((cellIDeachSwitch,nSwitchesThisCell*[indcell]))
        strEachCell.append(cellStr)
        meanRespEachCond[indcell] = zScoreData['meanRespEachCond']
        pValueMod[indcell] = zScoreData['pValueMod']

        # -- Evaluate if consistent change across switches --
        modDirEachCell = np.diff(zScoreData['meanRespEachSwitch'],axis=1)>0
        sumModDir = np.sum(modDirEachCell)
        if (sumModDir==0) | (sumModDir==nSwitchesThisCell):
            consistentMod[indcell]=True

    strAllAnimals = '-'.join(animalsNames)
    outputDir = settings.PROCESSED_REVERSAL_PATH%('all')
    if not os.path.exists(outputDir):
        print 'Creating output directory: %s'%(outputDir)
        os.makedirs(outputDir)
    outputFileName = os.path.join(outputDir,'summary_mod_%s_%s.npz'%(strAllAnimals,lockedTo))
    print 'Saving summary to %s'%outputFileName
    np.savez(outputFileName,strEachCell=strEachCell,meanRespEachSwitch=meanRespEachSwitch,
             pValueEachSwitch=pValueEachSwitch,typeEachSwitch=typeEachSwitch,
             typeEachSwitchLabels=zScoreData['typeEachSwitchLabels'],
             cellIDeachSwitch=cellIDeachSwitch,
             eachCondLabel=zScoreData['eachCondLabel'],
             meanRespEachCond=meanRespEachCond,
             pValueMod=pValueMod, consistentMod=consistentMod)
    
def print_summary_modulation(animalName):
    '''
    Print summary of modulation.
    '''
    lockedTo = 'SoundOn'
    
    dataPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
    dataDir = os.path.join(dataPath,'zscores_modulation_%s'%lockedTo)
    outputFileName = os.path.join(dataDir,'zscore_mod_All_%s_%s.npz'%(animalName,lockedTo))
    mdata = np.load(outputFileName)

    dataDir = os.path.join(dataPath,'zscores_response_%s'%lockedTo)
    outputFileName = os.path.join(dataDir,'zscore_resp_All_%s_%s.npz'%(animalName,lockedTo))
    rdata = np.load(outputFileName)
    
    print 'Cells with different sound response depending on action out of %d:'%len(mdata['pValueMod'])
    print 'From all blocks merged (p<0.05): %0.1f%%'%( 100*np.mean(mdata['pValueMod']<0.05) )
    print 'Consistent across switches (signif and not-signif): %0.1f%%'%( 100*np.mean(mdata['consistentMod']) )
    consistSignif = (mdata['pValueMod']<0.05) & mdata['consistentMod']
    print 'Consistent and significant on merged blocks (p<0.05): %0.1f%%'%( 100*np.mean(consistSignif) )
    print ''
    print 'Percentage of significantly modulated switches (n=%d): %0.1f%%'%(len(mdata['pValueEachSwitch']),
                                                                            100*np.mean(mdata['pValueEachSwitch']))
    print ''
    print 'Responsive to mid freq: %0.1f%% (%d)'%( 100*np.mean(rdata['responsiveMidFreq']),
                                                   sum(rdata['responsiveMidFreq']) )
    respConsistSignif = (mdata['pValueMod']<0.05) & mdata['consistentMod'] & rdata['responsiveMidFreq']
    responsiveTotal = np.sum(rdata['responsiveMidFreq'])
    print 'Responsive, consistent and significant (p<0.05): %0.1f%% (%d)'%( 100*np.sum(respConsistSignif)/responsiveTotal,
                                                                            np.sum(respConsistSignif) )
    print ''
    #signifCells = mdata['pValueMod']<0.05
    highBoundStrongerThanLowBound = np.diff(mdata['meanRespEachCond'],axis=1)>0
    signifRight = np.sum(~highBoundStrongerThanLowBound[respConsistSignif])
    signifLeft = np.sum(highBoundStrongerThanLowBound[respConsistSignif])
    respRight = np.sum(~highBoundStrongerThanLowBound[rdata['responsiveMidFreq']])
    respLeft = np.sum(highBoundStrongerThanLowBound[rdata['responsiveMidFreq']])
    print 'Right vs Left (responsive): %d/%d   signifConsistResp: %d/%d'%(respRight,respLeft,signifRight,signifLeft)


def load_summary_modulation_one_animal(animalName):
    '''
    Load summary information.
    '''
    lockedTo = 'SoundOn'
    
    dataPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
    dataDir = os.path.join(dataPath,'zscores_modulation_%s'%lockedTo)
    outputFileName = os.path.join(dataDir,'zscore_mod_All_%s_%s.npz'%(animalName,lockedTo))
    mdata = np.load(outputFileName)

    dataDir = os.path.join(dataPath,'zscores_response_%s'%lockedTo)
    outputFileName = os.path.join(dataDir,'zscore_resp_All_%s_%s.npz'%(animalName,lockedTo))
    rdata = np.load(outputFileName)

    return(rdata,mdata)
    

def merge_summary_modulation(animalsName):
    '''
    OBSOLETE (NOT NEEDED ANYMORE): summary of responses and modulation should
                        be calculated for many animals at once. See save_summary_...
    Load summary information.
    '''
    fieldsToMerge = ['typeEachSwitch', 'consistentMod', 'pValueMod',
                     'meanRespEachSwitch', 'pValueEachSwitch', 'meanRespEachCond','strEachCell']
    fieldsToAdd = ['cellIDeachSwitch']
    fieldsToAppend = []
    fieldsToKeep = ['typeEachSwitchLabels','eachCondLabel']
    allFields = fieldsToMerge+fieldsToAdd+fieldsToAppend+fieldsToKeep
    allEmpty = [np.empty(0,dtype=int),np.empty(0,dtype=bool),np.empty(0),
                np.empty((0,2)),np.empty(0),np.empty((0,2)),
                [],[],[],
                np.empty(0)]
    allData = {}
    allData['responsiveMidFreq'] = np.empty(0,dtype=bool)
    nCellsSoFar = 0
    for field,initValue in zip(allFields,allEmpty):
        allData[field] = initValue
    for animalName in animalsName:
        (rdata,mdata) = load_summary_modulation_one_animal(animalName)
        for field in fieldsToMerge:
            allData[field] = np.concatenate((allData[field],mdata[field]))
        for field in fieldsToAdd:
            shiftedValues = mdata[field]+nCellsSoFar
            allData[field] = np.concatenate((allData[field],shiftedValues))
        for field in fieldsToAppend:
            allData[field].append(mdata[field])
        nCellsSoFar += len(mdata['strEachCell'])
        #for field in rfieldsToMerge:
        allData['responsiveMidFreq'] = np.concatenate((allData['responsiveMidFreq'],rdata['responsiveMidFreq']))
    for field in fieldsToKeep:
        allData[field] = mdata[field]
    return allData

def load_summary_modulation(animalsNames,lockedTo='SoundOn'):
    ''' 
    Load summary of modulation, responses and list of cells.
    '''
    cellDB = load_cells_database(animalsNames)
    strAllAnimals = '-'.join(animalsNames)
    dataDir = settings.PROCESSED_REVERSAL_PATH%('all')
    respFileName = os.path.join(dataDir,'summary_resp_%s_%s.npz'%(strAllAnimals,lockedTo))
    modFileName = os.path.join(dataDir,'summary_mod_%s_%s.npz'%(strAllAnimals,lockedTo))
    respData = np.load(respFileName)
    modData = np.load(modFileName)
    return (modData,respData,cellDB)

def plot_summary_modulation_TEMP(animalsNames,lockedTo='SoundOn'):
    '''
    TEMP FUNCTION:
    This version is too restrictive. It is here just for illustration.
    '''
    from scipy import stats
    (mdata,rdata,cellDB) = load_summary_modulation(animalsNames,lockedTo=lockedTo)
    respConsistSignif = (mdata['pValueMod']<0.05) & mdata['consistentMod'] & rdata['responsiveMidFreq']

    dataToPlot = mdata['meanRespEachCond'][rdata['responsiveMidFreq']]
    pValToPlot = mdata['pValueMod'][rdata['responsiveMidFreq']]
    pValToPlot[~mdata['consistentMod'][rdata['responsiveMidFreq']]] = 1

    # --- Plot results ---
    from extracellpy import extraplots
    plt.clf()
    plt.setp(plt.gcf(),facecolor='w')
    nBins = 32 #14
    plt.subplot(1,2,1)
    extraplots.plot_index_histogram(dataToPlot[:,0],
                                    dataToPlot[:,1],
                                    pValue=pValToPlot,nBins=nBins)
    plt.subplot(1,2,2)
    extraplots.plot_scatter(dataToPlot[:,0],
                            dataToPlot[:,1],
                            pValue=pValToPlot)

    (tstatw,pval) = stats.wilcoxon(dataToPlot[:,0],dataToPlot[:,1]) # paired test
    print 'p-value = %0.4f'%pval

    
#def show_raster_figures(animalName,lockedTo,cellIndexList):
def show_raster_figures(cellDB,lockedTo='SoundOn'):
    '''
    Show raster figured (previously saved as PNG) for each cell in cellDB

    Use for example:
    cellIndexList = np.argsort(pValueMod)
    show_raster_figures(cellDB[cellIndexList])
    '''
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    figFilenameFormat = 'raster_%s_%s_T%dc%d.png'

    #for indcell in cellIndexList:
    for onecell in cellDB:
        dataModule = 'allcells_%s'%(onecell.animalName)
        allcells = __import__(dataModule)
        reload(allcells)

        dataPath = settings.PROCESSED_REVERSAL_PATH%(onecell.animalName)
        figDir = os.path.join(dataPath,'rasters_lockedTo%s_4cond'%lockedTo)
        figFilename = figFilenameFormat%(onecell.animalName,onecell.ephysSession,
                                         onecell.tetrode,onecell.cluster)
        img=mpimg.imread(os.path.join(figDir,figFilename))
        plt.imshow(img,aspect='equal',interpolation='bilinear')
        #print pValueMod[indcell]
        plt.axis('off')
        plt.draw()
        plt.show()
        plt.waitforbuttonpress()

