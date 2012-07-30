'''
Module for analyzing frequency-tuning sessions.
'''

from extracellpy import settings
reload(settings) # Force reload
from extracellpy import loadneuralynx
reload(loadneuralynx) # Force reload
from extracellpy import loadbehavior
from extracellpy import spikesanalysis
from extracellpy import extrafuncs
import os, sys
import numpy as np
import matplotlib.pyplot as plt

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

BEHAVIORPATH = settings.BEHAVIOR_PATH
EPHYSPATH = settings.EPHYS_PATH

bitTRIALIND = 3      # TrialIndicator (bitID starting from 0)
bitPHOTOSTIMIND = 4  # PhotoStimIndicator (bitID starting from 0)
bitTARGETIND = 5     # TargetIndicator (bitID starting from 0)

RASTER_MARKERSIZE = 1
FONTSIZE = 14

def plot_frequency_tuning_raster(spikeTimesFromEventOnset,trialIndexForEachSpike,freqEachTrial,timeRange):
    xLims = 1e3*timeRange
    yLims = [0,trialIndexForEachSpike[-1]-1]
    xLabel = 'Time from sound onset (ms)'
    plt.clf()
    plt.gcf().set_facecolor('w')
    ax = plt.axes([0.2,0.2,0.7,0.7])
    pR = plt.plot(1e3*spikeTimesFromEventOnset,trialIndexForEachSpike,'.k')
    pR[0].set_markersize(RASTER_MARKERSIZE)
    plt.xlim(xLims)
    plt.ylim(yLims)
    plt.xlabel(xLabel,fontsize=FONTSIZE)
    #ax.set_yticklabels('')
    plt.ylabel('Trial',fontsize=FONTSIZE)

    # -- Draw background for each frequency --
    #FreqEachTrial = behavData['SoundFreq'][trialsOfInterest]
    #possibleFreqs = np.unique(behavData['SoundFreq'])
    possibleFreqs = np.unique(freqEachTrial)
    lastTrialEachFreq = np.flatnonzero(np.diff(freqEachTrial))
    firstTrialEachFreq = np.hstack((0,lastTrialEachFreq+1))
    lastTrialEachFreq = np.hstack((lastTrialEachFreq, len(freqEachTrial)))
    # NOTE: indexes are from 0 to N+1 (python style)
    nFreq = len(firstTrialEachFreq)
    plt.hold(True)
    possibleColors = ['0.95','0.9']
    for indf in range(nFreq):
        xpos = np.hstack((xLims,xLims[::-1]))
        ypos = np.hstack((np.tile(firstTrialEachFreq[indf]-0.5,2),np.tile(lastTrialEachFreq[indf]+0.5,2)))
        thisColor = possibleColors[np.mod(indf,2)]
        bg = plt.fill(xpos,ypos,color=thisColor,ec='none')
        ht = plt.text(0.99*xLims[1],firstTrialEachFreq[indf],'%d'%possibleFreqs[indf],color='b',
                  backgroundcolor='w',horizontalalignment='right')
    plt.hold(False)
    plt.ylim(yLims)

    plt.draw()
    plt.show()



def align_data(onemu):
    rasterDeltaT = 0.1e-3            # sec
    timeRange = np.array([-0.2,0.6])
    animalName = onemu.animalName
    behavSession = onemu.behavSession
    ephysSession = onemu.ephysSession
    tetrode = onemu.tetrode
    clusters = onemu.clusters
    
    # -- Load events from Neuralynx --
    dataDir = os.path.join(settings.EPHYS_PATH,'%s/%s/'%(animalName,ephysSession))
    clustersDir = os.path.join(settings.EPHYS_PATH,'%s/%s_kk/'%(animalName,ephysSession))
    eventsFile = os.path.join(dataDir,'Events.nev')
    events = loadneuralynx.DataEvents(eventsFile)
    trialEvents = (events.valueTTL & (1<<bitTRIALIND)) != 0
    trialStartTimeNL = 1e-6*events.timestamps[trialEvents]
    targetEvents = (events.valueTTL & (1<<bitTARGETIND)) != 0
    targetTimeNL = 1e-6*events.timestamps[targetEvents]
    #targetFreqInd = events.valueTTL[trialEvents]>>8  # 

    # -- Load events from behavior --
    behavDataDir = os.path.join(settings.BEHAVIOR_PATH,'%s/'%(animalName))
    behavFileName = 'data_saja_tuningcurve_santiago_%s_%s.h5'%(animalName,behavSession)
    behavFile = os.path.join(behavDataDir,behavFileName)
    behavData = loadbehavior.TuningBehaviorData(behavFile)
    behavData.extract_event_times()
    behavData.align_to_ephys(trialStartTimeNL)

    # FIXME: add line to check that number of trials of ephys & behavior are consistent
    
    # -- Check is alignment ephys/behavior is correct --
    #behavData.check_clock_drift()
    #waitforbuttonpress()

    # -- Remove first empty trial from data --
    behavData['nTrials'] = behavData['nTrials']-1
    nTrials = behavData['nTrials']
    behavData.trialStartTime = behavData.trialStartTime[1:]
    behavData.targetOnsetTime = behavData.targetOnsetTime[1:]
    behavData['SoundFreq'] = behavData['SoundFreq'][1:]
    # -- Remove incomplete trial from ephys data --
    #targetTimeNL = targetTimeNL[:nTrials]
    behavData.trialStartTimeEphys = behavData.trialStartTimeEphys[1:] # Remove first empty trial

    trialsOfInterest = np.argsort(behavData['SoundFreq'])

    eventOfInterest = behavData.targetOnsetTime[trialsOfInterest] - \
                      behavData.trialStartTime[trialsOfInterest] + \
                      behavData.trialStartTimeEphys[trialsOfInterest]
    timeVec = np.arange(timeRange[0],timeRange[-1]+rasterDeltaT,rasterDeltaT)

    freqEachTrial = behavData['SoundFreq'][trialsOfInterest]

    # -- Load spikes --
    tetrodeFile = os.path.join(dataDir,'TT%d.ntt'%tetrode)
    dataTT = loadneuralynx.DataTetrode(tetrodeFile)
    dataTT.timestamps = dataTT.timestamps.astype(np.float64)*1e-6  # in sec

    # -- Load clusters if required --
    #if (clustersEachTetrode is not None) and clustersEachTetrode.has_key(tetrode):
    
    ########## BUG: clustersEachTetrode is not defined anywhere !!! #############
    
    if len(clusters)>0:
        clustersFile = os.path.join(clustersDir,'TT%d.clu.1'%tetrode)
        dataTT.set_clusters(clustersFile)
        spikeInds = extrafuncs.ismember(dataTT.clusters,clustersEachTetrode[tetrode])
        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
          spikesanalysis.eventlocked_spiketimes(dataTT.timestamps[spikeInds],eventOfInterest,timeRange)
    else:
        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
          spikesanalysis.eventlocked_spiketimes(dataTT.timestamps,eventOfInterest,timeRange)
    
    return (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial,freqEachTrial,timeRange)

def find_trials_each_freq(freqEachTrial):
    possibleFreq = np.unique(freqEachTrial)
    trialsEachFreq = []
    for indf,freq in enumerate(possibleFreq):
        trialsEachFreq.append(np.flatnonzero(freqEachTrial==freq))
    return (possibleFreq,trialsEachFreq)

def estimate_frequency_tuning(spikeTimesFromEventOnset,indexLimitsEachTrial,freqEachTrial,responseRange=None):
    if responseRange is None:
        responseRange = np.array([0,0.150])
    else:
        responseRange = np.array(responseRange)
    baselineRange = responseRange-0.200
    (possibleFreq,trialsEachFreq) = find_trials_each_freq(freqEachTrial)
    meanRespEachFreq = np.empty(len(possibleFreq))
    semRespEachFreq = np.empty(len(possibleFreq))
    # -- Response --
    for indf,freq in enumerate(possibleFreq):
        theseTrials = trialsEachFreq[indf]
        nSpikes=spikesanalysis.count_spikes_in_range(spikeTimesFromEventOnset,
                                                     indexLimitsEachTrial[:,theseTrials],
                                                     responseRange)
        spikesPerSec = nSpikes/np.diff(responseRange)
        meanRespEachFreq[indf] = np.mean(spikesPerSec)
        semRespEachFreq[indf] = np.std(spikesPerSec)/np.sqrt(len(theseTrials))
    # -- Baseline --
    nSpikesBaseline=spikesanalysis.count_spikes_in_range(spikeTimesFromEventOnset,
                                                         indexLimitsEachTrial,
                                                         baselineRange)
    spikesPerSecBaseline = nSpikesBaseline/np.diff(baselineRange)
    meanBaseline = np.mean(spikesPerSecBaseline)
    semBaseline = np.mean(spikesPerSecBaseline)
    return (possibleFreq,meanRespEachFreq,semRespEachFreq,meanBaseline,semBaseline)

def plot_frequency_tuning(possibleFreq,meanRespEachFreq,semRespEachFreq=None,
                          meanBaseline=None,semBaseline=None):
    pcolor='k'
    plt.clf()
    plt.plot(np.log10(possibleFreq),meanRespEachFreq,'o-',ms=5,mew=2,mec=pcolor,mfc='w',color=pcolor)
    if meanBaseline is not None:
        plt.axhline(meanBaseline,ls='--',color='k')
    plt.draw()
    plt.show()
    
if __name__ == "__main__":
    '''
    animalName = 'saja125'
    # -- Load list of cells --
    sys.path.append(settings.CELL_LIST_PATH)
    dataModule = 'alltuning_%s'%(animalName)
    allMU = __import__(dataModule)
    reload(allMU)
    muDB = allMU.muDB
    '''
    '''
    animalName   = 'saja125'
    ephysSession = '2012-01-30_14-31-32'
    behavSession = '20120130a'
    tetrodes = [2]
    clustersEachTetrode = None
    '''
    from extracellpy import celldatabase
    reload(celldatabase)

    muDB = celldatabase.MultiUnitDatabase()
    animalName   = 'saja100'
    #ephysSession = '2011-11-27_14-33-15'; behavSession = '20111127a'
    #clustersEachTetrode = {7:[13]}
    #clustersEachTetrode = {7:range(2,14)}
    ephysSession = '2011-12-05_16-35-53'; behavSession = '20111205a'
    clustersEachTetrode = {8:range(2,30)}
    #clustersEachTetrode = {8:[3]}
    
    for tetrode,clusters in sorted(clustersEachTetrode.items()):
        oneCell = celldatabase.MultiUnitInfo(animalName = animalName,
                                             ephysSession = ephysSession,
                                             behavSession = behavSession,
                                             tetrode = tetrode,
                                             clusters = clusters)
        muDB.append(oneCell)

    outputPath = settings.PROCESSED_REVERSAL_PATH%(animalName)
    outputDir = os.path.join(outputPath,'freqtuning')

    CASE = 2
    if CASE==1:

        for indmu,onemu in enumerate(muDB):
            #muStr = '%s_%s_T%dmu'%(onemu.animalName, onemu.ephysSession,onemu.tetrode)
            titleString = '%s [%s] T%d'%(onemu.animalName,onemu.ephysSession,onemu.tetrode)
            (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial,
                freqEachTrial,timeRange) = align_data(onemu)

            print titleString,
            plot_frequency_tuning_raster(spikeTimesFromEventOnset,trialIndexForEachSpike,
                                         freqEachTrial,timeRange)

            import matplotlib.pyplot as plt
            plt.title(titleString,fontsize=FONTSIZE-2)
            plt.draw()
            plt.show()
            #plt.waitforbuttonpress()
        
        '''
        plt.gcf().set_size_inches((8,6))
        figFormat = 'png' #'png' #'pdf' #'svg'
        figName = 'freqtuning_%s_T%d.%s'%(onemu.behavSession,onemu.tetrode,figFormat)
        #plt.savefig(os.path.join(outputDir,figName),format=figFormat)
        plt.savefig(os.path.join('/tmp/',figName),format=figFormat)
        print '... figure saved.'
        '''

    elif CASE==2:
        #rRange=[0,0.025]
        rRange=[0.025,0.150]
        #rRange=[0,0.150]
        onemu = muDB[0]
        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial,
           freqEachTrial,timeRange) = align_data(onemu)
        (possibleFreq,meanRespEachFreq,semRespEachFreq,meanBaseline,semBaseline) = \
          estimate_frequency_tuning(spikeTimesFromEventOnset,
                                    indexLimitsEachTrial,freqEachTrial,
                                    responseRange=rRange)
        plt.figure(1)
        plot_frequency_tuning_raster(spikeTimesFromEventOnset,trialIndexForEachSpike,freqEachTrial,timeRange)
        plt.figure(2)
        plot_frequency_tuning(possibleFreq,meanRespEachFreq,meanBaseline=meanBaseline)
