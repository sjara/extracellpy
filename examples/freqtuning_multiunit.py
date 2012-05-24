''' 
Raster for frequency tuning.

Santiago Jaramillo - 2012-03-08
'''

from extracellpy import settings
from extracellpy import loadneuralynx, spikesanalysis, loadbehavior
from extracellpy import extraplots
from extracellpy import extrafuncs
import os, sys
import numpy as np
from pylab import *

### THESE NEED TO BE DEFINE BEFORE CALLING THIS SCRIPT ###
#animalName   = 'saja125'
#ephysSession = '2012-03-06_17-16-57'
#behavSession = '20120306a'
#tetrodes = [1,2,3,4,5,6,7,8]

markerSize = 2
fontSize = 14

rasterDeltaT = 0.1e-3            # sec
timeRange = np.array([-0.2,0.6])

# -- Load events from Neuralynx --
dataDir = os.path.join(settings.EPHYS_PATH,'%s/%s/'%(animalName,ephysSession))
clustersDir = os.path.join(settings.EPHYS_PATH,'%s/%s_kk/'%(animalName,ephysSession))
eventsFile = os.path.join(dataDir,'Events.nev')
events = loadneuralynx.DataEvents(eventsFile)
bitTRIALIND = 3  # bitID starting from 0
trialEvents = (events.valueTTL & (1<<bitTRIALIND)) != 0
trialStartTimeNL = 1e-6*events.timestamps[trialEvents]
bitTARGET  = 5 
targetEvents = (events.valueTTL & (1<<bitTARGET)) != 0
targetTimeNL = 1e-6*events.timestamps[targetEvents]
#targetFreqInd = events.valueTTL[trialEvents]>>8  # 

# -- Load events from behavior --
behavDataDir = os.path.join(settings.BEHAVIOR_PATH,'%s/'%(animalName))
behavFileName = 'data_saja_tuningcurve_santiago_%s_%s.h5'%(animalName,behavSession)
behavFile = os.path.join(behavDataDir,behavFileName)
behavData = loadbehavior.TuningBehaviorData(behavFile)
behavData.extract_event_times()
behavData.align_to_ephys(trialStartTimeNL)

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


trialsOfInterest = argsort(behavData['SoundFreq'])

xLabel = 'Time from sound onset (ms)'
eventOfInterest = behavData.targetOnsetTime[trialsOfInterest] - \
                  behavData.trialStartTime[trialsOfInterest] + \
                  behavData.trialStartTimeEphys[trialsOfInterest]
'''
eventOfInterest = behavData.trialStartTimeEphys[trialsOfInterest]
'''

timeVec = np.arange(timeRange[0],timeRange[-1]+rasterDeltaT,rasterDeltaT)

siteIDvec = range(len(tetrodes))
siteID = 0
print("Use keys '<' and '>' to advance, 'q' to quit")
while True:
    tetrode = tetrodes[siteID]
    # -- Load spikes --
    tetrodeFile = os.path.join(dataDir,'TT%d.ntt'%tetrode)
    dataTT = loadneuralynx.DataTetrode(tetrodeFile)
    dataTT.timestamps = dataTT.timestamps.astype(np.float64)*1e-6  # in sec

    # -- Load clusters if required --
    if (clustersEachTetrode is not None) & clustersEachTetrode.has_key(tetrode):
        clustersFile = os.path.join(clustersDir,'TT%d.clu.1'%tetrode)
        dataTT.set_clusters(clustersFile)
        spikeInds = extrafuncs.ismember(dataTT.clusters,clustersEachTetrode[tetrode])
        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
          spikesanalysis.eventlocked_spiketimes(dataTT.timestamps[spikeInds],eventOfInterest,timeRange)
    else:
        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
          spikesanalysis.eventlocked_spiketimes(dataTT.timestamps,eventOfInterest,timeRange)

    clf()
    gcf().set_facecolor('w')
    ax = axes([0.2,0.2,0.7,0.7])
    pR = plot(1e3*spikeTimesFromEventOnset,trialIndexForEachSpike,'.k')
    pR[0].set_markersize(markerSize)
    xlim(1e3*timeRange)
    ylim([0,len(eventOfInterest)-1])
    titleString = '%s [%s] T%d'%(animalName,ephysSession,tetrode)
    title(titleString,fontsize=fontSize-2)
    print titleString

    xlabel(xLabel,fontsize=fontSize)
    #ax.set_yticklabels('')
    ylabel('Trial',fontsize=fontSize)

    # -- Draw background for each frequency --
    FreqEachTrial = behavData['SoundFreq'][trialsOfInterest]
    possibleFreqs = unique(behavData['SoundFreq'])
    lastTrialEachFreq = flatnonzero(diff(FreqEachTrial))
    firstTrialEachFreq = hstack((0,lastTrialEachFreq+1))
    lastTrialEachFreq = hstack((lastTrialEachFreq, len(FreqEachTrial)))
    # NOTE: indexes are from 0 to N+1 (python style)
    nFreq = len(firstTrialEachFreq)
    hold(True)
    possibleColors = ['0.95','0.9']
    for indf in range(nFreq):
        xpos = hstack((xlim(),xlim()[::-1]))
        ypos = hstack((tile(firstTrialEachFreq[indf]-0.5,2),tile(lastTrialEachFreq[indf]+0.5,2)))
        thisColor = possibleColors[mod(indf,2)]
        bg = fill(xpos,ypos,color=thisColor,ec='none')
        ht = text(0.99*xlim()[1],firstTrialEachFreq[indf],'%d'%possibleFreqs[indf],color='b',
                  backgroundcolor='w',horizontalalignment='right')
    hold(False)
    ylim([0,nTrials])

    draw()
    show()
    try:
        siteID = extraplots.loop_keycontrol(gcf(),siteIDvec,siteID)
    except StopIteration:
        break

