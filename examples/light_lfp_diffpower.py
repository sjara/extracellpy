# -*- coding: utf-8 -*-
'''
Look at LFP evoked by ChR2 stimulation.
'''

from extracellpy import settings
from extracellpy import extraplots
from extracellpy import loadneuralynx
#from scipy.signal import lfilter, bessel
import numpy as np
import os
from pylab import *
reload(loadneuralynx)

'''
### IT NEEDS: ###
ephysSession =  '2012-03-06_18-23-51' # (L2) Diff power every 20 trials (x4)
electrodeList = [5,6]
trialRange = [range(0,20),range(20,40),range(40,60),range(60,80)]
rangeLabel = ['1.0','1.2','1.5','2.0']
timeRange = [0.390,0.500]
'''

OUTPUTDIR = '/tmp/'
FONTSIZE = 16

lightOnset = 0.4 # sec
AmpGain = 1245  ### NOT CORRECT
ADMaxValue = 32767
ADBitVolts = 0.0000000306
print('The amplitude units are arbitrary (not read from file yet)')

#timeRange = np.array([0.380,0.650])
#timeRange = np.array([0.300,0.650])
#timeRange = np.array([-1,2])
#timeRange = np.array([-1,2])

dataDir = os.path.join(settings.EPHYS_PATH,'%s/%s/'%(animalName,ephysSession))
eventsFile = os.path.join(dataDir,'Events.nev')
events = loadneuralynx.DataEvents(eventsFile)

bitTRIALIND = 3      # TrialIndicator (bitID starting from 0)
bitPHOTOSTIMIND = 4  # PhotoStimIndicator (bitID starting from 0)
bitTARGETIND = 5     # TargetIndicator (bitID starting from 0)

(trialOnset,trialOffset) = events.find_bit_changes(bitTRIALIND)
#(photoOnset,photoOffset) = events.find_bit_changes(bitPHOTOSTIMIND)
trialOnsetTime = 1e-6*events.timestamps[trialOnset]

clf()
electrodeIDvec = range(len(electrodeList))
electrodeID = 0
print("Use keys '<' and '>' to advance, 'q' to quit")
while True:
    electrode = electrodeList[electrodeID]
    #ax = axes([0.1,0.2,0.8,0.6])
    cla()
    
    fileName = 'CSC%d.ncs'%electrode
    contFile = os.path.join(dataDir,fileName)
    dataLFP = loadneuralynx.DataCont(contFile)
    if not len(dataLFP.samples):
        continue
    (lockedLFP,timeVec) = dataLFP.lock_to_event(trialOnsetTime,timeRange)
    baselineSamples = np.arange(*np.searchsorted(timeVec,baselineRange))

    hold(True)
    hpall=[]
    for powerind in range(len(trialRange)):
        meanLFP = np.mean(lockedLFP[trialRange[powerind]],axis=0)
        #meanLFP = np.median(lockedLFP[trialRange[powerind]],axis=0)
        meanBaseline = np.mean(meanLFP[baselineSamples])
        
        if 0:
            # -- Filter signal (high-pass) --
            highpassFreq = 10 # Hz
            Wn = highpassFreq*2/dataLFP.samplingRate # Normalized
            [b,a]=bessel(4, Wn, btype='high')
            meanLFP = lfilter(b,a,meanLFP)

        #subplot(2,2,powerind+1)
        hp, = plot(1e3*(timeVec-lightOnset),1e6*ADBitVolts*(meanLFP-meanBaseline),'-')
        hpall.append(hp)

    xlabel('Time from photo-stimulation onset',fontsize=FONTSIZE)
    #ax.set_xticks(range(-60,250,20))
    #xlim([-40,200])
    xlim(1e3*np.array(timeRange)-400)
    legend(rangeLabel)
    title(fileName)
    hold(False)
    draw()
    show()
    try:
        electrodeID = extraplots.loop_keycontrol(gcf(),electrodeIDvec,electrodeID)
    except StopIteration:
        break


# -- Save figure --
if 0:
    gcf().set_frameon(False)
    gcf().set_size_inches((12,4))
    figFormat = 'svg' # 'svg' #'pdf' #'svg' #
    savefig(os.path.join(OUTPUTDIR,'lightLFP.%s'%(figFormat)),format=figFormat)
    gcf().set_frameon(True)
    #draw(); show()
