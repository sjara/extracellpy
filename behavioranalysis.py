#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Tools for analyzing behavioral data.
'''

import settings
from pylab import *
import numpy as np
import datetime
import sys, os
import cPickle as pickle
import gzip

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

BEHAVIORPATH = settings.BEHAVIOR_PATH


def parse_isodate(dateStr):
    '''Convert ISO-formatted date (e.g., 2011-12-31) to python date.    
    '''
    dateElems = dateStr.split('-')
    dateElems = [int(x) for x in dateElems]
    return datetime.date(*dateElems)

def reshape_by_session(allBehavData,maxNtrials=1500):
    '''Reshape behavior data from many sessions, keeping only valid trials
    and organizing in masked arrays where each column is one session'''
    nSessions = len(unique(allBehavData['sessionID']))
    allBehavDataPerSession = {}
    mask = np.ones((maxNtrials,nSessions),dtype=bool)
    for key in allBehavData.keys():
        typeThisKey = allBehavData[key].dtype
        allBehavDataPerSession[key] = np.empty((maxNtrials,nSessions),
                                               dtype=typeThisKey)
    for inds in range(nSessions):
        validTrialsThisSession = (allBehavData['sessionID']==inds) & \
                                 (~allBehavData['early'].astype(bool))
        nValidThisSession = sum(validTrialsThisSession)
        if nValidThisSession>maxNtrials:
            raise ValueError('Session has more trials than max allowed')
        for key in allBehavData.keys():
            allBehavDataPerSession[key][:nValidThisSession,inds] = \
                allBehavData[key][validTrialsThisSession]   
            mask[:nValidThisSession,inds]=False
        allBehavDataPerSession['mask']=mask
    return allBehavDataPerSession
    

def reshape_by_session_OLD(allBehavData):
    '''Reshape behavior data from many sessions, keeping only valid trials
    and organizing in masked arrays where each column is one session'''
    maxNtrials = 1500
    nSessions = unique(allBehavData['sessionID'])[-1]
    allBehavDataPerSession = {}
    for key in allBehavData.keys():
        typeThisKey = allBehavData[key].dtype
        allBehavDataPerSession[key] = ma.masked_array(np.zeros((maxNtrials,nSessions)),
                                                      dtype=typeThisKey)
    for inds in range(nSessions):
        validTrialsThisSession = (allBehavData['sessionID']==inds) & \
                                 (~allBehavData['early'].astype(bool))
        nValidThisSession = sum(validTrialsThisSession)
        if nValidThisSession>maxNtrials:
            raise ValueError('Session has more trials than max allowed')
        for key in allBehavData.keys():
            allBehavDataPerSession[key][:nValidThisSession,inds] = \
                allBehavData[key][validTrialsThisSession]   
            allBehavDataPerSession[key][nValidThisSession:,inds]=ma.masked
    return allBehavDataPerSession
    

def save_many_sessions_reversal(animalNames,datesRange,outputFile,paramsToSave,
                                dtypes=[],dataClass=None,datesList=[],compression=None,
                                saveDate=False):
    '''Save selected parameters from a set of sessions.
    Example parameters:
    animalNames = ['saja065','saja099','saja101']
    datesRange = ['2011-03-14','2011-03-18']
    paramsToSave = ['TargetFreq','CurrentBlock']
    outputFile = '/tmp/allBehavData.pkl'

    It returns:
    allBehavData (a dict with all session concatenated) and two additional
    keys, animalID and sessionID

    It saves the dict allBehavData serialized (using pickle).
    '''
    if not dataClass:
        from extracellpy import loadbehavior
        dataClass = loadbehavior.ReversalBehaviorData
    if datesList:
        allDates = datesList
    else:
        datesLims = [parse_isodate(dateStr) for dateStr in datesRange]
        allDates = [datesLims[0]+datetime.timedelta(n) \
                        for n in range((datesLims[-1]-datesLims[0]).days+1)]

    nAnimals = len(animalNames)
    dataFileFormat = 'data_saja_reversal_santiago_%s_%s.h5'

    allBehavData = {}
    allBehavData['sessionID'] = np.empty(0,dtype='i2')
    allBehavData['animalID'] = np.empty(0,dtype='i1')
    sessionArray = []
    inds = 0
    if not dtypes:
        for indp,key in enumerate(paramsToSave):
            pass
            # FIXME: Need to find dtype for each
            #dtypes.append(behavData[key].dtype)
    for key,dtype in zip(paramsToSave,dtypes):
        allBehavData[key] = np.empty(0,dtype=dtype)
    for inda in range(nAnimals):
        animalName = animalNames[inda]
        for indd,oneDate in enumerate(allDates):
            behavSession = oneDate.strftime('%Y%m%da')
            behavDataDir = os.path.join(BEHAVIORPATH,'%s/'%(animalName))
            behavFileName = dataFileFormat%(animalName,behavSession)
            behavFile = os.path.join(behavDataDir,behavFileName)
            try:
                behavData = dataClass(behavFile)
            except IOError:
                print('No file named %s'%(behavFile))
                continue
            print 'Loaded %s'%behavFile
            nTrials = behavData['nTrials']
            sessionArray.append(behavSession)
            behavData.extract_event_times() # This is needed (it shouldn't be)
            behavData.find_trials_each_type()
            behavData.find_trial_index_each_block()
            for key,dtype in zip(paramsToSave,dtypes):
                if key in behavData:
                    allBehavData[key] = np.concatenate((allBehavData[key],
                                                        behavData[key]))
                elif hasattr(behavData,key):
                    thisdtype = allBehavData[key].dtype
                    allBehavData[key] = np.concatenate((allBehavData[key],
                                                        getattr(behavData,key)))
                else:
                    raise ValueError('Parameter not present in behavior data.')
                allBehavData[key] = allBehavData[key].astype(dtype)
            allBehavData['sessionID'] = np.concatenate((allBehavData['sessionID'],
                                                        np.tile(inds,nTrials)))
            allBehavData['animalID'] = np.concatenate((allBehavData['animalID'],
                                                        np.tile(inda,nTrials)))
            inds += 1
    if saveDate:
        allBehavData['sessionName'] = np.array(sessionArray)
    if outputFile:
        save_dict_as_HDF5(outputFile,allBehavData,compression=compression)
    return allBehavData

def save_dict_as_HDF5(outputFile,allBehavData,compression=None):
    import h5py
    output = h5py.File(outputFile,'w')
    try:
        for key,val in allBehavData.iteritems():
            dset = output.create_dataset(key, data=val, dtype=val.dtype,
                                         compression=compression)
    except:
        output.close()
        raise
    output.close()

def load_dict_from_HDF5(dataFile):
    '''Use this function only for testing.
    The efficient way is to load only what is needed
    '''
    import h5py
    data = h5py.File(dataFile,'r')
    try:
        allBehavData = {}
        for key in data.keys():
            allBehavData[key] = data[key][...]
    except:
        data.close()
        raise
    data.close()
    return allBehavData

def save_many_sessions_reversal_OLD(animalNames,datesRange,outputFile,paramsToSave,
                                attribToSave=[],dataClass=None,datesList=[]):
    '''Save selected parameters from a set of sessions.
    Example parameters:
    animalNames = ['saja065','saja099','saja101']
    datesRange = ['2011-03-14','2011-03-18']
    paramsToSave = ['TargetFreq','CurrentBlock']
    outputFile = '/tmp/allBehavData.pkl'

    It returns:
    allBehavData (a dict with all session concatenated) and two additional
    keys, animalID and sessionID

    It saves the dict allBehavData serialized (using pickle).
    '''
    if not dataClass:
        from extracellpy import loadbehavior
        dataClass = loadbehavior.ReversalBehaviorData
    if datesList:
        allDates = datesList
    else:
        datesLims = [parse_isodate(dateStr) for dateStr in datesRange]
        allDates = [datesLims[0]+datetime.timedelta(n) \
                        for n in range((datesLims[-1]-datesLims[0]).days+1)]

    nAnimals = len(animalNames)
    dataFileFormat = 'data_saja_reversal_santiago_%s_%s.h5'

    allBehavData = {}
    allBehavData['sessionID'] = np.empty(0,dtype=int)
    allBehavData['animalID'] = np.empty(0,dtype=int)
    inds = 0
    for key in (paramsToSave+attribToSave):
        allBehavData[key] = np.empty(0)
    for inda in range(nAnimals):
        animalName = animalNames[inda]
        for indd,oneDate in enumerate(allDates):
            behavSession = oneDate.strftime('%Y%m%da')
            behavDataDir = os.path.join(BEHAVIORPATH,'%s/'%(animalName))
            behavFileName = dataFileFormat%(animalName,behavSession)
            behavFile = os.path.join(behavDataDir,behavFileName)
            try:
                behavData = dataClass(behavFile)
            except IOError:
                print('No file named %s'%(behavFile))
                continue
            print 'Loaded %s'%behavFileName
            nTrials = behavData['nTrials']
            behavData.extract_event_times() # This is needed (it shouldn't be)
            behavData.find_trials_each_type()
            behavData.find_trial_index_each_block()
            for key in paramsToSave:
                thisdtype = allBehavData[key].dtype
                allBehavData[key] = np.concatenate((allBehavData[key].astype(thisdtype),
                                                    behavData[key]))
            for attrib in attribToSave:
                thisdtype = allBehavData[key].dtype
                allBehavData[attrib] = np.concatenate((allBehavData[attrib].astype(thisdtype),
                                                       getattr(behavData,attrib)))
            allBehavData['sessionID'] = np.concatenate((allBehavData['sessionID'],
                                                        np.tile(inds,nTrials)))
            allBehavData['animalID'] = np.concatenate((allBehavData['animalID'],
                                                        np.tile(inda,nTrials)))
            inds += 1
    #outputFileName = 'allBehavData.pkl'
    #os.path.join(outputDir,outputFileName)
    if outputFile:
        output = gzip.open(outputFile, 'wb')
        pickle.dump(allBehavData, output)
        output.close()
    return allBehavData


def merge_sessions(allBehavData,keysToKeep):
    '''Create a dict of 1D arrays with data from many sessions
    allBehavData: dict (nParams) of lists (nSessions) of arrays (nTrials)
    fieldsToKeep: keys to include in output array
    '''
    mergedBehavData = {}
    for key in keysToKeep:
        mergedBehavData[key] = np.concatenate(allBehavData[key])
    return mergedBehavData

def merge_animals(allBehavData,keysToKeep):
    '''Create a dict of 1D arrays with data from many sessions
    many animals.
    allBehavData: dict (nParams) of lists (nAnimals) of 
                  lists(nSessions) of arrays (nTrials)
    fieldsToKeep: keys to include in output array
    '''
    mergedBehavData = {}
    for key in keysToKeep:
        mergedOneAnimal = []
        for allSessionsOneAnimal in allBehavData[key]:
            mergedOneAnimal.appen(np.concatenate(allSessionsOneAnimal))
        mergedBehavData[key] = np.concatenate(mergedOneAnimal)
    return mergedBehavData
    

def find_trials_each_type(psyCurveParameter,psyCurveParameterPossibleValues,\
                          currentBlock,currentBlockPossibleValues,validTrials=[]):
    '''
    Parameters
    ----------
    psyCurveParameter
    psyCurveParameterPossibleValues
    currentBlock
    currentBlockPossibleValues
    validTrials

    Returns
    -------
    trialsEachType: [nTrials,nParamValues,nBlocks] (boolean)

    See also:
    extracellpy.sessionanalysis.trials_by_condition()
    '''
    nTrials = len(psyCurveParameter)
    nValues = len(psyCurveParameterPossibleValues)
    nBlocks = len(currentBlockPossibleValues)
    if(not len(validTrials)):
        ValidTrials = ones(nTrials,dtype=bool)

    trialsEachBlock = zeros((nTrials,nBlocks),dtype=bool)
    for indb,blockID in enumerate(currentBlockPossibleValues):
        trialsEachBlock[:,indb] = (currentBlock==blockID)
    trialsEachType = zeros((nTrials,nValues,nBlocks),dtype=bool)
    for indval,paramValue in enumerate(psyCurveParameterPossibleValues):
        trialsThisValue = (psyCurveParameter==paramValue)
        for indb in range(nBlocks):
            trialsEachType[:,indval,indb] = trialsThisValue & \
                                            trialsEachBlock[:,indb] & validTrials
        
    return trialsEachType
                           


def OLD_find_trials_each_type(PsyCurveParameter,PsyCurveParameterPossibleValues,\
                          CurrentBlock,CurrentBlockLabels,\
                          CatchTrials=[],ValidTrials=[]):
    '''Find trials of each of four types:
    short-delay common, short-delay rare, long-delay common, long-delay rare
    
    It returns (TrialsEachType,TrialTypeStrings), where:
     TrialsEachType has shape (Ntrials,Npossibleparamsvalues,4)
     TrialTypeStrings is a list of 4 strings indicating each trial type
    '''
    Ntrials = len(PsyCurveParameter)
    Nvalues = len(PsyCurveParameterPossibleValues)

    if(not len(CatchTrials)):
        CatchTrials = zeros(Ntrials,dtype=bool)
    if(not len(ValidTrials)):
        ValidTrials = ones(Ntrials,dtype=bool)

    TrialTypeStrings = ['Short common','Short rare','Long common','Long rare']
    TrialsEachType = zeros((Ntrials,Nvalues,len(TrialTypeStrings)),dtype=bool)

    for indvalue in range(Nvalues):
        TrialsThisParamValue = (PsyCurveParameter==PsyCurveParameterPossibleValues[indvalue])
        # -- Short trials on short block --
        TrialsEachType[:,indvalue,0] = TrialsThisParamValue & ValidTrials & \
            (CurrentBlock==CurrentBlockLabels['ShortDelay']) & ~CatchTrials
        # -- Short trials on long block --
        TrialsEachType[:,indvalue,1] = TrialsThisParamValue & ValidTrials & \
            (CurrentBlock==CurrentBlockLabels['LongDelay']) & CatchTrials;
        # -- Long trials on long block --
        TrialsEachType[:,indvalue,2] = TrialsThisParamValue & ValidTrials & \
            (CurrentBlock==CurrentBlockLabels['LongDelay']) & ~CatchTrials;
        # -- Long trials on short block --
        TrialsEachType[:,indvalue,3] = TrialsThisParamValue & ValidTrials & \
            (CurrentBlock==CurrentBlockLabels['ShortDelay']) & CatchTrials;

    return (TrialsEachType,TrialTypeStrings)


def serialize(arrayOfStructs,fieldToSerialize):
    '''Take an array of structs and extract an array with one of the
       items in it.
    '''
    serializedArray = numpy.zeros(len(arrayOfStructs))
    for ind, oneStruct in enumerate(arrayOfStructs):
        serializedArray[ind]=oneStruct.__getattribute__(fieldToSerialize)
    return serializedArray


def plot_summary(behavData):
    '''Show summary of performance.
    Input is an object created by loadbehavior.ReversalBehaviorData()
    '''
    import matplotlib.pyplot as plt
    behavData.extract_event_times()    # In case it hasn't
    behavData.find_trials_each_type()  # In case it hasn't

    percentEarly = 100*sum(behavData.early)/float(behavData['nTrials'])

    g1 = behavData.lowFreqs & behavData.leftReward
    g2 = behavData.lowFreqs & behavData.rightReward
    g3 = behavData.highFreqs & behavData.leftReward
    g4 = behavData.highFreqs & behavData.rightReward

    trialsEachCond = np.array([g1,g2,g3,g4])
    validTrialsEachCond = trialsEachCond & (~behavData.early)
    correctTrialsEachCond = validTrialsEachCond & behavData.correct
    percentEarly = 100*sum(behavData.early)/float(len(behavData.early))

    ### WARNING: this gives a different number of valid trials than matlab!
    ###          check if first empty trial is included or not

    nCorrectEachCond = np.sum(correctTrialsEachCond,axis=1)
    nValidEachCond = np.sum(validTrialsEachCond,axis=1)

    perfEachCond = nCorrectEachCond/nValidEachCond.astype(float)

    freqValues = [behavData['FreqLow'][-1],
                  behavData['FreqMid'][-1],
                  behavData['FreqMid'][-1],
                  behavData['FreqHigh'][-1]]
    freqLabels = ['%0.1f'%(x/1000) for x in freqValues]

    # --- Plot results ---
    xPos = [0,1,3,4]
    #plt.clf()
    #ax = plt.axes()
    ax = plt.gca()
    ax.set_xlim([-1,5])
    ax.set_ylim([0,100])
    plt.hold(True)
    hline50 = plt.axhline(50,linestyle=':',color='k',zorder=-1)
    hline75 = plt.axhline(75,linestyle=':',color='k',zorder=-1)
    hbars = plt.bar(xPos,100*perfEachCond,align='center',fc=[0.8,0.8,0.8],ec='k')
    #htrials(indtype) = text(BarXpos(indtype),10,num2str(NtrialsEachCond(indc)));
    for thispos,thistext in zip(xPos,nValidEachCond):
        plt.text(thispos,10,str(thistext),ha='center')
    ax.set_ylabel('% correct',fontsize=12);
    ax.set_xticks(xPos)
    ax.set_xticklabels(freqLabels)

    titleStr = '%s [%s] %s\n'%(behavData['Subject'],behavData['Date'],behavData['HostName'])
    titleStr += '%d valid, %d%% early'%(sum(nValidEachCond),round(percentEarly))
    ax.set_title(titleStr,fontweight='bold')

    #plt.gcf().set_size_inches(4,5)

    plt.draw()
    plt.show()

    #plt.savefig(FIGDIR+'figreport.svg',format='svg')
    

def weibull(xval,alpha,beta):
    '''Weibull function
    alpha: bias
    beta: related to slope
    NOTE: this function isnot symmetric
    '''
    return 1 - np.exp(-pow(xval/alpha,beta))

def gaussianCDF(xval,mu,sigma):
    from scipy.stats import norm
    return norm.cdf(xval,mu,sigma)

def logistic(xval,alpha,beta):
    return 1/(1+np.exp(-(xval-alpha)/beta))

def psychfun(xval,alpha,beta,lamb,gamma):
    '''Psychometric function that allowing arbitrary asymptotes.
    alpha: bias
    beta : related to slope
    lamb : lapse term (up)
    gamma: lapse term (down)
    '''
    #return gamma + (1-gamma-lamb)*weibull(xval,alpha,beta)
    #return gamma + (1-gamma-lamb)*gaussianCDF(xval,alpha,beta)
    return gamma + (1-gamma-lamb)*logistic(xval,alpha,beta)
    
def negloglikepsych(theta,xVec,nTrials,nHits):
    '''
    theta: [alpha,beta,gamma,lamb]
    '''
    psychVals = psychfun(xVec,*theta)
    tmpval = psychVals[0]
    # -- Avoid saturation ---
    psychVals[psychVals<0.000001] = 0.000001
    psychVals[psychVals>0.999999] = 0.999999
    # -- Likelihood (only terms that depend on theta) --
    #np.log(comb(nTrials,nHits)) # Extra term (indep of theta)
    likeEach = nHits*np.log(psychVals) + \
               (nTrials-nHits)*np.log(1-psychVals)
    totalLikelihood = -np.sum(likeEach)
    # FIXME: The following tries to avoid the psych function to be
    #        outside [0,1], but affects a lot the function's smoothness
    totalLikelihood = totalLikelihood+100*(theta[3]<0)
    totalLikelihood = totalLikelihood+100*(theta[2]<0)
    #print psychVals
    #print totalLikelihood
    return totalLikelihood

def plot_psychcurve_fit_varsize(xValues,nTrials,nHits,curveParams=[],color='k'):
    '''
    Plot average performance for each value and fitted curve.
    '''
    yValues = nHits.astype(float)/nTrials
    xRange = xValues[-1]-xValues[1]
    mSize = nTrials.astype(float)/np.max(nTrials)
    for ind in range(len(xValues)):
        hp=plot(xValues[ind],yValues[ind],'o',mfc=color)
        setp(hp,mec=color,markersize=16*(mSize[ind]**0.2))
        hold(True)
    if curveParams:
        fitxval = np.linspace(xValues[0]-0.1*xRange,xValues[-1]+0.1*xRange,40)
        fityval = psychfun(fitxval,*curveParams)
        hp = plot(fitxval,fityval,'-',linewidth=2,color=color)
    ylim([-0.1,1.1])
    hline = axhline(0.5)
    setp(hline,linestyle=':',color='k')
    #grid(True)
    #hold(False)
    
def plot_psychcurve_fit(xValues,nTrials,nHits,curveParams=[],color='k'):
    '''
    Plot average performance for each value and fitted curve.
    '''
    solidXvalues = np.flatnonzero((nTrials/sum(nTrials).astype(float))>(1.0/len(nTrials)))
    yValues = nHits.astype(float)/nTrials
    xRange = xValues[-1]-xValues[1]
    hfit = []
    if len(curveParams):
        fitxval = np.linspace(xValues[0]-0.1*xRange,xValues[-1]+0.1*xRange,40)
        fityval = psychfun(fitxval,*curveParams)
        hfit = plot(fitxval,100*fityval,'-',linewidth=2,color=color)
        hold(True)
    hp = []
    for ind in range(len(xValues)):
        htemp,=plot(xValues[ind],100*yValues[ind],'o',mfc=color)
        hp.append(htemp)
        hold(True)
    setp(hp,mec=color,mfc='w',mew=2,markersize=6)
    for solid in solidXvalues:
        setp(hp[solid],mfc=color,markersize=8)
    #ylim([-0.1,1.1])
    ylim([-10,110])
    #hline = axhline(0.5)
    #setp(hline,linestyle=':',color='k')
    return (hp,hfit)

def fit_psychometric(nTrials,nHits,xValues,theta0):
    '''Maximum likelihood fit of psychometric function to behavior data.
    See also module "psignifit".

    Parameters
    ----------
    nTrials: vector with number of trials on each x-value
    nHits: vector with number of hits on each x-value
    xValues: vector with values of independent variable
    theta0: initial guess for parameters
            [alpha,beta,gamma,lamb]
            alpha: bias
            beta : related to slope
            gamma: lapse term (down)
            lamb : lapse term (up)
    Returns
    -------
    theta_hat: best estimate of parameters given data
    '''
    from scipy.optimize import fmin
    thetahat = fmin(negloglikepsych,theta0,args=[xValues,nTrials,nHits])
    ###,maxiter=1e6)
    return thetahat


if __name__ == "__main__":

    CASE = 4
    if CASE==1:
        '''Testing find_trials_each_type'''
        PsyCurveParameter = numpy.array([1,3,1,3])
        PsyCurveParameterPossibleValues = numpy.array([1,2,3])
        CurrentBlock = numpy.array([1,1,2,2])
        CurrentBlockLabels = {'ShortDelay': 1, 'LongDelay': 2}


        (TrialsEachType,TrialTypeStrings) = find_trials_each_type(PsyCurveParameter,PsyCurveParameterPossibleValues,\
                                CurrentBlock,CurrentBlockLabels)

        print "TrialsEachType:",TrialsEachType
        print "TrialTypeStrings:",TrialTypeStrings

    elif CASE==2:
        import loadbehavior
        bfile = '/var/data/BControlData/Data/santiago/saja100/data_saja_reversal_santiago_saja100_20111025a.h5'
        behavData = loadbehavior.ReversalBehaviorData(bfile)
        plot_summary(behavData)
    elif CASE==3:
        xValues = np.linspace(0,4,7)
        nHits = np.array([  1, 1, 2, 4, 26, 37, 37])
        nTrials = np.array([10, 10, 10, 20, 30, 40, 40])
        theta0 = [2,10,0.1,0.1]
        thetahat = fit_psychometric(nTrials,nHits,xValues,theta0)

        fitxval = linspace(xValues[0],xValues[-1],40)
        fityval = psychfun(fitxval,*thetahat)

        import matplotlib.pyplot as plt        
        plt.clf()
        plt.plot(xValues,nHits.astype(float)/nTrials,'o')
        plt.hold(True)
        plt.plot(fitxval,fityval,'-r',linewidth=2)
        plt.hold(False)
        plt.ylim([0,1])
        hline = axhline(0.5)
        setp(hline,linestyle=':',color='k')
        plt.draw()
        plt.show()
    elif CASE==4:
        animalNames = ['saja099','saja101']
        datesRange = ['2011-03-14','2011-03-15']
        outputFile = '/tmp/allBehavData.pkl'
        paramsToSave = ['TargetFreq','CurrentBlock']
        attribToSave = ['early']
        allBehav = save_many_sessions_reversal(animalNames,datesRange,
                                               outputFile,paramsToSave,
                                               attribToSave)
