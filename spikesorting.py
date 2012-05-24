'''
Methods and classes for spike sorting and creating reports.
'''

import matplotlib.pyplot as plt
from extracellpy import settings
reload(settings)
from extracellpy import loadneuralynx
from extracellpy import extraplots
import numpy as np
import os

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

BEHAVIORPATH = settings.BEHAVIOR_PATH
SAMPLES_PER_SPIKE = 32
N_CHANNELS = 4

def plot_isi_loghist(timeStamps,nBins=350,fontsize=8):
    '''
    Plot histogram of inter-spike interval (in msec, log scale)

    Parameters
    ----------
    timeStamps : array (assumed to be integers in microsec)
    '''
    fontsizeLegend = fontsize
    xLims = [1e-1,1e4]
    ax = plt.gca()
    ISI = np.diff(timeStamps)
    if np.any(ISI<0):
        raise 'Times of events are not ordered (or there is at least one repeated).'
    if len(ISI)==0:  # Hack in case there is only one spike
        ISI = np.array(1)
    #if len(timeStamps)<2:
    #    return (0,0,0) ### FIXME: what to do when only one spike?
    logISI = np.log10(ISI)
    [ISIhistogram,ISIbinsLog] = np.histogram(logISI,bins=nBins)
    ISIbins = 1e-3*(10**ISIbinsLog[:-1]) # Conversion to msec
    percentViolation = 100*np.mean(ISI<1e3) # Assumes ISI in usec
    percentViolation2 = 100*np.mean(ISI<2e3) # Assumes ISI in usec
    
    hp, = plt.semilogx(ISIbins,ISIhistogram,color='k')
    #plt.ylabel('Cluster %d'%SelectedCluster)
    plt.setp(hp,lw=0.5,color='k')
    yLims = plt.ylim()
    plt.xlim(xLims)
    plt.text(0.15,0.85*yLims[-1],'N=%d'%len(timeStamps),fontsize=fontsizeLegend,va='top')
    #plt.text(0.15,0.7*yLims[-1],'%0.2f%%'%percentViolation,fontsize=fontsizeLegend)
    plt.text(0.15,0.6*yLims[-1],'%0.2f%%\n%0.2f%%'%(percentViolation,percentViolation2),
             fontsize=fontsizeLegend,va='top')
    #'VerticalAlignment','top','HorizontalAlignment','left','FontSize',FontSizeAxes);
    ax.xaxis.grid(True)
    ax.yaxis.grid(False)
    plt.xlabel('Interspike interval (ms)')
    ax.set_yticks(plt.ylim())
    extraplots.set_ticks_fontsize(ax,fontsize)
    return (hp,ISIhistogram,ISIbins)

def plot_events_in_time(timeStamps,nBins=50,fontsize=8):
    '''
    Plot histogram of inter-spike interval (in msec, log scale)

    Parameters
    ----------
    timeStamps : array (assumed to be integers in microsec)
    '''
    ax = plt.gca()
    timeBinEdges = np.linspace(timeStamps[0],timeStamps[-1],nBins) # in microsec
    # FIXME: xLimits depend on the time of the first spike (not of recording)
    (nEvents,binEdges) = np.histogram(timeStamps,bins=timeBinEdges)
    hp, = plt.plot(1e-6/60 * (binEdges-timeStamps[0]),np.r_[nEvents,0],drawstyle='steps-post')
    plt.setp(hp,lw=1,color='k')
    plt.xlabel('Time (min)')
    plt.axis('tight')
    ax.set_yticks(plt.ylim())
    extraplots.boxoff(ax)
    extraplots.set_ticks_fontsize(ax,fontsize)
    return hp

def plot_waveforms(waveforms,ntraces=40,fontsize=8):
    '''
    Plot waveforms given array of shape (nChannels,nSamplesPerSpike,nSpikes)
    '''
    (nChannels,nSamplesPerSpike,nSpikes) = waveforms.shape
    meanWaveforms = np.mean(waveforms,axis=2)
    scalebarSize = meanWaveforms.max()
    
    spikesToPlot = np.random.randint(nSpikes,size=ntraces)
    xRange = np.arange(nSamplesPerSpike)
    for indc in range(nChannels):
        newXrange = xRange+indc*(nSamplesPerSpike+2)
        wavesToPlot = waveforms[indc,:,spikesToPlot].T
        plt.plot(newXrange,wavesToPlot,color='k',lw=0.4,clip_on=False)
        plt.hold(True)
        plt.plot(newXrange,meanWaveforms[indc,:],color='0.75',lw=1.5,clip_on=False)
    plt.plot(2*[-7],[0,scalebarSize],color='0.5',lw=2)
    percentOfMax = 100*(scalebarSize/2**15)
    plt.text(-10,scalebarSize/2,'%d%%\nmax'%np.round(percentOfMax),
             ha='right',va='center',ma='center',fontsize=fontsize)
    plt.hold(False)
    plt.axis('off')

def plot_projections(waveforms,npoints=200):
    (nChannels,nSamplesPerSpike,nSpikes) = waveforms.shape
    spikesToPlot = np.random.randint(nSpikes,size=npoints)
    peaks = np.max(waveforms[:,:,spikesToPlot],axis=1)
    plt.plot(peaks[0,:],peaks[1,:],'.k',ms=0.5)
    plt.hold(True)
    plt.plot(-peaks[2,:],peaks[3,:],'.k',ms=0.5)
    plt.plot(0,0,'+',color='0.5')
    plt.hold(False)
    plt.axis('off')
    
   
class ClusterReportFromData(object):
    def __init__(self,dataTT,outputDir=None,showfig=True,nrows=12):
        self.dataTT = dataTT
        self.nSpikes = 0
        self.clustersList = []
        self.nClusters = 0
        self.spikesEachCluster = [] # Bool
        #self.fig = plt.figure(fignum)
        self.fig = plt.gcf()
        self.nRows = nrows
        self.set_parameters() # From dataTT

        self.plot_report(showfig)
        if outputDir is not None:
            self.save_report(outputDir)
    def set_parameters(self):
        self.nSpikes = len(self.dataTT.timestamps)
        self.clustersList = np.unique(self.dataTT.clusters)
        self.nClusters = len(self.clustersList)
        self.find_spikes_each_cluster()
    def __str__(self):
        return '%d clusters'%(self.nClusters)
    def find_spikes_each_cluster(self):
        self.spikesEachCluster = np.empty((self.nClusters,self.nSpikes),dtype=bool)
        for indc,clusterID in enumerate(self.clustersList):
            self.spikesEachCluster[indc,:] = (self.dataTT.clusters==clusterID)
    def plot_report(self,showfig=False):
        #plt.figure(self.fig)
        self.fig.clf()
        self.fig.set_facecolor('w')
        nCols = 3
        nRows = self.nRows
        #for indc,clusterID in enumerate(self.clustersList[:2]):
        for indc,clusterID in enumerate(self.clustersList):
            #print('Preparing cluster %d'%clusterID)
            tsThisCluster = self.dataTT.timestamps[self.spikesEachCluster[indc,:]]
            wavesThisCluster = self.dataTT.samples[:,:,self.spikesEachCluster[indc,:]]
            # -- Plot ISI histogram --
            plt.subplot(self.nRows,nCols,indc*nCols+1)
            plot_isi_loghist(tsThisCluster)
            if indc<(self.nClusters-1):
                plt.xlabel('')
                plt.gca().set_xticklabels('')
            plt.ylabel('c%d'%clusterID)
            # -- Plot events in time --
            plt.subplot(2*self.nRows,nCols,2*(indc*nCols)+6)
            plot_events_in_time(tsThisCluster)
            if indc<(self.nClusters-1):
                plt.xlabel('')
                plt.gca().set_xticklabels('')
            # -- Plot projections --
            plt.subplot(2*self.nRows,nCols,2*(indc*nCols)+3)
            plot_projections(wavesThisCluster)  
            # -- Plot waveforms --
            plt.subplot(self.nRows,nCols,indc*nCols+2)
            plot_waveforms(wavesThisCluster)
        figTitle = self.get_title()
        plt.figtext(0.5,0.92, figTitle,ha='center',fontweight='bold',fontsize=10)
        if showfig:
            #plt.draw()
            plt.show()
    def get_title(self):
        return ''
    def get_default_filename(self,figformat):
        return 'clusterReport.%s'%(figformat)
    def save_report(self,outputdir,filename=None,figformat=None):
        self.fig.set_size_inches((8.5,11))
        if figformat is None:
            figformat = 'png' #'png' #'pdf' #'svg'
        if filename is None:
            filename = self.get_default_filename(figformat)
        fullFileName = os.path.join(outputdir,filename)
        print 'Saving figure to %s'%fullFileName
        self.fig.savefig(fullFileName,format=figformat)
        #plt.close(self.fig)
        ###def closefig(self):


class ClusterReportTetrode(ClusterReportFromData):
    def __init__(self,animalName,ephysSession,tetrode,outputDir=None,showfig=False,nrows=12):
        self.animalName = animalName
        self.ephysSession = ephysSession
        self.tetrode = tetrode
        self.dataDir = ''
        self.clustersFile = ''
        self.tetrodeFile = ''
        #self.dataTT = []

        self.load_data()
        super(ClusterReportTetrode, self).__init__(self.dataTT,outputDir=outputDir,
                                                   showfig=showfig,nrows=nrows)
    def load_data(self):
        self.dataDir = os.path.join(settings.EPHYS_PATH,'%s/%s/'%(self.animalName,self.ephysSession))
        clustersDir = os.path.join(settings.EPHYS_PATH,'%s/%s_kk/'%(self.animalName,self.ephysSession))
        self.tetrodeFile = os.path.join(self.dataDir,'TT%d.ntt'%self.tetrode)
        print 'Loading data %s'%(self.tetrodeFile)
        dataTT = loadneuralynx.DataTetrode(self.tetrodeFile,readWaves=True)
        #dataTT.timestamps = dataTT.timestamps.astype(np.float64)*1e-6  # in sec
        ### The following line is not needed anymore (not done when loading data)
        #dataTT.samples = dataTT.samples.reshape((N_CHANNELS,SAMPLES_PER_SPIKE,-1),order='F')
        # -- Load clusters --
        self.clustersFile = os.path.join(clustersDir,'TT%d.clu.1'%self.tetrode)
        dataTT.set_clusters(self.clustersFile)
        self.dataTT = dataTT
    def get_title(self):
        return self.dataDir+' (T%d)'%self.tetrode  #tetrodeFile
    def __str__(self):
        return '%s  %s  T%d\n%d clusters'%(self.animalName,self.ephysSession,self.tetrode,self.nClusters)
    def get_default_filename(self,figformat):
        return '%s_%s_T%02d.%s'%(self.animalName,self.ephysSession,self.tetrode,figformat)

def save_all_reports(animalName,ephysSession,tetrodes,outputDir):
    if not os.path.exists(outputDir):
        print 'Creating output directory: %s'%(outputDir)
        os.makedirs(outputDir)
    for onetetrode in tetrodes:
        sreport = ClusterReportTetrode(animalName,ephysSession,onetetrode)
        sreport.save_report(outputDir)

def merge_kk_clusters(animalName,ephysSession,tetrode,clustersToMerge,reportDir=None):
    dataDir = os.path.join(settings.EPHYS_PATH,'%s/%s_kk/'%(animalName,ephysSession))
    #reportDir = os.path.join(settings.EPHYS_PATH,'%s/%s_reportkk/'%(animalName,ephysSession))
    if reportDir is None:
        reportDir = os.path.join(settings.PROCESSED_REVERSAL_PATH,settings.CLUSTERS_REPORTS_DIR)
    fileName = 'TT%d.clu.1'%(tetrode)
    fullFileName = os.path.join(dataDir,fileName)
    backupFileName = os.path.join(dataDir,fileName+'.orig')
    # --- Make backup of original cluster file ---
    print 'Making backup to %s'%backupFileName
    os.system('rsync -a %s %s'%(fullFileName,backupFileName))
    # --- Load cluster data, replace and resave ---
    clusterData = np.fromfile(fullFileName,dtype='int32',sep='\n')
    indNoiseSpike = np.flatnonzero(clusterData==1)[0]
    clusterData[clusterData==clustersToMerge[1]] = clustersToMerge[0]
    clusterData[indNoiseSpike] = clustersToMerge[1]
    clusterData.tofile(fullFileName,sep='\n',format='%d')
    # -- Create report --
    print 'Creating report in %s'%reportDir
    ClusterReportTetrode(animalName,ephysSession,tetrode,reportDir)


if __name__ == "__main__":
    CASE = 1
    if CASE==1:
        animalName   = 'saja125'
        ephysSession = '2012-01-31_14-37-44'
        tetrode = 6
        sreport = ClusterReportTetrode(animalName,ephysSession,tetrode,'/tmp/reports')
        #sreport.save_report('/tmp/reports/')
        #sreport.closefig()
    elif CASE==2:
        animalName   = 'saja125'
        ephysSession = '2012-04-23_16-10-15'
        #save_all_reports(animalName,ephysSession,np.arange(1,8+1),'/var/data/neuralynx/saja125_processed/cluster_reports')
        save_all_reports(animalName,ephysSession,[2],'/tmp/reports')
    elif CASE==3:
        animalName   = 'saja125'
        ephysSession = '2012-04-23_16-10-15'
        tetrode = 2
        #merge_kk_clusters(animalName,ephysSession,tetrode,[2,5],reportDir='/tmp/reports')
        #merge_kk_clusters(animalName,ephysSession,tetrode,[2,10],reportDir='/tmp/reports')

        
'''
animalName   = 'saja125'
ephysSession = '2012-02-07_14-18-20'
tetrode = 2
sreport = ClusterReportTetrode(animalName,ephysSession,tetrode,'/tmp/reports')
'''
