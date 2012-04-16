# -*- coding: utf-8 -*-
'''
Some statistics functions.

binofit: Parameter estimates and confidence intervals for binomial data.
'''

import scipy.stats
import numpy as np
import sys
#from np.core.multiarray import *   # FUNCTIONS: where
#from np.core.fromnumeric import *  # FUNCTIONS: nonzero
#from matplotlib.mlab import *          # FUNCTIONS: find
#from np.lib.shape_base import *     # FUNCTIONS: vstack

def binofit(xArray,nArray,alpha):
    '''Parameter estimates and confidence intervals for binomial data.
    (p,ci) = binofit(x,N,alpha)
    '''
    if isinstance(xArray,np.ndarray):
        origShape = xArray.shape
        Psuccess = np.empty(xArray.size)
        ConfIntervals = np.empty((xArray.size,2))
        if not isinstance(alpha,np.ndarray):
            alpha = alpha*np.ones(origShape)
        for inde in range(xArray.size):
            Psuccess[inde],ConfIntervals[inde,:] = binofit_scalar(xArray.flat[inde],
                                                                nArray.flat[inde],
                                                                alpha.flat[inde])
        return(Psuccess.reshape(origShape),ConfIntervals.reshape(origShape+(2,)))
    else:
        return binofit_scalar(xArray,nArray,alpha)
        

def binofit_scalar(x,n,alpha):
    '''Parameter estimates and confidence intervals for binomial data.
    (p,ci) = binofit(x,N,alpha)

    Source: Matlab's binofit.m
    Reference:
      [1]  Johnson, Norman L., Kotz, Samuel, & Kemp, Adrienne W.,
      "Univariate Discrete Distributions, Second Edition", Wiley
      1992 p. 124-130.
    http://books.google.com/books?id=JchiadWLnykC&printsec=frontcover&dq=Univariate+Discrete+Distributions#PPA131,M1

    Re-written by Santiago Jaramillo - 2008.10.06
    '''
    
    if n<1:
        Psuccess = np.NaN
        ConfIntervals = (np.NaN,np.NaN)
    else:
        Psuccess = float(x)/n
        nu1 = 2*x
        nu2 = 2*(n-x+1);
        F = scipy.stats.f.ppf(alpha/2,nu1,nu2)
        lb  = (nu1*F)/(nu2 + nu1*F)
        if x==0: lb=0
        nu1 = 2*(x+1)
        nu2 = 2*(n-x)
        F = scipy.stats.f.ppf(1-alpha/2,nu1,nu2)
        ub = (nu1*F)/(nu2 + nu1*F)
        if x==n: ub=1
        ConfIntervals = (lb,ub);
    return (Psuccess,ConfIntervals)


def twobinotest(Ntrials,Nhits,Nsets=int(1e4)):
    '''Using resampling (is it bootstrapping, jackknifing or permutation test?)
       to find p-value for comparing two binomial distributions.
       Ntrials and Nhits have to be numpy arrays with two elements each.

       See also scipy.stats.fisher_exact()
    '''
    print 'Using %d bootstrap sets. Calculating...'%Nsets,
    sys.stdout.flush()
    
    # -- Assume the first one is smaller than the second one --
    # I don't think that assumption is needed (2009.03.26)
    DiffMean = np.diff(Nhits.astype(float)/Ntrials)
    DiffMeanPos = abs(DiffMean)
    DiffMeanNeg = -DiffMeanPos
    
    NtrialsTotal = Ntrials.sum()
    Samples = np.zeros(NtrialsTotal,dtype='int')
    Samples[0:Nhits.sum()] = 1
    
    NhitsBS = np.zeros((Nsets,2),dtype='int')
    for inds in range(int(Nsets)):
        NewOrder = np.random.randint(NtrialsTotal,size=(NtrialsTotal))
        NewHits0 = Samples[NewOrder[0:Ntrials[0]]].sum()
        NewHits1 = Samples[NewOrder[Ntrials[0]+1:NtrialsTotal]].sum()
        NhitsBS[inds,:] = np.array([NewHits0,NewHits1])
    
    MeansBS = NhitsBS.astype('float') / Ntrials[np.newaxis]
    DiffMeansBS = np.diff(MeansBS)
    
    #Nbeyondvalue = sum(DiffMeansBS>=DiffMean)
    Nbeyondvalue = sum( (DiffMeansBS>=DiffMeanPos) | (DiffMeansBS<=DiffMeanNeg) )
    pValue = Nbeyondvalue.astype('float')/Nsets
    
    print 'Done!  pvalue = %0.6f'%pValue
    return pValue
    #raise(Exception('DEBUG'))

if __name__ == "__main__":
    '''Testing the functions on this module'''
    
    # -- Test twobinotest --
    Ntrials = np.array([100,200])
    Nhits = np.array([70,120])

    pValue = twobinotest(Ntrials,Nhits)
    
    '''
NewMeans = NewSums ./ repmat([Ntrials1,Ntrials2],Nsets,1);
    
    
    for indset=1:Nsets
    %%%RandVec = rand(1,NtrialsTotal);
    %%%[ignore,NewOrder] = sort(RandVec);
    NewOrder = ceil(NtrialsTotal*rand(NtrialsTotal,1));

  NewSum1 = sum(AllSamples(NewOrder(1:Ntrials1)));
  NewSum2 = sum(AllSamples(NewOrder(Ntrials1+1:NtrialsTotal)));
  NewSums(indset,:) = [NewSum1,NewSum2];
end

NewMeans = NewSums ./ repmat([Ntrials1,Ntrials2],Nsets,1);

NewMeansDiff = NewMeans(:,1)-NewMeans(:,2);

hist(NewMeansDiff,30);

Nbeyondvalue = sum(NewMeansDiff>=DiffValueToTest);

pValue = Nbeyondvalue/Nsets;
fprintf('pValue = %0.6f \n',pValue);
'''
    
    
    
