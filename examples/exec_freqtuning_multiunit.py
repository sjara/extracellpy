'''
Execute for freqtuning_multiunit for a particular dataset.

Santiago Jaramillo - 2012-03-08
'''

clustersEachTetrode = None

SESSIONID = 1
if SESSIONID==0:
    animalName   = 'saja125'
    ephysSession = '2012-03-06_17-16-57'
    behavSession = '20120306a'
    tetrodes = [1,2,3,4,5,6,7,8]
elif SESSIONID==1:
    animalName   = 'saja125'
    ephysSession = '2012-04-19_15-12-22'
    behavSession = '20120419a'
    tetrodes = [1,2,3,4,5,6,7,8]
    clustersEachTetrode = {1:[2,3,4,5,7,8,10,11],2:[3,4,5,6,8,9,10,11],
                           3:[3,4,5,7,9,10], 4:[2,3,4,6,8,9,10,11],
                           5:[2,3,5,7,8,9], 6:[2,3,4,5,6,7,8,9,11],
                           7:[3,4,5,7,9,10,11,12], 8:[2,6,7,9,10,11]}
    #clustersEachTetrode = {6:[11]}
elif SESSIONID==2:
    animalName   = 'saja125'
    ephysSession = '2012-04-19_16-19-50'
    behavSession = '20120419b'
    tetrodes = [1,2,3,4,5,6,7,8]
    clustersEachTetrode = {1:[2,5,6,7,8,10,12],2:[2,3,4,5,6,7,9,10,11,12],
                           3:[3,5,7,8,11],4:[2,4,5,6,7,8,9,10],
                           5:[2,3,4,5,7,9],6:[2,3,5,7,8,9,10,11,12],
                           7:[2,3,4,6,7,8,9,10],8:[2,3,4,7,8,10,11]}
    #clustersEachTetrode = {6:[7]}


execfile('freqtuning_multiunit.py')
