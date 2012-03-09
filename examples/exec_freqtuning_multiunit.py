'''
Execute for freqtuning_multiunit for a particular dataset.

Santiago Jaramillo - 2012-03-08
'''

SESSIONID = 0
if SESSIONID==0:
    animalName   = 'saja125'
    ephysSession = '2012-03-06_17-16-57'
    behavSession = '20120306a'
    tetrodes = [1,2,3,4,5,6,7,8]
elif SESSIONID==1:
    pass

execfile('freqtuning_multiunit.py')
