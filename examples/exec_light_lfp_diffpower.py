
'''
Look at LFP evoked by ChR2 stimulation.
'''


animalName = 'saja125'
#ephysSession =  '2012-01-31_16-39-09' # (L2) Diff power every 10 trials
#trialRange = [range(0,10),range(10,20),range(20,30),range(30,40)]
#ephysSession =  '2012-02-02_17-16-59' # (L7) Diff power every 100 trials (x2)
#trialRange = [range(0,80),range(80,160)]
#ephysSession =  '2012-02-03_18-32-48' # (L2) Diff power every 60 trials (x3)
#trialRange = [range(0,60),range(60,120),range(120,180)]
#ephysSession =  '2012-02-04_13-16-14' # (L2) Diff power every 60 trials (x3)
#trialRange = [range(0,60),range(60,120),range(120,180)]
#ephysSession =  '2012-02-05_17-34-15' # (L7) Diff power every 60 trials (x3)
#trialRange = [range(0,60),range(60,120),range(120,180)]

ephysSession =  '2012-03-06_18-23-51' # (L2) Diff power every 20 trials (x4)
electrodeList = [5,6]
trialRange = [range(0,20),range(20,40),range(40,60),range(60,80)]
rangeLabel = ['1.0','1.2','1.5','2.0']
timeRange = [0.360,0.500] # In sec (light onset at 0.4sec)
baselineRange = [0.398,0.400] # In sec

execfile('light_lfp_diffpower.py')
