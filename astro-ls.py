#! /bin/python

import sys
import numpy as np
from astropy.stats import LombScargle
import matplotlib.pyplot as plt

def read_residuals(resfile):
    it = 0
    skipped = 0
    res_ts = []
    with open(resfile) as fin:
        for line in fin.readlines():
            l = line.split()
            if not 'o' in l:
                l = [ float(x) for x in l ]
                res_ts.append({'mjd': l[0], 'x': l[1], 'y': l[3], 'z': l[5]})
            else:
                # print 'Ommiting outlier line: \"{}\"'.format(line.strip())
                skipped+=1
            it+=1
    print 'Read {}, skipped {} of them.'.format(it, skipped)
    return res_ts


resl = read_residuals(sys.argv[1])

t = np.array([ j['mjd'] for j in resl ])
x = np.array([ j['x'] for j in resl ])
y = np.array([ j['y'] for j in resl ])
z = np.array([ j['z'] for j in resl ])
# t as day difference
min_mjd = t[0]
t = t - min_mjd
# random waves
#nin = len(t)
#phi = 0.5 * np.pi
#A   = 0.5
#w   = 1.0
#x = A * np.sin(w*t+phi)

frequency, power = LombScargle(t, x).autopower()
print 'Component x:'
print '\tStarting frequency: {} or a Period of {} days'.format(frequency[0], 1/frequency[0])
print '\tEnding frequency:   {} or a Period of {} days'.format(frequency[len(frequency)-1], 1/frequency[len(frequency)-1])
mval, midx = max([(v,i) for i,v in enumerate(power)])
print '\tMaximum power at frequency {} -- or a period of {} days --(value {})'.format(frequency[midx], 1/frequency[midx], mval)
print '\tDelta Freq          {}'.format(frequency[1]-frequency[0])
print '\tT0:                 {}'.format(t[0])
print '\tTmax:               {}'.format(t[len(t)-1])
print '\tAverage value:      {}'.format(np.average(x))
axis_x = plt.subplot(3, 1, 1)
plt.plot(frequency, power)
plt.title('Lomb-Scargle Periodogram (astropy)')
plt.ylabel('North')

axis_y = plt.subplot(3, 1, 2)
frequency, power = LombScargle(t, y).autopower()
print 'Component y:'
print '\tStarting frequency: {} or a Period of {} days'.format(frequency[0], 1/frequency[0])
print '\tEnding frequency:   {} or a Period of {} days'.format(frequency[len(frequency)-1], 1/frequency[len(frequency)-1])
mval, midx = max([(v,i) for i,v in enumerate(power)])
print '\tMaximum power at frequency {} -- or a period of {} days --(value {})'.format(frequency[midx], 1/frequency[midx], mval)
print '\tDelta Freq          {}'.format(frequency[1]-frequency[0])
print '\tT0:                 {}'.format(t[0])
print '\tTmax:               {}'.format(t[len(t)-1])
print '\tAverage value:      {}'.format(np.average(y))
plt.plot(frequency, power)
plt.ylabel('East')

axis_z = plt.subplot(3, 1, 3)
frequency, power = LombScargle(t, z).autopower()
print 'Component z:'
print '\tStarting frequency: {} or a Period of {} days'.format(frequency[0], 1/frequency[0])
print '\tEnding frequency:   {} or a Period of {} days'.format(frequency[len(frequency)-1], 1/frequency[len(frequency)-1])
mval, midx = max([(v,i) for i,v in enumerate(power)])
print '\tMaximum power at frequency {} -- or a period of {} days --(value {})'.format(frequency[midx], 1/frequency[midx], mval)
print '\tDelta Freq          {}'.format(frequency[1]-frequency[0])
print '\tT0:                 {}'.format(t[0])
print '\tTmax:               {}'.format(t[len(t)-1])
print '\tAverage value:      {}'.format(np.average(z))
plt.plot(frequency, power)
plt.ylabel('Up')

## all done! show
plt.show()
