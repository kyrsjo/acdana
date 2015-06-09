#!/usr/bin/env python
# Read a t3p log file, extract the 

import sys

doPlot = False

ts  = []
tst = []

for fn in sys.argv[1:]:
    if fn == "PLOT":
        doPlot = True
        continue
    ifile = open(fn,'r')
    for line in ifile:
        if line[:4] == "+ P ":
            #print line
            ls = line.split()
            ts.append(int(ls[2]))
            tst.append(float(ls[5]))

#Doubles removal
# i = 0
# ts_p = None
# while True:
#     if i >= len(ts):
#         break
#     if ts[i] >= 1700 and ts[i] <= 1800:
#         print i, ts[i], len(ts),ts_p

#     if ts_p == None:
#         ts_p = ts[i]
#     elif ts_p > ts[i]:
#         print "GAP!"
        
#     i += 1

import numpy as np
print "Mean:", np.mean(tst), "steps/sec"

if doPlot:
    import matplotlib.pyplot as plt
    plt.plot(ts,tst)
    plt.show()
