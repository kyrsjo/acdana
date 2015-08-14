#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

from Wakefile import *

print 'Usage:'
print 'plotWake.py wakefile1(--name1(--cutlen1)) wakefile2(--name2(--cutlen2)) etc.'

wakes = []
imps  = []
envs  = []
names = []
maxS  = []

for arg in sys.argv[1:]:
    args = arg.split('--')
    if len(args) == 1:
        wakes.append( loadWakeFile(arg) )
        names.append( arg )
        maxS.append(None)
    elif len(args) == 2:
        wakes.append( loadWakeFile(args[0]) )
        names.append( args[1] )
        maxS.append(None)
    elif len(args) == 3:
        wakes.append( loadWakeFile(args[0]) )
        names.append( args[1] )
        maxS.append(float(args[2]))
        for w in wakes[-1]:
            w.cropToS(maxS[-1])
    else:
        print "ERROR"
        exit
    imp = []
    env = []
    for w in wakes[-1]:
        imp.append(ImpedanceSpectrum(w))
        env.append(Envelope(w))
    imps.append(imp)
    envs.append(env)

if len(wakes) == 0:
    print "Please provide at least one wakefile"
    exit (1)

for (w,i,e,n) in zip(wakes,imps,envs,names):
    print "Plotting '"+n+"'"
    plt.figure(1)
    wl = plt.plot(w[0].s, w[0].V, label=n)
    
    wlc = wl[0].get_color()

    plt.figure(2)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.real(i[0].Z[:i[0].goodIdx]), color=wlc, ls='-', label=n)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.imag(i[0].Z[:i[0].goodIdx]), color=wlc, ls='--')
    
    
    plt.figure(3)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.abs(i[0].Z[:i[0].goodIdx]), label=n)
    #print i[0].f[:i[0].goodIdx]
    
    plt.figure(4)
    plt.semilogy(e[0].s,e[0].Venv,label=n)

plt.figure(1)
plt.legend()
plt.xlabel("s [m]")
plt.ylabel("V [V/pC]")

plt.figure(2)
plt.legend()
plt.xlabel("f [GHz]")
plt.ylabel("Z [$\Omega$]")


plt.figure(3)
plt.legend()
plt.xlabel("f [GHz]")
plt.ylabel("|Z| [$\Omega$]")

plt.figure(4)
plt.legend()
plt.xlabel("s [m]")
plt.ylabel("V [V/pC]")

plt.show()
