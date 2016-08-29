#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

from Wakefile import *

print 'Usage:'
print 'plotWake.py wakefile1(--name1(--cutlen1(--scaleFactor1))) wakefile2(--name2(--cutlen2(--scaleFactor2))) (--trans) etc.'

doTrans = False

wakes       = []
imps        = []
envs        = []
wakes_trans = []
imps_trans  = []
envs_trans  = []
names       = []
maxS        = []
scaleFactor = []

for arg in sys.argv[1:]:
    if arg=="--trans":
        doTrans=True
        if not len(wakes)==0:
            print "'--trans' keyword should come first"
            exit(1)
        continue
    args = arg.split('--')
    if len(args) == 1:
        wakes.append( loadWakeFile(arg) )
        names.append( arg )
        maxS.append(None)
        scaleFactor.append(None)
    elif len(args) == 2:
        wakes.append( loadWakeFile(args[0]) )
        names.append( args[1] )
        maxS.append(None)
        scaleFactor.append(None)
    elif len(args) == 3:
        wakes.append( loadWakeFile(args[0]) )
        names.append( args[1] )
        maxS.append(float(args[2]))
        for w in wakes[-1]:
            w.cropToS(maxS[-1])
        scaleFactor.append(None)
    elif len(args) == 4:
        wakes.append( loadWakeFile(args[0]) )
        names.append( args[1] )
        maxS.append(float(args[2]))
        for w in wakes[-1]:
            w.cropToS(maxS[-1])
        scaleFactor.append(float(args[3]))
        for w in wakes[-1]:
            w.scaleV(scaleFactor[-1])
    else:
        print "ERROR: Unexpected number of arguments (", len(args), ") when splitting '"+arg+"' by '--'. PEBKAC?"
        exit(1)
    
    imp = []
    env = []
    for w in wakes[-1]:
        imp.append(ImpedanceSpectrum(w))
        env.append(Envelope(w))
    imps.append(imp)
    envs.append(env)

    if doTrans:
        vx  = []
        imp = []
        env = []

        #Calculate transverse x-component, assuming single wake and antisymmetry around x=0
        # Code copied from old "TRANSSYM" mode
        assert len(wakes[-1])==1
        wl = wakes[-1][0]

        dx = wl.x
        s  = wl.s
        dVzDx = np.zeros_like(s)
        V_x = np.zeros_like(s)
        for si in xrange(1,len(s)):
            dVzDx[si] = wl.V[si]/dx
            if si%2==0:
                V_x[si] = (s[si]-s[si-2])*(dVzDx[si-2]+4*dVzDx[si-1]+dVzDx[si])/6.0 + V_x[si-2]
            else:
                V_x[si] = 0.5*(s[si]-s[si-1])*(dVzDx[si-1] + dVzDx[si]) + V_x[si-1]
        wx = (Wake(s,V_x,wl.I,0.0,wl.y))
        
        wakes_trans.append((wx,))
        imps_trans.append((ImpedanceSpectrum(wx),))
        envs_trans.append((Envelope(wx),))
        
    else:
        wakes_trans.append(None)
        imps_trans.append(None)
        envs_trans.append(None)
        
if len(wakes) == 0:
    print "Please provide at least one wakefile"
    exit (1)

for (w,i,e,n, wt,it,et) in zip(wakes,imps,envs,names, wakes_trans,imps_trans,envs_trans):
    print "Plotting '"+n+"'"

    #Longitudinal wake
    plt.figure(1)
    wl = plt.plot(w[0].s, w[0].V/w[0].Q, label=n, ls='-')
    wlc = wl[0].get_color()
    Iscale = max(w[0].V/w[0].Q)/max(w[0].I)
    plt.plot(w[0].s, w[0].I*Iscale, ls='--',color=wlc)

    #Longitudinal wake spectrum (re/im)
    plt.figure(2)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.real(i[0].Z[:i[0].goodIdx]), color=wlc, ls='-', label=n)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.imag(i[0].Z[:i[0].goodIdx]), color=wlc, ls='--')
    
    #Longitudinal wake spectrum (abs)
    plt.figure(3)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.abs(i[0].Z[:i[0].goodIdx]), label=n,color=wlc,ls="-")
    #print i[0].f[:i[0].goodIdx]

    #Longitudinal wake envelope
    plt.figure(4)
    plt.plot(e[0].s,e[0].Venv/w[0].Q,label=n,color=wlc,ls="-")

    #Longitudinal wake (log y)
    plt.figure(5)
    plt.semilogy(e[0].s,e[0].Venv/w[0].Q,label=n,color=wlc,ls="-")
    
    #Normalized spectrum
    # plt.figure(6)
    # freqRatioN = i[0].f[:i[0].goodIdx] / 11.2455e3 #f/f_rev
    # plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.real(i[0].Z[:i[0].goodIdx]) / freqRatioN, color=wlc, ls='-', label=n)
    # plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.imag(i[0].Z[:i[0].goodIdx] / freqRatioN), color=wlc, ls='--')

    if doTrans:
        #Transverse wake
        plt.figure(11)
        plt.plot(wt[0].s, wt[0].V/wt[0].Q,label=n,ls="-")

        #Transverse wake spectrum (re/im)
        plt.figure(12)
        plt.plot(it[0].f[:i[0].goodIdx]/1e9, np.real(it[0].Z[:i[0].goodIdx]), color=wlc, ls='-', label=n)
        plt.plot(it[0].f[:i[0].goodIdx]/1e9, np.imag(it[0].Z[:i[0].goodIdx]), color=wlc, ls='--')

        #Transverse wake spectrum (abs)
        plt.figure(13)
        plt.plot(it[0].f[:i[0].goodIdx]/1e9, np.abs(it[0].Z[:i[0].goodIdx]), label=n,color=wlc,ls="-")
        #print i[0].f[:i[0].goodIdx]

        #Transverse wake envelope
        plt.figure(14)
        plt.plot(et[0].s,et[0].Venv/wt[0].Q,label=n,color=wlc,ls="-")
        
        #Transverse wake envelope (log y)
        plt.figure(15)
        plt.semilogy(et[0].s,et[0].Venv/wt[0].Q,label=n,color=wlc,ls="-")
        
        #Normalized transverse spectrum
        # plt.figure(6)
        # freqRatioN = it[0].f[:i[0].goodIdx] / 11.2455e3 #f/f_rev
        # plt.plot(it[0].f[:i[0].goodIdx]/1e9, np.real(it[0].Z[:i[0].goodIdx]) / freqRatioN, color=wlc, ls='-', label=n)
        # plt.plot(it[0].f[:i[0].goodIdx]/1e9, np.imag(it[0].Z[:i[0].goodIdx] / freqRatioN), color=wlc, ls='--')

        
plt.figure(1)
plt.legend()
plt.xlabel("s [m]")
plt.ylabel("$V_z$ [V/pC]")

plt.figure(2)
plt.legend(loc=0)
plt.xlabel("f [GHz]")
plt.ylabel("$Z_z$ [$\Omega$]")

plt.figure(3)
plt.legend()
plt.xlabel("f [GHz]")
plt.ylabel("$|Z_z|$ [$\Omega$]")

plt.figure(4)
plt.legend()
plt.xlabel("s [m]")
plt.ylabel("$V_z$ [V/pC]")

plt.figure(5)
plt.legend()
plt.xlabel("s [m]")
plt.ylabel("$V_z$ [V/pC]")


# plt.figure(6)
# plt.legend()
# plt.xlabel("f [GHz]")
# plt.ylabel("Z [$\Omega$]/($f/f_{rev}$)")

if doTrans:
    plt.figure(11)
    plt.legend()
    plt.xlabel("s [m]")
    plt.ylabel("$V_x$ [V/pC]")
    
    plt.figure(12)
    plt.legend()
    plt.xlabel("f [GHz]")
    plt.ylabel("$Z_x$ [$\Omega$]")
    
    plt.figure(13)
    plt.legend()
    plt.xlabel("f [GHz]")
    plt.ylabel("$|Z_x|$ [$\Omega$]")
    
    plt.figure(14)
    plt.legend()
    plt.xlabel("s [m]")
    plt.ylabel("$V_x$ [V/pC]")

    plt.figure(15)
    plt.legend()
    plt.xlabel("s [m]")
    plt.ylabel("$V_x$ [V/pC]")

plt.show()
