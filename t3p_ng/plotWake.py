#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

from Wakefile import *

print 'Usage:'
print 'plotWake.py wakefile1(--name1(--cutlen1(--scaleFactor1))) wakefile2(--name2(--cutlen2(--scaleFactor2))) (--trans{sym|double}) (--window={rect|rectN|tri|welch|hanning|hamming|blackmanHarris|flatTop}) etc.'
print "Note that the windows affect the wakes calculated AFTER it's specification, or until the the next window. The default window is 'rect' (aka 'no window')"

doTrans     = False
doTransSym  = False
doTransDouble = False

wakefiles   = []
imps        = []
envs        = []
wakes_trans = []
imps_trans  = []
envs_trans  = []
names       = []
maxS        = []
scaleFactor = []

window = None

doExport = False

for arg in sys.argv[1:]:
    print
    if arg.startswith("--trans"):
        doTrans=True
        if not len(wakefiles)==0:
            print "'--trans{sym|double}' keyword should come first"
            exit(1)
        if arg == "--transsym":
            doTransSym=True
        elif arg == "--transdouble":
            doTransDouble=True
        else:
            print "Expected either transsym or transdouble"
            exit(1)
        continue
    elif arg.startswith("--window="):
        winString = arg.split("=")[1]
        print "Window type set = '" + winString + "'"
        if winString.startswith("rect"):
            if winString=="rect":
                window=ImpedanceSpectrum.window_rect
            else:
                ImpedanceSpectrum.window_rectN_cut=int(winString[4:])
                window=ImpedanceSpectrum.window_rectN
        elif winString=="tri":
            window=ImpedanceSpectrum.window_tri
        elif winString=="welch":
            window=ImpedanceSpectrum.window_welch
        elif winString=="hanning":
            window=ImpedanceSpectrum.window_hanning
        elif winString=="hamming":
            window=ImpedanceSpectrum.window_hamming
        elif winString=="blackmanHarris":
            window=ImpedanceSpectrum.window_blackmanHarris
        elif winString=="flatTop":
            window=ImpedanceSpectrum.window_flatTop
        else:
            print "Cannot recognize window name"
            exit(1)
        continue
    elif arg.startswith("--export"):
        doExport = True
        continue
    
    args = arg.split('--')
    if len(args) == 1:
        wakefiles.append( WakeFile(arg) )
        names.append( arg )
        maxS.append(None)
        scaleFactor.append(None)
    elif len(args) == 2:
        wakefiles.append( WakeFile(args[0]) )
        names.append( args[1] )
        maxS.append(None)
        scaleFactor.append(None)
    elif len(args) == 3:
        wakefiles.append( WakeFile(args[0]) )
        names.append( args[1] )
        maxS.append(float(args[2]))
        wakefiles[-1].cropToS(maxS[-1])
        scaleFactor.append(None)
    elif len(args) == 4:
        wakefiles.append( WakeFile(args[0]) )
        names.append( args[1] )
        maxS.append(float(args[2]))
        wakefiles[-1].cropToS(maxS[-1])
        scaleFactor.append(float(args[3]))
        for w in wakefiles[-1].wakes:
            w.scaleV(scaleFactor[-1])
    else:
        print "ERROR: Unexpected number of arguments (", len(args), ") when splitting '"+arg+"' by '--'. PEBKAC?"
        exit(1)
    
    imp = []
    env = []
    print wakefiles
    for w in wakefiles[-1].wakes:
        imp.append(ImpedanceSpectrum(w,window))
        env.append(Envelope(w))
    imps.append(imp)
    envs.append(env)
    
    if doTransDouble:
        assert len(wakefiles[-1].wakes) == 2
        #Calculate the "radial" wake from r1 to r2
        w1 = wakefiles[-1].wakes[0]
        w2 = wakefiles[-1].wakes[1]
        assert len(w1.s)==len(w2.s)
        dx=w2.x-w1.x
        dy=w2.y-w1.y
        dr = np.sqrt(dx**2 + dy**2)

        print "Calculating transverse wakes along vector with angle=", np.arctan2(dy,dx)*180.0/np.pi, "degrees, dr=",dr,"[m]"

        dVzDx = np.zeros_like(w1.s)
        V_x = np.zeros_like(w1.s)
        for si in xrange(1,len(w1.s)):
            dVzDx[si] = (w2.V[si]-w1.V[si])/dr
            if si%2==0:
                V_x[si] = (w1.s[si]-w1.s[si-2])*(dVzDx[si-2]+4*dVzDx[si-1]+dVzDx[si])/6.0 + V_x[si-2]
            else:
                V_x[si] = 0.5*(w1.s[si]-w1.s[si-1])*(dVzDx[si-1] + dVzDx[si]) + V_x[si-1]
        wx = (Wake(w1.s,V_x,w1.I,w1.x+dx/2,w2.y+dy/2))
        
        wakes_trans.append((wx,))
        imps_trans.append((ImpedanceSpectrum(wx,window),))
        envs_trans.append((Envelope(wx),))

        
    elif doTransSym:
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
        imps_trans.append((ImpedanceSpectrum(wx,window),))
        envs_trans.append((Envelope(wx),))
        
    else:
        wakes_trans.append(None)
        imps_trans.append(None)
        envs_trans.append(None)
        
if len(wakefiles) == 0:
    print "Please provide at least one wakefile"
    exit (1)

for (w,i,e,n, wt,it,et) in zip(wakefiles,imps,envs,names, wakes_trans,imps_trans,envs_trans):
    print "Plotting '"+n+"'"

    #Longitudinal wake
    plt.figure(1)
    wl = plt.plot(w.s, w.wakes[0].V/w.wakes[0].Q, label=n, ls='-')
    wlc = wl[0].get_color()
    Iscale = max(w.wakes[0].V/w.wakes[0].Q)/max(w.wakes[0].I)
    plt.plot(w.s, w.wakes[0].I*Iscale, ls='--',color=wlc)

    #Longitudinal wake spectrum (re/im)
    plt.figure(2)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.real(i[0].Z[:i[0].goodIdx]), color=wlc, ls='-', label=n)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.imag(i[0].Z[:i[0].goodIdx]), color=wlc, ls='--')
    
    #Longitudinal voltage spectrum (re/im)
    plt.figure(3)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.real(i[0].V_FFT[:i[0].goodIdx]), color=wlc, ls='-', label=n)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.imag(i[0].V_FFT[:i[0].goodIdx]), color=wlc, ls='--')

    #Longitudinal bunch spectrum (re/im)
    plt.figure(4)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.real(i[0].I_FFT[:i[0].goodIdx]), color=wlc, ls='-', label=n)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.imag(i[0].I_FFT[:i[0].goodIdx]), color=wlc, ls='--')

    
    #Longitudinal wake spectrum (abs)
    plt.figure(5)
    plt.plot(i[0].f[:i[0].goodIdx]/1e9, np.abs(i[0].Z[:i[0].goodIdx]), label=n,color=wlc,ls="-")
    #print i[0].f[:i[0].goodIdx]

    #Longitudinal wake envelope
    plt.figure(6)
    plt.plot(e[0].s,e[0].Venv/w.wakes[0].Q,label=n,color=wlc,ls="-")

    #Longitudinal wake (log y)
    plt.figure(7)
    plt.semilogy(e[0].s,e[0].Venv/w.wakes[0].Q,label=n,color=wlc,ls="-")
    
    #Normalized spectrum
    # plt.figure(8)
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
plt.title("Longitudinal wake")

plt.figure(2)
plt.legend(loc=0)
plt.xlabel("f [GHz]")
plt.ylabel("$Z_z$ [$\Omega$]")
plt.title("Longitudinal wake spectrum ( re (-) / im (--) )")

plt.figure(3)
plt.legend(loc=0)
plt.xlabel("f [GHz]")
#plt.ylabel("$Z_z$ [$\Omega$]")
plt.title("Longitudinal voltage spectrum ( re (-) / im (--)")

plt.figure(4)
plt.legend(loc=0)
plt.xlabel("f [GHz]")
#plt.ylabel("$Z_z$ [$\Omega$]")
plt.title("Longitudinal bunch spectrum ( re (-) / im (--)")

plt.figure(5)
plt.legend()
plt.xlabel("f [GHz]")
plt.ylabel("$|Z_z|$ [$\Omega$]")
plt.title("Longitudinal wake spectrum (abs)")

plt.figure(6)
plt.legend()
plt.xlabel("s [m]")
plt.ylabel("$V_z$ [V/pC]")
plt.title("Longitudinal wake envelope")

plt.figure(7)
plt.legend()
plt.xlabel("s [m]")
plt.ylabel("$V_z$ [V/pC]")
plt.title("Longitudinal wake (log y)")

# plt.figure(8)
# plt.legend()
# plt.xlabel("f [GHz]")
# plt.ylabel("Z [$\Omega$]/($f/f_{rev}$)")
# plt.title("Normalized spectrum")

if doTrans:
    plt.figure(11)
    plt.legend()
    plt.xlabel("s [m]")
    plt.ylabel("$V_x$ [V/pC]")
    plt.title("Transverse wake")
  
    plt.figure(12)
    plt.legend()
    plt.xlabel("f [GHz]")
    plt.ylabel("$Z_x$ [$\Omega$]")
    plt.title("Transverse wake spectrum (re (-) / im (--)")
        
    plt.figure(13)
    plt.legend()
    plt.xlabel("f [GHz]")
    plt.ylabel("$|Z_x|$ [$\Omega$]")
    plt.title("Transverse wake spectrum (abs)")
    
    plt.figure(14)
    plt.legend()
    plt.xlabel("s [m]")
    plt.ylabel("$V_x$ [V/pC]")
    plt.title("Transverse wake envelope")
    
    plt.figure(15)
    plt.legend()
    plt.xlabel("s [m]")
    plt.ylabel("$V_x$ [V/pC]")
    plt.title("Transverse wake envelope (log y)")

if doExport:
    assert len(wakefiles)==1, "Please use only one single wake file for exporting"
    if doTrans:
        ofname = wakefiles[0].fname + "-transwake"
        print "Exporting '"+ ofname + "'"
        ofile = open(ofname,'w')
        ofile.write("# s[m] Wr[V/pC]\n")
        for (s,V) in zip(wt[0].s, wt[0].V/wt[0].Q):
            ofile.write("%f %f\n" % (s,-V))
        ofile.close()
        
plt.show()
