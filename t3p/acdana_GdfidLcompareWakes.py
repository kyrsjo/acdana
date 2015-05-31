#!/usr/bin/env python

import numpy as np
import scipy.integrate as sciInt
import scipy.signal as sciSig
import matplotlib.pyplot as plt

import scipy.optimize as sciOpt

import WakefieldReader as WR
import sys


gdfidl = WR.GdfidlData(sOffset=0.0, WxScale=None)

#Compare Wz and Wx
# peakGdfidL_z = max(gdfidl.Wz)
# peakGdfidL_x = max(gdfidl.Wx)
# scale = peakGdfidL_z/peakGdfidL_x
# print "Scaling: %g" % (scale)
# gdfidl.Wx *= scale

# plt.plot(gdfidl.Wz_s,gdfidl.Wz, label="Wz")
# plt.plot(gdfidl.Wx_s,gdfidl.Wx, label="Wx")
# plt.legend()
# plt.show()

#Plot Wake and current
FAC = 7.12227E+01
print "Q_total = %g [C], from Wx, FAC=%g" % (sciInt.simps(gdfidl.Wx_I/FAC, x=gdfidl.Wx_I_s),FAC)
FAC= 1.30091E+01
print "Q_total = %g [C], from Wz, FAC=%g" % (sciInt.simps(gdfidl.Wz_I/FAC, x=gdfidl.Wz_I_s),FAC)

plt.figure(1)
plt.plot(gdfidl.Wx_s, gdfidl.Wx)
plt.plot(gdfidl.Wx_I_s, gdfidl.Wx_I)
plt.xlabel("s [m]")
plt.ylabel("$W_x$ [V/pC/m/mm]")

plt.figure(2)
#plt.plot(gdfidl.ZzF,np.abs(gdfidl.Zz))
plt.plot(gdfidl.ZxF/1e9,np.abs(gdfidl.Zx))
plt.xlabel("Frequency [GHz]")
plt.ylabel("$|Z_x|$ [$\Omega/m/mm$]")

plt.figure(3)
plt.plot(gdfidl.ZxF/1e9,np.real(gdfidl.Zx), label="Real")
plt.plot(gdfidl.ZxF/1e9,np.imag(gdfidl.Zx), label="Imag")
plt.legend(loc=0)
plt.xlabel("Frequency [GHz]")
plt.ylabel("$|Z_x|$ [$\Omega/m/mm$]")

plt.figure(4)
plt.semilogy(gdfidl.Wx_s, np.abs(gdfidl.Wx))
#plt.plot(gdfidl.Wx_I_s, gdfidl.Wx_I)
plt.xlabel("s [m]")
plt.ylabel("$W_x$ [V/pC/m/mm]")

plt.figure(5)
sPeaks_up = []
wPeaks_up = []
sPeaks_dn = []
wPeaks_dn = []
for i in xrange(1,len(gdfidl.Wx_s)-1):
    if gdfidl.Wx[i-1] < gdfidl.Wx[i] and gdfidl.Wx[i] > gdfidl.Wx[i+1]:
        sPeaks_up.append(gdfidl.Wx_s[i])
        wPeaks_up.append(gdfidl.Wx[i])
    elif gdfidl.Wx[i-1] > gdfidl.Wx[i] and gdfidl.Wx[i] < gdfidl.Wx[i+1]:
        sPeaks_dn.append(gdfidl.Wx_s[i])
        wPeaks_dn.append(gdfidl.Wx[i])
wEnvelope = ( np.interp(gdfidl.Wx_s, sPeaks_up, wPeaks_up) - np.interp(gdfidl.Wx_s, sPeaks_dn, wPeaks_dn) ) / 2.0
plt.semilogy(gdfidl.Wx_s,wEnvelope)
plt.xlabel("s [m]")
plt.ylabel("Envelope of $V_x$ [V/pC/m/mm]")
plt.xlim(0,30)

plt.figure(6)
#s = gdfidl.Wx_s
Z = gdfidl.Zx
#ks =  np.fft.fftfreq(s.size, d=s[1]-s[0])
f = gdfidl.ZxF
ks = f/3e8
fitfunc = lambda p,k : p[2] * np.sqrt( (p[0]**2 * (4*p[1]**2-1)) / \
           (16*np.pi**2 *(k**2*p[0]**2+(k**2-p[0]**2)**2*p[1]**2 ) ) )
        
def makeFit(startStop, p_start):
    print "Fitting..."
    print "\t startStop =", startStop #start and stop frequencies [1/m]
    print "\t p_start =", p_start #f0 [1/m], Q, amplitude [V/pC/mm/m]
    (kstart,kstop) = startStop
    (kidxstart, kidxstop) = (0,-1)
    for kidx in xrange(len(ks)):
        if kstart < ks[kidx] and kidxstart == 0:
            kidxstart = kidx
        elif kstop < ks[kidx] and kidxstop == -1:
            kidxstop = kidx
            break
    print "\t Fit Cut indexes: ", kidxstart,kidxstop, " => ", ks[kidxstart], ks[kidxstop], "1/m =",\
        ks[kidxstart]*3e8/1e9, ks[kidxstop]*3e8/1e9, "GHz"

    errfunc = lambda p,k, ref : fitfunc(p,k) - ref
    (p_fitted,fit_converged) = sciOpt.leastsq(errfunc, p_start,\
             args=(ks[kidxstart:kidxstop], np.abs(Z[kidxstart:kidxstop]*3e8/1e12)))
    assert fit_converged
    #print "\t p_fitted =", p_fitted
    print "\t f0 =", p_fitted[0], "1/m =", p_fitted[0]*3e8/1e9, "GHz, Q =", p_fitted[1], ", A =", p_fitted[2], "V/pC/mm/m"
    return (p_fitted, p_start, (kidxstart,kidxstop))


fitParams = []
fitParams.append(makeFit((55,56.8),(55.9,10,170)))
fitParams.append(makeFit((54,57.6),(55.9,10,170)))
fitParams.append(makeFit((53,58.5),(55.9,10,170)))
fitParams.append(makeFit((52,59.4),(55.9,10,170)))
# fitParams.append(makeFit((52,60),(55.9,10,170)))
# fitParams.append(makeFit((52,60),(55.9,10,170)))


#print len(ks), len(Z)
plt.plot(ks, np.abs(Z*3e8/1e12), label="data")
for (fit,idx) in zip(fitParams, xrange(len(fitParams))):
    p_fitted = fit[0]
    p_start = fit[1]
    (kidxstart, kidxstop) = fit[2]
    # plt.plot(ks[kidxstart:kidxstop], fitfunc(p_fitted,ks[kidxstart:kidxstop]), label="fitted")
    ksfittedDense = np.linspace(ks[kidxstart],ks[kidxstop-1],100)
    #plt.plot(ksfittedDense, fitfunc(p_fitted, ksfittedDense), label="fitted (dense)")
    plt.plot(ksfittedDense, fitfunc(p_fitted, ksfittedDense), label="fit #%i"%(idx,))
    #plt.plot(ks[kidxstart:kidxstop], fitfunc(p_start,ks[kidxstart:kidxstop]))
plt.legend()
plt.xlabel("ks [1/m]")
plt.ylabel("$Z_k$ [Vm / pC/m/mm]")

plt.figure(7)
plt.plot(f/1e9, np.abs(Z), label="data", lw=10.0, color='k',zorder=-1)
for (fit,idx) in zip(fitParams, xrange(len(fitParams))):
    p_fitted_freq = fit[0]
    p_fitted_freq[0]*=3e8
    p_fitted_freq[2]*=1e12
    p_start = fit[1]
    (fidxstart, fidxstop) = fit[2]
    fFittedDense = np.linspace(f[fidxstart],f[fidxstop-1],100)
    #print fFittedDense, fitfunc(p_fitted_freq, fFittedDense)
    plt.plot(fFittedDense/1e9, fitfunc(p_fitted_freq, fFittedDense), label="fit #%i"%(idx,), lw=3.0,zorder=len(fitParams)-idx)
fFittedDense = np.linspace(10e9,20e9,300)
plt.plot(fFittedDense/1e9, fitfunc((16.91e9, 11.1, 125.0e12), fFittedDense), label="1 meter", lw=1.5,zorder=len(fitParams))
plt.legend(loc=0)
plt.xlabel("f [GHz]")
plt.ylabel("$|Z_x|$ [V/pC/m/mm]")


#plt.plot(



plt.show()
