#!/usr/bin/env python

import numpy as np
import scipy.integrate as sciInt
import scipy.signal as sciSig
import scipy.interpolate as sciIp
import scipy.optimize as sciOpt

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True

##########################################################
### COMMENT THIS OUT IF NOT MAKING PUBLICATION PLOTS! ####
# golden_mean = (np.sqrt(5)-1.0)/2.0      # Aesthetic ratio
# fig_width = 138.2795/72.27*1.5#3.25  # width in inches 
# fig_height = fig_width*golden_mean#*0.89      # height in inches
# fig_size =  [fig_width,fig_height]
# params = {'backend': 'pdf',
#           'axes.labelsize': 10,
#           'text.fontsize': 10,
#           'legend.fontsize': 10,
#           'xtick.labelsize': 8,
#           'ytick.labelsize': 8,
#           'text.usetex': True,
#           'figure.figsize': fig_size,
#           'font.family' : 'serif',
#           'font.serif' : ['Times']}
# import matplotlib
# matplotlib.rcParams.update(params)
##########################################################
#Thesis:
from matplotlib import rc
matplotlib.rc('font',**{'family':'serif','serif':['Times'],'size':10})

GOLDEN = (1.0+np.sqrt(5.0))/2.0

paperWidth  = 17.0-(2.2+2.0) #cm
paperHeigth = paperWidth/GOLDEN

FIGSIZE = (paperWidth*0.393700787,(paperHeigth-2.0)*0.393700787);
scaleFig = 1
FIGSIZE = map(lambda x: x*scaleFig, FIGSIZE)
print "FIGSIZE =", FIGSIZE
DPI = 300
##########################################################

import WakefieldReader as WR

import sys
import os
import re

#Read arguments and files
def USAGE():
    print "Usage: acdana_wakefieldFFTplot.py {LONG(-GdfidL) | TRANS(-GdfidL) | TRANSSYM(-GdfidL)} {t3p_results/OUTPUT/wakefield.out | auto | file1--title,file2--title,...} len_structure[m] beam_offset[mm]"
    print "len_structure and beam_offset is used for normalizing the wakes,"
    print "use beam_offset only in the case of transverse wakes."
    print "Set them to -1 to switch off normalization"
    print ""
    print "Alternative usage:"
    print "acdana_wakefieldFFTplot.py PLOTFILE(-GdfidL) plotfile--title plotfile--title ..." 
    print "This runMode plotfiles to specify file paths and normalizations"
    #exit(1)
    raise

if len(sys.argv) == 1:
    USAGE()

runMode = ""
if   sys.argv[1] == "LONG":
    runMode = "LONG"
    len_structure = float(sys.argv[3])
    gdfidl = None
    transsym = False
elif sys.argv[1] == "LONG-GdfidL":
    runMode = "LONG"
    len_structure = float(sys.argv[3])
    gdfidl = WR.GdfidlData()
    transsym = False
elif sys.argv[1] == "TRANS":
    runMode = "TRANS"
    len_structure = float(sys.argv[3])
    beam_offset = float(sys.argv[4])
    gdfidl = None
    transsym = False
elif sys.argv[1] == "TRANS-GdfidL":
    runMode = "TRANS"
    len_structure = float(sys.argv[3])
    beam_offset = float(sys.argv[4])
    gdfidl = WR.GdfidlData()
    transsym = False
elif sys.argv[1] == "TRANSSYM":
    runMode = "TRANS"
    len_structure = float(sys.argv[3])
    beam_offset = float(sys.argv[4])
    gdfidl = None
    transsym = True
elif sys.argv[1] == "TRANSSYM-GdfidL":
    runMode = "TRANS"
    len_structure = float(sys.argv[3])
    beam_offset = float(sys.argv[4])
    gdfidl = WR.GdfidlData()
    transsym = True
elif sys.argv[1] == "PLOTFILE":
    runMode = "PLOTFILE"
    gdfidl = None
elif sys.argv[1] == "PLOTFILE-GdfidL":
    runMode = "PLOTFILE"
    gdfidl = WR.GdfidlData()
else:
    print "Expected runMode to be one of  'LONG', 'TRANS', 'TRANSSYM' or 'PLOTFILE'"
    USAGE()

if runMode != "PLOTFILE":
    try:
        fname = sys.argv[2]
        if fname == "auto":
            fname = "t3p_results/OUTPUT/wakefield.out"
            data = [WR.WakefieldReader(fname),]
            if runMode == "TRANS":
                data[-1].offset_x = beam_offset
            data[-1].offset_s = 0.0
            data[-1].len_structure = len_structure
            data[-1].extraNorm = 1.0
            data[-1].normalizeShift()
            data[-1].transsym = transsym
            data[-1].fitParams = None

        elif "," in fname:
            fname = fname.split(",")
            data = []
            for f in fname:
                if f == "":
                    continue
                (name,title) = f.split("--")
                data.append(WR.WakefieldReader(name))
                data[-1].title = title
                data[-1].offset_x = beam_offset
                data[-1].offset_s = 0.0
                data[-1].len_structure = len_structure
                data[-1].extraNorm = 1.0
                data[-1].normalizeShift()
                data[-1].transsym = transsym
                data[-1].fitParams = None

        else:
            data = [WR.WakefieldReader(fname),]
            if runMode == "TRANS":
                data[-1].offset_x = beam_offset
            data[-1].offset_s = 0.0
            data[-1].len_structure = len_structure
            data[-1].extraNorm = 1.0
            data[-1].normalizeShift()
            data[-1].transsym = transsym
            data[-1].fitParams = None
    except:
        print "problem parsing fname"
        print
        USAGE()
elif runMode == "PLOTFILE":

    fname = sys.argv[2:]
    data = []
    runMode = None
    for f in fname:
        if f == "":
            continue
        (name,title) = f.split("--")
        try:
            print "Reading '" + name + "' name..."
            plotfile = open(name,'r')
            plotfileLines = plotfile.readlines()
            plotfile.close()
            
            fmode = re.match("type\s*=\s*(LONG|TRANSSYM|TRANS)",plotfileLines[0]).group(1)
            print "\t fmode = '" + fmode + "'"
            if fmode == "TRANSSYM":
                fmode = "TRANS"
                transsym = True
            else:
                transsym = False
            if not runMode:
                runMode = fmode
            else:
                assert fmode == runMode

            wakeFileNameRel = re.match("data\s*=\s*'(\S+)'",\
                                           plotfileLines[1]).group(1)
            print "\t wakeFileNameRel = '" + wakeFileNameRel + "'"
            wakeFileName = os.path.join(os.path.dirname(name),wakeFileNameRel)
            print "\t wakeFileName = '" + wakeFileName + "'"

            len_structure = float(re.match("len_structure\s*=\s*(\d+.\d*)",\
                                               plotfileLines[2]).group(1))
            print "\t len_structure = '" + str(len_structure) + "' [m]"
            
            offset_x = float(re.match("offset_x\s*=\s*([+-]?\d+.\d*)",\
                                          plotfileLines[3]).group(1))
            print "\t offset_x = '" + str(offset_x) + "' [mm]"
            
            offset_s = float(re.match("offset_s\s*=\s*([+-]?\d+.\d*)",\
                                          plotfileLines[4]).group(1))
            print "\t offset_s = '" + str(offset_s) + "' [m]"

            extraNorm = float(re.match("extraNorm\s*=\s*([+-]?\d+.\d*)",\
                                           plotfileLines[5]).group(1))
            print "\t extraNorm = '" + str(extraNorm) + "' [m]"
 
            try:
                cut_s = float(re.match("cut_s\s*=\s*([+-]?\d+.\d*)",\
                                          plotfileLines[6]).group(1))
            except:
                cut_s = None
            print "\t cut_s = '" + str(cut_s) + "'"
            
            DMS = "\d+.?\d*" #Decimal match string
            fitParams_thisfile = []
            for lineNum in xrange(7,len(plotfileLines)):
                ZfitMatch = re.match("\s*Zfit\s*=\s*\[\s*("+DMS+"),\s*(" + DMS + ")\s*\]\s*\[\s*("+DMS+"),\s*("+DMS+"),\s*("+DMS+")\s*\]\s*",plotfileLines[lineNum])
                if ZfitMatch != None:
                    ZfitMatch = ZfitMatch.groups()
                    ZfitMatch = map(float, ZfitMatch)
                    fitParams_thisfile.append([[ZfitMatch[0], ZfitMatch[1]], [ZfitMatch[2], ZfitMatch[3], ZfitMatch[4]]])
                    print "Got fitParams range=[" + str(ZfitMatch[0]) + ", " + str(ZfitMatch[1]) + "],"
                    print "\t initial parameters f0 = " + str(ZfitMatch[2]) + " m^-1, Q = " + \
                        str(ZfitMatch[3]) + ", A = " + str(ZfitMatch[4]) + " V/pC/mm/m"
                    
            data.append(WR.WakefieldReader(wakeFileName, cut_s))
            data[-1].transsym = transsym
            data[-1].len_structure = len_structure
            data[-1].offset_x = offset_x
            data[-1].offset_s = offset_s
            data[-1].extraNorm = extraNorm
            data[-1].title = title
            data[-1].normalizeShift()
            data[-1].fitParams = fitParams_thisfile

        except:
            print
            print "**************************************************"
            print 
            print "Tried to parse plotfile '%s'" %(name,)
            print "Expected syntax: (line numbers matter!)"
            print "type = {LONG|TRANS|TRANSSYM}"
            print "data = 'filename'"
            print "len_structure = float (-1 to disable normalizing)"
            print "offset_x = float (Only used in case of TRANS, set to -1 to disable normalizing) [m]"
            print "offset_s = float (How much to shift the s vector backwards, set to Nsigma*sigma [m]"
            print "extraNorm = float (Extra normalization factor, usefull if wrong in output. Set to 1.0 to disable.)"
            print "cut_s = float (optional) (don't use parts of the wake with s > this value NOT IMPLEMENTED)"
            print "The data filename should be relative to the folder which the plotfile is located in."
            USAGE()
            

##################################################################################
if runMode == "LONG":
    #Prepare data for plotting
    lenNormalized = False
    for d in data:
        if d.len_structure != -1.0:
            if not lenNormalized:
                assert d == data[0], "First wasn't normalized"
                lenNormalized = True
        else:
            assert lenNormalized == False, "Got non-normalized after normalized"

        for wl in d.wakefieldsLongitudinal:
            if d.title == None:
                #wl.plotLab = str(wl.x) +", " + str(wl.y)
                wl.plotLab = "(%e, %e)" % (wl.x, wl.y)
            else:
                #wl.plotLab = d.title + ": (" + str(wl.x) +", " + str(wl.y) + ")"
                wl.plotLab = "%s : (%.4g, %.4g) [mm]" % (d.title, wl.x*1e3, wl.y*1e3)

            wl.makeFFT()
        
            #Integrate total charge
            #print "Q_total = %g [C], %s" % (sciInt.simps(wl.I, x=wl.s), wl.plotLab)

    #Plotting
    plt.figure(1,figsize=FIGSIZE,dpi=DPI)
    for d in data:
        for wl in d.wakefieldsLongitudinal:
            plt.plot(wl.s,wl.W, label=wl.plotLab)
    if gdfidl:
        plt.plot(gdfidl.Wz_s, gdfidl.Wz, label="GdfidL")
    plt.legend(loc=0)
    plt.xlabel("s [m]")
    if lenNormalized:
        plt.ylabel("$W_z$ [V/pC/m]")
    else:
        plt.ylabel("$W_z$ [V/pC]")
#    plt.savefig("wakeplot_Wz-vs-S.png")

    plt.figure(2,figsize=FIGSIZE,dpi=DPI)
    for d in data:
        for wl in d.wakefieldsLongitudinal:
            plt.plot(wl.s,wl.I,label=wl.plotLab)
    #plt.plot(wl.s, 1e-12*np.exp(-(wl.s-2.5e-3*5)**2/(2*2.5e-3**2))/(2.5e-3*np.sqrt(2*np.pi)))
    plt.legend(loc=0)
    plt.xlabel("s [m]")
    plt.ylabel("Bunch current [C/m]")
#    plt.savefig("wakeplot_I-vs-S.png")
    # analytic_I = 1e-12*np.exp(-(wl.s-2.5e-3*5)**2/(2*2.5e-3**2))/(2.5e-3*np.sqrt(2*np.pi))
    # print "Analytic Qtot = ", sciInt.simps(analytic_I,x=wl.s)
    # print "Analytic \int |I|^2 dx =", sciInt.simps(analytic_I**2,x=wl.s)

    plt.figure(3,figsize=FIGSIZE,dpi=DPI)
    for d in data:
        for wl in d.wakefieldsLongitudinal:
            plt.plot(wl.s,np.abs(sciSig.hilbert(wl.W)), label=wl.plotLab)
            #plt.plot(wl.s,abs(wl.W), label=wl.plotLab)
    plt.legend(loc=0)
    plt.xlabel("s [m]")
    if lenNormalized:
        plt.ylabel("$|\text{hilbert}|(W_z)$ [V/pC/m]")
    else:
        plt.ylabel("$|\text{hilbert}|(W_z)$ [V/pC]")
#    plt.savefig("wakeplot_envWz-vs-S.png")

    plt.figure(4,figsize=FIGSIZE,dpi=DPI)
    for d in data:
        for wl in d.wakefieldsLongitudinal:
            #plt.loglog(wl.f[0:wl.safeLen], np.abs(wl.W_FFT[0:wl.safeLen]),label=wl.plotLab)
            #plt.semilogy(wl.f[0:wl.safeLen]/1e9, np.abs(wl.W_FFT[0:wl.safeLen]),label=wl.plotLab)
            #Calculate impedance = FFT(Wx) / FFT(I)
            # 1. Integrate total charge for normalization
            #  (so that total charge is 1pC, as in the units for the wake)
            Qtot = sciInt.simps(wl.I, x=wl.s)
            print "Q_total = %g [C], %s" % (Qtot, wl.plotLab)
            wl.I /= Qtot
            # 2. Get the FFT of I
            #  (already done)
            # 3. Find frequency cutoff for I
            #  (don't try to calculate impedance where the current is only noise)
            I0 = np.abs(wl.I_FFT[0])
            cutLen = wl.safeLen
            for i in xrange(wl.safeLen):
                if np.abs(wl.I_FFT[i]) < I0*np.exp(-6):
                    cutLen = i
                    print "Cutting at f=%g [GHz] corresponding to I0*exp(-3), cutLen/safeLen=%g %%, cutLen=%i" % (wl.f[i]/1e9, 100*float(cutLen)/float(wl.safeLen), cutLen)
                    break
            # 4. Calculate the impedance!
            Z = wl.W_FFT[:wl.safeLen]/(3e8*wl.I_FFT[:wl.safeLen])
            
            plt.plot(wl.f[:cutLen]/1e9,np.abs(Z[:cutLen]), label=wl.plotLab)         
    plt.legend(loc=0)
    plt.xlabel("Frequency [GHz]")
    if lenNormalized:
        plt.ylabel("$|Z_z|$ [Ohm/m]")
    else:
        plt.ylabel("$|Z_z|$ [Ohm]")
#    plt.savefig("wakeplot_Z-vs-F_log.png")

    # plt.figure(5,figsize=FIGSIZE,dpi=DPI)
    # for d in data:
    #     for wl in d.wakefieldsLongitudinal:
    #         #plt.loglog(wl.f[0:wl.safeLen], abs(wl.I_FFT[0:wl.safeLen]),label=wl.plotLab)
    #         plt.semilogy(wl.f[0:wl.safeLen]/1e9, abs(wl.I_FFT[0:wl.safeLen]),label=wl.plotLab)
    # plt.legend()
    # plt.xlabel("Frequency [GHz]")
    # plt.ylabel("Bunch current spectrum")
    # plt.savefig("wakeplot_IFFT-vs-F_log.png")
    

    plt.show()

#######################################################################
elif runMode == "TRANS":
    lenNormalized = False
    offsetNormalized = False
    minS = 1000000000.0;

    for d in data:
        #Prepare data for plotting
        assert ((not d.transsym) and len(d.wakefieldsLongitudinal) == 2) \
            or (d.transsym and len(d.wakefieldsLongitudinal) == 1)
        s0 = len(d.wakefieldsLongitudinal[0].s)
        y0 = d.wakefieldsLongitudinal[0].y

        if d.len_structure != -1.0:
            if not lenNormalized:
                assert d == data[0], "First wasn't normalized"
                lenNormalized = True
        else:
            assert lenNormalized == False, "Got non-normalized after normalized"

        for w in d.wakefieldsLongitudinal:
            assert w.y == y0
            assert len(w.s) == s0
        if d.title == None:
            if gdfidl:
                plotLabel = "T3P"
            else:
                plotLabel = ""
        else:
            plotLabel = d.title

        #Sort data in order of increasing x
        wl = d.wakefieldsLongitudinal
        if not d.transsym:
            if wl[0].x > wl[1].x:
                wltmp = wl[1]
                wl[1] = wl[0]
                wl[0] = wltmp
            assert wl[0].x < wl[1].x

        #Calculate the transverse wake
        print "Calculating transverse wake..."
        # W_x(s) = \int_s0^s \frac{\partial W_z(s')}{\partial x} ds'
        if d.transsym:
            dx = wl[0].x
        else:
            dx = wl[1].x-wl[0].x
        s = wl[0].s
        dWzDx = np.zeros_like(s)
        W_x   = np.zeros_like(s)
        for si in xrange(1,s0):
            if d.transsym:
                dWzDx[si] = wl[0].W[si]/dx
            else:
                assert wl[0].s[si] == wl[1].s[si]
                #Just use a normal derivative formula..
                dWzDx[si] = (wl[1].W[si]-wl[0].W[si])/dx

            #Integrate wrt. s
            #Trapezoidal rule (BUG??!)
            #W_x1[si] = 0.5*(s[si]-s[si-1])*(dWzDs1[si]+dWzDs1[si-1])
            
            #Simpson's rule from the beginning (SLOOOOW)
            #W_x[si] = sciInt.simps(dWzDx[:si+1], x=s[:si+1])
            
            #Simpson's rule in a clever way (assuming equal intervals...)
            if si%2==0:
                W_x[si] = (s[si]-s[si-2])*(dWzDx[si-2]+4*dWzDx[si-1]+dWzDx[si])/6.0 + W_x[si-2]
            else:
                W_x[si] = 0.5*(s[si]-s[si-1])*(dWzDx[si-1] + dWzDx[si]) + W_x[si-1]

        #Debug
        # print "ds = ", s[1]-s[0], " len(s)=", len(s)
        # if len(s) == 4670:
        #     W_x1 /= 2


        if d.offset_x != -1.0:
            if not offsetNormalized:
                assert d == data[0], "First wasn't normalized"
                offsetNormalized = True
            print "Offset-normalizing with", d.offset_x
            W_x /= d.offset_x
        else:
            assert offsetNormalized == False, "Got non-normalized after normalized"
        
        if s[0] < minS:
            minS = s[0];

        #RAW DATA
        plt.figure(1,figsize=FIGSIZE,dpi=DPI)
        wxline = plt.plot(s, wl[0].W, label=plotLabel)
        wxcolor = wxline[0].get_color() #Used when plotting more than one line pr data,
                                        # to get same color scheme everywhere
        
        #SIGNAL
        plt.figure(2,figsize=FIGSIZE,dpi=DPI)
        plt.plot(s, W_x, label=plotLabel)
        
        #ENVELOPE (peak jumper)
        # plt.figure(8,figsize=FIGSIZE,dpi=DPI)
        # sPeaks = []
        # wPeaks = []
        # wAbs = np.abs(W_x)
        # for i in xrange(1,len(s)-1):
        #     if wAbs[i-1] < wAbs[i] and wAbs[i] > wAbs[i+1]:
        #         sPeaks.append(s[i])
        #         wPeaks.append(wAbs[i])
        # plt.semilogy(sPeaks,wPeaks,label=plotLabel)

        #Envelope (double peak jumper, difference between peaks)
        plt.figure(9,figsize=FIGSIZE,dpi=DPI)
        sPeaks_up = []
        wPeaks_up = []
        sPeaks_dn = []
        wPeaks_dn = []
        #wAbs = np.abs(W_x)
        for i in xrange(1,len(s)-1):
            if W_x[i-1] < W_x[i] and W_x[i] > W_x[i+1]:
                sPeaks_up.append(s[i])
                wPeaks_up.append(W_x[i])
            elif W_x[i-1] > W_x[i] and W_x[i] < W_x[i+1]:
                sPeaks_dn.append(s[i])
                wPeaks_dn.append(W_x[i])
        wEnvelope = ( np.interp(s, sPeaks_up, wPeaks_up) - np.interp(s, sPeaks_dn, wPeaks_dn) ) / 2.0
        plt.semilogy(s,wEnvelope,label=plotLabel)
        
        #Test plot for double peak jumper
        # plt.figure(10,figsize=FIGSIZE,dpi=DPI)
        # plt.plot(s,W_x)
        # plt.plot(s,np.interp(s, sPeaks_up, wPeaks_up))
        # plt.plot(s,np.interp(s, sPeaks_dn, wPeaks_dn))
        

        #ENVELOPE
        plt.figure(3,figsize=FIGSIZE,dpi=DPI)
        def smoothListGaussian(list,strippedXs=False,degree=200):
            #Function borrowed from
            #http://www.swharden.com/blog/2008-11-17-linear-data-smoothing-in-python/
            print "Gaussian smoothing, degree=", degree
            window=degree*2-1  
            weight=np.array([1.0]*window)  
            weightGauss=[]  
            for i in range(window):  
                i=i-degree+1  
                frac=i/float(window)  
                gauss=1/(np.exp((4*(frac))**2))  
                weightGauss.append(gauss)  
            weight=np.array(weightGauss)*weight  
            smoothed=[0.0]*(len(list)-window)  
            for i in range(len(smoothed)):  
                smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)  
            print "Smoother done."
            return smoothed
        wakeEnvelope = np.abs(sciSig.hilbert(W_x))
        plt.semilogy(s, wakeEnvelope, label=plotLabel, color=wxcolor)
        #smoothDegree = 301;
        #smoothed = smoothListGaussian([0.0]*(smoothDegree)+list(wakeEnvelope)+[0.0]*(smoothDegree-1), degree=smoothDegree)
        #plt.semilogy(s, smoothed, color=wxcolor, ls="--")
        
        sFirstBunch = 0.5e-9*3e8
        idxSfirstBunch = None
        #plt.axvline(x=0.5e-9*3e8, ls='--', color="black")
        for i in xrange(len(s)):
            if s[i] >= sFirstBunch:
                idxSfirstBunch = i
                break
        #smoothedWakeValue = smoothed[idxSfirstBunch]
        print "sFirstBunch = ", sFirstBunch,\
              "\n\t s[idxSfirstBunch]=", s[idxSfirstBunch],\
              "\n\t idxSfirstBunch =", idxSfirstBunch
        print "\t wake value =", wakeEnvelope[idxSfirstBunch]#, \
        #      "\n\t smoothed wake value =", smoothed[idxSfirstBunch], "[V/pC(/m/mm)]"

        # plt.figure(4,figsize=FIGSIZE,dpi=DPI)
        # plt.semilogy(s, np.abs(W_x), label=plotLabel)

        #IMPEDANCE
        plt.figure(5,figsize=FIGSIZE,dpi=DPI)
        # 0. Make FFT of Wx signal
        t = s/3e8
        dt = t[1]-t[0] #Assume uniform s/t-step
        Wx_FFT = 1e12*t[-1]*np.fft.rfft(W_x)/(len(W_x)-1) #1e12 to scale pC^-1 -> C^-1
        f = np.fft.fftfreq(s.size, d=dt)
        if f[len(Wx_FFT)] >= 0.0:
            safeLen = len(Wx_FFT)
        else:
            safeLen = len(Wx_FFT)-1
        #Calculate impedance = FFT(Wx) / FFT(I)
        # 1. Integrate total charge for normalization
        #  (so that total charge is 1pC, as in the units for the wake)
        Qtot = sciInt.simps(d.wakefieldsLongitudinal[0].I, x=s)
        print "Q_total = %g [C]" % (Qtot,)
        d.wakefieldsLongitudinal[0].I /= Qtot
        # 2. Get the FFT of I
        d.wakefieldsLongitudinal[0].makeFFT()
        I_FFT = d.wakefieldsLongitudinal[0].I_FFT
        assert safeLen == d.wakefieldsLongitudinal[0].safeLen
        # 3. Find frequency cutoff for I
        #  (don't try to calculate impedance where the current is only noise)
        I0 = np.abs(I_FFT[0])
        cutLen = safeLen
        for i in xrange(safeLen):
            if np.abs(I_FFT[i]) < I0*np.exp(-3):
                cutLen = i
                print "Cutting at f=%g [GHz] corresponding to I0*exp(-3), cutLen/safeLen=%g %%, cutLen=%i" % (f[i]/1e9, 100*float(cutLen)/float(safeLen), cutLen)
                break
        # 4. Calculate the impedance!
        Z = Wx_FFT[:safeLen]/(3e8*I_FFT[:safeLen])
        #plt.plot(f[:safeLen]/1e9,np.abs(Wx_FFT[:safeLen]), label=plotLabel)
        ZplotLine = plt.plot(f[:cutLen]/1e9,np.abs(Z[:cutLen]), label=plotLabel)
        ZplotLineColor = ZplotLine[0].get_color()
        
        # (real and complex)
        plt.figure(10,figsize=FIGSIZE,dpi=DPI)
        if plotLabel == "":
            plt.plot(f[:cutLen]/1e9,np.real(Z[:cutLen]), label="Real", ls="-")
            plt.plot(f[:cutLen]/1e9,np.imag(Z[:cutLen]), label="Imag", ls="--")
        else:
            plt.plot(f[:cutLen]/1e9,np.real(Z[:cutLen]),
                     label=plotLabel+"-real", ls="-", color=ZplotLineColor)
            plt.plot(f[:cutLen]/1e9,np.imag(Z[:cutLen]),
                     label=plotLabel+"-imag", ls="--", color=ZplotLineColor)

        # Try to fit the Q/f0/amplitude
        print "FITTING Z PEAK of '" + plotLabel + "'"
        plt.figure(6,figsize=FIGSIZE,dpi=DPI)
        
        ks =  np.fft.fftfreq(s.size, d=s[1]-s[0])
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
            print "\t Fit Cut indexes: ", kidxstart,kidxstop, " => ", ks[kidxstart], ks[kidxstop]

            errfunc = lambda p,k, ref : fitfunc(p,k) - ref
            (p_fitted,fit_converged) = sciOpt.leastsq(errfunc, p_start,\
                    args=(ks[kidxstart:kidxstop], np.abs(Z[kidxstart:kidxstop]*3e8/1e12)))
            assert fit_converged
            #print "\t p_fitted =", p_fitted
            print "\t f0 =", p_fitted[0], "1/m =", p_fitted[0]*3e8/1e9, "GHz, Q =", p_fitted[1], ", A =", p_fitted[2], "V/pC/mm/m"
            return (p_fitted, p_start, (kidxstart,kidxstop))

        # SETUP FOR CLIC_G FIRST CELL
        #fitParams = []
        #fitParams.append(makeFit((49.0,63.0), (56.0,10,150)))
        #fitParams.append(makeFit((71.0,74.0), (73.0,30,30)))
        #fitParams.append(makeFit((87.0,95.0), (91.0,10,50)))
        #fitParams.append(makeFit((128.0,135.0), (132.0,10,50)))
        # SETUP FOR CLIC_502 FIRST CELL
        #fitParams = []
        #fitParams.append(makeFit((45.0,58.0), (50.0,10,100)))
        #fitParams.append(makeFit((73,83), (78,10,30)))
        # SETUP FOR CLIC_502 LAST CELL
        #fitParams = []
        #fitParams.append(makeFit((40.0,64.0), (50.0,10,100)))
        #fitParams.append(makeFit((73,85), (78,10,30)))

        # Get from file
        fitParams = []
        freqSum=0.0; freqSum2=0.0
        QSum=0.0; QSum2=0.0
        ASum=0.0; ASum2=0.0
        if d.fitParams:
            for fp in d.fitParams:
                fitParams.append(makeFit(fp[0], fp[1]))
                
                freqSum  += fitParams[-1][0][0]*3e8/1e9
                freqSum2 += (fitParams[-1][0][0]*3e8/1e9)**2
                
                QSum  += fitParams[-1][0][1]
                QSum2 += fitParams[-1][0][1]**2
                
                ASum  += fitParams[-1][0][2]
                ASum2 += fitParams[-1][0][2]**2
        if len(fitParams) > 0:
            print "Freq. average =", freqSum/float(len(fitParams)), " +- ", \
                np.sqrt( freqSum2/float(len(fitParams)) - (freqSum/float(len(fitParams)))**2 )
            print "Q     average =", QSum/float(len(fitParams)), " +- ", \
                np.sqrt( QSum2/float(len(fitParams)) - (QSum/float(len(fitParams)))**2 )
            print "Ampl. average =", ASum/float(len(fitParams)), " +- ", \
                np.sqrt( ASum2/float(len(fitParams)) - (ASum/float(len(fitParams)))**2 )
        
        plt.plot(ks[:cutLen], np.abs(Z[:cutLen])*3e8/1e12, label="Data")
        for fit in fitParams:
            p_fitted = fit[0]
            p_start = fit[1]
            (kidxstart, kidxstop) = fit[2]
            plt.plot(ks[kidxstart:kidxstop], fitfunc(p_fitted,ks[kidxstart:kidxstop]), '+')#, label="fitted")
            ksfittedDense = np.linspace(ks[kidxstart],ks[kidxstop-1],100)
            plt.plot(ksfittedDense, fitfunc(p_fitted, ksfittedDense), label="fitted (dense)")
            #plt.plot(ks[kidxstart:kidxstop], fitfunc(p_start,ks[kidxstart:kidxstop]))
        plt.legend(loc=0)
        plt.xlabel("ks [1/m]")
        plt.ylabel("$Z_k$ [Vm / pC/m/mm]")
        #plt.subplots_adjust(left=0.15, right=0.96, bottom=0.19, top=0.97)
#        plt.savefig("ImpedanceFitting.pdf", transparent=True)
        print "END FITTING Z PEAK"
        
        # #Bunch
        # plt.figure(20,figsize=FIGSIZE,dpi=DPI)
        # plt.plot(s,d.wakefieldsLongitudinal[0].I)
        # Icalculated = np.exp(-s**2/(2*(2e-3)**2))/(2e-3 * np.sqrt(2*np.pi))
        # plt.plot(s,Icalculated,'r--')
        # plt.xlabel("s [m]")
        # plt.ylabel("I/Qtot")
        # Qtot = sciInt.simps(d.wakefieldsLongitudinal[0].I, x=s)
        # #print "Q_norm (file) = %g [C]" % (sciInt.simps(d.wakefieldsLongitudinal[0].I, x=s),)
        # #print "Q_norm (expr) = %g [C]" % (sciInt.simps(Icalculated, x=s),)
        
        # calcIFFT = np.exp(-np.pi**2 * ks**2 * 2 * (2e-3)**2)
        # calcIFFT *= np.abs(I_FFT[0])/calcIFFT[0]
        # calcIFFT2  = np.exp(-ks**2/(2*( 1.0/(np.pi**2*4*(2e-3)**2) ) ))
        # calcIFFT2 *= np.abs(I_FFT[0])/calcIFFT2[0]

        
        # plt.figure(21,figsize=FIGSIZE,dpi=DPI)
        # plt.plot(ks[:safeLen], np.abs(I_FFT[:safeLen]))
        # plt.plot(ks[:safeLen], calcIFFT[:safeLen], 'r--')
        # plt.plot(ks[:safeLen], calcIFFT2[:safeLen], 'k-')
        # plt.axvline(ks[cutLen])
        # plt.axhline(np.abs(I_FFT[0])*np.exp(-3))
        # plt.xlabel("ks [1/m]")
        # plt.ylabel("I_FFT")
        
        # plt.figure(22,figsize=FIGSIZE,dpi=DPI)
        # plt.plot(ks[:safeLen], np.real(I_FFT[:safeLen]), label="Real")
        # plt.plot(ks[:safeLen], np.imag(I_FFT[:safeLen]), label="Imag")
        # plt.legend()
        # plt.xlabel("ks [1/m]")
        # plt.ylabel("I_FFT")

        # plt.figure(23,figsize=FIGSIZE,dpi=DPI)
        # plt.plot(ks,f)
        # plt.figure(24,figsize=FIGSIZE,dpi=DPI)
        # plt.plot(f/ks)

        #NOMINAL WAKE
        plt.figure(7,figsize=FIGSIZE,dpi=DPI)
        W_nom = np.fft.irfft(3e8*Z[:cutLen]*I_FFT[0], len(s))
        print "Using dirac driving pulse with I_FFT[0] =", I_FFT[0]
        plt.plot(s+d.offset_s,W_nom, label=plotLabel+"-reconstructed", color=wxcolor,ls='--')
        plt.plot(s,W_x, color=wxcolor, label=plotLabel + "-data")

        #Plot wake from fitting
        wakeFunc = lambda p, s: -p[2]*np.exp(-2*np.pi*p[0]*s/(2*p[1]))*\
            np.sin(2*np.pi*p[0]*s*np.sqrt(1-(4*p[1])**-2))
        s_anaWake = np.linspace(0,s.max(),5000)
        anaWake = np.zeros_like(s_anaWake)
        for fit in fitParams:
            p_fitted = fit[0]
            anaWake += wakeFunc([p_fitted[0],p_fitted[1], p_fitted[2]],s_anaWake)
        plt.plot(s_anaWake,anaWake, label='From fitted impedance', ls="..")
        
    ## End loop over datafiles
    
    plt.figure(1,figsize=FIGSIZE,dpi=DPI)
    if gdfidl:
        if minS > gdfidl.Wx_s[0]:
            minS = gdfidl.Wx_s[0]
    plt.xlabel("s [m]")
    unitString = "V/pC"
    if lenNormalized:
        unitString += "/m"
    plt.ylabel("$V_z$ ["+unitString+"]")
    plt.legend(loc=0)
    (got_smin,got_smax) = plt.xlim()
    plt.xlim(minS,got_smax)
    plt.subplots_adjust(left=0.15,bottom=0.15,right=0.98,top=0.97)
#    plt.savefig("wakeplot_Vz.png")
    plt.savefig("wakeplot_Vz.pdf")
    plt.xlim(minS,0.5)
#    plt.savefig("wakeplot_Vz_zoom.png")
    plt.savefig("wakeplot_Vz_zoom.pdf")

    plt.figure(2,figsize=FIGSIZE,dpi=DPI)
    if gdfidl:
        gdfidlWakeLine = plt.plot(gdfidl.Wx_s, gdfidl.Wx, 'k', label="GdfidL")
        gdfidlColor = gdfidlWakeLine[0].get_color()
    plt.xlabel("s [m]")
    unitString = "V/pC"
    if lenNormalized:
        unitString += "/m"
    if offsetNormalized:
        unitString += "/mm"
    plt.ylabel("$V_x$ [" + unitString + "]")
    plt.legend(loc=0)
    plt.xlim(minS,got_smax)
    plt.subplots_adjust(left=0.15,bottom=0.15,right=0.98,top=0.98)
#    plt.savefig("wakeplot_Vx.png")
    plt.savefig("wakeplot_Vx.pdf")
    plt.xlim(minS,0.5)
#    plt.savefig("wakeplot_Wx_zoom.png")
    plt.savefig("wakeplot_Vx_zoom.pdf")

    
    plt.figure(3,figsize=FIGSIZE,dpi=DPI)
    if gdfidl:
        iclFac=1.0
        plt.semilogy(gdfidl.Wx_s[:int(len(gdfidl.Wx_s)*iclFac)], np.abs(sciSig.hilbert(gdfidl.Wx))[:int(len(gdfidl.Wx_s)*iclFac)], label="GdfidL", color=gdfidlColor)
    plt.xlabel("s [m]")
    plt.ylabel("$|\mathrm{hilbert}(V_x)|$ [" + unitString + "]")
    plt.legend(loc=0)
    plt.xlim(minS,got_smax)
#    plt.savefig("wakeplot_envHilbertWx-vs-S.png")

    # plt.figure(4,figsize=FIGSIZE,dpi=DPI)
    # if gdfidl:
    #     plt.semilogy(gdfidl.Wx_s, np.abs(gdfidl.Wx), label="GdfidL")
    # plt.xlabel("s [m]")
    # plt.ylabel("$|W_x|$ [" + unitString + "]")
    # plt.legend()
    # plt.xlim(minS,got_smax)
    # plt.savefig("wakeplot_envAbsWx-vs-S.png")
    
    plt.figure(5,figsize=FIGSIZE,dpi=DPI)
    if gdfidl:
        plt.plot(gdfidl.ZxF/1e9, np.abs(gdfidl.Zx), label="GdfidL", color=gdfidlColor)
    unitString = "\Omega\mathrm{"
    if lenNormalized:
        unitString += "/m"
    if offsetNormalized:
        unitString += "/mm"
    unitString += "}"
    plt.xlabel("f [GHz]")
    plt.ylabel("$|Z_x|$ [$"+unitString+"$]")
    plt.legend(loc=0)
    plt.subplots_adjust(left=0.15,bottom=0.15,right=0.98,top=0.96)
    plt.savefig("impedance.pdf")
#    plt.savefig("impedance.png")

    plt.figure(7,figsize=FIGSIZE,dpi=DPI)
    if gdfidl:
        plt.plot(gdfidl.Wx_s, gdfidl.Wx, label="GdfidL", color=gdfidlColor)
    plt.xlabel("s [m]")
    unitString = "V/pC"
    if lenNormalized:
        unitString += "/m"
    if offsetNormalized:
        unitString += "/mm"
    plt.ylabel("$W_x$ [" + unitString + "]")
    plt.legend(loc=0)
    plt.xlim(minS,got_smax)
#    plt.savefig("wakeplot_WxNOM-vs-S.png")

    # plt.figure(8,figsize=FIGSIZE,dpi=DPI)
    # plt.xlabel("s [m]")
    # unitString = "V/pC"
    # if lenNormalized:
    #     unitString += "/m"
    # if offsetNormalized:
    #     unitString += "/mm"
    # plt.ylabel("Envelope of $V_x$ [" + unitString + "]")
    # plt.legend(loc=0)

    plt.figure(9,figsize=FIGSIZE,dpi=DPI)
    if gdfidl:
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
        wEnvelope = ( np.interp(gdfidl.Wx_s, sPeaks_up, wPeaks_up) - 
                      np.interp(gdfidl.Wx_s, sPeaks_dn, wPeaks_dn) ) / 2.0
        plt.semilogy(gdfidl.Wx_s,wEnvelope, label='GdfidL', color=gdfidlColor)
    plt.xlabel("s [m]")
    unitString = "V/pC"
    if lenNormalized:
        unitString += "/m"
    if offsetNormalized:
        unitString += "/mm"
    plt.ylabel("Envelope of $V_x$ [" + unitString + "]")
    plt.xlim(minS,got_smax)
    plt.legend(loc=0)
    plt.subplots_adjust(left=0.15,bottom=0.15,right=0.98,top=0.96)
    plt.savefig("envelope.pdf")
#    plt.savefig("envelope.png")
    plt.xlim(minS,0.5)
    plt.savefig("envelope_zoom.pdf")
#    plt.savefig("envelope_zoom.png")


    plt.figure(10,figsize=FIGSIZE,dpi=DPI)
    if gdfidl:
        plt.plot(gdfidl.ZxF/1e9, np.real(gdfidl.Zx), label="GdfidL-real", color=gdfidlColor)
        plt.plot(gdfidl.ZxF/1e9, np.imag(gdfidl.Zx), label="GdfidL-imag", color=gdfidlColor)
    unitString = "\Omega\mathrm{"
    if lenNormalized:
        unitString += "/m"
    if offsetNormalized:
        unitString += "/mm"
    unitString += "}"
    plt.xlabel("f [GHz]")
    plt.ylabel("$Z_x$ [$"+unitString+"$]")
    plt.legend(loc=0)
    plt.subplots_adjust(left=0.15,bottom=0.15,right=0.98,top=0.96)
    plt.savefig("impedance_ReIm.pdf")
#    plt.savefig("impedance_ReIm.png")
    
    plt.show()
