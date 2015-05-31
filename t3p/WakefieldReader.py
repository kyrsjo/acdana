import numpy as np
import re

floatReg = "-?\d+.\d+e[-\+]\d+"

class WakefieldLongitudinal:

    #Position to read wakefield
    x = None
    y = None
    
    #Raw data
    s = None #s [m]
    W = None #W_long(s) [V/pC]
    I = None #I_bunch(s) [C/m]

    #FFT results
    c = 3e8 #[m/s]
    t = None
    dt = None
    f = None
    W_FFT = None
    I_FFT = None
    safeLen = None #Safe length to plot, where we have FFT and positive frequencies
    
    def __init__(self, dataParse, cut_s=None):
        "Initialize from a text string including comment lines"
        self.s = []
        self.W = []
        self.I = []

        print "CUT S = ", cut_s

        lcounter = 0
        for line in dataParse.splitlines():
            if lcounter == 1:
                #Parse x,y
                print line
                mo = re.match("#\((" + floatReg + "),("+floatReg+")\)",line)
                self.x = float(mo.group(1))
                self.y = float(mo.group(2))
            elif lcounter > 2:
                #Get s,W,I data points
                (s,W,I) = map(float, line.split())
                self.s.append(s)
                self.W.append(W)
                self.I.append(I)
                if cut_s != None and self.s[-1] > cut_s:
                    #Trim the end of the arrays
                    self.s = self.s[:-1]
                    self.W = self.W[:-1]
                    self.I = self.I[:-1]
                    break
            lcounter += 1


        self.s = np.array(self.s, 'float')
        self.W = np.array(self.W, 'float')
        self.I = np.array(self.I, 'float')

    def makeFFT(self):
        self.t = self.s/self.c
        dt = self.t[1] - self.t[0] #assume uniform s-step => t-step
        # for i in xrange(1,len(self.t)):
        #     print self.t[i]-self.t[i-1] #Check that assumption
        self.W_FFT = 1e12*self.t[-1]*np.fft.rfft(self.W)/(len(self.W)-1)
        self.I_FFT = self.t[-1]*np.fft.rfft(self.I)/(len(self.I)-1)
        self.f = np.fft.fftfreq(self.s.size, d=dt)

        if self.f[len(self.W_FFT)] >= 0.0:
            self.safeLen = len(self.W_FFT)
        else:
            self.safeLen = len(self.W_FFT)-1

class WakefieldReader:
    title = None
    len_structure = None
    offset_x = None
    offset_s = None
    transsym = None
    extraNorm = None
    cut_s = None

    wakefieldsLongitudinal = None
    
    def __init__(self, dataFile, cut_s=None):
        dataFilePtr = open(dataFile, 'r')
        self.cut_s = cut_s
        print "CUTS =", self.cut_s

        #Loop over x,y positions
        self.wakefieldsLongitudinal = []
        lastWasComment = True
        data = ""
        for line in dataFilePtr.readlines():
            if lastWasComment:
                if line[0] != "#":
                    #End of a header section
                    lastWasComment = False
            else:
                if line[0] == "#":
                    #Start of a new header section
                    self.wakefieldsLongitudinal.append(WakefieldLongitudinal(data, self.cut_s))
                    data = ""
                    lastWasComment = True
            data += line
        #Get the last block
        self.wakefieldsLongitudinal.append(WakefieldLongitudinal(data, self.cut_s))
        dataFilePtr.close()
        print "WakefieldReader read ", len(self.wakefieldsLongitudinal), " longitudinal wakefields"
    def normalizeShift(self):
        print "Normalizing and shifting WakefieldReader with title='" + str(self.title) + "'"
        for wl in self.wakefieldsLongitudinal:
            if self.len_structure != -1.0:
                print "Length-normalizing with len_structure=", self.len_structure," and extraNorm=", self.extraNorm
                wl.W /= self.len_structure*self.extraNorm
            if self.offset_s != 0.0:
                print "Shifting with backwards with s_offset=", self.offset_s
                wl.s -= self.offset_s
                

class GdfidlData:
    #folder = "/home/kyrre/PhD/TD-500GeV-T3P/baseline/CLIC_G_WakesAndCells/"
    folder = "/afs/cern.ch/project/CLICopti/CLIC_G/TimeDomain/CLIC_G_WakesAndCells/"
    basename = "clicg_first_47m_50micr_meshplane_"
    #basename = "clicg_last_47m_50micr_meshplane_"
    WzFilename   = folder + basename + "Wq_AT_XY.0001"
    WxFilename   = folder + basename + "WXq_AT_XY.0001"
    ReZzFilename  = folder + "clicg_first_47m_50micr_meshplane_ReZ_AT_XY.0001"
    ImZzFilename  = folder + "clicg_first_47m_50micr_meshplane_ImZ_AT_XY.0001"
    ReZxFilename = folder + "clicg_first_47m_50micr_meshplane_ReZx_AT_XY.0001"
    ImZxFilename = folder + "clicg_first_47m_50micr_meshplane_ImZx_AT_XY.0001"

    #To put bunch center at the same point as for T3P
    #sOffset = 5*2.5e-3
    #sOffset = 5*1.2e-3
    sOffset = 0.0 #Do this on T3P data instead
    #Beam offset x structure length (50cells + beampipe) x charge symmetry factor
    # (simulated with 1 pC on H-wall)
    #WxScale = 0.5 * (50+2)*8.332e-3 * 2 #[mm x m]
    WxScale = 0.5 * 50*8.332e-3 * 2 #[mm x m]
    ZxScale = WxScale

    Wz   = None
    Wz_s = None
    Wz_I = None
    Wz_I_s = None

    Zz = None
    ZzF = None

    Wx   = None
    Wx_s = None
    Wx_I = None
    Wx_I_s = None

    Zx = None
    ZxF = None

    def __init__(self, sOffset=None, WxScale=None, ZxScale=None):

        #Read wakes
        self.Wz = []
        self.Wz_s = []
        self.Wz_I = []
        self.Wz_I_s = []
        self.readWake(self.WzFilename,self.Wz, self.Wz_s, self.Wz_I, self.Wz_I_s)
        self.Wz = np.asarray(self.Wz)
        self.Wz_s = np.asarray(self.Wz_s)
        self.Wz_I = np.asarray(self.Wz_I)
        self.Wz_I_s = np.asarray(self.Wz_I_s)

        self.Wx = []
        self.Wx_s = []
        self.Wx_I = []
        self.Wx_I_s = []
        self.readWake(self.WxFilename,self.Wx, self.Wx_s, self.Wx_I, self.Wx_I_s)
        self.Wx = np.asarray(self.Wx)
        self.Wx_s = np.asarray(self.Wx_s)
        self.Wx_I = np.asarray(self.Wx_I)
        self.Wx_I_s = np.asarray(self.Wx_I_s)

        if sOffset != None:
            self.sOffset = sOffset
        if sOffset != 0.0:
            print "GdfidL: Shifting s with sOffset %g" %(self.sOffset,)
            self.Wz_s   += self.sOffset
            self.Wz_I_s += self.sOffset
            self.Wx_s   += self.sOffset
            self.Wx_I_s += self.sOffset

        if WxScale != None:
            self.WxScale = WxScale
        print "GdfidL: Scaling with WxScale %g" %(self.WxScale,)
        self.Wx /= self.WxScale

        #Read impedances
        ZzRe = []
        ZzReF = []
        self.readImpedance(self.ReZzFilename, ZzReF, ZzRe)
        ZzIm = []
        ZzImF = []
        self.readImpedance(self.ImZzFilename, ZzImF, ZzIm)

        ZxRe = []
        ZxReF = []
        self.readImpedance(self.ReZxFilename, ZxReF, ZxRe)
        ZxIm = []
        ZxImF = []
        self.readImpedance(self.ImZxFilename, ZxImF, ZxIm)

        def combineComplex(fRe,fIm, ZRe,ZIm):
            assert len(fRe) == len(fIm), "len(fRe)=%g, len(fIm)=%g" % (len(fRe),len(fIm))
            Z = []
            
            for i in xrange(len(fRe)):
                assert fRe[i] == fIm[i]
                Z.append(complex(ZRe[i],ZIm[i]))
                #print fRe[i], complex(ZRe[i],ZIm[i])
            Z = np.asarray(Z)
            return Z

        self.Zz = combineComplex(ZzReF, ZzImF, ZzRe, ZzIm)
        self.ZzF = np.array(ZzReF)
        #del ZzRe, ZzIm, ZzReF, ZzImF

        self.Zx = combineComplex(ZxReF, ZxImF, ZxRe, ZxIm)
        self.ZxF = np.array(ZxReF)
        #del ZxRe, ZxIm, ZxReF, ZxImF

        if ZxScale != None:
            self.ZxScale = ZxScale
        if ZxScale != 1.0:
            print "GdfidL: Scaling with ZxScale %g" %(self.ZxScale,)
            self.Zx /= self.ZxScale


    def readWake(self,fname, W, s, I, I_s):
        print "GdfidlData::readWake() reading '" + fname +"' ..."
        ifile = open(fname,'r')

        dumpOn = False
        dumpOnNext = False
        wakeNotCharge = True
        for line in ifile:
            if dumpOnNext:
                dumpOn = True
                dumpOnNext = False
            elif line == " # BEGIN_DATA\n":
                print "Wake data begins"
                dumpOnNext = True
            elif dumpOn and (line == " # BEGIN_CHARGE\n" or line == " # END_CHARGE\n"):
                print "Data ends"
                dumpOn = False
            elif line == " % linecolor= 4\n":
                print "Charge data begins"
                dumpOnNext = True
                wakeNotCharge = False
            if dumpOn:
                ls = line.split()
                if wakeNotCharge:
                    W.append(float(ls[1]))
                    s.append(float(ls[0]))
                else:
                    I.append(float(ls[1]))
                    I_s.append(float(ls[0]))
        
        ifile.close()

    def readImpedance(self, fname, f,Z):
        print "GdfidlData::readImpedance() reading '" + fname + "' ..."
        ifile = open(fname, 'r')

        dumpOn = False
        dumpOnNext = False
        for line in ifile:
            if dumpOnNext:
                dumpOn = True
                dumpOnNext = False
            elif line == " # start of data\n":
                print "Impedance data begins"
                dumpOnNext = True
            if dumpOn:
                ls = line.split()
                f.append(float(ls[0]))
                Z.append(float(ls[1]))
        ifile.close()
