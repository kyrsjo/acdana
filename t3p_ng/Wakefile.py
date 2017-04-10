#Classes for loading a wake file

import numpy as np
import scipy.integrate as sciInt

class Wake:
    s = None # [m]

    V = None # [V]
    I = None # [C/m]

    x = None # [m]
    y = None # [m]

    Q = None # [pC]

    def __init__(self,s,V,I,x,y):
        self.s = np.asarray(s)
        self.V = -1.0*np.asarray(V) #To conform to convention W(0+) > 0 -- ACE3P uses inverted definition.
        self.I = np.asarray(I)
        self.x = x
        self.y = y
        assert len(s) == len(V)
        assert len(s) == len(I)
        
        self.Q = sciInt.simps(self.I, x=self.s)*1e12
        print "Q =", self.Q, "[pC]"
        print "maxS =", self.s[-1], "[m]"
        
    def cropToS(self,maxS):
        imax = 0
        print "Wakefield::cropToS(): Cutting to s=",maxS
        for S in self.s:
            imax = imax+1
            if S > maxS:
                break
        if imax < len(self.s):
            print "CUTTING!"
            self.s = self.s[:imax]
            self.V = self.V[:imax]
            self.I = self.I[:imax]
    
    def scaleV(self,scaleFactor):
        self.V *= scaleFactor

class ImpedanceSpectrum:
    wake = None
    
    V_FFT = None # [V/frequency]
    I_FFT = None # [A/frequency]
    
    Z = None # [Ohm]

    t = None
    f = None
    
    cutIdx = None
    goodIdx = None

    windowFunc = None
    
    c = 299792458 #[m/s]

    def __init__(self, wake, windowFunc=None):
        self.wake = wake

        self.t = self.wake.s/self.c
        if windowFunc==None:
            windowFunc = self.window_rect
        self.windowFunc = windowFunc
        window = self.windowFunc(len(self.wake.V))
        
        self.V_FFT = self.t[-1]*np.fft.rfft(self.wake.V*window)/(len(self.wake.V)-1)
        self.I_FFT = self.c * self.t[-1]*np.fft.rfft(self.wake.I*window)/(len(self.wake.I)-1)
        
        self.f = np.fft.fftfreq(self.t.size, d=self.t[1]-self.t[0])

        if not self.f[len(self.V_FFT)] >= 0.0:
            self.cutIdx = len(self.V_FFT)
        else:
            self.cutIdx = len(self.V_FFT)-1
        self.f = self.f[:self.cutIdx]
            
        I0 = np.abs(self.I_FFT[0])
        self.goodIdx = self.cutIdx
        for i in xrange(self.cutIdx):
            if np.abs(self.I_FFT[i])/I0 < np.exp(-6):
                self.goodIdx = i
                print "Cutting at f=%g [GHz] corresponding to I0*exp(-3), cutLen/cutIdx=%g %%, goodLen=%i" % (self.f[i]/1e9, 100*float(self.goodIdx)/float(self.cutIdx), self.goodIdx)
                break

        self.Z = self.V_FFT/self.I_FFT

    # Window functions from
    # https://en.wikipedia.org/wiki/Window_function
    
    @staticmethod
    def window_rect(N):
        return 1.0

    window_rectN_cut = 0
    @staticmethod
    def window_rectN(N):
        w = np.ones(N)
        print N,ImpedanceSpectrum.window_rectN_cut, w
        w[-ImpedanceSpectrum.window_rectN_cut:] = 0.0
        print w
        return w
        
    @staticmethod
    def window_tri(N):
        n = np.arange(N)
        return 1.0 - np.abs(n-(N-2)/2.0) / ((N-1)/2.0)

    @staticmethod
    def window_welch(N):
        n = np.arange(N)
        return 1.0 - (np.abs(n-(N-2)/2.0) / ((N-1)/2.0))**2

    @staticmethod
    def window_hanning(N):
        n = np.arange(N)
        return 0.5*(1.0 - np.cos(2*np.pi*n/(N-1)))

    @staticmethod
    def window_hamming(N):
        n = np.arange(N)
        return (25.0/46.0) - ((46.0-25.0)/46.0)*np.cos(2*np.pi*n/float(N-1))

    @staticmethod
    def window_blackmanHarris(N):
        n = np.arange(N)
        return 0.35875 - 0.48829*np.cos(2*np.pi*n/float(N-1)) + 0.14128*np.cos(4*np.pi*n/float(N-1)) - 0.01168*np.cos(6*np.pi*n/float(N-1))

    @staticmethod
    def window_flatTop(N):
        n = np.arange(N)
        return 1 - 1.93*np.cos(2*np.pi*n/float(N-1)) + 1.29*np.cos(4*np.pi*n/float(N-1)) - 0.338*np.cos(6*np.pi*n/float(N-1)) + 0.028*np.cos(8*np.pi*n/float(N-1))
    
class Envelope:
    wake = None
    s = None
    Venv = None
    
    s_peaksUp = None
    V_peaksUp = None
    s_peaksDn = None
    V_peaksDn = None

    def __init__(self,wake):
        self.wake = wake
        
        self.s = self.wake.s

        self.s_peaksUp = []
        self.V_peaksUp = []
        self.s_peaksDn = []
        self.V_peaksDn = []
        
        #Double peak jumper envelope finder
        for i in xrange(1,len(self.s)-1):
            if self.wake.V[i-1] < self.wake.V[i] and self.wake.V[i] > self.wake.V[i+1]:
                self.s_peaksUp.append(self.s[i])
                self.V_peaksUp.append(self.wake.V[i])
            elif self.wake.V[i-1] > self.wake.V[i] and self.wake.V[i] < self.wake.V[i+1]:
                self.s_peaksDn.append(self.s[i])
                self.V_peaksDn.append(self.wake.V[i])
        self.Venv = ( np.interp(self.s, self.s_peaksUp, self.V_peaksUp) - np.interp(self.s, self.s_peaksDn, self.V_peaksDn) ) / 2.0

class WakeFile:
    s = []
    fname = []
    wakes = []
    def __init__(self,fname):
        print "WakeFile.__init__ reading '"+fname+"'"
        self.fname=fname
        self.wakes = []
        
        #Temps for loading each wake
        s = []
        V = []
        I = []
        x = None
        y = None

        wfile = open(fname,'r')

        inCommentBlock = True
        commentBlockLineCounter = 0
        for l in wfile.xreadlines():
            if l[0] == '#':
                if not inCommentBlock:
                    inCommentBlock = True
                    commentBlockLineCounter = 0

                    if len(s) >= 0:
                        assert len(s)==len(V) and len(s)==len(I)
                        #Add the previous wake to the list
                        self.wakes.append(Wake(s,V,I,x,y))
                        #reset the reader
                        s = []
                        V = []
                        I = []
                        x = None
                        y = None

                commentBlockLineCounter = commentBlockLineCounter + 1

                if commentBlockLineCounter == 2:
                    (x,y) = map(float,l[2:-2].split(','))
                    print "reading point x=%g, y=%g"%(x,y)
            else:
                if inCommentBlock:
                    assert commentBlockLineCounter == 4 or commentBlockLineCounter == 3
                    inCommentBlock = False
                    commentBlockLineCounter=0

                # Slow and ineff, but whatevs.
                sVI = l.split()
                s.append(float(sVI[0]))
                V.append(float(sVI[1]))
                I.append(float(sVI[2]))

        wfile.close()
        
        #Add the last wake to the list
        self.wakes.append(Wake(s,V,I,x,y))

        self.s = self.wakes[0].s
        #Check that s is the same every time
        for w in self.wakes:
            assert np.array_equal(w.s,self.s)
    def cropToS(self,maxS):
        for w in self.wakes:
            w.cropToS(maxS)

        self.s = self.wakes[0].s
        #Check that s is the same every time
        for w in self.wakes:
            assert np.array_equal(w.s,self.s)
