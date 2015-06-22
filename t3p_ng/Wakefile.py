#Classes for loading a wake file

import numpy as np
import scipy.integrate as sciInt

class Wake:
    s = None # [m]

    V = None # [V/pC]
    I = None # [C/m]

    x = None # [m]
    y = None # [m]

    Q = None # [C]

    def __init__(self,s,V,I,x,y):
        self.s = np.asarray(s)
        self.V = np.asarray(V)
        self.I = np.asarray(I)
        self.x = x
        self.y = y
        assert len(s) == len(V)
        assert len(s) == len(I)
        
        self.Q = sciInt.simps(self.I, x=self.s)*1e12
        print "Q = ", self.Q, "[pC]"

class ImpedanceSpectrum:
    wake = None
    
    V_FFT = None
    I_FFT = None
    
    Z = None

    t = None
    f = None
    
    cutIdx = None
    goodIdx = None
    
    c = 299792458 #[m/s]
    
    def __init__(self,wake):
        self.wake = wake
        
        self.t = self.wake.s/self.c

        self.V_FFT = 1e12*self.t[-1]*np.fft.rfft(self.wake.V)/(len(self.wake.V)-1) #1e12 to scale pC^-1 -> C^-1
        self.I_FFT = self.t[-1]*np.fft.rfft(self.wake.I)/(len(self.wake.I)-1)
        
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

        self.Z = self.V_FFT/(3e8*self.I_FFT)

def loadWakeFile(fname):
    wakes = []
    
    #loadStrings = ""
    s = []
    V = []
    I = []
    x = None
    y = None
    
    wfile = open(fname,'r')
    
    inCommentBlock = True
    commentBlockLineCounter = 0
    for l in wfile:
        if l[0] == '#':
            if not inCommentBlock:
                if commentBlockLineCounter < 0:
                    wakes.append(Wake(s,V,I,x,y))
                    s = []
                    V = []
                    I = []
                    x = None
                    y = None

                inCommentBlock = True
                commentBlockLineCounter = 0
            
            commentBlockLineCounter = commentBlockLineCounter + 1
                
            if commentBlockLineCounter == 2:
                (x,y) = map(float,l[2:-2].split(','))
                print "reading:",x,y
        else:
            if inCommentBlock:
                assert commentBlockLineCounter == 4 or commentBlockLineCounter == 3
                inCommentBlock = False
    
            # Slow and ineff, but whatevs.
            sVI = l.split()
            s.append(float(sVI[0]))
            V.append(float(sVI[1]))
            I.append(float(sVI[2]))
            
    #last one
    wakes.append(Wake(s,V,I,x,y))
    s = []
    V = []
    I = []
    x = None
    y = None
    wfile.close()
    
    return wakes


# import sys
# wakes = loadWakeFile(sys.argv[1])
# print wakes
