import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import os

jobs = []
titles = []
useDB = False

for arg in sys.argv[1:]:
    
    if arg == "--dB":
        useDB = True
        continue
    
    argsplit = arg.split("--")
    if len(argsplit) == 2:
        jobs  .append(argsplit[0])
        titles.append(argsplit[1])
    else:
        jobs.append(arg)
        titles.append(arg)
if len(jobs) == 0:
    jobs.append("./")
    titles.append("")

class SparamFileData:
    
    numModes = None
    
    freq = None
    S    = None
    Sabs = None
    Sang = None
    Sdet = None
    
    def __init__(self,folder):
        print "Now loading folder '"+folder+"'"
        
        SparamFile = open(os.path.join(folder,"SParameter.out"),'r')

        print SparamFile.readline(),
        print SparamFile.readline(),
        print SparamFile.readline(),
        headerLine = SparamFile.readline()
        print headerLine,
        print

        #TODO: Number of reads (in header) also depends on this...
        numModes = int(np.sqrt( len(headerLine.split())-1 ))
        print "Number of modes:", numModes

        freq  = []
        S     = []
        Sabs  = []
        Sang  = []
        Sdet  = []

        while True:
            line = SparamFile.readline()
            if line == "":
                break
    
            matchObj = re.match("(\S+)"+" \( ?(\S+),  ?(\S+)\)"*numModes**2, line)
            #print matchObj.groups()
            if matchObj == None:
                print "No matchobj?"
                print line,
                exit
    
            ls = line.split()
            freq.append(float(matchObj.group(1)))
    
            S.append   (np.zeros((numModes,numModes),dtype=complex))
            #Sabs.append(np.zeros((numModes,numModes),dtype=float))
            for i in xrange(numModes):
                for j in xrange(numModes):
                    matchIdx = 2+(j+i*numModes)*2
                    real = float(matchObj.group(matchIdx))
                    imag = float(matchObj.group(matchIdx+1))
                    #print i,j, complex(real,imag)
                    S   [-1][i,j] = complex(real,imag)
                    #Sabs[-1][i,j] = np.abs(S[-1][i,j])
            Sabs.append(np.abs(S[-1]))
            Sang.append(np.angle(S[-1]))
            
            Sdet.append(np.linalg.det(S[-1]))
            
            #Store the data
            self.numModes = numModes
            
            self.freq = np.asarray(freq)
            self.S    = np.asarray(S)
            self.Sabs = np.asarray(Sabs)
            self.Sang = np.asarray(Sang)
            self.Sdet = np.asarray(Sdet)

#TODO: The numModes should be the same for all data in order to compare

simData = []
for job,title,idx in zip(jobs,titles,xrange(len(jobs))):
    simData.append(SparamFileData(job))
    
    for i in xrange(simData[-1].numModes):
        for j in xrange(i,simData[-1].numModes):
            plt.figure(1)
            if useDB:
                plt.plot(simData[-1].freq,10*np.log10(simData[-1].Sabs[:,i,j]), label="$S_{"+str(i)+","+str(j)+"}$ (" + title + ")")
            else:
                plt.plot(simData[-1].freq,simData[-1].Sabs[:,i,j], label="$S_{"+str(i)+","+str(j)+"}$ (" + title + ")")
            
            plt.figure(2)
            plt.plot(simData[-1].freq,simData[-1].Sang[:,i,j], label="$S_{"+str(i)+","+str(j)+"}$ (" + title + ")")

            plt.figure(3)
            plt.plot(simData[-1].freq,np.abs(simData[-1].Sdet), label=title)
            
            plt.figure(100+i+j*simData[-1].numModes)
            if useDB:
                plt.plot(simData[-1].freq,np.log10(simData[-1].Sabs[:,i,j]), label=title)
            else:
                plt.plot(simData[-1].freq,simData[-1].Sabs[:,i,j], label=title)

            plt.figure(300+i+j*simData[-1].numModes)
            plt.plot(simData[-1].freq,simData[-1].Sang[:,i,j], label=title)
            
            if idx == len(jobs)-1:
                plt.figure(100+i+j*simData[-1].numModes)
                plt.title("$S_{"+str(i)+","+str(j)+"}$")
                if useDB:
                    plt.ylabel("|S| [dB]")
                else:
                    plt.ylabel("|S|")
                plt.xlabel("Frequency [Hz]")
                plt.legend(loc=0)

                plt.figure(300+i+j*simData[-1].numModes)
                plt.title("$S_{"+str(i)+","+str(j)+"}$")
                plt.ylabel("angle(S) [radians]")
                plt.xlabel("Frequency [Hz]")
                plt.legend(loc=0)
                plt.yticks((-np.pi,-np.pi/2,0,np.pi/2,np.pi),
                           ("$-\pi$","$-\pi/2$",0,"$\pi/2$","$\pi$"))
                
    if idx == len(jobs)-1:
        plt.figure(1)
        plt.legend(loc=0)
        if useDB:
            plt.ylabel("|S| [dB]")
        else:
            plt.ylabel("|S|")
        plt.xlabel("Frequency [Hz]")

        plt.figure(2)
        plt.legend(loc=0)
        plt.ylabel("angle(S) [radians]")
        plt.xlabel("Frequency [Hz]")
        plt.yticks((-np.pi,-np.pi/2,0,np.pi/2,np.pi),
                   ("$-\pi$","$-\pi/2$",0,"$\pi/2$","$\pi$"))
        
        plt.figure(3)        
        plt.axhline(1.0,ls="--", color="r")
        plt.ylabel("|det(S)|")
        plt.xlabel("Frequency [Hz]")
        plt.legend(loc=0)

plt.show()
