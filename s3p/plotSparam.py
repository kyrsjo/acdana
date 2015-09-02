import numpy as np
import matplotlib.pyplot as plt
import re
import sys

jobs = []

for arg in sys.argv[1:]:
    jobs.append(arg)
if len(jobs) == 0:
    jobs.append("./")
print jobs

SparamFile = open("SParameter.out",'r')

print SparamFile.readline(),
print SparamFile.readline(),
print SparamFile.readline(),
headerLine = SparamFile.readline()
print headerLine,
print

#TODO: Number of reads also depends on this...
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
    
# for i in xrange(len(freq)):
#     print freq[i]
#     print S[i]
#     print Sabs[i]
#     print Sang[i]
#     print abs(np.linalg.det(S[i])) # S should be unitary, so this should be 1.0
#     print

S    = np.asarray(S)
Sabs = np.asarray(Sabs)
Sang = np.asarray(Sang)
Sdet = np.asarray(Sdet)

plt.figure()
for i in xrange(numModes):
    for j in xrange(i,numModes):
        plt.plot(freq,Sabs[:,i,j], label="$S_{"+str(i)+","+str(j)+"}$")
plt.legend()
plt.xlabel("|S|")
plt.ylabel("Frequency [Hz]")

plt.figure()
for i in xrange(numModes):
    for j in xrange(i,numModes):
        plt.plot(freq,Sang[:,i,j], label="$S_{"+str(i)+","+str(j)+"}$")
plt.legend()
plt.xlabel("angle(S) [radians]")
plt.ylabel("Frequency [Hz]")

plt.figure()
plt.plot(freq,np.abs(Sdet))
plt.axhline(1.0,ls="--", color="r")
plt.ylabel("|det(S)|")
plt.xlabel("Frequency [Hz]")


plt.show()
