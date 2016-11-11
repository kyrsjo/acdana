#!/usr/bin/env python

import netCDF4 as ncdf
import sys
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print "Usage: findBadBC.py filename.ncdf"
    exit(1)
mfile_name = sys.argv[1]

#Get the root group in the file
mfile=ncdf.Dataset(mfile_name,'r')
print mfile

print mfile.dimensions["ncoords"]
#a = len(mfile.dimensions["ncoords"])

dim_tetexterior = len(mfile.dimensions["tetexterior"])
dim_tetexteriorsize = len(mfile.dimensions["tetexteriorsize"])

assert dim_tetexteriorsize == 9

#Find real surface tris that are alone. They should have SIDESET != -1

lonely_tris = [] #Each element is a sorted array of 3 indexes+SIDESET
doubleBC_tris = []
def checkTri(tri):
    tri[0].sort()
    #assert tri[0][0] < tri[0][1] < tri[0][2]
    #print tri

    #Is this TRI already known?
    i_delete = -1
    for (t,i) in zip(lonely_tris,xrange(len(lonely_tris))):
        if (t[0]==tri[0]).all():
            print "MATCH!", t, tri, i
            i_delete = i
            
            #BAD: at least one of them are not a true internal surf!
            if t[1] != -1:
                doubleBC_tris.append(t)
                print "BAD", t
            if tri[1] != -1:
                doubleBC_tris.append(tri)
                print "BAD", tri
                
            break
    if i_delete == -1:
        lonely_tris.append(tri)
    else:
        del lonely_tris[i_delete]
        print len (lonely_tris)
    #print i
    

# Columns in table tetrahedon_exterior (guessed):
# BLOCK corner1 corner2 corner3 corner4 mid1 mid2 mid3 SIDESET
for i in xrange(dim_tetexterior):
    if i%100==0:
        print float(i)/dim_tetexterior*100
    #print mfile.variables["tetrahedron_exterior"][i,:]
    
    tet = mfile.variables["tetrahedron_exterior"][i,:]

    checkTri((np.array( (tet[1],tet[2],tet[3]) ), tet[8]) )
    #checkTri((np.array( (tet[1],tet[2],tet[4]) ), tet[8]) )
    #checkTri((np.array( (tet[1],tet[3],tet[4]) ), tet[8]) )
    #checkTri((np.array( (tet[2],tet[3],tet[4]) ), tet[8]) )
    

X = []
Y = []
Z = []

for t in lonely_tris:
    if t[1] == -1:
        print "LONELY/NO BC:", t
        X.append(mfile.variables["coords"][t[0][0],0])
        Y.append(mfile.variables["coords"][t[0][0],1])
        Z.append(mfile.variables["coords"][t[0][0],2])

        X.append(mfile.variables["coords"][t[0][1],0])
        Y.append(mfile.variables["coords"][t[0][1],1])
        Z.append(mfile.variables["coords"][t[0][1],2])

        X.append(mfile.variables["coords"][t[0][2],0])
        Y.append(mfile.variables["coords"][t[0][2],1])
        Z.append(mfile.variables["coords"][t[0][2],2])
        
fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")

ax.scatter(X,Y,Z)
plt.show()
