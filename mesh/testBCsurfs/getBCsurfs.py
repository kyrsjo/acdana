#!/usr/bin/env

# Small script to understand the ordering of the exterior tets array

import netCDF4 as ncdf
import sys
import numpy as np

if len(sys.argv) != 2:
    print "Usage: getBCsurfs.py filename.ncdf"
    exit(1)
mfile_name = sys.argv[1]

#Get the root group in the file
mfile=ncdf.Dataset(mfile_name,'r')
print mfile

dim_tetexterior = len(mfile.dimensions["tetexterior"])
dim_tetexteriorsize = len(mfile.dimensions["tetexteriorsize"])
assert dim_tetexteriorsize == 9

X = []
Y = []
Z = []

for i in xrange(dim_tetexterior):
    if i%1000==0:
        print "Progress:",float(i)/dim_tetexterior*100,"%"
    
    tet = mfile.variables["tetrahedron_exterior"][i,:]

    tri = []
    tri.append( (tet[1],tet[2],tet[3]) )
    tri.append( (tet[1],tet[2],tet[4]) )
    tri.append( (tet[1],tet[3],tet[4]) )
    tri.append( (tet[2],tet[3],tet[4]) )

    SS = []
    SS.append( tet[5] )
    SS.append( tet[6] )
    SS.append( tet[7] )
    SS.append( tet[8] )

    for i in xrange(len(SS)):
        if SS[i] != -1:
            X.append(mfile.variables["coords"][tri[i][0],0])
            Y.append(mfile.variables["coords"][tri[i][0],1])
            Z.append(mfile.variables["coords"][tri[i][0],2])

            X.append(mfile.variables["coords"][tri[i][1],0])
            Y.append(mfile.variables["coords"][tri[i][1],1])
            Z.append(mfile.variables["coords"][tri[i][1],2])

            X.append(mfile.variables["coords"][tri[i][2],0])
            Y.append(mfile.variables["coords"][tri[i][2],1])
            Z.append(mfile.variables["coords"][tri[i][2],2])

CSVfile = open(mfile_name+"_surf.csv",'w')
CSVfile.write("x,y,z,s\n")
for i in xrange(len(X)):
    CSVfile.write("%f,%f,%f,1\n"%(X[i],Y[i],Z[i]))
CSVfile.close()

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")

ax.scatter(X,Y,Z)
plt.show()
