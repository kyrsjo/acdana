#!/usr/bin/env python

import netCDF4 as ncdf
import sys
import numpy as np

#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print "Usage: findBadBC.py filename_in.ncdf filename_out"
    exit(1)
mfile_in_name  = sys.argv[1]
mfile_out_name = sys.argv[2]


#Get the root group in the file
mfile_in=ncdf.Dataset(mfile_in_name,'r')
print "Got the following mesh:"
print mfile_in

#Open the output file
mfile_out = ncdf.Dataset(mfile_out_name, "w", format=mfile_in.data_model)
#Copy unchanged data
for d in mfile_in.dimensions:
    print "Adding dimension '"+str(d)+"' size",len(mfile_in.dimensions[d])
    mfile_out.createDimension(str(d),len(mfile_in.dimensions[d]))

for v in mfile_in.variables:
    print "Copying data table '"+v+"', type =",mfile_in.variables[v].datatype, ", dimensions =",mfile_in.variables[v].dimensions
    print mfile_in.variables[v].datatype
    print mfile_in.variables[v].dimensions
    mfile_out.createVariable(v,mfile_in.variables[v].datatype,mfile_in.variables[v].dimensions)
    if v == "coords":
        #Do some transformation on this one, so don't copy the data just yet...
        continue
    mfile_out.variables[v][:,:] = mfile_in.variables[v][:,:]

###########################################################
#TRANSFORMATION PARAMETERS!!!

ri = 5e-3
r0 = 101e-3/2.0
r1 = 111.0e-3/2.0

y_shift = 4.99e-3

#Exaggarated parmeters for demo
#ri = 10e-3
#r0 = 50e-3/2.0
#r1 = 111.0e-3/2.0
#y_shift = 25e-3

assert y_shift <= (r1-r0)

z0_up = -326.86e-3    # End of beam pipe
z1_up =  -26.86e-3    # Transition to central region

z0_dn = 442.607799e-3 # End of beam pipe
z1_dn = 142.607799e-3 # Transition to central region
###########################################################

#Transformation 1: Narrower at the pipe ends
def scaleR(z,r, z0,z1):
    assert r <= r1*(1+1e-10) # We should not be scaling anything outside of this region 

    #Inner boundary of scaled region
    r_inner = ri
    #Outer boundary of scaled region -- a function of z
    #Linearly interpolate from (r0,z0) to (r1,z1)
    r_outer = (1-(z-z0)/(z1-z0))*r0 + ((z-z0)/(z1-z0))*r1
    
    if r < r_inner:
        R = r
    else:
        part = (r-r_inner)/(r1-ri)
        R = (1-part)*r_inner + part*r_outer
    
    S = R/r
    return S
#Transformation 2: Shift y
def shiftY(z,r):
    if z < z1_up:
        z0 = z0_up
        z1 = z1_up

        #Inner boundary of scaled region
        r_inner = ri
        #Outer boundary of scaled region -- a function of z
        #Linearly interpolate from (r0,z0) to (r1,z1)
        r_outer = (1-(z-z0)/(z1-z0))*r0 + ((z-z0)/(z1-z0))*r1
        
        assert r >= 0
        if r < r_inner:
            yshift_scale = 0
        elif r < r_outer:
            #Linearly interpolate between (r_inner,0) and (r_outer,1.0)
            yshift_scale = (1-(r-r_inner)/(r_outer-r_inner))*0.0 + ((r-r_inner)/(r_outer-r_inner))*1.0
        else:
            yshift_scale = 1.0
            
        #Linearly interpolate between (z0_up,0) and (z1_up,y_shift)
        yshift = (1-(z-z0_up)/(z1_up-z0_up))*0 + (z-z0_up)/(z1_up-z0_up)*y_shift

        return yshift_scale * yshift
    elif z>z1_dn:
        z0 = z0_dn
        z1 = z1_dn
        
        #Inner boundary of scaled region
        r_inner = ri
        #Outer boundary of scaled region -- a function of z
        #Linearly interpolate from (r0,z0) to (r1,z1)
        r_outer = (1-(z-z0)/(z1-z0))*r0 + ((z-z0)/(z1-z0))*r1

        assert r >= 0
        if r < r_inner:
            yshift_scale = 0
        elif r < r_outer:
            #Linearly interpolate between (r_inner,0) and (r_outer,1.0)
            yshift_scale = (1-(r-r_inner)/(r_outer-r_inner))*0.0 + ((r-r_inner)/(r_outer-r_inner))*1.0
        else:
            yshift_scale = 1.0
        
        #Linearly interpolate between (z1_dn,-y_shift) and (z0_dn,0)
        yshift = (1-(z-z1_dn)/(z0_dn-z1_dn))*(-y_shift) + (z-z1_dn)/(z0_dn-z1_dn)*0
        return yshift_scale * yshift
    else:

        #Inner boundary of scaled region
        r_inner = ri
        #Outer boundary of scaled region
        r_outer = r1

        assert r >= 0
        if r < r_inner:
            yshift_scale = 0
        elif r < r_outer:
            #Linearly interpolate between (r_inner,0) and (r_outer,1.0)
            yshift_scale = (1-(r-r_inner)/(r_outer-r_inner))*0.0 + ((r-r_inner)/(r_outer-r_inner))*1.0
        else:
            yshift_scale = 1.0
        
        #Linearly interpolate between (z1_up,y_shift) and (z1_dn,-y_shift)
        yshift = (1-(z-z1_up)/(z1_dn-z1_up))*y_shift + ((z-z1_up)/(z1_dn-z1_up))*(-y_shift)
        return yshift_scale * yshift
#Loop over the coordinates, apply transformation as needed
ncoords = len(mfile_in.dimensions["ncoords"])

coords_in = mfile_in.variables["coords"]
coords_out = mfile_out.variables["coords"]

tenPercent = ncoords / 10
for i in xrange(ncoords):
    if i % tenPercent == 0:
        print (i*100.0/ncoords),"%"
        
    x,y,z = coords_in[i,:]
    #print i, x,y,z
    
    #r = np.sqrt(x**2+y**2)
    #Z_pre[i] = z
    #R_pre[i] = r
    
    if z < z1_up: #In upstream end
        r = np.sqrt(x**2+y**2)
        S = scaleR(z,r, z0_up,z1_up)
        x *= S
        y *= S
    elif z > z1_dn: #In downstream end
        r = np.sqrt(x**2+y**2)
        S = scaleR(z,r, z0_dn,z1_dn)
        x *= S
        y *= S
    
    r = np.sqrt(x**2+y**2)    
    y += shiftY(z,r)
    
    
    #r = np.sqrt(x**2+y**2)
    #Z_post[i] = z
    #R_post[i] = r
    
    coords_out[i,:] = [x,y,z]

mfile_out.close()    

