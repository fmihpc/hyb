# This file is part of the HYB simulation platform.
#
# Copyright 2014- Finnish Meteorological Institute
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#!/usr/bin/python
#
# Example script for pyhc module. Creates fieldline trace for magnetic field in vtk files.
# Use visit to view the trace files.
#
# NOTE!!!: edit sys.path.append path bellow to the directory containing _pyhc.so.
#

import sys
sys.path.append('/home/user/pyhc')
import pyhc
from pylab import *
#from numpy import *
import dircache


format = 'vtk'     # values: vtk or m (for matlab).
vector = 'B'     # values: rhov, v, B, B0, B1, E, S, K, j, 
ds = 40000       
maxsteps = 400 
planetary_boundary = 2440e3
initialpts = column_stack((linspace(2500e3, 4500e3, 3),zeros(3),zeros(3)))  # Default initial points for tracing 

xmax = [0,0,0]
xmin = [0,0,0]


if vector=='v' or vector=='rhov' or vector=='B' or vector=='j' or vector=='E'  or vector=='S' or vector=='K':
	field_x = vector + 'x'
	field_y = vector + 'y'
	field_z = vector + 'z'
else:
	field_x = vector[0] + 'x' + vector[1]
	field_y = vector[0] + 'y' + vector[1]
	field_z = vector[0] + 'z' + vector[1]


def stop(pnt):
    "Check if tracing should be stopped"
    global ds,maxsteps,xmax,xmin,planetary_boundary,initialpts

    if pnt[0] < xmin[0] or pnt[0] > xmax[0]:
        return True
    if pnt[1] < xmin[1] or pnt[1] > xmax[1]:
        return True
    if pnt[2] < xmin[2] or pnt[2] > xmax[2]:
        return True
    return False
    if planetary_boundary!=0  and (sqrt(dot(pnt,pnt))<=planetary_boundary):
        return True
    return False


def getbox(g):
    "Simulation box boundaries"
    global ds,maxsteps,xmax,xmin,planetary_boundary,initialpts
    xmin[0] = g.xmin0()
    xmin[1] = g.xmin1()
    xmin[2] = g.xmin2()
    xmax[0] = g.xmax0()
    xmax[1] = g.xmax1()
    xmax[2] = g.xmax2()


def fieldtrace(g,bgn):
    "B field line trace with midpoint method"
    global ds,maxsteps,xmax,xmin,planetary_boundary,initialpts

    F  = zeros(3)
    x  = zeros((maxsteps,3))
    x[0,0] = bgn[0]
    x[0,1] = bgn[1]
    x[0,2] = bgn[2]

    n = 0
    for n in range(maxsteps-1):
        F[0] =  g.intpol(x[n,0],x[n,1],x[n,2],field_x)
        F[1] =  g.intpol(x[n,0],x[n,1],x[n,2],field_y)
        F[2] =  g.intpol(x[n,0],x[n,1],x[n,2],field_z)
        
        if dot(F,F) == 0:
            break

        F = F/sqrt(dot(F,F))
        mp = x[n,:] + 0.5*F*ds

        if stop(mp):
            break

        F[0] =  g.intpol(mp[0],mp[1],mp[2],field_x)
        F[1] =  g.intpol(mp[0],mp[1],mp[2],field_y)
        F[2] =  g.intpol(mp[0],mp[1],mp[2],field_z)

        if dot(F,F) == 0:
            break

        F = F/sqrt(dot(F,F))
        x[n+1,0] = x[n,0] + ds*F[0]
        x[n+1,1] = x[n,1] + ds*F[1]
        x[n+1,2] = x[n,2] + ds*F[2]

        if stop(x[n+1,:]):
            break
    x=x[:n]
    return x
        

def traceHC(fname):
    "Trace hcfile"
    global ds,maxsteps,xmax,xmin,planetary_boundary,initialpts

    g = pyhc.PGrid()
    g.open(fname)
    if not g.isok():
        print 'Error: could not open hcfile'
        sys.exit(0)
    
    try:
        outf = open('mag_'+fname+'.'+format,'w')
    except IOError:
        print 'Error: could not open ', 'mag_'+fname+'.'+format
        sys.exit(-1)
    print "Creating ", 'mag_'+fname+'.'+format

    getbox(g)
    corns = array([[xmin[0], xmin[1], xmin[2]]])
    corns = append(corns,[[xmin[0], xmax[1], xmin[2]]],axis=0)
    corns = append(corns,[[xmax[0], xmin[1], xmin[2]]],axis=0)
    corns = append(corns,[[xmax[0], xmax[1], xmin[2]]],axis=0)
    corns = append(corns,[[xmin[0], xmin[1], xmax[2]]],axis=0)
    corns = append(corns,[[xmin[0], xmax[1], xmax[2]]],axis=0)
    corns = append(corns,[[xmax[0], xmin[1], xmax[2]]],axis=0)
    corns = append(corns,[[xmax[0], xmax[1], xmax[2]]],axis=0)

    lines = []
    for r in range(initialpts.shape[0]):
        xs = fieldtrace(g,initialpts[r,:])
        lines.append(xs)

    ds=-ds
    for r in range(initialpts.shape[0]):
        xs =  fieldtrace(g,initialpts[r,:])
        lines.append(xs)

    if format == 'vtk':
        writeVTKLine(outf,corns,lines)
    elif format == 'm':
        writeMfile(outf,corns,lines)


def writeMfile(outf,corns,lines):
    "Makes matlab m file of trace."
    global ds,maxsteps,xmax,xmin,planetary_boundary,starts
#    print >> outf, 'x y z id'  # uncomment this and rename file extension to '.3D' to view with visit
    for c in corns:
        print >> outf, c[0], c[1], c[2], -1
    for i in range(len(lines)):
        for x in lines[i]:
            print >> outf,  x[0], x[1], x[2], i%len(initialpts)  # mod used to get same id for backward and forward trace.


def writeVTKLine(fobj,corners,lines):
    "Write a vtk file containing traces. (See www.vtk.org/VTK/img/file-formats.pdf)"
    global initialpts

    print >> fobj, '# vtk DataFile Version 2.0'
    print >> fobj, 'pyhc magtracer test'
    print >> fobj, 'ASCII'
    print >> fobj, 'DATASET UNSTRUCTURED_GRID'

    cornnum = corners.shape[0]
    ptsnum = cornnum
    for line in lines:
        ptsnum = ptsnum + line.shape[0]
    print >> fobj, 'POINTS ',ptsnum,' double'

    for x in corners:
        print >> fobj, x[0], x[1], x[2]
    for line in lines:
        for x in line:
            print >> fobj, x[0], x[1], x[2]

    cellnum = 2*cornnum
    for line in lines:
        cellnum = cellnum + line.shape[0] + 1

    print >> fobj, 'CELLS', cornnum+len(lines), cellnum
    for i in range(cornnum): # gridcorners
        print >> fobj, 1, i

    n = 8
    for line in lines:
        print >> fobj, line.shape[0],
        for i in range(line.shape[0]):
            print >> fobj, n,
            n = n+1
        print >> fobj

    print >> fobj, 'CELL_TYPES', cornnum+len(lines)
    for i in range(cornnum): # gridcorners
        print >> fobj, 1
    for line in lines:
        print >> fobj, 4

    print >> fobj, 'POINT_DATA', ptsnum
    print >> fobj, 'SCALARS trace double 1'
    print >> fobj, 'LOOKUP_TABLE default'
    for i in range(ptsnum):
        print >> fobj, 1

    print >> fobj, 'SCALARS id double 1'
    print >> fobj, 'LOOKUP_TABLE default'
    for i in range(len(lines)):
        for x in lines[i]:
            print >> fobj, i%len(initialpts)  # mod used to get same id for backward and forward trace.


def loadInitialPoints(fn):
    "Load initial points for tracing. (file format: ascii txt three numbers per line indicating coordinate)"
    global initialpts
    try: 
        initialpts = loadtxt(fn)
    except IOError:
        print 'Error: cannot read initial points file'
        exit(-1)


def usage():
    print 'usage: ', sys.argv[0], ' [-f initial_points_file] hcfile_prefix'
    sys.exit(-1)


#### main program ####

if len(sys.argv) == 4 and sys.argv[1] == '-f':
    loadInitialPoints(sys.argv[2])
    prefix = sys.argv[3]
    prelen = len(prefix)
elif len(sys.argv) == 2:
    prefix = sys.argv[1]
    prelen = len(prefix)
else:
    usage()

files = dircache.listdir('./')
for file in files:
    if file[-3:] == '.hc' and file[:prelen] == prefix:
        traceHC(file)


