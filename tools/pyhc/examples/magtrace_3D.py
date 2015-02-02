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
# Example script for pyhc module. 
# Simple version of fieldline trace for magnetic fields. Creates 3D format files.
# Use visit to view the trace files.
#
# NOTE!!!: edit sys.path.append path bellow.
#

import sys
sys.path.append('/home/user/pyhc')
import pyhc
from pylab import *
import dircache

ds = 40000
maxsteps = 400
planetary_boundary = 2440e3
starts = linspace(4000e3, -4000e3, 30)   # Starting point for traces

xmax = [0,0,0]
xmin = [0,0,0]

def stop(pnt):
    "Check if tracing should be stopped"
    global ds,maxsteps,xmax,xmin,planetary_boundary,starts

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
    global ds,maxsteps,xmax,xmin,planetary_boundary,starts
    xmin[0] = g.xmin0()
    xmin[1] = g.xmin1()
    xmin[2] = g.xmin2()
    xmax[0] = g.xmax0()
    xmax[1] = g.xmax1()
    xmax[2] = g.xmax2()


def fieldtrace(g,bgn):
    "B field line trace with midpoint method"
    global ds,maxsteps,xmax,xmin,planetary_boundary,starts

    F  = zeros(3)
    x  = zeros((maxsteps,3))
    x[0,0] = bgn[0]
    x[0,1] = bgn[1]
    x[0,2] = bgn[2]

    n = 0
    for n in range(maxsteps-1):
        F[0] =  g.intpol(x[n,0],x[n,1],x[n,2],'Bx')
        F[1] =  g.intpol(x[n,0],x[n,1],x[n,2],'By')
        F[2] =  g.intpol(x[n,0],x[n,1],x[n,2],'Bz')
        
        if dot(F,F) == 0:
            break

        F = F/sqrt(dot(F,F))
        mp = x[n,:] + 0.5*F*ds

        if stop(mp):
            break

        F[0] =  g.intpol(mp[0],mp[1],mp[2],'Bx')
        F[1] =  g.intpol(mp[0],mp[1],mp[2],'By')
        F[2] =  g.intpol(mp[0],mp[1],mp[2],'Bz')

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
    "Makes traces and writes trace to 3D file"
    global ds,maxsteps,xmax,xmin,planetary_boundary,starts

    g = pyhc.PGrid()
    g.open(fname)
    if not g.isok():
        print 'Error: could not open hcfile'
        sys.exit(0)
    
    try:
        outf = open('mag_'+fname+'.3D','w')
    except IOError:
        print 'Error: could not open ', 'mag_'+fname+'.3D'
        sys.exit(-1)

    print "Creating ", 'mag_'+fname+'.3D'
    getbox(g)
    print >> outf, 'x y z num'

    print >> outf, xmin[0], xmin[1], xmin[2]
    print >> outf, xmin[0], xmax[1], xmin[2]
    print >> outf, xmax[0], xmin[1], xmin[2]
    print >> outf, xmax[0], xmax[1], xmin[2]
    print >> outf, xmin[0], xmin[1], xmax[2]
    print >> outf, xmin[0], xmax[1], xmax[2]
    print >> outf, xmax[0], xmin[1], xmax[2]
    print >> outf, xmax[0], xmax[1], xmax[2]

    for start in starts:
        xs = fieldtrace(g,[start,0,0])
        for x in xs:
            print >> outf, x[0], x[1], x[2]
    ds=-ds
    for start in starts:
        xs = fieldtrace(g,[start,0,0])
        for x in xs:
            print >> outf, x[0], x[1], x[2]


## main program

if len(sys.argv) != 2:
    print 'usage: ', sys.argv[0], 'hcfile_prefix'
    sys.exit(-1)
prefix = sys.argv[1]
prelen = len(prefix)

files = dircache.listdir('./')
for file in files:
    if file[-3:] == '.hc' and file[:prelen] == prefix:
        traceHC(file)







