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
# vx at nose diffusion box 
#
# NOTE: edit sys.path.append path bellow.
#

import sys
sys.path.append('/home/user/pyhc')
import pyhc
from pylab import *
import dircache

dx=20e3
start=9948e3 
end=9548e3
points = linspace(start,end,int(abs(start-end)/dx)+1)   # Starting point for traces
Bvals = []


def plotBLine(fname):
    g = pyhc.PGrid()
    g.open(fname)
    if not g.isok():
        print 'Error: could not open hcfile'
        sys.exit(0)

    for x in points:
        Bvals.append(g.zintpol(x,0.,0.,'Bz'))
    plot(points,Bvals,'ro-')
    hold(True)
    axvline(9738e3)
    hold(False)
    draw()
    show()


if len(sys.argv) != 2:
    print 'usage: ', sys.argv[0], 'hcfile_prefix'
    sys.exit(-1)
prefix = sys.argv[1]
prelen = len(prefix)

#ion()
files = dircache.listdir('./')
for file in files:
    if file[-3:] == '.hc' and file[:prelen] == prefix:
        Bvals = []
        plotBLine(file)







