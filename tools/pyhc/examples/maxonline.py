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
# maximum of magnetic field at nose on different timesteps. 
#
# NOTE: edit sys.path.append path bellow.
#

import sys
sys.path.append('/home/user/pyhc')
import pyhc
from pylab import *
import dircache

dx=6051e3/6.
start=-3.*6051e3
end=3.*6051e3
points = linspace(start,end,int(abs(start-end)/dx))   # Starting point for traces


def getMax(fname):
    g = pyhc.PGrid()
    g.open(fname)
    if not g.isok():
        print 'Error: could not open hcfile'
        sys.exit(0)
    B = array([])
    for x in points:
        B = append(B, g.zintpol(x,0.,0.,'Bz'))
    return max(abs(B))


if len(sys.argv) != 2:
    print 'usage: ', sys.argv[0], 'hcfile_prefix'
    sys.exit(-1)
prefix = sys.argv[1]
prelen = len(prefix)

ion()
maxBs=[]
files = dircache.listdir('./')
for file in files:
    if file[-3:] == '.hc' and file[:prelen] == prefix:
        m=getMax(file)
        maxBs.append(m)
        print m


plot(range(len(maxBs)),maxBs,'ro-')
title('Max B on line depending on timestep')
show()







