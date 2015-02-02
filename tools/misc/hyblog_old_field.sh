#!/bin/bash
# hybridcode field log plotter (uses gnuplot)
# USAGE: hyblog_old_field.sh [number of variable]
# Old name was of this script was: hlog_field

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

# Logfile
LOGFILENAME=field.log

# Template gnuplot command file
TEMPCMDFILENAME=tempAIOJEW6ferui.gnuplotcmd

# Check if the file exist
if [ ! -e $LOGFILENAME ]; then
 echo "ERROR: field log file ($LOGFILENAME) does not exists";
 exit -1;
fi  

# Select variable
if [ -z "$1" ]; then
 PS3='Choose variable: '
 select variable in "time" "avgBx" "avgBy" "avgBz" "avgB" "max(B)" "avg(div(B))" "max(div(B))" "max(dx*div(B)/B)" "energy(B^2/2*mu0)" "cutE" "cutRhoQ" "cutUe"; do break; done;
 # Check selection
 if [ -z "$variable" ]; then
  echo "ERROR: selected variable does not exist ($REPLY)"; exit -1;
 fi
 column=$REPLY
else
 column=$1
fi

# Variable's name
varname=$(cut -f $column $LOGFILENAME | sed q)

# Generate template gnuplot command file
echo "plot \"$LOGFILENAME\" using 1:$column title \"$varname\" w l; pause -1;" >$TEMPCMDFILENAME

# Run gnuplot
gnuplot $TEMPCMDFILENAME

# Remove template command file
rm $TEMPCMDFILENAME
