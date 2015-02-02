#!/bin/bash
# hybridcode particle log plotter (uses gnuplot)
# USAGE: hyblog_old_particles.sh [number of population / all] [number of variable]
# Old name was of this script was: hlog_particles

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
LOGFILENAME=particles.log

# Template gnuplot command file
TEMPCMDFILENAME=tempAIOJEW6ferui.gnuplotcmd

# Check if the file exist
if [ ! -e $LOGFILENAME ]; then
 echo "ERROR: particle log file ($LOGFILENAME) does not exists";
 exit -1;
fi    

# Variables per population
VARS_PER_POP=20

# Number of populations
npops=$(cut -f 3 -d ' ' $LOGFILENAME | sed -n '1p')

# Check structure of the file header
datacol=$((npops + 3))
checkname=$(eval "cut -f 2 -d ' ' $LOGFILENAME | sed -n '${datacol}p'")
if [ "$checkname" != "DATA_COLUMNS" ]; then
 echo "ERROR: bad logfile header structure"
 exit -1
fi  

# Number of header lines
nheader=$((npops+4));

# Number of data columns in the file
ndatacols=$(eval "cut -f 3 -d ' ' $LOGFILENAME | sed -n '${datacol}p'")

# check columns
this_should_be_zero=$(((ndatacols - 1)%npops));
this_should_be_VARS_PER_POP=$(((ndatacols - 1)/npops));
if [ "$this_should_be_zero" != "0" ] || [ "$this_should_be_VARS_PER_POP" != "$VARS_PER_POP" ]; then
 echo "ERROR: bad logfile column structure"
 exit -1
fi

# List particle populations into a string
poplist=""
nn=3
while [ $nn -le $((npops+2)) ]; do
 popname=$(eval "sed -n '${nn}p' $LOGFILENAME | cut -c3-999")
 poplist="$poplist \"($popname)\""
 nn=$((nn+1))
done
poplist="$poplist \"all\""

# Select particle population
if [ -z "$1" ]; then
 PS3='Choose population (name mass charge): '
 eval "select population in $poplist; do break; done;"
 # Check selection
 if [ -z "$population" ]; then
  echo "ERROR: selected population does not exist ($REPLY)"; exit -1;
 fi
 selected_pop=$REPLY
else
 # Population given in a command line
 if [ "$1" != "all" ]; then
  population="command line"
  selected_pop=$1
 else
  population="all"
  selected_pop=$((npops+1))
 fi
fi

# Select variable
if [ -z "$2" ]; then
 PS3='Choose variable: '
 select variable in "Particles" "Macroparticles" "SplitRate" "JoinRate" "avgVx" "avgVy" "avgVy" "avgV" "cutV" "KineticEnergy" "EscapeRate" "EscapeRateMomX" "EscapeRateMomY" "EscapeRateMomZ" "EscapeRateKinEn" "ImpactRate" "ImpactRateMomX" "ImpactRateMomY" "ImpactRateMomZ" "ImpactRateKinEn"; do break; done;
 # Check selection
 if [ -z "$variable" ]; then
  echo "ERROR: selected variable does not exist ($REPLY)"; exit -1;
 fi
 selected_var=$REPLY
else
 # Variable given in a command line
 variable="command line"
 selected_var=$2
fi

if [ "$population" != "all" ]; then
 # Column number
 column=$((1 + (selected_pop-1)*VARS_PER_POP + selected_var))
 # Variable's name
 varname=$(eval "cut -f $column $LOGFILENAME | sed -n '${nheader}p' | cut -c6-999")
else
 nn=1
 while [ $nn -le $npops ]; do
  column[nn]=$((1 + (nn-1)*VARS_PER_POP + selected_var))
  varname[nn]=$(eval "cut -f ${column[nn]} $LOGFILENAME | sed -n '${nheader}p' | cut -c6-999")
  nn=$((nn+1))
 done 
fi

# Construct plotting command for gnuplot
if [ "$population" == "all" ]; then
 plotcmd="plot ";
 nn=1
 while [ $nn -le $npops ]; do
  plotcmd="$plotcmd \"$LOGFILENAME\" using 1:${column[nn]} title \"${varname[nn]}\" w l"
  if [ $nn -le $((npops-1)) ]; then
   plotcmd="$plotcmd, "
  fi
  nn=$((nn+1))
 done
 varname="multiple"
else
 plotcmd="plot \"$LOGFILENAME\" using 1:$column title \"$varname\" w l"
fi

# Plotting
echo
echo "Plotting variable \"$varname\" from $LOGFILENAME"
echo

# Generate template gnuplot command file
echo "set xlabel \"Time [s]\"; set ylabel \"$varname\"; $plotcmd; pause -1;" >$TEMPCMDFILENAME

# Run gnuplot
gnuplot $TEMPCMDFILENAME

# Remove template command file
rm $TEMPCMDFILENAME
