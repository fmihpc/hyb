#!/bin/bash
# hybridcode log pdf plotter (uses gnuplot)
# USAGE: hyblog_old_pdf.sh
# Old name was of this script was: hlog_pdf

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

# Set gnuplot command
if [ "$(hostname)" = "jumbo" -o "$(hostname)" = "sambo" -o "$(hostname)" = "mortti" -o "$(hostname)" = "vertti" ]; then
 GNUPLOT=/usr/local/bin/gnuplot
else
 GNUPLOT=gnuplot
fi

# gnuplot terminal
GNUPLOT_TERMINAL="set terminal postscript landscape color \"Times-Roman\" 6; ";

# Text coordinates
TEXT_RIGHTLIM=0.79
TEXT_TOPLIM=1.05
TEXT_BOTTOMLIM=-0.05

ALLPSFILES=""


# FIELD PLOTTING


# Logfile
LOGFILENAME=field.log

# Check if the file exist
if [ ! -e $LOGFILENAME ]; then
 echo "ERROR: field log file ($LOGFILENAME) does not exists";
 exit -1;
fi  

echo "Generating fieldlog plots..";

OUTFILENAME=log_field;

parse_fieldlog_column()
{
# Template gnuplot command file
TEMPCMDFILENAME=$1

# Select variable
if [ -z "$2" ]; then
 PS3='Choose variable: '
 select variable in "time" "avgBx" "avgBy" "avgBz" "avgB" "max(B)" "avg(div(B))" "max(div(B))" "max(dx*div(B)/B)" "energy(B^2/2*mu0)" "cutE" "cutRhoQ" "cutUe"; do break; done;
 # Check selection
 if [ -z "$variable" ]; then
  echo "ERROR: selected variable does not exist ($REPLY)"; exit -1;
 fi
 column=$REPLY
else
 column=$2
fi

# Variable's name
varname=$(cut -f $column $LOGFILENAME | sed q)

# gnuplot commands
echo "set xlabel \"Time [s]\"; plot \"$LOGFILENAME\" using 1:$column title \"$varname\" w l;" >>$TEMPCMDFILENAME
}

TEMPCOMMANDFILE=tempAIOJEW6ferui.gnuplotcmd;

# Construct gnuplot command file
echo $GNUPLOT_TERMINAL >$TEMPCOMMANDFILE;
echo "set output '$OUTFILENAME.ps'; " >>$TEMPCOMMANDFILE;
echo "set multiplot; set size 0.25,0.25; " >>$TEMPCOMMANDFILE;

echo "set label \"HYBRIDCODE FIELD LOGS\" at screen 0, screen $TEXT_TOPLIM font \"Times-Roman,12\";" >>$TEMPCOMMANDFILE;
echo "set label \"$(pwd)\" at screen 0, screen $TEXT_BOTTOMLIM font \"Times-Roman,14\";" >>$TEMPCOMMANDFILE;
echo "set label \"$(date -R)\" at screen $TEXT_RIGHTLIM, screen $TEXT_BOTTOMLIM font \"Times-Roman,12\";" >>$TEMPCOMMANDFILE;

echo "set origin 0,0.75; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 2;
echo "set origin 0,0.5; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 3;
echo "set origin 0,0.25; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 4;
echo "set origin 0,0; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 5;

echo "set origin 0.25,0.75; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 6;
echo "set origin 0.25,0.5; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 7;
echo "set origin 0.25,0.25; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 8;
echo "set origin 0.25,0; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 9;

echo "set origin 0.5,0.75; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 10;
echo "set origin 0.5,0.5; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 11;
echo "set origin 0.5,0.25; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 12;
echo "set origin 0.5,0; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 13;

#echo "set origin 0.75,0; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 14;
#echo "set origin 0.75,0.25; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 15;
#echo "set origin 0.75,0.5; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 16;
#echo "set origin 0.75,0.75; " >>$TEMPCOMMANDFILE; parse_fieldlog_column $TEMPCOMMANDFILE 17;

echo "unset multiplot;" >>$TEMPCOMMANDFILE;

# Run gnuplot
$GNUPLOT $TEMPCOMMANDFILE >>/dev/null 2>>/dev/null;

# Generate pdf
#ps2pdf $OUTFILENAME.ps;

# Remove command file
rm $TEMPCOMMANDFILE;

ALLPSFILES="$ALLPSFILES $OUTFILENAME.ps ";


# PARTICLE PLOTTING


# Logfile
LOGFILENAME=particles.log

# Check if the file exist
if [ ! -e $LOGFILENAME ]; then
 echo "ERROR: field log file ($LOGFILENAME) does not exists";
 exit -1;
fi  

echo -n "Generating particlelog plots.. population: ";

popinfos=""

npops=-1

parse_particlelog_column()
{
# Template gnuplot command file
TEMPCMDFILENAME=$1

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
 echo "";
 echo "ERROR: bad logfile column structure"
 exit -1
fi

# List particle populations into a string
poplist=""
nn=3
while [ $nn -le $((npops+2)) ]; do
 popname=$(eval "sed -n '${nn}p' $LOGFILENAME | cut -c3-999")
 popinfos[$nn]=$popname
 poplist="$poplist \"($popname)\""
 nn=$((nn+1))
done
poplist="$poplist \"all\""

# Select particle population
if [ "$2" != "all" ]; then
 population="command line"
 selected_pop=$2
else
 population="all"
 selected_pop=$((npops+1))
fi

# Select variable
variable="command line"
selected_var=$3

if [ "$population" != "all" ]; then
 # Column number
 column=$((1 + (selected_pop-1)*VARS_PER_POP + selected_var))
 # Variable's name
 varname=$(eval "cut -f $column $LOGFILENAME | sed -n '${nheader}p' | cut -c6-999")
 varname_tempA=$(echo $varname | cut -f 1 -d '(')
 varname_tempB=$(echo $varname | cut -f 2 -d ')') 
 varname=$varname_tempA$varname_tempB
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
else
 plotcmd="plot \"$LOGFILENAME\" using 1:$column title \"$varname\" w l"
fi

# Generate template gnuplot command file
echo "set xlabel \"Time [s]\"; $plotcmd; " >>$TEMPCMDFILENAME
}

echo -n "all";

TEMPCOMMANDFILE=tempAIOJEW6ferui.gnuplotcmd;

# Construct gnuplot command file
OUTFILENAME=log_population_all_a;
echo $GNUPLOT_TERMINAL >$TEMPCOMMANDFILE;
echo "set output '$OUTFILENAME.ps'; " >>$TEMPCOMMANDFILE;
echo "set multiplot; set size 0.25,0.25; " >>$TEMPCOMMANDFILE;

# Function call to get the population infos
parse_particlelog_column "/dev/null" all 1;

# All population informations
allpopinfos="POPULATIONS\n\npopname  particlemass [amu] particlecharge [e]\n";
mm=1
while [ $mm -le $npops ]; do
 allpopinfos="$allpopinfos\n$mm. ${popinfos[$((2 + mm))]}";
 mm=$((mm+1))
done

echo "set label \"HYBRIDCODE PARTICLE LOGS (all populations)\" at screen 0, screen $TEXT_TOPLIM font \"Times-Roman,12\";" >>$TEMPCOMMANDFILE;
echo "set label \"$allpopinfos\" at screen $TEXT_RIGHTLIM, screen $TEXT_TOPLIM font \"Times-Roman,10\";" >>$TEMPCOMMANDFILE;
echo "set label \"$(pwd)\" at screen 0, screen $TEXT_BOTTOMLIM font \"Times-Roman,14\";" >>$TEMPCOMMANDFILE;
echo "set label \"$(date -R)\" at screen $TEXT_RIGHTLIM, screen $TEXT_BOTTOMLIM font \"Times-Roman,12\";" >>$TEMPCOMMANDFILE;

echo "set origin 0,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column "$TEMPCOMMANDFILE" all 1;
echo "set origin 0,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column "$TEMPCOMMANDFILE" all 2;
echo "set origin 0,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column "$TEMPCOMMANDFILE" all 3;
echo "set origin 0,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column "$TEMPCOMMANDFILE" all 4;

echo "set origin 0.25,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 5;
echo "set origin 0.25,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 6;
echo "set origin 0.25,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 7;
echo "set origin 0.25,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 8;

echo "set origin 0.5,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 9;
echo "set origin 0.5,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 10;
#echo "set origin 0.5,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 11;
#echo "set origin 0.5,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 16;

#echo "set origin 0.75,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 13;
#echo "set origin 0.75,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 14;
#echo "set origin 0.75,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 15;
#echo "set origin 0.75,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 16;

echo "unset multiplot;" >>$TEMPCOMMANDFILE;

# Run gnuplot
$GNUPLOT $TEMPCOMMANDFILE >>/dev/null 2>>/dev/null;

# Generate pdf
#ps2pdf $OUTFILENAME.ps;

# Remove template command file
rm $TEMPCMDFILENAME;

ALLPSFILES="$ALLPSFILES $OUTFILENAME.ps ";

# Second page

# Construct gnuplot command file
OUTFILENAME=log_population_all_b;
echo $GNUPLOT_TERMINAL >$TEMPCOMMANDFILE;
echo "set output '$OUTFILENAME.ps'; " >>$TEMPCOMMANDFILE;
echo "set multiplot; set size 0.25,0.25; " >>$TEMPCOMMANDFILE;

# Function call to get the population infos
parse_particlelog_column "/dev/null" all 1;

# All population informations
allpopinfos="POPULATIONS\n\npopname  particlemass [amu] particlecharge [e]\n";
mm=1
while [ $mm -le $npops ]; do
 allpopinfos="$allpopinfos\n$mm. ${popinfos[$((2 + mm))]}";
 mm=$((mm+1))
done

echo "set label \"HYBRIDCODE PARTICLE LOGS (all populations)\" at screen 0, screen $TEXT_TOPLIM font \"Times-Roman,12\";" >>$TEMPCOMMANDFILE;
echo "set label \"$allpopinfos\" at screen $TEXT_RIGHTLIM, screen $TEXT_TOPLIM font \"Times-Roman,10\";" >>$TEMPCOMMANDFILE;
echo "set label \"$(pwd)\" at screen 0, screen $TEXT_BOTTOMLIM font \"Times-Roman,14\";" >>$TEMPCOMMANDFILE;
echo "set label \"$(date -R)\" at screen $TEXT_RIGHTLIM, screen $TEXT_BOTTOMLIM font \"Times-Roman,12\";" >>$TEMPCOMMANDFILE;

echo "set origin 0,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column "$TEMPCOMMANDFILE" all 11;
echo "set origin 0,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column "$TEMPCOMMANDFILE" all 12;
echo "set origin 0,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column "$TEMPCOMMANDFILE" all 13;
echo "set origin 0,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column "$TEMPCOMMANDFILE" all 14;

echo "set origin 0.25,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 15;
echo "set origin 0.25,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 16;
echo "set origin 0.25,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 17;
echo "set origin 0.25,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 18;

echo "set origin 0.5,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 19;
echo "set origin 0.5,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 20;
#echo "set origin 0.5,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 11;
#echo "set origin 0.5,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 16;

#echo "set origin 0.75,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 13;
#echo "set origin 0.75,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 14;
#echo "set origin 0.75,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 15;
#echo "set origin 0.75,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE all 16;

echo "unset multiplot;" >>$TEMPCOMMANDFILE;

# Run gnuplot
$GNUPLOT $TEMPCOMMANDFILE >>/dev/null 2>>/dev/null;

# Generate pdf
#ps2pdf $OUTFILENAME.ps;

# Remove template command file
rm $TEMPCMDFILENAME;

ALLPSFILES="$ALLPSFILES $OUTFILENAME.ps ";

# Populations
ii=1
while [ $ii -le $npops ]; do
 echo -n " $ii";

 # Function call to get the population infos
 parse_particlelog_column "/dev/null" all 1;

 # Population informations
 populationinfo="POPULATIONS\n\npopname  particlemass [amu] particlecharge [e]\n";
 populationinfo="$populationinfo\n$ii. ${popinfos[$((2 + ii))]}";

 OUTFILENAME="log_population_${ii}_a";
 echo $GNUPLOT_TERMINAL >$TEMPCOMMANDFILE;
 echo "set output '$OUTFILENAME.ps'; " >>$TEMPCOMMANDFILE;
 echo "set multiplot; set size 0.25,0.25; " >>$TEMPCOMMANDFILE;

 echo "set label \"HYBRIDCODE PARTICLE LOGS\" at screen 0, screen $TEXT_TOPLIM font \"Times-Roman,12\";" >>$TEMPCOMMANDFILE;
 echo "set label \"$populationinfo\" at screen $TEXT_RIGHTLIM, screen $TEXT_TOPLIM font \"Times-Roman,10\";" >>$TEMPCOMMANDFILE;
 echo "set label \"$(pwd)\" at screen 0, screen $TEXT_BOTTOMLIM font \"Times-Roman,14\";" >>$TEMPCOMMANDFILE;
 echo "set label \"$(date -R)\" at screen $TEXT_RIGHTLIM, screen $TEXT_BOTTOMLIM font \"Times-Roman,12\";" >>$TEMPCOMMANDFILE;

 echo "set origin 0,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 1;
 echo "set origin 0,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 2;
 echo "set origin 0,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 3;
 echo "set origin 0,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 4;

 echo "set origin 0.25,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 5;
 echo "set origin 0.25,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 6;
 echo "set origin 0.25,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 7;
 echo "set origin 0.25,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 8;

 echo "set origin 0.5,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 9;
 echo "set origin 0.5,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 10;
 #echo "set origin 0.5,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 11;
 #echo "set origin 0.5,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 16;

 #echo "set origin 0.75,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 13;
 #echo "set origin 0.75,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 14;
 #echo "set origin 0.75,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 15;
 #echo "set origin 0.75,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 16;

 echo "unset multiplot;" >>$TEMPCOMMANDFILE;

 $GNUPLOT $TEMPCOMMANDFILE >>/dev/null 2>>/dev/null;
 #ps2pdf $OUTFILENAME.ps;
 rm $TEMPCOMMANDFILE;
 ALLPSFILES="$ALLPSFILES $OUTFILENAME.ps "
 
 OUTFILENAME="log_population_${ii}_b";
 echo $GNUPLOT_TERMINAL >$TEMPCOMMANDFILE;
 echo "set output '$OUTFILENAME.ps'; " >>$TEMPCOMMANDFILE;
 echo "set multiplot; set size 0.25,0.25; " >>$TEMPCOMMANDFILE;

 echo "set label \"HYBRIDCODE PARTICLE LOGS\" at screen 0, screen $TEXT_TOPLIM font \"Times-Roman,12\";" >>$TEMPCOMMANDFILE;
 echo "set label \"$populationinfo\" at screen $TEXT_RIGHTLIM, screen $TEXT_TOPLIM font \"Times-Roman,10\";" >>$TEMPCOMMANDFILE;
 echo "set label \"$(pwd)\" at screen 0, screen $TEXT_BOTTOMLIM font \"Times-Roman,14\";" >>$TEMPCOMMANDFILE;
 echo "set label \"$(date -R)\" at screen $TEXT_RIGHTLIM, screen $TEXT_BOTTOMLIM font \"Times-Roman,12\";" >>$TEMPCOMMANDFILE;

 echo "set origin 0,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 11;
 echo "set origin 0,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 12;
 echo "set origin 0,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 13;
 echo "set origin 0,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 14;

 echo "set origin 0.25,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 15;
 echo "set origin 0.25,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 16;
 echo "set origin 0.25,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 17;
 echo "set origin 0.25,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 18;

 echo "set origin 0.5,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 19;
 echo "set origin 0.5,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 20;
 #echo "set origin 0.5,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 11;
 #echo "set origin 0.5,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 16;

 #echo "set origin 0.75,0; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 13;
 #echo "set origin 0.75,0.25; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 14;
 #echo "set origin 0.75,0.5; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 15;
 #echo "set origin 0.75,0.75; " >>$TEMPCOMMANDFILE; parse_particlelog_column $TEMPCOMMANDFILE $ii 16;

 echo "unset multiplot;" >>$TEMPCOMMANDFILE;

 $GNUPLOT $TEMPCOMMANDFILE >>/dev/null 2>>/dev/null;
 #ps2pdf $OUTFILENAME.ps;
 rm $TEMPCOMMANDFILE;
 ALLPSFILES="$ALLPSFILES $OUTFILENAME.ps "

 ii=$((ii+1))
done

echo "";

# Generate all-in-one pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=hybridlogs.pdf $ALLPSFILES;

# Remove postscript files
rm $ALLPSFILES;

exit 0;
