#!/bin/bash

# This file is part of the HYB simulation platform.
#
# Copyright 2018- Aalto University
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

CMDFILE=gnuplot.cmd

logfile=$1
uselogscale=0

# multiple dirs
if [ "$#" -le "0" ]; then
 echo "USAGE: hyblog_pdf.sh field.log lin/log [dir1/] [dir2/] [dir3/] [...]"
 exit 1
elif [ "$#" -eq "2" ]; then
 sdirs=./
 logfile_full=$logfile
else
 sdirs=$(echo $@ | cut -d " " -f 3-)
 logfile_full=$(echo $sdirs |cut -f 1 -d ' ')$logfile
fi

if [ "$2" == "lin" ]; then
 uselogscale=0
elif [ "$2" == "log" ]; then
 uselogscale=1
 else
 echo "USAGE: hyblog_pdf.sh field.log lin/log [dir1/] [dir2/] [dir3/] [...]"
 exit 1
fi

# parse header
str1=$(sed -n '2p' $logfile_full | cut -f 2 -d ' ')
str2=$(sed -n '4p' $logfile_full | cut -f 2 -d ' ')
if [ "$str1" = "columns" ]; then
 echo "field log";
 cols=$(sed -n '2p' $logfile_full | cut -c13-)
 nc=2
elif [ "$str2" = "columns" ]; then
 echo "particle log";
 cols=$(sed -n '4p' $logfile_full | cut -c13-)
 nc=4
else
 echo "bad file structure";
 exit -1;
fi

# header lines
hl=$((n+cols-1))

psfiles=$(echo $logfile | cut -f 1 -d .)
pdffile=$(echo $logfile | cut -f 1 -d .).pdf

echo "plotting $cols data columns from $logfile_full"

echo "set terminal postscript landscape color \"Times-Roman\" 12;"  >$CMDFILE
echo set xlabel \"Time [s]\" >>$CMDFILE

# construct gnuplot command file
i=1
while [ $i -le $cols ]; do
 
 # output file
 echo -n set output \"$psfiles >>$CMDFILE
 if [ $i -lt 10 ]; then
  echo -n 000 >>$CMDFILE
 elif [ $i -lt 100 ]; then
  echo -n 00 >>$CMDFILE
 elif [ $i -lt 1000 ]; then
  echo -n 0 >>$CMDFILE
 fi
 echo ${i}.ps\" >>$CMDFILE
 
 # title
 ntitle=$((nc+i))
 title=$(sed -n ''${ntitle}p'' $logfile_full | cut -c3-)
 echo set ylabel \"$title\" >>$CMDFILE
 if [ $uselogscale == 1 ] && [[ "$title" =~ ^("05. avg(|B|) [T]"|"06. max(|B|) [T]"|"10. energy(sum(dV*B^2/2*mu0)) [J]"|"07. avg(|V|) [m/s]"|"08. Kinetic energy [J]")$ ]] ; then
  echo set logscale y >>$CMDFILE
  echo set format y \"10^{%L}\" >>$CMDFILE
 else
  echo unset logscale y >>$CMDFILE
  echo unset format >>$CMDFILE
 fi
 
 # plot command
 #echo plot \"$logfile_full\" using 1:$i title \"\" w l >>$CMDFILE
 echo -n  plot >>$CMDFILE
 comma=''
 for sdir in $sdirs ; do
  echo -n $comma\"$sdir$logfile\" using 1:$i w l >>$CMDFILE
  comma=','
 done
 echo >>$CMDFILE
 
 i=$((i+1))
done
  
gnuplot $CMDFILE >>/dev/null 2>>/dev/null
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$pdffile $psfiles*.ps
rm $psfiles*.ps $CMDFILE
