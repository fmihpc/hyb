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

# multiple dirs
if [ "$#" -le "0" ]; then
 echo "USAGE: hyblog_pdf_all.sh lin/log [dir1/] [dir2/] [dir3/] [...]"
elif [ "$#" -eq "1" ]; then
 sdirs=""
 firstdir=./
else
 sdirs=$(echo $@ | cut -d " " -f 2-)
 firstdir=$2
fi

if ! [[ "$1" =~ ^("lin"|"log")$ ]]; then
 echo "USAGE: hyblog_pdf_all.sh lin/log [dir1/] [dir2/] [dir3/] [...]"
 exit 1
fi

for i in $(/bin/ls ${firstdir}pop*.log)
do
 i=$(basename $i)
 echo "plotting $i"
 hyblog_pdf.sh $i $1 $sdirs
done
echo "plotting field.log"
hyblog_pdf.sh field.log $1 $sdirs
