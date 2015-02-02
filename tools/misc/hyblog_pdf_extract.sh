#!/bin/bash
# Old name was of this script was: extractlogpdfs.sh

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

if ! hash pdftk 2>/dev/null; then
 echo "pdftk not found"
 exit -1
fi
if [ "$#" != "3" ]; then
 echo "USAGE: hyblog_pdf_extract.sh colN \"pop*.pdf\" out.pdf"
 exit -1
fi
if [ -e "$3" ]; then
 echo "$3 already exists"
 exit -1
fi
colN=$1
files=$(ls $2)
out=$3
i=1
pagestr=""
for x in $files; do
 p=temp_abc421_page$i.pdf
 pdftk A=$x cat A$colN output $p
 pagestr="$pagestr $p"
 i=$((i+1))
done
pdftk $pagestr cat output $out
rm $pagestr
