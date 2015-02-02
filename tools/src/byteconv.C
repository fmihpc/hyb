/** This file is part of the HYB simulation platform.
 *
 *  Copyright 2014- Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef __GNUC__
#  pragma implementation "byteconv.H"
#endif

#include "byteconv.H"
#include <iostream>
using namespace std;

extern void doabort();

#ifdef LITTLE_ENDIAN
void ByteConversion(int sz, unsigned char *x, int n)
{
	int i,ind,j;
	unsigned int b[16];
	if (sz > 16) {cerr << "*** ByteConversion:: internal error\n"; doabort();}
	for (i=ind=0; i<n; i++,ind+=sz) {
		for (j=0; j<sz; j++) b[j] = x[ind+j];
		for (j=0; j<sz; j++) x[ind+sz-1-j] = b[j];
	}
}
#endif
