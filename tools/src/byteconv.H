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

#ifndef BYTECONV_H

#ifdef __GNUC__
#  pragma interface
#endif

#include <cmath>	// this defines LITTLE_ENDIAN on some systems, in which case preprocessor warns
#include <cstdlib>	// (also this can behave the same)
// about doubly defined symbol. Including it BEFORE the next if...endif ensures that this does not happen.

#include <iostream>
using namespace std;

#if (defined(__i386) || defined(__alpha)) && !defined(LITTLE_ENDIAN)
#  define LITTLE_ENDIAN
#endif

#define WRITE_BUFLEN 1000

#ifdef LITTLE_ENDIAN
extern void ByteConversion(int sz, unsigned char *x, int n);
#else
inline void ByteConversion(int sz, unsigned char *x, int n) {}
#endif
#define ByteConversion_input ByteConversion

#if defined(_CRAY1) && !defined(_CRAYIEEE)

extern "C" int IEG2CRAY(const int& type, const int& num,
						const float foreign[], const int& bitoff,
						double crayarr[], const int& stride=1, char craych[]=0);

extern "C" int CRAY2IEG(const int& type, const int& num,
						float foreign[], const int& bitoff,
						const double crayarr[], const int& stride=1, char craych[]=0);

extern void doabort()
#ifdef __GNUC__
	__attribute__((noreturn))
#endif
;

inline void WriteFloatsToFile(ostream& o, const float xf[], int n)
{
	float foreign[WRITE_BUFLEN];
	if (n > WRITE_BUFLEN) {cerr << "*** n > " << WRITE_BUFLEN << " in WriteDoubles\n"; doabort();}
	const int ierr = CRAY2IEG(2,n,foreign,0,(const double*)xf);
	if (ierr < 0) cerr << "*** CRAY2IEG/float: parameter error\n";
	o.write((char*)foreign,4*n);
}

inline void WriteDoublesToFile(ostream& o, const double x[], int n)
{
	float foreign[WRITE_BUFLEN];
	if (n > WRITE_BUFLEN) {cerr << "*** n > " << WRITE_BUFLEN << " in WriteDoubles\n"; doabort();}
	const int ierr = CRAY2IEG(8,n,foreign,0,x);
	if (ierr < 0) cerr << "*** CRAY2IEG/double: parameter error\n";
	o.write((char*)foreign,8*n);
}

inline void ReadFloatsFromFile(istream& o, float xf[], int n)
{
	float foreign[WRITE_BUFLEN];
	if (n > WRITE_BUFLEN) {cerr << "*** n > " << WRITE_BUFLEN << " in ReadDoubles\n"; doabort();}
	o.read((char*)foreign,4*n);
	const int ierr = IEG2CRAY(2,n,foreign,0,(double *)xf);
	if (ierr < 0) cerr << "*** IEG2CRAY/float: parameter error\n";
}

inline void ReadDoublesFromFile(istream& o, double x[], int n)
{
	float foreign[WRITE_BUFLEN];
	if (n > WRITE_BUFLEN) {cerr << "*** n > " << WRITE_BUFLEN << " in ReadDoubles\n"; doabort();}
	o.read((char*)foreign,8*n);
	const int ierr = IEG2CRAY(8,n,foreign,0,x);
	if (ierr < 0) cerr << "*** IEG2CRAY/double: parameter error\n";
}

#else

inline void WriteFloatsToFile(ostream& o, const float xf[], int n) {o.write((const char*)xf,sizeof(float)*n);}
inline void WriteDoublesToFile(ostream& o, const double x[], int n) {o.write((const char*)x,sizeof(double)*n);}
inline void ReadFloatsFromFile(istream& o, float xf[], int n) {o.read((char*)xf,sizeof(float)*n);}
inline void ReadDoublesFromFile(istream& o, double x[], int n) {o.read((char*)x,sizeof(double)*n);}

#endif

#define BYTECONV_H

#endif
