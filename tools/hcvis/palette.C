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
#  pragma implementation "palette.H"
#endif

#include "palette.H"
#include <cstdio>
#include <iostream>
#include <cmath>
using namespace std;

#define MAX_N 512

void Tpalette::default_init()
{
	int i;
	N = 256;
	pal = new unsigned char [3*256];
	for (i=0; i<256; i++) {
		pal[3*i] = pal[3*i+1] = pal[3*i+2] = i;
	}
}

void Tpalette::invert()
{
	int i,j;
	unsigned char tmp;
#define SWAP(x,y) tmp=x; x=y; y=tmp
	for (i=0; i<N/2; i++) {
		j = N-1-i;
		SWAP(pal[3*i],pal[3*j]);
		SWAP(pal[3*i+1],pal[3*j+1]);
		SWAP(pal[3*i+2],pal[3*j+2]);
	}
#undef SWAP
	pal[3*(N-1)] = pal[3*(N-2)] = pal[3*(N-3)];
	pal[3*(N-1)+1] = pal[3*(N-2)+1] = pal[3*(N-3)+1];
	pal[3*(N-1)+2] = pal[3*(N-2)+2] = pal[3*(N-3)+2];
}

void Tpalette::band(int n)
{
	int i;
	for (i=0; i<N; i++) {
		const double s = sin(2*3.14159265358979323846*i*double(n)/double(N));
		const double f = 1 - 0.5*s*s;
		pal[3*i  ] = (unsigned char)(f*pal[3*i  ]);
		pal[3*i+1] = (unsigned char)(f*pal[3*i+1]);
		pal[3*i+2] = (unsigned char)(f*pal[3*i+2]);
	}
}

void Tpalette::dim()
{
	int i;
	for (i=0; i<3*N; i++) pal[i] = (unsigned char)(0.85*pal[i]);
}

void Tpalette::load(const char *fn)
{
	int i;
	const int RGBorder = 0;
	unsigned char pal1[MAX_N*3];
	FILE *const fp = fopen(fn,"r");
	if (!fp) {
		cerr << "*** Tpalette::Tpalette: Cannot open file \"" << fn << "\"\n";
		return;
	}
	for (i=0; i<MAX_N; i++) {
		pal1[3*i  ] = fgetc(fp);
		pal1[3*i+1] = fgetc(fp);
		pal1[3*i+2] = fgetc(fp);
		if (feof(fp)) break;
	}
	N = i-1;
	pal = new unsigned char [3*N];
	for (i=0; i<N; i++) {
		if (RGBorder) {
			pal[3*i  ] = pal1[3*i  ];
			pal[3*i+1] = pal1[3*i+1];
			pal[3*i+2] = pal1[3*i+2];
		} else {
			pal[3*i  ] = pal1[i];
			pal[3*i+1] = pal1[i+N];
			pal[3*i+2] = pal1[i+2*N];
		}
	}
	fclose(fp);
}
