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
#  pragma implementation "vec.H"
#endif

#include "vec.H"
#include <math.h>
#define sqr(x) ((x)*(x))

smallnat Tvec::veclen = 0;


#define VECTOR Tvec

void genbase(const VECTOR n[3], VECTOR t1[3], VECTOR t2[3], smallnat vlen)
// Generate vlen orthonormal bases (n,t1,t2). Input: unit vectors n. Output: unit vectors t1,t2.
// The generated basis (n,t1,t2) will be right-handed. The input vectors n MUST be normalized before call.
// Also, all 3 components of n must be set. In 2D, set n[2] to zero. In 1D, set n[1]=n[2]=0.
// In 2D, only t1 is sensible and has t1[2]=0. t2 is computed but is z-directional.
/* Algorithm:
   (1) Choose e as one of ex,ey,ez (corresponding to the smallest component of n)
   (2) Form the cross product t1 = e x n
   (3) Normalize t1: t1=t1/|t1|
   (4) Form the cross product t2 = n x t1
In 2D: nz=0, d[v] = 2 always (e=ez)
In 1D: d[v] = 2 always also. t1 will be +-ey, t2 will be +-ez.
   */
{
	VECTOR::setlen(vlen);
	smallnat v;
	smallnat d[VECLEN];		// will be either 0,1,2 depending which element of n3 is smallest
	VECSTATIC VECTOR e[3],invt1len;
	VDIRS;
	for (v=0; v<VL; v++) {
		const real nx = fabs(n[0](v));
		const real ny = fabs(n[1](v));
		const real nz = fabs(n[2](v));
		if (nx < ny) {
			d[v] = (nx < nz) ? 0 : 2;
		} else {
			d[v] = (ny < nz) ? 1 : 2;
		}
		e[0][v] = e[1][v] = e[2][v] = 0;
	}
	VDIRS;
	for (v=0; v<VL; v++) e[d[v]][v] = 1;
	VDIRS;
	for (v=0; v<VL; v++) {
		t1[0][v] = -n[1](v)*e[2](v) + n[2](v)*e[1](v);
		t1[1][v] = -n[2](v)*e[0](v) + n[0](v)*e[2](v);
		t1[2][v] = -n[0](v)*e[1](v) + n[1](v)*e[0](v);
		invt1len[v] = 1.0/sqrt(sqr(t1[0](v)) + sqr(t1[1](v)) + sqr(t1[2](v)));
	}
	VDIRS;
	for (v=0; v<VL; v++) {
		t1[0][v]*= invt1len(v);
		t1[1][v]*= invt1len(v);
		t1[2][v]*= invt1len(v);
	}
	VDIRS;
	for (v=0; v<VL; v++) {
		t2[0][v] = n[1](v)*t1[2](v) - n[2](v)*t1[1](v);
		t2[1][v] = n[2](v)*t1[0](v) - n[0](v)*t1[2](v);
		t2[2][v] = n[0](v)*t1[1](v) - n[1](v)*t1[0](v);
	}
}

#undef VECTOR
