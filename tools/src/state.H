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

#ifndef STATE_H

#ifdef __GNUC__
#  pragma interface
#endif

#define USE_BMHD 0		/* if 1, include BMHD.H (Boris-MHD), otherwise MHD.H. Tested everywhere where (B)MHD.H is included. */

#include "gengrid.H"
#include "vec.H"
#include <math.h>

class params {
private:
	static real gam;				// adiabatic constant
	static real gamminus1;			// gam - 1
	static real inv_gamminus1;		// 1/(gam-1)
	static real mu0;
	static real sqrtmu0;			// sqrt(mu0)
	static real invmu0;				// 1/mu0
	static real invsqrtmu0;			// 1/sqrt(mu0)
	static real c;					// (possibly reduced, cf. Boris-correction) speed of light
	static real invc2;				// 1/c^2
	static real eps0;				// 1/(mu0*c^2)
public:
	static void setgamma(real g) {gam=g; gamminus1=gam-1; inv_gamminus1=1/gamminus1;}
	static void setmu0(real m) {mu0=m; sqrtmu0=sqrt(mu0); invmu0=1/mu0; invsqrtmu0=sqrt(invmu0);}
	static void setc(real c1) {c=c1; invc2=1.0/(c*c); eps0=1.0/(mu0*c*c);}
	static real gamma() {return gam;}
	static real gamma_minus_1() {return gamminus1;}
	static real inv_gamma_minus_1() {return inv_gamminus1;}
	static real Mu0() {return mu0;}
	static real SqrtMu0() {return sqrtmu0;}
	static real InvMu0() {return invmu0;}
	static real InvSqrtMu0() {return invsqrtmu0;}
	static real LightSpeed() {return c;}
	static real InvC2() {return invc2;}
	static real Epsilon0() {return eps0;}
};

template <smallnat sz, smallnat sz2> class Tstate {
private:
	int len;
#	ifdef TRANSPOSED_STORAGE
	real data[sz+sz2][VECLEN];
	real& uref(int i, smallnat v=0) {
		assert(0<=i && i<sz+sz2);
		assert(v<VECLEN);
		return data[i][v];
	}
	real u(int i, smallnat v=0) const {
		assert(0<=i && i<sz+sz2);
		assert(v<VECLEN);
		return data[i][v];
	}
#	else
	real data[VECLEN][sz];
	real& uref(int i, smallnat v=0) {
		assert(0<=i && i<sz+sz2);
		assert(v<VECLEN);
		return data[v][i];
	}
	real u(int i, smallnat v=0) const {
		assert(0<=i && i<sz+sz2);
		assert(v<VECLEN);
		return data[v][i];
	}
#	endif
public:
	Tstate() {len=1;}
	void setlength(int L) {len=L;}
	int length() const {return len;}
	real& operator()(int i, smallnat v) {return uref(i,v);}
	real& operator[](int i) {return uref(i,0);}
	real operator()(int i) const {return u(i,0);}
	static int size() {return sz;}
	static void permutation(smallnat dim, smallnat perm[sz]);
};

// should_also_be_const is used in MHD.C, Euler.C
#ifndef should_also_be_const
#  if defined(__GNUC__)
//   Newer g++ vers (2.7.2) has problems with this!
#    define should_also_be_const
#  else
#    define should_also_be_const const
#  endif
#endif

#define STATE_H

#endif
