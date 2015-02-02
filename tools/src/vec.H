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

#ifndef VEC_H

#ifdef __GNUC__
#  pragma interface
#endif

#include "realptr.H"	// to get real, smallnat and VECLEN

class Tvec {
private:
	real data[VECLEN];
	static smallnat veclen;
public:
	Tvec() {}
	static smallnat vl() {return veclen;}
	static void setlen(smallnat vlen) {veclen = vlen;}
	real& operator[](smallnat i) {return data[i];}
	real operator()(smallnat i) const {return data[i];}
};

class Ttrivialvec {
private:
	real data;
public:
	static smallnat vl() {return 1;}
	static void setlen(smallnat) {}
	real& operator[](smallnat) {return data;}
	real operator()(smallnat) const {return data;}
	operator real() const {return data;}
};

#define VL (VECTOR::vl())

#ifdef _CRAY1
#  define VECSTATIC static
#else
#  define VECSTATIC
#endif

extern void genbase(const Tvec n[3], Tvec t1[3], Tvec t2[3], smallnat vlen);

#define VEC_H

#endif
