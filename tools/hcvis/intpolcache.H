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

#ifndef INTPOLCACHE_H

#ifdef __GNUC__
#  pragma interface
#endif

#include "realptr.H"		// to get real

#define INTPOL_CACHE_LEN 1223

class TIntpolCache {
private:
	int nstores, nreads, nhits;
	real Dx, x0, y0;
	int nbits;
	struct Tnode {
		Tnode *next;
		int ix, iy;
		real u;
	};
	Tnode *lists[INTPOL_CACHE_LEN];
	void quantize(real x, real y, int& ix, int& iy) const {
		ix = int((1<<nbits)*((x-x0)/Dx) + 0.5);
		iy = int((1<<nbits)*((y-y0)/Dx) + 0.5);
	}
	int hashkey(int x, int y) const {return (unsigned(x) + (unsigned(y) << 15)) % INTPOL_CACHE_LEN;}
public:
	TIntpolCache(real Dx1, real x01, real y01);
	void store(real u, real x, real y);
	bool read(real& u, real x, real y);
	real hitratio() const {return nreads ? nhits/real(nreads) : 0;}
	~TIntpolCache();
};

#define INTPOLCACHE_H

#endif
