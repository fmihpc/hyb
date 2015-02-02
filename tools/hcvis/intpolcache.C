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
#  pragma implementation "intpolcache.H"
#endif

#include "intpolcache.H"

TIntpolCache::TIntpolCache(real Dx1, real x01, real y01)
	: nstores(0), nreads(0), nhits(0), Dx(Dx1), x0(x01), y0(y01), nbits(16)
{
	int i;
	for (i=0; i<INTPOL_CACHE_LEN; i++) lists[i] = 0;
}

void TIntpolCache::store(real u, real x, real y)
{
	int ix,iy;
	quantize(x,y, ix,iy);
	const int k = hashkey(ix,iy);
	Tnode *const n = new Tnode;
	n->next = lists[k];
	n->ix = ix;
	n->iy = iy;
	n->u = u;
	lists[k] = n;
	nstores++;
}

bool TIntpolCache::read(real& u, real x, real y)
{
	int ix,iy;
	quantize(x,y, ix,iy);
	const int k = hashkey(ix,iy);
	if (k < 0 || k >= INTPOL_CACHE_LEN) {
		cerr << "*** k = " << k << "\n";
	}
	Tnode *p;
	bool retval = false;
	for (p=lists[k]; p; p=p->next)
		if (p->ix == ix && p->iy == iy) {
			u = p->u;
			nhits++;
			retval = true;
		}
	nreads++;
	return retval;
}

TIntpolCache::~TIntpolCache()
{
	int i;
	Tnode *p,*q;
	for (i=0; i<INTPOL_CACHE_LEN; i++) {
		for (p=lists[i]; p; p=q) {
			q = p->next;
			delete p;
		}
	}	
}

