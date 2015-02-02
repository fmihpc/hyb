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
#  pragma implementation "cellcache.H"
#endif

#include "cellcache.H"
#include <stdlib.h>		// to get exit()

TCellCache::TCellCache()
{
	int i;
	for (i=0; i<CELLCACHE_HASHTABLE_SIZE; i++) hash[i] = 0;
	n = 0;
	maxn = CELLCACHE_HASHTABLE_SIZE;
	first = 0;
	last = 0;
	dx = 0.01;
	invdx = 1.0/dx;
	nsuccessful_probes = nfailed_probes = 0;
}

void TCellCache::quantize(const Tdimvec& X, smallnat dim, long iq[3]) const
{
	smallnat d;
	for (d=0; d<dim; d++) iq[d] = (long)(X(d)*invdx*65536.0 + 0.5);
}

int TCellCache::hashfunction(const long iq[3], smallnat dim) const
{
	unsigned long x = iq[0] << 7;
	if (dim >= 2) {
		x+= iq[1];
		if (dim == 3)
			x+= iq[2] >> 7;
	}
	return int(x % CELLCACHE_HASHTABLE_SIZE);
}

void TCellCache::clear()
{
	int i;
	Tnode *p;
	for (i=0; i<CELLCACHE_HASHTABLE_SIZE; i++) {
		while (hash[i]) {
			p = hash[i];
			hash[i] = p->next;
			delete p;
		}
	}
	n = 0;
	first = 0;
	last = 0;
	nsuccessful_probes = nfailed_probes = 0;
}

void TCellCache::printstats(ostream& o)
{
	long n = nfailed_probes + nsuccessful_probes;
	if (n==0) n = 1;
	o << "TCellCache: " << 100.0*double(nsuccessful_probes)/double(n) << "% success rate";
}

void TCellCache::store(const Tdimvec& X, smallnat dim, const real u[max_ncd], smallnat ncd)
{
	long iq[3];
	quantize(X,dim,iq);
	const int i = hashfunction(iq,dim);
	Tnode *const p = new Tnode;
	p->iq[0] = iq[0];
	p->iq[1] = iq[1];
	p->iq[2] = iq[2];
	smallnat a;
	for (a=0; a<ncd; a++)
		p->u[a] = u[a];
	p->next = hash[i];
	hash[i] = p;
	p->next2 = 0;
	if (last) last->next2 = p;
	if (!first) first = p;
	last = p;
	n++;
	if (n > maxn) {
		// Delete the least-frequently stored node
		const int j = hashfunction(first->iq,dim);
		Tnode *q = hash[j];
		if (q == first) {
			hash[j] = q->next;
			first = q->next2;
			delete q;
			n--;
		} else {
			while (q && q->next != first) q = q->next;
			if (!q) {
				cerr << "*** Internal error in cellcache.C::store()\n";
				exit(1);
			}
			// now q->next == first
			Tnode *const del = q->next;
			q->next = del->next;
			first = del->next2;
			delete del;
		}
	}
}

bool TCellCache::probe(const Tdimvec& X, smallnat dim, real u[max_ncd], smallnat ncd)
{
	long iq[3];
	quantize(X,dim,iq);
	const int i = hashfunction(iq,dim);
	Tnode *p;
	smallnat d,a;
	for (p=hash[i]; p; p=p->next) {
		bool found = true;
		for (d=0; d<dim; d++)
			if (iq[d] != p->iq[d]) {
				found = false;
				break;
			}
		if (found) {
			for (a=0; a<ncd; a++) u[a] = p->u[a];
			nsuccessful_probes++;
			return true;
		}
	}
	nfailed_probes++;
	return false;
}
