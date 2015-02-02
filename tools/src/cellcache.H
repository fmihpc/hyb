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

#ifndef CELLCACHE_H

#ifdef __GNUC__
#  pragma interface
#endif

#include "realptr.H"		// to get real and smallnat
#include <iostream>
using namespace std;
#include "BC.H"				// to get Tdimvec
#define max_ncd 1000
#define CELLCACHE_HASHTABLE_SIZE 131

class TCellCache {
private:
	struct Tnode {
		long iq[3];
		real u[max_ncd];
		Tnode *next;		// pointer to next in synonym chain
		Tnode *next2;		// pointer to next temporally added node
	};
	Tnode *hash[CELLCACHE_HASHTABLE_SIZE];
	int n;					// number of entries in hash table
	int maxn;				// maximum allowed n
	real dx, invdx;
	Tnode *first;			// next2-chained list, in temporal order
	Tnode *last;			// pointer to last element in next2-chainged first-list
	long nsuccessful_probes, nfailed_probes;		// for statistics
	void quantize(const Tdimvec& X, smallnat dim, long iq[3]) const;
	int hashfunction(const long iq[3], smallnat dim) const;
public:
	TCellCache();
	void set_dx(real dx1) {dx=dx1; invdx=1.0/dx;}
	void store(const Tdimvec& X, smallnat dim, const real u[max_ncd], smallnat ncd);
	bool probe(const Tdimvec& X, smallnat dim, real u[max_ncd], smallnat ncd);
	void printstats(ostream& o);
	void clear();
	~TCellCache() {clear();}
};

#define CELLCACHE_H

#endif
