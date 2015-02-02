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
#  pragma implementation "cache.H"
#endif

#include <stdlib.h>
#include "cache.H"
#include "shmemintf.H"

#define CACHE_SIZE 7000

Tcache Dcache;

void Tcache::init()
{
	cache_size = CACHE_SIZE;
	cache = new TCacheWord [2*CACHE_SIZE];
	cache = (TCacheWord *)(CACHE_SIZE*(long(cache + CACHE_SIZE)/CACHE_SIZE));
	words_allocated = 0;
	verbose_flag = false;
}

TCacheWord *Tcache::alloc(int sz, const char *msg)
// allocate sz bytes from the cache, aligning at TCacheWord boundary
{
	const int words_to_allocate = (sz+sizeof(TCacheWord)-1)/sizeof(TCacheWord);
	if (words_to_allocate < CACHE_SIZE - words_allocated) {
		// we can satisfy the alloc request
		TCacheWord *const result = &cache[words_allocated];
		words_allocated+= words_to_allocate;
		if (verbose_flag)
			ATOMIC_OUTPUT(
				o << "Tcache::alloc ";
				if (msg) o << msg << ' ';
				o << words_to_allocate << ": "
				  << words_allocated << "/" << CACHE_SIZE << "\n";
			);
		return result;
	} else {
		cerr << "*** Tcache::alloc(" << sz << "): cache overflow\n";
		exit(1);
		return 0;
	}
}
