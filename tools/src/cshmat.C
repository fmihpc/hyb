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
#  pragma implementation "cshmat.H"
#endif

#include "cshmat.H"
#include <iostream>	// for cerr
#include <cstdlib>		// for exit()
#include <cstring>		// for memcpy()
#if USE_SHMEM
#  include <cmalloc>	// for shmalloc, shfree
#endif
using namespace std;

#if USE_SHMEM
// defined in sharr.C
extern void *sharedmalloc(size_t);
extern void sharedfree(void *);
#endif

template <class T1, class T2>
TCachedSharedMatrix<T1,T2>::TCachedSharedMatrix()
{
	ptr = 0;
	nrows = 0; ncols1 = ncols2 = ncols = 0;
	cache_size = 0;
	cache = 0;
	n_cells_in_cache = 0;
}

template <class T1, class T2>
void TCachedSharedMatrix<T1,T2>::init(TGridIndex nrows1, TCompIndex ncols_1, TCompIndex ncols_2)
{
	if (ptr && nrows > 0 && ncols > 0) {
#		if USE_SHMEM
		cerr << "TCachedSharedMatrix::init warning: Calling shfree(ptr)\n";
		sharedfree(ptr);
#		else
		delete [] ptr;
#		endif
	}
//	cache_size = 31;
	cache_size = 64;
	/* cache_size 31 = 1 + 6 + 4*6. In 3D, one cell (1) has up to 6 immediate neighbours.
	   In the worst case, each neighbour is denser than the first cell, and thus has 4 interfaces. */
	nrows = nrows1;
	ncols1 = (ncols_1*sizeof(T1)+sizeof(TCacheWord)-1)/sizeof(TCacheWord);
	ncols2 = (ncols_2*sizeof(T2)+sizeof(TCacheWord)-1)/sizeof(TCacheWord);
	ncols = ncols1 + ncols2;
#	if USE_SHMEM
	nrows = DivNpes(nrows + Npes-1)*Npes;
	const TGridIndex nrows_local = DivNpes(nrows);
//	shmembarrierall();		// not necessary?
	ptr = (TCacheWord *)sharedmalloc(nrows_local*ncols*sizeof(TCacheWord));
#	else
	ptr = new TCacheWord [nrows*ncols];
	memset(ptr,0,nrows*ncols*sizeof(TCacheWord));
#	endif
	if (Dcache.verboseflag() && mype==0)
		cout << "TCachedSharedMatrix: "
			 << ncols1*sizeof(TCacheWord) << "+" << ncols2*sizeof(TCacheWord) << "="
			 << ncols*sizeof(TCacheWord) << " bytes of memory per cell\n" << flush;
	if (!ptr) {
		cerr << "*** TCachedSharedMatrix(" << nrows1 << "," << ncols_1 << "," << ncols_2 << "): out of memory\n";
		exit(1);
	}
	if (!cache) {
		const int nbytes1 = cache_size*ncols*sizeof(TCacheWord);
		const int nbytes2 = cache_size*sizeof(TCacheWord);
		const int nbytes = nbytes1 + nbytes2;
		if (sizeof(TCacheWord) < sizeof(TGridIndex)) {
			cerr << "*** TCachedSharedMatrix::init: sizeof(TCacheWord)=" << sizeof(TCacheWord)
				 << ", sizeof(TGridIndex)=" << sizeof(TGridIndex) << "\n";
			exit(1);
		}
		if (Dcache.bytes_free() >= nbytes) {
			cache = (TCacheWord *)Dcache.alloc(nbytes1,"cell cache");
			indices_in_cache = (TGridIndex *)Dcache.alloc(nbytes2,"cache tags");
			cache_was_allocated_using_new = false;
		} else {
			cache = new TCacheWord [cache_size*ncols];
			indices_in_cache = new TGridIndex [cache_size];
			cache_was_allocated_using_new = true;
			cerr << "TCachedSharedMatrix::init warning: Out of Dcache, using operator new\n";
		}
	}
	n_cells_in_cache = 0;
}

long int n_loads = 0;
static long int n_nonhits = 0;

template <class T1, class T2>
TCachedSharedMatrix<T1,T2>::~TCachedSharedMatrix()
{
	if (Dcache.verboseflag()) {
#		ifdef _UNICOS
		static int n_loads_int, n_nonhits_int;
		n_loads_int = n_loads;
		n_nonhits_int = n_nonhits;
		shmemsumtoall(&n_loads_int);
		shmemsumtoall(&n_nonhits_int);
		n_loads = n_loads_int; n_nonhits = n_nonhits_int;
#		endif
		if (n_loads > 0 && mype==0) {
			cerr << "TCachedSharedMatrix: " << n_loads << " refs, " << n_loads-n_nonhits << " hits, "
				 << 100.0*(n_loads-n_nonhits)/double(n_loads) << " % hit rate\n";
		}
	}
	if (cache_was_allocated_using_new) {delete [] cache; delete [] indices_in_cache; cache = 0;}
#	if USE_SHMEM
	sharedfree(ptr);
#	else
	delete [] ptr;
#	endif
}

template <class T1, class T2>
TCacheIndex TCachedSharedMatrix<T1,T2>::load(TGridIndex cell) const
{
	assert(cell >= 0);
#	ifndef _UNICOS
	n_loads++;
#	endif
	// Check if it is already in cache
	TCacheIndex ci;
	for (ci=0; ci<n_cells_in_cache; ci++)
		if (indices_in_cache[ci] == cell)
#			if AVOID_MEMCPY_HACK
			return -ci-1;
#			else
			return ci;
#			endif
	if (n_cells_in_cache < cache_size) {
#if USE_LOCALITY_TESTS
		if (pe(cell) == mype) {
#			if AVOID_MEMCPY_HACK
#ifndef _UNICOS
			n_nonhits++;
#endif
			return cell;
#			else
			memcpy(&cache[n_cells_in_cache*ncols],&ptr[iL(cell)*ncols],ncols*sizeof(TCacheWord));
#			endif
		} else
#endif	/* USE_LOCALITY_TESTS */
			shmemget(&cache[n_cells_in_cache*ncols],&ptr[iL(cell)*ncols],ncols,pe(cell));
			
		indices_in_cache[n_cells_in_cache] = cell;	// mark the cache tag
#ifndef _UNICOS
		n_nonhits++;
#endif
#if HAVE_MUTABLE
#		if AVOID_MEMCPY_HACK
		return -(n_cells_in_cache++)-1;
#		else
		return n_cells_in_cache++;
#		endif
#else
		const smallnat retval = n_cells_in_cache;
		const smallnat newval = n_cells_in_cache + 1;
		*((smallnat *)(&n_cells_in_cache)) = newval;
#		if AVOID_MEMCPY_HACK
		return -retval-1;
#		else
		return retval;
#		endif
#endif
	} else {
		cerr << "*** TCachedSharedMatrix::load: Cache overflow. Increase cache_size (" << cache_size << ") in cshmat.C.\n";
		abort();
		return 0;
	}
}

#if HAS_TEMPLATE_INIT_SYNTAX
template class TCachedSharedMatrix<real,TGridIndex>;
#elif defined(__sgi)
#  pragma instantiate TCachedSharedMatrix<real,TGridIndex>
#endif
