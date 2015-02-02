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

#ifndef CSHMAT_H

#ifdef __GNUC__
#  pragma interface
#endif

#include "realptr.H"
#include "cache.H"
#include "shmemintf.H"
#include <stdlib.h>		// for abort()

#ifndef AVOID_MEMCPY_HACK
#  define AVOID_MEMCPY_HACK 1
#endif

#if USE_SHMEM
#  ifndef USE_LOCALITY_TESTS
#    if SIMULATE_SLOW_SHMEM
#      define USE_LOCALITY_TESTS 0
#    else
#      define USE_LOCALITY_TESTS 1
#    endif
#  endif
#else
// If we are not using SHMEM, turn on locality tests in all cases:
// (otherwise doesn't always work (?))
#  define USE_LOCALITY_TESTS 1
#endif

/*
 * If AVOID_MEMCPY_HACK is defined as 1, TCacheIndex is a signed quantity
 * (TGridIndex) encoded as follows: Values 0,1,2,... correspond to direct local memory
 * cell indices (they are not in the cache). Values ...,-2,-1 map to values 0,1,... in the
 * cache (by formula cache_index = -n-1).
 */

#if AVOID_MEMCPY_HACK
#  define TCacheIndex TGridIndex
#else
   typedef smallnat TCacheIndex;
#endif
typedef smallnat TCompIndex;

#if defined(__GNUC__) || defined(_UNICOS) || defined(__DECCXX)
#  define HAVE_MUTABLE 1		// Whether our C++ has the 'mutable' keyword. Mutable members can be changed by const functions.
#else
#  define HAVE_MUTABLE 0		// If we don't have 'mutable', we have to use a backdoor way (pointer type casts)
#endif

/*
 * A TCachedSharedMatrix<T1,T2>(nrows,ncols_1,ncols_2) is actually two matrices, one with
 * component type T1 (nrows x ncols_1), the other one with component type T2 (nrows x ncols_2).
 * 
 */

template <class T1, class T2>
class TCachedSharedMatrix {
private:
	TCacheWord *ptr;
	TGridIndex nrows;
	TCompIndex ncols1,ncols2,ncols;
	smallnat cache_size;		// constant, set in init(), currently equal to 31
	TCacheWord *cache;			// cache_size*ncols TCacheWord items
	TGridIndex *indices_in_cache;	// cache_size TGridIndex items
#if HAVE_MUTABLE
	mutable
#endif
	smallnat n_cells_in_cache;
	bool cache_was_allocated_using_new;
#if USE_SHMEM
	int pe(TGridIndex cell) const {return ModNpes(cell);}
#else
	int pe(TGridIndex) const {return 0;}
#endif
	TGridIndex iL(TGridIndex cell) const {return DivNpes(cell);}
public:
	TCachedSharedMatrix();
	void init(TGridIndex nrows1, TCompIndex ncols_1, TCompIndex ncols_2);
	TCachedSharedMatrix(TGridIndex nrows1, TCompIndex ncols_1, TCompIndex ncols_2) {ptr=0; init(nrows1,ncols_1,ncols_2);}
	void clearcache() const {
#if HAVE_MUTABLE
		n_cells_in_cache = 0;
#else
		*((smallnat *)(&n_cells_in_cache)) = 0;
#endif
	}
	TCacheIndex load(TGridIndex cell) const;
	void popcache() const {
		assert(n_cells_in_cache > 0);
#		if HAVE_MUTABLE
		n_cells_in_cache--;
#		else
		*((smallnat *)(&n_cells_in_cache)) = n_cells_in_cache - 1;
#		endif
	}

	T1 M1C(TCacheIndex cell, TCompIndex comp) const {
#		if AVOID_MEMCPY_HACK
		return *((T1*)(cell < 0 ? &cache[(-cell-1)*ncols] : &ptr[iL(cell)*ncols]) + comp);
#		else
		return *((T1*)(&cache[cell*ncols]) + comp);
#		endif
	}
	T1 M1(TGridIndex cell, TCompIndex comp) const {
		assert(cell >= 0);
#if !SIMULATE_SLOW_SHMEM || USE_LOCALITY_TESTS
		if (pe(cell) == mype)
			return *((T1*)(&ptr[iL(cell)*ncols]) + comp);
		else
#endif
			return shmemget((T1*)(&ptr[iL(cell)*ncols]) + comp, pe(cell));
	}
	T1 M1local(TGridIndex cell, TCompIndex comp) const {
		assert(cell >= 0);
		assert(pe(cell) == mype);
		return *((T1*)(&ptr[iL(cell)*ncols]) + comp);
	}
	void M1set(TGridIndex cell, TGridIndex comp, T1 val) {
		assert(cell >= 0);
#		if USE_SHMEM
		shmemput((T1*)(&ptr[iL(cell)*ncols]) + comp, val, pe(cell));
#		else
		*((T1*)(&ptr[cell*ncols]) + comp) = val;
#		endif
	}
	void M1setlocal(TGridIndex cell, TGridIndex comp, T1 val) {
		assert(pe(cell) == mype);
		*((T1*)(&ptr[iL(cell)*ncols]) + comp) = val;
	}
	void M1set_multi(TGridIndex cell, smallnat comp0, smallnat ncomps, const T1 src[], int src_stride) {
#		if USE_SHMEM
		shmemiput((T1*)(&ptr[iL(cell)*ncols]) + comp0, src, 1, src_stride, ncomps, pe(cell));
#		else
		smallnat comp,q;
		for (comp=comp0,q=0; q<ncomps; comp++,q++)
			*((T1*)(&ptr[cell*ncols]) + comp) = src[q*src_stride];
#		endif
	}

	T2 M2C(TCacheIndex cell, TGridIndex comp) const {
#		if AVOID_MEMCPY_HACK
		return *((T2*)(cell < 0 ? &cache[(-cell-1)*ncols + ncols1] : &ptr[iL(cell)*ncols + ncols1]) + comp);
#		else
		return *((T2*)(&cache[cell*ncols + ncols1]) + comp);
#		endif
	}
	T2 M2(TGridIndex cell, TGridIndex comp) const {
		assert(cell >= 0);
#if !SIMULATE_SLOW_SHMEM || USE_LOCALITY_TESTS
		if (pe(cell) == mype)
			return *((T2*)(&ptr[iL(cell)*ncols + ncols1]) + comp);
		else
#endif
			return shmemget((T2*)(&ptr[iL(cell)*ncols + ncols1]) + comp, pe(cell));
	}
	T2 M2local(TGridIndex cell, TGridIndex comp) const {
		assert(cell >= 0);
		assert(pe(cell) == mype);
		return *((T2*)(&ptr[iL(cell)*ncols + ncols1]) + comp);
	}
	void M2set(TGridIndex cell, TGridIndex comp, T2 val) {
		assert(cell >= 0);
#		if USE_SHMEM
		shmemput((T2*)(&ptr[iL(cell)*ncols + ncols1]) + comp, val, pe(cell));
#		else
		*((T2*)(&ptr[cell*ncols + ncols1]) + comp) = val;
#		endif
	}
	void M2setlocal(TGridIndex cell, TGridIndex comp, T2 val) {
		assert(cell >= 0);
		assert(pe(cell) == mype);
		*((T2*)(&ptr[iL(cell)*ncols + ncols1]) + comp) = val;
	}
	T2 M2atomically_read_and_set(TGridIndex cell, TGridIndex comp, T2 val) {
		assert(cell >= 0);
		return shmemswap((T2*)(&ptr[iL(cell)*ncols + ncols1]) + comp, val, pe(cell));
	}
#	if USE_SHMEM
	void M2setlock(TGridIndex cell, TGridIndex comp) {
		assert(cell >= 0);
		int cnt = 0;
		while (shmemswap((T2*)(&ptr[iL(cell)*ncols + ncols1]) + comp, T2(1), pe(cell)) && cnt < 1000000) cnt++;
		if (cnt >= 1000000) {cerr << "*** M2setlock: deadlock\n"; abort();}
	}
	void M2unsetlock(TGridIndex cell, TGridIndex comp) {
		assert(cell >= 0);
		int cnt = 0;
		if (shmemswap((T2*)(&ptr[iL(cell)*ncols + ncols1]) + comp, T2(0), pe(cell)) != T2(1)) {
			cerr << "*** M2unsetlock: old value was not 1\n";
			abort();
		}
	}
#	endif
	~TCachedSharedMatrix();
};

#define CSHMAT_H

#endif
