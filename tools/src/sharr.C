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
#  pragma implementation "sharr.H"
#endif

#include "sharr.H"
#if USE_SHMEM
#  include <malloc.h>		// for shmalloc, shfree
#ifdef __sgi
#  define shmalloc_check(x)		// there is no shmalloc_check on IRIX, thus define it away
#endif
#endif
#include <stdlib.h>			// for memset

#if USE_SHMEM

const int shmalloc_padding = 0 /*10000*/;

// Wrappers for shmalloc, shfree, made extern since also called from cshmat.C
void *sharedmalloc(size_t sz)
{
	void *result;
	void *result1 = shmalloc(sz+2*shmalloc_padding*sizeof(int));
	result = (int *)result1 + shmalloc_padding;
	// We zero the result arrays. They have to be some legal floating point values
	// so that we can call applyBC() even before initializing the grid with values.
	memset(result1,0,sz+2*shmalloc_padding*sizeof(int));
	shmalloc_check(0);
	shmembarrierall();
	return result;
}

void sharedfree(void *ptr)
{
	shfree((int *)ptr - shmalloc_padding);
	shmalloc_check(0);
}
#endif

template <class T>
TSharedArray<T>::TSharedArray()
	: chunksize(1)
{
	ptr = 0;
	len = -1;
}

extern void doabort();

#if CHUNKSIZE_ALWAYS_ONE && !USE_CSHMAT
# error CHUNKSIZE_ALWAYS_ONE=1 and USE_CSHMAT=0 is not possible combination
#endif

template <class T>
void TSharedArray<T>::init(TGridIndex n, TGridIndex chunksize1)
{
#	if CHUNKSIZE_ALWAYS_ONE
	if (chunksize1 != 1) {
		cerr << "*** TSharedArray<T>::init: CHUNKSIZE_ALWAYS_ONE is defined but chunksize=" << chunksize1 << "\n";
		doabort();
	}
#	endif
	chunksize = chunksize1;
//	if (chunksize != 1) cerr << "################# CHUNKSIZE=" << chunksize << "\n";
#if USE_SHMEM
	if (ptr && len > 0) {
		cerr << "TSharedArray<T>::init warning: Calling shfree(ptr)\n";
		sharedfree(ptr);
	}
	const TGridIndex bs = chunksize*Npes;
	len = ((n+bs-1)/bs)*bs;
	locallen = DivNpes(len);
	ptr = (T *)sharedmalloc(locallen*sizeof(T));
	if (!ptr) {
		cerr << "*** TSharedArray<T>::init: shmalloc(" << locallen*sizeof(T) << ") returned NULL, Out of memory\n";
		cerr << "    n=" << n << ", chunksize=" << chunksize << ", Npes=" << Npes << ", len=" << len << "\n" << flush;
		doabort();
		exit(1);
	}
#else
	if (ptr && len > 0) delete [] ptr;
	ptr = new T [n];
	memset(ptr,0,n*sizeof(T));
	len = n;
#endif
}

template <class T>
void TSharedArray<T>::zero() {
	if (!ptr) {cerr << "*** TSharedArray<T>::zero: Array not yet inited\n"; return;}
#if USE_SHMEM
	memset(ptr,0,sizeof(T)*locallen);
	shmembarrierall();
#else
	memset(ptr,0,sizeof(T)*len);
#endif
}

template <class T>
TSharedArray<T>::~TSharedArray() {
#if USE_SHMEM
	shfree(ptr);
#else
	delete [] ptr;
#endif
}

#if USE_SHMEM

template <class T>
void TSharedArray<T>::gather(TGridIndex i0,const TGridIndexVector& iv, TGridIndex mult,
							 T* result, TGridIndex result_stride) const
{
	smallnat v;
	const int n = iv.length();
	if (sizeof(T) == 4) {
		static short indexbuff[VECLEN];
		static T resultbuff[VECLEN];
#		pragma _CRI cache_align indexbuff,resultbuff
		unroll4;
		for (v=0; v<n; v++) indexbuff[v] = (iv(v) > 0 ? iv(v) : 0)*mult;	// allow for iv(v)==NOINDEX without causing problems
		// check that all map to same PE
		const int pe_i0 = pe(i0+indexbuff[0]);
		for (v=0; v<n; v++) if (pe(i0+indexbuff[v]) != pe_i0) {cerr << "*** TSharedArray::gather32 error\n"; exit(1);}
		// localize indexbuff
		for (v=0; v<n; v++) indexbuff[v] = iL(i0+indexbuff[v]) - iL(i0);
		shmem_ixget32(resultbuff,ptr+iL(i0),indexbuff,n,pe_i0);
		unroll4;
		for (v=0; v<n; v++) result[v*result_stride] = resultbuff[v];
	} else {
		static long indexbuff[VECLEN];
		static T resultbuff[VECLEN];
#		pragma _CRI cache_align indexbuff,resultbuff
		unroll4;
		for (v=0; v<n; v++) indexbuff[v] = (iv(v) > 0 ? iv(v) : 0)*mult;	// allow for iv(v)==NOINDEX without causing problems
		// check that all map to same PE
		const int pe_i0 = pe(i0+indexbuff[0]);
		for (v=0; v<n; v++) if (pe(i0+indexbuff[v]) != pe_i0) {cerr << "*** TSharedArray::gather64 error\n"; exit(1);}
		// localize indexbuff
		for (v=0; v<n; v++) indexbuff[v] = iL(i0+indexbuff[v]) - iL(i0);
		shmem_ixget64(resultbuff,ptr+iL(i0),indexbuff,n,pe_i0);
		unroll4;
		for (v=0; v<n; v++) result[v*result_stride] = resultbuff[v];
	}
}

template <class T>
void TSharedArray<T>::scatter(TGridIndex i0, const TGridIndexVector& iv, TGridIndex mult,
							  const T* source, TGridIndex source_stride)
{
	smallnat v;
	const int n = iv.length();
	if (sizeof(T) == 4) {
		static short indexbuff[VECLEN];
		static T sourcebuff[VECLEN];
#		pragma _CRI cache_align indexbuff,sourcebuff
		unroll4;
		for (v=0; v<n; v++) indexbuff[v] = iv(v)*mult;
		// check that all map to same PE
		const int pe_i0 = pe(i0+indexbuff[0]);
		for (v=0; v<n; v++) if (pe(i0+indexbuff[v]) != pe_i0) {cerr << "*** TSharedArray::scatter32 error\n"; exit(1);}
		// localize indexbuff
		for (v=0; v<n; v++) indexbuff[v] = iL(i0+indexbuff[v]) - iL(i0);
		unroll4;
		for (v=0; v<n; v++) sourcebuff[v] = source[v*source_stride];
		shmem_ixput32(ptr+iL(i0),sourcebuff,indexbuff,n,pe_i0);
	} else {
		static long indexbuff[VECLEN];
		static T sourcebuff[VECLEN];
#		pragma _CRI cache_align indexbuff,sourcebuff
		unroll4;
		for (v=0; v<n; v++) indexbuff[v] = iv(v)*mult;
		// check that all map to same PE
		const int pe_i0 = pe(i0+indexbuff[0]);
		for (v=0; v<n; v++) if (pe(i0+indexbuff[v]) != pe_i0) {cerr << "*** TSharedArray::scatter64 error\n"; exit(1);}
		// localize indexbuff
		for (v=0; v<n; v++) indexbuff[v] = iL(i0+indexbuff[v]) - iL(i0);
		unroll4;
		for (v=0; v<n; v++) sourcebuff[v] = source[v*source_stride];
		shmem_ixput64(ptr+iL(i0),sourcebuff,indexbuff,n,pe_i0);
	}
}

#endif

#if HAS_TEMPLATE_INIT_SYNTAX
template class TSharedArray<TGridIndex>;
template class TSharedArray<real>;
template class TSharedArray<TVoidPtr>;
#elif defined(__sgi)
#  pragma instantiate TSharedArray<TGridIndex>
#  pragma instantiate TSharedArray<real>
#  pragma instantiate TSharedArray<TVoidPtr>
#endif
