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

/*

  Shared array class
  ------------------

  Initializer: Must know: Array length. (We support only 1D arrays here.(Do we, see below?))
  Must find out: Number of PEs.
  Must decide: Placement pattern. User-specifiable.

  i <--> (iL,pe) where i is the user-visible (global) index.

  Allocate memory using shmalloc on each PE.

  Provide methods to transfer data to/from the shared array
  and arrays in local memory. These methods must support nonunit stride
  for both source and destination (shmem_iget, shmem_iput routines).

  Unit stride in the user-visible global index does not necessarily translate
  to unit stride in the local index iL. It depends on the layout.
  
  Must also take into account the possibility that the access crosses
  PE boundary.

  Finally we must do both integer and real arrays, so make the shared array
  class into a template class.

UGLINESS: Must know the physical layout for the local buffer tables in HCgrid.C.
Now compiler takes care of this. Other solution would be to use a flat buffer table
in between.

In HCgrid.C, we need the pool matrix, not 1D array.
Anyway we must get rid of the & operators as in
&pool.M(0,1) - &pool.M(0,0);

since these would be global addresses.

 */

#ifndef SHARR_H

#ifdef __GNUC__
#  pragma interface
#endif
#include "realptr.H"
#include "shmemintf.H"

// Warning: If you define USE_CSHMAT=1 in makeflags.inc (for a non-SHMEM machine,
// for SHMEM_machines USE_CSHMAT is the default), you should also define
// CHUNKSIZE_ALWAYS_ONE=1 there.
#ifndef CHUNKSIZE_ALWAYS_ONE
#  if USE_SHMEM
#    define CHUNKSIZE_ALWAYS_ONE 1
#  else
#    define CHUNKSIZE_ALWAYS_ONE 0
#  endif
#endif

// () indexing works for arbitrary element.
// If you are sure that the index corresponds to local PE, you can use [] indexing,
// which is faster because it does not call SHMEM but accesses the element directly.
// Likewise, you can use lput() instead of put() if the address is local.
// Locality is checked (asserted) if DEBUG is defined.
// Index i is local if and only if (i % Npes) == mype

template <class T> class TSharedArray {
private:
	T *ptr;					// Pointer to the data vector. Symmetric pointer for SHMEM (via shmalloc)
	TGridIndex len;			// Total length of data vector (across all PEs). May be increased if USE_SHMEM.
	TGridIndex chunksize;	// chunksize consecutive data items reside on the same PE
#if USE_SHMEM
	TGridIndex locallen;		// Number of data items on one PE (always holds: len = locallen*Npes)
#	if CHUNKSIZE_ALWAYS_ONE
	int pe(TGridIndex i) const {return ModNpes(i);}
	TGridIndex iL(TGridIndex i) const {return DivNpes(i);}
	TGridIndex compose(TGridIndex iL, int pe) const {return iL*Npes + pe;}
#	else
	int pe(TGridIndex i) const {return ModNpes(i/chunksize);}
	TGridIndex iL(TGridIndex i) const {return (i/(Npes*chunksize))*chunksize + (i % chunksize);}
	TGridIndex compose(TGridIndex iL, int pe) const {return (iL/chunksize)*(Npes*chunksize) + (iL % chunksize) + chunksize*pe;}
#	endif
#endif
public:
	TSharedArray();
	bool isempty() const {return ptr==0;}
	void init(TGridIndex n, TGridIndex chunksize1=1);
	void zero();		// fill the entire array with zeros, init must be called first, of course
	TSharedArray(TGridIndex n, TGridIndex chunksize1=1) {ptr = 0; len = -1; init(n,chunksize1);}
#if USE_SHMEM
	T operator()(TGridIndex i) const {return shmemget(ptr+iL(i),pe(i));}
	T operator[](TGridIndex i) const {assert(pe(i)==mype); return ptr[iL(i)];}
	void put(TGridIndex i, T val) {shmemput(ptr+iL(i),val,pe(i));}
	void lput(TGridIndex i, T val) {assert(pe(i)==mype); ptr[iL(i)] = val;}
	T atomically_read_and_set(TGridIndex i, T val) {return shmemswap(ptr+iL(i),val,pe(i));}
#else
	T operator()(TGridIndex i) const {return ptr[i];}
	T operator[](TGridIndex i) const {return ptr[i];}
	void put(TGridIndex i, T val) {ptr[i] = val;}
	void lput(TGridIndex i, T val) {ptr[i] = val;}
	T atomically_read_and_set(TGridIndex i, T val) {const T result = ptr[i]; ptr[i] = val; return result;}
#endif
	// transfer of L-length block starting from i (unit stride) to/from array
	void get(TGridIndex i, TGridIndex L, T* result, TGridIndex result_stride=1) const;
	void put(TGridIndex i, TGridIndex L, const T* source, TGridIndex source_stride=1);
	// transfer of L-length block starting from i (arbitrary stride) to/from array
	void get(TGridIndex i, TGridIndex L, TGridIndex stride, T* result, TGridIndex result_stride=1) const;
	void put(TGridIndex i, TGridIndex L, TGridIndex stride, const T* source, TGridIndex source_stride=1);
	// gather/scatter
	void gather(TGridIndex i0, const TGridIndexVector& iv, TGridIndex mult, T* result, TGridIndex result_stride=1) const;
	void scatter(TGridIndex i0, const TGridIndexVector& iv, TGridIndex mult, const T* source, TGridIndex source_stride=1);
	TGridIndex length() const {return len;}
	~TSharedArray();
};

#if USE_SHMEM

template <class T>
inline void TSharedArray<T>::get(TGridIndex i, TGridIndex L, TGridIndex stride, T* result, TGridIndex result_stride) const
{
	// Check that all indices i:i+L-1 map to the same PE
	const int pe_i = pe(i);
	TGridIndex j;
	for (j=i; j<i+L; j++) if (pe(j) != pe_i) {cerr << "*** TSharedArray::get error!\n"; exit(1);}
	shmemiget(result,ptr+iL(i),result_stride,stride,L,pe(i));
}

template <class T>
inline void TSharedArray<T>::put(TGridIndex i, TGridIndex L, TGridIndex stride, const T* source, TGridIndex source_stride)
{
	// Check that all indices i:i+L-1 map to the same PE
	const int pe_i = pe(i);
	TGridIndex j;
	for (j=i; j<i+L; j++) if (pe(j) != pe_i) {cerr << "*** TSharedArray::put error!\n"; exit(1);}
	shmemiput(ptr+iL(i),source,stride,source_stride,L,pe(i));
}

#else

template <class T>
inline void TSharedArray<T>::get(TGridIndex i, TGridIndex L, TGridIndex stride, T* result, TGridIndex result_stride) const
{
	TGridIndex p;
	VDIRS; unroll4;
	for (p=0; p<L; p++) result[p*result_stride] = ptr[i+p*stride];
}

template <class T>
inline void TSharedArray<T>::put(TGridIndex i, TGridIndex L, TGridIndex stride, const T* source, TGridIndex source_stride)
{
	TGridIndex p;
	VDIRS; unroll4;
	for (p=0; p<L; p++) ptr[i+p*stride] = source[p*source_stride];
}

template <class T>
inline void TSharedArray<T>::gather(TGridIndex i0,const TGridIndexVector& iv, TGridIndex mult,
									T* result, TGridIndex result_stride) const
{
	smallnat v;
	const smallnat n = iv.length();
	TGridIndex jd;
	VDIRS; unroll4;
	for (v=0,jd=0; v<n; v++,jd+=result_stride) result[jd] = ptr[i0+iv(v)*mult];
}

template <class T>
inline void TSharedArray<T>::scatter(TGridIndex i0, const TGridIndexVector& iv, TGridIndex mult,
									 const T* source, TGridIndex source_stride)
{
	smallnat v;
	const smallnat n = iv.length();
	TGridIndex js;
	VDIRS; unroll4;
	for (v=0,js=0; v<n; v++,js+=source_stride) ptr[i0+iv(v)*mult] = source[js];
}

#endif

template <class T>
inline void TSharedArray<T>::get(TGridIndex i, TGridIndex L, T* result, TGridIndex result_stride) const
{
	get(i,L,1,result,result_stride);
}

template <class T>
inline void TSharedArray<T>::put(TGridIndex i, TGridIndex L, const T* source, TGridIndex source_stride)
{
	put(i,L,1,source,source_stride);
}

// Useful typedefs. Also, gcc compiler has sometimes problems without.
typedef TSharedArray<TGridIndex> TIndexTable;
typedef void *TVoidPtr;

#define SHARR_H

#endif
