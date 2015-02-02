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

#ifndef SHMEMINTF_H

#ifdef __GNUC__
#  pragma interface
#endif

#include <cstring>		// to get memcpy
#include <iostream>
using namespace std;

#ifndef USE_SHMEM
#define USE_SHMEM 0
#endif

#ifndef BARRIER_DEBUG
#  ifdef DEBUG
#    define BARRIER_DEBUG 1
#  else
#    define BARRIER_DEBUG 0
#  endif
#endif

#define ROOT_PE 0

#if USE_SHMEM

#define real real_hack
#ifdef __sgi
#  include <mpp/shmem.h>
#else
/*
 * For Cray T3E, starting from 3.1.0.0 compiler version (see CC -V), the system include
 * file mpp/shmem.h seems to be C++ compatible. Thus we don't include our own version shmem.H
 * any more. If you have an CC version, you may have to include shmem.H instead, but be aware
 * that it is potentially unsafe. The safest way is to copy mpp/shmem.h to shmem.H and to fix
 * the few formal things associated with complex numbers that cause an incompatibility with C++.
 */
//#  include "shmem.H"
#include <mpp/shmem.h>
#endif
#undef real

#ifdef FIXED_LOG2_NPES

   extern short int mype;
#  define log2_Npes FIXED_LOG2_NPES
#  define Npes (1 << FIXED_LOG2_NPES)
#  define Npes_minus_1 (Npes-1)
#  define DivNpes(x) ((x) >> FIXED_LOG2_NPES)
#  define ModNpes(x) ((x) & Npes_minus_1)

#else

   extern short int mype, Npes;		// making them short int is slightly faster

#  ifndef NPES_ONLY_POWER_OF_TWO
#    define NPES_ONLY_POWER_OF_TWO 0
#  endif

#  if NPES_ONLY_POWER_OF_TWO
     extern short int log2_Npes, Npes_minus_1;
#    define DivNpes(x) ((x) >> log2_Npes)
#    define ModNpes(x) ((x) & Npes_minus_1)
#  else
#    define DivNpes(x) ((x)/Npes)
#    define ModNpes(x) ((x)%Npes)
#  endif

#endif

#include <sstream>
#define ATOMIC_OUTPUT(code) {        \
      strstream o;                   \
	  o << "PE " << mype << ": ";    \
      code;                          \
      o << ends;                     \
      char *buff = o.str();          \
      cerr << buff;                  \
      delete [] buff;                \
   }

inline int shmemmype() {return shmem_my_pe();}
inline int shmemnpes() {return shmem_n_pes();}

#ifdef __sgi

#if BARRIER_DEBUG
#  define shmembarrierall() ATOMIC_OUTPUT(o << "Barrier call from " << __FILE__ << ":" << __LINE__ << "\n"); shmem_barrier_all()
#else
#  define shmembarrierall() shmem_barrier_all()
#endif

#else

extern "C" void barrier(void);
#if BARRIER_DEBUG
#  define shmembarrierall() ATOMIC_OUTPUT(o << "Barrier call from " << __FILE__ << ":" << __LINE__ << "\n"); barrier()
#else
#  define shmembarrierall() barrier()
#endif

#endif

extern void startpes();

inline void shmemquiet() {shmem_quiet();}

inline float shmemget(const float *src, int pe) {return shmem_float_g(src,pe);}
inline double shmemget(const double *src, int pe) {return shmem_double_g(src,pe);}
inline short shmemget(const short *src, int pe) {return shmem_short_g(src,pe);}
inline int shmemget(const int *src, int pe) {return shmem_int_g(src,pe);}
inline long shmemget(const long *src, int pe) {return shmem_long_g(src,pe);}
inline void *shmemget(void *const *src, int pe) {return (void *)shmem_long_g((long *)src,pe);}

inline void shmemget(float *target, const float *src, int len, int pe) {shmem_float_get(target,src,len,pe);}
inline void shmemget(double *target, const double *src, int len, int pe) {shmem_double_get(target,src,len,pe);}
inline void shmemget(short *target, const short *src, int len, int pe) {shmem_short_get(target,src,len,pe);}
inline void shmemget(int *target, const int *src, int len, int pe) {shmem_int_get(target,src,len,pe);}
inline void shmemget(long *target, const long *src, int len, int pe) {shmem_long_get(target,src,len,pe);}

inline void shmemput(float *target, float val, int pe) {shmem_float_p(target,val,pe);}
inline void shmemput(double *target, double val, int pe) {shmem_double_p(target,val,pe);}
inline void shmemput(short *target, short val, int pe) {shmem_short_p(target,val,pe);}
inline void shmemput(int *target, int val, int pe) {shmem_int_p(target,val,pe);}
inline void shmemput(long *target, long val, int pe) {shmem_long_p(target,val,pe);}
inline void shmemput(void **target, void *val, int pe) {shmem_long_p((long*)target,long(val),pe);}

inline void shmemput(float *target, const float *src, int len, int pe) {shmem_float_put(target,src,len,pe);}
inline void shmemput(double *target, const double *src, int len, int pe) {shmem_double_put(target,src,len,pe);}
inline void shmemput(short *target, const short *src, int len, int pe) {shmem_short_put(target,src,len,pe);}
inline void shmemput(int *target, const int *src, int len, int pe) {shmem_int_put(target,src,len,pe);}
inline void shmemput(long *target, const long *src, int len, int pe) {shmem_long_put(target,src,len,pe);}

/*
 * SGI doesn't seem to have shmem_float_iget, shmem_double_iget, etc., even though the man page lists them.
 * The only environment supported for libsma is -64, so we can assume that we are in 64 bit mode.
 * In that mode, sizeof(short)=2, sizeof(int)=4, sizeof(long)=8. Thus, float and int correspond to 32
 * while double and long correspond to 64. This leaves us with the problem of what to do with 'short int'.
 * We implement shmemiget(short*,..) and shmemiput(short*,..) for __sgi by explicit looping. This is slow,
 * but it should work.
 */

#ifdef __sgi
#define shmem_float_iget shmem_iget32
#define shmem_double_iget shmem_iget64
#define shmem_int_iget shmem_iget32
#define shmem_long_iget shmem_iget64
#define shmem_float_iput shmem_iput32
#define shmem_double_iput shmem_iput64
#define shmem_int_iput shmem_iput32
#define shmem_long_iput shmem_iput64
#endif

inline void shmemiget(float *target, const float *src, int target_stride, int src_stride, int len, int pe) {
	shmem_float_iget(target,src,target_stride,src_stride,len,pe);
}
inline void shmemiget(double *target, const double *src, int target_stride, int src_stride, int len, int pe) {
	shmem_double_iget(target,src,target_stride,src_stride,len,pe);
}
inline void shmemiget(short *target, const short *src, int target_stride, int src_stride, int len, int pe) {
#	ifdef __sgi
	int i;
	for (i=0; i<len; i++) target[i*target_stride] = shmemget(src+i*src_stride, pe);
#	else
	shmem_short_iget(target,src,target_stride,src_stride,len,pe);
#	endif
}
inline void shmemiget(int *target, const int *src, int target_stride, int src_stride, int len, int pe) {
	shmem_int_iget(target,src,target_stride,src_stride,len,pe);
}
inline void shmemiget(long *target, const long *src, int target_stride, int src_stride, int len, int pe) {
	shmem_long_iget(target,src,target_stride,src_stride,len,pe);
}
inline void shmemiget(void **target, void *const *src, int target_stride, int src_stride, int len, int pe) {
	shmem_long_iget((long*)target,(const long*)src,target_stride,src_stride,len,pe);
}

inline void shmemiput(float *target, const float *src, int target_stride, int src_stride, int len, int pe) {
	shmem_float_iput(target,src,target_stride,src_stride,len,pe);
}
inline void shmemiput(double *target, const double *src, int target_stride, int src_stride, int len, int pe) {
	shmem_double_iput(target,src,target_stride,src_stride,len,pe);
}
inline void shmemiput(short *target, const short *src, int target_stride, int src_stride, int len, int pe) {
	shmem_short_iput(target,src,target_stride,src_stride,len,pe);
}
inline void shmemiput(int *target, const int *src, int target_stride, int src_stride, int len, int pe) {
	shmem_int_iput(target,src,target_stride,src_stride,len,pe);
}
inline void shmemiput(long *target, const long *src, int target_stride, int src_stride, int len, int pe) {
	shmem_long_iput(target,src,target_stride,src_stride,len,pe);
}
inline void shmemiput(void **target, void *const *src, int target_stride, int src_stride, int len, int pe) {
	shmem_long_iput((long*)target,(const long*)src, target_stride,src_stride,len,pe);
}

inline short shmemfinc(short *target, int pe) {return shmem_short_finc(target,pe);}

inline void shmemsumtoall(short*trg, short*src, int sz, int pestart, int log_pestride, int pesize, short *pwrk, long *psync)
{
	shmem_short_sum_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}

inline void shmemsumtoall(int*trg, int*src, int sz, int pestart, int log_pestride, int pesize, int *pwrk, long *psync)
{
	shmem_int_sum_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}

inline void shmemsumtoall(float*trg, float*src, int sz, int pestart, int log_pestride, int pesize, float *pwrk, long *psync)
{
	shmem_float_sum_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}

inline void shmemsumtoall(double*trg, double*src, int sz, int pestart, int log_pestride, int pesize, double *pwrk, long *psync)
{
	shmem_double_sum_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}


inline void shmemmaxtoall(short*trg, short*src, int sz, int pestart, int log_pestride, int pesize, short *pwrk, long *psync)
{
	shmem_short_max_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}

// On IRIX we may require also shmemmaxtoall for smallnat type, which is unsigned short.
// This assumes that no smallnat uses the highest bit, which would make the corresponding signed short value negative
inline void shmemmaxtoall(unsigned short *trg, unsigned short *src,
						  int sz, int pestart, int log_pestride, int pesize, unsigned short *pwrk, long *psync)
{
	shmem_short_max_to_all((short *)trg,(short *)src, sz, pestart,log_pestride,pesize, (short *)pwrk,psync);
}

inline void shmemmaxtoall(int*trg, int*src, int sz, int pestart, int log_pestride, int pesize, int *pwrk, long *psync)
{
	shmem_int_max_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}

inline void shmemmaxtoall(float*trg, float*src, int sz, int pestart, int log_pestride, int pesize, float *pwrk, long *psync)
{
	shmem_float_max_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}

inline void shmemmaxtoall(double*trg, double*src, int sz, int pestart, int log_pestride, int pesize, double *pwrk, long *psync)
{
	shmem_double_max_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}


inline void shmemmintoall(short*trg, short*src, int sz, int pestart, int log_pestride, int pesize, short *pwrk, long *psync)
{
	shmem_short_min_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}

inline void shmemmintoall(int*trg, int*src, int sz, int pestart, int log_pestride, int pesize, int *pwrk, long *psync)
{
	shmem_int_min_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}

inline void shmemmintoall(float*trg, float*src, int sz, int pestart, int log_pestride, int pesize, float *pwrk, long *psync)
{
	shmem_float_min_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}

inline void shmemmintoall(double*trg, double*src, int sz, int pestart, int log_pestride, int pesize, double *pwrk, long *psync)
{
	shmem_double_min_to_all(trg,src, sz, pestart,log_pestride,pesize, pwrk,psync);
}

inline void shmemsetlock(long *lockptr) {shmem_set_lock(lockptr);}
inline void shmemclearlock(long *lockptr) {shmem_clear_lock(lockptr);}

template <class T>
inline void shmemsumtoall(T *ptr)
// Replace *ptr by the all-PE-wide sum of *ptr 's. ptr must be symmetric pointer (same on all PEs)
{
	static T pWrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE+1];
	static long pSync[_SHMEM_REDUCE_SYNC_SIZE];
	static T glob;
	pragma("_CRI cache_align pWrk,pSync");
	int q;
	for (q=0; q<_SHMEM_REDUCE_SYNC_SIZE; q++) pSync[q] = _SHMEM_SYNC_VALUE;
	shmembarrierall();		// barrier is always needed before reduction routine
	shmemsumtoall(&glob,ptr,1, 0,0,Npes, pWrk,pSync);
	*ptr = glob;	// copy result to *ptr
}

template <class T>
inline void shmemmaxtoall(T *ptr)
// Replace *ptr by the all-PE-wide maximum of *ptr 's. ptr must be symmetric pointer (same on all PEs)
{
	static T pWrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE+1];
	static long pSync[_SHMEM_REDUCE_SYNC_SIZE];
	static T glob;
	pragma("_CRI cache_align pWrk,pSync");
	int q;
	for (q=0; q<_SHMEM_REDUCE_SYNC_SIZE; q++) pSync[q] = _SHMEM_SYNC_VALUE;
	shmembarrierall();		// barrier is always needed before reduction routine
	shmemmaxtoall(&glob,ptr,1, 0,0,Npes, pWrk,pSync);
	*ptr = glob;	// copy result to *ptr
}

template <class T>
inline void shmemmintoall(T *ptr)
// Replace *ptr by the all-PE-wide minimum of *ptr 's. ptr must be symmetric pointer (same on all PEs)
{
	static T pWrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE+1];
	static long pSync[_SHMEM_REDUCE_SYNC_SIZE];
	static T glob;
	pragma("_CRI cache_align pWrk,pSync");
	int q;
	for (q=0; q<_SHMEM_REDUCE_SYNC_SIZE; q++) pSync[q] = _SHMEM_SYNC_VALUE;
	shmembarrierall();		// barrier is always needed before reduction routine
	shmemmintoall(&glob,ptr,1, 0,0,Npes, pWrk,pSync);
	*ptr = glob;	// copy result to *ptr
}

inline short shmemswap(short* target, short value, int pe) {return shmem_short_swap(target,value,pe);}
inline int shmemswap(int* target, int value, int pe) {return shmem_int_swap(target,value,pe);}
inline long shmemswap(long* target, long value, int pe) {return shmem_long_swap(target,value,pe);}
inline float shmemswap(float *target, float value, int pe) {return shmem_float_swap(target,value,pe);}
inline double shmemswap(double *target, double value, int pe) {return shmem_double_swap(target,value,pe);}
inline void *shmemswap(void **target, void *value, int pe) {return (void *)shmem_long_swap((long *)target,(long)value,pe);}

#else	/* USE_SHMEM == 0 */

#define mype 0
#define Npes 1
#define ATOMIC_OUTPUT(code) {ostream& o = cerr; code;}

#define DivNpes(x) (x)
#define ModNpes(x) 0

inline int shmemnpes() {return 1;}
inline int shmemmype() {return 0;}

#if BARRIER_DEBUG
#  define shmembarrierall() ATOMIC_OUTPUT(o << "Barrier call from " << __FILE__ << ":" << __LINE__ << "\n")
#else
#  define shmembarrierall()
#endif

inline void startpes() {}

#ifndef SIMULATE_SLOW_SHMEM
#  define SIMULATE_SLOW_SHMEM 0
#endif

#if SIMULATE_SLOW_SHMEM

extern float shmemget(const float *src, int pe);
extern double shmemget(const double *src, int pe);
extern short shmemget(const short *src, int pe);
extern int shmemget(const int *src, int pe);
extern long shmemget(const long *src, int pe);

extern void shmemget(float *target, const float *src, int len, int pe);
extern void shmemget(double *target, const double *src, int len, int pe);
extern void shmemget(short *target, const short *src, int len, int pe);
extern void shmemget(int *target, const int *src, int len, int pe);
extern void shmemget(long *target, const long *src, int len, int pe);

extern void shmemiget(float *target, const float *src, int target_stride, int src_stride, int len, int pe);
extern void shmemiget(double *target, const double *src, int target_stride, int src_stride, int len, int pe);
extern void shmemiget(short *target, const short *src, int target_stride, int src_stride, int len, int pe);
extern void shmemiget(int *target, const int *src, int target_stride, int src_stride, int len, int pe);
extern void shmemiget(long *target, const long *src, int target_stride, int src_stride, int len, int pe);

#else

inline void shmemquiet() {}

inline float shmemget(const float *src, int) {return *src;}
inline double shmemget(const double *src, int) {return *src;}
inline short shmemget(const short *src, int) {return *src;}
inline int shmemget(const int *src, int) {return *src;}
inline long shmemget(const long *src, int) {return *src;}

inline void shmemget(float *target, const float *src, int len, int) {
	int i;
#	ifdef _UNICOS
#	pragma _CRI unroll 4
#	endif
	for (i=0; i<len; i++) target[i] = src[i];
}

inline void shmemget(double *target, const double *src, int len, int) {
	int i;
#	ifdef _UNICOS
#	pragma _CRI unroll 4
//#	pragma _CRI cache_bypass src
#	endif
	for (i=0; i<len; i++) target[i] = src[i];
}

inline void shmemget(short *target, const short *src, int len, int) {
	int i;
#	ifdef _UNICOS
#	pragma _CRI unroll 4
#	endif
	for (i=0; i<len; i++) target[i] = src[i];
}

inline void shmemget(int *target, const int *src, int len, int) {
	int i;
#	ifdef _UNICOS
#	pragma _CRI unroll 4
//#	pragma _CRI cache_bypass src
#	endif
	for (i=0; i<len; i++) target[i] = src[i];
}

inline void shmemget(long *target, const long *src, int len, int) {
	int i;
#	ifdef _UNICOS
#	pragma _CRI unroll 4
//#	pragma _CRI cache_bypass src
#	endif
	for (i=0; i<len; i++) target[i] = src[i];
}

#endif

extern void shmemiget(float *target, const float *src, int target_stride, int src_stride, int len, int pe);
extern void shmemiget(double *target, const double *src, int target_stride, int src_stride, int len, int pe);
extern void shmemiget(short *target, const short *src, int target_stride, int src_stride, int len, int pe);
extern void shmemiget(int *target, const int *src, int target_stride, int src_stride, int len, int pe);
extern void shmemiget(long *target, const long *src, int target_stride, int src_stride, int len, int pe);

inline void shmemput(float *target, float val, int) {*target = val;}
inline void shmemput(double *target, double val, int) {*target = val;}
inline void shmemput(short *target, short val, int) {*target = val;}
inline void shmemput(int *target, int val, int) {*target = val;}
inline void shmemput(long *target, long val, int) {*target = val;}

inline void shmemput(float *target, const float *src, int len, int) {memcpy(target,src,len*sizeof(float));}
inline void shmemput(double *target, const double *src, int len, int) {memcpy(target,src,len*sizeof(double));}
inline void shmemput(short *target, const short *src, int len, int) {memcpy(target,src,len*sizeof(short));}
inline void shmemput(int *target, const int *src, int len, int) {memcpy(target,src,len*sizeof(int));}
inline void shmemput(long *target, const long *src, int len, int) {memcpy(target,src,len*sizeof(long));}

inline void shmemiput(float *target, const float *src, int target_stride, int src_stride, int len, int) {
	int i;
	for (i=0; i<len; i++) target[i*target_stride] = src[i*src_stride];
}

inline void shmemiput(double *target, const double *src, int target_stride, int src_stride, int len, int) {
	int i;
	for (i=0; i<len; i++) target[i*target_stride] = src[i*src_stride];
}

inline void shmemiput(short *target, const short *src, int target_stride, int src_stride, int len, int) {
	int i;
	for (i=0; i<len; i++) target[i*target_stride] = src[i*src_stride];
}

inline void shmemiput(int *target, const int *src, int target_stride, int src_stride, int len, int) {
	int i;
	for (i=0; i<len; i++) target[i*target_stride] = src[i*src_stride];
}

inline void shmemiput(long *target, const long *src, int target_stride, int src_stride, int len, int) {
	int i;
	for (i=0; i<len; i++) target[i*target_stride] = src[i*src_stride];
}

inline short shmemfinc(short *target, int) {return (*target)++;}

inline void shmemsumtoall(short*trg, short*src, int, int,int,int, short *, long *) {*trg=*src;}
inline void shmemsumtoall(int*trg, int*src, int, int,int,int, int *, long *) {*trg=*src;}
inline void shmemsumtoall(float*trg, float*src, int, int,int,int, float *, long *) {*trg=*src;}
inline void shmemsumtoall(double*trg, double*src, int, int,int,int, double *, long *) {*trg=*src;}

inline void shmemmaxtoall(short*trg, short*src, int, int,int,int, short *, long *) {*trg=*src;}
inline void shmemmaxtoall(int*trg, int*src, int, int,int,int, int *, long *) {*trg=*src;}
inline void shmemmaxtoall(float*trg, float*src, int, int,int,int, float *, long *) {*trg=*src;}
inline void shmemmaxtoall(double*trg, double*src, int, int,int,int, double *, long *) {*trg=*src;}

inline void shmemmintoall(short*trg, short*src, int, int,int,int, short *, long *) {*trg=*src;}
inline void shmemmintoall(int*trg, int*src, int, int,int,int, int *, long *) {*trg=*src;}
inline void shmemmintoall(float*trg, float*src, int, int,int,int, float *, long *) {*trg=*src;}
inline void shmemmintoall(double*trg, double*src, int, int,int,int, double *, long *) {*trg=*src;}

inline short shmemswap(short *target, short val, int) {const short result = *target; *target = val; return result;}
inline int shmemswap(int *target, int val, int) {const int result = *target; *target = val; return result;}
inline long shmemswap(long *target, long val, int) {const long result = *target; *target = val; return result;}
inline float shmemswap(float *target, float val, int) {const float result = *target; *target = val; return result;}
inline double shmemswap(double *target, double val, int) {const double result = *target; *target = val; return result;}
inline void *shmemswap(void **target, void *val, int) {void *const result = *target; *target = val; return result;}

inline void shmemsetlock(long *) {}
inline void shmemclearlock(long *) {}

#define shmemsumtoall(ptr)
#define shmemmaxtoall(ptr)
#define shmemmintoall(ptr)

#endif	/* USE_SHMEM */

class TBarrierOnFunctionReturn {
public:
	TBarrierOnFunctionReturn() {}
	~TBarrierOnFunctionReturn() {shmembarrierall();}
};

/* Usage: Write SINGLE_PE; as a first thing in a function that must
 * be executed by a single PE only (PE 0). Other PEs wait for PE 0 to finish
 * the function.
 *
 * Implementation: SINGLE_PE declares a special object whose destructor
 * barriers. Then, if the current PE is NOT PE 0, it returns immediately.
 * This return will call the destructor, i.e. barrier. If the current PE IS PE 0,
 * the function body is executed, and a barrier is encountered upon return.
 *
 * More advanced usage:
 * SINGLE_PE need not be the first statement. Statements before it
 * are executed by all processors.
 *
 * Using SINGLE_PE in an inner {} block is confusing and not recommended.
 */

#define SINGLE_PE TBarrierOnFunctionReturn BarrierOnFunctionReturn; if (shmemmype() != ROOT_PE) return;

#ifndef MAX_NPES
#  if USE_SHMEM
#    define MAX_NPES 256
#  else
#    define MAX_NPES 1
#  endif
#endif

#define SHMEMINTF_H

#endif
