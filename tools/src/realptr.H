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

#ifndef REALPTR_H

#ifdef __GNUC__
#  pragma interface
#endif

/*
  =======================================================================
  The class TRealPtr is defined.
  The macro "real" is defined also (real=double).
  If REALARRAY_DEBUG or DEBUG is defined, length checks are performed.

  The array memory space must be allocated by the user either directly
  or using operator new, for example:

# define N 100
  real araw[3*2*4*N];
  TRealPtr a(araw,3*2*4*N);

  TRealPtr stores the pointer, and the length in debug mode.
  =======================================================================

*/

#include <iostream>
#include <sstream>
#include <cstring>
using namespace std;

typedef double real;
typedef double hareal;		// high-accuracy real used in Tanaka's Riemann-solver
//typedef long double hareal;		// high-accuracy real used in Tanaka's Riemann-solver
typedef int TGridIndex;

#define should_be_const const

#ifdef DEBUG
#  define REALARRAY_DEBUG
#endif

//#define ASSERTIONS

#ifdef DEBUG
#  define ASSERTIONS
#endif

#ifdef ASSERTIONS
   extern void assert_action(const char *file, int lineno);
#  define assert(x) if (!(x)) assert_action(__FILE__,__LINE__)
#else
#  undef assert
#  define assert(x)
#endif

#if defined(__GNUC__) || defined(__DECCXX) || defined(_UNICOS)
#  define HAS_TEMPLATE_INIT_SYNTAX 1
#else
#  define HAS_TEMPLATE_INIT_SYNTAX 0
#endif

#ifdef _UNICOS
#  define pragma(x) _Pragma(x)
#else
#  define pragma(x)
#endif

#define VECLEN 8

#ifdef _UNICOS
typedef int smallnat;
#else
typedef unsigned short smallnat;	// smallnat = small natural number, can be unsigned
#endif

#if VECLEN <= 128
#  if VECLEN <= 64
#    define VDIRS pragma("_CRI ivdep"); pragma("_CRI shortloop"); pragma("_CRI prefervector")
#  else
#    define VDIRS pragma("_CRI ivdep"); pragma("_CRI shortloop128"); pragma("_CRI prefervector")
#  endif
#else
#  define VDIRS pragma("_CRI ivdep"); pragma("_CRI prefervector")
#endif
#  define unroll2
#  define unroll3
#  define unroll4
#define VLOOP(iv) VDIRS; for (v=0; v<(iv).length(); v++)
#define VLOOPN(iv,n) VDIRS; for (v=0; v<(n); v++)

extern real do_some_work(real x);

class TGridIndexVector {
private:
	TGridIndex tab[VECLEN];
	smallnat len;
	void invalidate() {
#ifdef DEBUG
		memset(tab,0xFF,VECLEN*sizeof(TGridIndex));
#endif
	}
public:
	TGridIndexVector() {len=0; invalidate();}
	void setlength(smallnat n) {assert(n<=VECLEN); len=n;}
	TGridIndexVector(int n) {setlength(n); invalidate();}
	smallnat length() const {return len;}
	bool isempty() const {return len==0;}
	TGridIndex& operator[](int i) {assert(0<=i && i<VECLEN); return tab[i];}
	TGridIndex operator()(int i) const {assert(0<=i && i<len); return tab[i];}
};

#ifdef REALARRAY_DEBUG

class TRealPtr {
private:
	real *ptr;
	int len;
//	void assert_action() const;
//	void assert(bool flag) const {if (!flag) assert_action();}
public:
	TRealPtr() {ptr=0; len = -1;}
	TRealPtr(real *p) {ptr=p; len = -1;}
	TRealPtr(real *p, int n) {ptr=p; len=n;}
	void init(int n);
	TRealPtr(int n) {ptr = new real [n]; len=n;}
	real& operator[](int i) const {assert(0<=i && (len < 0 || i<len)); return ptr[i];}
	TRealPtr& operator=(const TRealPtr& rp) {ptr=rp.ptr; len=rp.len; return *this;}
	bool operator!() const {return (ptr == 0);}
	void invalidate() {if (len>0) memset(ptr,0xFF,len*sizeof(real));}
	int length() const {return len;}
	~TRealPtr() {if (ptr && len > 0) delete [] ptr;}
};

#define TCacheBasedRealPtr TRealPtr

#else

class TRealPtr {
private:
	real *ptr;
public:
	TRealPtr() {ptr=0;}
	TRealPtr(real *p) {ptr=p;}
	TRealPtr(real *p, int) {ptr=p;}
	void init(int n);
	real& operator[](int i) const {return ptr[i];}
	TRealPtr& operator=(const TRealPtr& rp) {ptr=rp.ptr; return *this;}
	bool operator!() const {return (ptr == 0);}
	int length() const {return -1;}
	void invalidate() {}
	~TRealPtr() {if (ptr) delete [] ptr;}
};

class TCacheBasedRealPtr {
private:
	real *ptr;
public:
	TCacheBasedRealPtr() {ptr=0;}
	TCacheBasedRealPtr(real *p) {ptr=p;}
	TCacheBasedRealPtr(real *p, int) {ptr=p;}
	void init(int n);
	real& operator[](int i) const {return ptr[i];}
	TCacheBasedRealPtr& operator=(const TCacheBasedRealPtr& rp) {ptr=rp.ptr; return *this;}
	bool operator!() const {return (ptr == 0);}
	int length() const {return -1;}
	void invalidate() {}
	~TCacheBasedRealPtr() {}
};

#endif

#ifndef RECORD_FLOPS
#  define RECORD_FLOPS 1
#endif

extern double flops;            // The global flops counter. Incremented by one by every FLOP, if RECORD_FLOPS is 1

/* We assume that one division is one flop.
   On Cray C90, one sqrt() is 10 additions, 8 multiplications and 3 reciprocals.
   One inverse (1/x) on Cray is 1 reciprocal and 2 mults ==> 6 of the 8 mults
   are dummy. That is, the true numbers for the algorithm used by Cray C90 are:
   10 additions, 2 multiplications, 3 reciprocals, that is 15 FLOPS altogether.
   For Cray T3E, one sqrt() is 49 operations (reported by pat), 47 if vectorized,
   and 35 if the -lmfastv library is also used.
   One reciprocal is only one operation.
*/
#define flops_div 1
#define flops_sqrt 49
#define flops_log flops_sqrt

#if RECORD_FLOPS
#  define recflops(f) flops+=(f)
#  define clear_flops() flops=0
#else
#  define recflops(f)
#  define clear_flops()
#endif

#define REALPTR_H

#endif
