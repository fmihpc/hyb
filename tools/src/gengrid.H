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

#ifndef GENGRID_H

#ifdef __GNUC__
#  pragma interface
#endif

#if defined(__GNUC__) && __GNUC__ < 3
#  define min(x,y) ((x) <? (y))
#  define max(x,y) ((x) >? (y))
//#  define sqr(x) ({__typeof(x) _arg = (x); _arg*_arg;})
#  define sqr(x) ((x)*(x))
#else
#  include <algorithm>
   using std::max;
   using std::min;
//   template <class T> inline T min(T x, T y) {return (x < y) ? x : y;}
//   template <class T> inline T max(T x, T y) {return (x > y) ? x : y;}
   template <class T> inline T sqr(T x) {return x*x;}
#endif

#include "realptr.H"
#include "cache.H"

#define warn(s) \
    { \
		 static int count = 0; \
		 static int cnt_interval = 1; \
	     if (count++ % cnt_interval == 0) { \
			 cerr << "*** (" << cnt_interval << ") Warning at line " << __LINE__ << " in file " << __FILE__ \
				  << ": " << (s) << "\n"; \
			 cnt_interval*= 10; \
		 } \
	}

#define MAXDIM 3
#define MAXNEI (1 << (MAXDIM-1))		// maximum number of occuring maxnei

#define disp(x) cout << #x " = " << x << "\n"

class machine {
private:
	static real eps;				// machine epsilon
	static real tny;				// smallest possible real number
public:
	static void init();
//	machine() {init();}
	static real epsilon() {return eps;}
	static real tiny() {return tny;}
};

// Class TGenericGrid is the dimension-independent base class,
// whereas Tgrid, TCartesianGrid, THCgrid are templates over dim.
// TGenericGrid holds the buffer tables and so can be used for interfacing
// with the solver routines. Thus the solvers need not be templates over
// dim, which saves us some code triplication.

class TGenericGrid
{
private:
	TCacheBasedRealPtr fluxtab,areatab,utab,celltab,surftab,nntab;
protected:
	smallnat ncd;					// Number of Cell Data items
	smallnat ncd_saved;				// Number of saved Cell Data items, default=ncd. Affects streamsave() only.
	smallnat nsd;					// Number of Surf Data items
	smallnat nfq;					// Number of Flux Quantities
	bool dirty;
	void invalidateFT() {fluxtab.invalidate();}
	void invalidateAT() {areatab.invalidate();}
	void invalidateUT() {utab.invalidate();}
	void invalidateCT() {celltab.invalidate();}
	void invalidateST() {surftab.invalidate();}
	void invalidateNN() {nntab.invalidate();}
	void invalidate() {
		invalidateFT(); invalidateUT(); invalidateCT(); invalidateST(); invalidateAT(); invalidateNN();
	}
public:
	TGenericGrid() {dirty=true; machine::init();}
	void init(smallnat ncd1, smallnat nsd1, smallnat nfq1);
	void set_saved_ncd(smallnat ncd_saved1) {ncd_saved = ncd_saved1;}
	TGenericGrid(smallnat ncd1, smallnat nsd1, smallnat nfq1) {init(ncd1,nsd1,nfq1);}
	// Access functions to buffer tables
	real& FT(smallnat d, smallnat dir, smallnat s, smallnat v=0) const;
	real& AT(smallnat d, smallnat dir, smallnat v=0) const;
	real& UT(smallnat dir, smallnat c, smallnat v=0) const;
	real& CT(smallnat c, smallnat v=0) const;
	real& ST(smallnat k, smallnat s, smallnat v=0) const;
	real& NN(smallnat v=0) const;
	~TGenericGrid();
};

inline real& TGenericGrid::FT
    (smallnat d, smallnat dir, smallnat s, smallnat v) const
{
	assert(d < MAXDIM);
	assert(dir==0 || dir==1);
	assert(s < nsd);
	assert(v < VECLEN);
#	ifdef VECTOR_MACHINE
	// [MAXDIM][2][nsd][VECLEN]
	return fluxtab[((d*2+dir)*nsd+s)*VECLEN+v];
#	else
	// [VECLEN][MAXDIM][2][nsd]
	return fluxtab[((v*MAXDIM+d)*2+dir)*nsd+s];
#	endif
}

inline real& TGenericGrid::AT
    (smallnat d, smallnat dir, smallnat v) const
{
	assert(d < MAXDIM);
	assert(dir==0 || dir==1);
	assert(v < VECLEN);
#	ifdef VECTOR_MACHINE
	// [MAXDIM][2][VECLEN]
	return areatab[(d*2+dir)*VECLEN+v];
#	else
	// [VECLEN][MAXDIM][2]
	return areatab[(v*MAXDIM+d)*2+dir];
#	endif
}

inline real& TGenericGrid::UT(smallnat dir, smallnat c, smallnat v) const
{
	assert(dir==0 || dir==1);
	assert(c < ncd);
	assert(v < VECLEN);
#	ifdef VECTOR_MACHINE
	// [2][ncd][VECLEN]
	return utab[(dir*ncd+c)*VECLEN+v];
#	else
	// [VECLEN][2][ncd]
	return utab[(v*2+dir)*ncd+c];
#	endif
}

inline real& TGenericGrid::CT(smallnat c, smallnat v) const
{
	assert(c < ncd);
	assert(v < VECLEN);
#	ifdef VECTOR_MACHINE
	// [ncd][VECLEN]
	return celltab[c*VECLEN+v];
#	else
	// [VECLEN][ncd]
	return celltab[v*ncd+c];
#	endif
}

inline real& TGenericGrid::ST(smallnat k, smallnat s, smallnat v) const
{
	assert(s < nsd);
	assert(v < VECLEN);
#	ifdef VECTOR_MACHINE
	// [MAXNEI][nsd][VECLEN]
	return surftab[(k*nsd+s)*VECLEN+v];
#	else
	// [MAXNEI][VECLEN][nsd]
	return surftab[(k*VECLEN+v)*nsd+s];
#	endif
}

inline real& TGenericGrid::NN(smallnat v) const
{
	return nntab[v];
}

#define GENGRID_H

#endif
