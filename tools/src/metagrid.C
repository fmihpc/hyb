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
#  pragma implementation "metagrid.H"
#endif

#include <fstream>
#include <cstring>
#include "metagrid.H"
#include "fileheader.H"
using namespace std;

void Tmetagrid::streamload(istream& i)
{
	dirty = true;
	const bool use_mapping = false;
	Theader h;
	i >> h;
	if (!h.good()) return;
	dim = smallnat(h.getint("dim"));
	const smallnat ncd = smallnat(h.getint("ncd"));
	const smallnat nsd = smallnat(h.getint("nsd"));
	if (dim == 0 || ncd == 0) return;	// insensible data

//	const TGridIndex maxnc = TGridIndex(h.getint("maxnc"));
	// If merge_load is false, we are probably calling it from "A".
	// In this case it is desirable not to leave any freepool between heap1 and heap2
	// to conserve memory, and it would be unnecessary because we are not going to
	// add any grid cells.
	const TGridIndex maxnc = TGridIndex(h.getint("ncells")) + TGridIndex(h.getint("nablocks"));

	const int n1 = h.getint("n1");
	const int n2 = h.getint("n2");
	const int n3 = h.getint("n3");
	const int ntot = n1*n2*n3;
	real xmin[3], xmax[3], dx;
	xmin[0] = h.getreal("xmin0");
	xmin[1] = h.getreal("xmin1");
	xmin[2] = h.getreal("xmin2");
	xmax[0] = h.getreal("xmax0");
	xmax[1] = h.getreal("xmax1");
	xmax[2] = h.getreal("xmax2");
	dx = h.getreal("dx");
	c1 = 0;
	h1 = 0;
#	if MAXDIM >= 2
	c2 = 0;
	h2 = 0;
#	endif
#	if MAXDIM >= 3
	c3 = 0;
	h3 = 0;
#	endif
	if (!strcmp(h.getstr("type"),"hc"))
		hc = true;
	else if (!strcmp(h.getstr("type"),"cartesian"))
		hc = false;
	else {
		cerr << "*** Type is " << h.getstr("type") << ", not hc nor cartesian,\n";
		return;
	}
	switch (dim) {
	case 1:
		if (hc) {
			h1 = new THCgrid<1>(ncd,nsd,max(ncd,nsd),maxnc,use_mapping);
			h1->setbox(xmin[0],xmax[0]);
			h1->set_fast_regular();
			h1->regular(n1);
			if (!h1->streamload(i,&h)) {
				cerr << "*** (dim=1), hc load failed\n";
				h1 = 0;
				return;
			}
		} else {
			c1 = new TCartesianGrid<1>(ncd,nsd,max(ncd,nsd),ntot,use_mapping);
			c1->setbox(xmin[0],xmax[0]);
			c1->set_fast_regular();
			c1->regular(n1);
			if (!c1->streamload(i,&h)) {
				cerr << "*** (dim=1), cartesian load failed\n";
				c1 = 0;
				return;
			}
		}
		break;
#if MAXDIM >= 2
	case 2:
		if (hc) {
			h2 = new THCgrid<2>(ncd,nsd,max(ncd,nsd),maxnc,use_mapping);
			h2->setbox(xmin[0],xmax[0], xmin[1],xmax[1]);
			h2->set_fast_regular();
			h2->regular(n1,n2);
			if (!h2->streamload(i,&h)) {
				cerr << "*** (dim=2), hc load failed\n";
				h2 = 0;
				return;
			}
		} else {
			c2 = new TCartesianGrid<2>(ncd,nsd,max(ncd,nsd),ntot,use_mapping);
			c2->setbox(xmin[0],xmax[0], xmin[1],xmax[1]);
			c2->set_fast_regular();
			c2->regular(n1,n2);
			if (!c2->streamload(i,&h)) {
				cerr << "*** (dim=2), cartesian load failed\n";
				c2 = 0;
				return;
			}
		}
		break;
#endif
#if MAXDIM >= 3
		case 3:
			if (hc) {
				h3 = new THCgrid<3>(ncd,nsd,max(ncd,nsd),maxnc,use_mapping);
				h3->setbox(xmin[0],xmax[0], xmin[1],xmax[1], xmin[2],xmax[2]);
				h3->set_fast_regular();
				h3->regular(n1,n2,n3);
				if (!h3->streamload(i,&h)) {
					cerr << "*** (dim=3), hc load failed\n";
					h3 = 0;
					return;
				}
			} else {
				c3 = new TCartesianGrid<3>(ncd,nsd,max(ncd,nsd),ntot,use_mapping);
				c3->setbox(xmin[0],xmax[0], xmin[1],xmax[1], xmin[2],xmax[2]);
				c3->set_fast_regular();
				c3->regular(n1,n2,n3);
				if (!c3->streamload(i,&h)) {
					cerr << "*** (dim=3), cartesian load failed\n";
					c3 = 0;
					return;
				}
			}
			break;
#endif
		}
	dirty = false;
}

bool Tmetagrid::regular(smallnat dim1, smallnat ncd1, smallnat nsd1, real dx, const real xmin[3], const real xmax[3], bool hcflag)
{
	dim = dim1;
	dirty = false;
	if (hcflag) {
		cerr << "*** Tmetagrid::regular: hcflag=true not yet implemented, ignored\n";
	}
	hc = false;
	switch (dim) {
	case 1:
		c1 = new TCartesianGrid<1>(ncd1,nsd1,max(ncd1,nsd1));
		c1->setbox(xmin[0],xmax[0]);
		c1->regular(dx);
		break;
#if MAXDIM >= 2
	case 2:
		c2 = new TCartesianGrid<2>(ncd1,nsd1,max(ncd1,nsd1));
		c2->setbox(xmin[0],xmax[0], xmin[1],xmax[1]);
		c2->regular(dx);
		break;
#endif
#if MAXDIM >= 3
	case 3:
		c3 = new TCartesianGrid<3>(ncd1,nsd1,max(ncd1,nsd1));
		c3->setbox(xmin[0],xmax[0], xmin[1],xmax[1], xmin[2],xmax[2]);
		c3->regular(dx);
		break;
#endif
	default:
		cerr << "*** Tmetagrid::regular: bad dimension=" << dim1 << "\n";
		dirty = true;
	}
	return !dirty;
}

#if MAXDIM >= 2
#  define IF2(x) x
#else
#  define IF2(x)
#endif
#if MAXDIM >= 3
#  define IF3(x) x
#else
#  define IF3(x)
#endif

#define CASES(fname,code)\
if (dirty) {warn("Tmetagrid::" fname " called on dirty object"); return;}\
switch (dim) {\
case 1: if (hc) h1->code; else c1->code; break;\
IF2(case 2: if (hc) h2->code; else c2->code; break;)\
IF3(case 3: if (hc) h3->code; else c3->code; break;)\
}

#define RESULTCASES(fname,code)\
if (dirty) {warn("Tmetagrid::" fname " called on dirty object"); return result;}\
switch (dim) {\
case 1: result = hc ? h1->code : c1->code; break;\
IF2(case 2: result = hc ? h2->code : c2->code; break;)\
IF3(case 3: result = hc ? h3->code : c3->code; break;)\
}

void Tmetagrid::getbox(real xmin[3], real xmax[3]) const {CASES("getbox",getbox(xmin,xmax));}

void Tmetagrid::get_exterior_box(real xmin[3], real xmax[3]) const {
	CASES("get_exterior_box",get_exterior_box(xmin,xmax));
}

void Tmetagrid::scale(real scaling)
{
	CASES("scale",scale(scaling));
}

void Tmetagrid::centroid(TGridIndex i, Tdimvec& X) const {CASES("centroid",centroid(i,X));}

bool Tmetagrid::intpol(const Tdimvec& X, int order, bool quiet) {
	bool result=false;
	RESULTCASES("intpol",intpol(X,order,quiet));
	return result;
}

void Tmetagrid::vintpol(const real X[3][VECLEN], smallnat vlen, int order) {CASES("vintpol",vintpol(X,vlen,order));}

TGridIndex Tmetagrid::find(const Tdimvec& X) const {
	TGridIndex result = NOINDEX;
	RESULTCASES("find",find(X));
	return result;
}

TGridIndex Tmetagrid::first() const
{
	TGridIndex result = NOINDEX;
	RESULTCASES("first",first());
	return result;
}

bool Tmetagrid::isover(TGridIndex i) const
{
	bool result = true;
	RESULTCASES("isover",isover(i));
	return result;
}

TGridIndex Tmetagrid::next(TGridIndex i) const
{
	TGridIndex result = i;
	RESULTCASES("next",next(i));
	return result;
}

int Tmetagrid::Ncelldata() const
{
	int result = 0;
	RESULTCASES("Ncelldata",Ncelldata());
	return result;
}

smallnat Tmetagrid::Nneighbours(TGridIndex i, smallnat d, smallnat dir) const
{
	smallnat result = 0;
	RESULTCASES("Nneighbours",Nneighbours(i,d,dir));
	return result;
}

bool Tmetagrid::save(const char *fn, int precision, ios::openmode om) const
{
	bool result = false;
	RESULTCASES("save",save(fn,precision,om));
	return result;
}

TRealFormat Tmetagrid::realformat(TRealFormat rf)
{
	TRealFormat result = REALFORMAT_INQUIRE;
	RESULTCASES("realformat",realformat(rf));
	return result;
}

void Tmetagrid::set_remove_gaps(bool flag)
{
	CASES("set_remove_gaps",set_remove_gaps(flag));
}

void Tmetagrid::set_intpol_cacheing(bool flag)
{
	CASES("set_intpol_cacheing",set_intpol_cacheing(flag));
}

void Tmetagrid::clear_intpol_cache()
{
	CASES("clear_intpol_cache",clear_intpol_cache());
}

void Tmetagrid::print_corner_cache_stats(ostream& o)
{
	CASES("print_corner_cache_stats",print_corner_cache_stats(o));
}

real Tmetagrid::MinimumGridSpacing() const
{
	real result = -1;
	RESULTCASES("MinimumGridSpacing",MinimumGridSpacing());
	return result;
}

real Tmetagrid::BasegridCellSize() const
{
	real result = -1;
	RESULTCASES("BasegridCellSize",BasegridCellSize());
	return result;
}

bool Tmetagrid::isactive(TGridIndex i) const
{
	bool result = false;
	RESULTCASES("isactive",isactive(i));
	return result;
}

bool Tmetagrid::isleaf(TGridIndex i) const
{
	bool result = false;
	RESULTCASES("isleaf",isleaf(i));
	return result;
}

TCellType Tmetagrid::celltype(TGridIndex i) const
{
	TCellType result = INTERIOR_CELL;
	RESULTCASES("celltype",celltype(i));
	return result;
}

real Tmetagrid::cellsize(TGridIndex i) const
{
	real result = 0.0;
	RESULTCASES("cellsize",cellsize(i));
	return result;
}

static real dummy_real = 0;

real& Tmetagrid::CT(smallnat c, smallnat v)
{
	if (dirty) {warn("Tmetagrid::CT called on dirty object"); return dummy_real;}
	switch (dim) {
	case 1:
		return hc ? h1->CT(c,v) : c1->CT(c,v);
#if MAXDIM >= 2
	case 2:
		return hc ? h2->CT(c,v) : c2->CT(c,v);
#endif
#if MAXDIM >= 3
	case 3:
		return hc ? h3->CT(c,v) : c3->CT(c,v);
#endif
	}
	return dummy_real;
}

ostream& operator<<(ostream& o, const Tmetagrid& G)
{
	switch (G.dim) {
	case 1:
		if (G.hc) o << *G.h1; else o << *G.c1;
		break;
#if MAXDIM >=2
	case 2:
		if (G.hc) o << *G.h2; else o << *G.c2;
		break;
#endif
#if MAXDIM >= 3
	case 3:
		if (G.hc) o << *G.h3; else o << *G.c3;
		break;
#endif
	}
	return o;
}

Tmetagrid& Tmetagrid::operator<<(TGridIndex i)
{
	switch (dim) {
	case 1: if (hc) h1->putcell(i); else c1->putcell(i); break;
#if MAXDIM >= 2
	case 2: if (hc) h2->putcell(i); else c2->putcell(i); break;
#endif
#if MAXDIM >= 3
	case 3: if (hc) h3->putcell(i); else c3->putcell(i); break;
#endif
	}
	return *this;
}

Tmetagrid& Tmetagrid::operator>(TGridIndex i)
{
	switch (dim) {
	case 1: if (hc) {h1->clearcache(); h1->fgetcell(i);} else c1->fgetcell(i); break;
#if MAXDIM >= 2
	case 2: if (hc) {h2->clearcache(); h2->fgetcell(i);} else c2->fgetcell(i); break;
#endif
#if MAXDIM >= 3
	case 3: if (hc) {h3->clearcache(); h3->fgetcell(i);} else c3->fgetcell(i); break;
#endif
	}
	return *this;
}

Tmetagrid& Tmetagrid::operator>>(TGridIndex i)
{
	switch (dim) {
	case 1: if (hc) {h1->clearcache(); h1->getcell(i,true);} else c1->getcell(i,true); break;
#if MAXDIM >= 2
	case 2: if (hc) {h2->clearcache(); h2->getcell(i,true);} else c2->getcell(i,true); break;
#endif
#if MAXDIM >= 3
	case 3: if (hc) {h3->clearcache(); h3->getcell(i,true);} else c3->getcell(i,true); break;
#endif
	}
	return *this;
}

void Tmetagrid::dealloc()
{
	if (dirty) return;	// without doing anything
	switch (dim) {
	case 1: if (hc) delete h1; else delete c1; break;
#if MAXDIM >= 2
	case 2: if (hc) delete h2; else delete c2; break;
#endif
#if MAXDIM >= 3
	case 3: if (hc) delete h3; else delete c3; break;
#endif
	}
}

