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

#ifndef METAGRID_H

#ifdef __GNUC__
#  pragma interface
#endif

#include <fstream>
#include "grid.H"
#include "cartgrid.H"
#include "HCgrid.H"
using namespace std;

class Tmetagrid {
private:
	bool dirty;
	smallnat dim;
	bool hc;
	TCartesianGrid<1> *c1;
	THCgrid<1> *h1;
#	if MAXDIM >= 2
	TCartesianGrid<2> *c2;
	THCgrid<2> *h2;
#	endif
#	if MAXDIM >= 3
	TCartesianGrid<3> *c3;
	THCgrid<3> *h3;
#	endif
	void streamload(istream& i);
	void load(const char *fn) {ifstream i(fn); if (i.good()) streamload(i);}
	void dealloc();
public:
	Tmetagrid() {dirty=true;}
	Tmetagrid(const char *fn) {dirty=true; load(fn); /*if (dirty) cerr << "*** Tmetagrid constructor failed.\n";*/}
	bool regular(smallnat dim1, smallnat ncd1, smallnat nsd1, real dx, const real xmin[3], const real xmax[3], bool hcflag=false);
	Tmetagrid(smallnat dim1, smallnat ncd1, smallnat nsd1, real dx, const real xmin[3], const real xmax[3], bool hcflag=false) {
		dirty=true;
		regular(dim1,ncd1,nsd1,dx,xmin,xmax,hcflag);
	}
	bool good() const {return !dirty;}
	smallnat dimension() const {return dim;}
	void getbox(real xmin[3], real xmax[3]) const;
	void get_exterior_box(real xmin[3], real xmax[3]) const;
	void scale(real scaling);
	void centroid(TGridIndex i, Tdimvec& X) const;
	bool iscartesian() const {return !hc;}
	bool ishc() const {return hc;}
	TGridIndex first() const;
	bool isover(TGridIndex) const;
	TGridIndex next(TGridIndex) const;
	bool intpol(const Tdimvec& X, int order=0, bool quiet=false);
	void vintpol(const real X[3][VECLEN], smallnat vlen, int order=0);
	TGridIndex find(const Tdimvec& X) const;
	int Ncelldata() const;
	real MinimumGridSpacing() const;
	real BasegridCellSize() const;
	bool isactive(TGridIndex i) const;
	bool isleaf(TGridIndex i) const;
	TCellType celltype(TGridIndex i) const;
	real cellsize(TGridIndex i) const;
	real& CT(smallnat c, smallnat v=0);
	TRealFormat realformat(TRealFormat rf=REALFORMAT_INQUIRE);
	void set_remove_gaps(bool flag=true);
	void set_intpol_cacheing(bool flag=true);
	void clear_intpol_cache();
	void print_corner_cache_stats(ostream& o);
	smallnat Nneighbours(TGridIndex i, smallnat d, smallnat dir) const;
	bool save(const char *fn, int precision=10, ios::openmode om=ios::out) const;
	// Syntactic sugar for save/load
	friend ostream& operator<<(ostream& o, const Tmetagrid& G);
	friend istream& operator>>(istream& i, Tmetagrid& G) {G.streamload(i); return i;}
	// Sugar for putcell/getcell
	Tmetagrid& operator<<(TGridIndex i);
	Tmetagrid& operator>(TGridIndex i);
	Tmetagrid& operator>>(TGridIndex i);
	~Tmetagrid() {dealloc();}
};

#define METAGRID_H

#endif
