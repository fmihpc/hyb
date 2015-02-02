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
#  pragma implementation "grid.H"
#endif

#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstdio>
#include "grid.H"
#include "state.H"
using namespace std;

/*
  Interface to grid data
  ----------------------

  Buffer tables
  -------------

  CV(c,v)				Cell table		Filled by getcell/vgetcell
  FT(d,dir,s,v)         Flux table      Filled by getcell/vgetcell, if FTflag==true
    FT values are the flux intensities. They must be multiplied by areas (AT) to get the total fluxes.
  AT(d,dir,v)           Surf area table Filled by getcell/vgetcell, if FTflag==true
  NT(d,dir,v)           Number of neighbours table  Filled by getcell/vgetcell, if FTflag==true
  ST(k,s,v)				Surf table		Needed by putsurf/vputsurf. Filled by (v)getsurf, if STflag==true
  UT(dir,c,v)			U-table			Filled by getsurf/vgetsurf
  AT2(v)                2nd area table  Filled by getsurf/vgetsurf
  NN(v)                 Neighbour index table Filled by getsurf/vgetsurf

  c: in 0..ncd-1, v: in 0..VECLEN-1, s: in 0..nsd-1,
  d: in 0..dim-1, dir: either 0 or 1, nei: 0..maxnei-1

  c labels cell data items
  v is vector index (last argument with default 0 in most cases)
  s labels surf data items
  d,dir,nei label dimension, direction and "which surface".
  i is grid index (labels cells in the grid)

  Functions
  ---------

  j = getneighbour(i, d,dir)
      Return neighbouring cell in dimension d (0,1,2), direction dir (0=left,1=right).
	  Return NOINDEX if exceeds domain.
	  Notice that in case of HC grid, j may not be leaf cell even when i is!
  
  getcell(i, FTflag=true, v=0)
      Transfer cell i data to CT(c,v) for all c.
	  If FTflag is true, also fill FT(d,dir,s,v) and AT(d,dir,v)
	    for all d,dir,nei,s;
	    that is, transfer the surf data of all touching surfs in flux table FT,
		the corresponding surf areas (area table AT),
		and the number of neighbours for surf (d,dir).
      Getcell is used to retrieve all data related to a single cell.

  putcell(i,v=0)
      The reverse of getcell. Transfer data from CV(c,v) (for all c)
	  into the grid. FT is not used.
	  Putcell is used to set the cell data for a single cell.

  getsurf(i,ss,STflag=false, v=0)   (TSurfSpec ss)
      Transfer surf (i,ss) left and right cell data
	  to UT(dir,c,v) (dir=0,1). A surface is specified by the pair (i,ss)
	  where i defines a cell and ss contains (d,dir,nei).
	  Getsurf is used to retrieve surf's touching cell pair data.
	  Also sets the surface area in AT2() and the neighbour index to NN().
	  If STflag==true, also ST(0,s,v) is filled.

  putsurf(i,ss, v=0)
      The reverse of getsurf, except that UT is not used,
	  and that putsurf uses ST(k,s,v) whereas getsurf produces ST(0,s,v).
	  Putsurf is used to set the surf data for a single surf.

  Vector versions
  ---------------
  Trivial vector versions are provided in this file. They call the scalar
  versions in a loop. They should be overwritten for better performance in
  derived classes, especially for vector machines.

  vgetneighbour(jv, iv,d,dir)
      All jv,iv,d,dir are vectors. d,dir are smallnat vectors, jv,iv TGridIndexVectors.
  vgetneighbour_s(jv, iv,d,dir)
      The same as vgetneighbour, but d,dir are scalars, not vectors.
  vgetcell(iv, FTflag=true)
      Same as getcell but for many cells simultaneously. Fills CT(c,v) for all v,
	  and FT, AT and NT if FTflag==true.
  vputcell(iv)
      Set cell data for many cells simultaneously. Needs CT(c,v) for all v.
  vgetsurc(iv,sv)
      Transfer to ST(k,s,v), UT(dir,c,v), AT2(v), NN(v).
  vputsurf(iv,sv)
      Transfer from ST(k,s,v) into grid.
*/


#ifndef PROMOTE_INLINING
template <smallnat dim>
void Tgrid<dim>::vnext1(TGridIndexVector& iv, TGridIndex start, TGridIndex stride, Ttraversal trav) const
{
	TGridIndex i;
	smallnat j;
	for (i=start,j=0; !isover(i) && j<VECLEN; i+=stride)
		if (trav==TRAV_ALL || isactive_vfn(i)) iv[j++] = i;
	iv.setlength(j);
	warn("using Tgrid::vnext1");
}
#endif

template <smallnat dim>
void Tgrid<dim>::getbox(real Xmin[3], real Xmax[3]) const
{
	smallnat d;
	for (d=0; d<dim; d++) {
		Xmin[d] = xmin[d];
		Xmax[d] = xmax[d];
	}
}

template <smallnat dim>
void Tgrid<dim>::get_exterior_box(real Xmin[3], real Xmax[3]) const
{
	getbox(Xmin,Xmax);
	smallnat d;
	for (d=0; d<dim; d++) {
		Xmin[d]-= NB*dx;
		Xmax[d]+= NB*dx;
	}
}

template <smallnat dim>
void Tgrid<dim>::writeinfo(ostream& o, TGridIndex i) const
{
	o << "[i=" << i << " ";
	switch (celltype(i)) {
	case INTERIOR_CELL: o << "interior"; break;
	case GHOST_CELL: o << "ghost"; break;
	case DEAD_CELL: o << "dead"; break;
	case REMOVED_CELL: o << "removed"; break;
	}
	Tdimvec X;
	centroid(i,X);
	o << " X=" << X(0);
	if (dim >= 2) o << " Y=" << X(1);
	if (dim == 3) o << " Z=" << X(2);
	o << ", lev=" << level_vfn(i) << "]";
}

template <smallnat dim>
void Tgrid<dim>::default_map(should_be_const real Xin[MAXDIM][VECLEN], real Xout[MAXDIM][VECLEN], smallnat vlen)
{
	smallnat v,d;
	for (d=0; d<dim; d++) {
		VDIRS;
		for (v=0; v<vlen; v++) Xout[d][v] = Xin[d][v];
	}
}

template <smallnat dim>
void Tgrid<dim>::init(smallnat ncd1, smallnat nsd1, smallnat nfq1, TGridIndex maxnc1, bool use_mapping1)
{
	TGenericGrid::init(ncd1,nsd1,nfq1);
	use_mapping = use_mapping1;
	merge_load = false;
	remove_gaps_when_saving = true;
	parallel_io = true;		// unless the machine has SHMEM and Npes>1, this has no effect
	nsd_priv = use_mapping1 ? 1+MAXDIM : 1;			// Number of surf data items, private
	// 1+MAXDIM: surf area (1 item), plus normal vector (MAXDIM items)
	// 1: just the surf area (1 item)
	ncd_priv = use_mapping1 ? 2 : 1;				// Number of cell data items, private
	maxnei = 1;
	output_realformat = REALFORMAT_DOUBLE;
	Nbcs = 0;
	setbox_called = false;
	regular_called = false;
	do_fast_regular = false;
	cache_intpols = false;
	cell_prepare_hook = 0;
	source_term_hook = 0;
	
	maxnc = DivNpes(maxnc1+Npes-1)*Npes;
	if (ncd1 > max_ncd) {cerr << "*** Tgrid::Tgrid: ncd=" << ncd1 << " is too large, maximum is " << max_ncd << "\n"; exit(1);}
	if (nsd1 > max_ncd) {cerr << "*** Tgrid::Tgrid: nsd=" << nsd1 << " is too large, maximum is " << max_ncd << "\n"; exit(1);}
	clen = ncd1 + dim*(nsd1+nsd_priv) + ncd_priv;
	xmin[0] = xmin[1] = xmin[2] = 0;
	xmax[0] = xmax[1] = xmax[2] = 1;
	smallnat d;
	for (d=0; d<MAXDIM; d++) {
		box_bctypes[d][0] = box_bctypes[d][1] = NEUMANN_BC;
		box_bcfuncs[d][0] = box_bcfuncs[d][1] = 0;
	}
	map = &default_map;
}

template <smallnat dim>
TRealFormat Tgrid<dim>::realformat(TRealFormat rf)
{
	if (rf != REALFORMAT_INQUIRE) output_realformat = rf;
	return output_realformat;
}

template <smallnat dim>
void Tgrid<dim>::setmap(void(*map1)(should_be_const real[MAXDIM][VECLEN], real [MAXDIM][VECLEN], smallnat))
{
	if (!use_mapping) cerr << "*** Tgrid::setmap warning: Calling setup while use_mapping=false in constructor\n";
	map = map1;
}

template <smallnat dim>
void Tgrid<dim>::setbox(real x1, real x2)
{
	if (regular_called) {
		cerr << "*** Tgrid<" << dim << ">::setbox(): regular() already called\n";
	}
	xmin[0] = x1;
	xmax[0] = x2;
	setbox_called=true;
}

template <smallnat dim>
void Tgrid<dim>::setbox(real x1, real x2, real y1, real y2)
{
	if (regular_called) {
		cerr << "*** Tgrid<" << dim << ">::setbox(): regular() already called\n";
	}
	xmin[0] = x1;
	xmax[0] = x2;
	xmin[1] = y1;
	xmax[1] = y2;
	setbox_called=true;
}

template <smallnat dim>
void Tgrid<dim>::setbox(real x1, real x2, real y1, real y2, real z1, real z2)
{
	if (regular_called) {
		cerr << "*** Tgrid<" << dim << ">::setbox(): regular() already called\n";
	}
	xmin[0] = x1;
	xmax[0] = x2;
	xmin[1] = y1;
	xmax[1] = y2;
	xmin[2] = z1;
	xmax[2] = z2;
	setbox_called = true;
}

#ifdef _CRAYT3E
#include <intrinsics.h>
#endif

real our_variable1, our_variable2;

template <smallnat dim>
void Tgrid<dim>::regular(real dx1)
{
	smallnat d;
	dx = dx1; cornercache.set_dx(dx);
	for (d=0; d<dim; d++) {
		n123[d] = int(floor((xmax[d] -xmin[d])/dx + 1e-5));
		xmax[d] = xmin[d] + n123[d]*dx;
		n123[d]+= 2*NB;
	}
	n1 = n123[0]; n2 = n123[1]; n3 = n123[2];
	ntot = n1*n2*n3;
	
	// Add box boundary conditions
	for (d=0; d<dim; d++) {
		const real deltax = 10*machine::epsilon()*max(fabs(xmin[d]),fabs(xmax[d]));
		const real xval_low = xmin[d] - deltax;
		const real xval_high = xmax[d] + deltax;
		addBC(TPlaneBoundaryGeometry(dim,d,xval_low,false), box_bctypes[d][0], true, box_bcfuncs[d][0]);
		addBC(TPlaneBoundaryGeometry(dim,d,xval_high,true), box_bctypes[d][1], true, box_bcfuncs[d][1]);
	}
	regular_called = true;
	if (do_fast_regular) return;
	setup_after_regular();
	shmembarrierall();

	bool were_ghosts;
	if (mype == ROOT_PE) mature_cells(0,ntot-1, were_ghosts);
	shmembarrierall();
}

template <smallnat dim>
void Tgrid<dim>::set_box_BC(smallnat d, smallnat dir, TBCtype type, TEvaluatorFunctionPtr ev)
{
	if (regular_called) cerr << "Tgrid<" << dim << ">::set_box_BC: regular() has already been called, no effect\n";
	box_bctypes[d][dir] = type;
	box_bcfuncs[d][dir] = ev;
}

template <smallnat dim>
void Tgrid<dim>::regular(int ntot1)
{
	if (ntot1 < 1 || ntot1 > maxnc) {
		cerr << "*** Tgrid::regular(" << ntot1 << "): arg out of range 1.." << maxnc << "\n";
		exit(1);
	}
	if (dim == 1) {
		const real dx1 = (xmax[0] - xmin[0])/(ntot1 - 2*NB);
		regular(dx1);
	} else if (dim == 2) {
		const real dx1 = sqrt( (xmax[0]-xmin[0])*(xmax[1]-xmin[1])/ntot1 );
		int i1 = int(floor((xmax[0]-xmin[0])/dx1)) - 2*NB;
		int i2 = int(floor((xmax[1]-xmin[1])/dx1)) - 2*NB;
		regular(i1,i2);
	} else {
		const real dx1 = pow( (xmax[0]-xmin[0])*(xmax[1]-xmin[1])*(xmax[2]-xmin[2])/ntot1, 1.0/3.0 );
		int i1 = int(floor((xmax[0]-xmin[0])/dx1)) - 2*NB;
		int i2 = int(floor((xmax[1]-xmin[1])/dx1)) - 2*NB;
		int i3 = int(floor((xmax[2]-xmin[2])/dx1)) - 2*NB;
		regular(i1,i2,i3);
	}
}

template <smallnat dim>
void Tgrid<dim>::regular(int i1, int i2)
{
	if (i1 < 1 || i2 < 1 || i1*i2 > maxnc) {
		cerr << "*** Tgrid::regular(" << i1 << "," << i2 << "): bad args, max product is " << maxnc << "\n";
		exit(1);
	}
	if (dim != 2) {
		cerr << "*** Tgrid::regular(" << i1 << "," << i2 << "): dim=" << dim << " should be 2\n";
		exit(1);
	}
	const real dx1 = (xmax[0] - xmin[0])/(i1 - 2*NB);
	const real dy1 = (xmax[1] - xmin[1])/(i2 - 2*NB);
	regular(max(dx1,dy1));
}

template <smallnat dim>
void Tgrid<dim>::regular(int i1, int i2, int i3)
{
	if (i1 < 1 || i2 < 1 || i3 < 1 || i1*i2*i3 > maxnc) {
		cerr << "*** Tgrid::regular(" << i1 << "," << i2 << "," << i3 << "):\n"
			 << "    bad args, max product is " << maxnc << "\n";
		exit(1);
	}
	if (dim != 3) {
		cerr << "*** Tgrid::regular(" << i1 << "," << i2 << "," << i3 << "):\n"
			 << "    dim=" << dim << " should be 3\n";
		exit(1);
	}
	const real dx1 = (xmax[0] - xmin[0])/(i1 - 2*NB);
	const real dy1 = (xmax[1] - xmin[1])/(i2 - 2*NB);
	const real dz1 = (xmax[2] - xmin[2])/(i3 - 2*NB);
	regular(max(max(dx1,dy1),dz1));
}

template <smallnat dim>
void Tgrid<dim>::scale(real scaling)
// Multiply dx,xmin,xmax by scaling
{
	dx*= scaling; cornercache.set_dx(dx);
	smallnat d;
	for (d=0; d<dim; d++) {
		xmin[d]*= scaling;
		xmax[d]*= scaling;
	}
}

template <smallnat dim>
Tgrid<dim>::~Tgrid()
{
	if (!dirty) {
//		if (use_mapping) delete [] physsize.RawPtr();
//		delete [] invvol.RawPtr();
//		delete [] (char*)cellinfo;
		//delete [] cellinfo;
	}
}

template <smallnat dim>
void Tgrid<dim>::vgetcell(const TGridIndexVector& iv, bool FTflag)
{
	smallnat v;
	VLOOP(iv) getcell(iv(v),FTflag,v);
}

template <smallnat dim>
void Tgrid<dim>::vgetcell_allow_noindex(const TGridIndexVector& iv)
{
	smallnat v;
	VLOOP(iv) if (iv(v) != NOINDEX) getcell(iv(v),false,v);
}

template <smallnat dim>
void Tgrid<dim>::vgetnormal(const TGridIndexVector& iv, const TSurfSpecVector& sv, Tvec n3[3])
{
	smallnat v;
	real n[3];
	VLOOP(iv) {
		getnormal(iv(v),sv(v), n);
		n3[0][v] = n[0]; n3[1][v] = n[1]; n3[2][v] = n[2];
	}
}

template <smallnat dim>
void Tgrid<dim>::vputcell(const TGridIndexVector& iv)
{
	smallnat v;
	VLOOP(iv) putcell(iv(v),v);
}

template <smallnat dim>
void Tgrid<dim>::vgetsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv, bool STflag)
{
	smallnat v;
	VLOOP(iv) getsurf(iv(v),sv(v),STflag,v);
}

template <smallnat dim>
void Tgrid<dim>::vputsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv)
{
	smallnat v;
	VLOOP(iv) putsurf(iv(v),sv(v),v);
}

template <smallnat dim>
bool Tgrid<dim>::save1(ostream& o) const
{
	o << "#\n";
	o << "dim = " << dim << "\n";
	o << "ncd = " << ncd_saved << "\n";
	o << "nsd = " << 0 /*nsd*/ << "\n";
	o << "maxnc = " << maxnc << "\n";
	o << "realformat = ";
	switch (output_realformat) {
	case REALFORMAT_DOUBLE: o << "double"; break;
	case REALFORMAT_FLOAT: o << "float"; break;
	case REALFORMAT_ASCII: o << "ascii"; break;
	case REALFORMAT_INQUIRE: o << "INQUIRE"; break;
	}
	o << "\n";
	o << "n1 = " << n1 << "\n";
	o << "n2 = " << n2 << "\n";
	o << "n3 = " << n3 << "\n";
	const int oldprec = o.precision();
	o.precision(16);
	o << "xmin0 = " << xmin[0] << "\n";
	o << "xmin1 = " << xmin[1] << "\n";
	o << "xmin2 = " << xmin[2] << "\n";
	o << "xmax0 = " << xmax[0] << "\n";
	o << "xmax1 = " << xmax[1] << "\n";
	o << "xmax2 = " << xmax[2] << "\n";
	o << "dx = " << dx << "\n";
	o.precision(oldprec);
	return o.good();
}

template <smallnat dim>
bool Tgrid<dim>::load(const char *fn)
{
	ifstream f(fn);
	if (!f.good()) return false;
	return streamload(f);
}

template <smallnat dim>
bool Tgrid<dim>::save(const char *origfn, int precision, ios::openmode om) const
{
	bool retval = false;
	if (Npes > 1 && parallel_io) {
		char *fn = new char [strlen(origfn) + 20];
		sprintf(fn,"%s.head",origfn);
		ofstream f(fn,ios::out);		// The header file must be rewritten even if om == ios::app
		f.precision(precision);
		if (f.good()) {
			retval = streamsave(f,origfn);
			f << flush;
		}
		delete [] fn;
	} else {
		ofstream f(origfn,om);
		f.precision(precision);
		if (f.good()) {
			retval = streamsave(f,0 /*origfn*/);
			f << flush;
		}
		f.close();
	}
	return retval;
}

template <smallnat dim>
void Tgrid<dim>::corner_deplist
    (const Tdimvec& X, TGridIndex depcells[8], real weights[8], smallnat& Ndepcells) const
{
	TGridIndex i;
	smallnat d,ch;
#ifndef nchildren
	const smallnat nchildren = (1 << dim);
#endif
	Tdimvec X1;
	int maxn = 0;
	for (d=0; d<dim; d++) maxn = max(maxn, n123[d]);
	const real deltax = 100*machine::epsilon()*dx*maxn;
	static const int dir[8][3] =
	{{-1,-1,-1},
	 {+1,-1,-1},
	 {-1,+1,-1},
	 {+1,+1,-1},
	 {-1,-1,+1},
	 {+1,-1,+1},
	 {-1,+1,+1},
	 {+1,+1,+1}};
	smallnat c = 0;
	for (ch=0; ch<nchildren; ch++) {
		for (d=0; d<dim; d++) {
			X1[d] = X(d) + dir[ch][d]*deltax;
		}
		i = find(X1);
		if (i != NOINDEX && celltype(i) != DEAD_CELL && celltype(i) != REMOVED_CELL) {
			depcells[c] = i;
			weights[c] = 1.0/cellsize_vfn(i);
			c++;
		}
	}
	Ndepcells = c;
	real weightsum = 0.0;
	for (ch=0; ch<c; ch++) weightsum+= weights[ch];
	if (weightsum != 0) {
		const real invweightsum = 1.0/weightsum;
		for (ch=0; ch<c; ch++) weights[ch]*= invweightsum;
	}
}

template <smallnat dim>
bool Tgrid<dim>::intpol(const Tdimvec& X, int order, bool quiet)
{
	TGridIndex i = find(X);
	if (i == NOINDEX) {
		if (!quiet) warn("Tgrid::intpol: point is outside domain");
		return false;
	}
	clearcache();
	if (order == 0) {
		// Zeroth-order interpolation (just return contents of cell i)
		getcell(i,false);
	} else if (order == 1) {
		// Linear interpolation
		smallnat c,d,n,ch;
		TGridIndex depcells[8][8];
		real weights[8][8];
		smallnat Ndepcells[8];
#ifndef nchildren
		const smallnat nchildren = (1 << dim);
#endif
		const real deltax = cellsize_vfn(i);
		const real invdeltax = 1.0/deltax;
		const real halfdeltax = 0.5*deltax;
		Tdimvec Xc,Xlow;
		centroid(i,Xc);
		for (d=0; d<dim; d++) Xlow[d] = Xc(d) - halfdeltax;
		static const int dir[8][3] =
		{{-1,-1,-1},
		 {+1,-1,-1},
		 {-1,+1,-1},
		 {+1,+1,-1},
		 {-1,-1,+1},
		 {+1,-1,+1},
		 {-1,+1,+1},
		 {+1,+1,+1}};
		Tdimvec X1;
		real ucorner[max_ncd], result[max_ncd];
		real t[3];
		for (d=0; d<dim; d++) {
			t[d] = (X(d) - Xlow(d))*invdeltax;
			assert(-1e5*machine::epsilon()<=t[d] && t[d]<=1+1e5*machine::epsilon());
		}
		for (c=0; c<ncd; c++) result[c] = 0.0;
		real cornerweightsum = 0;
		for (ch=0; ch<nchildren; ch++) {
			real cornerweight = 1.0;
			for (d=0; d<dim; d++) {
				cornerweight*= (dir[ch][d] > 0) ? t[d] : 1-t[d];
			}
			if (fabs(cornerweight) > 1e-10) {
				for (d=0; d<dim; d++) {
					X1[d] = Xc(d) + dir[ch][d]*halfdeltax;
				}
				bool found_in_cache = false;
				if (cache_intpols) found_in_cache = cornercache.probe(X1,dim,ucorner,ncd);
				if (!found_in_cache) {
					corner_deplist(X1,depcells[ch],weights[ch],Ndepcells[ch]);
					for (c=0; c<ncd; c++) ucorner[c] = 0.0;
					for (n=0; n<Ndepcells[ch]; n++) {
						const real W = weights[ch][n];
						getcell(depcells[ch][n],false);	// get cornerdep cell contents in CT(c)
						for (c=0; c<ncd; c++) ucorner[c]+= W*CT(c);
					}
					if (cache_intpols) cornercache.store(X1,dim,ucorner,ncd);
				}
				for (c=0; c<ncd; c++) result[c]+= cornerweight*ucorner[c];
				cornerweightsum+= cornerweight;
			}
		}
		assert(fabs(cornerweightsum-1) < 1e-6);
		for (c=0; c<ncd; c++) CT(c) = result[c];
	} else {
		warn("Tgrid::intpol: only order=0 and order=1 currently supported");
		// Use zeroth-order interpolation in this case
		getcell(i,false);
	}
	// result is in CT(c,0)
	return true;	// successful return
}

template <smallnat dim>
void Tgrid<dim>::vintpol(const real X[MAXDIM][VECLEN], smallnat vlen, int order)
{
	smallnat v,d,c;
	Tdimvec X1;
	// Do v=0 later, therefore start from 1 here
	for (v=1; v<vlen; v++) {
		for (d=0; d<dim; d++) X1[d] = X[d][v];
		intpol(X1,order);
		for (c=0; c<ncd; c++) CT(c,v) = CT(c,0);
	}
	// Do v=0 component only now so that correct result remains in CT(c,0)
	for (d=0; d<dim; d++) X1[d] = X[d][0];
	intpol(X1,order);
}

template <smallnat dim>
void Tgrid<dim>::addBC(const TBoundaryGeometry& geom, TBCtype type, bool isbox, TEvaluatorFunctionPtr ev)
{
	if (Nbcs >= MAX_N_BC) {
		cerr << "*** Tgrid<" << dim << ">::addBC: Too many (" << Nbcs << " boundary conditions\n";
		exit(1);
	}
	boundgeoms[Nbcs] = geom.copy();
	bound_types[Nbcs] = type;
	bound_funcs[Nbcs] = ev;
	if (type != DIRICHLET_BC && type != CUSTOM_BC && ev != 0) {
		warn("Non-null evaluator ignored since neither DIRICHLET_BC nor CUSTOM_BC");
		bound_funcs[Nbcs] = 0;
	}
	is_box_bc[Nbcs] = isbox;
	Nbcs++;
}

template <smallnat dim>
void Tgrid<dim>::mature_cells(TGridIndex bottom, TGridIndex top, bool& were_ghosts)
/*
- N cells are added "simultaneously" where N=1 for recoarsening,
  N=8 (in 3D) for subdivision, and N=N for constructor phase
- Those cells which are inside some boundary body are flagged
  as ghost cells, others are flagged as interior cells.
  The index of the applying boundary body is recorded in each cell.
- Those that were flagged as ghost cells are scanned.
  If no neighbour is interior cell, the cell is turned from ghost to dead.
- Ghost cells are scanned and inserted in the ghost cell list, if a
  ghost cell list is in use.
*/
{
	TGridIndex i,j;
	int k;
	smallnat d,dir;
	were_ghosts = false;
	for (i=bottom; i<=top; i++) {
		if (!isleaf_vfn(i)) {cerr << "*** mature_cells: not leaf cell\n"; doabort();}
		bool isghost = false;
		Tdimvec X;
		centroid(i,X);
		for (k=0; k<Nbcs; k++)
			if (boundgeoms[k]->IsInsideBody(X)) {
				isghost = true;
				were_ghosts = true;
				break;
			}
		set_celltype(i, isghost ? GHOST_CELL : INTERIOR_CELL);
		if (isghost) set_BCindex(i,k);
	}
	for (i=bottom; i<=top; i++) {
		if (celltype(i) == GHOST_CELL) {
			bool neighbour_is_interior = false;
			for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
				j = dogetneighbour(i,d,dir);
				if (j == NOINDEX) continue;
				if (isleaf_vfn(j)) {
					if (celltype(j) == INTERIOR_CELL) {
						neighbour_is_interior = true;
						goto ende;
					} else {
					}
				} else {
					TGridIndex j1;
					const TGridIndex k = firstchild_vfn(j);
					smallnat ch;
					for (ch=0; ch<(1<<dim); ch++) {
						if (GetBit(d,ch) == dir) continue;
						j1 = k + ch;
						if (celltype(j1) == INTERIOR_CELL) {
							neighbour_is_interior = true;
							goto ende;
						} else {
						}
					}
				}
			}
		ende:
			if (!neighbour_is_interior) {
				set_celltype(i,DEAD_CELL);
			}
		}
	}
	if (cell_prepare_hook)
		for (i=bottom; i<=top; i++) if (celltype(i) != DEAD_CELL) (*cell_prepare_hook)(*this,i);
}

template <smallnat dim>
void Tgrid<dim>::reprepare_all_cells()
{
	if (!cell_prepare_hook) return;
	TGridIndex i;
	const TGridIndex ncells = Ncells_vfn();
	int cnt = 0;
	for (i=0; i<ncells; i++)
		if ((celltype(i) == INTERIOR_CELL || celltype(i) == GHOST_CELL) && isleaf_vfn(i)) {
			(*cell_prepare_hook)(*this,i);
			cnt++;
		}
	cout << "reprepare_all_cells: called cell_prepare_hook " << cnt << " times (ncells=" << ncells << ")\n";
}

template <smallnat dim>
void Tgrid<dim>::applyBC(TGridIndex i, real t)
{
	clearcache();
	const TBCtype type = boundary_type(i);
	if (type == DIRICHLET_BC) {
		fgetcell(i);
		const TEvaluatorFunctionPtr ev = boundary_func(i);
		if (ev) {
			smallnat d,a;
			Tdimvec X1;
			real X[MAXDIM][VECLEN];
			centroid(i,X1);
			for (d=0; d<dim; d++) X[d][0] = X1(d);
			for (d=dim; d<MAXDIM; d++) X[d][0] = 0;
			for (a=0; a<ncd; a++) CT(a) = 0;
			TGridIndexVector iv(1);
			iv[0] = i;
			(*ev)(*this,iv,t,X,smallnat(1),0);
			putcell(i);
		}
		return;
	}
	if (type == CUSTOM_BC) {
		const TEvaluatorFunctionPtr ev = boundary_func(i);
		if (ev) {
			smallnat d,a;
			Tdimvec X1;
			real X[MAXDIM][VECLEN];
			centroid(i,X1);
			for (d=0; d<dim; d++) X[d][0] = X1(d);
			for (d=dim; d<MAXDIM; d++) X[d][0] = 0;
			for (a=0; a<ncd; a++) CT(a) = 0;
			TGridIndexVector iv(1);
			iv[0] = i;
			(*ev)(*this,iv,t,X,smallnat(1),0/*u*/);
			putcell(i);
		} else
			warn("CUSTOM_BC without evaluator, not modifying ghost");
		return;
	}

	// i is a ghost cell, u will be the average of neighbouring interior cells
	real u[max_ncd],u1[max_ncd];
	smallnat a,d,dir;
	TGridIndex j;
	for (a=0; a<ncd; a++) u[a] = u1[a] = 0;
	real norm;
	norm = 0.0;
	for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
		j = dogetneighbour(i,d,dir);
		if (j == NOINDEX) continue;
		const int result1 = child_interior_average(j,u1);
		if (result1 > 0) {
			norm = norm + 1;
			for (a=0; a<ncd; a++) u[a]+= u1[a];
		}
	}
	if (norm == 0) {
		// We consider this a normal case now. Ghosts that lose contact with interior cells
		// because of subdivision are detected here, not in subdivide() as earlier.
		// In subdivide() we did this by calling mature_all_children for all neighbours,
		// but this causes problems in parallel adaptation.
		set_celltype(i,DEAD_CELL);
		return;
	}
	norm = 1.0/norm;
	for (a=0; a<ncd; a++) u[a]*= norm;
	// Previously CUSTOM_BC branch was here (it used u)
	
	if (type == WALL_BC) {
		Tdimvec Xc;
		centroid(i,Xc);
		real nx,ny,nz;
		boundary_geometry(i)->Normal(Xc, nx,ny,nz);
		Tdimvec n;
		smallnat d1;
		for (d1=0; d1<MAXDIM; d1++) n[d1] = 0;
		n[0] = nx;
		if (dim >= 2) n[1] = ny;
		if (dim >= 3) n[2] = nz;
		if (use_mapping) {
			real Xin[MAXDIM][VECLEN], Xout[MAXDIM][VECLEN];
			real XMAX = 0;
			for (d1=0; d1<dim; d1++) XMAX = max(XMAX, max(fabs(xmin[d1]), fabs(xmax[d1])));
			const real epsilon = sqrt(machine::epsilon())*XMAX;
			switch (dim) {
			case 1: break;
			case 2:
		    {
				Tdimvec t;
				// t = ez x n
				t[0] = -n(1);
				t[1] = n(0);
				for (d1=0; d1<dim; d1++) {
					Xin[d1][0] = Xc(d1);
					Xin[d1][1] = Xc(d1) + epsilon*t(d1);
				}
				(*map)(Xin,Xout,2);
				for (d1=0; d1<dim; d1++) t[d1] = Xout[d1][1] - Xout[d1][0];
				// now t is the mapped unnormalized tangent vector
				// n = t x ez
				n[0] = t(1);
				n[1] = -t(0);
				n.normalize();
			}
			break;
			case 3:
	 	    {
				Tvec N[3], T1[3], T2[3];
				for (d1=0; d1<3; d1++) N[d1][0] = n(d1);
				genbase(N,T1,T2,smallnat(1));
				for (d1=0; d1<dim; d1++) {
					Xin[d1][0] = Xc(d1);
					Xin[d1][1] = Xc(d1) + epsilon*T1[d1](0);
					Xin[d1][2] = Xc(d1) + epsilon*T2[d1](0);
				}
				(*map)(Xin,Xout,3);
				Tdimvec t1,t2;
				for (d1=0; d1<dim; d1++) {
					t1[d1] = Xout[d1][1] - Xout[d1][0];
					t2[d1] = Xout[d1][2] - Xout[d1][0];
				}
				n[0] = t1(1)*t2(2) - t1(2)*t2(1);
				n[1] = t1(2)*t2(0) - t1(0)*t2(2);
				n[2] = t1(0)*t2(1) - t1(1)*t2(0);
				n.normalize();
			}
			break;
			}
		}
		if (ncd >= 5) {
			// We want zero normal velocity. Do this by negating the normal component of p=rhov.
			// That is, we replace p by p - 2*(p.n)n.
			// Notice that this does not change p^2 ==> no need to call toPrimitive/fromPrimitive.
			real p_dot_n = 0;
			for (d1=0; d1<dim; d1++) p_dot_n+= u[1+d1]*n(d1);
			for (d1=0; d1<dim; d1++) u[1+d1] -= 2*p_dot_n*n(d1);
			if (ncd >= 8) {
				// In case of MHD, we want a perfectly conducting boundary ==> zero normal B.
				// As for p above, we replace B by B - 2*(B.n)n.
				real B_dot_n = 0;
				for (d1=0; d1<dim; d1++) B_dot_n+= u[5+d1]*n(d1);
				for (d1=0; d1<dim; d1++) u[5+d1] -= 2*B_dot_n*n(d1);
			}
		}
	}
	for (a=0; a<ncd; a++) CT(a) = u[a];
	putcell(i);
}


void BlockWorkDivide(TGridIndex& a, TGridIndex& b)
// Given an integer range a..b of iterations, distribute the iterations as evenly as possible
// to Npes available processors, where the current process-ID is mype (0<=mype<Npes).
// Npes,mype are external int variables which are initialized in gengrid.C:machine::init()
// (or, they are #defined as 1 and 0, if USE_SHMEM is 0.)
// a and b will be modified as the new range for processor mype.
// Calling this routine once in each PE is enough to reduce the global range a..b
// to the local range.
{
	if (b < a) {
		return;
	}
	const TGridIndex n = b-a+1;		// number of iterations (amount of work) to distribute
	const double n1 = double(n)/double(Npes);
	const TGridIndex a1 = a + TGridIndex(mype*n1 + 0.5);
	const TGridIndex b1 = a + TGridIndex((mype+1)*n1 + 0.5) - 1;
	// Properties:
	// (1) a1next = a1/.mype->mype+1 = a + round((mype+1)*n1) = b1+1 by definition,
	//     thus, the range will be always covered without 'holes'
	// (2) b1last = a + round(Npes*n1) - 1 = a + round(n) - 1 = a+n-1 = a+b-a+1-1 = b, q.e.d.
	// (3) If Npes==1 (and mype==0), a1=a, b1 = a + round(n1) - 1 = a + n - 1 = a+b-a+1-1 = b, q.e.d.
	a = a1;
	b = b1;
}


#if HAS_TEMPLATE_INIT_SYNTAX
template class Tgrid<1>;
#if MAXDIM >= 2
template class Tgrid<2>;
#endif
#if MAXDIM >= 3
template class Tgrid<3>;
#endif
#endif


