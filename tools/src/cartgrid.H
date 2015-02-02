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

#ifndef CARTGRID_H

#ifdef __GNUC__
#  pragma interface
#endif

#include "grid.H"
#include "sharr.H"

#ifndef INSTANTIATE_CARTESIAN
#  if defined(_CRAYT3E) || defined(__linux)
#    define INSTANTIATE_CARTESIAN 0
#  else
#    define INSTANTIATE_CARTESIAN 1
#  endif
#endif
// If set to 0, disable Cartesian instantiation and appl/tst2cart and speed up build time somewhat.
// By default we disable Cartesian instantiation for Linux and Cray-T3E to reduce compilation time.
// We enable it on SGI and Cray-C90 because on SGI compilation time doesn't matter and for C90
// the Cartesian case might really be faster in execution.
// In any case, one can override these defaults in makeflags.inc.

template <smallnat dim>
class TCartesianGrid : public Tgrid<dim>
{
protected:
	using Tgrid<dim>::dx;
	using Tgrid<dim>::maxnc;
	using Tgrid<dim>::ntot;
	using Tgrid<dim>::n123;
	using Tgrid<dim>::n1;
	using Tgrid<dim>::n2;
	using Tgrid<dim>::n3;
	using Tgrid<dim>::xmin;
	using Tgrid<dim>::ncd;
	using Tgrid<dim>::ncd_saved;
	using Tgrid<dim>::nsd;
	using Tgrid<dim>::output_realformat;
	using Tgrid<dim>::compose;
	using Tgrid<dim>::decompose;
	using Tgrid<dim>::save1;
	using Tgrid<dim>::surfdataoffset;
	using Tgrid<dim>::areadataoffset;
	using Tgrid<dim>::use_mapping;
	using Tgrid<dim>::bound_funcs;
	using Tgrid<dim>::bound_types;
	using Tgrid<dim>::boundgeoms;
	using Tgrid<dim>::clen;
	using Tgrid<dim>::cornercache;
private:
	TSharedArray<real> /*TRealPtr*/ data;
	TSharedArray<TGridIndex> cellinfotab;
	TRealPtr invvol;			// 1.0/cell_volume for every cell
	TRealPtr physsize;			// The physical (mapped) cell size for every cell
	void vdecompose(const TGridIndexVector& iv, int ijk[3][VECLEN]) const;
	void vcompose(TGridIndexVector& iv, should_be_const int ijk[3][VECLEN], int vlen) const;
	void vnext1(TGridIndexVector& iv, TGridIndex start, TGridIndex stride, Ttraversal trav=TRAV_ALL) const;
#	ifdef TRANSPOSED_STORAGE
	TGridIndex flatindex(int cell, int comp) const {return comp*maxnc + cell;}
	TGridIndex compstep() const {return maxnc;}		// flatindex distance between consecutive comp elements
#	else
	TGridIndex flatindex(int cell, int comp) const {return cell*clen + comp;}
	TGridIndex compstep() const {return 1;}
#	endif
	real M(int cell, int comp) const {return data(flatindex(cell,comp));}
	void Mset(int cell, int comp, real val) {data.put(flatindex(cell,comp),val);}
protected:
	void setup_after_regular();
public:
	using Tgrid<dim>::first;
	using Tgrid<dim>::CT;
	using Tgrid<dim>::writeinfo;
	using Tgrid<dim>::ST;
	using Tgrid<dim>::AT;
	using Tgrid<dim>::FT;
	using Tgrid<dim>::UT;
	using Tgrid<dim>::NN;
	TCartesianGrid() : Tgrid<dim>() {
		n1=n2=n3=ntot=n123[0]=n123[1]=n123[2]=1; dx=0.01;
	}
	void init(smallnat ncd1, smallnat nsd1, smallnat nfq1, TGridIndex maxnc1, bool use_mapping1);
	TCartesianGrid(smallnat ncd1, smallnat nsd1, smallnat nfq1, TGridIndex maxnc1=10000, bool use_mapping1=false) {
		init(ncd1,nsd1,nfq1,maxnc1,use_mapping1);
	}
	void PEcoherency(int rootpe=0);
	// Some dummy functions
	TGridIndex parent(TGridIndex) const {return NOINDEX;}
	TGridIndex level(TGridIndex) const {return 0;}
	TGridIndex child(TGridIndex, smallnat) const {return NOINDEX;}
	bool can_subdivide() const {return false;}
	bool children_have_nonleaf_neighbours(TGridIndex) const {return false;}
	TGridIndex subdivide(TGridIndex, bool =true) {return NOINDEX;}
	bool recoarsen(TGridIndex, bool =true) {return false;}
	bool all_children_are_leaf(TGridIndex) const {return false;}
	void v_all_children_are_leaf(const TGridIndexVector&, TGridIndexVector&) const {}
	bool atomically_test_and_set_mark(TGridIndex) {return false;}
	// For iterating
	TGridIndex next(TGridIndex i) const {return i+Npes;}
	TGridIndex next_1PE(TGridIndex i) const {return i+1;}
	bool isover(TGridIndex i) const {return i>= ntot;}
	TGridIndex Ncells() const {return ntot;}
	virtual TGridIndex Ncells_vfn() const {return ntot;}
	TGridIndex allocated_mem_as_units_of_cell() const {return ntot;}
	// Setting cell volume for each cell
	void setvolume(TGridIndex i, real V) {invvol[i] = 1.0/V;}
	void setinvvolume(TGridIndex i, real invV) {invvol[i] = invV;}
	real cellvolume(TGridIndex i) const {return 1.0/invvol[i];}
	real invcellvolume(TGridIndex i) const {return invvol[i];}
	// setphyssize is relevant only if mapping is in use
	void setphyssize(TGridIndex i, real s) {if (use_mapping) physsize[i] = s;}
	real cellsize_phys(TGridIndex i) const {
		return use_mapping ? physsize[i] : cellsize_vfn(i);
	}

	TCellInfoType cellinfo(TGridIndex i) const {return cellinfotab(i);}
	void set_cellinfo(TGridIndex i, TCellInfoType t) {cellinfotab.put(i, t);}
	bool may_subdivide(TGridIndex i) const {return TCellInfo::may_subdivide(cellinfotab(i));}
	bool may_recoarsen(TGridIndex i) const {return TCellInfo::may_recoarsen(cellinfotab(i));}
	void set_norecoars(TGridIndex i) {cellinfotab.put(i, TCellInfo::set_norecoars(cellinfotab(i)));}
	TCellType celltype(TGridIndex i) const {return TCellInfo::celltype(cellinfotab(i));}
	void set_celltype(TGridIndex i, TCellType ct) {cellinfotab.put(i, TCellInfo::set_celltype(cellinfotab(i),ct));}
	smallnat timeclass(TGridIndex i) const {return TCellInfo::timeclass(cellinfotab(i));}
	void set_timeclass(TGridIndex i, smallnat tc) {cellinfotab.put(i, TCellInfo::set_timeclass(cellinfotab(i),tc));}
	bool isdense(TGridIndex i, smallnat d, smallnat dir) const {return TCellInfo::isdense(cellinfotab(i),d,dir);}
	bool any_neighbour_dense(TGridIndex i) const {return TCellInfo::any_neighbour_dense(cellinfotab(i));}
	bool isindomain(TGridIndex i) const {return celltype(i)==INTERIOR_CELL;}
	TBCtype boundary_type(TGridIndex i) const {return bound_types[TCellInfo::BCindex(cellinfotab(i))];}
	TEvaluatorFunctionPtr boundary_func(TGridIndex i) const {return bound_funcs[TCellInfo::BCindex(cellinfotab(i))];}
	TBoundaryGeometry *boundary_geometry(TGridIndex i) const {return boundgeoms[TCellInfo::BCindex(cellinfotab(i))];}
	smallnat BCindex(TGridIndex i) const {return TCellInfo::BCindex(cellinfotab(i));}
	void set_BCindex(TGridIndex i, smallnat bcind) {cellinfotab.put(i, TCellInfo::set_BCindex(cellinfotab(i),bcind));}
	bool has_calculable_right_surf(TGridIndex i, smallnat d, smallnat& nnei, smallnat& tci) const;

    // Access functions
	real cellsize(TGridIndex i) const;
	real cellsize_vfn(TGridIndex i) const {return cellsize(i);}
	TGridIndex level_vfn(TGridIndex) const {return 0;}
	TGridIndex Nlevels() const {return 1;}
	smallnat Nneighbours_vfn(TGridIndex,smallnat,smallnat) const {return 1;}
	TGridIndex firstchild_vfn(TGridIndex) const {return NOINDEX;}
	smallnat Nneighbours(TGridIndex, smallnat, smallnat) const {return 1;}
	void update_cached(TGridIndex) {}
	void update_cached_neighbours(TGridIndex) {}
	bool isleaf(TGridIndex) const {return true;}
	bool isleaf_vfn(TGridIndex) const {return true;}
	bool isactive(TGridIndex i) const {return isleaf(i) && isindomain(i);}
	bool isactive_vfn(TGridIndex i) const {return isactive(i);}
	int child_interior_average(TGridIndex i, real uave[max_ncd]);
	void centroid(TGridIndex i, Tdimvec& X) const;
	smallnat order(TGridIndex i) const;
	real MinimumGridSpacing() const {return dx;}
	// Scalar access functions:
	TGridIndex dogetneighbour(TGridIndex i, smallnat d, smallnat dir) const;
	TGridIndex getneighbour(TGridIndex i, smallnat d, smallnat dir) const {
		return dogetneighbour(i,d,dir);
	}
	void vdogetneighbour(TGridIndexVector& result,
						 const TGridIndexVector& iv,
						 const smallnat d[VECLEN], const smallnat dir[VECLEN]) const;
	void vgetneighbour_dir0(TGridIndexVector& result,
							const TGridIndexVector& iv,
							const smallnat d[VECLEN]) const;
	void vdogetneighbour_s(TGridIndexVector& result,
						   const TGridIndexVector& iv,
						   smallnat d, smallnat dir) const;
	void vgetneighbour(TGridIndexVector& result,
					   const TGridIndexVector& iv,
					   const smallnat d[VECLEN], const smallnat dir[VECLEN]) const {
		vdogetneighbour(result,iv,d,dir);
	}
	void vgetneighbour_s(TGridIndexVector& result,
						 const TGridIndexVector& iv,
						 smallnat d, smallnat dir) const {
		vdogetneighbour_s(result,iv,d,dir);
	}
	void getcell(TGridIndex i, bool FTflag=true, smallnat v=0);
	void getcell(TGridIndex i, smallnat min_d, smallnat max_d, smallnat v=0);
	void fgetcell(TGridIndex i, smallnat v=0);
	void getnormal(TGridIndex i, TSurfSpec ss, real n3[3]);
	void putcell(TGridIndex i, smallnat v=0);
	void getsurf(TGridIndex i, TSurfSpec ss, bool STflag=false, smallnat v=0);
	void putsurf(TGridIndex i, TSurfSpec ss, smallnat v=0);
#if 0
	// "Syntactic sugar"
	TCartesianGrid& operator>(TGridIndex i) {fgetcell(i); return *this;}
	TCartesianGrid& operator>(const TSurfDef& sd) {getsurf(sd.i,sd.ss); return *this;}
	TCartesianGrid& operator>(const TGridIndexVector& iv) {vgetcell(iv,false); return *this;}
	TCartesianGrid& operator>(const TSurfDefVector& sdv) {vgetsurf(sdv.iv,sdv.sv); return *this;}
#endif
	// Faster versions of vector access functions that override the defaults:
	void vputcell(const TGridIndexVector& iv);
	void vgetcell(const TGridIndexVector& iv, bool FTflag=true);
	void vgetcell_allow_noindex(const TGridIndexVector& iv);
	void vputsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv);
	void vgetsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv, bool STflag=false);
	bool streamload(istream& i, const Theader *hp=0);
	bool streamsave(ostream& o, const char *filename_base=0) const;
	void write_meminfo(ostream&) const {}
	TGridIndex find(const Tdimvec& X) const;
	void forbid_recoarsen() {}
	void clearcache() {}
	bool iseven(TGridIndex i) const;
	int neighbour_timeclass_fix(TGridIndex i);
};

template <smallnat dim>
inline void TCartesianGrid<dim>::getcell(TGridIndex i, bool FTflag, smallnat v)
{
	if (FTflag)
		getcell(i,0,dim-1,v);
	else
		fgetcell(i,v);
}

template <smallnat dim>
inline real TCartesianGrid<dim>::cellsize(TGridIndex) const
{
	return dx;
}

template <smallnat dim>
inline int TCartesianGrid<dim>::child_interior_average(TGridIndex i, real uave[max_ncd])
{
	if (celltype(i) == INTERIOR_CELL) {
		smallnat a;
		fgetcell(i,0);
		for (a=0; a<ncd; a++) uave[a] = CT(a);
		return 1;
	} else
		return 0;
}

template <smallnat dim>
inline void TCartesianGrid<dim>::centroid(TGridIndex i, Tdimvec& X) const
{
	int ijk[3];
	decompose(i, ijk);
	smallnat d;
	for (d=0; d<dim; d++)
		X[d] = xmin[d] + (ijk[d] - NB + 0.5)*dx;
}

template <smallnat dim>
inline smallnat TCartesianGrid<dim>::order(TGridIndex i) const
{
	int ijk[3];
	smallnat result;
	decompose(i, ijk);
	switch (dim) {
	case 1: result = smallnat(ijk[0] % 3); break;
	case 2: result = smallnat((ijk[0] + 2*ijk[1]) % 5); break;
	case 3: result = smallnat((ijk[0] + 2*ijk[1] + 3*ijk[2]) % 7); break;
	}
	return result;
}

template <smallnat dim>
inline bool TCartesianGrid<dim>::has_calculable_right_surf(TGridIndex i, smallnat d, smallnat& nnei, smallnat& tci) const
// Return true if either (a) j is interior cell or (b) j is ghost but i is interior, otherwise false.
// Here j is the right neighbour of i.
// If returning true, set nnei to 1.
// If returning false, nnei is not modified.
{
	TGridIndex j = getneighbour(i,d,1);
	if (j == NOINDEX) return false;
	const TCellInfoType jinfo = cellinfo(j);
	const TCellType jtype = TCellInfo::celltype(jinfo);
	const TCellInfoType iinfo = cellinfo(i);
	if (jtype != INTERIOR_CELL && (jtype != GHOST_CELL || TCellInfo::celltype(iinfo) != INTERIOR_CELL)) return false;
	nnei = 1;
	tci = TCellInfo::timeclass(iinfo);
	return true;
}

template <smallnat dim>
inline bool TCartesianGrid<dim>::iseven(TGridIndex i) const
{
	int ijk[3];
	decompose(i, ijk);
	smallnat d;
	TGridIndex sumindex = 0;
	for (d=0; d<dim; d++) sumindex+= ijk[d];
	return (sumindex & 0x1) == 0;
}

#define CARTGRID_H

#endif
