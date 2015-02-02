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

#ifndef HCGRID_H

#ifdef __GNUC__
#  pragma interface
#endif

#include "mempool.H"

/*
 * 4+dim integer fields:
 *   Pointer to first child, or NOINDEX if leaf cell
 *   Pointer to parent cell, or NOINDEX if basegrid cell
 *   Child order (0..2^dim-1), or 0 if basegrid cell
 *   Cell level (0 for basegrid cell, 1 for its children, etc.)
 *   Dim pointers to a-blocks, or NOINDEX
 *   If CACHE_NEIGHBOURS is defined as 1, then we have 2*dim additional
 *   fields (total is then 4+3*dim).
 */

#define nchildren (1 << dim)

template <smallnat dim>
class THCgrid : public Tgrid<dim> {
protected:
	using Tgrid<dim>::dx;
	using Tgrid<dim>::n1;
	using Tgrid<dim>::n2;
	using Tgrid<dim>::n3;
	using Tgrid<dim>::ntot;
	using Tgrid<dim>::n123;
	using Tgrid<dim>::xmin;
	using Tgrid<dim>::xmax;
	using Tgrid<dim>::maxnc;
	using Tgrid<dim>::clen;
	using Tgrid<dim>::maxnei;
	using Tgrid<dim>::output_realformat;
	using Tgrid<dim>::nsd_priv;
	using Tgrid<dim>::cornercache;
	using Tgrid<dim>::merge_load;
	using Tgrid<dim>::remove_gaps_when_saving;
	using Tgrid<dim>::use_mapping;
	using Tgrid<dim>::parallel_io;
	using Tgrid<dim>::boundgeoms;
	using Tgrid<dim>::bound_types;
	using Tgrid<dim>::bound_funcs;
	using Tgrid<dim>::map;
	using Tgrid<dim>::ncd;
	using Tgrid<dim>::ncd_saved;
	using Tgrid<dim>::nsd;
	using Tgrid<dim>::compose;
	using Tgrid<dim>::decompose;
	using Tgrid<dim>::save1;
	using Tgrid<dim>::surfdataoffset;
	using Tgrid<dim>::areadataoffset;
	using Tgrid<dim>::normaldataoffset;
	using Tgrid<dim>::cell_prepare_hook;
	using Tgrid<dim>::volumedataoffset;
	using Tgrid<dim>::physsizedataoffset;
private:
	Tmempool pool;
	TSharedArray<TVoidPtr> ghostptrs;
#if USE_SHMEM
	TSharedArray<TGridIndex> ghostptr_PEs;
#endif
	int find_k_table[MAXDIM][2][1 << MAXDIM];
	void vnext1(TGridIndexVector& iv, TGridIndex start, TGridIndex stride, Ttraversal trav=TRAV_ALL) const;
	/*
	 * For freelist saving to work, parent index must be 0.
	 * Also, if you change indices here you must check args to pool.streamload/streamsave.
	 */
	// Integer fields in IM(:,:):
	// 0 - parent, 1 - firstchild, 2 - childorder, 3 - level, 4 - cellinfo,
	// 5,6,7 - a_block pointers
	// If PARALLEL_ADAPTATION, then 5th field is the "marked" flag and 6,7,8 are the a_block pointers
	int parent_offset() const {return 0;}
	int firstchild_offset() const {return 1;}
	int childorder_offset() const {return 2;}
	int level_offset() const {return 3;}
	int cellinfo_offset() const {return 4;}
#	if PARALLEL_ADAPTATION
	int mark_offset() const {return 5;}
	void markset(TGridIndex i) {pool.IMset(i,mark_offset(),1);}
	void markclear(TGridIndex i) {pool.IMset(i,mark_offset(),0);}
public:
	/* NOTICE! If you modify ablockpointer_offset, you MUST do the corresponding change in
	   mempool.C:streamload() also! */
	int ablockpointer_offset() const {return 6;}
	bool ismarked(TGridIndex i) const {return bool(pool.IM(i,mark_offset()));}
	bool atomically_test_and_set_mark(TGridIndex i) {return pool.IMatomically_read_and_set(i,mark_offset(),1);}
#	else
public:
	/* NOTICE! If you modify ablockpointer_offset, you MUST do the corresponding change in
	   mempool.C:streamload() also! */
	int ablockpointer_offset() const {return 5;}
	void markclear(TGridIndex) {}
	bool ismarked(TGridIndex) {return false;}
	bool atomically_test_and_set_mark(TGridIndex) {return false;}
#	endif
#if USE_CSHMAT
public:
	TGridIndex parentC(TCacheIndex iC) const {return pool.IMC(iC,parent_offset());}
	TGridIndex levelC(TCacheIndex iC) const {return pool.IMC(iC,level_offset());}
	TGridIndex childC(TCacheIndex iC, smallnat ch) const {return firstchildC(iC)+ch;}
private:
	TGridIndex firstchildC(TCacheIndex iC) const {return pool.IMC(iC,firstchild_offset());}
	TGridIndex childorderC(TCacheIndex iC) const {return pool.IMC(iC,childorder_offset());}
	TGridIndex a_blockC(TCacheIndex iC, smallnat d) const {return pool.IMC(iC,ablockpointer_offset()+d);}
#endif
public:
	TGridIndex parent(TGridIndex i) const {return pool.IM(i,parent_offset());}
	TGridIndex level(TGridIndex i) const {return pool.IM(i,level_offset());}
	TGridIndex level_vfn(TGridIndex i) const {return level(i);}
	TGridIndex Nlevels() const {return 20;}		// CHANGE THIS!!!
	TGridIndex child(TGridIndex i, smallnat ch) const {return firstchild(i)+ch;}
	TGridIndex firstchild(TGridIndex i) const {return pool.IM(i,firstchild_offset());}
	TGridIndex firstchild_vfn(TGridIndex i) const {return firstchild(i);}
private:
	TGridIndex childorder(TGridIndex i) const {return pool.IM(i,childorder_offset());}
	TGridIndex a_block(TGridIndex i, smallnat d) const {return pool.IM(i,ablockpointer_offset()+d);}
	void parentset(TGridIndex i, TGridIndex val) {pool.IMset(i,parent_offset(),val);}
	void firstchildset(TGridIndex i, TGridIndex val) {pool.IMset(i,firstchild_offset(),val);}
	void childorderset(TGridIndex i, TGridIndex val) {pool.IMset(i,childorder_offset(),val);}
	void levelset(TGridIndex i, TGridIndex val) {pool.IMset(i,level_offset(),val);}
	void a_block_set(TGridIndex i, smallnat d, TGridIndex val) {pool.IMset(i,ablockpointer_offset()+d,val);}
#if CACHE_NEIGHBOURS
#	if USE_CSHMAT
	TGridIndex neighbourC(TCacheIndex iC, smallnat d, smallnat dir) const {
		return pool.IMC(iC,ablockpointer_offset()+dim+2*d+dir);
	}
#	endif
	TGridIndex neighbour(TGridIndex i, smallnat d, smallnat dir) const {
		return pool.IM(i,ablockpointer_offset()+dim+2*d+dir);
	}
	void neighbour_set(TGridIndex i, smallnat d, smallnat dir, TGridIndex val) {
		pool.IMset(i,ablockpointer_offset()+dim+2*d+dir,val);
	}
#endif
	void mature_all_children(TGridIndex);
	unsigned char EncodeCellInfo(TGridIndex) const;
	TCellInfoType cellinfotab(TGridIndex i) const {return TCellInfoType(pool.IM(i,cellinfo_offset()));}
#if USE_CSHMAT
	TCellInfoType cellinfotabC(TCacheIndex iC) const {return TCellInfoType(pool.IMC(iC,cellinfo_offset()));}
#endif
	void cellinfotab_put(TGridIndex i, TCellInfoType val) {pool.IMset(i,cellinfo_offset(),TGridIndex(val)); shmemquiet();}
protected:
	virtual void setup_after_regular();
	void update_cached(TGridIndex i, smallnat d, smallnat dir);
	void update_cached(TGridIndex i);
	int ablock_surfdataoffset(smallnat k) const {return k*(nsd+nsd_priv)+nsd_priv;}
	int ablock_areadataoffset(smallnat k) const {return k*(nsd+nsd_priv);}
	int ablock_normaldataoffset(smallnat k) const {return k*(nsd+nsd_priv) + 1;}
public:
	using Tgrid<dim>::first;
	using Tgrid<dim>::CT;
	using Tgrid<dim>::writeinfo;
	using Tgrid<dim>::vfirst;
	using Tgrid<dim>::vnext;
	using Tgrid<dim>::ST;
	using Tgrid<dim>::AT;
	using Tgrid<dim>::FT;
	using Tgrid<dim>::UT;
	using Tgrid<dim>::NN;
	using Tgrid<dim>::applyBC;
	using Tgrid<dim>::mature_cells;
	using Tgrid<dim>::first_1PE;
	THCgrid() {}
	void init(smallnat ncd1, smallnat nsd1, smallnat nfq1, TGridIndex maxnc1, bool use_mapping1);
	THCgrid(smallnat ncd1, smallnat nsd1, smallnat nfq1, TGridIndex maxnc1, bool use_mapping1) {
		init(ncd1,nsd1,nfq1,maxnc1,use_mapping1);
	}
	void PEcoherency(int rootpe=ROOT_PE);
	bool can_subdivide() const {return true;}
	bool has_nonleaf_neighbours(TGridIndex) const;
	bool children_have_nonleaf_neighbours(TGridIndex) const;
	// For iterating
	TGridIndex next(TGridIndex i) const {return i+Npes;}
	TGridIndex next_1PE(TGridIndex i) const {return i+1;}
	bool isover(TGridIndex i) const {return i>= pool.Ncells_fast();}	// assume: PEcoherency persist always when isover is needed
	TGridIndex Ncells() const {return pool.Ncells();}
	virtual TGridIndex Ncells_vfn() const {return pool.Ncells();}
	TGridIndex allocated_mem_as_units_of_cell() const {return pool.allocated_mem_as_units_of_cell();}
	// Access functions
#if USE_CSHMAT
	real cellsizeC(TCacheIndex iC) const {return dx/(1 << levelC(iC));}
#endif
	real cellsize(TGridIndex i) const {return dx/(1 << level(i));}
	real cellsize_vfn(TGridIndex i) const {return cellsize(i);}
	smallnat Nneighbours(TGridIndex i, smallnat d, smallnat dir) const {
		return isdense(i,d,dir) ? (1 << (dim-1)) : 1;
	}
#if USE_CSHMAT
	smallnat NneighboursC(TCacheIndex iC, smallnat d, smallnat dir) const {
		return isdenseC(iC,d,dir) ? (1 << (dim-1)) : 1;
	}
#endif
	smallnat Nneighbours_vfn(TGridIndex i, smallnat d, smallnat dir) const {return Nneighbours(i,d,dir);}
#if USE_CSHMAT
	bool isleafC(TCacheIndex iC) const {return firstchildC(iC) == NOINDEX;}
#endif
	bool isleaf(TGridIndex i) const {return firstchild(i) == NOINDEX;}
	bool isleaf_vfn(TGridIndex i) const {return isleaf(i);}
	bool isactive(TGridIndex i) const {return isleaf(i) && isindomain(i);}
	bool isactive_vfn(TGridIndex i) const {return isactive(i);}
	smallnat MaximumLevel() const;
	real MinimumGridSpacing() const;
    // Setting cell volume for each cell
	void setvolume(TGridIndex i, real V) {pool.Mset(i,volumedataoffset(), 1.0/V);}
	void setinvvolume(TGridIndex i, real invV) {pool.Mset(i,volumedataoffset(), invV);}
	real cellvolume(TGridIndex i) const {return 1.0/pool.M(i,volumedataoffset());}
	real invcellvolume(TGridIndex i) const {return pool.M(i,volumedataoffset());}
	// setphyssize is relevant only if mapping is in use
	void setphyssize(TGridIndex i, real s) {if (use_mapping) pool.Mset(i,physsizedataoffset(), s);}
	real cellsize_phys(TGridIndex i) const {
		return use_mapping ? pool.M(i,physsizedataoffset()) : cellsize_vfn(i);
	}

	TVoidPtr ghostptr(TGridIndex i) const {return ghostptrs(i);}
	void set_ghostptr(TGridIndex i, TVoidPtr val) {ghostptrs.put(i,val);}
#	if USE_SHMEM
	int ghostptr_PE(TGridIndex i) const {return int(ghostptr_PEs(i));}
	void set_ghostptr_PE(TGridIndex i, int pe) {ghostptr_PEs.put(i,pe);}
#	else
	int ghostptr_PE(TGridIndex) const {return 0;}
	void set_ghostptr_PE(TGridIndex,int) {}
#	endif

	TCellInfoType cellinfo(TGridIndex i) const {return cellinfotab(i);}
	void set_cellinfo(TGridIndex i, TCellInfoType t) {cellinfotab_put(i, t);}
	bool may_subdivide(TGridIndex i) const {return TCellInfo::may_subdivide(cellinfotab(i));}
	bool may_recoarsen(TGridIndex i) const {return TCellInfo::may_recoarsen(cellinfotab(i));}
	void set_norecoars(TGridIndex i) {cellinfotab_put(i, TCellInfo::set_norecoars(cellinfotab(i)));}
	TCellType celltype(TGridIndex i) const {return TCellInfo::celltype(cellinfotab(i));}
	void set_celltype(TGridIndex i, TCellType ct) {cellinfotab_put(i, TCellInfo::set_celltype(cellinfotab(i),ct));}
	smallnat timeclass(TGridIndex i) const {return TCellInfo::timeclass(cellinfotab(i));}
	void set_timeclass(TGridIndex i, smallnat tc) {cellinfotab_put(i, TCellInfo::set_timeclass(cellinfotab(i),tc));}
	bool isdense(TGridIndex i, smallnat d, smallnat dir) const {return TCellInfo::isdense(cellinfotab(i),d,dir);}
	bool any_neighbour_dense(TGridIndex i) const {return TCellInfo::any_neighbour_dense(cellinfotab(i));}
	bool isindomain(TGridIndex i) const {return celltype(i)==INTERIOR_CELL;}
	TBCtype boundary_type(TGridIndex i) const {return bound_types[TCellInfo::BCindex(cellinfotab(i))];}
	TEvaluatorFunctionPtr boundary_func(TGridIndex i) const {return bound_funcs[TCellInfo::BCindex(cellinfotab(i))];}
	TBoundaryGeometry *boundary_geometry(TGridIndex i) const {return boundgeoms[TCellInfo::BCindex(cellinfotab(i))];}
	smallnat BCindex(TGridIndex i) const {return TCellInfo::BCindex(cellinfotab(i));}
	void set_BCindex(TGridIndex i, smallnat bcind) {cellinfotab_put(i, TCellInfo::set_BCindex(cellinfotab(i),bcind));}
	bool has_calculable_right_surf(TGridIndex i, smallnat d, smallnat& nnei, smallnat& tci) const;
#if USE_CSHMAT
	TCellInfoType cellinfoC(TCacheIndex iC) const {return cellinfotabC(iC);}
	bool isdenseC(TCacheIndex iC, smallnat d, smallnat dir) const {return TCellInfo::isdense(cellinfotabC(iC),d,dir);}
#endif

    // Scalar access functions:
	TGridIndex dogetneighbour(TGridIndex i, smallnat d, smallnat dir) const;
	void getcell(TGridIndex i, bool FTflag=true, smallnat v=0);
	void getcell(TGridIndex i, smallnat min_d, smallnat max_d, smallnat v=0);
	void fgetcell(TGridIndex i, smallnat v=0);
	void getnormal(TGridIndex i, TSurfSpec ss, real n3[3]);
	void putcell(TGridIndex i, smallnat v=0);
	void getsurf(TGridIndex i, TSurfSpec ss, bool STflag=false, smallnat v=0);
	void putsurf(TGridIndex i, TSurfSpec ss, smallnat v=0);
	// Faster versions of vector access functions that override the defaults:
#if !USE_CSHMAT
	void vgetcell(const TGridIndexVector& iv, bool FTflag=true);
	void vgetcell_allow_noindex(const TGridIndexVector& iv);
	void vgetnormal(const TGridIndexVector& iv, const TSurfSpecVector& sv, Tvec n3[3]);
	void vputcell(const TGridIndexVector& iv);
	void vgetsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv, bool STflag=false);
	void vputsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv);
#endif
#if CACHE_NEIGHBOURS
#	if USE_CSHMAT
	TGridIndex getneighbourC(TCacheIndex iC, smallnat d, smallnat dir) const {
		return neighbourC(iC,d,dir);
	}
#	endif
	TGridIndex getneighbour(TGridIndex i, smallnat d, smallnat dir) const {
		return neighbour(i,d,dir);
	}
	void vgetneighbour(TGridIndexVector& result,
					   const TGridIndexVector& iv,
					   const smallnat d[VECLEN], const smallnat dir[VECLEN]) const;
	void vgetneighbour_s(TGridIndexVector& result,
						 const TGridIndexVector& iv,
						 smallnat d, smallnat dir) const;
#else
	TGridIndex getneighbour(TGridIndex i, smallnat d, smallnat dir) const {
		return dogetneighbour(i,d,dir);
	}
	void vgetneighbour(TGridIndexVector& result,
					   const TGridIndexVector& iv,
					   const smallnat d[VECLEN], const smallnat dir[VECLEN]) const;
	void vgetneighbour_s(TGridIndexVector& result,
						 const TGridIndexVector& iv,
						 smallnat d, smallnat dir) const;
#endif
	int child_interior_average(TGridIndex i, real uave[max_ncd]);
    // recursive average of interior cells, may use CT as workspace, returns result in uave, return number of contributing cells
	void centroid(TGridIndex i, Tdimvec& X) const;
	smallnat order(TGridIndex i) const;
	bool streamload(istream& i, const Theader *hp=0);
	bool streamsave(ostream& o, const char *filename_base=0) const;
	void write_meminfo(ostream& o) const {o << pool;}
	TGridIndex find(const Tdimvec& X) const;
	// HC specific functions (subdivide is actually dummy in cartgrid)
	TGridIndex subdivide(TGridIndex i, real t);
	bool recoarsen(TGridIndex i, real t, bool may_recurse=true);
	void check();
	int find_k(smallnat d, smallnat dir, smallnat ch) const;
	smallnat find_k_fast(smallnat d, smallnat dir, smallnat ch) const;
#if USE_CSHMAT
	bool all_children_are_leafC(TCacheIndex iC) const;
#endif
	bool all_children_are_leaf(TGridIndex i) const;
	void v_all_children_are_leaf(const TGridIndexVector& iv, TGridIndexVector& result) const;
	void forbid_recoarsen();
	void clearcache() {pool.clearcache();}
	bool iseven(TGridIndex i) const;
	int neighbour_timeclass_fix(TGridIndex i);
	virtual ~THCgrid() {}
};

template <smallnat dim>
inline void THCgrid<dim>::getcell(TGridIndex i, bool FTflag, smallnat v)
{
	if (FTflag)
		getcell(i,0,dim-1,v);
	else
		fgetcell(i,v);
}

#if CACHE_NEIGHBOURS

template <smallnat dim>
inline void THCgrid<dim>::vgetneighbour(
	TGridIndexVector& result,
	const TGridIndexVector& iv,
	const smallnat d[VECLEN], const smallnat dir[VECLEN]) const
{
	smallnat v;
	VLOOP(iv) result[v] = neighbour(iv(v),d[v],dir[v]);
	result.setlength(iv.length());
}

template <smallnat dim>
inline void THCgrid<dim>::vgetneighbour_s(
	TGridIndexVector& result,
	const TGridIndexVector& iv,
	smallnat d, smallnat dir) const
{
	smallnat v;
	VLOOP(iv) result[v] = neighbour(iv(v),d,dir);
	result.setlength(iv.length());
}

#else

template <smallnat dim>
inline void THCgrid<dim>::vgetneighbour(TGridIndexVector& result,
										const TGridIndexVector& iv,
										const smallnat d[VECLEN], const smallnat dir[VECLEN]) const
{
	smallnat v;
	VLOOP(iv) result[v] = dogetneighbour(iv(v),d[v],dir[v]);
}

template <smallnat dim>
inline void THCgrid<dim>::vgetneighbour_s(TGridIndexVector& result,
										  const TGridIndexVector& iv,
										  smallnat d, smallnat dir) const
{
	smallnat v;
	VLOOP(iv) result[v] = dogetneighbour(iv(v),d,dir);
}

#endif	/* CACHE_NEIGHBOURS */

#if USE_CSHMAT
template <smallnat dim>
inline bool THCgrid<dim>::all_children_are_leafC(TCacheIndex iC) const
{
	smallnat ch;
	if (isleafC(iC)) return false;
	for (ch=0; ch<nchildren; ch++)
		if (!isleafC(pool.load(childC(iC,ch)))) return false;
	return true;
}
#endif

template <smallnat dim>
inline bool THCgrid<dim>::all_children_are_leaf(TGridIndex i) const
{
	smallnat ch;
	if (isleaf(i)) return false;
	for (ch=0; ch<nchildren; ch++)
		if (!isleaf(child(i,ch))) return false;
	return true;
}

template <smallnat dim>
inline bool THCgrid<dim>::has_calculable_right_surf(TGridIndex i, smallnat d, smallnat& nnei, smallnat& tci) const
// Return true if either (a) j is interior cell or (b) j is ghost but i is interior, otherwise false.
// Here j is the right neighbour of i.
// If returning true, set nnei to the number of right neighbours of i.
// If returning false, nnei,tci are not modified.
{
#	if USE_CSHMAT
	const TCacheIndex iC = pool.load(i);
	TGridIndex j = getneighbourC(iC,d,1);
	if (j == NOINDEX) {
		return false;
	}
	const TCacheIndex jC = pool.load(j);
	const TCellInfoType jinfo = cellinfoC(jC);
	const TCellType jtype = TCellInfo::celltype(jinfo);
	const TCellInfoType iinfo = cellinfoC(iC);
	const TCellType itype = TCellInfo::celltype(iinfo);
	assert(!(itype==INTERIOR_CELL && jtype==DEAD_CELL));
	if (jtype != INTERIOR_CELL && (jtype != GHOST_CELL || itype != INTERIOR_CELL)) {
		return false;
	}
	nnei = NneighboursC(iC,d,1);
#	else
	TGridIndex j = getneighbour(i,d,1);
	if (j == NOINDEX) return false;
	const TCellInfoType jinfo = cellinfo(j);
	const TCellType jtype = TCellInfo::celltype(jinfo);
	const TCellInfoType iinfo = cellinfo(i);
	if (jtype != INTERIOR_CELL && (jtype != GHOST_CELL || TCellInfo::celltype(iinfo) != INTERIOR_CELL)) return false;
	nnei = Nneighbours(i,d,1);
#	endif
	tci = TCellInfo::timeclass(iinfo);
	return true;
}

#define HCGRID_H

#endif
