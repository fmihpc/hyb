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

#ifndef GRID_H

#ifdef __GNUC__
#  pragma interface
#endif

#include "BC.H"
#include "fileheader.H"
#include "vec.H"
#include "shmemintf.H"
#include "sharr.H"
#include "cellcache.H"

#ifndef DEBUG
#  define PROMOTE_INLINING
#endif

#ifndef PARALLEL_ADAPTATION
#  define PARALLEL_ADAPTATION 1
#endif

#define NOINDEX (-1234567)

extern void doabort()
#ifdef __GNUC__
	__attribute__((noreturn))
#endif
;

inline void vselect_nonnegative(TGridIndexVector& result_ind, const TGridIndexVector& src, const TGridIndexVector& ind)
{
	smallnat v,v2;
	VDIRS;
	for (v=v2=0; v<ind.length(); v++)
		if (src(v) >= 0) result_ind[v2++] = ind(v);
	result_ind.setlength(v2);
}

class TSurfSpec {
private:
	// d: 0,1,2 ==> 2 bits enough
	// dir: 0,1 ==> 1 bit enough
	// k: 0,1,2,3 ==> 2 bits enough
#ifdef _CRAY1
	smallnat d0,dir0,k0;
#else
	smallnat d0:2,dir0:1,k0:2;
#endif
public:
	TSurfSpec() {}
	TSurfSpec(smallnat d1, smallnat dir1, smallnat k1) {d0=d1; dir0=dir1; k0=k1;}
	smallnat d() const {return d0;}
	smallnat dir() const {return dir0;}
	smallnat k() const {return k0;}
	void set_d(smallnat d1) {d0=d1;}
	void set_dir(smallnat dir1) {dir0=dir1;}
	void set_k(smallnat k1) {k0=k1;}
};

enum TCellType {INTERIOR_CELL, GHOST_CELL, DEAD_CELL, REMOVED_CELL};

typedef TGridIndex TCellInfoType;	// Integer type larger enough to hold 31 bits

/*
 * We do TCellInfo by explicit bit manipulations rather than using struct bit-fields.
 * The primary reason is that on Cray T3E, sizeof(any struct) is at least 8 bytes,
 * which does not fit in TGridIndex (4 bytes).
 * TCellInfo contains only static member functions for manipulating the individual bits.
 * The data is passed as an explicit argument to these functions.
 * The bit layout is as follows:
 *
 *                | <-- timeclass ----> | |<-------bcindex----->| |  | denseflags| |type
 * 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
 *                                                                ^  ^
 *                                                                |  |
 *                                                                |  nosubdiv
 *                                                                norecoars
 *
 * The type field is the cell type, and takes four possible values (enum TCellType above).
 * We allocate 3 bits for it anyway to have some room for future expansion.
 * denseflags (neighbourdenseflags) has 6 bits for 6 directions in 3D (d=0,1,2, dir=0,1).
 * Bits inside neighbourdenseflags are numbered as: b=2*d+dir, (b=0..5). Bit ON means that the
 * neighbour in that direction exists and is dense. Direction dir=0 is left and dir=1 is right.
 * Coordinate d=0 is X, d=1 is Y and d=2 is Z.
 * Bit 9 is the nosubdiv bit, when it is on, subdivision of that cell is forbidden.
 * Bit 10 is the norecoars bit, when it is on, recoarsening of that cell is forbidden.
 * The norecoars bit is used in body refinement.
 * Bits 11-18 (8 bits) is the bcindex, which is the index of the boundary condition for this cell.
 * This field is meaningful only if the cell is a ghost cell (celltype=GHOST_CELL).
 * Since it contains 8 bits, the number of different boundary conditions (MAX_N_BCS below)
 * is 256.
 * The last field (bits 19-26, 8 bits) is allocated to the timeclass field. Thus, the maximum
 * timeclass is 256, which is enough for any purpose (since the time step is 2^(-timeclass)).
 * Bits 27-31 are currently not used.
 *
 * A zero is a useful default bit pattern.
 * It corresponds to an interior cell with no refined neighbours and with zero timeclass
 * and no subdivision/recoarsening forbidding flags.
 */

#define BCINDEX_BITS 8

class TCellInfo {
public:
	static TCellType celltype(TCellInfoType i) {return TCellType(i & 0x7);}
	static TCellInfoType set_celltype(TCellInfoType i, TCellType ct) {return (i & (~7)) | TCellInfoType(ct);}
	static bool isdense(TCellInfoType i, smallnat d, smallnat dir) {
		const smallnat b = 3 + 2*d + dir;
		return (i & (TCellInfoType(1) << b)) != 0;
	}
	static bool any_neighbour_dense(TCellInfoType i) {return (i & 0x01F8) != 0;}
	static TCellInfoType clear_dense_flags(TCellInfoType i) {return i & (~0x01F8);}
	static TCellInfoType set_dense(TCellInfoType i, smallnat d, smallnat dir, bool flag=true) {
		const smallnat b = 3 + 2*d + dir;
		if (flag) return i | (TCellInfoType(1) << b); else return i & (~(TCellInfoType(1) << b));
	}
	static bool may_subdivide(TCellInfoType i) {return (i & (TCellInfoType(1) << 9)) == 0;}
	static TCellInfoType set_nosubdiv(TCellInfoType i, bool flag=true) {
		if (flag) return i | (TCellInfoType(1) << 9); else return i & (~(TCellInfoType(1) << 9));
	}
	static bool may_recoarsen(TCellInfoType i) {return (i & (TCellInfoType(1) << 10)) == 0;}
	static TCellInfoType set_norecoars(TCellInfoType i, bool flag=true) {
		if (flag) return i | (TCellInfoType(1) << 10); else return i & (~(TCellInfoType(1) << 10));
	}
	static smallnat BCindex(TCellInfoType i) {return smallnat((i & 0x0007F800) >> 11);}
	static TCellInfoType set_BCindex(TCellInfoType i, smallnat bcindex) {
		return (i & (~0x0007F800)) | ((bcindex & 0xF) << 11);
	}
	static smallnat timeclass(TCellInfoType i) {return smallnat((i & 0x07F80000) >> 19);}
	static TCellInfoType set_timeclass(TCellInfoType i, smallnat timeclass) {
		return (i & (~0x07F80000)) | ((timeclass & 0xF) << 19);
	}
};

class TSurfSpecVector {
private:
#ifdef VECTOR_MACHINE
	smallnat d0[VECLEN], dir0[VECLEN], k0[VECLEN];
#else
	TSurfSpec tab[VECLEN];
#endif
	int len;
public:
	TSurfSpecVector() {len=0;}
	void setlength(int n) {assert(0<=n && n<=VECLEN); len=n;}
	TSurfSpecVector(int n) {setlength(n);}
	int length() const {return len;}
#ifdef VECTOR_MACHINE
	smallnat d(int i) const {assert(0<=i && i<len); return d0[i];}
	smallnat dir(int i) const {assert(0<=i && i<len); return dir0[i];}
	smallnat k(int i) const {assert(0<=i && i<len); return k0[i];}
	void set_d(int i, smallnat d) {assert(0<=i && i<VECLEN); d0[i] = d;}
	void set_dir(int i, smallnat dir) {assert(0<=i && i<VECLEN); dir0[i] = dir;}
	void set_k(int i, smallnat k) {assert(0<=i && i<VECLEN); k0[i] = k;}
	TSurfSpec operator()(int i) const {assert(0<=i && i<len); return TSurfSpec(d0[i],dir0[i],k0[i]);}
#else
	smallnat d(int i) const {assert(0<=i && i<len); return tab[i].d();}
	smallnat dir(int i) const {assert(0<=i && i<len); return tab[i].dir();}
	smallnat k(int i) const {assert(0<=i && i<len); return tab[i].k();}
	void set_d(int i, smallnat d) {assert(0<=i && i<VECLEN); tab[i].set_d(d);}
	void set_dir(int i, smallnat dir) {assert(0<=i && i<VECLEN); tab[i].set_dir(dir);}
	void set_k(int i, smallnat k) {assert(0<=i && i<VECLEN); tab[i].set_k(k);}
	TSurfSpec operator()(int i) const {assert(0<=i && i<len); return tab[i];}
#endif
};

#define NB 1

enum TRealFormat {REALFORMAT_INQUIRE, REALFORMAT_DOUBLE, REALFORMAT_FLOAT, REALFORMAT_ASCII};

struct TSurfDef {
	TGridIndex i;
	TSurfSpec& ss;
	TSurfDef(TGridIndex i1, TSurfSpec& ss1) : i(i1), ss(ss1) {}
};

struct TSurfDefVector {
	TGridIndexVector& iv;
	TSurfSpecVector& sv;
	TSurfDefVector(TGridIndexVector& iv1, TSurfSpecVector& sv1) : iv(iv1), sv(sv1) {}
};

inline TSurfDef operator,(TGridIndex i, TSurfSpec& ss) {return TSurfDef(i,ss);}
inline TSurfDefVector operator,(TGridIndexVector& iv, TSurfSpecVector& sv) {return TSurfDefVector(iv,sv);}

enum Ttraversal {TRAV_ALL, TRAV_LEAF_ONLY, TRAV_DOMAIN_ONLY, TRAV_ACTIVE_ONLY, TRAV_LEAFGHOST_ONLY};

//#define MAX_N_BC 64		/* the width of TCellInfo::bcindex bit field MUST agree with this! it must be 6, if MAX_N_BC is 64 */
#define MAX_N_BC (1 << BCINDEX_BITS)

template <smallnat dim>
class Tgrid : public TGenericGrid
{
public:
	typedef void (*TSourceTerm)(Tgrid<dim>& g, should_be_const real X[MAXDIM][VECLEN], real src[max_ncd][VECLEN], smallnat vlen);
private:
	bool regular_called, setbox_called, do_fast_regular;
	bool cache_intpols;				// set to true only if the grid is not modified, like in hcvis and A
	static void default_map(should_be_const real Xin[MAXDIM][VECLEN], real Xout[MAXDIM][VECLEN], smallnat vlen);
protected:
	TCellCache cornercache;		// in use if cache_intpols is true
	int nsd_priv;				// number of "private" surface data items (area, +possibly normal vec)
	int ncd_priv;				// number of "private" cell data items (volume, +possibly physical cell size)
	bool merge_load;			// if true, old mempool structure is not destroyed on streamload
	bool remove_gaps_when_saving;	// save some space of output files by not storing deleted-in-between cells, usually true
	bool use_mapping;			// whether analytic mapping is in use, if it is, normal vectors are included in nsd_priv
	bool parallel_io;			// whether every PE writes its own file for SHMEM machine(s) with Npes>1, do hccomb later to concat
//	TSharedArray<TGridIndex> cellinfotab;
	TGridIndex maxnc;			// Maximum number of cells
	int clen;					// = ncd + dim*(nsd+nsd_priv) + ncd_priv
	smallnat maxnei;			// must be overwritten by all THCgrid constructors
	real xmin[3], xmax[3];		// box size
	int n1,n2,n3;
	int n123[3];
	int ntot;
	real dx;
	TRealFormat output_realformat;
	// Boundary condition stuff:
	int Nbcs;
	TBoundaryGeometry *boundgeoms[MAX_N_BC];
	TBCtype bound_types[MAX_N_BC];
	TEvaluatorFunctionPtr bound_funcs[MAX_N_BC];
	TBCtype box_bctypes[MAXDIM][2];
	TEvaluatorFunctionPtr box_bcfuncs[MAXDIM][2];
	bool is_box_bc[MAX_N_BC];
	void (*map)(should_be_const real Xin[MAXDIM][VECLEN], real Xout[MAXDIM][VECLEN], smallnat vlen);
	void (*cell_prepare_hook)(Tgrid<dim>& g, TGridIndex i);
	TSourceTerm source_term_hook;
	bool save1(ostream& o) const;
#ifdef PROMOTE_INLINING
	virtual void vnext1(TGridIndexVector& iv, TGridIndex start, TGridIndex stride, Ttraversal trav=TRAV_ALL) const=0;
#else
	virtual void vnext1(TGridIndexVector& iv, TGridIndex start, TGridIndex stride, Ttraversal trav=TRAV_ALL) const;
#endif
	TGridIndex compose(const int ijk[3]) const;
	void decompose(TGridIndex p, int ijk[3]) const;
	TGridIndex compose(int i, int j) const {return i*n2+j;}
	TGridIndex compose(int i, int j, int k) const {return (i*n2+j)*n3+k;}
	void corner_deplist(const Tdimvec& X,
						TGridIndex depcells[8], real weights[8], smallnat& Ndepcells) const;
	virtual void setup_after_regular() {}
	// access functions to nsd_priv data items
	int areadataoffset(smallnat d) const {return ncd + d*(nsd+nsd_priv);}				// the first priv item is the area
	int normaldataoffset(smallnat d) const {return ncd + d*(nsd+nsd_priv) + 1;}			// the the normal vector, if any
	int surfdataoffset(smallnat d) const {return ncd + d*(nsd+nsd_priv) + nsd_priv;}	// "public" surf data is beyond "private"
	// access functions to ncd_priv data items
	int volumedataoffset() const {return ncd + dim*(nsd+nsd_priv);}
	int physsizedataoffset() const {return volumedataoffset()+1;}
public:
	Tgrid() {}
	void init(smallnat ncd1, smallnat nsd1, smallnat nfq1, TGridIndex maxnc1, bool use_mapping1);
	Tgrid(smallnat ncd1, smallnat nsd1, smallnat nfq1, TGridIndex maxnc1, bool use_mapping1) {
		init(ncd1,nsd1,nfq1,maxnc1,use_mapping1);
	}
	virtual void PEcoherency(int rootpe=0) = 0;
	// Initializers
	void setbox(real x1, real x2);
	void setbox(real x1, real x2, real y1, real y2);
	void setbox(real x1, real x2, real y1, real y2, real z1, real z2);
	void setmap(void(*map1)(should_be_const real[MAXDIM][VECLEN], real [MAXDIM][VECLEN], smallnat));
	bool mapping_in_use() const {return use_mapping;}
	void set_cell_prepare_func(void (*hook)(Tgrid<dim>&, TGridIndex)) {cell_prepare_hook = hook;}
	void set_source_term(TSourceTerm hook) {source_term_hook = hook;}
	TSourceTerm sourceterm() const {return source_term_hook;}
	// Initializers
	void regular(real dx1);
	void regular(int ntot1);
	void regular(int i1, int i2);
	void regular(int i1, int i2, int i3);
	// Scaling of all dimensions, affects dx,xmin,xmax (box size)
	void scale(real scaling);
	// Using set_fast_regular prevents calling update_cached for every basegrid cell.
	// This is used by Tmetagrid::streamload(), which needs to call regular() only to set
	// n1,n2,n3, ntot etc. The speedup in e.g., vis/A, is quite substantial.
	void set_fast_regular(bool flag=true) {do_fast_regular=flag;}
	// Query function, returns true for HCgrid, false for cartgrid
	virtual bool can_subdivide() const = 0;
	// Setting real format
	TRealFormat realformat(TRealFormat rf=REALFORMAT_INQUIRE);
	// Setting some flags
	void set_merge_load(bool flag=true) {merge_load=flag;}
	void set_remove_gaps(bool flag=true) {remove_gaps_when_saving=flag;}
	void set_intpol_cacheing(bool flag=true) {cache_intpols=flag; if (!flag) cornercache.clear();}
	void clear_intpol_cache() {cornercache.clear();}
	void print_corner_cache_stats(ostream& o) {cornercache.printstats(o);}
	// Grid index iteration
	TGridIndex first() const {return mype;}
	TGridIndex first_1PE() const {return 0;}
	virtual TGridIndex next(TGridIndex) const = 0;
	virtual TGridIndex next_1PE(TGridIndex) const = 0;
	void vfirst(TGridIndexVector& iv, Ttraversal trav=TRAV_ALL) const {vnext1(iv,mype,Npes,trav);}
	void vfirst_1PE(TGridIndexVector& iv, Ttraversal trav=TRAV_ALL) const {vnext1(iv,0,1,trav);}
	void vnext(TGridIndexVector& iv, Ttraversal trav=TRAV_ALL) const {
		if (iv.isempty()) return;
		vnext1(iv,iv(iv.length()-1)+Npes,Npes,trav);
	}
	void vnext_1PE(TGridIndexVector& iv, Ttraversal trav=TRAV_ALL) const {
		if (iv.isempty()) return;
		vnext1(iv,iv(iv.length()-1)+1,1,trav);
	}
	// Helper functions
	smallnat maxneighbours() const {return maxnei;}
	void getbox(real Xmin[3], real Xmax[3]) const;
	void get_exterior_box(real Xmin[3], real Xmax[3]) const;
	smallnat dimension() const {return dim;}
	static smallnat Dimension() {return dim;}
	smallnat Ncelldata() const {return ncd;}
	smallnat Nfluxquant() const {return nfq;}
	TGridIndex MaxNCells() const {return maxnc;}

	virtual TCellType celltype(TGridIndex i) const = 0;
	virtual void set_celltype(TGridIndex i, TCellType ct) = 0;
	virtual smallnat BCindex(TGridIndex i) const = 0;
	virtual void set_BCindex(TGridIndex i, smallnat bcind) = 0;
	virtual TBCtype boundary_type(TGridIndex i) const = 0;
	virtual TEvaluatorFunctionPtr boundary_func(TGridIndex i) const = 0;
	virtual TBoundaryGeometry *boundary_geometry(TGridIndex i) const = 0;

	void set_box_BC(smallnat d, smallnat dir, TBCtype type, TEvaluatorFunctionPtr ev=0);
	void applyBC(TGridIndex i, real t);
#ifndef PROMOTE_INLINING
	virtual TGridIndex Ncells() const = 0;
	virtual bool isover(TGridIndex) const = 0;
	virtual bool isleaf(TGridIndex i) const = 0;
	virtual real cellsize(TGridIndex i) const = 0;
	virtual smallnat Nneighbours(TGridIndex i, smallnat d, smallnat dir) const = 0;
#endif
	virtual TGridIndex Ncells_vfn() const = 0;
	virtual bool isactive_vfn(TGridIndex) const = 0;
	virtual real cellsize_vfn(TGridIndex i) const = 0;
	virtual smallnat Nneighbours_vfn(TGridIndex i, smallnat d, smallnat dir) const = 0;
	virtual TGridIndex level_vfn(TGridIndex i) const = 0;
	virtual TGridIndex Nlevels() const = 0;
	virtual bool isleaf_vfn(TGridIndex i) const = 0;
	virtual TGridIndex firstchild_vfn(TGridIndex i) const = 0;
	virtual int child_interior_average(TGridIndex i, real uave[max_ncd]) = 0;	// uses CT as workspace, thus cannot be const
	virtual void centroid(TGridIndex i, Tdimvec& X) const = 0;
	void writeinfo(ostream& o, TGridIndex i) const;
	virtual real MinimumGridSpacing() const = 0;
	real BasegridCellSize() const {return dx;}
	// Scalar access functions. Last argument (v) is usually zero in user calls.
	virtual TGridIndex dogetneighbour(TGridIndex i, smallnat d, smallnat dir) const = 0;
	virtual void getcell(TGridIndex i, bool FTflag=true, smallnat v=0) = 0;
    // --> fills CT(c,v) [and FT(d,dir,s,v), AT(d,dir,v) if FTflag]
	virtual void getcell(TGridIndex i, smallnat min_d, smallnat max_d, smallnat v=0) = 0;
	// --> same as ordinary getcell but fills FT,AT only for d=min_d..max_d
	virtual void fgetcell(TGridIndex i, smallnat v=0) = 0;
	// --> fills CT(c,v)
	virtual void getnormal(TGridIndex i, TSurfSpec s, real n3[3]) = 0;
	virtual void vgetnormal(const TGridIndexVector& iv, const TSurfSpecVector& sv, Tvec n3[3]);
	virtual void putcell(TGridIndex i, smallnat v=0) = 0;					// <-- requires CT(c,v)
	virtual void getsurf(TGridIndex i, TSurfSpec s, bool STflag=false, smallnat v=0) = 0;
    // --> fills UT(dir,c,v), and ST(0,s,v) if STflag==true
	virtual void putsurf(TGridIndex i, TSurfSpec s, smallnat v=0) = 0;	// <-- requires ST(k,s,v)
	// Vector access functions. These are not pure, default versions call scalar accessors.
	// Therefore we hade the v argument also for scalar case.
	virtual void vgetcell(const TGridIndexVector& iv, bool FTflag=true);
    // --> fills CT(c,v) [and FT(d,dir,s,v), AT(d,dir,v) if FTflag]
	virtual void vgetcell_allow_noindex(const TGridIndexVector& iv);		// Same as vgetcell but allow and ignore NOINDEX
	virtual void vputcell(const TGridIndexVector& iv);						// <-- requires CT(c,v) [v=0..iv.length()-1]
	virtual void vgetsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv, bool STflag=false);
	// --> fills UT(dir,c,v), and ST(0,s,v) if STflag==true
	virtual void vputsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv);	// <-- requires ST(k,s,v)
	void addBC(const TBoundaryGeometry& geom, TBCtype type, bool isbox=false, TEvaluatorFunctionPtr ev=0);
	bool is_box_ghost(TGridIndex i) const {return celltype(i)==GHOST_CELL && is_box_bc[BCindex(i)];}
	void mature_cells(TGridIndex bottom, TGridIndex top, bool& were_ghosts);
	void reprepare_all_cells();		// call when you want to recompute background magnetic field (calls cell_prepare_hook)
	virtual bool streamload(istream& i, const Theader *hp=0) = 0;
	virtual bool streamsave(ostream& o, const char *filename_base=0) const = 0;
	bool load(const char *fn);
	bool save(const char *fn, int precision=10, ios::openmode om=ios::out) const;
	bool save_append(const char *fn, int precision=10) const {return save(fn,precision,ios::app);}
	void set_parallel_IO(bool flag=true) {parallel_io = flag;}
	virtual TGridIndex find(const Tdimvec& X) const = 0;
	bool intpol(const Tdimvec& X, int order=0, bool quiet=false);
	void vintpol(const real X[MAXDIM][VECLEN], smallnat vlen, int order=0);
	real intpol(const Tdimvec& X, int order, smallnat comp) {intpol(X,order); return CT(comp);}
	// Syntactic sugar for grid input/output
	friend ostream& operator<<(ostream& o, const Tgrid& g) {g.streamsave(o); return o << flush;}
	friend istream& operator>>(istream& i, Tgrid& g) {g.streamload(i); return i;}
	// Function to write memory info to stream
	virtual void write_meminfo(ostream& o) const = 0;
	// Other
	virtual void forbid_recoarsen() = 0;
	virtual void clearcache() = 0;
	virtual int ablockpointer_offset() const {return 0;}	// not very important, needed in gridsimu.C:run() maxmem stuff
	virtual ~Tgrid();
};

template <smallnat dim>
inline TGridIndex Tgrid<dim>::compose(const int ijk[3]) const
{
	if (dim == 1)
		return ijk[0];
	else if (dim == 2)
		return ijk[0]*n2 + ijk[1];
	else
		return (ijk[0]*n2 + ijk[1])*n3 + ijk[2];
}

template <smallnat dim>
inline void Tgrid<dim>::decompose(TGridIndex p, int ijk[3]) const
{
	if (dim == 1) {
		ijk[0] = p;
	} else if (dim == 2) {
		// p = i*n2 + j
		const int j = p % n2;
		const int i = (p-j)/n2;
		ijk[0] = i;
		ijk[1] = j;
	} else {
		// 3D: p = (i*n2 + j)*n3 + k
		// n = i*n2 + j
		const int k = p % n3;
		const int n = (p-k)/n3;
		const int j = n % n2;
		const int i = n/n2;
		ijk[0] = i;
		ijk[1] = j;
		ijk[2] = k;
	}
}

extern void BlockWorkDivide(TGridIndex& a, TGridIndex& b);
/* Given an integer range a..b of iterations, distribute the iterations as evenly as possible
 * to Npes available processors, where the current process-ID is mype (0<=mype<Npes).
 * Npes,mype are external int variables which are initialized in gengrid.C:machine::init()
 * (or, they are #defined as 1 and 0, if USE_SHMEM is not defined.)
 * a and b will be modified as the new range for processor mype.
 * Calling this routine once in each PE is enough to reduce the global range a..b
 * to the local range.
 */


// --------------- Some bit manipulation functions -----------------

inline bool GetBit(smallnat d, unsigned ch)
// Return the d'th bit of ch
{
	return ((ch & (1 << d)) != 0);
}

inline unsigned InvertBit(smallnat d, unsigned ch)
// Return ch with d'th bit inverted
{
	return (ch ^ (1 << d));
}

inline unsigned SetBit(smallnat d, unsigned ch)
// Return ch with d'th bit set to 1
{
	return (ch | (1 << d));
}

inline unsigned ClearBit(smallnat d, unsigned ch)
// Return ch with d'th bit set to 0
{
	return (ch & ~(1 << d));
}

inline unsigned InsertBit(smallnat n, smallnat b, unsigned ch)
// Make space for n'th bit in bit field ch and set the bit to value b
{
	const unsigned highmask = (~0) << n;
	const unsigned lowmask = ~highmask;
	return ((ch & highmask) << 1) | (b << n) | (ch & lowmask) ;
}

inline unsigned NeighbourChild(smallnat d, smallnat dir, smallnat k)
{
	return InsertBit(d,!dir,k);
}

#define GRID_H

#endif
