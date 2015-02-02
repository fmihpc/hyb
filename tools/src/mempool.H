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

#ifndef MEMPOOL_H

#ifdef __GNUC__
#  pragma interface
#endif

#include "grid.H"
#include "fileheader.H"

#ifndef USE_CSHMAT
#  if USE_SHMEM
#    define USE_CSHMAT 1
#  else
#    define USE_CSHMAT 0
#  endif
#endif

#ifndef CACHE_NEIGHBOURS
#  define CACHE_NEIGHBOURS 1
#endif

// Actually se need only cshmat.H if USE_CSHMAT is 1, but we avoid
// conditional inclusion because it makes makdep.dat system-dependent
// and we are able to generate it only on gcc/g++ compiler!
#include "cshmat.H"
#include "sharr.H"

extern real NotANumber;		// filled by any Tmempool constructor

class Tmempool {
private:
	bool merge_load;
	bool abort_on_memfinito;
#if USE_CSHMAT
	TCachedSharedMatrix<real,TGridIndex> data;
#else
	TSharedArray<real> rdata;
	TSharedArray<TGridIndex> idata;
#endif
	TGridIndex maxnc;
	// These items must be handled by PEcoherency():
	// ---------------------------------------------
	TGridIndex heap1;
	TGridIndex freelist1;
	TGridIndex freepool;
	TGridIndex heap2;
	TGridIndex freelist2;
	TGridIndex Nfreelist1, Nfreelist2;
	// ---------------------------------------------
	int clen_r;				// number of reals in one clen-object
	int clen_i;				// number of TGridIndex'es in one clen-object
	int c1;					// number of clen-objects in heap1 blocks (in practice 2,4 or 8)
	int c2;					// number of clen-objects in heap2 blocks (in practice 0,1 or 2)
	bool dirty;				// cleared when alloc_base has been called
	TRealFormat loaded_realformat;
	void memfinito();
	void WriteCells(ostream *optr,
					const TIndexTable& newindex, TGridIndex Nnonremoved,
					const int savedreals[], int Nsavedreals,
					int parent_ptr_index, int child_ptr_index,
					TRealFormat realformat, int n_int_bytes, bool parallel_IO) const;
	TGridIndex list_length(TGridIndex list) const;
public:
	Tmempool() {dirty=true;}
	void init(TGridIndex maxnc1,
			 int clen_r1, int clen_i1,
			 int c1_1, int c2_1);
	Tmempool(TGridIndex maxnc1,
			 int clen_r1, int clen_i1,
			 int c1_1, int c2_1)
	{
		init(maxnc1,clen_r1,clen_i1,c1_1,c2_1);
	}
#if !USE_CSHMAT
	TSharedArray<real>& M() {return rdata;}
#endif
	void PEcoherency(int rootpe=ROOT_PE);
	TGridIndex get_heap1() const {
#		if PARALLEL_ADAPTATION
		return shmemget(&heap1,ROOT_PE);
#		else
		return heap1;
#		endif
	}
	void alloc_base(int n);
	void set_merge_load(bool flag=true) {merge_load=flag;}
	TGridIndex alloc_heap1();
	TGridIndex alloc_heap2();
	void dealloc_heap1(TGridIndex);
	void dealloc_heap2(TGridIndex);
	TGridIndex freecells() const {
#		if PARALLEL_ADAPTATION
		return shmemget(&heap2,ROOT_PE) - shmemget(&freepool,ROOT_PE);
#		else
		return heap2 - freepool;
#		endif
	}
	TGridIndex Ncells() const {
#		if PARALLEL_ADAPTATION
		return shmemget(&freepool,ROOT_PE);
#		else
		return freepool;
#		endif
	}
	TGridIndex Ncells_fast() const {		// You can call this if PEcoherency is known to persist
		return freepool;
	}
	TGridIndex allocated_mem_as_units_of_cell() const {
#		if PARALLEL_ADAPTATION
		return maxnc - shmemget(&heap2,ROOT_PE) + shmemget(&freepool,ROOT_PE)
			- c1*shmemget(&Nfreelist1,ROOT_PE) - c2*shmemget(&Nfreelist2,ROOT_PE);
#		else
		return (maxnc-heap2) + freepool - c1*Nfreelist1 - c2*Nfreelist2;
#		endif
	}
#	if USE_CSHMAT
	TCacheIndex load(TGridIndex i) const {return data.load(i);}
	void popcache() const {data.popcache();}
	void clearcache() const {data.clearcache();}
	real MC(TCacheIndex cell, smallnat comp) const {return data.M1C(cell,comp);}
	TGridIndex IMC(TCacheIndex cell, smallnat comp) const {return data.M2C(cell,comp);}
	real M(TGridIndex cell, smallnat comp) const {return data.M1(cell,comp);}
	real Mlocal(TGridIndex cell, smallnat comp) const {return data.M1local(cell,comp);}
	void Mset(TGridIndex cell, smallnat comp, real val) {data.M1set(cell,comp,val);}
	void Msetlocal(TGridIndex cell, smallnat comp, real val) {data.M1setlocal(cell,comp,val);}
	void Mset_multi(TGridIndex cell, smallnat comp0, smallnat ncomps, const real src[], int src_stride) {
		data.M1set_multi(cell,comp0,ncomps,src,src_stride);
	}
	TGridIndex IM(TGridIndex cell, smallnat comp) const {return data.M2(cell,comp);}
	TGridIndex IMlocal(TGridIndex cell, smallnat comp) const {return data.M2local(cell,comp);}
	void IMset(TGridIndex cell, smallnat comp, TGridIndex val) {data.M2set(cell,comp,val);}
	TGridIndex IMatomically_read_and_set(TGridIndex cell, smallnat comp, TGridIndex val) {
		return data.M2atomically_read_and_set(cell,comp,val);
	}
	void IMsetlocal(TGridIndex cell, smallnat comp, TGridIndex val) {data.M2setlocal(cell,comp,val);}
#	if USE_SHMEM
	void IMsetlock(TGridIndex cell, smallnat comp) {data.M2setlock(cell,comp);}
	void IMunsetlock(TGridIndex cell, smallnat comp) {data.M2unsetlock(cell,comp);}
#	endif
#	else
	void clearcache() const {}
#	ifdef TRANSPOSED_STORAGE
	TGridIndex flatMindex(TGridIndex cell, int comp) const {return comp*maxnc + cell;}
	TGridIndex flatIMindex(TGridIndex cell, int comp) const {return comp*maxnc + cell;}
	TGridIndex compstep() const {return maxnc;}		// flatindex distance between consecutive comp elements
	TGridIndex indexMstep() const {return 1;}		// flatindex distance between consecutive cell elements
	TGridIndex indexIMstep() const {return 1;}
#	else
	TGridIndex flatMindex(TGridIndex cell, int comp) const {return cell*clen_r + comp;}
	TGridIndex flatIMindex(TGridIndex cell, int comp) const {return cell*clen_i + comp;}
	TGridIndex compstep() const {return 1;}		// see above
	TGridIndex indexMstep() const {return clen_r;}
	TGridIndex indexIMstep() const {return clen_i;}
#	endif
	real M(TGridIndex cell, int comp) const {return rdata(flatMindex(cell,comp));}
	real Mlocal(TGridIndex cell, int comp) const {return rdata[ flatMindex(cell,comp) ];}
	void Mset(TGridIndex cell, int comp, real val) {rdata.put(flatMindex(cell,comp),val);}
	void Msetlocal(TGridIndex cell, int comp, real val) {rdata.lput(flatMindex(cell,comp),val);}
	real M(TGridIndex flatindex) const {return rdata(flatindex);}
	void Mset(TGridIndex flatindex, real val) {rdata.put(flatindex,val);}
	TGridIndex IM(TGridIndex cell, int comp) const {return idata(flatIMindex(cell,comp));}
	TGridIndex IMlocal(TGridIndex cell, int comp) const {return idata[ flatMindex(cell,comp) ];}
	void IMset(TGridIndex cell, int comp, TGridIndex val) {idata.put(flatIMindex(cell,comp),val);}
	void IMsetlocal(TGridIndex cell, int comp, TGridIndex val) {idata.lput(flatIMindex(cell,comp),val);}
	TGridIndex IMatomically_read_and_set(TGridIndex cell, smallnat comp, TGridIndex val) {
		return idata.atomically_read_and_set(flatIMindex(cell,comp),val);
	}
#	endif
	bool streamsave(ostream& o,
					const TIndexTable& newindex, TGridIndex Nnonremoved,
					TRealFormat realformat,
					int parent_ptr_index, int child_ptr_index,
					const int savedreals[], int Nsavedreals,
					const int savedints_leaf[], int Nsavedints_leaf,
					const int savedints_nonleaf[], int Nsavedints_nonleaf,
					const char *filename_base=0
					) const;
	bool streamload(istream& o, const Theader& h,
					int parent_ptr_index, int child_ptr_index);
	TRealFormat LoadedRealformat() const {return loaded_realformat;}
	friend ostream& operator<<(ostream& o, const Tmempool& mp);
	~Tmempool();
};

extern ostream& operator<<(ostream& o, const Tmempool& mp);

#define MEMPOOL_H

#endif
