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
#  pragma implementation "mempool.H"
#endif

#include <cstring>
#include <cmath>
#include <fstream>
#include <cstdio>
#if defined(_UNICOS) && !defined(_CRAYIEEE)
#  include <float.h>
#endif
using namespace std;
#include "mempool.H"
#include "byteconv.H"

extern void doabort();
real NotANumber;		// filled by the first Tmempool constructor

void Tmempool::init(TGridIndex maxnc1,
					int clen_r1, int clen_i1,
					int c1_1, int c2_1)
{
	abort_on_memfinito = false;
	dirty = true;
	maxnc = maxnc1;
	clen_r = clen_r1;
	clen_i = clen_i1;
	c1 = c1_1;
	c2 = c2_1;
#if USE_CSHMAT
	data.init(maxnc,clen_r,clen_i);
#else
	rdata.init(maxnc*clen_r,clen_r);
	idata.init(maxnc*clen_i,clen_i);
#endif
	heap1 = 0;
	heap2 = maxnc;
	freepool = 0;
	freelist1 = freelist2 = NOINDEX;
	Nfreelist1 = Nfreelist2 = 0;
	loaded_realformat = REALFORMAT_ASCII;
	/* Fill the NotANumber external variable.
	 * Filling a float or double with -1 (0xFF) bytes gives a NaN in all
	 * known arithmetics: CRAY PVP, little-endian IEEE, big-endian IEEE.
	 * (Actually, in CRAY PVP it gives floating point exception if any FP-operation
	 * is attempted on NotANumber, since Cray does not have NaNs.)
	 */
	memset(&NotANumber,0xFF,sizeof(real));
	shmembarrierall();
}

void Tmempool::PEcoherency(int rootpe)
{
	if (shmemmype() != rootpe) {
		heap1 = shmemget(&heap1,rootpe);
		freelist1 = shmemget(&freelist1,rootpe);
		freepool = shmemget(&freepool,rootpe);
		heap2 = shmemget(&heap2,rootpe);
		freelist2 = shmemget(&freelist2,rootpe);
		Nfreelist1 = shmemget(&Nfreelist1,rootpe);
		Nfreelist2 = shmemget(&Nfreelist2,rootpe);
	}
}

void Tmempool::alloc_base(int n)
{
	if (dirty) {
		heap1 = freepool = n;
		if (heap1 > maxnc) {
			cerr << "*** Tmempool::alloc_base(" << n << "): out of memory, max n is " << maxnc << "\n";
			doabort();
		}
		dirty = false;
	} else {
		cerr << "*** Tmempool::alloc_base called second time\n";
		doabort();
	}
	shmembarrierall();	// Make sure that ROOT_PE has completed alloc_base before others try to access heap1, freepool using SHMEM
}

void Tmempool::memfinito()
{
	cerr << "*** Tmempool: out of memory\n";
	doabort();
}

#if PARALLEL_ADAPTATION
static long mempool_lock = 0;
#define enter_critical() shmemsetlock(&mempool_lock)
#define exit_critical() shmemclearlock(&mempool_lock)
#else
#define enter_critical()
#define exit_critical()
#endif

TGridIndex Tmempool::alloc_heap1()
// Allocate c1 heap items from heap1.
// If there is no more memory, abort with an error message if abort_on_memfinito is true,
// otherwise return NOINDEX.
{
	assert(!dirty);
	TGridIndex result;
	enter_critical();
#		if PARALLEL_ADAPTATION
		result = shmemget(&freepool,ROOT_PE);
		shmemput(&freepool, result + c1, ROOT_PE);
		if (result + c1 > shmemget(&heap2, ROOT_PE)) {
#		else
		result = freepool;
		freepool+= c1;
		if (freepool > heap2) {
#		endif
			if (abort_on_memfinito)
				memfinito();
			else {
				warn("Tmempool::alloc_heap1: out of memory");
#				if PARALLEL_ADAPTATION
				shmemput(&freepool,result,ROOT_PE);
#				else
				freepool = result;
#				endif
				result = NOINDEX;
			}
		}
	exit_critical();
	return result;
}

void Tmempool::dealloc_heap1(TGridIndex i)
{
	assert(!dirty);
#	if !PARALLEL_ADAPTATION
	assert(i>=heap1 && i<freepool && (i-heap1)%c1 == 0);
#	endif
	enter_critical();
#	if PARALLEL_ADAPTATION
	IMset(i,0, shmemget(&freelist1,ROOT_PE));
	shmemput(&freelist1,i,ROOT_PE);
	shmemput(&Nfreelist1, shmemget(&Nfreelist1,ROOT_PE) + 1, ROOT_PE);
#	else
	IMset(i,0, freelist1);
	freelist1 = i;
	Nfreelist1++;
#	endif
	exit_critical();
}

TGridIndex Tmempool::alloc_heap2()
{
	assert(!dirty);
	TGridIndex result;
	enter_critical();
#		if USE_SHMEM
		TGridIndex hp2 = shmemget(&heap2,ROOT_PE);
		shmemput(&heap2, hp2 - c2, ROOT_PE);
		if (shmemget(&freepool,ROOT_PE) > hp2-c2) {
#		else
		heap2-= c2;
		if (freepool > heap2) {
#		endif
				memfinito();
		}
#		if PARALLEL_ADAPTATION && USE_SHMEM
		result = hp2 - c2;
#		else
		result = heap2;
#		endif
	exit_critical();
	return result;
}

void Tmempool::dealloc_heap2(TGridIndex i)
{
	assert(!dirty);
	assert(i < maxnc);
#	if !PARALLEL_ADAPTATION
	assert(i >= heap2);
	assert((i-heap2) % c2 == 0);
#	endif
	enter_critical();
#	if PARALLEL_ADAPTATION
	TGridIndex fl2 = shmemget(&freelist2,ROOT_PE);
	IMset(i,0, fl2);
	shmemput(&freelist2, i, ROOT_PE);
	shmemput(&Nfreelist2, shmemget(&Nfreelist2,ROOT_PE) + 1, ROOT_PE);
#	else
	IMset(i,0, freelist2);
	freelist2 = i;
	Nfreelist2++;
#	endif
	exit_critical();
}

ostream& operator<<(ostream& o, const Tmempool& mp)
{
	// Write 5 numbers: heap1_free, heap1_alloc, freepool, heap2_alloc, heap2_free
	o << mp.Nfreelist1;
	o << ' ' << mp.freepool - mp.heap1 - mp.Nfreelist1;
	o << ' ' << mp.heap2 - mp.freepool;
	o << ' ' << mp.maxnc - mp.heap2 - mp.Nfreelist2;
	o << ' ' << mp.Nfreelist2;
	o << '\n';
	return o;
}

static int NumberOfBytes(TGridIndex x)
// Take care that e.g. NumberOfBytes(50000) is 3, not 2, because the sign bit must
// also be included. That's why we compute (nbits+8)/8, otherwise it would be
// (nbits+7)/8.
{
	if (x < 0) x = -x;
	const int maxnbits = 8*sizeof(x);
	int i,nbits=0;
	for (i=0; i<maxnbits; i++)
		if (((x >> i) & 0x1) != 0) nbits = i+1;
	return ((nbits+8)/8);
}

inline void WriteInt(ostream& o, const TIndexTable& newindex,
					 TGridIndex x,
					 int n_int_bytes, TRealFormat realformat)
{
	TGridIndex xx = (x == NOINDEX || newindex.isempty()) ? x : newindex(x);
#ifdef DEBUG
	if (!newindex.isempty() && x != NOINDEX) assert(xx != NOINDEX);
#endif
	if (xx < 0) xx = -1;		// the external representation of NOINDEX is always -1
	if (realformat == REALFORMAT_ASCII)
		o << ' ' << xx;
	else {
		int i;
		for (i=0; i<n_int_bytes; i++)
			o.put((xx >> 8*i) & 0xFF);
	}
}

inline void ReadInt(istream& o, TGridIndex& x, int n_int_bytes, TRealFormat realformat)
{
	if (realformat == REALFORMAT_ASCII)
		o >> x;
	else {
		// We have to be careful: TGridIndex is a signed integer.
		// If the one stored in the file has fewer bytes, we have to
		// explicitly sign-extend, in case it is negative.
		// Use a temporary unsigned long during the shifts.
		// --------------------------------------------------------
		// NOTE: This works only if two's complement representation
		// is in use.
		// --------------------------------------------------------
		assert(sizeof(TGridIndex) <= sizeof(long));
		int i;
		unsigned long ux = 0;
		for (i=0; i<n_int_bytes; i++)
			ux|= ((unsigned long)o.get() << 8*i);
		if (ux & (1 << (8*n_int_bytes-1))) {
			// The number is negative. Form a OR-mask which contains
			// zeros on those bit positions we already have read in.
			const long mask = (-1) << (8*n_int_bytes);
			ux|= (unsigned long)mask;
		}
		x = TGridIndex(ux);
	}
	if (x < 0) x = NOINDEX;		// The external representation of NOINDEX is always -1, but still NOINDEX may be != -1
}


inline void WriteDoubles(ostream& o, double x[], int n, TRealFormat realformat)
{
	int i;
	float *xf = (float *)x;
	switch (realformat) {
	case REALFORMAT_ASCII:
		for (i=0; i<n; i++)
			o << ' ' << x[i];
		break;
	case REALFORMAT_FLOAT:
		if (sizeof(float) != sizeof(double))
			for (i=0; i<n; i++) xf[i] = float(x[i]);
		ByteConversion(sizeof(float),(unsigned char*)xf,n);
		WriteFloatsToFile(o,xf,n);
		break;
	case REALFORMAT_DOUBLE:
		ByteConversion(sizeof(double),(unsigned char*)x,n);
		WriteDoublesToFile(o,x,n);
		break;
	case REALFORMAT_INQUIRE:
		break;
	}
}

inline void ReadDoubles(istream& o, double x[], int n, TRealFormat realformat)
{
	int i;
	float *xf = (float*)x;
	switch (realformat) {
	case REALFORMAT_ASCII:
		for (i=0; i<n; i++)
			o >> x[i];
		break;
	case REALFORMAT_FLOAT:
		if (sizeof(float) != sizeof(double)) {
			float *xf = new float [n];
			o.read((char*)xf,sizeof(float)*n);
			ByteConversion_input(sizeof(float),(unsigned char*)xf,n);
			for (i=0; i<n; i++) x[i] = xf[i];
			delete [] xf;
		} else {
			ReadFloatsFromFile(o,xf,n);
			ByteConversion_input(sizeof(double),(unsigned char*)xf,n);
		}
		break;
	case REALFORMAT_DOUBLE:
		ReadDoublesFromFile(o,x,n);
		ByteConversion_input(sizeof(double),(unsigned char*)x,n);
		break;
	case REALFORMAT_INQUIRE:
		break;
	}
}

inline void Tmempool::WriteCells(ostream *optr,
								 const TIndexTable& newindex, TGridIndex Nnonremoved,
								 const int savedreals[], int Nsavedreals,
								 int parent_ptr_index,
								 int child_ptr_index,
								 TRealFormat realformat,
								 int n_int_bytes,
								 bool parallel_IO) const
{
	TGridIndex i;
	smallnat a;
	double *dbuff = new double [Nsavedreals];
	char *FFFF = new char [Nsavedreals*sizeof(double)];
	memset(FFFF,0xFF,Nsavedreals*sizeof(double));
	TGridIndex Nsaved = 0;
	TGridIndex i1=0,i2=freepool-1;
	if (parallel_IO) {
		BlockWorkDivide(i1,i2);
	}
	for (i=i1; i<=i2; i++) {
		if (!newindex.isempty() && newindex(i) == NOINDEX) continue;
		// Process cell i
		// If it's leaf cell ('L'), zero leaf cell ('Z') or non-leaf ('N')
		Nsaved++;
		const bool isleaf = (IM(i,child_ptr_index) == NOINDEX);
		if (isleaf) {
			pragma("_CRI ivdep"); pragma("_CRI shortloop");
			for (a=0; a<Nsavedreals; a++)
				dbuff[a] = M(i,savedreals[a]);
			// If NaN's, zero them before writing to file. Eases postprocessing.
			// Dead cells always contains NaNs, and will again when read from file.
			if (!memcmp(dbuff,FFFF,sizeof(double)*Nsavedreals)) {
				pragma("_CRI ivdep"); pragma("_CRI shortloop");
				for (a=0; a<Nsavedreals; a++) dbuff[a] = 0;
			}
			bool allzeros = true;
			pragma("_CRI ivdep"); pragma("_CRI shortloop");
			for (a=0; a<Nsavedreals; a++) if (dbuff[a] != 0) allzeros = false;
			optr->put(allzeros ? 'Z' : 'L');
			WriteInt(*optr,newindex,IM(i,parent_ptr_index),n_int_bytes,realformat);
			if (!allzeros) WriteDoubles(*optr,dbuff,Nsavedreals,realformat);
		} else {
			// Is not leaf cell, write pointer to first child
			optr->put('N');
			WriteInt(*optr,newindex,IM(i,parent_ptr_index),n_int_bytes,realformat);
			WriteInt(*optr,newindex,IM(i,child_ptr_index),n_int_bytes,realformat);
		}
		if (realformat == REALFORMAT_ASCII) *optr << '\n';
	}
	if (Nsaved != Nnonremoved && Npes==1)
		cerr << "Tmempool::WriteCells: Nsaved=" << Nsaved << ", Nnonremoved=" << Nnonremoved << "\n";
	delete [] FFFF;
	delete [] dbuff;
}

bool Tmempool::streamsave(ostream& o,
						  const TIndexTable& newindex, TGridIndex Nnonremoved,
						  TRealFormat realformat,
						  int parent_ptr_index, int child_ptr_index,
						  const int savedreals[], int Nsavedreals,
						  const int savedints_leaf[], int Nsavedints_leaf,
						  const int savedints_nonleaf[], int Nsavedints_nonleaf,
						  const char *filename_base
						  ) const
{
	const bool parallel_IO = (filename_base != 0);
	assert(!dirty);
	if (Nsavedints_leaf != 1 || savedints_leaf[0] != parent_ptr_index) {
		cerr << "*** Tmempool::streamsave: Nsavedints_leaf=" << Nsavedints_leaf
			 << ", savedints_leaf[0]=" << savedints_leaf[0] << "\n";
		doabort();
	}
	if (Nsavedints_nonleaf != 2 || savedints_nonleaf[0] != parent_ptr_index || savedints_nonleaf[1] != child_ptr_index) {
		cerr << "*** Tmempool::streamsave: Nsavedints_nonleaf=" << Nsavedints_nonleaf
			 << ", savedints_nonleaf[0]=" << savedints_nonleaf[0] << "\n";
		doabort();
	}
#	if defined(_CRAY1) && defined(_CRAYIEEE)
	// T90 cannot produce 32-bit IEEE floats
	if (realformat == REALFORMAT_FLOAT) {
		warn("IEEE Cray PVP system does not support REALFORMAT_FLOAT");
		realformat = REALFORMAT_DOUBLE;
	}
#	endif
	if (!newindex.isempty() && freepool - Nnonremoved != Nfreelist1*c1 && mype==0) {
		cerr << "Tmempool::streamsave warning: freepool=" << freepool << ", Nnonremoved=" << Nnonremoved
			 << ", Nfreelist1=" << Nfreelist1 << "\n";
		cerr << "  (" << freepool - Nnonremoved - Nfreelist1*c1 << " too short freelist1)\n";
	}
	const int n_int_bytes = NumberOfBytes(maxnc);
	if (mype==ROOT_PE) {
		o << "ncells = " << Nnonremoved /*freepool*/ << "\n";
		o << "nablocks = 0\n";
		o << "bytes_per_int = " << n_int_bytes << "\n";
		o << "freelist1 = " << (!newindex.isempty() ? NOINDEX : freelist1) << "\n";
		// If newindex is in use, freelist1 will be saved as NOINDEX, because we do not save removed cells
		o << "freelist2 = " << NOINDEX << "\n";
		// Freelist2 will be saved as NOINDEX, because we do not save a-blocks at all
		o << "eoh\n";
	}
	ostream *optr = &o;
	ofstream o2;
	if (Npes > 1 && parallel_IO && filename_base != 0) {
		char *fn = new char [strlen(filename_base) + 20];
		sprintf(fn,"%s.%.4d",filename_base,mype);
		o2.open(fn);
		o2.precision(o.precision());
		optr = &o2;
		delete [] fn;
	} else {
		if (mype != ROOT_PE) return true;
	}
	switch (realformat) {
	case REALFORMAT_ASCII:
		WriteCells(optr,newindex,Nnonremoved,savedreals,Nsavedreals,
				   parent_ptr_index,child_ptr_index,REALFORMAT_ASCII,n_int_bytes, parallel_IO);
		break;
	case REALFORMAT_DOUBLE:
		WriteCells(optr,newindex,Nnonremoved,savedreals,Nsavedreals,
				   parent_ptr_index,child_ptr_index,REALFORMAT_DOUBLE,n_int_bytes, parallel_IO);
		break;
	case REALFORMAT_FLOAT:
		WriteCells(optr,newindex,Nnonremoved,savedreals,Nsavedreals,
				   parent_ptr_index,child_ptr_index,REALFORMAT_FLOAT,n_int_bytes, parallel_IO);
		break;
	case REALFORMAT_INQUIRE:
		break;
	}
	return (mype==ROOT_PE) ? o.good() : true;
}

TGridIndex Tmempool::list_length(TGridIndex list) const
{
	TGridIndex result = 0, p = list;
	while (p >= 0) {result++; p=IM(p,0);}
	return result;
}

bool Tmempool::streamload(istream& o, const Theader& h,
						  int parent_ptr_index, int child_ptr_index)
{
	if (strcmp(h.getstr("type"),"hc")) {
		cerr << "*** Tmempool::streamload: is not HC grid file\n";
		return false;
	}
	TGridIndex maxnc_file = TGridIndex(h.getint("maxnc"));
	freepool = TGridIndex(h.getint("ncells"));
	if (merge_load) {
		if (maxnc_file > maxnc) {
			if (freepool > maxnc) {
				cerr << "*** Tmempool::streamload: File ncells=" << freepool << " too large - load cancelled\n";
				return false;
			} else {
				warn("File maxnc is larger, memory shortage possible later");
				maxnc_file = maxnc;
			}
		}
	} else
		maxnc = maxnc_file;
	heap1 = TGridIndex(h.getint("n1"))*TGridIndex(h.getint("n2"))*TGridIndex(h.getint("n3"));
	heap2 = maxnc;
	freelist1 = TGridIndex(h.getint("freelist1"));
	freelist2 = TGridIndex(h.getint("freelist2"));
	assert(freelist2 == NOINDEX);
	Nfreelist1 = list_length(freelist1);
	Nfreelist2 = list_length(freelist2);
	assert(Nfreelist2 == 0);
	const int n_int_bytes = int(h.getint("bytes_per_int"));
	char *realformat_str = h.getstr("realformat");
	TRealFormat realformat=REALFORMAT_DOUBLE;
	if (!strcmp(realformat_str,"float"))
		realformat = REALFORMAT_FLOAT;
	else if (!strcmp(realformat_str,"double"))
		realformat = REALFORMAT_DOUBLE;
	else if (!strcmp(realformat_str,"ascii"))
		realformat = REALFORMAT_ASCII;
	loaded_realformat = realformat;
#	if defined(_CRAY1) && defined(_CRAYIEEE)
	if (realformat == REALFORMAT_FLOAT) {
		cerr << "*** IEEE Cray PVP system cannot read realformat=float data!\n";
		doabort();
	}
#	endif
	const int ncd = int(h.getint("ncd"));
	const int nsd = int(h.getint("nsd"));
	const int dim = int(h.getint("dim"));
	if (!merge_load) {

		clen_r = ncd + dim*nsd;
		/*
		 * WARNING : Actually we should write clen_i = ablockpointer_offset() + dim,
		 * but ablockpointer_offset is defined inside class THCgrid and unaccessible here.
		 * Moving the definition outside THCgrid would be ugly, too, since then the other
		 * offset inliners should also be moved, otherwise they are spread in two places.
		 * PLEASE REMEMBER: If you have to add more integer fields to THCgrid so that
		 * ablockpointer_offset() ever changes, you MUST make the corresponding changes also here.
		 */
		clen_i =
#if PARALLEL_ADAPTATION
			6
#else
			5
#endif
			+ dim;
		c1 = (1 << dim);
		c2 = ((1 << (dim-1))*nsd > ncd+dim*nsd) ? 2 : 1;
	}
	const int Nsavedreals = ncd + dim*nsd;
	TGridIndex i;
	smallnat a;
	char buff[1001];	// used only if realformat=REALFORMAT_ASCII
	double *dbuff = new double [Nsavedreals];
	for (i=0; i<freepool; i++) {
		// Process cell i
		// If it's leaf cell ('L'), zero leaf cell ('Z'), or non-leaf cell ('N') ...
		const char LZN = o.get();
		// Set NOINDEX to every int field first.
		for (a=0; a<clen_i; a++) IMset(i,a,NOINDEX);
		if (LZN == 'L' || LZN == 'Z') {
			// Is leaf cell, read pointer to parent (int), and savedreals
			TGridIndex tmp;
			ReadInt(o,tmp,n_int_bytes,realformat);
			IMset(i,parent_ptr_index, tmp);
			// Invalidate the data first. Do it for consistency. These values should never be accessed.
			// Note that we use clen_r, not Nsavedreals, these may be different if merge_load is true.
			for (a=0; a<clen_r; a++) Mset(i,a, NotANumber);
			if (LZN == 'L') {
				ReadDoubles(o,dbuff,Nsavedreals,realformat);
				pragma("_CRI ivdep"); pragma("_CRI shortloop");
				for (a=0; a<Nsavedreals; a++)
					Mset(i,a /*savedreals[a]*/, dbuff[a]);
			}
			// Set child pointer to null
			IMset(i,child_ptr_index, NOINDEX);
		} else if (LZN == 'N') {
			// Is not leaf cell, read pointers to parent and first child
			TGridIndex tmp;
			ReadInt(o,tmp,n_int_bytes,realformat);
			IMset(i,parent_ptr_index,tmp);
			ReadInt(o,tmp,n_int_bytes,realformat);
			IMset(i,child_ptr_index,tmp);
			// Invalidate the data. Do it for consistency. These values should never be accessed.
			// Note that we use clen_r, not Nsavedreals, these may be different if merge_load is true.
			for (a=0; a<clen_r; a++) Mset(i,a, NotANumber);
		} else {
			cerr << "*** Syntax error in " << (realformat==REALFORMAT_ASCII ? "data" : "binary") << " part of HC file:\n";
			cerr << "    Flagbyte is neither 'L', 'Z' nor 'N' but '" << LZN << "'.\n";
			return false;
		}
		if (realformat == REALFORMAT_ASCII) o.getline(buff,1000);	// pass '\n'
	}
	delete [] dbuff;
	return o.good();
}

Tmempool::~Tmempool()
{
//	delete [] idata;
//	delete [] rdata.RawPtr();
}
