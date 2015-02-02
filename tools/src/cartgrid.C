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
#  pragma implementation "cartgrid.H"
#endif

#include "cartgrid.H"
#include "fileheader.H"
#include "byteconv.H"
#include <stdlib.h>
#include <math.h>
#include <string.h>

template <smallnat dim>
void TCartesianGrid<dim>::init(smallnat ncd1, smallnat nsd1, smallnat nfq1, TGridIndex maxnc1, bool use_mapping1)
{
	init(ncd1,nsd1,nfq1,maxnc1,use_mapping1);
	n1=n2=n3=ntot=n123[0]=n123[1]=n123[2]=1; dx=0.01;
	cornercache.set_dx(dx);
	data.init(maxnc*clen,clen);
	invvol.init(maxnc);
	if (use_mapping) physsize.init(maxnc);

	cellinfotab.init(maxnc);
	TGridIndex i;
	pragma("_CRI ivdep");
	for (i=0; i<maxnc; i++)
		cellinfotab.put(i, 0);		// 0 is a suitable default cell type, see grid.H
}

template <smallnat dim> void TCartesianGrid<dim>::PEcoherency(int) {}

template <smallnat dim>
void TCartesianGrid<dim>::setup_after_regular()
{
	const real area = pow(dx,real(dim-1));
	const real invV = 1.0/pow(dx,real(dim));
	TGridIndex i;
	smallnat d;
	for (i=0; i<ntot; i++) {
		setinvvolume(i,invV);
		if (use_mapping) setphyssize(i,dx);
		for (d=0; d<dim; d++) Mset(i,areadataoffset(d),area);
	}
}

template <smallnat dim>
void TCartesianGrid<dim>::vnext1(TGridIndexVector& iv, TGridIndex start, TGridIndex stride, Ttraversal trav) const
{
	TGridIndex i;
	smallnat j;
#ifdef _CRAY1
	const TGridIndex buflen = 128;
#else
	const TGridIndex buflen = 13;
#endif
	if (trav != TRAV_ALL && trav != TRAV_LEAF_ONLY) {
		int nresult = 0, k, n;
		i = start;
		TGridIndex buf[buflen];
		while (nresult < VECLEN && i<ntot) {
			n = min(buflen,ntot-i);
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			for (j=0,k=0; j<n; i+=stride,j++)
				if (isindomain(i))
					buf[k++] = i;
			n = min(k,VECLEN-nresult);
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			for (j=0; j<n; nresult++,j++) iv[nresult] = buf[j];
		}
		iv.setlength(nresult);
	} else if (trav == TRAV_LEAFGHOST_ONLY) {
		TGridIndex nresult = 0, k, n;
		i = start;
		TGridIndex buf[buflen];
		while (nresult < VECLEN && i<ntot) {
			n = min(buflen,ntot-i);
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			for (j=0,k=0; j<n; i+=stride,j++)
				if (isleaf(i) && celltype(i) == GHOST_CELL)
					buf[k++] = i;
			n = min(k,TGridIndex(VECLEN)-nresult);
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			for (j=0; j<n; nresult++,j++) iv[nresult] = buf[j];
		}
		iv.setlength(nresult);
	} else {
		const TGridIndex n = min(TGridIndex(VECLEN),ntot-start);
		const TGridIndex stop = start + n;
		VDIRS;
		for (i=start,j=0; i<stop; i+=stride,j++) iv[j] = i;
		iv.setlength(n);
	}
}

#ifdef _CRAY1
#  define INLINE inline
#else
#  define INLINE
#endif

template <smallnat dim>
INLINE void TCartesianGrid<dim>::vdecompose(const TGridIndexVector& iv, int ijk[3][VECLEN]) const
{
	smallnat v;
	const int n = iv.length();
	if (dim == 1) {
		VLOOPN(iv,n) ijk[0][v] = iv(v);
	} else if (dim == 2) {
		// p = i*n2 + j
		VLOOPN(iv,n) {
			const int p = iv(v);
			const int j = p % n2;
			const int i = (p-j)/n2;
			ijk[0][v] = i;
			ijk[1][v] = j;
		}
	} else {
		// 3D: p = (i*n2 + j)*n3 + k
		// n = i*n2 + j
		VLOOPN(iv,n) {
			const int p = iv(v);
			const int k = p % n3;
			const int n = (p-k)/n3;
			const int j = n % n2;
			const int i = n/n2;
			ijk[0][v] = i;
			ijk[1][v] = j;
			ijk[2][v] = k;
		}
	}
}

template <smallnat dim>
INLINE void TCartesianGrid<dim>::vcompose(TGridIndexVector& iv, should_be_const int ijk[3][VECLEN], int vlen) const
{
	smallnat v;
	if (dim == 1) {
		VDIRS;
		for (v=0; v<vlen; v++) iv[v] = ijk[0][v];
	} else if (dim == 2) {
		VDIRS;
		for (v=0; v<vlen; v++) iv[v] = ijk[0][v]*n2 + ijk[1][v];
	} else {
		VDIRS;
		for (v=0; v<vlen; v++) iv[v] = (ijk[0][v]*n2 + ijk[1][v])*n3 + ijk[2][v];
	}
	iv.setlength(vlen);
}

template <smallnat dim>
TGridIndex TCartesianGrid<dim>::dogetneighbour(TGridIndex p, smallnat d, smallnat dir) const
{
	int ijk[3];
	decompose(p, ijk);
	if (dir) {
		ijk[d]++;
		if (ijk[d] >= n123[d]) {
			/*if (probeflag)*/ return NOINDEX;
		}
	} else {
		ijk[d]--;
		if (ijk[d] < 0) {
			/*if (probeflag)*/ return NOINDEX;
		}
	}
	return compose(ijk);
}

void vgetneighbour_error(smallnat dim, int ijk[3][VECLEN], TGridIndex error_v, smallnat d)
{
	// reverse the change to report correct error message
	if (ijk[d][error_v] < 0)
		ijk[d][error_v]++;
	else
		ijk[d][error_v]--;
	cerr << "*** TCartesianGrid::vgetneighbour_s:\n";
	cerr << "    (";
	smallnat d1;
	for (d1=0; d1<dim; d1++) {
		cerr << ijk[d1][error_v];
		if (d1<dim-1) cerr << ",";
	}
	cerr << ") has no left/right neighbour in dim " << d << "\n";
	doabort();
}

template <smallnat dim>
void TCartesianGrid<dim>::vdogetneighbour(
	TGridIndexVector& result,
	const TGridIndexVector& iv,
	const smallnat d[VECLEN], const smallnat dir[VECLEN]) const
{
	smallnat v;
	int ijk[3][VECLEN];
	vdecompose(iv, ijk);
	const smallnat n = iv.length();
	const bool probeflag = true;
	if (probeflag) {
		int errormask[VECLEN];		// Cray cannot vectorized if bool is used (bool is same as char for Cray)
		VLOOPN(iv,n) {
			const smallnat dv = d[v];
			if (dir[v]) {
				ijk[dv][v]++;
				errormask[v] = (ijk[dv][v] >= n123[dv]);
			} else {
				ijk[dv][v]--;
				errormask[v] = (ijk[dv] < 0);
			}
		}
		vcompose(result,ijk,iv.length());
		VLOOPN(iv,n) if (errormask[v]) result[v] = NOINDEX;
	} else {
		int error_v = NOINDEX;
		VLOOPN(iv,n) {
			const smallnat dv = d[v];
			if (dir[v]) {
				ijk[dv][v]++;
				if (ijk[dv][v] >= n123[dv]) error_v = v;
			} else {
				ijk[dv][v]--;
				if (ijk[dv] < 0) error_v = v;
			}
		}
		if (error_v != NOINDEX) vgetneighbour_error(dim,ijk,error_v,d[error_v]);
		vcompose(result,ijk,iv.length());
	}
}

template <smallnat dim>
INLINE void TCartesianGrid<dim>::vgetneighbour_dir0(
	TGridIndexVector& result,
	const TGridIndexVector& iv,
	const smallnat d[VECLEN]) const
// Special (optimized) case of vgetneighbour for dir==0
{
	smallnat v;
	int ijk[3][VECLEN];
	vdecompose(iv, ijk);
	const smallnat n = iv.length();
	const bool probeflag = true;
	if (probeflag) {
		int errormask[VECLEN];
		VLOOPN(iv,n) {
			const smallnat dv = d[v];
			ijk[dv][v]--;
			errormask[v] = (ijk[dv] < 0);
		}
		vcompose(result,ijk,iv.length());
		VLOOPN(iv,n) if (errormask[v]) result[v] = NOINDEX;
	} else {
		int error_v = NOINDEX;
		VLOOPN(iv,n) {
			const smallnat dv = d[v];
			ijk[dv][v]--;
			if (ijk[dv] < 0) error_v = v;
		}
		if (error_v != NOINDEX) vgetneighbour_error(dim,ijk,error_v,d[error_v]);
		vcompose(result,ijk,iv.length());
	}
}

template <smallnat dim>
void TCartesianGrid<dim>::vdogetneighbour_s(
	TGridIndexVector& result,
	const TGridIndexVector& iv,
	smallnat d, smallnat dir) const
{
	smallnat v;
	int ijk[3][VECLEN];
	vdecompose(iv, ijk);
	const smallnat n = iv.length();
	const bool probeflag = true;
	if (probeflag) {
		int errormask[VECLEN];
		if (dir) {
			VLOOPN(iv,n) {
				ijk[d][v]++;
				errormask[v] = (ijk[d][v] >= n123[d]);
			}
		} else {
			VLOOPN(iv,n) {
				ijk[d][v]--;
				errormask[v] = (ijk[d] < 0);
			}
		}
		vcompose(result,ijk,iv.length());
		VLOOPN(iv,n) if (errormask[v]) result[v] = NOINDEX;
	} else {
		TGridIndex error_v = NOINDEX;
		if (dir) {
			VLOOPN(iv,n) {
				ijk[d][v]++;
				if (ijk[d][v] >= n123[d]) error_v = v;
			}
		} else {
			VLOOPN(iv,n) {
				ijk[d][v]--;
				if (ijk[d] < 0) error_v = v;
			}
		}
		if (error_v != NOINDEX) vgetneighbour_error(dim,ijk,error_v,d);
		vcompose(result,ijk,iv.length());
	}
}

template <smallnat dim>
void TCartesianGrid<dim>::putcell(TGridIndex i, smallnat v)
{
	smallnat c;
	pragma("_CRI ivdep");
	for (c=0; c<ncd; c++) Mset(i,c,CT(c,v));
	assert(CT(0,v) > 0);
}

// Faster version of vputcell than the trivial one in grid.C
template <smallnat dim>
void TCartesianGrid<dim>::vputcell(const TGridIndexVector& iv)
{
	smallnat c,v;
#	ifdef VECTOR_MACHINE
	for (c=0; c<ncd; c++) {
		VLOOP(iv) {
			Mset(iv(v),c,CT(c,v));
			assert(CT(0,v) > 0);
			assert(CT(3,v) == 0);
		}
	}
#	else
	VLOOP(iv) {
		const TGridIndex i = iv(v);
		for (c=0; c<ncd; c++) Mset(i,c,CT(c,v));
	}
#	endif
}

template <smallnat dim>
void TCartesianGrid<dim>::getcell(TGridIndex i, smallnat min_d, smallnat max_d, smallnat v)
{
	smallnat c,s,d;
	TGridIndex j;
	pragma("_CRI ivdep");
	for (c=0; c<ncd; c++) CT(c,v) = M(i,c);
	for (d=min_d; d<=max_d; d++) {
		const int offset = surfdataoffset(d);
		const int aoffset = areadataoffset(d);
		pragma("_CRI ivdep");
		for (s=0; s<nsd; s++) FT(d,1,s,v) = M(i,offset+s);
		AT(d,1,v) = M(i,aoffset);
		assert(AT(d,1,v) > 0);
		j = getneighbour(i,d,0);
		pragma("_CRI ivdep");
		for (s=0; s<nsd; s++) FT(d,0,s,v) = M(j,offset+s);
		AT(d,0,v) = M(j,aoffset);
		assert(AT(d,0,v) > 0);
	}
}

template <smallnat dim>
void TCartesianGrid<dim>::fgetcell(TGridIndex i, smallnat v)
{
	smallnat c;
	pragma("_CRI ivdep");
	for (c=0; c<ncd; c++) CT(c,v) = M(i,c);
}

template <smallnat dim>
void TCartesianGrid<dim>::getnormal(TGridIndex i, TSurfSpec s, real n3[3])
{
	smallnat d1;
	const smallnat d = s.d();
	const smallnat dir = s.dir();
	TGridIndex j;
	for (d1=0; d1<3; d1++) n3[d1] = 0;
	if (dir == 1) {
	   for (d1=0; d1<dim; d1++) n3[d1] = M(i,this->normaldataoffset(s.d())+d1);
	} else {
		j = getneighbour(i,d,0);
	   for (d1=0; d1<dim; d1++) n3[d1] = M(j,this->normaldataoffset(s.d())+d1);
	}
}

// Faster version of vgetcell than the trivial one in grid.C
template <smallnat dim>
void TCartesianGrid<dim>::vgetcell(const TGridIndexVector& iv, bool FTflag)
{
	smallnat v,c,s,d;
	TGridIndex i,j;
#	ifdef VECTOR_MACHINE
	const int n = iv.length();
	for (c=0; c<ncd; c++) {
		VLOOPN(iv,n) {
			CT(c,v) = M(iv(v),c);
			assert(CT(0,v) > 0);
			assert(CT(3,v) == 0);
		}
	}
#	else
	VLOOP(iv) {
		i = iv(v);
		for (c=0; c<ncd; c++) CT(c,v) = M(i,c);
	}
#	endif
	if (!FTflag) return;
#	ifdef VECTOR_MACHINE
	for (d=0; d<dim; d++) {
		const int offset = surfdataoffset(d);
		const int aoffset = areadataoffset(d);
		for (s=0; s<nsd; s++) {
			VLOOPN(iv,n) {
				FT(d,1,s,v) = M(iv(v),offset+s);
			}
		}
		VLOOPN(iv,n) {
			AT(d,1,v) = M(iv(v),aoffset);
			assert(AT(d,1,v) > 0);
		}
		TGridIndexVector jv;
		vgetneighbour_s(jv,iv,d,0);
		for (s=0; s<nsd; s++) {
			VLOOPN(iv,n) {
				FT(d,0,s,v) = M(jv(v),offset+s);
				assert(jv(v) != NOINDEX);
			}
		}
		VLOOPN(iv,n) {
			AT(d,0,v) = M(jv(v),aoffset);
			assert(AT(d,0,v) > 0);
		}
	}
#	else
	VLOOP(iv) {
		i = iv(v);
		for (d=0; d<dim; d++) {
			const int offset = surfdataoffset(d);
			const int aoffset = areadataoffset(d);
			for (s=0; s<nsd; s++) FT(d,1,s,v) = M(i,offset+s);
			AT(d,1,v) = M(i,aoffset);
			assert(AT(d,1,v) > 0);
			j = getneighbour(i,d,0);
			for (s=0; s<nsd; s++) FT(d,0,s,v) = M(j,offset+s);
			AT(d,0,v) = M(j,aoffset);
			assert(AT(d,0,v) > 0);
		}
	}
#	endif
}

template <smallnat dim>
void TCartesianGrid<dim>::vgetcell_allow_noindex(const TGridIndexVector& iv)
{
	smallnat v,c;
	TGridIndex i;
#	ifdef VECTOR_MACHINE
	const int n = iv.length();
	for (c=0; c<ncd; c++) {
		VLOOPN(iv,n) {
			if (iv(v) != NOINDEX) CT(c,v) = M(iv(v),c);
		}
	}
#	else
	VLOOP(iv) {
		i = iv(v);
		if (i == NOINDEX)
			for (c=0; c<ncd; c++) CT(c,v) = M(i,c);
	}
#	endif
}

template <smallnat dim>
void TCartesianGrid<dim>::putsurf(TGridIndex i, TSurfSpec ss, smallnat v)
{
	smallnat s;
	const int offset = surfdataoffset(ss.d());
	const TGridIndex j = (ss.dir() ? i : getneighbour(i,ss.d(),0));
	pragma("_CRI ivdep");
	for (s=0; s<nsd; s++) Mset(j,offset+s, ST(0,s,v));
}


// Faster version of vputsurf than the trivial one in grid.C
template <smallnat dim>
void TCartesianGrid<dim>::vputsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv)
{
	smallnat v,s;
	TGridIndex i,j;
	TSurfSpec ss;
#	ifdef VECTOR_MACHINE
	smallnat svd[VECLEN];
	const smallnat n = iv.length();
	TGridIndexVector nei,offset(n),jvec(n);
	long int zero = 0;
	VLOOPN(iv,n) {
		svd[v] = sv.d(v);
	}
	vgetneighbour_dir0(nei,iv,svd);
	VLOOPN(iv,n) {
		jvec[v] = sv.dir(v) ? iv(v) : nei(v);
		offset[v] = surfdataoffset(svd[v]);
	}

	// A bit cryptic version. More readable is the if 0... endif below.
	TGridIndexVector s0offset(n);
	VLOOPN(iv,n) s0offset[v] = flatindex(jvec(v),offset(v));
	const TGridIndex soffset = compstep() /*&M(0,1) - &M(0,0)*/;
	for (s=0; s<nsd; s++) {
		VLOOPN(iv,n) data.put(s0offset(v) + soffset*s, ST(0,s,v));
	}
#	else
	VLOOP(iv) {
		i = iv(v);
		ss = sv(v);
		const int offset = surfdataoffset(ss.d());
		j = (ss.dir() ? i : getneighbour(i,ss.d(),0));
		for (s=0; s<nsd; s++) Mset(j,offset+s, ST(0,s,v));
	}
#	endif
}

template <smallnat dim>
void TCartesianGrid<dim>::getsurf(TGridIndex i, TSurfSpec ss, bool STflag, smallnat v)
{
	smallnat c;
	TGridIndex i1,i2;
	if (ss.dir()) {
		i1 = i;
		NN() = i2 = getneighbour(i,ss.d(),1);
	} else {
		i2 = i;
		NN() = i1 = getneighbour(i,ss.d(),0);
	}
	pragma("_CRI ivdep");
	for (c=0; c<ncd; c++) {
		UT(0,c,v) = M(i1,c);
		UT(1,c,v) = M(i2,c);
	}
	assert(UT(0,0,v) > 0);
	assert(UT(1,0,v) > 0);
	if (!STflag) return;
	smallnat s;
	const int offset = surfdataoffset(ss.d());
	pragma("_CRI ivdep");
	for (s=0; s<nsd; s++) ST(0,s,v) = M(i1,offset+s);
}

template <smallnat dim>
void TCartesianGrid<dim>::vgetsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv, bool STflag)
{
	const smallnat n = iv.length();
	TGridIndexVector i1(n),i2(n),nei;
	smallnat v,c;
	smallnat svd[VECLEN],svdir[VECLEN];
	VLOOPN(iv,n) {
		svd[v] = sv.d(v);
		svdir[v] = sv.dir(v);
	}
	vgetneighbour(nei,iv,svd,svdir);
	VLOOPN(iv,n) {
		if (svdir[v]) {
			i1[v] = iv(v);
			i2[v] = nei(v);
		} else {
			i2[v] = iv(v);
			i1[v] = nei(v);
		}
	}
	VLOOPN(iv,n) NN(v) = nei(v);
	i1.setlength(n); i2.setlength(n);
	for (c=0; c<ncd; c++) {
		VLOOPN(iv,n) {
			UT(0,c,v) = M(i1(v),c);
			UT(1,c,v) = M(i2(v),c);
		}
	}
#ifdef ASSERTIONS
	VLOOP(iv) {
		assert(UT(0,0,v) > 0);
		assert(UT(1,0,v) > 0);
	}
#endif
	if (!STflag) return;
	smallnat s;
	TGridIndexVector offset(n);
	VLOOPN(iv,n) offset[v] = surfdataoffset(svd[v]);
	for (s=0; s<nsd; s++) {
		VLOOPN(iv,n) ST(0,s,v) = M(i1(v),offset(v)+s);
	}
}

template <smallnat dim>
bool TCartesianGrid<dim>::streamload(istream& i, const Theader *hp1)
{
	Theader h;
	const Theader *hp;
	if (hp1)
		hp = hp1;
	else {
		i >> h;
		hp = &h;
	}
	if (strcmp(hp->getstr("type"),"cartesian")) {
		cerr << "*** TCartesianGrid::streamload: not a Cartesian grid file\n";
		return false;
	}
	char *const rfstr = hp->getstr("realformat");
	TRealFormat rf = REALFORMAT_ASCII;
	if (!strcmp(rfstr,"ascii"))
		rf = REALFORMAT_ASCII;
	else if (!strcmp(rfstr,"float"))
		rf = REALFORMAT_FLOAT;
	else if (!strcmp(rfstr,"double"))
		rf = REALFORMAT_DOUBLE;
	else {
		cerr << "*** TCartesianGrid::streamload: bad realformat = '" << rfstr << "'\n";
		return false;
	}
	TGridIndex j;
	smallnat a;
	switch (rf) {
	case REALFORMAT_DOUBLE:
		{
			double *buff = new double [ncd];
			for (j=0; j<ntot; j++) {
				ReadDoublesFromFile(i,buff,ncd);
				ByteConversion_input(sizeof(double),(unsigned char*)buff,ncd);
				for (a=0; a<ncd; a++) Mset(j,a, buff[a]);
			}
			delete [] buff;
		}
		break;
	case REALFORMAT_FLOAT:
		{
			float *buff = new float [ncd];
			for (j=0; j<ntot; j++) {
				ReadFloatsFromFile(i,buff,ncd);
				ByteConversion_input(sizeof(float),(unsigned char*)buff,ncd);
				for (a=0; a<ncd; a++) Mset(j,a, buff[a]);
			}
			delete [] buff;
		}
		break;
	case REALFORMAT_ASCII:
		for (j=0; j<ntot; j++) {
			for (a=0; a<ncd; a++) {real tmp; i >> tmp; Mset(j,a,tmp);}
		}
		break;
	case REALFORMAT_INQUIRE:
		break;
	}
	return i.good();
}

template <smallnat dim>
bool TCartesianGrid<dim>::streamsave(ostream& o, const char*) const
{
	save1(o);
	o << "type = cartesian\n";
	o << "eoh\n";
	TGridIndex i0,i,j;
	smallnat a;
	switch (output_realformat) {
	case REALFORMAT_DOUBLE:
		{
			double *buff = new double [VECLEN*ncd_saved];
			for (i0=0; i0<ntot; i0+= VECLEN) {
				const TGridIndex jlimit = min(ntot-i0,TGridIndex(VECLEN));
				for (a=0; a<ncd_saved; a++) {
					// We have a risk of bank conflicts if ncd_saved==8 for instance, but still
					// it is better to vectorize the long loop rather than a 8-length loop.
					VDIRS;
					for (j=0,i=i0; j<jlimit; j++,i++)
						buff[j*ncd_saved+a] = M(i,a);
				}
				ByteConversion(sizeof(double),(unsigned char*)buff,ncd_saved*jlimit);
				WriteDoublesToFile(o,buff,ncd_saved*jlimit);
			}
			delete [] buff;
		}
		break;
	case REALFORMAT_FLOAT:
		{
			float *buff = new float [ncd_saved*VECLEN];
			for (i0=0; i0<ntot; i0+= VECLEN) {
				const TGridIndex jlimit = min(ntot-i0,TGridIndex(VECLEN));
				for (a=0; a<ncd_saved; a++) {
					// We have a risk of bank conflicts if ncd_saved==8 for instance, but still
					// it is better to vectorize the long loop rather than a 8-length loop.
					VDIRS;
					for (j=0,i=i0; j<jlimit; j++,i++)
						buff[j*ncd_saved+a] = M(i,a);
				}
				ByteConversion(sizeof(float),(unsigned char*)buff,ncd_saved*jlimit);
				WriteFloatsToFile(o,buff,ncd_saved*jlimit);
			}
			delete [] buff;
		}
		break;
	case REALFORMAT_ASCII:
		for (i=0; i<ntot; i++) {
			for (a=0; a<ncd_saved; a++) o << ' ' << M(i,a);
			o << '\n';
		}
		break;
	case REALFORMAT_INQUIRE:
		break;
	}
	return o.good();
}

template <smallnat dim>
TGridIndex TCartesianGrid<dim>::find(const Tdimvec& X) const
{
	int ijk[3];
	smallnat d;
	for (d=0; d<dim; d++) {
		ijk[d] = TGridIndex(floor((X(d) - xmin[d])/dx)) + NB;
		if (ijk[d] < 0 || ijk[d] >= n123[d]) return NOINDEX;
	}
	return compose(ijk);
}

template <smallnat dim>
int TCartesianGrid<dim>::neighbour_timeclass_fix(TGridIndex i)
{
	smallnat d,dir;
	TGridIndex j;
	const unsigned short itc = timeclass(i);
	unsigned short jtc;
	unsigned short maxjtc = 0;
	for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
		j = getneighbour(i,d,dir);
		if (celltype(j) != INTERIOR_CELL && celltype(j) != GHOST_CELL) continue;
		jtc = timeclass(j);
		maxjtc = max(maxjtc, jtc);
	}
	if (maxjtc > itc+1) {
		set_timeclass(i, maxjtc-1);
		return maxjtc - itc;
	}
	return 0;
}

#if HAS_TEMPLATE_INIT_SYNTAX
  template class TCartesianGrid<1>;
# if MAXDIM >= 2
    template class TCartesianGrid<2>;
# endif
# if MAXDIM >= 3
    template class TCartesianGrid<3>;
# endif
#endif

