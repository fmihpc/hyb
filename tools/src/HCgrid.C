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
#  pragma implementation "HCgrid.H"
#endif

#include <fstream>
#include <cmath>
#include <cstdio>
#include "HCgrid.H"
#include "fileheader.H"

template <smallnat dim>
void THCgrid<dim>::init(smallnat ncd1, smallnat nsd1, smallnat nfq1, TGridIndex maxnc1, bool use_mapping1)
{
	Tgrid<dim>::init(ncd1,nsd1,nfq1,maxnc1,use_mapping1);
	pool.init(maxnc1,
			  ncd1+dim*(nsd1+nsd_priv)+(use_mapping1 ? 2 : 1),	// clen_r
			  ablockpointer_offset()+dim						// clen_i
#if CACHE_NEIGHBOURS
			  +2*dim
#endif
			  ,
			  1<<dim,																	// c1
			  ((1<<(dim-1))*(nsd1+nsd_priv) > ncd1+dim*(nsd1+nsd_priv)) ? 2 : 1);		// c2
	ghostptrs.init(maxnc1);
	ghostptrs.zero();
#	if USE_SHMEM
	ghostptr_PEs.init(maxnc1);
	ghostptr_PEs.zero();
#	endif
	n1 = n2 =n3 = ntot = n123[0] = n123[1] = n123[2] = 1;
	dx = 0.01; cornercache.set_dx(dx);
	maxnei = (1 << (dim-1));
	// Initialize find_k_table, for function find_k_fast to work
	smallnat d,dir,ch;
	for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) for (ch=0; ch<nchildren; ch++)
		find_k_table[d][dir][ch] = find_k(d,dir,ch);
	TGridIndex i;
	pragma("_CRI ivdep");
	for (i=0; i<maxnc; i++)
		cellinfotab_put(i, 0);		// 0 is a suitable default cell type, see grid.H
}

template <smallnat dim>
void THCgrid<dim>::PEcoherency(int rootpe)
{
	pool.PEcoherency(rootpe);
	if (shmemmype() != rootpe) {
		// stuff set by regular():
		n1 = shmemget(&n1,rootpe);
		n2 = shmemget(&n2,rootpe);
		n3 = shmemget(&n3,rootpe);
		shmemget(&xmin[0],&xmin[0],3,rootpe);
		shmemget(&xmax[0],&xmax[0],3,rootpe);
		shmemget(&n123[0],&n123[0],3,rootpe);
		ntot = shmemget(&ntot,rootpe);
		dx = shmemget(&dx,rootpe); cornercache.set_dx(dx);
	}
}

template <smallnat dim>
void THCgrid<dim>::setup_after_regular()
{
	TGridIndex i;
	smallnat d;
	// NOTICE! Must call pool.alloc_base() already here, so that isover(i) produces correct result below!
	pool.alloc_base(ntot);
	VDIRS;
	for (i=first(); !isover(i); i=next(i)) {
		firstchildset(i,NOINDEX);
		parentset(i,NOINDEX);
		childorderset(i,0);
		levelset(i,0);
		markclear(i);
	}
#	ifdef VECTOR_MACHINE
	for (d=0; d<dim; d++) {
		VDIRS;
		for (i=first(); !isover(i); i=next(i)) a_block_set(i,d,NOINDEX);
	}
#	else
	for (i=first(); !isover(i); i=next(i)) for (d=0; d<dim; d++) a_block_set(i,d,NOINDEX);
#	endif
	shmembarrierall();
	if (mype == ROOT_PE) {
		for (i=first_1PE(); !isover(i); i=next_1PE(i)) update_cached(i);
	}
	shmembarrierall();
}

template <smallnat dim>
void THCgrid<dim>::update_cached(TGridIndex i, smallnat d, smallnat dir)
// Update neighbourdenseflags in cellinfo[i] in given direction (d,dir).
// If CAHCE_NEIGHBOURS is defined, update also neighbour pointers.
{
	real Xin[MAXDIM][VECLEN], Xout[MAXDIM][VECLEN];
	Tdimvec Xc;
	centroid(i,Xc);
	const real halfdx = 0.5*cellsize(i);
	const TGridIndex j = dogetneighbour(i,d,dir);
	// Notice that it is important that we use dogetneighbour rather than getneighbour here!
	if (merge_load) {
		// If merge_load is false, we just set the dense flags (see below)
#if CACHE_NEIGHBOURS
		neighbour_set(i,d,dir, j);
#endif
		real V;

		if (use_mapping) {

#if VECLEN < 8
#error VECLEN must not be smaller than 8
#endif
			real DX=0;
			switch (dim) {
			case 1:
				// In 1D the "volume" of a cell is the distance between the endpoints
				Xin[0][0] = Xc(0) - halfdx;
				Xin[0][1] = Xc(0) + halfdx;
				(*map)(Xin,Xout,2);
				V = fabs(Xout[0][1] - Xout[0][0]);
				// The cell physical size is the same as the cell "volume"
				DX = V;
				break;
			case 2:
			{
				// In 2D the "volume" of a cell is the area of the quadrilateral,
				// which is equal to half the cross product of the cross-sectional vectors (HT!)
				Xin[0][0] = Xc(0) - halfdx;
				Xin[1][0] = Xc(1) - halfdx;
				Xin[0][1] = Xc(0) + halfdx;
				Xin[1][1] = Xc(1) - halfdx;
				Xin[0][2] = Xc(0) - halfdx;
				Xin[1][2] = Xc(1) + halfdx;
				Xin[0][3] = Xc(0) + halfdx;
				Xin[1][3] = Xc(1) + halfdx;
				(*map)(Xin,Xout,4);
				const real r23_x = Xout[0][1] - Xout[0][2];
				const real r23_y = Xout[1][1] - Xout[1][2];
				const real r41_x = Xout[0][3] - Xout[0][0];
				const real r41_y = Xout[1][3] - Xout[1][0];
				V = 0.5*fabs(r23_x * r41_y - r23_y * r41_x);
				// The cell physical size is the minimum of the 4 difference distances
				const real rsqr12 = sqr(Xout[0][1] - Xout[0][0]) + sqr(Xout[1][1] - Xout[1][0]);
				const real rsqr23 = sqr(Xout[0][2] - Xout[0][1]) + sqr(Xout[1][2] - Xout[1][1]);
				const real rsqr34 = sqr(Xout[0][3] - Xout[0][2]) + sqr(Xout[1][3] - Xout[1][2]);
				const real rsqr41 = sqr(Xout[0][0] - Xout[0][3]) + sqr(Xout[1][0] - Xout[1][3]);
				DX = sqrt( min(min(rsqr12,rsqr23),min(rsqr34,rsqr41)) );
			}
			break;
			case 3:
			default:
			{
				smallnat ch,d;
				for (ch=0; ch<8; ch++) for (d=0; d<3; d++) {
					Xin[d][ch] = Xc(d) + (GetBit(d,ch) ? -halfdx : +halfdx);
				}
				(*map)(Xin,Xout,8);
				const real onesixth = 1.0/6.0;
				Tdimvec a,b,c;
#			define tripleproduct(a1,a2, b1,b2, c1,c2) \
				(   a[0] = Xout[0][a2] - Xout[0][a1],\
					a[1] = Xout[1][a2] - Xout[1][a1],\
					a[2] = Xout[2][a2] - Xout[2][a1],\
					b[0] = Xout[0][b2] - Xout[0][b1],\
					b[1] = Xout[1][b2] - Xout[1][b1],\
					b[2] = Xout[2][b2] - Xout[2][b1],\
					c[0] = Xout[0][c2] - Xout[0][c1],\
					c[1] = Xout[1][c2] - Xout[1][c1],\
					c[2] = Xout[2][c2] - Xout[2][c1],\
					a[0]*(b[1]*c[2] - b[2]*c[1]) + a[1]*(b[2]*c[0] - b[0]*c[2]) + a[2]*(b[0]*c[1] - b[1]*c[0]) )
					const real V0 = onesixth * fabs(tripleproduct(1,0, 2,0, 4,0));
				const real V3 = onesixth * fabs(tripleproduct(2,3, 1,3, 7,3));
				const real V5 = onesixth * fabs(tripleproduct(4,5, 7,5, 1,5));
				const real V6 = onesixth * fabs(tripleproduct(7,6, 4,6, 2,6));
				const real Vtetra = onesixth * fabs(tripleproduct(4,1, 7,1, 2,1));
				V = V0 + V3 + V5 + V6 + Vtetra;
				// The cell physical size is the minimum of the 12 difference distances
				bool FirstTime = true;
				smallnat ch1,ch2;
				for (ch1=0; ch1<8; ch1++) for (ch2=ch1+1; ch2<8; ch2++) {
					const real DX1 = sqrt(
						sqr(Xout[0][ch2] - Xout[0][ch1]) +
						sqr(Xout[1][ch2] - Xout[1][ch1]) +
						sqr(Xout[2][ch2] - Xout[2][ch1]) );
					if (FirstTime) DX=DX1; else DX=min(DX,DX1);
					FirstTime = false;
				}
			}
			break;
			}
			setphyssize(i,DX);

		} else {

		switch (dim) {
		case 1:
			V = cellsize(i);
			break;
		case 2:
			V = sqr(cellsize(i));
			break;
		case 3:
			V = cellsize(i)*sqr(cellsize(i));
			break;
		}

		}		// use_mapping
		setinvvolume(i,1.0/V);
	}	// end of if (merge_load)
	if (!isleaf(i)) return;
	if (j == NOINDEX) {
		cellinfotab_put(i, TCellInfo::set_dense(cellinfotab(i), d,dir, false));
		return;
	}
	const bool j_leaf = isleaf(j);
	cellinfotab_put(i, TCellInfo::set_dense(cellinfotab(i), d,dir, !j_leaf));
	if (!merge_load) return;
	if (dir == 1) {
		real i_area;
		//
		// In case of subdivided right neighbour we make the approximation that
		// the areas of the small faces are equal, that is, we do not compute
		// the mapping of the central point. Improving this would require:
		// (1) mapping the central point, (2) knowing which k value corresponds
		// to which neighbour.

		if (use_mapping) {

			// Compute i_area and normal vector n
			Tdimvec n;
			smallnat d1;
			for (d1=0; d1<MAXDIM; d1++) n[d1] = 0;
			switch (dim) {
			case 1:
				i_area = 1.0;
				n[0] = dir ? +1 : -1;
				break;
			case 2:
		    {
				smallnat d1;
				for (d1=0; d1<2; d1++) {
					if (d1 == d) {
						Xin[d1][0] = Xin[d1][1] = Xc(d1) + (dir ? +halfdx : -halfdx);
					} else {
						Xin[d1][0] = Xc(d1) - halfdx;
						Xin[d1][1] = Xc(d1) + halfdx;
					}
				}
				(*map)(Xin,Xout,2);
				// Xout[:][0] and Xout[:][1] are now the two endpoints of the cell edge.
				// i_area will be the distance between these endpoints:
				i_area = sqrt( sqr(Xout[0][1] - Xout[0][0]) + sqr(Xout[1][1] - Xout[1][0]) );
				// A normal vector can be found by making a cross product with ez.
				// The direction is determined by requiring that n.(0.5*(Xout[:][0]+Xout[:][1]) - Xc) > 0.
				// First, compute n = ez x (Xout[:][1] - Xout[:][0]).
				n[0] = -(Xout[1][1] - Xout[1][0]);
				n[1] = (Xout[0][1] - Xout[0][0]);
				// Second, change n's sign if necessary
				const real avepoint_x = 0.5*(Xout[0][0] + Xout[0][1]);
				const real avepoint_y = 0.5*(Xout[1][0] + Xout[1][1]);
				if (n(0)*(avepoint_x - Xc(0)) + n(1)*(avepoint_y - Xc(1)) <= 0) {
					n[0] = -n(0);
					n[1] = -n(1);
				}
			}
    		break;
			case 3:
		    {
				smallnat d1;
				bool dim_1 = true;
				for (d1=0; d1<3; d1++) {
					if (d1 == d) {
						Xin[d1][0] = Xin[d1][1] = Xin[d1][2] = Xin[d1][3] = Xc(d1) + (dir ? +halfdx : -halfdx);
					} else {
						if (dim_1) {
							Xin[d1][0] = Xc(d1) - halfdx;
							Xin[d1][1] = Xc(d1) + halfdx;
							Xin[d1][2] = Xc(d1) - halfdx;
							Xin[d1][3] = Xc(d1) + halfdx;
						} else {
							Xin[d1][0] = Xc(d1) - halfdx;
							Xin[d1][1] = Xc(d1) - halfdx;
							Xin[d1][2] = Xc(d1) + halfdx;
							Xin[d1][3] = Xc(d1) + halfdx;
						}
						dim_1 = false;
					}
				}
				(*map)(Xin,Xout,4);
				Tdimvec r12, r30;
				for (d1=0; d1<3; d1++) {
					r12[d1] = Xout[d1][1] - Xout[d1][2];
					r30[d1] = Xout[d1][3] - Xout[d1][0];
				}
				// n = cross product of r12, r30
				n[0] = r12[1]*r30[2] - r12[2]*r30[1];
				n[1] = r12[2]*r30[0] - r12[0]*r30[2];
				n[2] = r12[0]*r30[1] - r12[1]*r30[0];
				// i_area is half the cross product of the crosser vectors
				i_area = 0.5*sqrt( sqr(n[0]) + sqr(n[1]) + sqr(n[2]) );
				// The normal vector is now ready, except for the sign (and normalization, done after switch stmt).
				// The direction is determined by requiring that n.(0.25*(Xout[:][0]+Xout[:][1]+Xout[:][2]+Xout[:][3]) - Xc) > 0.
				Tdimvec avevec;
				for (d1=0; d1<3; d1++)
					avevec[d1] = 0.25*(Xout[d1][0] + Xout[d1][1] + Xout[d1][2] + Xout[d1][3]);
				if (n(0)*(avevec(0) - Xc(0)) + n(1)*(avevec(1) - Xc(1)) + n(2)*(avevec(2) - Xc(2)) <= 0) {
					n[0] = -n(0);
					n[1] = -n(1);
					n[2] = -n(2);
				}
		    }
		    break;
			}
			// Normalize n
			real nnorm = 0;
			for (d1=0; d1<dim; d1++) nnorm+= sqr(n(d1));
			nnorm = 1.0/sqrt(nnorm);
			for (d1=0; d1<dim; d1++) n[d1]*= nnorm;
			// Transfer vector n into grid structure
			if (j_leaf || dim == 1) {
				for (d1=0; d1<dim; d1++) pool.Mset(i,normaldataoffset(d)+d1, n[d1]);
			} else {
				TGridIndex ablock = a_block(i,d);
				if (ablock == NOINDEX) {
					ablock = pool.alloc_heap2();
					a_block_set(i,d,ablock);
				}
				smallnat k;
				for (d1=0; d1<dim; d1++) for (k=0; k<maxnei; k++)
					pool.Mset(ablock,ablock_normaldataoffset(k)+d1, n[d1]);
			}

		} else {

			switch (dim) {
			case 1: i_area = 1.0; break;
			case 2: i_area = cellsize(i); break;
			case 3: i_area = sqr(cellsize(i)); break;
			}

		}	// use_mapping

		if (j_leaf || dim == 1)
			pool.Mset(i,areadataoffset(d), i_area);
		else {
			switch (dim) {
			case 1: break;
			case 2: i_area*= 0.5; break;
			case 3: i_area*= 0.25; break;
			}
			TGridIndex ablock = a_block(i,d);
			if (ablock == NOINDEX) {
				ablock = pool.alloc_heap2();
				a_block_set(i,d,ablock);
			}
			smallnat k;
			for (k=0; k<maxnei; k++)
				pool.Mset(ablock,ablock_areadataoffset(k), i_area);
		}
	}
}

template <smallnat dim>
void THCgrid<dim>::update_cached(TGridIndex i)
// Update cached cell quantities in all direction (d,dir)
{
	smallnat d,dir;
	for (d=0; d<dim; d++) for (dir=0; dir<2; dir++)
		update_cached(i,d,dir);
}

#define buflen 128

#define VNEXT1(f) \
	ndone = 0;    \
	i = start;    \
	while (ndone < VECLEN && !isover(i)) { \
		maxj = DivNpes(Ncells() - i);      \
		if (buflen < maxj) maxj = buflen;  \
		pragma("_CRI ivdep");              \
		pragma("_CRI shortloop128");       \
		for (j=0,k=0; j<maxj; j++,i+=Npes) \
			if (isleaf(i)) buf[k++] = i;   \
		maxj = k;                          \
		if (VECLEN-ndone < maxj) maxj = VECLEN-ndone;          \
		pragma("_CRI ivdep");                                  \
		pragma("_CRI shortloop128");                           \
		for (j=0; j<maxj; j++) iv[ndone++] = buf[j];           \
	}                                      \
	iv.setlength(ndone)

/*
 * Cray compiler seems to produce core dumps if the above macro is used
 * instead of writing the text explicitly as we did below. Apparently some of
 * the pragmas are somehow line-specific, since also it gives a vectorization
 * info about only one loop in case of VNEXT1 even though the body contains
 * two loops.
 */

template <smallnat dim>
void THCgrid<dim>::vnext1(TGridIndexVector& iv, TGridIndex start, TGridIndex stride, Ttraversal trav) const
{
	TGridIndex i,j;
#ifdef VECTOR_MACHINE
	TGridIndex k,ndone,maxj,buf[buflen];
#endif
	switch (trav) {
	case TRAV_ALL:
		unroll4;
		for (i=start,j=0; !isover(i) && j<VECLEN; i+=stride) {
			iv[j++] = i;
		}
		iv.setlength(j);
		break;
	case TRAV_LEAF_ONLY:
#ifdef VECTOR_MACHINE
		ndone = 0;
		i = start;
		while (ndone < VECLEN && !isover(i)) {
			maxj = DivNpes(Ncells() - i);
			if (buflen < maxj) maxj = buflen;
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			unroll3;
			for (j=0,k=0; j<maxj; j++,i+=stride)
				if (isleaf(i) && celltype(i)!=DEAD_CELL && celltype(i)!=REMOVED_CELL) buf[k++] = i;
			/* now buf contains k (0<=k<=buflen) valid indices */
			maxj = VECLEN-ndone;
			if (k < maxj) maxj = k;
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			unroll4;
			for (j=0; j<maxj; j++) iv[ndone++] = buf[j];
		}
		iv.setlength(ndone);
#else
		unroll4;
		for (i=start,j=0; !isover(i) && j<VECLEN; i+=stride) {
			if (isleaf(i) && celltype(i)!=DEAD_CELL && celltype(i)!=REMOVED_CELL) iv[j++] = i;
		}
		iv.setlength(j);
#endif
		break;
	case TRAV_DOMAIN_ONLY:
#ifdef VECTOR_MACHINE
		ndone = 0;
		i = start;
		while (ndone < VECLEN && !isover(i)) {
			maxj = DivNpes(Ncells() - i);
			if (buflen < maxj) maxj = buflen;
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			unroll3;
			for (j=0,k=0; j<maxj; j++,i+=stride)
				if (isindomain(i)) buf[k++] = i;
			/* now buf contains k (0<=k<=buflen) valid indices */
			maxj = VECLEN-ndone;
			if (k < maxj) maxj = k;
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			unroll4;
			for (j=0; j<maxj; j++) iv[ndone++] = buf[j];
		}
		iv.setlength(ndone);
#else
		unroll4;
		for (i=start,j=0; !isover(i) && j<VECLEN; i+=stride) {
			if (isindomain(i)) iv[j++] = i;
		}
		iv.setlength(j);
#endif
		break;
	case TRAV_ACTIVE_ONLY:
#ifdef VECTOR_MACHINE
		ndone = 0;
		i = start;
		while (ndone < VECLEN && !isover(i)) {
			maxj = DivNpes(Ncells() - i);
			if (buflen < maxj) maxj = buflen;
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			unroll3;
			for (j=0,k=0; j<maxj; j++,i+=stride)
				if (isactive(i)) buf[k++] = i;
			/* now buf contains k (0<=k<=buflen) valid indices */
			maxj = VECLEN-ndone;
			if (k < maxj) maxj = k;
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			unroll4;
			for (j=0; j<maxj; j++) iv[ndone++] = buf[j];
		}
		iv.setlength(ndone);
#else
		unroll4;
		for (i=start,j=0; !isover(i) && j<VECLEN; i+=stride) {
			if (isactive(i)) iv[j++] = i;
		}
		iv.setlength(j);
#endif
		break;
	case TRAV_LEAFGHOST_ONLY:
#ifdef VECTOR_MACHINE
		ndone = 0;
		i = start;
		while (ndone < VECLEN && !isover(i)) {
			maxj = DivNpes(Ncells() - i);
			if (buflen < maxj) maxj = buflen;
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			unroll3;
			for (j=0,k=0; j<maxj; j++,i+=stride)
				if (isleaf(i) && celltype(i)==GHOST_CELL) buf[k++] = i;
			/* now buf contains k (0<=k<=buflen) valid indices */
			maxj = VECLEN-ndone;
			if (k < maxj) maxj = k;
			pragma("_CRI ivdep");
			pragma("_CRI shortloop128");
			unroll4;
			for (j=0; j<maxj; j++) iv[ndone++] = buf[j];
		}
		iv.setlength(ndone);
#else
		unroll4;
		for (i=start,j=0; !isover(i) && j<VECLEN; i+=stride) {
			if (isleaf(i) && celltype(i)==GHOST_CELL) iv[j++] = i;
		}
		iv.setlength(j);
#endif
		break;
	}
}

#undef buflen

template <smallnat dim>
smallnat THCgrid<dim>::MaximumLevel() const
{
	smallnat v, maxlevel = 0;
	TGridIndexVector iv;
	for (vfirst(iv); !iv.isempty(); vnext(iv)) {
		const smallnat n = iv.length();
		VLOOPN(iv,n) maxlevel = max(maxlevel,smallnat(level(iv(v))));
	}
	return maxlevel;
}

template <smallnat dim>
real THCgrid<dim>::MinimumGridSpacing() const
{
	const smallnat maxlevel = MaximumLevel();
	return dx*pow(0.5,real(maxlevel));
}

template <smallnat dim>
TGridIndex THCgrid<dim>::dogetneighbour(TGridIndex i,
										smallnat d, smallnat dir) const
{
	pool.clearcache();
	TGridIndex pc = parent(i);
	assert(i != NOINDEX);
	if (pc == NOINDEX) {
		int ijk[3];
		decompose(i, ijk);
		if (dir) ijk[d]++; else ijk[d]--;
		if (ijk[d] < 0 || ijk[d] >= n123[d]) {
			 // No neighbour, we are at the boundary of the basegrid
			return NOINDEX;
		}
		return compose(ijk);
	} else {
		// Find child cell index ch such that child(pc,ch) == c. Exactly one should be found.
		smallnat ch1,ch=0,nfound=0;
#		if USE_CSHMAT
		const TCacheIndex pcC = pool.load(pc);
		assert(firstchildC(pcC) == firstchild(pc));
#		endif
		for (ch1=0; ch1<nchildren; ch1++) {
#			if USE_CSHMAT
			assert(childC(pcC,ch1) == child(pc,ch1));
			if (childC(pcC,ch1) == i)
#			else
			if (child(pc,ch1) == i)
#			endif
			{
				ch = ch1;
				nfound++;
			}
		}
		if (nfound != 1) {
			cerr << "*** THCgrid<" << dim << ">::dogetneighbour: found "
				 << nfound << " matching children\n";
			doabort();
		}
		assert(child(pc,ch) == i);
		// Let b be the d'th bit of ch
		const bool b = GetBit(d,ch);
		// If b != dir, the neighbour is found among pc's children, and its ch-index ch_n
		// is obtained by inverting the d'th bit of ch
		if (b != dir) {
			const unsigned ch_n = InvertBit(d,ch);
			assert(child(pc,ch_n) != parent(i));
			assert(child(pc,ch_n) != i);
#			if USE_CSHMAT
			const TGridIndex retval = childC(pcC,ch_n);
			assert(retval == child(pc,ch_n));
//			pool.popcache();
			return retval;
#			else
			return child(pc,ch_n);
#			endif
		} else {
#			if USE_CSHMAT
//			pool.popcache();
#			endif
			// Else b==dir, and the neighbour is not found at this level
			// We find the neighbour of the parent cell pc.
			TGridIndex pcneighbour = dogetneighbour(pc,d,dir);
			if (pcneighbour == NOINDEX) return NOINDEX;
			if (firstchild(pcneighbour) != NOINDEX) {
				// Now the looked-for neighbour is among pcneighbour's children.
				// Its ch-index ch_n is obtained by inverting the d'th bit of ch
				const int ch_n = InvertBit(d,ch);
				assert(child(pcneighbour,ch_n) != parent(i));
				assert(child(pcneighbour,ch_n) != i);
				return child(pcneighbour,ch_n);
			} else {
				assert(pcneighbour != parent(i));
				assert(pcneighbour != i);
				return pcneighbour;
			}
		}
	}
}

template <smallnat dim>
int THCgrid<dim>::find_k(smallnat d, smallnat dir, smallnat ch) const
// This function is not properly defined for all values (d,dir,ch).
// We signal these cases by returning -1.
{
	smallnat k;
	for (k=0; k<maxnei; k++) 
		if (NeighbourChild(d,dir,k) == ch) return k;
	return -1;
}

template <smallnat dim>
inline smallnat THCgrid<dim>::find_k_fast(smallnat d, smallnat dir, smallnat ch) const
{
	const int retval = find_k_table[d][dir][ch];
	assert(retval >= 0);
	return smallnat(retval);
}

template <smallnat dim>
void THCgrid<dim>::fgetcell(TGridIndex i, smallnat v)		// same as getcell(i,false,v), but maybe faster
{
#	if USE_CSHMAT
	smallnat c;
	const TCacheIndex iC = pool.load(i);
	for (c=0; c<ncd; c++) CT(c,v) = pool.MC(iC,c);
#	else
	pool.M().get(pool.flatMindex(i,0),ncd,pool.compstep(),&CT(0,v),&CT(1,v)-&CT(0,v));
#	endif
}

template <smallnat dim>
void THCgrid<dim>::getcell(TGridIndex i, smallnat min_d, smallnat max_d, smallnat v)
{
	smallnat s,d,k;
#if USE_CSHMAT
	smallnat c;
	const TCacheIndex iC = pool.load(i);
	pragma("_CRI ivdep");
	for (c=0; c<ncd; c++) CT(c,v) = pool.MC(iC,c);
	real Asum,A,surfflux1;
	for (d=min_d; d<=max_d; d++) {
		const int offset = surfdataoffset(d);
		const int aoffset = areadataoffset(d);

		// left: dir=0
		const TGridIndex jleft = getneighbourC(iC,d,0);
		assert(jleft != NOINDEX);
		const TCacheIndex jleftC = pool.load(jleft);
		if (firstchildC(jleftC) == NOINDEX) {
			if (dim > 1 && levelC(jleftC) < levelC(iC)) {
				const smallnat the_k = find_k_fast(d,1,childorderC(iC));
				TGridIndex ablock = a_blockC(jleftC,d);
				assert(ablock != NOINDEX);
				const TCacheIndex ablockC = pool.load(ablock);
				for (s=0; s<nsd; s++) {
					FT(d,0,s,v) = pool.MC(ablockC,ablock_surfdataoffset(the_k)+s);
					assert(finite(FT(d,0,s,v)));
				}
				AT(d,0,v) = pool.MC(ablockC,ablock_areadataoffset(the_k));
			} else {
#ifdef DEBUG
				if (dim > 1) assert(levelC(jleftC) == levelC(iC));
#endif
				for (s=0; s<nsd; s++) {
					FT(d,0,s,v) = pool.MC(jleftC,offset+s);
					assert(finite(FT(d,0,s,v)));
				}
				AT(d,0,v) = pool.MC(jleftC,aoffset);
			}
		} else {
			// left neighbour has children, it is denser than i
			// accumulate everything to k=0 flux component
			for (s=0; s<nsd; s++) FT(d,0,s,v) = 0;
			Asum = 0;
			for (k=0; k<maxnei; k++) {
				const TGridIndex leftchild = childC(jleftC,NeighbourChild(d,0,k));
				const TCacheIndex leftchildC = pool.load(leftchild);
				A = pool.MC(leftchildC,aoffset);
				Asum+= A;
				for (s=0; s<nsd; s++) {
					surfflux1 = pool.MC(leftchildC,offset+s);
					FT(d,0,s,v)+= A*surfflux1;
					assert(finite(FT(d,0,s,v)));
				}
			}
			const real invAsum = 1.0/Asum;
			for (s=0; s<nsd; s++) FT(d,0,s,v)*= invAsum;
			AT(d,0,v) = Asum;
			recflops(maxnei*(1+2*nsd)+flops_div+nsd);
		}

		// right: dir=1
		if (dim==1 || !isdenseC(iC,d,1) /*firstchildC(jrightC) == NOINDEX*/) {
			if (dim > 1) assert(firstchild(getneighbourC(iC,d,1))==NOINDEX);
			// right neighbour not denser than i, or in 1D ==> do not need a-blocks
			for (s=0; s<nsd; s++) {
				FT(d,1,s,v) = pool.MC(iC,offset+s);
				assert(finite(FT(d,1,s,v)));
			}
			AT(d,1,v) = pool.MC(iC,aoffset);
		} else {
			// right neighbour is denser than i, need to access i's a-blocks
			// accumulate everything to k=0 flux component
			const TGridIndex ablock = a_blockC(iC,d);
			const TCacheIndex ablockC = pool.load(ablock);
			for (s=0; s<nsd; s++) FT(d,1,s,v) = 0;
			Asum = 0;
			for (k=0; k<maxnei; k++) {
				A = pool.MC(ablockC,ablock_areadataoffset(k));
				Asum+= A;
				for (s=0; s<nsd; s++) {
					FT(d,1,s,v)+= A*pool.MC(ablockC,ablock_surfdataoffset(k)+s);
					assert(finite(FT(d,1,s,v)));
				}
			}
			const real invAsum = 1.0/Asum;
			for (s=0; s<nsd; s++) FT(d,1,s,v)*= invAsum;
			AT(d,1,v) = Asum;
			recflops(maxnei*(1+2*nsd)+flops_div+nsd);
		}
	}
#else	/* USE_CSHMAT */
	pool.M().get(pool.flatMindex(i,0),ncd,pool.compstep(),&CT(0,v),&CT(1,v)-&CT(0,v));
	static real FTbuff[max_ncd];
	pragma("_CRI cache_align FTbuff");
	real Asum,A;
	for (d=min_d; d<=max_d; d++) {
		const int offset = surfdataoffset(d);
		const int aoffset = areadataoffset(d);

		// left: dir=0
		const TGridIndex jleft = getneighbour(i,d,0);
		assert(jleft != NOINDEX);
		if (firstchild(jleft) == NOINDEX) {
			// left neighbour not denser than i
			if (dim > 1 && level(jleft) < level(i)) {
				// left neighbour jleft is bigger than i
				const smallnat the_k = find_k_fast(d,1,childorder(i));
				TGridIndex ablock = a_block(jleft,d);
				assert(ablock != NOINDEX);
				pool.M().get(pool.flatMindex(ablock,ablock_surfdataoffset(the_k)),nsd,pool.compstep(),
							 &FT(d,0,0,v),&FT(d,0,1,v)-&FT(d,0,0,v));
				AT(d,0,v) = pool.M(ablock,ablock_areadataoffset(the_k));
			} else {
				// left neighbour jleft is of the same size as i (or 1D)
				if (dim > 1) assert(level(jleft) == level(i));
				pool.M().get(pool.flatMindex(jleft,offset),nsd,pool.compstep(),
							 &FT(d,0,0,v),&FT(d,0,1,v)-&FT(d,0,0,v));
				AT(d,0,v) = pool.M(jleft,aoffset);
			}
		} else {
			// left neighbour has children, it is denser than i
			// accumulate everything to k=0 flux component
			for (s=0; s<nsd; s++) FT(d,0,s,v) = 0;
			Asum = 0;
			for (k=0; k<maxnei; k++) {
				const TGridIndex leftchild = child(jleft,NeighbourChild(d,0,k));
				pool.M().get(pool.flatMindex(leftchild,offset),nsd,pool.compstep(),FTbuff,1);
				A = pool.M(leftchild,aoffset);
				Asum+= A;
				for (s=0; s<nsd; s++) FT(d,0,s,v)+= A*FTbuff[s];
			}
			const real invAsum = 1.0/Asum;
			for (s=0; s<nsd; s++) FT(d,0,s,v)*= invAsum;
			AT(d,0,v) = Asum;
			recflops(maxnei*(1+2*nsd)+flops_div+nsd);
		}

		// right: dir=1
		TGridIndex jright = getneighbour(i,d,1);
		assert(jright != NOINDEX);
		if (dim==1 || firstchild(jright) == NOINDEX) {
			// right neighbour not denser than i, or in 1D ==> do not need a-blocks
			pool.M().get(pool.flatMindex(i,offset),nsd,pool.compstep(),&FT(d,1,0,v),&FT(d,1,1,v)-&FT(d,1,0,v));
			AT(d,1,v) = pool.M(i,aoffset);
		} else {
			// right neighbour is denser than i, need to access i's a-blocks
			// accumulate everything to k=0 flux component
			const TGridIndex ablock = a_block(i,d);
			for (s=0; s<nsd; s++) FT(d,1,s,v) = 0;
			Asum = 0;
			for (k=0; k<maxnei; k++) {
				pool.M().get(pool.flatMindex(ablock,ablock_surfdataoffset(k)),nsd,pool.compstep(),FTbuff,1);
				A = pool.M(ablock,ablock_areadataoffset(k));
				Asum+= A;
				for (s=0; s<nsd; s++) FT(d,1,s,v)+= A*FTbuff[s];
			}
			const real invAsum = 1.0/Asum;
			for (s=0; s<nsd; s++) FT(d,1,s,v)*= invAsum;
			AT(d,1,v) = Asum;
			recflops(maxnei*(1+2*nsd)+flops_div+nsd);
		}
	}
#endif	/* USE_CSHMAT */
}

template <smallnat dim>
void THCgrid<dim>::getnormal(TGridIndex i, TSurfSpec s, real n3[3])
{
	smallnat d1;
	const smallnat d = s.d();
	const smallnat k = s.k();
	const smallnat dir = s.dir();
	for (d1=0; d1<3; d1++) n3[d1] = 0;
#if USE_CSHMAT
	const TCacheIndex iC = pool.load(i);
	if (dir == 0) {
		const TGridIndex jleft = getneighbourC(iC,d,0);
		const TCacheIndex jleftC = pool.load(jleft);
		if (firstchildC(jleftC) == NOINDEX) {
			// left neighbour not denser than i
			if (dim > 1 && levelC(jleftC) < levelC(iC)) {
				// left neighbour jleft is bigger than i
				const smallnat the_k = find_k_fast(d,1,childorderC(iC));
				TGridIndex ablock = a_blockC(jleftC,d);
				assert(ablock != NOINDEX);
				const TCacheIndex ablockC = pool.load(ablock);
				for (d1=0; d1<dim; d1++) n3[d1] = pool.MC(ablockC,ablock_normaldataoffset(the_k)+d1);
			} else {
				// left neighbour jleft is of the same size as i (or 1D)
				if (dim > 1) assert(levelC(jleftC) == levelC(iC));
				for (d1=0; d1<dim; d1++) n3[d1] = pool.MC(jleftC,normaldataoffset(d)+d1);
			}
		} else {
			// left neighbour has children, it is denser than i
			const TGridIndex leftchild = childC(jleftC,NeighbourChild(d,0,k));
			const TCacheIndex leftchildC = pool.load(leftchild);
			for (d1=0; d1<dim; d1++) n3[d1] = pool.MC(leftchildC,normaldataoffset(d)+d1);
		}
	} else {
		// right: dir=1
//		const TGridIndex jright = getneighbourC(iC,d,1);
//		const TCacheIndex jrightC = pool.load(jright);
		if (dim==1 || !isdenseC(iC,d,1) /*firstchildC(jrightC) == NOINDEX*/) {
			// right neighbour not denser than i, or in 1D ==> do not need a-blocks
			for (d1=0; d1<dim; d1++) n3[d1] = pool.MC(iC,normaldataoffset(d)+d1);
		} else {
			// right neighbour is denser than i, need to access i's a-blocks
			const TGridIndex ablock = a_blockC(iC,d);
			const TCacheIndex ablockC = pool.load(ablock);
			for (d1=0; d1<dim; d1++) n3[d1] = pool.MC(ablockC,ablock_normaldataoffset(k)+d1);
		}
	}
#else	/* USE_CSHMAT */
	if (dir == 0) {
		const TGridIndex jleft = getneighbour(i,d,0);
		if (firstchild(jleft) == NOINDEX) {
			// left neighbour not denser than i
			if (dim > 1 && level(jleft) < level(i)) {
				// left neighbour jleft is bigger than i
				const smallnat the_k = find_k_fast(d,1,childorder(i));
				TGridIndex ablock = a_block(jleft,d);
				assert(ablock != NOINDEX);
				for (d1=0; d1<dim; d1++) n3[d1] = pool.M(ablock,ablock_normaldataoffset(the_k)+d1);
			} else {
				// left neighbour jleft is of the same size as i (or 1D)
				if (dim > 1) assert(level(jleft) == level(i));
				for (d1=0; d1<dim; d1++) n3[d1] = pool.M(jleft,normaldataoffset(d)+d1);
			}
		} else {
			// left neighbour has children, it is denser than i
			const TGridIndex leftchild = child(jleft,NeighbourChild(d,0,k));
			for (d1=0; d1<dim; d1++) n3[d1] = pool.M(leftchild,normaldataoffset(d)+d1);
		}
	} else {
		// right: dir=1
		TGridIndex jright = getneighbour(i,d,1);
		if (dim==1 || firstchild(jright) == NOINDEX) {
			// right neighbour not denser than i, or in 1D ==> do not need a-blocks
			for (d1=0; d1<dim; d1++) n3[d1] = pool.M(i,normaldataoffset(d)+d1);
		} else {
			// right neighbour is denser than i, need to access i's a-blocks
			const TGridIndex ablock = a_block(i,d);
			for (d1=0; d1<dim; d1++) n3[d1] = pool.M(ablock,ablock_normaldataoffset(k)+d1);
		}
	}
#endif	/* USE_CSHMAT */
}

template <smallnat dim>
void THCgrid<dim>::putcell(TGridIndex i, smallnat v)
{
#	if USE_CSHMAT
	pool.Mset_multi(i,0,ncd, &CT(0,v),&CT(1,v)-&CT(0,v));
#	else
	pool.M().put(pool.flatMindex(i,0),ncd,pool.compstep(),&CT(0,v),&CT(1,v)-&CT(0,v));
#	endif
}

template <smallnat dim>
void THCgrid<dim>::getsurf(TGridIndex i, TSurfSpec ss, bool STflag, smallnat v)
{
	TGridIndex i1,i2;
#if USE_CSHMAT
	smallnat c;
	TCacheIndex i1C,i2C;
	const TCacheIndex iC = pool.load(i);
	assert(firstchildC(iC) == NOINDEX);
	const smallnat ch = NeighbourChild(ss.d(),ss.dir(),ss.k());
	if (ss.dir()) {
		i1 = i;
		i1C = iC;
		i2 = getneighbourC(iC,ss.d(),1);
		i2C = pool.load(i2);
		if (firstchildC(i2C) != NOINDEX) {i2 = childC(i2C,ch); i2C = pool.load(i2);}
		NN() = i2;
		assert(firstchildC(i2C) == NOINDEX);
	} else {
		i2 = i;
		i2C = iC;
		i1 = getneighbourC(iC,ss.d(),0);
		i1C = pool.load(i1);
		if (firstchildC(i1C) != NOINDEX) {i1 = childC(i1C,ch); i1C = pool.load(i1);}
		NN() = i1;
		assert(firstchildC(i1C) == NOINDEX);
	}
	pragma("_CRI ivdep");
	for (c=0; c<ncd; c++) {
		UT(0,c,v) = pool.MC(i1C,c);
		UT(1,c,v) = pool.MC(i2C,c);
	}

	if (!STflag) return;
	smallnat s;
	const int offset = surfdataoffset(ss.d());
	assert(firstchildC(iC) == NOINDEX);
	if (ss.dir()) {

		// right, may need to access i's a-block
//		const TGridIndex jright = getneighbourC(iC,ss.d(),1);
//		const TCacheIndex jrightC = pool.load(jright);
		if (dim==1 || !isdenseC(iC,ss.d(),1) /*firstchildC(jrightC) == NOINDEX*/) {
			// jright is of the same size or bigger than i
			assert(ss.k()==0);
			pragma("_CRI ivdep");
			for (s=0; s<nsd; s++) ST(0,s,v) = pool.MC(iC,offset+s);
		} else {
			// jright is denser than i, need to use i's a-block
			TGridIndex ablock = a_blockC(iC,ss.d());
			assert(ablock != NOINDEX);
			if (ablock == NOINDEX) {
				cerr << "getsurf warning: allocating new a-block (" << i << "," << ss.d() << ")\n";
				ablock = pool.alloc_heap2();
				a_block_set(i,ss.d(),ablock);
			}
			const TCacheIndex ablockC = pool.load(ablock);
			pragma("_CRI ivdep");
			for (s=0; s<nsd; s++) ST(0,s,v) = pool.MC(ablockC,ablock_surfdataoffset(ss.k())+s);
		}
	} else {

		cerr << "*** Leftside getsurf with STflag==true not yet implemented when USE_CSHMAT==1\n";
		doabort();
		
	}
#else	/* USE_CSHMAT */
	assert(firstchild(i) == NOINDEX);
	const smallnat ch = NeighbourChild(ss.d(),ss.dir(),ss.k());
	if (ss.dir()) {
		i1 = i;
		i2 = getneighbour(i,ss.d(),1);
		if (firstchild(i2) != NOINDEX) i2 = child(i2,ch);
		NN() = i2;
		assert(firstchild(i2) == NOINDEX);
	} else {
		i2 = i;
		i1 = getneighbour(i,ss.d(),0);
		if (firstchild(i1) != NOINDEX) i1 = child(i1,ch);
		NN() = i1;
		assert(firstchild(i1) == NOINDEX);
	}
	pool.M().get(pool.flatMindex(i1,0),ncd,pool.compstep(),&UT(0,0,v),&UT(0,1,v)-&UT(0,0,v));
	pool.M().get(pool.flatMindex(i2,0),ncd,pool.compstep(),&UT(1,0,v),&UT(1,1,v)-&UT(1,0,v));
	if (!STflag) return;
	const int offset = surfdataoffset(ss.d());
	assert(firstchild(i) == NOINDEX);
	if (ss.dir()) {

		// right, may need to access i's a-block
		const TGridIndex jright = getneighbour(i,ss.d(),1);
		if (dim==1 || firstchild(jright) == NOINDEX) {
			// jright is of the same size or bigger than i
			assert(ss.k()==0);
			pool.M().get(pool.flatMindex(i,offset),nsd,pool.compstep(),&ST(0,0,v),&ST(0,1,v)-&ST(0,0,v));
		} else {
			// jright is denser than i, need to use i's a-block
			TGridIndex ablock = a_block(i,ss.d());
			if (ablock == NOINDEX) {
				cerr << "getsurf warning: allocating new a-block (" << i << "," << ss.d() << ")\n";
				ablock = pool.alloc_heap2();
				a_block_set(i,ss.d(),ablock);
			}
			assert(ablock != NOINDEX);
			pool.M().get(pool.flatMindex(ablock,ablock_surfdataoffset(ss.k())),nsd,pool.compstep(),
						 &ST(0,0,v),&ST(0,1,v)-&ST(0,0,v));
		}

	} else {

		smallnat s;
		// left, may also need a-blocks
		const TGridIndex jleft = getneighbour(i,ss.d(),0);
		if (dim==1 || firstchild(jleft) == NOINDEX) {
			if (dim > 1 && level(jleft) < level(i)) {
				// left neighbour is bigger, need to use jleft's a-blocks
				const smallnat the_k = find_k_fast(ss.d(),1,childorder(i));
				TGridIndex ablock = a_block(jleft,ss.d());
				assert(ablock != NOINDEX);
				pragma("_CRI ivdep");
				for (s=0; s<nsd; s++) ST(0,s,v) = pool.M(ablock,ablock_surfdataoffset(the_k)+s);
			} else {
				// left neighbour jleft is of the same size as i
				assert(level(jleft) == level(i));
				assert(ss.k()==0);
				pragma("_CRI ivdep");
				for (s=0; s<nsd; s++) ST(0,s,v) = pool.M(jleft,offset+s);
			}
		} else {
			const TGridIndex leftchild = child(jleft,NeighbourChild(ss.d(),0,ss.k()));
			pragma("_CRI ivdep");
			for (s=0; s<nsd; s++) ST(0,s,v) = pool.M(leftchild,offset+s);
		}
		
	}
#endif
}

template <smallnat dim>
void THCgrid<dim>::putsurf(TGridIndex i, TSurfSpec ss, smallnat v)
{
	const int offset = surfdataoffset(ss.d());
#if USE_CSHMAT
	smallnat s;
	const TCacheIndex iC = pool.load(i);
	assert(firstchildC(iC) == NOINDEX);
	if (ss.dir()) {

		// right, may need to access i's a-block
//		const TGridIndex jright = getneighbourC(iC,ss.d(),1);
//		const TCacheIndex jrightC = pool.load(jright);
		if (dim==1 || !isdenseC(iC,ss.d(),1) /*firstchildC(jrightC) == NOINDEX*/) {
			// jright is of the same size or bigger than i
			assert(ss.k()==0);
			pragma("_CRI ivdep");
			for (s=0; s<nsd; s++) pool.Mset(i,offset+s, ST(ss.k(),s,v));
			assert(finite(ST(0,0,v)));
		} else {
			// jright is denser than i, need to use i's a-block
			TGridIndex ablock = a_block(i,ss.d());
			if (ablock == NOINDEX) {
				cerr << "putsurf warning: allocating new a-block (" << i << "," << ss.d() << ")\n";
				ablock = pool.alloc_heap2();
				a_block_set(i,ss.d(),ablock);
			}
			assert(ablock != NOINDEX);
			pragma("_CRI ivdep");
			for (s=0; s<nsd; s++) pool.Mset(ablock,ablock_surfdataoffset(ss.k())+s, ST(ss.k(),s,v));
			assert(finite(ST(ss.k(),0,v)));
		}

	} else {
		// left, may also need a-blocks
		const TGridIndex jleft = getneighbour(i,ss.d(),0);
		const TCacheIndex jleftC = pool.load(jleft);
		if (dim==1 || firstchildC(jleftC) == NOINDEX) {
			if (dim > 1 && levelC(jleftC) < levelC(iC)) {
				// left neighbour is bigger, need to use jleft's a-blocks
				const smallnat the_k = find_k_fast(ss.d(),1,childorderC(iC));
				TGridIndex ablock = a_blockC(jleftC,ss.d());
				assert(ablock != NOINDEX);
				pragma("_CRI ivdep");
				for (s=0; s<nsd; s++) pool.Mset(ablock,ablock_surfdataoffset(the_k)+s, ST(ss.k(),s,v));
				assert(finite(ST(ss.k(),0,v)));
			} else {
				// left neighbour jlet is of the same size as i
				assert(levelC(jleftC) == levelC(iC));
				assert(ss.k()==0);
				pragma("_CRI ivdep");
				for (s=0; s<nsd; s++) pool.Mset(jleft,offset+s, ST(ss.k(),s,v));
				assert(finite(ST(ss.k(),0,v)));
			}
		} else {
			const TGridIndex leftchild = childC(jleftC,NeighbourChild(ss.d(),0,ss.k()));
			pragma("_CRI ivdep");
			for (s=0; s<nsd; s++) pool.Mset(leftchild,offset+s, ST(ss.k(),s,v));
			assert(finite(ST(ss.k(),0,v)));
		}
	}
#else	/* USE_CSHMAT */
	assert(firstchild(i) == NOINDEX);
	if (ss.dir()) {

		// right, may need to access i's a-block
		const TGridIndex jright = getneighbour(i,ss.d(),1);
		if (dim==1 || firstchild(jright) == NOINDEX) {
			// jright is of the same size or bigger than i
			assert(ss.k()==0);
			pool.M().put(pool.flatMindex(i,offset),nsd,pool.compstep(),&ST(ss.k(),0,v),&ST(ss.k(),1,v)-&ST(ss.k(),0,v));
		} else {
			// jright is denser than i, need to use i's a-block
			TGridIndex ablock = a_block(i,ss.d());
			if (ablock == NOINDEX) {
				cerr << "putsurf warning: allocating new a-block (" << i << "," << ss.d() << ")\n";
				ablock = pool.alloc_heap2();
				a_block_set(i,ss.d(),ablock);
			}
			assert(ablock != NOINDEX);
			pool.M().put(pool.flatMindex(ablock,ablock_surfdataoffset(ss.k())),nsd,pool.compstep(),
						 &ST(ss.k(),0,v),&ST(ss.k(),1,v)-&ST(ss.k(),0,v));
		}

	} else {
		// left, may also need a-blocks
		const TGridIndex jleft = getneighbour(i,ss.d(),0);
		if (dim==1 || firstchild(jleft) == NOINDEX) {
			if (dim > 1 && level(jleft) < level(i)) {
				// left neighbour is bigger, need to use jleft's a-blocks
				const smallnat the_k = find_k_fast(ss.d(),1,childorder(i));
				TGridIndex ablock = a_block(jleft,ss.d());
				assert(ablock != NOINDEX);
				pool.M().put(pool.flatMindex(ablock,ablock_surfdataoffset(the_k)),nsd,pool.compstep(),
							 &ST(ss.k(),0,v),&ST(ss.k(),1,v)-&ST(ss.k(),0,v));
			} else {
				// left neighbour jlet is of the same size as i
				assert(level(jleft) == level(i));
				assert(ss.k()==0);
				pool.M().put(pool.flatMindex(jleft,offset),nsd,pool.compstep(),&ST(ss.k(),0,v),&ST(ss.k(),1,v)-&ST(ss.k(),0,v));
			}
		} else {
			const TGridIndex leftchild = child(jleft,NeighbourChild(ss.d(),0,ss.k()));
			pool.M().put(pool.flatMindex(leftchild,offset),nsd,pool.compstep(),&ST(ss.k(),0,v),&ST(ss.k(),1,v)-&ST(ss.k(),0,v));
		}
	}
#endif
}

#if !USE_CSHMAT

// -------------------- Vectorized functions -----------------------

template <smallnat dim>
void THCgrid<dim>::vgetcell(const TGridIndexVector& iv, bool FTflag)
// See also: vgetcell_allow_noindex
{
#ifdef DEBUG
	invalidateCT();
#endif
	smallnat c,s,d,v,k;
	const smallnat n = iv.length();
	for (c=0; c<ncd; c++) {
//		VLOOPN(iv,n) CT(c,v) = pool.M(iv(v),c);
		pool.M().gather(pool.flatMindex(0,c),iv,pool.indexMstep(),&CT(c,0),&CT(c,1)-&CT(c,0));
	}
	if (!FTflag) return;

	static bool leftisdenser[VECLEN];
	static bool leftisbigger[VECLEN];
	static bool rightisdenser[VECLEN];
	static smallnat vind[VECLEN];
	static TGridIndexVector jleft,jright,V,VA,ablock,leftchild,the_k;
	static real A[VECLEN],Asum[VECLEN],invAsum[VECLEN];
	jleft.setlength(n);
	jright.setlength(n);
	V.setlength(n);
	VA.setlength(n);
	ablock.setlength(n);
	leftchild.setlength(n);
	the_k.setlength(n);
	pragma("_CRI cache_align leftisdenser,leftisbigger,rightisdenser,vind,jleft,jright,V,VA,ablock,leftchild,the_k");
	pragma("_CRI cache_align A,Asum,invAsum");
	const int sstep = pool.compstep();

	for (d=0; d<dim; d++) {

		const int offset = surfdataoffset(d);
		const int aoffset = areadataoffset(d);

		// left: dir=0
		vgetneighbour_s(jleft,iv,d,0);
		bool any_leftisdenser = false;
		VLOOPN(iv,n) {
			leftisdenser[v] = (firstchild(jleft(v)) != NOINDEX);
			if (leftisdenser[v]) any_leftisdenser = true;
			leftchild[v] = leftisdenser[v] ? child(jleft(v),NeighbourChild(d,0,0)) : jleft(v);
		}
		if (dim == 1) {
			for (s=0; s<nsd; s++) {
				VLOOPN(iv,n) FT(d,0,s,v) = pool.M(leftchild(v),offset+s);
			}
			VLOOPN(iv,n) AT(d,0,v) = pool.M(leftchild(v),aoffset);
		} else {
			// dim > 1

			// Left2 contains 0.45 memory conflicts apparently because find_k_fast table gets
			// multiple reads since several iv(v) can map to same childorder
			VLOOPN(iv,n) {
				leftisbigger[v] = !leftisdenser[v] && level(jleft(v)) < level(iv(v));
				ablock[v] = leftisbigger[v] ? a_block(jleft(v),d) : v /*NOINDEX*/;
				the_k[v] = leftisbigger[v] ? find_k_fast(d,1,childorder(iv(v))) : 0;
			}

			VLOOPN(iv,n)
				V[v] =
				leftisbigger[v]
				? pool.flatMindex(ablock(v),ablock_surfdataoffset(the_k(v)))
				: pool.flatMindex(leftchild(v),offset);
			for (s=0; s<nsd; s++) {
//				VLOOPN(iv,n) FT(d,0,s,v) = pool.M(sstep*s+V(v));
				pool.M().gather(sstep*s,V,1,&FT(d,0,s,0),&FT(d,0,s,1)-&FT(d,0,s,0));
			}
			VLOOPN(iv,n)
				V[v] =
				leftisbigger[v]
				? pool.flatMindex(ablock(v),ablock_areadataoffset(the_k(v)))
				: pool.flatMindex(leftchild(v),aoffset);
//			VLOOPN(iv,n) AT(d,0,v) = pool.M(V(v));
			pool.M().gather(0,V,1,&AT(d,0,0),&AT(d,0,1)-&AT(d,0,0));

			if (any_leftisdenser) {
				// leftisdenser[vind[p]] will be true
				smallnat p = 0;
				VLOOPN(iv,n) if (leftisdenser[v]) vind[p++] = v;
				const smallnat Nvind = p;	// length of vector vind
				VDIRS;
				for (p=0; p<Nvind; p++) Asum[p] = 0;
				for (s=0; s<nsd; s++) {
					VDIRS;
					for (p=0; p<Nvind; p++) FT(d,0,s,vind[p]) = 0;
				}
				for (k=0; k<maxnei; k++) {
					VLOOPN(iv,n)
						if (leftisdenser[v])
							leftchild[v] = child(jleft(v),NeighbourChild(d,0,k));
					VLOOPN(iv,n) {
						V[v] = leftisdenser[v] ? pool.flatMindex(leftchild(v),offset) : v;
						VA[v] = leftisdenser[v] ? pool.flatMindex(leftchild(v),aoffset) : v;
					}
					VDIRS;
					for (p=0; p<Nvind; p++) {
						A[p] = pool.M(VA(vind[p]));
						Asum[p]+= A[p];
					}
					for (s=0; s<nsd; s++) {
						VDIRS;
						for (p=0; p<Nvind; p++) FT(d,0,s,vind[p])+= A[p]*pool.M(sstep*s + V(vind[p]));
					}
				}
				VDIRS;
				for (p=0; p<Nvind; p++) {
					AT(d,0,vind[p]) = Asum[p];
					invAsum[p] = 1.0/Asum[p];
				}
				for (s=0; s<nsd; s++) {
					VDIRS;
					for (p=0; p<Nvind; p++) FT(d,0,s,vind[p])*= invAsum[p];
				}
				recflops(Nvind*(1+2*nsd+flops_div+nsd));
			}
		}

		// right: dir=1
		vgetneighbour_s(jright,iv,d,1);
		if (dim == 1) {
			// In 1D ==> do not need a-blocks
//			VLOOPN(iv,n) NT(d,1,v) = 1;
			for (s=0; s<nsd; s++) {
//				VLOOPN(iv,n) FT(d,1,s,v) = pool.M(iv(v),offset+s);
				pool.M().gather(pool.flatMindex(0,offset+s),iv,pool.indexMstep(),&FT(d,1,s,0),&FT(d,1,s,1)-&FT(d,1,s,0));
			}
//			VLOOPN(iv,n) AT(d,1,v) = pool.M(iv(v),aoffset);
			pool.M().gather(pool.flatMindex(0,aoffset),iv,pool.indexMstep(),&AT(d,1,0),&AT(d,1,1)-&AT(d,1,0));
		} else {
			// dim > 1
			bool any_rightisdenser = false;
			VLOOPN(iv,n) {
				rightisdenser[v] = (firstchild(jright(v)) != NOINDEX);
				if (rightisdenser[v]) any_rightisdenser = true;
				ablock[v] = rightisdenser[v] ? a_block(iv(v),d) : NOINDEX;
			}

			VLOOPN(iv,n)
				V[v] =
				rightisdenser[v]
				? pool.flatMindex(ablock[v],ablock_surfdataoffset(0))
				: pool.flatMindex(iv(v),offset);
			for (s=0; s<nsd; s++) {
//				VLOOPN(iv,n) FT(d,1,s,v) = pool.M(s*sstep + V(v));
				pool.M().gather(s*sstep,V,1,&FT(d,1,s,0),&FT(d,1,s,1)-&FT(d,1,s,0));
			}
			VLOOPN(iv,n)
				V[v] =
				rightisdenser[v]
				? pool.flatMindex(ablock[v],ablock_areadataoffset(0))
				: pool.flatMindex(iv(v),aoffset);
//			VLOOPN(iv,n) AT(d,1,v) = pool.M(V(v));
			pool.M().gather(0,V,1,&AT(d,1,0),&AT(d,1,1)-&AT(d,1,0));
			if (any_rightisdenser) {
				smallnat vind[VECLEN];
				// rightisdenser[vind[p]] will be true
				smallnat p = 0;
				VLOOPN(iv,n) if (rightisdenser[v]) vind[p++] = v;
				const smallnat Nvind = p;	// length of vector vind
				VDIRS;
				for (p=0; p<Nvind; p++) Asum[p] = 0;
				for (s=0; s<nsd; s++) {
					VDIRS;
					for (p=0; p<Nvind; p++) FT(d,1,s,vind[p]) = 0;
				}
				for (k=0; k<maxnei; k++) {
					VLOOPN(iv,n) {
						V[v] = pool.flatMindex(ablock(v),ablock_surfdataoffset(k));
						VA[v] = pool.flatMindex(ablock(v),ablock_areadataoffset(k));
					}
					VDIRS;
					for (p=0; p<Nvind; p++) {
						A[p] = pool.M(VA(vind[p]));
						Asum[p]+= A[p];
					}
					for (s=0; s<nsd; s++) {
						VDIRS;
						for (p=0; p<Nvind; p++) FT(d,1,s,vind[p])+= A[p]*pool.M(s*sstep + V(vind[p]));
					}
				}
				VDIRS;
				for (p=0; p<Nvind; p++) {
					invAsum[p] = 1.0/Asum[p];
					AT(d,1,vind[p]) = Asum[p];
				}
				for (s=0; s<nsd; s++) {
					VDIRS;
					for (p=0; p<Nvind; p++) FT(d,1,s,vind[p])*= invAsum[p];
				}
				recflops(Nvind*(1+2*nsd+flops_div+nsd));
			}
		}
	}
}

template <smallnat dim>
void THCgrid<dim>::vgetcell_allow_noindex(const TGridIndexVector& iv)
// Same as vgetcell, but allow NOINDEX values in index vector
// (these values are simply ignored). FTflag is always false in this routine,
// so this routine does not transfer the cell-surrounding fluxes.
{
#ifdef DEBUG
	invalidateCT();
#endif
	smallnat c,v;
	const smallnat n = iv.length();
	/*
	 * Added if (iv(v) != NOINDEX) here. Easier to write HCtools.C:vComputeRefinementIndex
	 * So that it is possible to do g > iv when some of the iv's may be NOINDEX.
	 * The corresponding elements in CT are simply not assigned then.
	 * (NOTE: It MIGHT be possible just to do max(iv(v),0), then the elements would
	 *  be assigned, were legal memory addresses, but dummy.)
	 */
	for (c=0; c<ncd; c++) {
#if USE_SHMEM
		pool.M().gather(pool.flatMindex(0,c),iv,pool.indexMstep(),&CT(c,0),&CT(c,1)-&CT(c,0));
#else
		VLOOPN(iv,n) if (iv(v) != NOINDEX) CT(c,v) = pool.M(iv(v),c);
#endif
	}
}

template <smallnat dim>
void THCgrid<dim>::vgetnormal(const TGridIndexVector& iv, const TSurfSpecVector& sv, Tvec n3v[3])
{
	smallnat v,d1;
	TGridIndex i;
	real n3[3];
	VLOOP(iv) {
		i = iv(v);
		const smallnat d = sv.d(v);
		const smallnat k = sv.k(v);
		const smallnat dir = sv.dir(v);
		for (d1=0; d1<3; d1++) n3[d1] = 0;
		if (dir == 0) {
			const TGridIndex jleft = getneighbour(i,d,0);
			if (firstchild(jleft) == NOINDEX) {
				// left neighbour not denser than i
				if (dim > 1 && level(jleft) < level(i)) {
					// left neighbour jleft is bigger than i
					const smallnat the_k = find_k_fast(d,1,childorder(i));
					TGridIndex ablock = a_block(jleft,d);
					assert(ablock != NOINDEX);
					for (d1=0; d1<dim; d1++) n3[d1] = pool.M(ablock,ablock_normaldataoffset(the_k)+d1);
				} else {
					// left neighbour jleft is of the same size as i (or 1D)
					if (dim > 1) assert(level(jleft) == level(i));
					for (d1=0; d1<dim; d1++) n3[d1] = pool.M(jleft,normaldataoffset(d)+d1);
				}
			} else {
				// left neighbour has children, it is denser than i
				const TGridIndex leftchild = child(jleft,NeighbourChild(d,0,k));
				for (d1=0; d1<dim; d1++) n3[d1] = pool.M(leftchild,normaldataoffset(d)+d1);
			}
		} else {
			// right: dir=1
			TGridIndex jright = getneighbour(i,d,1);
			if (dim==1 || firstchild(jright) == NOINDEX) {
				// right neighbour not denser than i, or in 1D ==> do not need a-blocks
				for (d1=0; d1<dim; d1++) n3[d1] = pool.M(i,normaldataoffset(d)+d1);
			} else {
				// right neighbour is denser than i, need to access i's a-blocks
				const TGridIndex ablock = a_block(i,d);
				for (d1=0; d1<dim; d1++) n3[d1] = pool.M(ablock,ablock_normaldataoffset(k)+d1);
			}
		}
		n3v[0][v] = n3[0]; n3v[1][v] = n3[1]; n3v[2][v] = n3[2];
	}
}

template <smallnat dim>
void THCgrid<dim>::vputcell(const TGridIndexVector& iv)
{
	smallnat c;
	// Non-cryptic version
	for (c=0; c<ncd; c++) {
//		VLOOP(iv) pool.Mset(iv(v),c, CT(c,v));
		pool.M().scatter(pool.flatMindex(0,c),iv,pool.indexMstep(),&CT(c,0),&CT(c,1)-&CT(c,0));
	}

}

template <smallnat dim>
void THCgrid<dim>::vgetsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv, bool STflag)
{
#ifdef DEBUG
	invalidateUT();
#endif
	smallnat c,v;
	TGridIndex i,j;
	const smallnat n = iv.length();
	TGridIndexVector i1(n),i2(n),inei(n);
	smallnat dv[VECLEN],dirv[VECLEN];
	VLOOPN(iv,n) {
		dv[v] = sv.d(v);
		dirv[v] = sv.dir(v);
	}
	vgetneighbour(inei,iv,dv,dirv);
	VLOOPN(iv,n) {
		i = iv(v);
		assert(firstchild(i) == NOINDEX);
		const smallnat ch = NeighbourChild(sv.d(v),sv.dir(v),sv.k(v));
		j = inei(v);
		if (firstchild(j) != NOINDEX) j = child(j,ch);
		assert(firstchild(j) == NOINDEX);
		NN(v) = j;
		i1[v] = sv.dir(v) ? i : j;
		i2[v] = sv.dir(v) ? j : i;
	}
	
	// Cryptic version
	TGridIndexVector V1(n),V2(n);
	VLOOPN(iv,n) {
		V1[v] = pool.flatMindex(i1(v),0);
		V2[v] = pool.flatMindex(i2(v),0);
	}
	const int cstep = pool.compstep() /*&pool.M(0,1) - &pool.M(0,0)*/;
	for (c=0; c<ncd; c++) {
		pool.M().gather(c*cstep,V1,1,&UT(0,c,0),&UT(0,c,1)-&UT(0,c,0));
		pool.M().gather(c*cstep,V2,1,&UT(1,c,0),&UT(1,c,1)-&UT(1,c,0));
	}

	if (!STflag) return;

	smallnat s;
	VLOOPN(iv,n) assert(sv.dir(v) == 1);
	TGridIndexVector jright;
	smallnat unityvec[VECLEN];
	VLOOPN(iv,n) {
		unityvec[v] = 1;
	}
	vgetneighbour(jright,iv,dv,unityvec);
	bool rightisdenser[VECLEN];
	TGridIndexVector ablock(n);
	VLOOPN(iv,n) {
		rightisdenser[v] = (firstchild(jright(v)) != NOINDEX);
		if (rightisdenser[v]) ablock[v] = a_block(iv(v),sv.d(v));
	}
	if (dim == 1) {
		for (s=0; s<nsd; s++) {
			VLOOPN(iv,n) {
				const int offset = surfdataoffset(sv.d(v));
				assert(firstchild(iv(v)) == NOINDEX);
				assert(sv.k(v)==0);
				ST(0,s,v) = pool.M(iv(v),offset+s);
			}
		}
	} else {
		// dim > 1
		TGridIndexVector V(n);
		VLOOPN(iv,n){
			const int offset = surfdataoffset(sv.d(v));
			assert(firstchild(iv(v)) == NOINDEX);
			V[v] =
				rightisdenser[v]
				? pool.flatMindex(ablock(v),ablock_surfdataoffset(sv.k(v)))
				: pool.flatMindex(iv(v),offset);
		}
		const int sstep = pool.compstep() /*&pool.M(0,1) - &pool.M(0,0)*/;
		for (s=0; s<nsd; s++) {
			pool.M().gather(sstep*s,V,1,&ST(0,s,0),&ST(0,s,1)-&ST(0,s,0));
		}
	}
	
}

template <smallnat dim>
void THCgrid<dim>::vputsurf(const TGridIndexVector& iv, const TSurfSpecVector& sv)
{
	smallnat s,v;
	const smallnat n = iv.length();
	VLOOPN(iv,n) assert(sv.dir(v) == 1);
	TGridIndexVector jright;
	smallnat dv[VECLEN], unityvec[VECLEN];
	VLOOPN(iv,n) {
		dv[v] = sv.d(v);
		unityvec[v] = 1;
	}
	vgetneighbour(jright,iv,dv,unityvec);
	bool rightisdenser[VECLEN];
	TGridIndexVector ablock(n);
	VLOOPN(iv,n) {
		rightisdenser[v] = (firstchild(jright(v)) != NOINDEX);
		if (rightisdenser[v]) ablock[v] = a_block(iv(v),sv.d(v));
	}
	if (dim == 1) {
		for (s=0; s<nsd; s++) {
			VLOOPN(iv,n) {
				const int offset = surfdataoffset(sv.d(v));
				assert(firstchild(iv(v)) == NOINDEX);
				assert(sv.k(v)==0);
				pool.Mset(iv(v),offset+s, ST(sv.k(v),s,v));
			}
		}
	} else {
		// dim > 1
		TGridIndexVector V(n);
		VLOOPN(iv,n){
			const int offset = surfdataoffset(sv.d(v));
			assert(firstchild(iv(v)) == NOINDEX);
			V[v] =
				rightisdenser[v]
				? pool.flatMindex(ablock(v),ablock_surfdataoffset(sv.k(v)))
				: pool.flatMindex(iv(v),offset);
		}
		const int sstep = pool.compstep() /*&pool.M(0,1) - &pool.M(0,0)*/;
		for (s=0; s<nsd; s++) {
			real buff[VECLEN];
			VLOOPN(iv,n) buff[v] = ST(sv.k(v),s,v);
			pool.M().scatter(sstep*s,V,1,buff,1);
		}
	}
}

#endif		/* !USE_CSHMAT */

template <smallnat dim>
void THCgrid<dim>::centroid(TGridIndex i, Tdimvec& X) const
{
	// dX will be X - X[basegridcell]
	// we accumulate dX on our way up from i to basegrid cell
	TGridIndex pc = i, parpc;
	Tdimvec dX;
	smallnat d;
	for (d=0; d<dim; d++) dX[d] = 0;
	while ((parpc=parent(pc)) != NOINDEX) {
		const unsigned ch = childorder(pc);
		pc = parpc;
		for (d=0; d<dim; d++) dX[d]+= (GetBit(d,ch) ? +0.25 : -0.25)*cellsize(pc);
	}
	int ijk[3];
	decompose(pc, ijk);
	for (d=0; d<dim; d++)
		X[d] = xmin[d] + (ijk[d] - NB + 0.5)*dx + dX[d];
}

template <smallnat dim>
smallnat THCgrid<dim>::order(TGridIndex i) const
{
	smallnat d;
	TGridIndex j = i,jparent;
	int shifter = 0;
	int ijk[3] = {0,0,0};
	while ((jparent=parent(j)) != NOINDEX) {
		const unsigned ch = childorder(j);
		for (d=0; d<dim; d++) ijk[d]+= (GetBit(d,ch) << shifter);
		j = jparent;
		shifter++;
	}
	int ijk_base[3];
	decompose(j, ijk_base);
	for (d=0; d<dim; d++) ijk[d]+= (ijk_base[d] << shifter);
	smallnat result;
	switch (dim) {
	case 1: result = smallnat(ijk[0] % 3); break;
	case 2: result = smallnat((ijk[0] + 2*ijk[1]) % 5); break;
	case 3: result = smallnat((ijk[0] + 2*ijk[1] + 3*ijk[2]) % 7); break;
	}
	return result;
}

template <smallnat dim>
int THCgrid<dim>::child_interior_average(TGridIndex i, real uave[max_ncd])
{
	smallnat a,ch;
	TGridIndex j;
	real uave1[max_ncd];
	assert(i != NOINDEX);
	if (isleaf(i)) {
		if (celltype(i) == INTERIOR_CELL) {
			fgetcell(i);
			for (a=0; a<ncd; a++) uave[a] = CT(a);
			return 1;
		} else
			return 0;
	} else {
		// has children
		for (a=0; a<ncd; a++) uave[a] = 0;
		int result = 0;
		for (ch=0; ch<nchildren; ch++) {
			j = child(i,ch);
			const int result1 = child_interior_average(j,uave1);
			result+= result1;
			if (result1 > 0) for (a=0; a<ncd; a++) uave[a]+= uave1[a];
		}
		if (result > 0) {
			const real norm = 1.0/result;
			for (a=0; a<ncd; a++) uave[a]*= norm;
		}
		return result;
	}
}

template <smallnat dim>
inline unsigned char THCgrid<dim>::EncodeCellInfo(TGridIndex i) const
{
	unsigned char ch;
	ch = (unsigned char)celltype(i);
	if (may_subdivide(i)) ch|= 0x80;
	if (may_recoarsen(i)) ch|= 0x40;
	return ch;
}

inline TCellInfoType DecodeCellInfo(TCellInfoType i, unsigned char ch)
{
	TCellInfoType result = TCellInfo::set_celltype(i, TCellType(ch & 0x3F));
	result = TCellInfo::set_nosubdiv(result, (ch & 0x80) == 0);
	result = TCellInfo::set_norecoars(result, (ch & 0x40) == 0);
	return result;
}

template <smallnat dim>
bool THCgrid<dim>::streamload(istream& o, const Theader *hp1)
{
	Theader h;
	const Theader *hp;
	if (hp1)
		hp = hp1;
	else {
		o >> h;
		hp = &h;
	}
	if (!o.good()) return false;
	pool.set_merge_load(merge_load);
	bool retval = pool.streamload(o,*hp,parent_offset(),firstchild_offset());
	// Load cellinfo vector
	TGridIndex i,j;
	const TGridIndex ncells = pool.Ncells();
	for (i=0; i<ncells; i++) {
		TCellInfoType info = 0;
		if (pool.LoadedRealformat() == REALFORMAT_ASCII) {
			unsigned int ch;
			o >> ch;
			info = TCellInfo::set_celltype(info, TCellType(ch-'0'));
			o >> ch;
			info = TCellInfo::set_nosubdiv(info, ch == '0');
			o >> ch;
			info = TCellInfo::set_norecoars(info, ch == '0');
			o >> ch;
			info = TCellInfo::set_BCindex(info, ch);
		} else {
			info = DecodeCellInfo(info, o.get());
			info = TCellInfo::set_BCindex(info, (unsigned short)o.get());
		}
		cellinfotab_put(i,info);
	}
	
	// Setup childorder fields and mark bits
	const TGridIndex heap1 = pool.get_heap1();
	for (i=0; i<heap1; i++) {
		childorderset(i,0);
		markclear(i);
	}
	for (i=heap1; i<ncells; i++) {
		childorderset(i, (i-heap1) % nchildren);
		markclear(i);
	}
	// Setup level fields
	for (i=0; i<heap1; i++) levelset(i,0);
	for (i=heap1; i<ncells; i++) {
		TGridIndex lev = 0;
		j = i;
		// potentially infinite while-loop, if tree structure is broken
		while (parent(j) != NOINDEX) {j = parent(j); lev++;}
		levelset(i,lev);
	}
	// Setup a-block fields
	smallnat d;
	for (i=0; i<ncells; i++) for (d=0; d<dim; d++) a_block_set(i,d,NOINDEX);

	// Update cached values
	clearcache();
	for (i=0; i<ncells; i++)
		if (celltype(i) != REMOVED_CELL) update_cached(i);
	if (!merge_load) return retval;
	if (cell_prepare_hook)
		for (i=0; i<ncells; i++)
			if (isleaf(i) && celltype(i)==INTERIOR_CELL) (*cell_prepare_hook)(*this,i);
	clearcache();	// seems to be needed
	return retval;
}

template <smallnat dim>
bool THCgrid<dim>::streamsave(ostream& o, const char *filename_base) const
{
	TIndexTable newindex;
	TGridIndex i,ni;
	const TGridIndex ncells = pool.Ncells();
	TGridIndex i1=0, i2=ncells-1, j;
	if (parallel_io && Npes > 1) BlockWorkDivide(i1,i2);
	if (remove_gaps_when_saving) {
		// Do not save REMOVED_CELLs. Renumber cells while saving.
		// newindex(i) will be the gap-removed index.
		newindex.init(ncells);
		if (Npes == 1 || !parallel_io) {
			for (i=i1,ni=0; i<=i2; i++)
				newindex.put(i, (celltype(i) != REMOVED_CELL) ? ni++ : NOINDEX);
#ifdef DEBUG
			for (i=0; i<ncells; i++) {
				if (newindex(i) == NOINDEX) continue;
				assert((newindex(i) - i) % nchildren == 0);
			}
#endif
		} else {
			const TGridIndex ncells_loc = i2 - i1 + 1;
			TGridIndex *newindex_loc = new TGridIndex [ncells_loc];
			static TGridIndex niL;
			niL = 0;
			for (i=i1,j=0; i<=i2; i++,j++)
				newindex_loc[j] = (celltype(i) != REMOVED_CELL) ? niL++ : NOINDEX;
			// niL is the number of non-removed cells in our local newindex_loc vector
			// sum niL's from previous PEs to niLoffset
			shmembarrierall();
			int pe;
			TGridIndex niLoffset = 0;
			for (pe=0; pe<mype; pe++) niLoffset+= shmemget(&niL,pe);	// sum only pe < mype
			// add niLoffset to our newindex_loc
			for (j=0; j<ncells_loc; j++) if (newindex_loc[j] != NOINDEX) newindex_loc[j]+= niLoffset;
			// sum all newindex lengths
			shmemsumtoall(&niL);
			ni = niL;
			// Transfer from newindex_loc to newindex
			for (i=i1,j=0; i<=i2; i++,j++) newindex.put(i, newindex_loc[j]);
			delete [] newindex_loc;
		}
	} else {
		ni = ncells;
	}

	if (mype==0) {
		save1(o);
		o << "type = hc\n";
	}
	int *savedreals = new int [ncd];
	bool succ;
	static const int zero_one[2] = {0,1};
	for (i=0; i<ncd_saved; i++) savedreals[i] = i;
	shmembarrierall();		// wait until newindex is completely ready before going to pool.streamsave
	succ = pool.streamsave(o,newindex,ni,output_realformat,
						   parent_offset(),firstchild_offset(),
						   savedreals,ncd_saved,
						   zero_one,1,		// store parent_ptr_index for leaf cells
						   zero_one,2,		// store parent_ptr_index and child_ptr_index for nonleaf cells
						   (parallel_io && Npes > 1) ? filename_base : 0
		);
	if (succ) {
		// Save cell type and bcindex (two bytes per cell)
		if (MAX_N_BC > 256) {
			cerr << "*** MAX_N_BC is defined > 256, BCindex field won't fit in one byte, file recovery won't work.\n";
		}
		if (output_realformat == REALFORMAT_ASCII) {
			// ASCII output for scalar machine
			if (Npes > 1 && parallel_io && filename_base != 0) {
				// MPP machine (T3E): each PE writes its own file. Filename is ${filename_base}.tail.nnnn
				char *fn = new char [strlen(filename_base) + 20];
				sprintf(fn,"%s.tail.%.4d",filename_base,mype);
				ofstream o1(fn);
				for (i=i1; i<=i2; i++) {
					if (remove_gaps_when_saving) {if (newindex(i) == NOINDEX) continue;}
					o1 << (unsigned int)(celltype(i))
					   << ' ' << (unsigned int)(may_subdivide(i))
					   << ' ' << (unsigned int)(may_recoarsen(i))
					   << ' ' << (unsigned int)(BCindex(i))
					   << '\n';
				}
				delete [] fn;
			} else {
				// Scalar machine (non-MPP, or filename_base not given). The normal case.
				// Also MPP goes this way if -no_parallel_io was given.
				if (mype == 0)
					for (i=i1; i<=i2; i++) {
						if (remove_gaps_when_saving) {
							if (newindex(i) == NOINDEX) {
								continue;
							}
						}
						o << (unsigned int)(celltype(i))
						  << ' ' << (unsigned int)(may_subdivide(i))
						  << ' ' << (unsigned int)(may_recoarsen(i))
						  << ' ' << (unsigned int)(BCindex(i))
						  << '\n';
					}
				shmembarrierall();
			}
		} else {
			// Binary output for scalar machine
			if (Npes > 1 && parallel_io && filename_base != 0) {
				// MPP machine (T3E): each PE writes its own file. Filename is ${filename_base}.tail.nnnn
				char *fn = new char [strlen(filename_base) + 20];
				sprintf(fn,"%s.tail.%.4d",filename_base,mype);
				ofstream o1(fn);
				for (i=i1; i<=i2; i++) {
					if (remove_gaps_when_saving) {if (newindex(i) == NOINDEX) continue;}
					o1.put(EncodeCellInfo(i));
					o1.put(char(BCindex(i)));
				}
				delete [] fn;
			} else {
				// Binary output for scalar machine (non-MPP, or filename_base not given). The normal case.
				// Also MPP goes this way if -no_parallel_io was given. Thus we have to reset i1,i2
				// which were modified by BlockWorkDivide above.
				if (mype == 0)
					for (i=i1; i<=i2; i++) {
						if (remove_gaps_when_saving) {if (newindex(i) == NOINDEX) continue;}
						o.put(EncodeCellInfo(i));
						o.put(char(BCindex(i)));
					}
				shmembarrierall();
			}
		}
		succ = o.good();
	}
	delete [] savedreals;
	return succ;
}

template <smallnat dim>
TGridIndex THCgrid<dim>::find(const Tdimvec& X) const
{
	int ijk[3];
	int d;	// NOTE! must be signed type, since we do d>=0 test in the downgoing for loop below!
	const real invdx = 1.0/dx;
	for (d=0; d<dim; d++) {
		ijk[d] = TGridIndex(floor((X(d) - xmin[d])*invdx)) + NB;
		if (ijk[d] < 0 || ijk[d] >= n123[d]) return NOINDEX;
	}
	TGridIndex i = compose(ijk);		// this is the basegrid cell
	Tdimvec Xc;
	smallnat ch;
	int cnt=0;
	real current_dx = dx;
	while (!isleaf(i)) {
		ch = 0;
		if (cnt==0) centroid(i,Xc);
		// Avoid calling centroid for cnt>0. Instead we update Xc explicitly below.
		for (d=dim-1; d>=0; d--) {
			ch<<= 1;
			if (X(d) > Xc[d]) ch|= 1;
		}
		assert(ch < (1 << dim));
		i = child(i,ch);
		for (d=0; d<dim; d++) Xc[d]+= (GetBit(smallnat(d),ch) ? 0.25 : -0.25)*current_dx;
		current_dx*= 0.5;
		cnt++;
	}
	return i;
}

template <smallnat dim>
void THCgrid<dim>::mature_all_children(TGridIndex i)
{
	bool were_ghosts;
	if (isleaf(i))
		mature_cells(i,i,were_ghosts);
	else {
		smallnat ch;
		for (ch=0; ch<nchildren; ch++)
			mature_all_children(firstchild(i)+ch);
	}
}

template <smallnat dim>
void THCgrid<dim>::forbid_recoarsen()
{
	TGridIndex i;
	for (i=first(); !isover(i); i=next(i)) {
		if (!isleaf(i)) set_norecoars(i);
	}
}

bool can_do_induced_ghost_subdivisions = true;	// set this to false during body-refinement

template <smallnat dim>
TGridIndex THCgrid<dim>::subdivide(TGridIndex i, real t)
{
	markclear(i);
	if (!may_subdivide(i)) return NOINDEX;
	if (pool.freecells() < 0.01*maxnc) {
		warn("subdivision cancelled: out of memory");
		return NOINDEX;
	}
	smallnat d,dir,ch,a;
	TGridIndex j;
	// If already subdivided, do nothing
	if (!isleaf(i)) {
//		warn("cell already subdivided");
		return firstchild(i);
	}
	// If any neighbour is larger, subdivide it first. May also recurse.
	const TGridIndex ilevel = level(i);
	for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
		j = getneighbour(i,d,dir);
		if (j != NOINDEX && level(j) < ilevel) {
#			if PARALLEL_ADAPTATION
			warn("would need remnant induced subdivision with PARALLEL_ADAPTATION");
			return NOINDEX;
#			else
			subdivide(j,t);
#			endif
		}
	}
	// If any neighbour is still larger, do nothing.
	// If any neighbour is of equal size and ghost, subdivide it first
	for (dir=0; dir<2; dir++) for (d=0; d<dim; d++) {
		j = getneighbour(i,d,dir);
		if (j == NOINDEX) continue;
		if (level(j) < ilevel) {
			warn("cell neighbour is larger, subdivision forbidden");
			return NOINDEX;
		}
		if (can_do_induced_ghost_subdivisions
			&& celltype(i) == INTERIOR_CELL
			&& level(j) == ilevel
			&& isleaf(j)
			&& celltype(j) == GHOST_CELL) {
#			if PARALLEL_ADAPTATION
			if (!ismarked(j)) {
				// If j is marked it will be subdivided later, and this is perfectly safe.
				// So we do nothing here if j is marked.
				warn("would need remnant induced ghost-subdivision with PARALLEL_ADAPTATION");
				cerr << "ind-ghost-subdiv:: ";
				cerr << "i: "; writeinfo(cerr,i); cerr << ", j: "; writeinfo(cerr,j); cerr << "\n";
				// We continue with subdivide anyway, it is not a fatal problem to have ghost and interior
				// of differing subdivision levels (and near corners, this might be even necessary
				// in some cases).
			}
#			else
			subdivide(j,t /*,false*/);
#			endif
		}
	}
	const bool parent_was_ghost = (celltype(i) == GHOST_CELL);
	// Deallocate any a-blocks
	for (d=0; d<dim; d++)
		if (a_block(i,d) != NOINDEX) {
//			clog << "normal deallocation of a-block (" << i << "," << d << ")\n";
			pool.dealloc_heap2(a_block(i,d));
			a_block_set(i,d, NOINDEX);
		}
	TGridIndex result = pool.alloc_heap1();
	if (result == NOINDEX) {
		warn("subdivide: out of memory (this should not happen!)");
		return NOINDEX;
	}
	firstchildset(i, result);
	// Copy cell data from parent
//	*this > i;
	fgetcell(i);
	for (ch=0; ch<nchildren; ch++) {
		const TGridIndex child = result+ch;
		parentset(child, i);
		firstchildset(child, NOINDEX);
		levelset(child, level(i)+1);
		childorderset(child, ch);
		markclear(child);
		for (d=0; d<dim; d++) a_block_set(child,d,NOINDEX);
		// Copy parent's cell data to child
//		*this << child;
		putcell(child);
	}
	// Invalidate the parent data. This is not absolutely necessary,
	// but we do it for consistency. These quantities should never be
	// accessed.
	for (a=0; a<ncd; a++) CT(a) = NotANumber;
//	*this << i;
	putcell(i);
	// Allocate a-blocks if necessary
	if (dim > 1) {
		for (d=0; d<dim; d++) {
			const TGridIndex jleft = getneighbour(i,d,0);
			if (jleft != NOINDEX && isleaf(jleft)) {
				// left neighbour is less dense than newly subdivided siblings,
				// need to allocate a-block for left cell
				a_block_set(jleft,d, pool.alloc_heap2());
			}
		}
	}

	// Update cached quantities
	// Pass through neighbours of the parent cell and update pointers
	// pointing back towards the neighbour (that's why !dir)
	for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
		// Can use getneighbour here, since we can assume that the grid
		// neighbour pointers are in correct shape for the parent cell i.
		// In case of any doubt, one could substitute dogetneighbour here.
		j = getneighbour(i,d,dir);
		if (j == NOINDEX) continue;
		// If j is leaf, process it, otherwise process its children.
		// (only those that are on i's side).
		if (isleaf(j))
			update_cached(j,d,!dir);
		else {
			for (ch=0; ch<nchildren; ch++) {
				if (GetBit(d,ch) == !dir) update_cached(child(j,ch),d,!dir);
			}
		}
	}
	// Update the newly created cells in all directions.
	for (ch=0; ch<nchildren; ch++) update_cached(result+ch);
	bool were_ghosts;
	mature_cells(result,result+nchildren-1, were_ghosts);		// Mature the newly created children
	// Why do we need the following ?
#if !PARALLEL_ADAPTATION
	if (were_ghosts) {
		for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
			j = dogetneighbour(i,d,dir);
			if (j == NOINDEX) continue;
			mature_all_children(j);
		}
	}
#endif

	if (parent_was_ghost) {
		for (i=result; i<result+nchildren; i++)
			if (celltype(i) == INTERIOR_CELL) {
				for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
					j = dogetneighbour(i,d,dir);
					if (j == NOINDEX) continue;
					if (celltype(j) == DEAD_CELL) {set_celltype(j,GHOST_CELL); shmemquiet();}
					if (celltype(j) == GHOST_CELL && level(j) <= ilevel && !ismarked(j)) {
						warn("(harmless) Rare induced ghost subdivision");
						subdivide(j,t);
					}
				}
			}
	}
	for (i=result; i<result+nchildren; i++)
		if (celltype(i) == GHOST_CELL) applyBC(i,t);
	return result;
}

template <smallnat dim>
bool THCgrid<dim>::has_nonleaf_neighbours(TGridIndex i) const
{
	smallnat d,dir;
	TGridIndex j;
	for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
		j = getneighbour(i,d,dir);
		if (j == NOINDEX) continue;
		if (!isleaf(j)) return true;
	}
	return false;
}

template <smallnat dim>
bool THCgrid<dim>::children_have_nonleaf_neighbours(TGridIndex i) const
{
	smallnat ch;
	if (isleaf(i)) return false;
	for (ch=0; ch<nchildren; ch++) if (has_nonleaf_neighbours(child(i,ch))) return true;
	return false;
}

template <smallnat dim>
bool THCgrid<dim>::recoarsen(TGridIndex i, real t, bool may_recurse)
{
	markclear(i);
	if (!may_recoarsen(i)) return false;
	smallnat ch,d,a,dir;
	TGridIndex j;
#if !PARALLEL_ADAPTATION
	if (may_recurse && celltype(i) == GHOST_CELL) return false;
#endif
	if (isleaf(i)) {
		warn("recoarsening leaf cell, doing nothing");
		return false;
	}
	if (celltype(i) == DEAD_CELL || celltype(i) == REMOVED_CELL) {
		warn("recoarsening dead/removed cell, doing nothing");
		return false;
	}
	// Check that all children are leaf cells.
	// And, if any child has a nonleaf neighbour, do nothing.
	pragma("_CRI ivdep");
	for (ch=0; ch<nchildren; ch++) {
		j = child(i,ch);
		if (!isleaf(j)) {
			warn("(harmless) recoarsening a cell with a grandchild, doing nothing");
			return false;
		}
		if (has_nonleaf_neighbours(j)) {
#if PARALLEL_ADAPTATION
			warn("trying to recoarsen a cell whose child has a nonleaf neighbour");
#endif
			return false;
		}
	}
	// Now we know that recoarsening of i can be done, i.e.
	// children of i can be deallocated.
	if (may_recurse && celltype(i) == INTERIOR_CELL) {
		for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
			j = getneighbour(i,d,dir);
			if (j != NOINDEX && celltype(j)==GHOST_CELL && !isleaf(j)) {
#				if PARALLEL_ADAPTATION
				if (!ismarked(j)) {
					// If j is marked it will be recoarsened later, and this is perfectly safe.
					// So we do nothing here if j is marked.
					warn("would need remnant induced ghost-recoarsening with PARALLEL_ADAPTATION");
					cerr << "i: "; writeinfo(cerr,i); cerr << ", j: "; writeinfo(cerr,j); cerr << "\n";
					// we proceed anyway to recoarsen i. j is left unrecoarsened, but this should not be a very severe problem.
				}
#				else
				recoarsen(j,t,false);
#				endif
			}
		}
	}
	// Neighbours of i are either leaf cells, or have children, but no grandchildren.
	// If a left neighbour of i is a leaf cell, it has a-blocks which must be deleted.
	for (d=0; d<dim; d++) {
		j = getneighbour(i,d,0);
		if (j != NOINDEX && isleaf(j) && a_block(j,d) != NOINDEX) {
			pool.dealloc_heap2(a_block(j,d));
			a_block_set(j,d,NOINDEX);
		}
	}
	// If a right neighbour of i is non-leaf, a-blocks must be created for i.
	for (d=0; d<dim; d++) {
		j = getneighbour(i,d,1);
		if (j != NOINDEX && !isleaf(j)) a_block_set(i,d, pool.alloc_heap2());
	}
	// Synthesize parent value in u. Exclude dead cells from the averaging process.
	// Mark the children as removed at the same time.
	real *u = new real [ncd];
	pragma("_CRI shortloop"); pragma("_CRI ivdep");
	for (a=0; a<ncd; a++) u[a] = 0;
	real norm = 0;
	for (ch=0; ch<nchildren; ch++) {
		j = child(i,ch);
		if (celltype(j) != DEAD_CELL && celltype(j) != REMOVED_CELL) {
//			*this > j;
			fgetcell(j);
			pragma("_CRI shortloop"); pragma("_CRI ivdep");
			for (a=0; a<ncd; a++) u[a]+= CT(a);
			norm+= 1.0;
		}
		set_celltype(j,REMOVED_CELL);
		// Invalidate the child data. This is not absolutely necessary, but we do it for
		// consistency.
		for (a=0; a<ncd; a++) CT(a) = NotANumber;
//		*this << j;
		putcell(j);
	}
	if (norm != 0) norm = 1.0/norm;
	pragma("_CRI shortloop"); pragma("_CRI ivdep");
	for (a=0; a<ncd; a++) u[a]*= norm;
	// Put parent value in place
	pragma("_CRI shortloop"); pragma("_CRI ivdep");
	for(a=0; a<ncd; a++) CT(a) = u[a];
//	*this << i;
	putcell(i);
	delete [] u;
	const TGridIndex first_child = firstchild(i);
//	firstchildref(i) = NOINDEX;
	firstchildset(i,NOINDEX);
	// Update cached quantities, similar to what we did in subdivide() above.
	// Pass through neighbours of the parent cell and update pointers
	// pointing back towards the neighbour (that's why !dir)
	for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
		// Must use dogetneighbour here
		j = dogetneighbour(i,d,dir);
		if (j == NOINDEX) continue;
		if (isleaf(j))
			update_cached(j,d,!dir);
		else {
			for (ch=0; ch<nchildren; ch++) {
				if (GetBit(d,ch) == !dir) update_cached(child(j,ch),d,!dir);
			}
		}
	}
	// Update i in all directions.
	update_cached(i);
	pool.dealloc_heap1(first_child);
	bool were_ghosts;
	mature_cells(i,i, were_ghosts);		// Mature the recoarsened cell
	// The following was originally to turn those ghosts to dead that had lost contact
	// with interior cells due to subdivision. We now let these ghosts hang around
	// until they are detected in grid.C:applyBC. The same comments apply to
	// recoarsen() below.
#if !PARALLEL_ADAPTATION
	if (were_ghosts) {
		for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
			j = dogetneighbour(i,d,dir);
			if (j == NOINDEX) continue;
			mature_all_children(j);
		}
	}
#endif
	if (celltype(i) == GHOST_CELL) applyBC(i,t);
	return true;
}

template <smallnat dim>
void THCgrid<dim>::check()
// (1) Check that a-blocks are where they should be
// (2) Check that getneighbour works properly
{
	TGridIndex i;
	smallnat d,dir=1;
	if (dim == 1) return;	// no a-blocks in 1D grids
	for (i=first(); !isover(i); i=next(i)) {
		if (!isleaf(i) || celltype(i)==REMOVED_CELL) continue;
		for (d=0; d<dim; d++) {
			const TGridIndex j = getneighbour(i,d,dir);
			if (j == NOINDEX) continue;
			if (!isleaf(j)) {
				// right neighbour is denser
				if (a_block(i,d) == NOINDEX) {
					clog << "*** THCgrid<" << dim << ">::check: missing a-block (" << i << "," << d << ")\n";
				} else {
				}
			} else {
				if (a_block(i,d) != NOINDEX) {
					clog << "THCgrid<" << dim << ">::check: superfluous a-block (" << i << "," << d << ")\n";
				} else {
				}
			}
		}
	}
	// (2) Check that getneighbour works properly.
	// The maximum of the coordinate distances between cell i and neighbouring cell
	// should be equal to half the sum of cell sizes.
	for (i=first(); !isover(i); i=next(i)) {
		if (!isleaf(i) || celltype(i)==REMOVED_CELL) continue;
		Tdimvec Xi,Xj;
		centroid(i,Xi);
		for (d=0; d<dim; d++) for (dir=0; dir<2; dir++) {
			const TGridIndex j = getneighbour(i,d,dir);
			if (j == NOINDEX || celltype(j) == REMOVED_CELL) continue;
			centroid(j,Xj);
			smallnat d1;
			real maxdiff = 0; smallnat maxdir = 0;
			for (d1=0; d1<dim; d1++)
				if (fabs(Xi(d1) - Xj(d1)) > maxdiff) {
					maxdiff = fabs(Xi(d1) - Xj(d1));
					maxdir = d1;
				}
			if (fabs(maxdiff - 0.5*(cellsize(i) + cellsize(j))) > 1e-10) {
				clog << "*** THCgrid<" << dim << ">::check: cell centroid diff="
					 << maxdiff << ", should be " << 0.5*(cellsize(i) + cellsize(j)) << "\n";
				clog << "    Xi = (" << Xi(0) << ", " << Xi(1) << ")\n";
				clog << "    Xj = (" << Xj(0) << ", " << Xj(1) << ")\n";
			} else {
			}
		}
	}
}

template <smallnat dim>
void THCgrid<dim>::v_all_children_are_leaf(const TGridIndexVector& iv, TGridIndexVector& result) const
// Gathers those indices in result which (1) do have children, (2) but do NOT have grandchildren.
{
	smallnat v,ch,nresult=0;
	bool bresult[VECLEN], is_not_leaf[VECLEN];
	const smallnat n = iv.length();
	// At first: does not have children ==> bresult=false, otherwise true
	VLOOPN(iv,n) bresult[v] = is_not_leaf[v] = !isleaf(iv(v));
	// Then turn those bresult[v] false for which some children has children
	for (ch=0; ch<nchildren; ch++) {
		VLOOPN(iv,n)
			if (is_not_leaf[v]) {
				if (!isleaf(child(iv(v),ch))) bresult[v] = false;
			}
	}
	VLOOPN(iv,n) if (bresult[v]) result[nresult++] = iv(v);
	result.setlength(nresult);
}

template <smallnat dim>
bool THCgrid<dim>::iseven(TGridIndex i) const
{
	if (parent(i) == NOINDEX) {
		int ijk[3];
		decompose(i, ijk);
		smallnat d;
		TGridIndex sumindex = 0;
		for (d=0; d<dim; d++) sumindex+= ijk[d];
		return (sumindex & 0x1) == 0;
	} else {
		const TGridIndex order = childorder(i);
		if (dim == 1) {
			return (order == 0);
		} else if (dim == 2) {
			return (order == 0 || order == 3);
		} else {
			return (order == 0 || order == 3 || order == 5 || order == 6);
		}
	}
}

template <smallnat dim>
int THCgrid<dim>::neighbour_timeclass_fix(TGridIndex i)
{
	smallnat d,dir,k;
	TGridIndex j;
	const unsigned short itc = timeclass(i);
	unsigned short jtc;
	unsigned short maxjtc = 0;
	for (d=0; d<dim; d++) for (dir=0; dir<2; dir++)
		if (isdense(i,d,dir)) {
			for (k=0; k<maxnei; k++) {
				j = child(getneighbour(i,d,dir), NeighbourChild(d,dir,k));
				if (celltype(j) == INTERIOR_CELL || celltype(j) == GHOST_CELL) {
					jtc = timeclass(j);
					maxjtc = max(maxjtc, jtc);
				}
			}
		} else {
			j = getneighbour(i,d,dir);
			if (j != NOINDEX && (celltype(j) == INTERIOR_CELL || celltype(j) == GHOST_CELL)) {
				jtc = timeclass(j);
				maxjtc = max(maxjtc, jtc);
			}
		}
	if (maxjtc > itc+1) {
		set_timeclass(i, maxjtc-1);
		return maxjtc - itc;
	}
	return 0;
}

#if HAS_TEMPLATE_INIT_SYNTAX
  template class THCgrid<1>;
# if MAXDIM >= 2
    template class THCgrid<2>;
# endif
# if MAXDIM >= 3
    template class THCgrid<3>;
# endif
#endif

