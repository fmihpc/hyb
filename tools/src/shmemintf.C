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
#  pragma implementation "shmemintf.H"
#endif

#include "shmemintf.H"

#if USE_SHMEM

short int mype = 0;
#ifndef FIXED_LOG2_NPES
short int Npes = 1;
#endif

#if NPES_ONLY_POWER_OF_TWO
#ifndef FIXED_LOG2_NPES
short int log2_Npes = 0, Npes_minus_1 = 0;
#endif
#include <stdlib.h>		// for exit()
#endif

void startpes()
{
	mype = shmem_my_pe();
#ifdef FIXED_LOG2_NPES
	if (shmem_n_pes() != Npes && mype==0) {
		cerr << "*** startpes(): FIXED_LOG2_NPES=" << FIXED_LOG2_NPES << " but trying to run it for " << shmem_n_pes() << " PEs\n";
		exit(1);
	}
#else
	Npes = shmem_n_pes();
#	if NPES_ONLY_POWER_OF_TWO
	log2_Npes = -1;
	int x = Npes;
	while (x != 0) {
		x>>= 1;
		log2_Npes++;
	}
	if (Npes != (1 << log2_Npes)) {
		cerr << "*** shmemintf.C:startpes(): NPES_ONLY_POWER_OF_TWO is 1 but Npes=" << Npes << "\n";
		exit(1);
	}
	Npes_minus_1 = Npes - 1;
#	endif
#endif
}

#else

#include <math.h>
#include <string.h>

double xxx = 0.3432e-2;

inline void simulate_shmem_get_delay() {/*xxx = sin(sin(sin(sin(xxx))));*/}

#if SIMULATE_SLOW_SHMEM

float shmemget(const float *src, int) {simulate_shmem_get_delay(); return *src;}
double shmemget(const double *src, int) {simulate_shmem_get_delay(); return *src;}
short shmemget(const short *src, int) {simulate_shmem_get_delay(); return *src;}
int shmemget(const int *src, int) {simulate_shmem_get_delay(); return *src;}
long shmemget(const long *src, int) {simulate_shmem_get_delay(); return *src;}

void shmemget(float *target, const float *src, int len, int) {
	simulate_shmem_get_delay();
	memcpy(target,src,len*sizeof(float));
}

void shmemget(double *target, const double *src, int len, int) {
	simulate_shmem_get_delay();
	memcpy(target,src,len*sizeof(double));
}

void shmemget(short *target, const short *src, int len, int) {
	simulate_shmem_get_delay();
	memcpy(target,src,len*sizeof(short));
}

void shmemget(int *target, const int *src, int len, int) {
	simulate_shmem_get_delay();
	memcpy(target,src,len*sizeof(int));
}

void shmemget(long *target, const long *src, int len, int) {
	simulate_shmem_get_delay();
	memcpy(target,src,len*sizeof(long));
}

#endif	/* SIMULATE_SLOW_SHMEM */

void shmemiget(float *target, const float *src, int target_stride, int src_stride, int len, int) {
	int i;
	for (i=0; i<len; i++) target[i*target_stride] = src[i*src_stride];
}

void shmemiget(double *target, const double *src, int target_stride, int src_stride, int len, int) {
	int i;
	for (i=0; i<len; i++) target[i*target_stride] = src[i*src_stride];
}

void shmemiget(short *target, const short *src, int target_stride, int src_stride, int len, int) {
	int i;
	for (i=0; i<len; i++) target[i*target_stride] = src[i*src_stride];
}

void shmemiget(int *target, const int *src, int target_stride, int src_stride, int len, int) {
	int i;
	for (i=0; i<len; i++) target[i*target_stride] = src[i*src_stride];
}

void shmemiget(long *target, const long *src, int target_stride, int src_stride, int len, int) {
	int i;
	for (i=0; i<len; i++) target[i*target_stride] = src[i*src_stride];
}

#endif


