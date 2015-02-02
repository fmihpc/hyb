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
#  pragma implementation "gengrid.H"
#endif

#include "gengrid.H"

static real FindMachineEpsilon()
{
	// We have to use the volatile variable a in order to prevent compiler
	// from assigning it in register. On 386 (Linux) at least, the accuracy
	// of FP registers is greater than the accuracy of double variables,
	// and this would lead to too small epsilon being found.
	real eps = 1;
	volatile real a;
	do {
		eps/= 2;
		a = 1 + eps;
	} while (a > 1);
	eps*= 2;
	return eps;
}

static real FindSmallestRealNumber()
{
	real tiny = 1;
	volatile real a;
	do {
		tiny*= 0.125;
		a = tiny*0.125;
	} while (a != 0);
	return tiny;
}

real machine::eps;
real machine::tny;

machine the_machine;	// create one instance so that init() gets called and eps,tiny computed

void machine::init()
{
	eps = FindMachineEpsilon();
	tny = FindSmallestRealNumber();
	// Outputting to clog here causes a core dump on DEC probably because iostream
	// gets its global ctors called only after this one...
	// on sumppu the message gets printed 5 times, which is weird... but hopefully
	// does not harm.
//	clog << "machine::init(): eps=" << eps << ", tiny=" << tny << "\n" << flush;
}


void TGenericGrid::init(smallnat ncd1, smallnat nsd1, smallnat nfq1)
{
	machine::init();
	ncd = ncd1;
	ncd_saved = ncd;
	nsd = nsd1;
	nfq = nfq1;

	fluxtab.init(MAXDIM*2*nsd*VECLEN);
	areatab.init(MAXDIM*2*VECLEN);

	utab.init(2*ncd*VECLEN);
	celltab.init(ncd*VECLEN);
	surftab.init(MAXNEI*nsd*VECLEN);
	nntab.init(VECLEN);

	dirty = false;
}

TGenericGrid::~TGenericGrid()
{
	if (!dirty) {
//		delete [] surftab.RawPtr();
//		delete [] celltab.RawPtr();
//		delete [] utab.RawPtr();
//		delete [] areatab.RawPtr();
//		delete [] fluxtab.RawPtr();
	}
	Dcache.global_reset();
}

