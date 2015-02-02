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
#  pragma implementation "realptr.H"
#endif

#include "realptr.H"
#include "cache.H"
#include <stdlib.h>
#include <math.h>

real do_some_work(real x1)
{
	real x = x1;
	if (x < 0) x = 0; if (x > 0.1) x = 0.1;
	const int n = 150;
	int i;
	for (i=0; i<n; i++) {
		x = sin(sin(x));
	}
	recflops((2*flops_sqrt)*n);
	return 1e-30*x;
}

void TRealPtr::init(int n)
{
#ifdef REALARRAY_DEBUG
	if (ptr && len > 0) delete [] ptr;
	len = n;
#else
	if (ptr) delete [] ptr;
#endif
	ptr = new real [n];
}

#ifndef REALARRAY_DEBUG
void TCacheBasedRealPtr::init(int n)
{
	if (ptr) delete [] ptr;
	if (Dcache.bytes_free() >= n*int(sizeof(real)))
		ptr = Dcache.alloc(n*sizeof(real),"buff");
	else {
		cerr << "TCacheBasedRealPtr::init(" << n << ") warning: Out of Dcache, resorting to operator new\n";
		cerr << "  (Dcache.bytes_free() = " << Dcache.bytes_free() << ", required " << n*int(sizeof(real)) << "\n";
		ptr = new real [n];
	}
}
#endif


double flops = 0.0;

void doabort() {abort();}

#ifdef ASSERTIONS
void assert_action(const char *file, int lineno) {
	cerr << file << ":" << lineno << ": failed assertion\n";
	cout << flush;
	doabort();
}
#endif


