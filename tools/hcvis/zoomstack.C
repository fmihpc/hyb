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
#  pragma implementation "zoomstack.H"
#endif
#include "realptr.H"		// to get real
#include "zoomstack.H"

void TZoomStack::push(real wmin1,real wmax1, real hmin1,real hmax1)
{
	if (n < MAX_ZOOM_STACK) {
		wmin[n] = wmin1;
		wmax[n] = wmax1;
		hmin[n] = hmin1;
		hmax[n] = hmax1;
		n++;
	}
}

void TZoomStack::pop(real& wmin1, real& wmax1, real& hmin1, real& hmax1)
{
	if (n > 0) {
		n--;
		wmin1 = wmin[n];
		wmax1 = wmax[n];
		hmin1 = hmin[n];
		hmax1 = hmax[n];
	} else {
		cerr << "*** TZoomStack: trying to pop from empty stack\n";
	}
}

