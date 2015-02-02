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

#ifndef ZOOMSTACK_H

#ifdef __GNUC__
#  pragma interface
#endif

#define MAX_ZOOM_STACK 100

class TZoomStack {
private:
	real wmin[MAX_ZOOM_STACK],wmax[MAX_ZOOM_STACK], hmin[MAX_ZOOM_STACK],hmax[MAX_ZOOM_STACK];
	int n;
public:
	TZoomStack() {n=0;}
	void push(real wmin1,real wmax1, real hmin1,real hmax1);
	void pop(real& wmin1, real& wmax1, real& hmin1, real& hmax1);
	void reset() {n=0;}
	int isempty() {return n <= 0;}
};

#define ZOOMSTACK_H

#endif
