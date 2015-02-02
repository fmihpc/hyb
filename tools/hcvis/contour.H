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

#ifndef CONTOUR_H

#ifdef __GNUC__
#  pragma interface
#endif
#include <GL/gl.h>

class TContourSpec {
private:
	int N;
	real *contour;
	void init(real u1, real u2, int n=10);
	void init(real u1, real u2, real du);
public:
	TContourSpec() {N=1; contour=new real[1]; contour[0]=0;}
	TContourSpec(real u1, real u2, int n=10) {init(u1,u2,n);}
	TContourSpec(int n) {init(0,1,n);}
	void setcontour(int i, real val);
	void renew(real u1, real u2, int n=10) {delete [] contour; init(u1,u2,n);}
	void renew(real u1, real u2, real du) {delete [] contour; init(u1,u2,du);}
	int Ncontours() const {return N;}
	real nthcontour(int n) const {return contour[n];}
	~TContourSpec() {delete [] contour;}
};

extern void GLContourTriangle(real x1, real y1, real z1,
							  real x2, real y2, real z2,
							  real x3, real y3, real z3,
							  real u1, real u2, real u3,
							  const TContourSpec& cs, int dim_select);

extern void GLContourTriangle_3D(real x1, real y1, real z1,
								 real x2, real y2, real z2,
								 real x3, real y3, real z3,
								 real u1, real u2, real u3,
								 const TContourSpec& cs);

extern void GLContourSquare(real xc, real yc, real zc, real dx, const real u[4],
							const bool udens[4], const real ufacecenters[4],
							const TContourSpec& cs,
							const real ex[3], const real ey[3],		// unit vectors along X,Y of square
							int dim_select);	// dim_select=-1: 3D, dim_select=0: select y,z; =1: select x,z; =2: select x,y

extern void GLIsosurfTetrahedron(real x1, real y1, real z1,
								 real x2, real y2, real z2,
								 real x3, real y3, real z3,
								 real x4, real y4, real z4,
								 real u1, real u2, real u3, real u4,
								 const TContourSpec& cs,
								 void (*ComputeGradient)(real,real,real, real&,real&,real&));

extern void GLIsosurfCube_nosubdiv(real xc, real yc, real zc,
								   real dx,
								   const real u[8],
								   const TContourSpec& cs,
								   void (*ComputeGradient)(real,real,real, real&,real&,real&)=0);

extern void GLIsosurfCube(real xc, real yc, real zc,
						  real dx,
						  const real u[8],
						  const bool udens[6], const real ufacecenters[6],
						  const TContourSpec& cs,
						  real (*Interpolate)(real,real,real),
						  void (*ComputeGradient)(real,real,real, real&,real&,real&)=0);

#define CONTOUR_H

#endif
