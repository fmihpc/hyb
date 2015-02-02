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
#  pragma implementation "contour.H"
#endif
#include "realptr.H"	// to get real
#include "contour.H"
#include <math.h>

inline int min(int x, int y) {return (x < y) ? x : y;}
inline int max(int x, int y) {return (x > y) ? x : y;}
inline real min(real x, real y) {return (x < y) ? x : y;}
inline real max(real x, real y) {return (x > y) ? x : y;}

void TContourSpec::init(real u1, real u2, int n)
{
	if (n < 2) {
		contour = new real [1];
		contour[0] = 0.5*(u1+u2);
		N = 1;
	} else {
		N = n;
		const real du = (u2-u1)/(n+1);
		contour = new real [N];
		int i;
		for (i=0; i<N; i++) contour[i] = u1 + (i+1)*du;
	}
}

void TContourSpec::init(real u1, real u2, real du)
{
	int i;
	if (u1 < 0 && u2 > 0) {
		const int Nnonneg = int(ceil(u2/du));
		const int Nneg = int(floor(-u1/du));
		N = Nnonneg + Nneg;
		contour = new real [N];
		for (i=0; i<Nneg; i++) contour[i] = -(Nneg-i)*du;
		for (i=0; i<Nnonneg; i++) contour[i+Nneg] = i*du;
	} else {
		real u;
		for (N=0,u=u1; u<u2; u+=du) N++;
		contour = new real [N];
		for (i=0; i<N; i++) contour[i] = u1 + i*du;
	}
}

void TContourSpec::setcontour(int i, real val)
{
	if (i < 0 || i >= N) {
		cerr << "*** TContourSpec::setcontour: index i=" << i << " out of range 0.." << N-1 << "\n";
		return;
	}
	contour[i] = val;
}

void GLContourTriangle(real x1, real y1, real z1,
					   real x2, real y2, real z2,
					   real x3, real y3, real z3,
					   real u1, real u2, real u3,
					   const TContourSpec& cs,
					   int dim_select)
{
	int i;
	const real min12 = min(u1,u2);
	const real min13 = min(u1,u3);
	const real min23 = min(u2,u3);
	const real max12 = max(u1,u2);
	const real max13 = max(u1,u3);
	const real max23 = max(u2,u3);
	const real umin = min(u1,min23);
	const real umax = max(u1,max23);
	real X[3],Y[3],Z[3];
	for (i=0; i<cs.Ncontours(); i++) {
		const real u = cs.nthcontour(i);
		if (u <= umin || u >= umax) continue;
		const int in_12 = min12 < u && u < max12;
		const int in_13 = min13 < u && u < max13;
		const int in_23 = min23 < u && u < max23;
		int ind = 0;
		if (in_12) {
			const real t = (u-min12)/(max12-min12);
			X[ind] = (u1 < u2) ? (1-t)*x1 + t*x2 : (1-t)*x2 + t*x1;
			Y[ind] = (u1 < u2) ? (1-t)*y1 + t*y2 : (1-t)*y2 + t*y1;
			Z[ind] = (u1 < u2) ? (1-t)*z1 + t*z2 : (1-t)*z2 + t*z1;
			ind++;
		}
		if (in_13) {
			const real t = (u-min13)/(max13-min13);
			X[ind] = (u1 < u3) ? (1-t)*x1 + t*x3 : (1-t)*x3 + t*x1;
			Y[ind] = (u1 < u3) ? (1-t)*y1 + t*y3 : (1-t)*y3 + t*y1;
			Z[ind] = (u1 < u3) ? (1-t)*z1 + t*z3 : (1-t)*z3 + t*z1;
			ind++;
		}
		if (in_23) {
			const real t = (u-min23)/(max23-min23);
			X[ind] = (u2 < u3) ? (1-t)*x2 + t*x3 : (1-t)*x3 + t*x2;
			Y[ind] = (u2 < u3) ? (1-t)*y2 + t*y3 : (1-t)*y3 + t*y2;
			Z[ind] = (u2 < u3) ? (1-t)*z2 + t*z3 : (1-t)*z3 + t*z2;
			ind++;
		}
		if (ind >= 2) {
			glBegin(GL_LINE_STRIP);
			switch (dim_select) {
			case -1:
				glVertex3f(X[0],Y[0],Z[0]);
				glVertex3f(X[1],Y[1],Z[1]);
				break;
			case 0:
				glVertex2f(Y[0],Z[0]);
				glVertex2f(Y[1],Z[1]);
				break;
			case 1:
				glVertex2f(X[0],Z[0]);
				glVertex2f(X[1],Z[1]);
				break;
			case 2:
				glVertex2f(X[0],Y[0]);
				glVertex2f(X[1],Y[1]);
				break;
			}
			glEnd();
		} 
		if (ind != 2) cerr << "contour warning: ind=" << ind << "\n";
	}
}

void GLContourTriangle_3D(real x1, real y1, real z1,
						  real x2, real y2, real z2,
						  real x3, real y3, real z3,
						  real u1, real u2, real u3,
						  const TContourSpec& cs)
{
	int i;
	const real min12 = min(u1,u2);
	const real min13 = min(u1,u3);
	const real min23 = min(u2,u3);
	const real max12 = max(u1,u2);
	const real max13 = max(u1,u3);
	const real max23 = max(u2,u3);
	const real umin = min(u1,min23);
	const real umax = max(u1,max23);
	real X[3],Y[3],Z[3];
	for (i=0; i<cs.Ncontours(); i++) {
		const real u = cs.nthcontour(i);
		if (u <= umin || u >= umax) continue;
		const int in_12 = min12 < u && u < max12;
		const int in_13 = min13 < u && u < max13;
		const int in_23 = min23 < u && u < max23;
		int ind = 0;
		if (in_12) {
			const real t = (u-min12)/(max12-min12);
			X[ind] = (u1 < u2) ? (1-t)*x1 + t*x2 : (1-t)*x2 + t*x1;
			Y[ind] = (u1 < u2) ? (1-t)*y1 + t*y2 : (1-t)*y2 + t*y1;
			Z[ind] = (u1 < u2) ? (1-t)*z1 + t*z2 : (1-t)*z2 + t*z1;
			ind++;
		}
		if (in_13) {
			const real t = (u-min13)/(max13-min13);
			X[ind] = (u1 < u3) ? (1-t)*x1 + t*x3 : (1-t)*x3 + t*x1;
			Y[ind] = (u1 < u3) ? (1-t)*y1 + t*y3 : (1-t)*y3 + t*y1;
			Z[ind] = (u1 < u3) ? (1-t)*z1 + t*z3 : (1-t)*z3 + t*z1;
			ind++;
		}
		if (in_23) {
			const real t = (u-min23)/(max23-min23);
			X[ind] = (u2 < u3) ? (1-t)*x2 + t*x3 : (1-t)*x3 + t*x2;
			Y[ind] = (u2 < u3) ? (1-t)*y2 + t*y3 : (1-t)*y3 + t*y2;
			Z[ind] = (u2 < u3) ? (1-t)*z2 + t*z3 : (1-t)*z3 + t*z2;
			ind++;
		}
		if (ind >= 2) {
			glBegin(GL_LINE_STRIP);
			glVertex3f(X[0],Y[0],Z[0]);
			glVertex3f(X[1],Y[1],Z[1]);
			glEnd();
		} 
		if (ind != 2) cerr << "contour warning: ind=" << ind << "\n";
	}
}

void GLContourSquare(real xc, real yc, real zc, real dx, const real u[4],
					 const bool udens[4], const real ufacecenters[4],
					 const TContourSpec& cs,
					 const real ex[3], const real ey[3],
					 int dim_select)
//
//  u2  --- u23 ---  u3
//   |       |        |
//  u02 --- uc  ---  u13
//   |       |        |
//  u0  --- u01 ---  u1
//
// udens, ufacecenters directions:
//
//       3
//       |
//  0  -----  1
//       |
//       2
{
	const real halfdx = 0.5*dx;
	const real u0 = u[0];
	const real u1 = u[1];
	const real u2 = u[2];
	const real u3 = u[3];
	const real u01 = ufacecenters[2];
	const real u02 = ufacecenters[0];
	const real u23 = ufacecenters[3];
	const real u13 = ufacecenters[1];
	const real uc = 0.25*(u0+u1+u2+u3);

	const real x0 = xc - halfdx*ex[0] - halfdx*ey[0];
	const real x1 = xc + halfdx*ex[0] - halfdx*ey[0];
	const real x2 = xc - halfdx*ex[0] + halfdx*ey[0];
	const real x3 = xc + halfdx*ex[0] + halfdx*ey[0];

	const real y0 = yc - halfdx*ex[1] - halfdx*ey[1];
	const real y1 = yc + halfdx*ex[1] - halfdx*ey[1];
	const real y2 = yc - halfdx*ex[1] + halfdx*ey[1];
	const real y3 = yc + halfdx*ex[1] + halfdx*ey[1];

	const real z0 = zc - halfdx*ex[2] - halfdx*ey[2];
	const real z1 = zc + halfdx*ex[2] - halfdx*ey[2];
	const real z2 = zc - halfdx*ex[2] + halfdx*ey[2];
	const real z3 = zc + halfdx*ex[2] + halfdx*ey[2];

	const real x01 = 0.5*(x0 + x1);
	const real y01 = 0.5*(y0 + y1);
	const real z01 = 0.5*(z0 + z1);

	const real x02 = 0.5*(x0 + x2);
	const real y02 = 0.5*(y0 + y2);
	const real z02 = 0.5*(z0 + z2);

	const real x13 = 0.5*(x1 + x3);
	const real y13 = 0.5*(y1 + y3);
	const real z13 = 0.5*(z1 + z3);

	const real x23 = 0.5*(x2 + x3);
	const real y23 = 0.5*(y2 + y3);
	const real z23 = 0.5*(z2 + z3);

#	define triangle(A,B,C) GLContourTriangle(x##A,y##A,z##A, x##B,y##B,z##B, x##C,y##C,z##C, u##A,u##B,u##C,cs,dim_select)
	if (udens[0]) {
		triangle(c,2,02);
		triangle(0,c,02);
	} else
		triangle(c,0,2);

	if (udens[2]) {
		triangle(0,c,01);
		triangle(01,c,1);
	} else
		triangle(c,0,1);

	if (udens[1]) {
		triangle(1,c,13);
		triangle(c,13,3);
	} else
		triangle(c,1,3);

	if (udens[3]) {
		triangle(c,3,23);
		triangle(c,23,2);
	} else
		triangle(c,2,3);
#	undef triangle
}

inline void normalize(real& nx, real& ny, real& nz)
{
	real nnorm = sqrt(nx*nx + ny*ny + nz*nz);
	nnorm = 1.0/nnorm;
	nx*= nnorm; ny*= nnorm; nz*= nnorm;
}

void GLIsosurfTetrahedron(real x1, real y1, real z1,
						  real x2, real y2, real z2,
						  real x3, real y3, real z3,
						  real x4, real y4, real z4,
						  real u1, real u2, real u3, real u4,
						  const TContourSpec& cs,
						  void (*ComputeGradient)(real,real,real, real&,real&,real&))
{
	int i;
	const real min12 = min(u1,u2);
	const real min34 = min(u3,u4);
	const real max12 = max(u1,u2);
	const real max34 = max(u3,u4);
	const real umin = min(min12,min34);
	const real umax = max(max12,max34);
	real X[6],Y[6],Z[6];
	for (i=0; i<cs.Ncontours(); i++) {
		const real u = cs.nthcontour(i);
		if (u <= umin || u >= umax) continue;
		const real min13 = min(u1,u3);
		const real min14 = min(u1,u4);
		const real min23 = min(u2,u3);
		const real min24 = min(u2,u4);
		const real max13 = max(u1,u3);
		const real max14 = max(u1,u4);
		const real max23 = max(u2,u3);
		const real max24 = max(u2,u4);
		const int in_12 = min12 < u && u < max12;
		const int in_13 = min13 < u && u < max13;
		const int in_14 = min14 < u && u < max14;
		const int in_23 = min23 < u && u < max23;
		const int in_24 = min24 < u && u < max24;
		const int in_34 = min34 < u && u < max34;
		int ind = 0;
		if (in_12) {
			const real t = (u-min12)/(max12-min12);
			X[ind] = (u1 < u2) ? (1-t)*x1 + t*x2 : (1-t)*x2 + t*x1;
			Y[ind] = (u1 < u2) ? (1-t)*y1 + t*y2 : (1-t)*y2 + t*y1;
			Z[ind] = (u1 < u2) ? (1-t)*z1 + t*z2 : (1-t)*z2 + t*z1;
			ind++;
		}
		if (in_13) {
			const real t = (u-min13)/(max13-min13);
			X[ind] = (u1 < u3) ? (1-t)*x1 + t*x3 : (1-t)*x3 + t*x1;
			Y[ind] = (u1 < u3) ? (1-t)*y1 + t*y3 : (1-t)*y3 + t*y1;
			Z[ind] = (u1 < u3) ? (1-t)*z1 + t*z3 : (1-t)*z3 + t*z1;
			ind++;
		}
		if (in_14) {
			const real t = (u-min14)/(max14-min14);
			X[ind] = (u1 < u4) ? (1-t)*x1 + t*x4 : (1-t)*x4 + t*x1;
			Y[ind] = (u1 < u4) ? (1-t)*y1 + t*y4 : (1-t)*y4 + t*y1;
			Z[ind] = (u1 < u4) ? (1-t)*z1 + t*z4 : (1-t)*z4 + t*z1;
			ind++;
		}
		if (in_23) {
			const real t = (u-min23)/(max23-min23);
			X[ind] = (u2 < u3) ? (1-t)*x2 + t*x3 : (1-t)*x3 + t*x2;
			Y[ind] = (u2 < u3) ? (1-t)*y2 + t*y3 : (1-t)*y3 + t*y2;
			Z[ind] = (u2 < u3) ? (1-t)*z2 + t*z3 : (1-t)*z3 + t*z2;
			ind++;
		}
		if (in_24) {
			const real t = (u-min24)/(max24-min24);
			X[ind] = (u2 < u4) ? (1-t)*x2 + t*x4 : (1-t)*x4 + t*x2;
			Y[ind] = (u2 < u4) ? (1-t)*y2 + t*y4 : (1-t)*y4 + t*y2;
			Z[ind] = (u2 < u4) ? (1-t)*z2 + t*z4 : (1-t)*z4 + t*z2;
			ind++;
		}
		if (in_34) {
			const real t = (u-min34)/(max34-min34);
			X[ind] = (u3 < u4) ? (1-t)*x3 + t*x4 : (1-t)*x4 + t*x3;
			Y[ind] = (u3 < u4) ? (1-t)*y3 + t*y4 : (1-t)*y4 + t*y3;
			Z[ind] = (u3 < u4) ? (1-t)*z3 + t*z4 : (1-t)*z4 + t*z3;
			ind++;
		}

#		define vert(n) \
			if (ComputeGradient) {(*ComputeGradient)(X[n],Y[n],Z[n],nx,ny,nz); nx=-nx; ny=-ny; nz=-nz;\
			normalize(nx,ny,nz); glNormal3f(nx,ny,nz);} glVertex3f(X[n],Y[n],Z[n])

		real nx,ny,nz;
		if (!ComputeGradient) {
			// No gradient calculator function provided, compute gradient by finite differencing the 4 tetrad points
			const real u1mod = u2 - u1;
			const real u2mod = u3 - u1;
			const real u3mod = u4 - u1;
			const real x1_ = x2 - x1;
			const real x2_ = x3 - x1;
			const real x3_ = x4 - x1;
			const real y1_ = y2 - y1;
			const real y2_ = y3 - y1;
			const real y3_ = y4 - y1;
			const real z1_ = z2 - z1;
			const real z2_ = z3 - z1;
			const real z3_ = z4 - z1;
			const real det = -(x3_*y2_*z1_) + x2_*y3_*z1_ + x3_*y1_*z2_ - x1_*y3_*z2_ - x2_*y1_*z3_ + x1_*y2_*z3_;
			const real invdet = det;		// actually, invdet=1.0/det, but since we normalize below, sign is enough!

			real nx = -(-(u3mod*y2_*z1_) + u2mod*y3_*z1_ + u3mod*y1_*z2_ - u1mod*y3_*z2_ - u2mod*y1_*z3_ + u1mod*y2_*z3_)*invdet;
			real ny = (-(u3mod*x2_*z1_) + u2mod*x3_*z1_ + u3mod*x1_*z2_ - u1mod*x3_*z2_ - u2mod*x1_*z3_ + u1mod*x2_*z3_)*invdet;
			real nz = -(-(u3mod*x2_*y1_) + u2mod*x3_*y1_ + u3mod*x1_*y2_ - u1mod*x3_*y2_ - u2mod*x1_*y3_ + u1mod*x2_*y3_)*invdet;

			normalize(nx,ny,nz);

			glNormal3f(nx,ny,nz);
		}

		if (ind == 3) {
			glBegin(GL_TRIANGLES);
			vert(0); vert(1); vert(2);
			glEnd();
		} else if (ind == 4) {
			glBegin(GL_QUADS);
			if (in_13 && in_14 && in_23 && in_24) {
				// 13, 14, 24, 23
				vert(0); vert(1); vert(3); vert(2);
			} else if (in_12 && in_14 && in_23 && in_34) {
				// 12, 14, 34, 23
//				vert(1); vert(0); vert(2); vert(3);
				vert(0); vert(1); vert(3); vert(2);
			} else if (in_12 && in_13 && in_24 && in_34) {
				// 12, 13, 34, 24
//				vert(0); vert(2); vert(3); vert(1);
				vert(0); vert(1); vert(3); vert(2);
			} else {
				cerr << "*** isosurf: impossible case\n";
			}
			glEnd();
		} else {
			cerr << "*** isosurf warning: ind=" << ind << "\n";
		}
#		undef vert
	}
}

void GLIsosurfCube_nosubdiv(real xc, real yc, real zc,
							real dx,
							const real u[8],
							const TContourSpec& cs,
							void (*ComputeGradient)(real,real,real, real&,real&,real&))
//
//         6  --------- 7
//       / |          / |
//      /  |         /  |
//    4  --|------  5   |
//    |    |        |   |
//    |    2 -------|-  3
//    |  /          |  /
//    | /           | /
//    0  ---------  1
{
	// Fast return if none of the isosurfs in cs applies
	real umin=u[0], umax=umin;
	int a,i;
	for (a=1; a<8; a++) {
		umin = min(umin,u[a]);
		umax = max(umax,u[a]);
	}
	int found = 0;
	for (i=0; i<cs.Ncontours(); i++) {
		const real u = cs.nthcontour(i);
		if (umin < u && u < umax) {found=1; break;}
	}
	if (!found) return;
	const real halfdx = 0.5*dx;
	const real x0 = xc - halfdx;
	const real x1 = xc + halfdx;
	const real x2 = x0;
	const real x3 = x1;
	const real x4 = x0;
	const real x5 = x1;
	const real x6 = x0;
	const real x7 = x1;
	const real y0 = yc - halfdx;
	const real y1 = y0;
	const real y2 = yc + halfdx;
	const real y3 = y2;
	const real y4 = y0;
	const real y5 = y0;
	const real y6 = y2;
	const real y7 = y2;
	const real z0 = zc - halfdx;
	const real z1 = z0;
	const real z2 = z0;
	const real z3 = z0;
	const real z4 = zc + halfdx;
	const real z5 = z4;
	const real z6 = z4;
	const real z7 = z4;
	// fc0,1,2,3,4,5: face centers
	const real xfc0 = x0;
	const real yfc0 = yc;
	const real zfc0 = zc;
	const real xfc1 = x1;
	const real yfc1 = yc;
	const real zfc1 = zc;
	const real xfc2 = xc;
	const real yfc2 = y0;
	const real zfc2 = zc;
	const real xfc3 = xc;
	const real yfc3 = y2;
	const real zfc3 = zc;
	const real xfc4 = xc;
	const real yfc4 = yc;
	const real zfc4 = z0;
	const real xfc5 = xc;
	const real yfc5 = yc;
	const real zfc5 = z4;
	// u values:
	const real u0 = u[0];
	const real u1 = u[1];
	const real u2 = u[2];
	const real u3 = u[3];
	const real u4 = u[4];
	const real u5 = u[5];
	const real u6 = u[6];
	const real u7 = u[7];
	const real ufc0 = 0.25*(u0+u2+u6+u4);
	const real ufc1 = 0.25*(u1+u3+u7+u5);
	const real ufc2 = 0.25*(u0+u1+u5+u4);
	const real ufc3 = 0.25*(u2+u3+u7+u6);
	const real ufc4 = 0.25*(u0+u1+u3+u2);
	const real ufc5 = 0.25*(u4+u5+u7+u6);
	const real uc = (ufc0+ufc1+ufc2+ufc3+ufc4+ufc5)/6.0;
#	define tetra(A,B,C,D) \
	GLIsosurfTetrahedron(x##A,y##A,z##A, x##B,y##B,z##B, x##C,y##C,z##C, x##D,y##D,z##D,\
						 u##A,u##B,u##C,u##D,cs,ComputeGradient)
	// 0: 0264
	tetra(0,2,fc0,c);
	tetra(2,6,fc0,c);
	tetra(6,4,fc0,c);
	tetra(4,0,fc0,c);
	// 1: 1375
	tetra(1,3,fc1,c);
	tetra(3,7,fc1,c);
	tetra(7,5,fc1,c);
	tetra(5,1,fc1,c);
	// 2: 0154
	tetra(0,1,fc2,c);
	tetra(1,5,fc2,c);
	tetra(5,4,fc2,c);
	tetra(4,0,fc2,c);
	// 3: 2376
	tetra(2,3,fc3,c);
	tetra(3,7,fc3,c);
	tetra(7,6,fc3,c);
	tetra(6,2,fc3,c);
	// 4: 0132
	tetra(0,1,fc4,c);
	tetra(1,3,fc4,c);
	tetra(3,2,fc4,c);
	tetra(2,0,fc4,c);
	// 5: 4576
	tetra(4,5,fc5,c);
	tetra(5,7,fc5,c);
	tetra(7,6,fc5,c);
	tetra(6,4,fc5,c);
#	undef tetra
}

void GLIsosurfCube(real xc, real yc, real zc,
				   real dx,
				   const real u[8],
				   const bool udens[6], const real ufacecenters[6],
				   const TContourSpec& cs,
				   real (*Interpolate)(real,real,real),
				   void (*ComputeGradient)(real,real,real, real&,real&,real&))
//
//         6  --------- 7
//       / |          / |
//      /  |         /  |
//    4  --|------  5   |
//    |    |        |   |
//    |    2 -------|-  3
//    |  /          |  /
//    | /           | /
//    0  ---------  1
{
	// Fast return if none of the isosurfs in cs applies
	real umin=u[0], umax=umin;
	int a,i;
	for (a=1; a<8; a++) {
		umin = min(umin,u[a]);
		umax = max(umax,u[a]);
	}
	int found = 0;
	for (i=0; i<cs.Ncontours(); i++) {
		const real u = cs.nthcontour(i);
		if (umin < u && u < umax) {found=1; break;}
	}
	if (!found) return;
	
	if (udens[0] || udens[1] || udens[2] || udens[3] || udens[4] || udens[5]) {
		int dirx,diry,dirz;
		const real qdx = 0.25*dx;
		const real hdx = 0.5*dx;
		real u1[8];
		real utab[3][3][3];
		utab[0][0][0] = u[0];
		utab[2][0][0] = u[1];
		utab[0][2][0] = u[2];
		utab[2][2][0] = u[3];
		utab[0][0][2] = u[4];
		utab[2][0][2] = u[5];
		utab[0][2][2] = u[6];
		utab[2][2][2] = u[7];
		utab[0][1][1] = udens[0] ? ufacecenters[0] : 0.25*(u[0]+u[2]+u[4]+u[6]);
		utab[2][1][1] = udens[1] ? ufacecenters[1] : 0.25*(u[1]+u[3]+u[5]+u[7]);
		utab[1][0][1] = udens[2] ? ufacecenters[2] : 0.25*(u[0]+u[1]+u[4]+u[5]);
		utab[1][2][1] = udens[3] ? ufacecenters[3] : 0.25*(u[2]+u[3]+u[6]+u[7]);
		utab[1][1][0] = udens[4] ? ufacecenters[4] : 0.25*(u[0]+u[1]+u[2]+u[3]);
		utab[1][1][2] = udens[5] ? ufacecenters[5] : 0.25*(u[4]+u[5]+u[6]+u[7]);
		// x=xmin:
		utab[0][0][1] = (udens[0] || udens[2]) ? (*Interpolate)(xc-hdx,yc-hdx,zc    ) : 0.5*(u[0] + u[4]);
		utab[0][2][1] = (udens[0] || udens[3]) ? (*Interpolate)(xc-hdx,yc+hdx,zc    ) : 0.5*(u[2] + u[6]);
		utab[0][1][0] = (udens[0] || udens[4]) ? (*Interpolate)(xc-hdx,yc,    zc-hdx) : 0.5*(u[0] + u[2]);
		utab[0][1][2] = (udens[0] || udens[5]) ? (*Interpolate)(xc-hdx,yc,    zc+hdx) : 0.5*(u[4] + u[6]);
		// x=xmax:
		utab[2][0][1] = (udens[1] || udens[2]) ? (*Interpolate)(xc+hdx,yc-hdx,zc    ) : 0.5*(u[1] + u[5]);
		utab[2][2][1] = (udens[1] || udens[3]) ? (*Interpolate)(xc+hdx,yc+hdx,zc    ) : 0.5*(u[3] + u[7]);
		utab[2][1][0] = (udens[1] || udens[4]) ? (*Interpolate)(xc+hdx,yc,    zc-hdx) : 0.5*(u[1] + u[3]);
		utab[2][1][2] = (udens[1] || udens[5]) ? (*Interpolate)(xc+hdx,yc,    zc+hdx) : 0.5*(u[5] + u[7]);
		// others:
		utab[1][0][0] = (udens[2] || udens[4]) ? (*Interpolate)(xc,yc-hdx,zc-hdx) : 0.5*(u[0] + u[1]);
		utab[1][2][0] = (udens[3] || udens[4]) ? (*Interpolate)(xc,yc+hdx,zc-hdx) : 0.5*(u[2] + u[3]);
		utab[1][0][2] = (udens[2] || udens[5]) ? (*Interpolate)(xc,yc-hdx,zc+hdx) : 0.5*(u[4] + u[5]);
		utab[1][2][2] = (udens[3] || udens[5]) ? (*Interpolate)(xc,yc+hdx,zc+hdx) : 0.5*(u[6] + u[7]);
		// centroid:
		utab[1][1][1] = (utab[0][1][1] + utab[2][1][1] + utab[1][0][1] + utab[1][2][1] + utab[1][1][0] + utab[1][1][2])/6.0;

		for (dirx=-1; dirx<=1; dirx+=2) for (diry=-1; diry<=1; diry+=2) for (dirz=-1; dirz<=1; dirz+=2) {
			const int i=(1+dirx)/2, j=(1+diry)/2, k=(1+dirz)/2;
			u1[0] = utab[i][j][k];
			u1[1] = utab[i+1][j][k];
			u1[2] = utab[i][j+1][k];
			u1[3] = utab[i+1][j+1][k];
			u1[4] = utab[i][j][k+1];
			u1[5] = utab[i+1][j][k+1];
			u1[6] = utab[i][j+1][k+1];
			u1[7] = utab[i+1][j+1][k+1];
			GLIsosurfCube_nosubdiv(xc+dirx*qdx,yc+diry*qdx,zc+dirz*qdx,hdx,u1,cs,ComputeGradient);
		}
	} else {
		GLIsosurfCube_nosubdiv(xc,yc,zc,dx,u,cs,ComputeGradient);
	}
}
