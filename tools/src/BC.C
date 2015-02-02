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
#  pragma implementation "BC.H"
#endif
#include "BC.H"
#include <math.h>
#include <stdlib.h>


real Tdimvec::OutputXScaling = 1;

ostream& operator<<(ostream& o, const Tdimvec& u) {
	unsigned short a;
	o << "#(";
	for (a=0; a<MAXDIM-1; a++) o << u(a)/Tdimvec::OutputXScaling << ',';
	return (o << u(MAXDIM-1)/Tdimvec::OutputXScaling << ')');
}

Tdimvec Tdimvec::operator+(const Tdimvec& b) const
{
	unsigned short d;
	Tdimvec result;
	for (d=0; d<MAXDIM; d++) result[d] = u[d] + b(d);
	return result;
}

Tdimvec Tdimvec::operator*(const Tdimvec& b) const
{
	unsigned short d;
	Tdimvec result;
	for (d=0; d<MAXDIM; d++) result[d] = u[d] * b(d);
	return result;
}

Tdimvec& Tdimvec::operator+=(const Tdimvec& b)
{
	unsigned short d;
	for (d=0; d<MAXDIM; d++) u[d]+= b(d);
	return *this;
}

Tdimvec operator*(real x, const Tdimvec& a)
{
	unsigned short d;
	Tdimvec result = 0;
	for (d=0; d<MAXDIM; d++) result[d] = x*a(d);
	return result;
}

bool TBoundaryGeometry::BarIntersection
    (const Tdimvec& /*r1*/, const Tdimvec& /*r2*/, Tdimvec& /*result*/) const
{
	cerr << "*** TBoundaryGeometry::BarIntersection called\n";
	cerr << "    This should never happen. You have derived from TBoundaryGeometry,\n";
	cerr << "    defined CanClip()=1 and forgotten to define BarIntersection()\n";
	exit(1);
	return 0;
}

bool TPlaneBoundaryGeometry::IsInsideBody(const Tdimvec& point) const
{
	const real x = point(d) - xval;
	const int retval = (neg_x_in_sim_domain ? (x >= 0) : (x < 0));
	return retval;
}

bool TInclinedPlaneBoundaryGeometry::IsInsideBody(const Tdimvec& point) const
{
	// return normal.(point - onepoint) >= 0
	unsigned short d;
	real dr, dotprod=0;
	for (d=0; d<dim; d++) {
		dr = point(d) - onepoint(d);
		dotprod+= normal(d)*dr;
	}
	const int retvalue = (dotprod >= 0);
//	cout << "incl.plane.IsInsideBody(" << point << ") = " << retvalue << "\n";
	return retvalue;
}

bool TSphericalBoundaryGeometry::IsInsideBody(const Tdimvec& point) const
{
	unsigned short d;
	real radSquared = 0;
	for (d=0; d<dim; d++) radSquared+= sqr(point(d) - center(d));
	const int isinside = (radSquared <= rSquared);
	return isinside;
}

bool TEllipsoidalBoundaryGeometry::IsInsideBody(const Tdimvec& point) const
{
	unsigned short d;
	real radSquared = 0;
	for (d=0; d<dim; d++) radSquared+= sqr((point(d) - center(d))/r(d));
	const int isinside = (radSquared <= 1);
	return isinside;
}

bool TRectangularBoxBoundaryGeometry::IsInsideBody(const Tdimvec& point) const
{
	unsigned short d;
	int isinside = 1;
	for (d=0; d<dim; d++) {
		if (!(xmin(d) <= point(d) && point(d) <= xmax(d))) {
			isinside = 0;
			break;
		}
	}
	return isinside;
}

bool TInclinedPlaneBoundaryGeometry::BarIntersection
    (const Tdimvec& r1, const Tdimvec& r2, Tdimvec& result) const
{
	// The equation of the plane is: normal.(result - onepoint) == 0
	// The bar is parameterized by:  result = r1 + t*u, where u = r2-r1
	unsigned short d;
	Tdimvec u, dr;
	for (d=0; d<dim; d++) {
		u[d] = r2(d) - r1(d);
		dr[d] = onepoint(d) - r1(d);
	}
	real n_dot_u = 0, n_dot_dr = 0;
	for (d=0; d<dim; d++) {
		n_dot_u+= normal(d)*u(d);
		n_dot_dr+= normal(d)*dr(d);
	}
	if (n_dot_u == 0) return 0;
	const real t = n_dot_dr/n_dot_u;
	if (!(0 <= t && t <= 1)) return 0;
	for (d=0; d<dim; d++)
		result[d] = r1(d) + t*u(d);
	return 1;
}

bool TSphericalBoundaryGeometry::BarIntersection
    (const Tdimvec& r1, const Tdimvec& r2, Tdimvec& result) const
{
	if (dim != 2 && dim != 3) {
		cerr << "*** TSphericalBoundaryGeometry::BarIntersection in dim=" << dim << " not implemented\n";
		return 0;
	}
	unsigned short d;
	Tdimvec u;
	for (d=0; d<dim; d++) u[d] = r2(d) - r1(d);
	real u2 = 0;
	for (d=0; d<dim; d++) u2+= sqr(u(d));
	const real a = u2;
	real b = 0;
	for (d=0; d<dim; d++) b+= 2*u(d)*(r1(d) - center(d));
	real r1_dot_r0 = 0;
	for (d=0; d<dim; d++) r1_dot_r0+= r1(d)*center(d);
	real r12 = 0, r02 = 0;
	for (d=0; d<dim; d++) {
		r12+= sqr(r1(d));
		r02+= sqr(center(d));
	}
	const real c = r12 + r02 - rSquared - 2*r1_dot_r0;
	const real D = sqr(b) - 4*a*c;
	// Now we want to solve a 2nd order algebraic equation a*t^2 + b*t + c == 0.
	// See Numerical Recipes p. 145. The discriminant is already computed in D.
	if (D < 0) return 0;	// If no real roots, there is no intersection.
	const int sgnb = (b >= 0 ? +1 : -1);
	const real q = -0.5*(b + sgnb*sqrt(D));
	// If division by zero would result, there is certainly no intersection.
	if (a == 0 || q == 0) return 0;
	// Now compute the two roots t1 and t2. Eqs. 5.5.4-5.5.5 in Numerical Recipes 1st ed.
	const real t1 = q/a;
	const real t2 = c/q;
	// Valid roots are in the range 0..1. If both roots are valid then the sphere
	// intersects the bar twice. We judge this case as a non-intersection, however.
	// (This is numerically sensible, if mathematically weird).
	if (0 <= t1 && t1 <= 1) {
		if (0 <= t2 && t2 <= 1) return 0;	// both valid ==> both invalid
		// Now only t1 is valid.
		for (d=0; d<dim; d++) result[d] = r1(d) + t1*u(d);
	} else {
		// t1 is invalid. Check is t2 is valid.
		if (!(0 <= t2 && t2 <= 1)) return 0;	// It was invalid
		// Now only t2 is valid.
		for (d=0; d<dim; d++) result[d] = r1(d) + t2*u(d);
	}
	return 1;
}

bool TEllipsoidalBoundaryGeometry::BarIntersection
     (const Tdimvec& r1, const Tdimvec& r2, Tdimvec& result) const
{
	if (dim != 2 && dim != 3) {
		cerr << "*** TEllipsoidalBoundaryGeometry::BarIntersection in dim=" << dim << " not implemented\n";
		return 0;
	}
	unsigned short d;
	Tdimvec u;
	for (d=0; d<dim; d++) u[d] = r2(d) - r1(d);
	real a = 0;
	for (d=0; d<dim; d++) a+= sqr(u(d)/r(d));
	real b = 0;
	for (d=0; d<dim; d++) b+= (u(d)/rSquared(d))*(r1(d) - center(d));
	b*= 2;
	real c = 0;
	for (d=0; d<dim; d++) c+= sqr((center(d) - r1(d))/r(d));
	c = -1 + c;
	/* The rest of the function is common with TSphericalBoundaryGeometry; it is the 2nd order solution */
	const real D = sqr(b) - 4*a*c;
	// Now we want to solve a 2nd order algebraic equation a*t^2 + b*t + c == 0.
	// See Numerical Recipes p. 145. The discriminant is already computed in D.
	if (D < 0) return 0;	// If no real roots, there is no intersection.
	const int sgnb = (b >= 0 ? +1 : -1);
	const real q = -0.5*(b + sgnb*sqrt(D));
	// If division by zero would result, there is certainly no intersection.
	if (a == 0 || q == 0) return 0;
	// Now compute the two roots t1 and t2. Eqs. 5.5.4-5.5.5 in Numerical Recipes 1st ed.
	const real t1 = q/a;
	const real t2 = c/q;
	// Valid roots are in the range 0..1. If both roots are valid then the sphere
	// intersects the bar twice. We judge this case as a non-intersection, however.
	// (This is numerically sensible, if mathematically weird).
	if (0 <= t1 && t1 <= 1) {
		if (0 <= t2 && t2 <= 1) return 0;	// both valid ==> both invalid
		// Now only t1 is valid.
		for (d=0; d<dim; d++) result[d] = r1(d) + t1*u(d);
	} else {
		// t1 is invalid. Check is t2 is valid.
		if (!(0 <= t2 && t2 <= 1)) return 0;	// It was invalid
		// Now only t2 is valid.
		for (d=0; d<dim; d++) result[d] = r1(d) + t2*u(d);
	}
	return 1;
}

void TPlaneBoundaryGeometry::Normal(const Tdimvec& /*point*/, real& nx, real& ny, real& nz) const
{
	if (neg_x_in_sim_domain) {
		if (d == 0) {
			nx = -1;
			ny = 0;
			nz = 0;
		} else if (d == 1) {
			nx = 0;
			ny = -1;
			nz = 0;
		} else {
			nx = 0;
			ny = 0;
			nz = -1;
		}
	} else {
		if (d == 0) {
			nx = 1;
			ny = 0;
			nz = 0;
		} else if (d == 1) {
			nx = 0;
			ny = 1;
			nz = 0;
		} else {
			nx = 0;
			ny = 0;
			nz = 1;
		}
	}
}

void TInclinedPlaneBoundaryGeometry::Normal(const Tdimvec& /*point*/, real& nx, real& ny, real& nz) const
{
	// 'normal' points from domain into body, but (nx,ny,nz) we want to point from body to domain
	nx = -normal(0);
	ny = -normal(1);
	nz = -normal(2);
}

void TSphericalBoundaryGeometry::Normal(const Tdimvec& point, real& nx, real& ny, real& nz) const
{
	nx = point(0) - center(0);
	ny = point(1) - center(1);
	if (dim >= 3)
		nz = point(2) - center(2);
	else
		nz = 0;
	const real C = 1/sqrt(sqr(nx) + sqr(ny) + sqr(nz));
	nx*= C; ny*= C; nz*= C;
}

void TEllipsoidalBoundaryGeometry::Normal(const Tdimvec& point, real& nx, real& ny, real& nz) const
{
	nx = (point(0) - center(0))/rSquared(0);
	ny = (point(1) - center(1))/rSquared(1);
	if (dim >= 3)
		nz = (point(2) - center(2))/rSquared(2);
	else
		nz = 0;
	const real C = 1/sqrt(sqr(nx) + sqr(ny) + sqr(nz));
	nx*= C; ny*= C; nz*= C;
}

void TRectangularBoxBoundaryGeometry::Normal(const Tdimvec& point, real& nx, real& ny, real& nz) const
{
	if (dim == 2) {
		// (cx,cy) is the box center
		const real cx = 0.5*(xmin(0) + xmax(0));
		const real cy = 0.5*(xmin(1) + xmax(1));
		// V1,V2,V3,V4 will be diff. vector between corners and center
		const real V1x = xmax(0) - cx;
		const real V1y = xmax(1) - cy;
		const real V2x = xmin(0) - cx;
		const real V2y = xmax(1) - cy;
		const real V4x = xmax(0) - cx;
		const real V4y = xmin(1) - cy;
		// V = point - c
		const real Vx = point(0) - cx;
		const real Vy = point(1) - cy;
		const int pos1 = (V1x*Vy - V1y*Vx >= 0);
		const int pos2 = (V2x*Vy - V2y*Vx >= 0);
		const int pos4 = (V4x*Vy - V4y*Vx >= 0);
		nz = 0;
		if (pos1) {
			// case is 1 or 2
			if (pos2) {
				// Case 2
				nx = -1;
				ny = 0;
			} else {
				// Case 1
				nx = 0;
				ny = 1;
			}
		} else {
			// case is 3 or 4
			if (pos4) {
				// Case 4
				nx = 1;
				ny = 0;
			} else {
				// Case 3
				nx = 0;
				ny = -1;
			}
		}
	} else if (dim == 3) {
		cerr << "*** TRectangularBoxBoundaryGeometry::Normal() not yet implemented in 3D\n";
		nx=1; ny=0; nz=0;
	} else {
		cerr << "*** TRectangularBoxBoundaryGeometry::Normal() called for dim=" << dim << "\n";
		nx=1; ny=0; nz=0;
	}
}




