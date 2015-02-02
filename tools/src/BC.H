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

#ifndef BC_H

#ifdef __GNUC__
#  pragma interface
#endif

#include "gengrid.H"
#include <math.h>

enum TBoundaryGeometryType {
	PLANAR_BOUNDARY, INCLINED_PLANAR_BOUNDARY,
	RECTANGULAR_BOX_BOUNDARY, SPHERICAL_BOUNDARY, ELLIPSOIDAL_BOUNDARY
};

// A small utility class for storing MAXDIM-vector
class Tdimvec {
private:
	real u[MAXDIM];
	void copyfrom(const Tdimvec& v) {
#if USE_MEMCPY
		memcpy(u,v.u,MAXDIM*sizeof(real));
#else
		u[0] = v.u[0];
#	if MAXDIM >= 2
		u[1] = v.u[1];
#		if MAXDIM == 3
		u[2] = v.u[2];
#		endif
#	endif
#endif
	}
public:
	Tdimvec() {}
	Tdimvec(real x) {							// assign the same value for all elements
		smallnat d;
		for (d=0; d<MAXDIM; d++) u[d] = x;
	}
#	if MAXDIM==2
	Tdimvec(real x, real y) {u[0] = x; u[1] = y;}
#	elif MAXDIM==3
	Tdimvec(real x, real y) {u[0] = x; u[1] = y; u[2] = 0;}
	Tdimvec(real x, real y, real z) {u[0] = x; u[1] = y; u[2] = z;}
#	endif
	static real OutputXScaling;						// Outputted cell coords are divided by this, default 1.
	Tdimvec(const Tdimvec& v) {copyfrom(v);}		// copy constructor
	Tdimvec& operator=(const Tdimvec& v) {copyfrom(v); return *this;}
	real& operator[](int a) {return u[a];}			// indexing, used on LHS of assignment stmt
	real operator()(int a) const {return u[a];}	// indexing, used on RHS of assignment stmt
	Tdimvec operator+(const Tdimvec& b) const;
	Tdimvec operator*(const Tdimvec& b) const;
	Tdimvec& operator+=(const Tdimvec& v);
	void normalize() {
		real norm=0; smallnat d;
		for (d=0; d<MAXDIM; d++) norm+= sqr(u[d]);
		norm = 1.0/sqrt(norm);
		for (d=0; d<MAXDIM; d++) u[d]*= norm;
	}
	friend ostream& operator<<(ostream& o, const Tdimvec& u);
	~Tdimvec() {}
};

extern Tdimvec operator*(real x, const Tdimvec& a);

extern ostream& operator<<(ostream& o, const Tdimvec& u);

/*
  TBoundaryGeometry is abstract base class (ABC) for an analytically defined boundary shape.
  The concrete classes must define:
  - at least one constructor
  - IsInsideBody() for testing whether a given point is outside the simulation domain (i.e., "inside the body")
  - type() which returns one of the TBoundaryGeometryType's above
  - CanClip() which tells whether BarIntersection() is defined or not
  - BarIntersection() [if CanClip()=1] for finding the intersection point of the boundary surface
    and a given bar (bar is a line between two points)
  - Normal() which returns a normal vector (from body into domain)
  - copy(). Just copy it from TPlaneBoundaryGeometry, for example.
*/

class TBoundaryGeometry {
protected:
	smallnat dim;			// dimension of the grid (1,2,3)
public:
	TBoundaryGeometry(smallnat dim1) {dim=dim1;}
//	int matches(const Tcell *c, TCellType *detailed_info_ptr=0) const;
	virtual bool IsInsideBody(const Tdimvec&) const = 0;
	// Test whether a single cell is at the boundary or not.
	virtual bool BarIntersection(const Tdimvec& r1, const Tdimvec& r2, Tdimvec& result) const;
	// find intersection of line from r1 to r2 and the boundary. Return true if does intersect, and return intersection point in r.
	// If the derived class has CanClip()==0, BarIntersection() need not be redefined.
	virtual TBoundaryGeometryType type() const = 0;
	virtual bool CanClip() const = 0;				// whether the boundary defines BarIntersection or not
	virtual void Normal(const Tdimvec& point, real& nx, real &ny, real& nz) const = 0;
	virtual TBoundaryGeometry* copy() const = 0;	// Just copy the TPlaneBoundaryGeometry definition to new classes.
	virtual ~TBoundaryGeometry() {}
};

// A planar boundary is a concrete example of abstract boundary.
class TPlaneBoundaryGeometry : public TBoundaryGeometry {
private:
	smallnat d;					// dimension which is perpendicular to the plane (3D: 0,1,2)
	real xval;					// the coordinate value defining the plane (X[d]=xval)
	bool neg_x_in_sim_domain;	// Boolean flag: true if x<xval[d] are inside simulation domain
public:
	// The constructor takes the number of dimensions (dim1), perpendicular dimension d (0,1,2)
	// and the coordinate value, as well as neg_x_in_sim_domain flag
	TPlaneBoundaryGeometry(smallnat dim1, smallnat d1, real xval1, bool flag) : TBoundaryGeometry(dim1)
	{
		d = d1;
		xval = xval1;
		neg_x_in_sim_domain = flag;
	}
	virtual bool IsInsideBody(const Tdimvec& point) const;
	TBoundaryGeometryType type() const {return PLANAR_BOUNDARY;}
	void Normal(const Tdimvec& point, real& nx, real& ny, real& nz) const;
	TBoundaryGeometry* copy() const {return new TPlaneBoundaryGeometry(*this);}
	bool CanClip() const {return false;}
	virtual ~TPlaneBoundaryGeometry() {}
};

// An inclined plane is a generalization of TPlaneBoundaryGeometry.
class TInclinedPlaneBoundaryGeometry : public TBoundaryGeometry {
private:
	Tdimvec normal;			// vector perpendicular to the plane
	Tdimvec onepoint;		// one point of the plane
	// Points r for which normal.(r - onepoint) >= 0 are "inside body" (i.e. not in simulation domain)
	// i.e., the normal points from the domain towards the body
public:
	TInclinedPlaneBoundaryGeometry(smallnat dim1, const Tdimvec& normal1, const Tdimvec& onepoint1)
		: TBoundaryGeometry(dim1), normal(normal1), onepoint(onepoint1) {}
	bool IsInsideBody(const Tdimvec& point) const;
	bool BarIntersection(const Tdimvec& r1, const Tdimvec& r2, Tdimvec& result) const;
	TBoundaryGeometryType type() const {return INCLINED_PLANAR_BOUNDARY;}
	void Normal(const Tdimvec& point, real& nx, real& ny, real& nz) const;
	TBoundaryGeometry* copy() const {return new TInclinedPlaneBoundaryGeometry(*this);}
	bool CanClip() const {return true;}
	virtual ~TInclinedPlaneBoundaryGeometry() {}
};

// A spherical body immersed in fluid
class TSphericalBoundaryGeometry : public TBoundaryGeometry {
private:
	real r;					// radius of the sphere
	real rSquared;			// sqr(r)
	Tdimvec center;			// center coordinates of the sphere
public:
	TSphericalBoundaryGeometry(real r1, smallnat dim1, const Tdimvec& center1=0)
		: TBoundaryGeometry(dim1), center(center1)
	{
		r = r1;
		rSquared = sqr(r);
	}
	bool IsInsideBody(const Tdimvec& point) const;
	bool BarIntersection(const Tdimvec& r1, const Tdimvec& r2, Tdimvec& result) const;
	TBoundaryGeometryType type() const {return SPHERICAL_BOUNDARY;}
	void Normal(const Tdimvec& point, real& nx, real& ny, real& nz) const;
	TBoundaryGeometry* copy() const {return new TSphericalBoundaryGeometry(*this);}
	bool CanClip() const {return true;}
	virtual ~TSphericalBoundaryGeometry() {}
};

// An coordinate-aligned ellipsoidal body immersed in fluid (generalization of TSphericalBoundaryGeometry)
class TEllipsoidalBoundaryGeometry : public TBoundaryGeometry {
private:
	Tdimvec r;			// the "radius" in different directions
	Tdimvec rSquared;
	Tdimvec center;		// center coordinates of the ellipsoid
public:
	TEllipsoidalBoundaryGeometry(const Tdimvec& r1, smallnat dim1, const Tdimvec& center1=0)
		: TBoundaryGeometry(dim1), r(r1), rSquared(r1*r1), center(center1) {}
	bool IsInsideBody(const Tdimvec& point) const;
	bool BarIntersection(const Tdimvec& r1, const Tdimvec& r2, Tdimvec& result) const;
	TBoundaryGeometryType type() const {return ELLIPSOIDAL_BOUNDARY;}
	void Normal(const Tdimvec& point, real& nx, real& ny, real& nz) const;
	TBoundaryGeometry* copy() const {return new TEllipsoidalBoundaryGeometry(*this);}
	bool CanClip() const {return true;}
	virtual ~TEllipsoidalBoundaryGeometry() {}
};

// A rectangular box immersed in fluid
class TRectangularBoxBoundaryGeometry : public TBoundaryGeometry {
private:
	Tdimvec xmin, xmax;
public:
	TRectangularBoxBoundaryGeometry(smallnat dim1, const Tdimvec& xmin1, const Tdimvec& xmax1)
		: TBoundaryGeometry(dim1), xmin(xmin1), xmax(xmax1) {}
	bool IsInsideBody(const Tdimvec& point) const;
	TBoundaryGeometryType type() const {return RECTANGULAR_BOX_BOUNDARY;}
	void Normal(const Tdimvec& point, real& nx, real& ny, real& nz) const;
	TBoundaryGeometry* copy() const {return new TRectangularBoxBoundaryGeometry(*this);}
	bool CanClip() const {return false;}
	virtual ~TRectangularBoxBoundaryGeometry() {}
};

enum TBCtype {DIRICHLET_BC, NEUMANN_BC, WALL_BC, CUSTOM_BC};

typedef void (*TEvaluatorFunctionPtr)(TGenericGrid& g, const TGridIndexVector& iv, real t, should_be_const real X[MAXDIM][VECLEN], smallnat n, const real *u);
// TEvaluatorFunctionPtr fp; (*fp)(g,i,t,X,n,0)
// should compute n ncd-tuples of dynamic variables, all at time t, for grid point i.
// The coordinates are passed in array X[MAXDIM][VECLEN].
// The function should write the results to g.CT(a,v), for a=0..ncd-1 and v=0..n-1.
// If the last argument (u) is non-null, it points to averaged data of nearest interior neighbours.
// In this case, n is always passed as 1.

#define BC_H

#endif
