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

#ifndef R3V_H
#define R3V_H

#include <cmath>
#include <iostream>
#include "definitions.h"
#include "random.h"

//! 3-vector of gridreal components for detectors - time critical
class Tgr3v
{
public:
    gridreal r[3];
    Tgr3v() {}
    Tgr3v(gridreal x, gridreal y, gridreal z) {
        r[0] = x;
        r[1] = y;
        r[2] = z;
    }
    Tgr3v(const gridreal rr[3]) {
        r[0] = rr[0];
        r[1] = rr[1];
        r[2] = rr[2];
    }
    gridreal& operator[](int k) {
        return r[k];
    }
    real operator()(int k) const {
        return r[k];
    }
    real magn2() const {
        return sqr(r[0]) + sqr(r[1]) + sqr(r[2]);
    }
    real magn() const {
        return sqrt(magn2());
    }
    Tgr3v operator*(real a) const {
        return Tgr3v(r[0]*a,r[1]*a,r[2]*a);
    }
    Tgr3v operator+(const Tgr3v& b) const {
        return Tgr3v(r[0]+b(0), r[1]+b(1), r[2]+b(2));
    }
    Tgr3v operator-(const Tgr3v& b) const {
        return Tgr3v(r[0]-b(0), r[1]-b(1), r[2]-b(2));
    }
    Tgr3v operator/(real a) const {
        const real inva=1.0/a;
        return Tgr3v(r[0]*inva,r[1]*inva,r[2]*inva);
    }
    Tgr3v& operator-=(const Tgr3v& b) {
        r[0]-= b(0);
        r[1]-= b(1);
        r[2]-= b(2);
        return *this;
    }
    Tgr3v operator-() const {
        return Tgr3v(-r[0],-r[1],-r[2]);
    }
    friend Tgr3v UnitVector(const Tgr3v& v) {
        const gridreal vmagn = v.magn();
        if (vmagn == 0) return Tgr3v(1,0,0);
        const gridreal inv_vmagn = 1.0/vmagn;
        Tgr3v u;
        for (int d=0; d<3; d++) u[d] = v(d)*inv_vmagn;
        return u;
    }
    friend Tgr3v Cross(const Tgr3v& a, const Tgr3v& b) {
        Tgr3v c;
        c[0] = a(1)*b(2) - a(2)*b(1);
        c[1] = a(2)*b(0) - a(0)*b(2);
        c[2] = a(0)*b(1) - a(1)*b(0);
        return c;
    }
    friend Tgr3v LinearlyIndependentVector(const Tgr3v& a) {
        Tgr3v u(0,0,0);
        if (fabs(a(0)) < fabs(a(1)))
            u[0] = 1;
        else
            u[1] = 1;
        return u;
    }
    friend Tgr3v RotateVector(const Tgr3v& a, real phi, const Tgr3v& ax); // Rotate vector a angle phi about axis vector ax
    friend Tgr3v UnitVector_robust(const Tgr3v& v, const Tgr3v& n);
    friend std::ostream& operator<<(std::ostream& o, const Tgr3v& x);
    std::string toString();
};

//! Multiply operator
inline Tgr3v operator*(fastreal a, const Tgr3v& b)
{
    return Tgr3v(a*b(0),a*b(1),a*b(2));
}

//! Dot product
inline fastreal dot(const Tgr3v& a, const Tgr3v& b)
{
    return a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
}

//! Magnitude squared
inline fastreal magn2(const Tgr3v a)
{
    return a.magn2();
}

//! Random perpendicular unit vector
inline Tgr3v RandomPerpendicularUnitVector(const Tgr3v& u)
{
    return UnitVector(RotateVector(Cross(u,LinearlyIndependentVector(u)),2*pi*uniformrnd(),u));
}

//! 3-vector of real components
class Tr3v
{
private:
    real r[3];
public:
    Tr3v() {}
    Tr3v(real x, real y, real z) {
        r[0] = x;
        r[1] = y;
        r[2] = z;
    }
    Tr3v(const gridreal rr[3]) {
        r[0] = rr[0];
        r[1] = rr[1];
        r[2] = rr[2];
    }
    real& operator[](int k) {
        return r[k];
    }
    real operator()(int k) const {
        return r[k];
    }
    real magn2() const {
        return sqr(r[0]) + sqr(r[1]) + sqr(r[2]);
    }
    real magn() const {
        return sqrt(magn2());
    }
    Tr3v operator*(real a) const {
        return Tr3v(r[0]*a,r[1]*a,r[2]*a);
    }
    Tr3v operator+(const Tr3v& b) const {
        return Tr3v(r[0]+b(0), r[1]+b(1), r[2]+b(2));
    }
    Tr3v operator-(const Tr3v& b) const {
        return Tr3v(r[0]-b(0), r[1]-b(1), r[2]-b(2));
    }
    Tr3v operator/(real a) const {
        const real inva=1.0/a;
        return Tr3v(r[0]*inva,r[1]*inva,r[2]*inva);
    }
    Tr3v& operator-=(const Tr3v& b) {
        r[0]-= b(0);
        r[1]-= b(1);
        r[2]-= b(2);
        return *this;
    }
    Tr3v operator-() const {
        return Tr3v(-r[0],-r[1],-r[2]);
    }
    friend Tr3v UnitVector(const Tr3v& v) {
        const real vmagn = v.magn();
        if (vmagn == 0) return Tr3v(1,0,0);
        const real inv_vmagn = 1.0/vmagn;
        int d;
        Tr3v u;
        for (d=0; d<3; d++) u[d] = v(d)*inv_vmagn;
        return u;
    }
    friend Tr3v Cross(const Tr3v& a, const Tr3v& b) {
        Tr3v c;
        c[0] = a(1)*b(2) - a(2)*b(1);
        c[1] = a(2)*b(0) - a(0)*b(2);
        c[2] = a(0)*b(1) - a(1)*b(0);
        return c;
    }
    friend Tr3v LinearlyIndependentVector(const Tr3v& a) {
        Tr3v u(0,0,0);
        if (fabs(a(0)) < fabs(a(1)))
            u[0] = 1;
        else
            u[1] = 1;
        return u;
    }
    friend Tr3v RotateVector(const Tr3v& a, real phi, const Tr3v& ax); // Rotate vector a angle phi about axis vector ax
    friend Tr3v UnitVector_robust(const Tr3v& v, const Tr3v& n);
    friend std::ostream& operator<<(std::ostream& o, const Tr3v& x);
    std::string toString();
};

//! Multiply operator
inline Tr3v operator*(real a, const Tr3v& b)
{
    return Tr3v(a*b(0),a*b(1),a*b(2));
}

//! Dot product
inline real dot(const Tr3v& a, const Tr3v& b)
{
    return a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
}

//! Magnitude squared
inline real magn2(const Tr3v a)
{
    return a.magn2();
}

//! Random perpendicular unit vector
inline Tr3v RandomPerpendicularUnitVector(const Tr3v& u)
{
    return UnitVector(RotateVector(Cross(u,LinearlyIndependentVector(u)),2*pi*uniformrnd(),u));
}

#endif

