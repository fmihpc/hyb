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

#include <sstream>
#include "vectors.h"
#include "params.h"

using namespace std;

//! Rotate vector a angle phi about axis vector ax
Tgr3v RotateVector(const Tgr3v& a, real phi, const Tgr3v& ax)
{
    const Tgr3v n = UnitVector(ax);
    const Tgr3v apar = n*dot(a,n);
    const Tgr3v b = Cross(n,a);
    // dot term to preserve components along axis
    Tgr3v r = (a-apar)*cos(phi) + b*sin(phi) + apar;
    return r;
}

//! Unit vector
Tgr3v UnitVector_robust(const Tgr3v& v, const Tgr3v& n)
{
    Tgr3v u;
    const real vmagn = v.magn();
    if (vmagn < 1e-7) {
        u = UnitVector(Cross(n,LinearlyIndependentVector(n)));
    } else {
        const real inv_vmagn = 1.0/vmagn;
        int d;
        for (d=0; d<3; d++) u[d] = v(d)*inv_vmagn;
    }
    return u;
}

//! Stream operator
ostream& operator<<(ostream& o, const Tgr3v& x)
{
    return o << '(' << x.r[0] << ',' << x.r[1] << ',' << x.r[2] << ')';
}

//! Returns string (x,y,z) in units of R_P
string Tgr3v::toString()
{
    stringstream ss;
    ss << "(";
    const int n = 3;
    for (int i = 0; i < n; ++i) {
        ss << r[i]/Params::R_P;
        if(i < n-1) {
            ss << ",";
        }
    }
    ss << ")";
    return ss.str();
}

//! Rotate vector a angle phi about axis vector ax
Tr3v RotateVector(const Tr3v& a, real phi, const Tr3v& ax)
{
    const Tr3v n = UnitVector(ax);
    const Tr3v apar = n*dot(a,n);
    const Tr3v b = Cross(n,a);
    // dot term to preserve components along axis
    Tr3v r = (a-apar)*cos(phi) + b*sin(phi) + apar;
    return r;
}

//! Unit vector
Tr3v UnitVector_robust(const Tr3v& v, const Tr3v& n)
{
    Tr3v u;
    const real vmagn = v.magn();
    if (vmagn < 1e-7) {
        u = UnitVector(Cross(n,LinearlyIndependentVector(n)));
    } else {
        const real inv_vmagn = 1.0/vmagn;
        int d;
        for (d=0; d<3; d++) u[d] = v(d)*inv_vmagn;
    }
    return u;
}

//! Stream operator
ostream& operator<<(ostream& o, const Tr3v& x)
{
    return o << '(' << x.r[0] << ',' << x.r[1] << ',' << x.r[2] << ')';
}

//! Returns string (x,y,z) in units of R_P
string Tr3v::toString()
{
    stringstream ss;
    ss << "(";
    const int n = 3;
    for (int i = 0; i < n; ++i) {
        ss << r[i]/Params::R_P;
        if(i < n-1) {
            ss << ",";
        }
    }
    ss << ")";
    return ss.str();
}

