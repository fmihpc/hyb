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

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

#include "definitions.h"
#include "params.h"

//! (SPHERICAL) Transformation: Position (x, y, z) -> (r, theta, phi)
inline void sph_transf_C2S_r(gridreal r[3])
{
    gridreal x = r[0];
    gridreal y = r[1];
    gridreal z = r[2];
    r[0] = sqrt(x*x + y*y + z*z);         // r = [r1, r2], r2 > r1 > 0
    r[1] = acos(z/sqrt(x*x + y*y + z*z)); // theta = [0, pi]
    r[2] = atan2(y,x);                    // phi =   [-pi, pi]
    //if(r[2] < 0)  r[2] = r[2] + 2*pi;
}

//! (SPHERICAL) Transformation: Position (x, y, z) -> (r, theta, phi)
inline void sph_transf_S2C_R(gridreal r[3])
{
    gridreal R     = r[0];
    gridreal theta = r[1];
    gridreal phi   = r[2];
    r[0] = R*sin(theta)*cos(phi);
    r[1] = R*sin(theta)*sin(phi);
    r[2] = R*cos(theta);
}

/** \brief (SPHERICAL) Transformation: Position Hybrid (x, y, z) -> (r, theta, phi)
 *
 * We obtain the formulas from the following assumptions (type 0 in conficuratin file):
 *
 *         y          theta - pi/2          z            phi
 *      ------- = ------------------- ,  ------- = --------------
 *      y0 + sph_dy   pi/2 - alpha_theta     z0 + sph_dz   pi - alpha_phi
 *
 * alpha_theta and alpha_phi are defined in .cfg file (normalized to pi).
 *
 * We obtain the formulas from the following assumptions (type 1 in conficuratin file):
 *
 *         y      theta - pi/2         z        phi
 *      ------- = -------------- , ------- = ---------
 *      y0 + sph_dy   pi/2 + dtheta    z0 + sph_dz   pi + dphi
 *
 * alpha_theta and alpha_phi are defined in .cfg file (normalized to pi).
 */
inline void sph_transf_H2S_R(gridreal r[3])
{
    gridreal x = r[0];
    gridreal y = r[1];
    gridreal z = r[2];
    if(Params::sph_BC_type == 0) {
        //r[1] = y/abs(Params::box_ymax + Params::sph_dy)*pi/2.0*(1.0 - 2.0*Params::sph_alpha_theta) + pi/2; // theta = [alpha_theta, pi - alpha_theta]
        //r[2] = z/abs(Params::box_zmax + Params::sph_dz)*pi*(1.0 - Params::sph_alpha_phi);                  // phi   = [-pi + alpha_phi, pi - alpha_phi]
        r[1] = y/abs(Params::box_ymax + Params::sph_dy)*(Params::sph_theta_max - pi/2) + pi/2;               // theta = [alpha_theta, pi - alpha_theta]
        r[2] = z/abs(Params::box_zmax + Params::sph_dz)*Params::sph_phi_max;                                 // phi   = [-pi + alpha_phi, pi - alpha_phi]
    } else if(Params::sph_BC_type == 1) {
        //r[1] = (y/abs(Params::box_ymax)+1.0)*pi/2.0; // theta = [0, pi]
        //r[2] = (z/abs(Params::box_zmax))*pi;         // phi   = [-pi, pi]
        //r[1] = y/abs(Params::box_ymax + Params::sph_dy)*(pi/2.0 + Params::sph_dtheta) + pi/2; // theta = [-dtheta, pi + dtheta]
        //r[2] = z/abs(Params::box_zmax + Params::sph_dz)*(pi + Params::sph_dphi);              // phi   = [-pi - dphi, pi + dphi]
        r[1] = y/abs(Params::box_ymax + Params::sph_dy)*(Params::sph_theta_max - pi/2 + Params::sph_dtheta) + pi/2; // theta = [-dtheta, pi + dtheta]
        r[2] = z/abs(Params::box_zmax + Params::sph_dz)*(Params::sph_phi_max + Params::sph_dphi);                   // phi   = [-pi - dphi, pi + dphi]
    }
}

//! (SPHERICAL) Transformation: Position (x, y, z) -> (r, theta, phi)
inline void sph_transf_S2H_R(gridreal r[3])
{
    gridreal R     = r[0];
    gridreal theta = r[1];
    gridreal phi   = r[2];
    if(Params::sph_BC_type == 0) {
        //r[1] = (2.0*theta/pi - 1.0)/(1.0 - 2.0*Params::sph_alpha_theta)*abs(Params::box_ymax + Params::sph_dy); // y = [-y0 + dy, y0 - dy]
        //r[2] = (phi/pi)*abs(Params::box_zmax + Params::sph_dz)/(1.0 - Params::sph_alpha_phi);                   // z = [-z0 + sph_dz, z0 - sph_dz]
        r[1] = (theta - pi/2)/(Params::sph_theta_max - pi/2)*abs(Params::box_ymax + Params::sph_dy);               // y = [-y0 + dy, y0 - dy]
        r[2] =  phi*abs(Params::box_zmax + Params::sph_dz)/Params::sph_phi_max;                                   // z = [-z0 + sph_dz, z0 - sph_dz]
    } else if(Params::sph_BC_type == 1) {
        //r[1] = (2.0*theta/pi-1.0)*abs(Params::box_ymax); // y = [-y0, y0]
        //r[2] = (phi/pi)*abs(Params::box_zmax);           // z = [-z0, z0]
        //r[1] = (theta - pi/2.0)/(pi/2.0 + Params::sph_dtheta)*abs(Params::box_ymax + Params::sph_dy); // y = [-y0 - dy, y0 + dy]
        //r[2] = phi/(pi + Params::sph_dphi)*abs(Params::box_zmax + Params::sph_dz);                    // z = [-z0 - sph_dz, z0 + sph_dz]
        r[1] = (theta - pi/2)/(Params::sph_theta_max - pi/2 + Params::sph_dtheta)*abs(Params::box_ymax + Params::sph_dy); // y = [-y0 - dy, y0 + dy]
        r[2] = phi/(Params::sph_phi_max + Params::sph_dphi)*abs(Params::box_zmax + Params::sph_dz);                       // z = [-z0 - sph_dz, z0 + sph_dz]
    }
}

// =================================================================================
// ================================== VELOCITIES ===================================
// =================================================================================


//! (SPHERICAL) Transformation: Velocities (Vr, Vtheta, Vphi) -> (Vx, Vy, Vz). Note! Here R[3] must be in spherical coordinates.
inline void sph_transf_S2C_V(gridreal r[3], gridreal v[3])
{
    gridreal R     = r[0];
    gridreal theta = r[1];
    gridreal phi   = r[2];
    gridreal vr_dot     = v[0];
    gridreal vtheta_dot = v[1];
    gridreal vphi_dot   = v[2];
    v[0] = vr_dot*sin(theta)*cos(phi) + vtheta_dot*R*cos(theta)*cos(phi) - vphi_dot*R*sin(theta)*sin(phi); // vx
    v[1] = vr_dot*sin(theta)*sin(phi) + vtheta_dot*R*cos(theta)*sin(phi) + vphi_dot*R*sin(theta)*cos(phi); // vy
    v[2] = vr_dot*cos(theta)          - vtheta_dot*R*sin(theta);                                           // vz
}

//! (SPHERICAL) Transformation: Velocities (Vx, Vy, Vz) -> (Vr, Vtheta, Vphi)
inline void sph_transf_C2S_V(gridreal r[3], gridreal v[3])
{
    gridreal R     = r[0];
    gridreal theta = r[1];
    gridreal phi   = r[2];
    gridreal vx = v[0];
    gridreal vy = v[1];
    gridreal vz = v[2];
    v[0] =  vx*sin(theta)*cos(phi) + vy*sin(theta)*sin(phi) + vz*cos(theta);  // vr
    v[1] =  vx*cos(theta)*cos(phi) + vy*cos(theta)*sin(phi) - vz*sin(theta);  // vtheta
    v[2] = -vx*sin(phi)            + vy*cos(phi);                             // vphi
    // Real spherical velocities
    v[0] =  v[0];               // vr_dot
    v[1] =  v[1]/R;             // vtheta_dot
    v[2] =  v[2]/R/sin(theta);  // vphi_dot
}

//! (SPHERICAL) Transformation: Velocities (Vr, Vtheta, Vphi) -> (Vx, Vy, Vz). This is just derivation from sph_transf_H2S_R.
inline void sph_transf_H2S_V(gridreal v[3])
{
    gridreal vh_theta = v[1];
    gridreal vh_phi   = v[2];
    if (Params::sph_BC_type == 0) {
        //v[1] = vh_theta/abs(Params::box_ymax + Params::sph_dy)*pi/2.0*(1.0 - 2.0*Params::sph_alpha_theta);
        //v[2] = vh_phi/abs(Params::box_zmax + Params::sph_dz)*pi*(1.0 - Params::sph_alpha_phi);
        v[1] = vh_theta/abs(Params::box_ymax + Params::sph_dy)*(Params::sph_theta_max - pi/2);
        v[2] = vh_phi/abs(Params::box_zmax + Params::sph_dz)*Params::sph_phi_max;
    }
    if (Params::sph_BC_type == 1) {
        //v[1] = vh_theta*pi/abs(Params::box_ymax)/2.0;
        //v[2] = vh_phi*pi/abs(Params::box_zmax);
        //v[1] = vh_theta/abs(Params::box_ymax + Params::sph_dy)*(pi/2.0 + Params::sph_dtheta);
        //v[2] = vh_phi/abs(Params::box_zmax + Params::sph_dz)*(pi + Params::sph_dphi);
        v[1] = vh_theta/abs(Params::box_ymax + Params::sph_dy)*(Params::sph_theta_max - pi/2 + Params::sph_dtheta);
        v[2] = vh_phi/abs(Params::box_zmax + Params::sph_dz)*(Params::sph_phi_max + Params::sph_dphi);
    }
}

//! (SPHERICAL) Transformation: Velocities (V_r, V_theta, V_phi) -> (Vh_r, Vh_theta, Vh_phi). This is just derivation from sph_transf_H2S_R.
inline void sph_transf_S2H_V(gridreal v[3])
{
    gridreal v_theta = v[1];
    gridreal v_phi   = v[2];
    if (Params::sph_BC_type == 0) {
        //v[1] = (2.0*v_theta/pi)/(1.0 - 2.0*Params::sph_alpha_theta)*abs(Params::box_ymax + Params::sph_dy);
        //v[2] = (v_phi/pi)*abs(Params::box_zmax + Params::sph_dz)/(1.0 - Params::sph_alpha_phi);
        v[1] =  v_theta/(Params::sph_theta_max - pi/2)*abs(Params::box_ymax + Params::sph_dy);
        v[2] =  v_phi*abs(Params::box_zmax + Params::sph_dz)/Params::sph_phi_max;
    } else if(Params::sph_BC_type == 1) {
        //v[1] = 2.0*v_theta*abs(Params::box_ymax)/pi;
        //v[2] = v_phi*abs(Params::box_zmax)/pi;
        //v[1] = v_theta/(pi/2.0 + Params::sph_dtheta)*abs(Params::box_ymax + Params::sph_dy);
        //v[2] = v_phi/(pi + Params::sph_dphi)*abs(Params::box_zmax + Params::sph_dz);
        v[1] = v_theta/(Params::sph_theta_max - pi/2 + Params::sph_dtheta)*abs(Params::box_ymax + Params::sph_dy);
        v[2] = v_phi/(Params::sph_phi_max + Params::sph_dphi)*abs(Params::box_zmax + Params::sph_dz);
    }
}

// =================================================================================
// ==================================== VECTORS ====================================
// =================================================================================

//! (SPHERICAL) Transformation: Vectors (Ar, Atheta, Aphi) -> (Ax, Ay, Az). Note! Here R[3] must be in spherical coordinates.
inline void sph_transf_S2C_A(gridreal r[3], real A[3])
{
    gridreal R     = r[0];
    gridreal theta = r[1];
    gridreal phi   = r[2];
    real A_r     = A[0];
    real A_theta = A[1];
    real A_phi   = A[2];
    A[0] = A_r*sin(theta)*cos(phi) + A_theta*cos(theta)*cos(phi) - A_phi*sin(phi); // Ax
    A[1] = A_r*sin(theta)*sin(phi) + A_theta*cos(theta)*sin(phi) + A_phi*cos(phi); // Ay
    A[2] = A_r*cos(theta)          - A_theta*sin(theta);                           // Az
}

//! (SPHERICAL) Transformation: Vectors (Ax, Ay, Az) -> (Ar, Atheta, Aphi)
inline void sph_transf_C2S_A(gridreal r[3], real A[3])
{
    gridreal R     = r[0];
    gridreal theta = r[1];
    gridreal phi   = r[2];
    real Ax = A[0];
    real Ay = A[1];
    real Az = A[2];
    A[0] =  Ax*sin(theta)*cos(phi) + Ay*sin(theta)*sin(phi) + Az*cos(theta);  // Ar
    A[1] =  Ax*cos(theta)*cos(phi) + Ay*cos(theta)*sin(phi) - Az*sin(theta);  // Atheta
    A[2] = -Ax*sin(phi)            + Ay*cos(phi);                             // Aphi
}

//! (SPHERICAL) Note! Here R[3] must be in spherical coordinates
inline void sph_transf_S2C_A1(gridreal r[3], gridreal A[3])
{
    gridreal R     = r[0];
    gridreal theta = r[1];
    gridreal phi   = r[2];
    gridreal A_r     = A[0];
    gridreal A_theta = A[1];
    gridreal A_phi   = A[2];
    A[0] = A_r*sin(theta)*cos(phi) + A_theta*cos(theta)*cos(phi) - A_phi*sin(phi); // Ax
    A[1] = A_r*sin(theta)*sin(phi) + A_theta*cos(theta)*sin(phi) + A_phi*cos(phi); // Ay
    A[2] = A_r*cos(theta)          - A_theta*sin(theta);                           // Az
}

//! (SPHERICAL) Transformation: Vectors (Ax, Ay, Az) -> (Ar, Atheta, Aphi)
inline void sph_transf_C2S_A1(gridreal r[3], gridreal A[3])
{
    gridreal R     = r[0];
    gridreal theta = r[1];
    gridreal phi   = r[2];
    gridreal Ax = A[0];
    gridreal Ay = A[1];
    gridreal Az = A[2];
    A[0] =  Ax*sin(theta)*cos(phi) + Ay*sin(theta)*sin(phi) + Az*cos(theta);  // Ar
    A[1] =  Ax*cos(theta)*cos(phi) + Ay*cos(theta)*sin(phi) - Az*sin(theta);  // Atheta
    A[2] = -Ax*sin(phi)            + Ay*cos(phi);                             // Aphi
}

//! (SPHERICAL) Transformation: VTK coordinates Hybrid (x, y, z) -> (r, theta, phi)
inline void sph_transf_H2S_VTK(gridreal coords[3])
{
    sph_transf_H2S_R(coords);
    gridreal r     = coords[0];
    gridreal theta = coords[1];
    gridreal phi   = coords[2];
    coords[0] = r*sin(theta)*cos(phi);
    coords[1] = r*sin(theta)*sin(phi);
    coords[2] = r*cos(theta);
}

#endif

#endif

