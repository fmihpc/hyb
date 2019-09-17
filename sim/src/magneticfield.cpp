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

#include "magneticfield.h"
#include "simulation.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "transformations.h"
#endif

using namespace std;

extern Tgrid g;
extern Params simuConfig;

//! Dummy constructor
VectorField::VectorField() { }

//! Dummy virtual destructor
VectorField::~VectorField() { }

//! Dummy implementation for a virtual interface function
void VectorField::getValue(const gridreal r[3],datareal V[3])
{
    WARNINGMSG("dummy implementation function called");
}

//! Check if the object is initialized
bool VectorField::isDefined()
{
    if(name.length() <= 0) {
        return false;
    } else {
        return true;
    }
}

//! Default constructor
MagneticFieldProfile::MagneticFieldProfile()
{
    // Make sure the class is not used without
    // proper initialization (this aborts execution
    // if getValue called).
    this->ptr = &MagneticFieldProfile::defaultFunction;
    resetParameters();
}

#define ELSEIF_MAGNETIC(func) else if(funcName.compare(#func) == 0) { this->ptr = &MagneticFieldProfile::func; setArgs_ ## func(); }

//! Constructor
MagneticFieldProfile::MagneticFieldProfile(string funcName,vector<real> args)
{
    this->name = funcName;
    this->ptr = &MagneticFieldProfile::defaultFunction;
    this->args = args;
    resetParameters();
    // MAGNETIC FIELD PROFILES
    if(funcName.compare("") == 0)  {
        ERRORMSG("empty magnetic field function name");
        doabort();
    }
    ELSEIF_MAGNETIC(constantB)
    ELSEIF_MAGNETIC(constantBx)
    ELSEIF_MAGNETIC(laminarFlowAroundSphereB)
    ELSEIF_MAGNETIC(laminarFlowAroundSphereBx)
    ELSEIF_MAGNETIC(dipoleB)
    ELSEIF_MAGNETIC(translateDipoleB)
    ELSEIF_MAGNETIC(generalDipoleB)
    ELSEIF_MAGNETIC(hemisphericDipoleB)
    ELSEIF_MAGNETIC(dipoleCuspB)
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    ELSEIF_MAGNETIC(sph_laminarFlowAroundSphereBx)
    ELSEIF_MAGNETIC(sph_laminarFlowAroundSphereBz)
    ELSEIF_MAGNETIC(sph_dipoleB)
    ELSEIF_MAGNETIC(sph_Br)
#endif
    else {
        ERRORMSG2("bad magnetic field function name",funcName);
        doabort();
    }
}

//! Destructor
MagneticFieldProfile::~MagneticFieldProfile() { }

//! Returns the magnetic field value at point r
void MagneticFieldProfile::getValue(const gridreal r[3],datareal B[3])
{
    (this->*ptr)(r,B);
}

//! Default function, which aborts the program if called
void MagneticFieldProfile::defaultFunction(const gridreal r[3],datareal B[3])
{
    ERRORMSG("function pointer not set");
    doabort();
}

//! Reset private class variables
void MagneticFieldProfile::resetParameters()
{
    // MAGNETIC FIELD PARAMETERS
    Bx = By = Bz = Btot = R = R3 = 0.0;
    sinTheta = cosTheta = sinPhi = cosPhi = xOrigin = yOrigin = zOrigin = 0.0;
    dipSurfB = dipSurfR = dipMomCoeff = dipRmin2 = hemiCoeffDip = hemiCoeffQuad = 0.0;
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    sph_B0 = sph_R0 = sph_d = 0.0;
#endif
}

// MAGNETIC FIELD PROFILES

/** \brief Set function arguments for homogeneous B field
 *
 * Homogenous (Bx,By,Bz) field.
 * Config file format: constantB Bx By Bz
 * Bx,By,Bz = components of the field
 */
void MagneticFieldProfile::setArgs_constantB()
{
    if(args.size() != 3) {
        ERRORMSG2("function takes three arguments",name);
        doabort();
    }
    Bx = args[0];
    By = args[1];
    Bz = args[2];
}

//! Constant homogeneous B field
void MagneticFieldProfile::constantB(const gridreal r[3],datareal B[3])
{
    B[0] += Bx;
    B[1] += By;
    B[2] += Bz;
}

/** \brief Set function arguments for homogeneous Bx field
 *
 * Homogenous (Bx,0,0) field.
 * Config file format: constantBx Bx
 * Bx = component of the field
 */
void MagneticFieldProfile::setArgs_constantBx()
{
    if(args.size() != 1) {
        ERRORMSG2("function takes one argument",name);
        doabort();
    }
    Bx = args[0];
}
//! Constant homogeneous Bx field
void MagneticFieldProfile::constantBx(const gridreal r[3],datareal B[3])
{
    B[0] += Bx;
}

/** \brief Set function arguments for general potential flow field around a sphere
 *
 * B = (Bx,By,Bz) when r >> R.
 * Config file format: laminarFlowAroundSphereBx R Bx By Bz
 * R = radius of the sphere, inside R the field is zero
 * Bx,By,Bz = components of the field
 */
void MagneticFieldProfile::setArgs_laminarFlowAroundSphereB()
{
    if(args.size() != 4) {
        ERRORMSG2("function takes four arguments",name);
        doabort();
    }
    R = args[0];
    Bx = args[1];
    By = args[2];
    Bz = args[3];
    R3 = cube(R);
    Btot = sqrt( sqr(Bx) + sqr(By) + sqr(Bz) );
    real Beq = sqrt( sqr(Bx) + sqr(By) );
    real phi, theta;
    if (Btot < 1e-25) {
        phi = 0;
        theta = 0;
        ERRORMSG2("Initial magfield is less than 1e-25 T",name);
    } else {
        if (fabs(Bx) < 1e-8*Btot) {
            phi = sign(By)*pi/2;
        } else {
            phi = atan( By/Bx) + ( 1-sign(Bx) )*pi/2;
        }
        if (Beq < 1e-8*Btot) {
            theta = ( 1 - sign(Bz) )/2*pi;
        } else {
            theta = acos(Bz/Btot);
        }
    }
    real theta_lat = pi/2 - theta;
    sinTheta = sin(theta_lat);
    cosTheta = cos(theta_lat);
    sinPhi = sin(phi);
    cosPhi = cos(phi);
}

//! General potential flow field around a sphere
void MagneticFieldProfile::laminarFlowAroundSphereB(const gridreal r[3],datareal B[3])
{
    /*
     Components of the magnetic field:
     B_x/B = -3/2*R^3*x^2/r^5 + 1 + 1/2*R^3/r^3
     B_y/B = -3/2*R^3*x*y/r^5
     B_z/B = -3/2*R^3*x*z/r^5

     Rotation matrix (theta_lat around y): (first this one, in a/b a is regular rotation and b is antirotation)
     x'         cos theta_lat     0   -/+sin theta_lat   x
     y'   =           0           1            0         y
     z'      +/-sin theta_lat     0      cos theta_lat   z

     Rotation matrix (phi around z):
     x'       cos phi   sin phi    0      x
     y'   =  -sin phi   cos phi    0      y
     z'          0         0       1      z
     */
    real rr = normvec(r);
    // Laminar solution is valid outside the radius R, otherwise don't add the laminar field
    if (rr <= R) {
        return;
    }
    // Some powers of r
    real r2 = sqr(rr);
    real r3 = r2*rr;
    real r5 = r3*r2;
    // Rotation of coordinates
    real x, x1, y, y1, z, z1, B_x1, B_y1, B_z1;
    x =  r[0]*cosPhi + r[1]*sinPhi;
    y = -r[0]*sinPhi + r[1]*cosPhi;
    z =  r[2];
    x1 = x*cosTheta + z*sinTheta;  // these need to be calculated with x, y, z
    y1 = y;
    z1 = -x*sinTheta + z*cosTheta;
    // laminar flow in X direction
    real coeff = -1.5 * Btot * R3 * x1 / r5;
    datareal B_x = coeff * x1 + Btot * (1 + 0.5 * R3 / r3);
    datareal B_y = coeff * y1;
    datareal B_z = coeff * z1;
    // Rotation of vector field (opposite direction)
    B_x1 = B_x*cosTheta - B_z*sinTheta;
    B_y1 = B_y;
    B_z1 = B_x*sinTheta + B_z*cosTheta;
    B[0] += B_x1*cosPhi - B_y1*sinPhi;
    B[1] += B_x1*sinPhi + B_y1*cosPhi;
    B[2] += B_z1;
}

/** \brief Set function arguments for potential flow field around a sphere in x-direction
 *
 * B = (Bx,0,0) when r >> R.
 * Config file format: laminarFlowAroundSphereBx R Bx
 * R = radius of the sphere, inside R the field is zero
 * Bx = x-component of the field
 */
void MagneticFieldProfile::setArgs_laminarFlowAroundSphereBx()
{
    if(args.size() != 2) {
        ERRORMSG2("function takes two arguments",name);
        doabort();
    }
    R = args[0];
    Bx = args[1];
    R3 = cube(R);
}

//! Potential flow field around a sphere in x-direction
void MagneticFieldProfile::laminarFlowAroundSphereBx(const gridreal r[3],datareal B[3])
{
    real rr = normvec(r);
    // Laminar solution is valid outside the radius R, otherwise don't add the laminar field
    if (rr <= R) {
        return;
    }
    // Some powers of r
    real r3 = pow (rr, 3);
    real r5 = pow (rr, 5);
    real coeff = -1.5*Bx*R3*r[0]/r5;
    B[0] += coeff*r[0] + Bx*( 1 + 0.5 * (R3/r3) );
    B[1] += coeff*r[1];
    B[2] += coeff*r[2];
}


/** \brief Set function arguments for a dipole field in the z-direction located at (0,0,0)
 *
 * Config file format: dipoleB surfB surfR Rmin
 * surfB = magnitude of the dipole field at the surface (equator), sign(surfB) gives the orientation
 * surfR = surface radius
 * Rmin  = the field is set to zero inside this radius (dipole field diverges as r->0)
 */
void MagneticFieldProfile::setArgs_dipoleB()
{
    if(args.size() != 3) {
        ERRORMSG2("function takes three arguments",name);
        doabort();
    }
    dipSurfB = args[0];
    dipSurfR = args[1];
    dipMomCoeff = 3.0*dipSurfB*cube(dipSurfR);
    dipRmin2 = sqr(args[2]);
}

//! Dipole field in the z-direction located at (0,0,0)
void MagneticFieldProfile::dipoleB(const gridreal r[3],datareal B[3])
{
    const gridreal x = r[0];
    const gridreal y = r[1];
    const gridreal z = r[2];
    const real r2 = sqr(x) + sqr(y) + sqr(z);
    if(r2 < dipRmin2) {
        return;
    }
    const real rr = sqrt(r2);
    const real r5 = sqr(r2)*rr;
    const real coeff = dipMomCoeff/r5;
    B[0] += coeff*x*z;
    B[1] += coeff*y*z;
    B[2] += coeff*(sqr(z) - r2/3.0);
}

/** \brief Set function arguments for dipole field in the z-direction located at (x0,y0,z0)
 *
 * Config file format: translateDipoleB surfB surfR Rmin x0 y0 z0
 * surfB = magnitude of the dipole field at the surface (equator), sign(surfB) gives the orientation
 * surfR = surface radius
 * Rmin  = the field is set to zero inside this radius (dipole field diverges as r->0)
 * x0    = x coordinate of the dipole
 * y0    = y coordinate of the dipole
 * z0    = z coordinate of the dipole
 */
void MagneticFieldProfile::setArgs_translateDipoleB()
{
    if(args.size() != 6) {
        ERRORMSG2("function takes six arguments",name);
        doabort();
    }
    dipSurfB = args[0];
    dipSurfR = args[1];
    dipMomCoeff = 3.0*dipSurfB*cube(dipSurfR);
    dipRmin2 = sqr(args[2]);
    xOrigin = args[3];
    yOrigin = args[4];
    zOrigin = args[5];
}

//! Dipole field in the z-direction located at (x0,y0,z0)
void MagneticFieldProfile::translateDipoleB(const gridreal r[3],datareal B[3])
{
    const gridreal x = r[0] - xOrigin;
    const gridreal y = r[1] - yOrigin;
    const gridreal z = r[2] - zOrigin;
    const real r2 = sqr(x) + sqr(y) + sqr(z);
    if(r2 < dipRmin2) {
        return;
    }
    const real rr = sqrt(r2);
    const real r5 = sqr(r2)*rr;
    const real coeff = dipMomCoeff/r5;
    B[0] += coeff*x*z;
    B[1] += coeff*y*z;
    B[2] += coeff*(sqr(z) - r2/3.0);
}

/** \brief Set function arguments for a rotated dipole field located at (x0,y0,z0)
 *
 * Config file format: generalDipoleB surfB surfR Rmin x0 y0 z0 theta phi
 * surfB = magnitude of the dipole field at the surface (equator), sign(surfB) gives the orientation
 * surfR = surface radius
 * Rmin  = the field is set to zero inside this radius (dipole field diverges as r->0)
 * x0    = x coordinate of the dipole
 * y0    = y coordinate of the dipole
 * z0    = z coordinate of the dipole
 * theta = rotation angle about the y-axis
 * phi   = rotation angle about the x-axis
 */
void MagneticFieldProfile::setArgs_generalDipoleB()
{
    if(args.size() != 8) {
        ERRORMSG2("function takes eight arguments",name);
        doabort();
    }
    dipSurfB = args[0];
    dipSurfR = args[1];
    dipMomCoeff = 3.0*dipSurfB*cube(dipSurfR);
    dipRmin2 = sqr(args[2]);
    xOrigin = args[3];
    yOrigin = args[4];
    zOrigin = args[5];
    real theta = args[6];
    if(theta >= -180 && theta <= 180) {
        theta *= pi/180.0;
    } else {
        ERRORMSG2("theta (rotation about y-axis) should be between -180 and 180 degrees",name);
        doabort();
    }
    if(fabs(theta) < 1e-10) {
        sinTheta = 0.0;
        cosTheta = 1.0;
    } else {
        sinTheta = sin(theta);
        cosTheta = cos(theta);
    }
    real phi = args[7];
    if(phi >= -180 && phi <= 180) {
        phi *= pi/180.0;
    } else {
        ERRORMSG2("phi (rotation about the x-axis) should be between -180 and 180 degrees",name);
        doabort();
    }
    if(fabs(phi) < 1e-10) {
        sinPhi = 0.0;
        cosPhi = 1.0;
    } else {
        sinPhi = sin(phi);
        cosPhi = cos(phi);
    }
}

//! Rotated dipole field located at (x0,y0,z0)
void MagneticFieldProfile::generalDipoleB(const gridreal r[3],datareal B[3])
{
    const gridreal x = r[0] - xOrigin;
    const gridreal y = r[1] - yOrigin;
    const gridreal z = r[2] - zOrigin;
    const real r2 = sqr(x) + sqr(y) + sqr(z);
    if(r2 < dipRmin2) {
        return;
    }
    const real rr = sqrt(r2);
    const real r5 = sqr(r2)*rr;
    const real coeff = dipMomCoeff/r5;
    B[0] += coeff*x*z;
    B[1] += coeff*y*z;
    B[2] += coeff*(sqr(z) - r2/3.0);
    //Rotation matrix R_y(theta) (=theta around y):
    // x'      cos theta     0   sin theta   x
    // y'   =      0         1      0        y
    // z'      -sin theta    0   cos theta   z
    //Rotation matrix R_x(phi) (=phi around x):
    // x'          1         0      0       x
    // y'   =      0   cos phi   -sin phi   y
    // z'          0   sin phi    cos phi   z
    // Rotation of coordinates
    // (xx,yy,zz) = R_y(theta) (x,y,z)
    const real xx = x*cosTheta + z*sinTheta;
    const real yy = y;
    const real zz = -x*sinTheta + z*cosTheta;
    // (xxx,yyy,zzz) = R_x(phi) (xx,yy,zz)
    const real xxx = xx;
    const real yyy = yy*cosPhi - zz*sinPhi;
    const real zzz = yy*sinPhi + zz*cosPhi;
    // Dipole in z-direction
    const datareal Bx = coeff * xxx * zzz;
    const datareal By = coeff * yyy * zzz;
    const datareal Bz = coeff * (sqr(zzz) - r2/3.0);
    // Rotation of field
    // B_ = R_x(-phi) B
    const real Bxx = Bx;
    const real Byy = By*cosPhi + Bz*sinPhi;
    const real Bzz = -By*sinPhi + Bz*cosPhi;
    // B = R_y(-theta) B_
    B[0] += Bxx*cosTheta - Bzz*sinTheta;
    B[1] += Byy;
    B[2] += Bxx*sinTheta + Bzz*cosTheta;
}

/** \brief Set function arguments for a hemispheric multipole field located at (x0,y0,z0)
 *
 * Config file format: hemisphericDipoleB hemiCoeffDip hemiCoeffQuad dipSurfB dipSurfR dipRmin x0 y0 z0
 * hemiCoeffDip  = dipole coefficient
 * hemiCoeffQuad = quadrupole coefficient
 * dipSurfB = magnitude of the dipole field at the surface (equator), sign(dipSurfB) gives the orientation
 * dipSurfR = surface radius
 * dipRmin = the field is set to zero inside this radius
 * x0    = x coordinate of the multipole
 * y0    = y coordinate of the multidipole
 * z0    = z coordinate of the multidipole
 */
void MagneticFieldProfile::setArgs_hemisphericDipoleB()
{
    if(args.size() != 8) {
        ERRORMSG2("function takes eight arguments",name);
        doabort();
    }
    hemiCoeffDip  = args[0];
    hemiCoeffQuad = args[1];
    dipSurfB = args[2];
    dipSurfR = args[3];
    dipMomCoeff = 3.0*dipSurfB*cube(dipSurfR);
    dipRmin2 = sqr(args[4]);
    xOrigin = args[5];
    yOrigin = args[6];
    zOrigin = args[7];
}

//! Hemispheric multipole field located at (x0,y0,z0)
void MagneticFieldProfile::hemisphericDipoleB(const gridreal r[3],datareal B[3])
{
    const gridreal x = r[0] - xOrigin;
    const gridreal y = r[1] - yOrigin;
    const gridreal z = r[2] - zOrigin;
    const real r2 = sqr(x) + sqr(y) + sqr(z);
    if(r2 < dipRmin2) {
        return;
    }
    const real rr = sqrt(r2);
    const real r5 = sqr(r2)*rr;
    const real coeff = dipMomCoeff/r5;
    B[0] += coeff*x*(hemiCoeffDip*z + 0.5*hemiCoeffQuad*dipSurfR*(5*sqr(z)/r2 - 1));
    B[1] += coeff*y*(hemiCoeffDip*z + 0.5*hemiCoeffQuad*dipSurfR*(5*sqr(z)/r2 - 1));
    B[2] += coeff*(hemiCoeffDip*(sqr(z) - r2/3.0) + 0.5*hemiCoeffQuad*dipSurfR*z*(5*sqr(z) - 3*r2)/r2);
}

/** \brief Set function arguments for dipole field in the x-direction located at (x0,y0,z0)
 *
 * Config file format: dipoleCuspB surfB surfR Rmin x0 y0 z0
 * surfB = magnitude of the dipole field at the surface (equator), sign(surfB) gives the orientation
 * surfR = surface radius
 * Rmin  = the field is set to zero inside this radius (dipole field diverges as r->0)
 * x0    = x coordinate of the dipole
 * y0    = y coordinate of the dipole
 * z0    = z coordinate of the dipole
 */
void MagneticFieldProfile::setArgs_dipoleCuspB()
{
    if(args.size() != 6) {
        ERRORMSG2("function takes six arguments",name);
        doabort();
    }
    dipSurfB = args[0];
    dipSurfR = args[1];
    dipMomCoeff = 3.0*dipSurfB*cube(dipSurfR);
    dipRmin2 = sqr(args[2]);
    xOrigin = args[3];
    yOrigin = args[4];
    zOrigin = args[5];
}

//! Dipole field in the x-direction located at (x0,y0,z0)
void MagneticFieldProfile::dipoleCuspB(const gridreal r[3],datareal B[3])
{
    const gridreal x = r[0] - xOrigin;
    const gridreal y = r[1] - yOrigin;
    const gridreal z = r[2] - zOrigin;
    const real r2 = sqr(x) + sqr(y) + sqr(z);
    if(r2 < dipRmin2) {
        return;
    }
    const real rr = sqrt(r2);
    const real r5 = sqr(r2)*rr;
    const real coeff = dipMomCoeff/r5;
    //B[0] += coeff*x*z;
    //B[1] += coeff*y*z;
    //B[2] += coeff*(sqr(z) - r2/3.0);
    B[0] += coeff*(sqr(x) - r2/3.0);
    B[1] += coeff*y*x;
    B[2] += coeff*z*x;
}

//! Add constant magnetic field in B
void addConstantMagneticField(const gridreal r[3], datareal B[3])
{
    // Go thru the constant magnetic field function calls (pointers actually)
    for(unsigned int i=0; i < Params::constantMagneticFieldProfile.size(); ++i) {
        Params::constantMagneticFieldProfile[i].getValue(r,B);
    }
}

//! Set initial magnetic field in B
void setInitialMagneticField(const gridreal r[3], datareal B[3])
{
    B[0] = 0;
    B[1] = 0;
    B[2] = 0;
    // Go thru the initial magnetic field function calls (pointers actually)
    for(unsigned int i=0; i < Params::initialMagneticFieldProfile.size(); ++i) {
        Params::initialMagneticFieldProfile[i].getValue(r,B);
    }
}

//! Initialize magnetic field profiles
void initializeMagneticField()
{
    MSGFUNCTIONCALL("initializeMagneticField");
    // Set initial magnetic field functions
    mainlog << "|--------------- INITIAL MAGNETIC FIELD PROFILES ---------------|\n";
    vector<string> tempNames;
    vector< vector<real> > tempArgs;
    bool nonZeroFuncs = simuConfig.getFunctionNamesAndArgs("initialMagneticFieldFUNC",tempNames,tempArgs);
    if(nonZeroFuncs == false || tempNames.size() <= 0) {
        mainlog << "none\n";
    } else if ( tempNames.size() != tempArgs.size() ) {
        ERRORMSG("internal error");
        doabort();
    } else {
        for(unsigned int i=0; i < tempNames.size(); ++i) {
            mainlog << tempNames[i] << " ";
            for(unsigned int j=0; j < tempArgs[i].size(); ++j) {
                mainlog << tempArgs[i][j] << " ";
            }
            mainlog << "\n";
        }
        Params::initialMagneticFieldProfile.clear();
        for(unsigned int i=0; i < tempNames.size(); ++i) {
            Params::initialMagneticFieldProfile.push_back( MagneticFieldProfile(tempNames[i],tempArgs[i]) );
            //g.set_B(Params::initialMagneticFieldProfile);
        }
        g.set_B(&setInitialMagneticField);
    }
    mainlog << "|---------------------------------------------------------------|\n\n";
    // Set constant magnetic field functions
    mainlog << "|--------------- CONSTANT MAGNETIC FIELD PROFILES ---------------|\n";
    tempNames.clear();
    tempArgs.clear();
    nonZeroFuncs = simuConfig.getFunctionNamesAndArgs("constantMagneticFieldFUNC",tempNames,tempArgs);
    if(nonZeroFuncs == false || tempNames.size() <= 0) {
        mainlog << "none\n";
    } else if ( tempNames.size() != tempArgs.size() ) {
        ERRORMSG("internal error");
        doabort();
    } else {
        for(unsigned int i=0; i < tempNames.size(); ++i) {
            mainlog << tempNames[i] << " ";
            for(unsigned int j=0; j < tempArgs[i].size(); ++j) {
                mainlog << tempArgs[i][j] << " ";
            }
            mainlog << "\n";
        }
        Params::constantMagneticFieldProfile.clear();
        for(unsigned int i=0; i < tempNames.size(); ++i) {
            Params::constantMagneticFieldProfile.push_back( MagneticFieldProfile(tempNames[i],tempArgs[i]) );
        }
    }
    mainlog << "|-----------------------------------------------------------------|\n";
    MSGFUNCTIONEND("initializeMagneticField");
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Set function arguments
void MagneticFieldProfile::setArgs_sph_laminarFlowAroundSphereBx()
{
    if(args.size() != 2) {
        ERRORMSG2("function takes two arguments",name);
        doabort();
    }
    R = args[0];
    Bx = args[1];
    R3 = cube(R);
}

//! (SPHERICAL) Spherical version of "laminarFlowAroundSphereBx"
void MagneticFieldProfile::sph_laminarFlowAroundSphereBx(const gridreal r[3],datareal B[3])
{
// Laminar solution is valid outside the radius R, otherwise don't add the laminar field
    if (r[0] <= R) {
        return;
    }
// We need to rename vector r in r_new
    gridreal r_new[3] = {r[0], r[1], r[2]};
    sph_transf_H2S_R(r_new);
    sph_transf_S2C_A(r_new, B);
    sph_transf_S2C_R(r_new);
    real rr = normvec(r_new);
//Propagation along x-axis
// Some powers of r
    real r3 = pow (rr, 3);
    real r5 = pow (rr, 5);
    real coeff = -1.5*Bx*R3*r_new[0]/r5;
    B[0] += coeff*r_new[0] + Bx*( 1 + 0.5 * (R3/r3) );
    B[1] += coeff*r_new[1];
    B[2] += coeff*r_new[2];
    sph_transf_C2S_r(r_new);
    sph_transf_C2S_A(r_new, B);
}

//! (SPHERICAL) Set function arguments
void MagneticFieldProfile::setArgs_sph_laminarFlowAroundSphereBz()
{
    if(args.size() != 2) {
        ERRORMSG2("function takes two arguments",name);
        doabort();
    }
    R = args[0];
    Bz = args[1];
    R3 = cube(R);
}

//! (SPHERICAL) Spherical version of "laminarFlowAroundSphereBx"
void MagneticFieldProfile::sph_laminarFlowAroundSphereBz(const gridreal r[3],datareal B[3])
{
// Laminar solution is valid outside the radius R, otherwise don't add the laminar field
    if (r[0] <= R) {
        return;
    }
// We need to rename vector r in r_new
    gridreal r_new[3] = {r[0], r[1], r[2]};
    sph_transf_H2S_R(r_new);
    sph_transf_S2C_A(r_new, B);
    sph_transf_S2C_R(r_new);
    real rr = normvec(r_new);
//Propagation along y-axis
// Some powers of r
    real r3 = pow (rr, 3);
    real r5 = pow (rr, 5);
    real coeff = -1.5*Bz*R3*r_new[2]/r5;
    B[0] += coeff*r_new[0];
    B[1] += coeff*r_new[1];
    B[2] += coeff*r_new[2] + Bz*( 1 + 0.5 * (R3/r3) );
    sph_transf_C2S_r(r_new);
    sph_transf_C2S_A(r_new, B);
}

// Spherical dipole field in any directions (x,y,z). The origin is located at (x=0,y=0,z=0).
// Config file format: generalDipoleB surfB rmin sph_dipole_dir
// surfB = magnitude of the dipole field at the surface (equator), sign(surfB) gives the orientation
// rmin = the field is set to zero inside this radius (dipole field diverges as r->0)
// sph_dipole_dir is the direction of magnetic moment 0-x, 1-y, 2-z
//! (SPHERICAL) Set function arguments
void MagneticFieldProfile::setArgs_sph_dipoleB()
{
    if(args.size() != 3) {
        ERRORMSG2("function takes three arguments",name);
        doabort();
    }
    dipSurfB = args[0];
    dipRmin2 = sqr(args[1]);
    sph_dipole_dir = args[2];
    // Dipole moment [A m^2]
    const real dipMom = dipSurfB*cube(Params::R_P)*4*pi/Params::mu_0;
    dipMomCoeff = 3.0*Params::mu_0*dipMom/(4.0*pi);
}

//! (SPHERICAL)
void MagneticFieldProfile::sph_dipoleB(const gridreal r[3],datareal B[3])
{
    gridreal r1[3] = {r[0], r[1], r[2]};
    datareal B_dipole[3] = {0.0 , 0.0, 0.0};
    sph_transf_H2S_R(r1);
    sph_transf_S2C_R(r1);
    real r2 = sqr(r1[0]) + sqr(r1[1]) + sqr(r1[2]);
    if(r2 < dipRmin2) {
        return;
    }
    real rr = sqrt(r2);
    real r5 = sqr(r2)*rr;
    real coeff = dipMomCoeff/r5;
    //Spherical version of Magnetic moment is along z axis
    //B[0] += coeff*2*cos(r1[1]);
    //B[1] += coeff*sin(r1[1]);
    //B[2] += 0.0;
    // Magnetic moment is along x axis
    if (sph_dipole_dir == 0) {
        B_dipole[0] += coeff*( sqr(r1[0]) - r2/3.0 );
        B_dipole[1] += coeff*r1[1]*r1[0];
        B_dipole[2] += coeff*r1[2]*r1[0];
    }
    // Magnetic moment is along y axis
    if (sph_dipole_dir == 1) {
        B_dipole[0] += coeff*r1[0]*r1[1];
        B_dipole[1] += coeff*( sqr(r1[1]) - r2/3.0 );
        B_dipole[2] += coeff*r1[2]*r1[1];
    }
    // Magnetic moment is along z axis
    if (sph_dipole_dir == 2) {
        B_dipole[0] += coeff*r1[0]*r1[2];
        B_dipole[1] += coeff*r1[1]*r1[2];
        B_dipole[2] += coeff*( sqr(r1[2]) - r2/3.0 );
    }
    sph_transf_C2S_r(r1);
    sph_transf_C2S_A(r1, B_dipole);
    B[0] += B_dipole[0];
    B[1] += B_dipole[1];
    B[2] += B_dipole[2];
}

//! (SPHERICAL) Set function arguments
void MagneticFieldProfile::setArgs_sph_Br()
{
    if(args.size() != 3) {
        ERRORMSG2("function takes three arguments",name);
        doabort();
    }
    sph_B0 = args[0]; // Magnetic field at the RO - radius
    sph_d  = args[1];
    sph_R0 = args[2]; // Internal radius of
}

//! (SPHERICAL) Br function: Br = B0*(R0/r)^d, B0 = magnetic field at the obstacle, R0 =  Radius of the obstacle, d - index of power
void MagneticFieldProfile::sph_Br(const gridreal r[3],datareal B[3])
{
    //Hybrid to spherical coordinate transformation
    gridreal r_sph[3] = {r[0],r[1],r[2]};
    sph_transf_H2S_R(r_sph);
    B[0] += sph_B0*pow(sph_R0/r_sph[0], sph_d)*pow(sin(r_sph[1]), 20*sph_d);
    B[1] += 0.0;
    B[2] += 0.0;
}

#endif

