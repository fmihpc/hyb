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
#include "atmosphere.h"
#include "simulation.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "transformations.h"
#endif

using namespace std;

extern Tgrid g;

//! Dummy constructor
ScalarField::ScalarField() { }

//! Dummy virtual destructor
ScalarField::~ScalarField() { }

//! Dummy implementation for a virtual interface function
real ScalarField::getValue(const gridreal[])
{
    WARNINGMSG("dummy implementation function called");
    return 0.0;
}

//! Check if the object is initialized
bool ScalarField::isDefined()
{
    if(name.length() <= 0) {
        return false;
    } else {
        return true;
    }
}

//! String summary
string ScalarField::toString()
{
    stringstream ss;
    ss << name << " ";
    for(unsigned int ii = 0; ii < args.size(); ++ii) {
        ss << args[ii] << " ";
    }
    return ss.str();
}

//! Default constructor
SpatialDistribution::SpatialDistribution()
{
    this->popid = 9999;
    // Make sure the class is not used without
    // proper initialization (this aborts execution
    // if getValue called).
    this->ptr = &SpatialDistribution::defaultFunction;
}

#define ELSEIF_DISTR_FUNC(func) else if(funcName.compare(#func) == 0) { this->ptr = &SpatialDistribution::func; }

//! Constructor
SpatialDistribution::SpatialDistribution(string funcName,vector<real> args,unsigned int popid)
{
    this->name = funcName;
    this->popid = popid;
    this->ptr = &SpatialDistribution::defaultFunction;
    this->args = args;
    if(this->popid > Params::pops.size()) {
        ERRORMSG2("given population id > max(id)",funcName);
        doabort();
    }
    // SPATIAL DISTRIBUTION PROFILES
    if(funcName.compare("") == 0) {
        ERRORMSG("empty spatial distribution function name");
        doabort();
    }
    ELSEIF_DISTR_FUNC(ionoConstantDayConstantNight)
    ELSEIF_DISTR_FUNC(ionoCosSzaDayConstantNight)
    ELSEIF_DISTR_FUNC(ionoCosSzaNoonToMidnight)
    ELSEIF_DISTR_FUNC(ionoLinearSzaDayConstantNight)
    ELSEIF_DISTR_FUNC(ionoLinearSzaNoonToMidnight)
    ELSEIF_DISTR_FUNC(ionoConstDayAndNight_anySolarDirection)
    ELSEIF_DISTR_FUNC(ionoCosSzaConstNight_anySolarDirection)
    ELSEIF_DISTR_FUNC(ionoCosSza_anySolarDirection)
    ELSEIF_DISTR_FUNC(ionoSmoothDistribution)
    ELSEIF_DISTR_FUNC(neutralDensityVenusHydrogen)
    ELSEIF_DISTR_FUNC(neutralDensityVenusOxygenHot)
    ELSEIF_DISTR_FUNC(photoionChamberlainTitan)
    ELSEIF_DISTR_FUNC(neutralDensityChamberlainT)
    ELSEIF_DISTR_FUNC(neutralDensityChamberlainH)
    ELSEIF_DISTR_FUNC(neutralDensityPowerLaw)
    ELSEIF_DISTR_FUNC(neutralDensityExponential)
    ELSEIF_DISTR_FUNC(ionizationConstant)
    ELSEIF_DISTR_FUNC(shadow)
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    ELSEIF_DISTR_FUNC(sph_ionoCosSzaDayConstantNight)
    ELSEIF_DISTR_FUNC(sph_neutralDensityVenusHydrogen)
    ELSEIF_DISTR_FUNC(sph_neutralDensityVenusOxygenHot)
    ELSEIF_DISTR_FUNC(sph_shadow)
#endif
    else {
        ERRORMSG2("bad spatial distribution function name",funcName);
        doabort();
    }
}

SpatialDistribution::~SpatialDistribution() { }

//! Returns the distribution value at point r
real SpatialDistribution::getValue(const gridreal r[])
{
    return (this->*ptr)(r);
}

//! Default function, which aborts the program if called
real SpatialDistribution::defaultFunction(const gridreal r[3])
{
    ERRORMSG("function pointer not set");
    doabort();
    return -1;
}

// SPATIAL DISTRIBUTION PROFILES

//! Ionospheric emission: Constant at all SZA
real SpatialDistribution::ionoConstantDayConstantNight(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 4) {
        ERRORMSG2("function takes four arguments",this->name);
        doabort();
    }
    const real daySideFactor = args[0];
    const real nightSideFactor = args[1];
    const real dummyScaleHeight = args[2];
    const real R = args[3];
    if (nightSideFactor < 0) {
        ERRORMSG2("nightSideFactor < 0",this->name);
        doabort();
    }
    if (daySideFactor < 0) {
        ERRORMSG2("daySideFactor < 0",this->name);
        doabort();
    }
    const fastreal rr = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
    const fastreal sza = acos(r[0]/rr);
    // No ions inside the radius
    if(rr < R) {
        return 0.0;
    }
    fastreal effdens;
    if(sza < pi/2) {
        effdens = daySideFactor;
    } else {
        effdens = nightSideFactor;
    }
    return effdens*exp(-(rr-R)/dummyScaleHeight);
}

//! Ionospheric emission: Goes from SZA=0 to 90 as cos(SZA) and SZA>90 is constant
real SpatialDistribution::ionoCosSzaDayConstantNight(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 4) {
        ERRORMSG2("function takes four arguments",this->name);
        doabort();
    }
    const real daySideFactor = args[0];
    const real nightSideFactor = args[1];
    const real dummyScaleHeight = args[2];
    const real R = args[3];
    if (nightSideFactor < 0) {
        ERRORMSG2("nightSideFactor < 0",this->name);
        doabort();
    }
    if (daySideFactor < 0) {
        ERRORMSG2("daySideFactor < 0",this->name);
        doabort();
    }
    const fastreal rr = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
    const fastreal sza = acos(r[0]/rr);
    // No ions inside the radius
    if(rr < R) {
        return 0.0;
    }
    fastreal effdens;
    if(sza < pi/2) {
        effdens = daySideFactor + (nightSideFactor - daySideFactor)*(1-cos(sza));
    } else {
        effdens = nightSideFactor;
    }
    return effdens*exp(-(rr-R)/dummyScaleHeight);
}

//! Ionospheric emission: cos(SZA) dependence at all SZA
real SpatialDistribution::ionoCosSzaNoonToMidnight(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 4) {
        ERRORMSG2("function takes four arguments",this->name);
        doabort();
    }
    const real daySideFactor = args[0];
    const real nightSideFactor = args[1];
    const real dummyScaleHeight = args[2];
    const real R = args[3];
    if (nightSideFactor < 0) {
        ERRORMSG2("nightSideFactor < 0",this->name);
        doabort();
    }
    if (daySideFactor < 0) {
        ERRORMSG2("daySideFactor < 0",this->name);
        doabort();
    }
    const fastreal rr = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
    const fastreal sza = acos(r[0]/rr);
    // No ions inside the radius
    if(rr < R) {
        return 0.0;
    }
    fastreal effdens = daySideFactor + (nightSideFactor - daySideFactor)*(1-cos(sza/2.0));
    return effdens*exp(-(rr-R)/dummyScaleHeight);
}

//! Ionospheric emission: Linear at SZA<90, constant at SZA>90
real SpatialDistribution::ionoLinearSzaDayConstantNight(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 4) {
        ERRORMSG2("function takes four arguments",this->name);
        doabort();
    }
    const real daySideFactor = args[0];
    const real nightSideFactor = args[1];
    const real dummyScaleHeight = args[2];
    const real R = args[3];
    if (nightSideFactor < 0) {
        ERRORMSG2("nightSideFactor < 0",this->name);
        doabort();
    }
    if (daySideFactor < 0) {
        ERRORMSG2("daySideFactor < 0",this->name);
        doabort();
    }
    const fastreal rr = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
    const fastreal sza = acos(r[0]/rr);
    // No ions inside the radius
    if(rr < R) {
        return 0.0;
    }
    fastreal effdens;
    if(sza < pi/2) {
        effdens = daySideFactor + (nightSideFactor - daySideFactor)*(sza/(pi/2));
    } else {
        effdens = nightSideFactor;
    }
    return effdens*exp(-(rr-R)/dummyScaleHeight);
}

//! Ionospheric emission: Linear at all SZA
real SpatialDistribution::ionoLinearSzaNoonToMidnight(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 4) {
        ERRORMSG2("function takes four arguments",this->name);
        doabort();
    }
    const real daySideFactor = args[0];
    const real nightSideFactor = args[1];
    const real dummyScaleHeight = args[2];
    const real R = args[3];
    if (nightSideFactor < 0) {
        ERRORMSG2("nightSideFactor < 0",this->name);
        doabort();
    }
    if (daySideFactor < 0) {
        ERRORMSG2("daySideFactor < 0",this->name);
        doabort();
    }
    const fastreal rr = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
    const fastreal sza = acos(r[0]/rr);
    // No ions inside the radius
    if(rr < R) {
        return 0.0;
    }
    fastreal effdens = daySideFactor + (nightSideFactor - daySideFactor)*(sza/pi);
    return effdens*exp(-(rr-R)/dummyScaleHeight);
}

//! Ionospheric emission with any solar direction: Goes from SZA=0 to 90 as cos(SZA) and SZA>90 is constant
real SpatialDistribution::ionoCosSzaConstNight_anySolarDirection(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes one argument",this->name);
        doabort();
    }
    const real nightSideFactor = args[0];
    if (ionoSmoothDistribution(r)<0.5) {
        return 0;
    }
    real cosSza=CosSZA(r);
    if (cosSza<0) {
        return nightSideFactor;
    } else {
        return nightSideFactor+cosSza*(1-nightSideFactor);
    }
}

//! Ionospheric emission with any solar direction: Constant at all SZA
real SpatialDistribution::ionoConstDayAndNight_anySolarDirection(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes one argument",this->name);
        doabort();
    }
    const real nightSideFactor = args[0];
    if (ionoSmoothDistribution(r)<0.5) {
        return 0;
    }
    if (CosSZA(r)<0) {
        return nightSideFactor;
    } else {
        return 1;
    }
}

//! Ionospheric emission with any solar direction: cos(SZA) dependence at all SZA
real SpatialDistribution::ionoCosSza_anySolarDirection(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes one argument",this->name);
        doabort();
    }
    const real nightSideFactor = args[0];
    if (ionoSmoothDistribution(r)<0.5) {
        return 0;
    }
    return nightSideFactor + (0.5*(1+CosSZA(r))*(1-nightSideFactor));
}

//! Ionospheric emission: Smoothing of any ionospheric distribution
real SpatialDistribution::ionoSmoothDistribution(const gridreal r[3])
{
    const real rr = sqr(r[0]) + sqr(r[1]) + sqr(r[2]);
    if (rr<sqr(Params::R_P)) return 0;
    static const real Rmin =
        min2(min3(fabs(Params::box_xmin),fabs(Params::box_ymin),fabs(Params::box_zmin)),
             min3(fabs(Params::box_xmax),fabs(Params::box_ymax),fabs(Params::box_zmax)));
    if (rr<sqr(Rmin)) return 1;
    else return 0;
}

//! Neutral corona profile: Venus thermal+hot hydrogen
real SpatialDistribution::neutralDensityVenusHydrogen(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes 1 argument",this->name);
        doabort();
    }
    const real zeroR = args[0];
    const fastreal rr = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
    // zero density below r
    if( rr < zeroR ) {
        return 0.0;
    }
    const fastreal sza = acos(r[0]/rr);
    const real exobaseR = Params::R_P + 170e3;
    // THERMAL HYDROGEN
    // zone: ZIP 1 = Dayside
    const fastreal n_0_zip1 = 1.32e11;
    const fastreal beta_zip1 = Params::G*Params::M_P*Params::m_p/(Params::k_B*285);
    // zone: ZIP 5 ~ Nightside
    const fastreal n_0_zip5 = 2.59e12;
    const fastreal beta_zip5 = Params::G*Params::M_P*Params::m_p/(Params::k_B*110);
    // densities
    const fastreal n_thermal_day = n_0_zip1*exp(-beta_zip1*(1/exobaseR - 1.0/rr));
    const fastreal n_thermal_night = n_0_zip5*exp(-beta_zip5*(1.0/exobaseR - 1.0/rr));
    // linear cos interpolation (noon -> midnight)
    const fastreal n_thermal = (1.0 - (sza/pi))*n_thermal_day + (sza/pi)*n_thermal_night;
    // HOT HYDROGEN
    // noon (sza = 0 deg)
    const fastreal a_1_noon = -6.2625e-5/1e3;
    const fastreal a_2_noon = 15.4817;
    const fastreal a_3_noon = 3.6414e4*1e3;
    // terminator (sza = 90 deg)
    const fastreal a_1_term = -8.4607e-5/1e3;
    const fastreal a_2_term = 15.9944;
    const fastreal a_3_term = 2.9743e4*1e3;
    // midnight (sza = 180 deg)
    const fastreal a_1_midn = -6.2309e-5/1e3;
    const fastreal a_2_midn = 15.2723;
    const fastreal a_3_midn = 4.3781e4*1e3;
    // densities
    const fastreal n_hot_noon = exp(a_1_noon*rr + a_2_noon + a_3_noon/rr);
    const fastreal n_hot_term = exp(a_1_term*rr + a_2_term + a_3_term/rr);
    const fastreal n_hot_midn = exp(a_1_midn*rr + a_2_midn + a_3_midn/rr);
    fastreal n_hot;
    // linear cos interpolation (noon -> terminator -> midnight)
    if(sza <= pi/2) { // dayside hemisphere
        n_hot = (1.0 - (sza/(pi/2)))*n_hot_noon + (sza/(pi/2))*n_hot_term;
    } else { // nightside hemisphere
        n_hot = (1.0 - (sza/(pi/2) - 1))*n_hot_term + (sza/(pi/2) - 1)*n_hot_midn;
    }
    return n_thermal + n_hot;
}

//! Neutral corona profile: Venus hot oxygen
real SpatialDistribution::neutralDensityVenusOxygenHot(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes 1 argument",this->name);
        doabort();
    }
    const real zeroR = args[0];
    const fastreal rr = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
    // zero density below r
    if( rr < zeroR ) {
        return 0.0;
    }
    const fastreal sza = acos(r[0]/rr);
    // ---- Hot neutral oxygen model from Gunell et al., Planetary and Space Science, 09/2005
    // dayside (originally from McElroy et al., 1982)
    const fastreal r_0_day = Params::R_V + 200e3;
    const fastreal beta_day = Params::G*Params::M_P*Params::m_O/(Params::k_B*6400);
    const fastreal n_0_day = 7.5e10;
    // nightside (originally from Nagy et al.,1981)
    const fastreal r_0_night = Params::R_V + 300e3;
    const fastreal beta_night = Params::G*Params::M_P*Params::m_O/(Params::k_B*4847);
    const fastreal n_0_night = 2e9;
    // densities
    fastreal n_day;
    fastreal n_night;
    // below exobase set density to zero
    if (rr < r_0_day) {
        n_day = 0;
    } else {
        n_day = n_0_day*exp(-beta_day*(1.0/r_0_day - 1.0/rr));
    }
    // below exobase set density to zero
    if (rr < r_0_night) {
        n_night = 0;
    } else {
        n_night = n_0_night*exp(-beta_night*(1.0/r_0_night - 1.0/rr));
    }
    // linear cos interpolation (noon -> midnight)
    return (1.0 - (sza/pi))*n_day + (sza/pi)*n_night;
}

//! Photoion profile: Titan
real SpatialDistribution::photoionChamberlainTitan(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 5) {
        ERRORMSG2("function takes five arguments (R_shadow,R_exo,n_exo,T_exo,ionizationrate)",this->name);
        doabort();
    }
    //R_shadow2 = (float) sqr(args[0]);
    //R_exo = args[1];
    //n_exo = args[2]; // #/m\B3
    //T_exo = args[3];
    //rate = args[4];
    if (args[1]<args[0]) {
        static bool errlogCnt = true;
        if (errlogCnt) {
            errorlog << "photoionChamberlainTitan: R_shadow > R_exo!" << endl;
            errlogCnt = false;
        }
    }
    const gridreal r2= sqr(r[0]) + sqr(r[1]) + sqr(r[2]);
    if (r2 < sqr(args[1])) return 0; //no emission below exobase.
    if (shadowTitan(r, args[0])==false) return 0;
    const real g_exo = Params::G*Params::M_P* Params::pops[this->popid]->m / (Params::k_B*args[3]);
    return  args[2]*args[4] * exp(g_exo*(1/sqrt(r2)-1/args[1]));
}

/** \brief Neutral density: Chamberlain distribution (temperature as an argument)
 *
 * n = sum_i n_i*exp( -G*M*m/(kB*T_i)*(1/r0 - 1/r) )
 *
 * arguments: neutralDensityChamberlainT r0 n1 T1 n2 T2 ... exobaseR
 */
real SpatialDistribution::neutralDensityChamberlainT(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    int N = this->args.size();
    if(N % 2 != 0 || N <= 0) {
        ERRORMSG2("function takes even amount arguments",this->name);
        doabort();
    }
    const real r0 = args[0];
    const real exobaseR = args[N-1];
    const fastreal rr = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
    // zero density below the exobase
    if( rr < exobaseR ) {
        return 0.0;
    }
    real ntot = 0.0;
    const real a = Params::G * Params::M_P * Params::pops[this->popid]->m / Params::k_B;
    for(int i=1; i<N-2; i+=2) {
        const real ni = args[i];
        const real Ti = args[i+1];
        const real beta = a/Ti;
        ntot += ni*exp( -beta*(1/r0 - 1/rr) );
    }
    return ntot;
}

/** \brief Neutral density: Chamberlain distribution (scale height as an argument)
 *
 * n = sum_i n_i*exp( -H_i*(1/r0 - 1/r) )
 *
 * arguments: neutralDensityChamberlainH r0 n1 H1 n2 H2 ... exobaseR
 */
real SpatialDistribution::neutralDensityChamberlainH(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    int N = this->args.size();
    if(N % 2 != 0 || N <= 0) {
        ERRORMSG2("function takes even amount arguments",this->name);
        doabort();
    }
    const real r0 = args[0];
    const real exobaseR = args[N-1];
    const fastreal rr = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
    // zero density below the exobase
    if( rr < exobaseR ) {
        return 0.0;
    }
    real ntot = 0.0;
    for(int i=1; i<N-2; i+=2) {
        const real ni = args[i];
        const real Hi = args[i+1];
        ntot += ni*exp( -Hi*(1/r0 - 1/rr) );
    }
    return ntot;
}

/** \brief Neutral density: power law
 *
 * n = sum_i n_i*(r0/r)^k_i
 *
 * arguments: neutralDensityPowerLaw r0 n1 k1 n2 k2 ... exobaseR
 */
real SpatialDistribution::neutralDensityPowerLaw(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    int N = this->args.size();
    if(N % 2 != 0 || N <= 0) {
        ERRORMSG2("function takes even amount arguments",this->name);
        doabort();
    }
    const real r0 = args[0];
    const real exobaseR = args[N-1];
    const fastreal rr = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
    // zero density below the exobase
    if( rr < exobaseR ) {
        return 0.0;
    }
    real ntot = 0.0;
    for(int i=1; i<N-2; i+=2) {
        const real ni = args[i];
        const real ki = args[i+1];
        ntot += ni*pow(r0/rr,ki);
    }
    return ntot;
}

/** \brief Neutral density: exponential
 *
 * n = sum_i n_i*exp( -(r-r0)/H_i )
 *
 * arguments: neutralDensityExponential r0 n1 H1 n2 H2 ... exobaseR
 */
real SpatialDistribution::neutralDensityExponential(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    int N = this->args.size();
    if(N % 2 != 0 || N <= 0) {
        ERRORMSG2("function takes even amount arguments",this->name);
        doabort();
    }
    const real r0 = args[0];
    const real exobaseR = args[N-1];
    const fastreal rr = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
    // zero density below the exobase
    if( rr < exobaseR ) {
        return 0.0;
    }
    real ntot = 0.0;
    for(int i=1; i<N-2; i+=2) {
        const real ni = args[i];
        const real Hi = args[i+1];
        ntot += ni*exp( -(rr-r0)/Hi );
    }
    return ntot;
}

//! Ionization profile: Constant everywhere
real SpatialDistribution::ionizationConstant(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes 1 argument",this->name);
        doabort();
    }
    const real ionizationRate = args[0];
    return ionizationRate;
}


//! Shadow profile: Optical shadow at r = shadowR
real SpatialDistribution::shadow(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes 1 argument",this->name);
        doabort();
    }
    const real shadowR = args[0];
    // return zero if in shadow
    if( r[0] < 0 && sqr(r[1])+sqr(r[2]) < sqr(shadowR) ) {
        return 0.0;
    } else {
        return 1.0;
    }
}

//! Shadow profile: Titan
bool SpatialDistribution::shadowTitan(const gridreal r[3], const real R_shadow)
{
    //alp=Saturn Local Time /24 *2*pi
    static const fastreal alp = Params::SaturnLocalTime/24*2*pi;
    static const fastreal ssl = Params::SubSolarLatitude/180*pi;
    // unit vector in the direction of the Sun - calculate only once
    static const gridreal e_Sun[] = {-sin(alp) * cos(ssl), -cos(alp) * cos(ssl), sin(ssl)};
    // distance to the Sun-Titan line
    gridreal dd[3];
    crossProduct(e_Sun,r,dd);
    // d=abs((rB-rA)x(rA-r))/abs(rB-rA), rA and rB on the line
    // here rA=[0,0,0] and rB=e_Sun
    //Plane through origin ax + by + cz = 0 perpendicular to e_Sun
    // is e[0]*x + e[1]*y + e[2]*z = 0
    if (sqr(dd[0]) + sqr(dd[1]) + sqr(dd[2]) < sqr(R_shadow) && (dotProduct(e_Sun,r) < 0)) return false;
    //shadow radius R_shadow (smaller than R_exo)
    else return true;
}

//! Default constructor
MultipleProductDistribution::MultipleProductDistribution()
{
    // Make sure the class is not used without
    // proper initialization (this aborts execution
    // if getValue called).
    this->ptr = &MultipleProductDistribution::defaultFunction;
}

//! Constructor
MultipleProductDistribution::MultipleProductDistribution(vector<SpatialDistribution> distFuncs)
{
    this->distFuncs = distFuncs;
    this->ptr = &MultipleProductDistribution::productFunction;
}

MultipleProductDistribution::~MultipleProductDistribution() { }

//! Wrapper function for the function pointer
real MultipleProductDistribution::getValue(const gridreal r[])
{
    return (this->*ptr)(r);
}

//! Returns value of the nth function in the vector
real MultipleProductDistribution::getValue(const gridreal r[],const unsigned int n)
{
    if(n < distFuncs.size()) {
        return distFuncs[n].getValue(r);
    } else {
        return -1.0;
    }
}

//! String summary
string MultipleProductDistribution::toString(string prefix,string delim)
{
    stringstream ss;
    for(unsigned int ii = 0; ii < distFuncs.size(); ii++) {
        ss << prefix << distFuncs[ii].toString() << delim;
    }
    return ss.str();
}

//! Default function, which aborts the program if called
real MultipleProductDistribution::defaultFunction(const gridreal r[3])
{
    ERRORMSG("function pointer not set");
    doabort();
    return -1;
}

//! Product function
real MultipleProductDistribution::productFunction(const gridreal r[3])
{
    real result = 1;
    for(unsigned int ii = 0; ii < distFuncs.size(); ii++) {
        result *= distFuncs[ii].getValue(r);
    }
    return result;
}

//! cos(SZA) at the point r
real CosSZA(const gridreal r[3])
{
    static const real SS_colat = 0.5*pi*(1 - (Params::SubSolarLatitude/180));
    static const real SS_phi = pi*(1.5 - (Params::SaturnLocalTime/12));
    static const gridreal SolarDir[3] = {
        static_cast<gridreal>(sin(SS_colat)*cos(SS_phi)),
        static_cast<gridreal>(sin(SS_colat)*sin(SS_phi)),
        static_cast<gridreal>(cos(SS_colat))
    };
    return dotProduct(SolarDir,r)/normvec(r);
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Ionospheric emission: Goes from SZA=0 to 90 as cos(SZA) and SZA>90 is constant
real SpatialDistribution::sph_ionoCosSzaDayConstantNight(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 4) {
        ERRORMSG2("function takes four arguments",this->name);
        doabort();
    }
    const real daySideFactor = args[0];
    const real nightSideFactor = args[1];
    const real dummyScaleHeight = args[2];
    const real R = args[3];
    fastreal sza;
    if (nightSideFactor < 0) {
        ERRORMSG2("nightSideFactor < 0",this->name);
        doabort();
    }
    if (daySideFactor < 0) {
        ERRORMSG2("daySideFactor < 0",this->name);
        doabort();
    }
    const fastreal rr = r[0];
    // We need to rename vector r in r_new
    gridreal r_new[3] = {r[0], r[1], r[2]};
    sph_transf_H2S_R(r_new);
    sph_transf_S2C_R(r_new);
    // Propagation along x-axis
    if (Params::sph_propagation_dir == 0) sza = acos(r_new[0]/rr);
    // Propagation along y-axis
    if (Params::sph_propagation_dir == 1) sza = acos(r_new[2]/rr);
    // No ions inside the radius
    if(rr < R) {
        return 0.0;
    }
    fastreal effdens;
    if(sza < pi/2) {
        effdens = daySideFactor + (nightSideFactor - daySideFactor)*(1-cos(sza));
    } else {
        effdens = nightSideFactor;
    }
    return effdens*exp(-(rr-R)/dummyScaleHeight);
}

//! (SPHERICAL) Neutral corona profile: Venus thermal+hot hydrogen
real SpatialDistribution::sph_neutralDensityVenusHydrogen(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes 1 argument",this->name);
        doabort();
    }
    const real exobaseR = args[0];
    fastreal sza = 0.0;
    const fastreal rr = r[0];
    // zero density below the exobase
    if( rr < exobaseR ) {
        return 0.0;
    }
    // We need to rename vector r in r_new
    gridreal r_new[3] = {r[0], r[1], r[2]};
    sph_transf_H2S_R(r_new);
    sph_transf_S2C_R(r_new);
    if (Params::sph_propagation_dir == 0) sza = acos(r_new[0]/rr);    //Propagation along x-axis
    else if (Params::sph_propagation_dir == 1) sza = acos(r_new[2]/rr);    //Propagation along y-axis
    // ---- Hot neutral hydrogen model from Gunell et al., Planetary and Space Science, 09/2005
    // (originally from Rodriguez et al. 1984)
    // zone: ZIP 1 = Dayside
    const fastreal n_0_zip1 = 1.32e11;
    const fastreal beta_zip1 = Params::G*Params::M_P*Params::m_p/(Params::k_B*285);
    // zone: ZIP 5 ~ Nightside
    const fastreal n_0_zip5 = 2.59e12;
    const fastreal beta_zip5 = Params::G*Params::M_P*Params::m_p/(Params::k_B*110);
    // densities
    const fastreal n_thermal_day = n_0_zip1*exp(-beta_zip1*(1/exobaseR - 1.0/rr));
    const fastreal n_thermal_night = n_0_zip5*exp(-beta_zip5*(1.0/exobaseR - 1.0/rr));
    // linear cos interpolation (noon -> midnight)
    const fastreal n_thermal = (1.0 - (sza/pi))*n_thermal_day + (sza/pi)*n_thermal_night;
    // ---- Hot neutral hydrogen model from Gunell et al., Planetary and Space Science, 09/2005
    // (fit to the analytical model by Rodriguez et al. 1984)
    // noon (sza = 0 deg)
    const fastreal a_1_noon = -6.2625e-5/1e3;
    const fastreal a_2_noon = 15.4817;
    const fastreal a_3_noon = 3.6414e4*1e3;
    // terminator (sza = 90 deg)
    const fastreal a_1_term = -8.4607e-5/1e3;
    const fastreal a_2_term = 15.9944;
    const fastreal a_3_term = 2.9743e4*1e3;
    // midnight (sza = 180 deg)
    const fastreal a_1_midn = -6.2309e-5/1e3;
    const fastreal a_2_midn = 15.2723;
    const fastreal a_3_midn = 4.3781e4*1e3;
    // densities
    const fastreal n_noon = exp(a_1_noon*rr + a_2_noon + a_3_noon/rr);
    const fastreal n_term = exp(a_1_term*rr + a_2_term + a_3_term/rr);
    const fastreal n_midn = exp(a_1_midn*rr + a_2_midn + a_3_midn/rr);
    // linear cos interpolation (noon -> terminator -> midnight)
    fastreal n_hot;
    if(sza <= pi/2) { // dayside hemisphere
        n_hot = (1.0 - (sza/(pi/2)))*n_noon + (sza/(pi/2))*n_term;
    } else { // nightside hemisphere
        n_hot = (1.0 - (sza/(pi/2) - 1))*n_term + (sza/(pi/2) - 1)*n_midn;
    }
    return n_thermal + n_hot;
}

//! (SPHERICAL) Neutral corona profile: Venus hot oxygen
real SpatialDistribution::sph_neutralDensityVenusOxygenHot(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes 1 argument",this->name);
        doabort();
    }
    const real exobaseR = args[0];
    fastreal sza = 0.0;
    const fastreal rr = r[0];
    // zero density below the exobase
    if( rr < exobaseR ) {
        return 0.0;
    }
    // We need to rename vector r in r_new
    gridreal r_new[3] = {r[0], r[1], r[2]};
    sph_transf_H2S_R(r_new);
    sph_transf_S2C_R(r_new);
    if (Params::sph_propagation_dir == 0) sza = acos(r_new[0]/rr);    //Propagation along x-axis
    else if(Params::sph_propagation_dir == 1) sza = acos(r_new[2]/rr);    //Propagation along y-axis
    // ---- Hot neutral oxygen model from Gunell et al., Planetary and Space Science, 09/2005
    // dayside (originally from McElroy et al., 1982)
    const fastreal r_0_day = Params::R_V + 200e3;
    const fastreal beta_day = Params::G*Params::M_P*Params::m_O/(Params::k_B*6400);
    const fastreal n_0_day = 7.5e10;
    // nightside (originally from Nagy et al.,1981)
    const fastreal r_0_night = Params::R_V + 300e3;
    const fastreal beta_night = Params::G*Params::M_P*Params::m_O/(Params::k_B*4847);
    const fastreal n_0_night = 2e9;
    // densities
    fastreal n_day;
    fastreal n_night;
    // below exobase set density to zero
    if (rr < r_0_day) {
        n_day = 0;
    } else {
        n_day = n_0_day*exp(-beta_day*(1.0/r_0_day - 1.0/rr));
    }
    // below exobase set density to zero
    if (rr < r_0_night) {
        n_night = 0;
    } else {
        n_night = n_0_night*exp(-beta_night*(1.0/r_0_night - 1.0/rr));
    }
    // linear cos interpolation (noon -> midnight)
    return (1.0 - (sza/pi))*n_day + (sza/pi)*n_night;
}

//! (SPHERICAL) Shadow profile: Optical shadow
real SpatialDistribution::sph_shadow(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes 1 argument",this->name);
        doabort();
    }
    const real shadowR = args[0];
    // return zero if in shadow
    gridreal r_new[3] = {r[0], r[1], r[2]};
    sph_transf_H2S_R(r_new);
    sph_transf_S2C_R(r_new);
    if(Params::sph_propagation_dir == 0) {
        if( r_new[0] < 0 && sqr(r_new[1])+sqr(r_new[2]) < sqr(shadowR) ) {
            return 0.0;
        } else {
            return 1.0;
        }
    } else { // (Params::sph_propagation_dir == 1)
        if( r_new[2] < 0 && sqr(r_new[0])+sqr(r_new[1]) < sqr(shadowR) ) {
            return 0.0;
        } else {
            return 1.0;
        }
    }
}

#endif

