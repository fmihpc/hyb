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
#  pragma implementation "variables.H"
#endif

#include "variables.H"
#include "constants.H"
#include "gridcache.H"
#include <string.h>
#include <math.h>

enum {VAR_RHO,
	  VAR_N,
	  VAR_RHOVX,VAR_RHOVY,VAR_RHOVZ,
          VAR_NV,
	  VAR_VX,VAR_VY,VAR_VZ,VAR_VR,
	  VAR_V, VAR_V2,
	  VAR_U, VAR_U1, VAR_UMINUSDIP, VAR_P, VAR_T, VAR_A, VAR_VA, VAR_VF,
	  VAR_MACH, VAR_ALFVENMACH, VAR_BETA, VAR_WORK,
	  VAR_FRACT_THERMAL, VAR_FRACT_KINETIC, VAR_FRACT_MAGNETIC,
	  VAR_BX,VAR_BY,VAR_BZ,
	  VAR_BX0,VAR_BY0,VAR_BZ0,
	  VAR_BX1,VAR_BY1,VAR_BZ1,
	  VAR_DIVB,VAR_DIVB_REL,
	  VAR_B,VAR_B0,VAR_BMINUSB0,VAR_B1, VAR_B2,VAR_B02,VAR_B12, VAR_LARMOR_P, VAR_LARMOR_P_OVER_DX,
	  VAR_LARMOR_E,
	  VAR_INERTLEN_P,VAR_INERTLEN_E,
	  VAR_GYRO_P,VAR_GYRO_E,
	  VAR_LOWERHYBRID,
	  VAR_PLASMAFREQ,
	  VAR_DEBYELEN,
	  VAR_EX, VAR_EY, VAR_EZ, VAR_E,
	  VAR_SX, VAR_SY, VAR_SZ, VAR_S,		// Poynting vector (1/mu0)*ExB
	  VAR_KX, VAR_KY, VAR_KZ, VAR_K,		// Energy flux vector K = (u+P-B^2/(2mu0))v + S
      VAR_JX,VAR_JY,VAR_JZ,VAR_J,VAR_J_OVER_ENV, VAR_JPAR};
#define VAR_LAST (VAR_JPAR+1)

const static char *varnames[VAR_LAST];
const static char *vardescr[VAR_LAST];

Tvariable::Tvariable()
{
	var = 0;
	logflag = false;
	pseudobackgroundflag = false;
	static bool FirstTime = true;
	if (FirstTime) {
		varnames[VAR_RHO] = "rho";
		varnames[VAR_N] = "n";
		varnames[VAR_RHOVX] = "rhovx";
		varnames[VAR_RHOVY] = "rhovy";
		varnames[VAR_RHOVZ] = "rhovz";
	        varnames[VAR_NV] = "nv";
		varnames[VAR_VX] = "vx";
		varnames[VAR_VY] = "vy";
		varnames[VAR_VZ] = "vz";
		varnames[VAR_VR] = "vr";
		varnames[VAR_V2] = "v2";
		varnames[VAR_V] = "v";
		varnames[VAR_U] = "U";
		varnames[VAR_U1] = "U1";
		varnames[VAR_UMINUSDIP] = "Uminusdip";
		varnames[VAR_P] = "P";
		varnames[VAR_T] = "T";
		varnames[VAR_A] = "a";
		varnames[VAR_VA] = "vA";
		varnames[VAR_VF] = "vf";
		varnames[VAR_MACH] = "Ms";
		varnames[VAR_ALFVENMACH] = "MA";
		varnames[VAR_BETA] = "beta";
		varnames[VAR_WORK] = "work";
		varnames[VAR_FRACT_THERMAL] = "Rth";
		varnames[VAR_FRACT_KINETIC] = "Rkin";
		varnames[VAR_FRACT_MAGNETIC] = "Rmag";
		varnames[VAR_BX] = "Bx";
		varnames[VAR_BY] = "By";
		varnames[VAR_BZ] = "Bz";
		varnames[VAR_BX0] = "Bx0";
		varnames[VAR_BY0] = "By0";
		varnames[VAR_BZ0] = "Bz0";
		varnames[VAR_BX1] = "Bx1";
		varnames[VAR_BY1] = "By1";
		varnames[VAR_BZ1] = "Bz1";
		varnames[VAR_B] = "B";
		varnames[VAR_B0] = "B0";
		varnames[VAR_BMINUSB0] = "BminusB0";
		varnames[VAR_B1] = "B1";
		varnames[VAR_B2] = "B2";
		varnames[VAR_B02] = "B02";
		varnames[VAR_B12] = "B12";
		varnames[VAR_DIVB] = "divB";
		varnames[VAR_DIVB_REL] = "divBrel";
		varnames[VAR_LARMOR_P] = "rLp";
		varnames[VAR_LARMOR_P_OVER_DX] = "rLp/dx";
		varnames[VAR_LARMOR_E] = "rLe";
		varnames[VAR_INERTLEN_P] = "c/wpi";
		varnames[VAR_INERTLEN_E] = "c/wpe";
		varnames[VAR_GYRO_P] = "fgyrop";
		varnames[VAR_GYRO_E] = "fgyroe";
		varnames[VAR_LOWERHYBRID] = "fLH";
		varnames[VAR_PLASMAFREQ] = "fpe";
		varnames[VAR_DEBYELEN] = "rDebye";
		varnames[VAR_JX] = "jx";
		varnames[VAR_JY] = "jy";
		varnames[VAR_JZ] = "jz";
		varnames[VAR_J] = "j";
		varnames[VAR_EX] = "Ex";
		varnames[VAR_EY] = "Ey";
		varnames[VAR_EZ] = "Ez";
		varnames[VAR_E] = "E";
		varnames[VAR_SX] = "Sx";
		varnames[VAR_SY] = "Sy";
		varnames[VAR_SZ] = "Sz";
		varnames[VAR_S] = "S";
		varnames[VAR_KX] = "Kx";
		varnames[VAR_KY] = "Ky";
		varnames[VAR_KZ] = "Kz";
		varnames[VAR_K] = "K";
		varnames[VAR_J_OVER_ENV] = "j/env";
		varnames[VAR_JPAR] = "jPar";
		vardescr[VAR_RHO] = "Mass density (kg/m^3), Primary variable";
		vardescr[VAR_N] = "Number density (1/m^3), Computed from rho as n=rho/mass (mass=particle mass)";
		vardescr[VAR_RHOVX] = "X-component of the momentum flux (kg/(m^2*s)), Primary variable";
		vardescr[VAR_RHOVY] = "Y-component of the momentum flux (kg/(m^2*s)), Primary variable";
		vardescr[VAR_RHOVZ] = "Z-component of the momentum flux (kg/(m^2*s)), Primary variable";
	        vardescr[VAR_NV] = "Magnitude of the particle flux (1/(m^2*s)), Computed as nv=sqrt(rhovx^2 + rhovy^2 + rhovz^2)/mass";
		vardescr[VAR_VX] = "X-component of the velocity (m/s), Computed as vx=rhovx/rho";
		vardescr[VAR_VY] = "Y-component of the velocity (m/s), Computed as vy=rhovy/rho";
		vardescr[VAR_VZ] = "Z-component of the velocity (m/s), Computed as vz=rhovz/rho";
		vardescr[VAR_VR] = "Radial component of the velocity (m/s)";
		vardescr[VAR_V2] = "Velocity squared, v2 = vx^2 + vy^2 + vz^2";
		vardescr[VAR_V] = "Velocity, v = sqrt(vx^2 + vy^2 + vz^2)";
		vardescr[VAR_U] = "Total energy density, thermal + kinetic + magnetic (J/m^3)";
		vardescr[VAR_U1] = "Total energy density excluding dipole field (J/m^3), Primary variable";
		vardescr[VAR_UMINUSDIP] = "Total energy density minus dipole field energy (J/m^3)";
		vardescr[VAR_P] = "Pressure (J/m^3), Computed as P = (gamma-1)*(U - 0.5*rho*v2 - 0.5*B12/mu0)";
		vardescr[VAR_T] = "Temperature (K), Computed as T = P/(kB*n)";
		vardescr[VAR_A] = "Sound speed (m/s), Computed as a = sqrt(gamma*P/rho)";
		vardescr[VAR_VA] = "Alfven speed (m/s), Computed as vA = B/sqrt(mu0*rho)";
		vardescr[VAR_VF] = "Fast mode speed (m/s), vf = sqrt(a^2 + vA^2)";
		vardescr[VAR_MACH] = "Sonic Mach number, Ms = v/a";
		vardescr[VAR_ALFVENMACH] = "Alfvenic Mach number, MA = v/vA";
		vardescr[VAR_BETA] = "Plasma beta parameter, beta = 2*mu0*P/B^2";
		vardescr[VAR_WORK] = "Computational load taking subcycling into account, work = vf/dx^(dim+1)\n"
			"where dx is the local grid spacing and dim is the grid dimensionality.\n"
			"Dimension: Cell updates per second per volume\n"
			"(or area, or length if 2D or 1D grid)";
		vardescr[VAR_FRACT_THERMAL] = "Thermal energy density / total energy density, Rth=P/((gamma-1)*U)";
		vardescr[VAR_FRACT_KINETIC] = "Kinetic energy density / total energy density, Rkin=0.5*rho*v2/U";
		vardescr[VAR_FRACT_MAGNETIC] = "Magnetic energy density / total energy density, Rmag=0.5*B12/(mu0*U)";
		vardescr[VAR_BX] = "X-component of the total magnetic field (T), Computed as Bx=Bx0+Bx1";
		vardescr[VAR_BY] = "Y-component of the total magnetic field (T), Computed as By=By0+By1";
		vardescr[VAR_BZ] = "Z-component of the total magnetic field (T), Computed as Bz=Bz0+Bz1";
		vardescr[VAR_BX0] = "X-component of the background magnetic field (T), Primary quantity";
		vardescr[VAR_BY0] = "Y-component of the background magnetic field (T), Primary quantity";
		vardescr[VAR_BZ0] = "Z-component of the background magnetic field (T), Primary quantity";
		vardescr[VAR_BX1] = "X-component of the variation magnetic field (T), Primary quantity";
		vardescr[VAR_BY1] = "Y-component of the variation magnetic field (T), Primary quantity";
		vardescr[VAR_BZ1] = "Z-component of the variation magnetic field (T), Primary quantity";
		vardescr[VAR_B] = "Magnitude of the total magnetic field (T), B = sqrt(Bx^2 + By^2 + Bz^2)";
		vardescr[VAR_B0] = "Magnitude of the background magnetic field (T),\n"
			"B0 = sqrt(Bx0^2 + By0^2 +Bz0^2)";
		vardescr[VAR_BMINUSB0] = "|B| - |B0| (T)";
		vardescr[VAR_B1] = "Magnitude of the variation magnetic field (T),\n"
			"B1 = sqrt(Bx1^2 + By1^2 + Bz1^2)";
		vardescr[VAR_B2] = "Square of the total magnetic field (T^2),\n"
			"B2 = Bx^2 + By^2 + Bz^2 = B^2";
		vardescr[VAR_B02] = "Square of the background magnetic field (T^2),\n"
			"B02 = Bx0^2 + By0^2 + Bz0^2 = B0^2";
		vardescr[VAR_B12] = "Square of the variation magnetic field (T^2),\n"
			"B12 = Bx1^2 + By1^2 + Bz1^2 = B1^2";
		vardescr[VAR_DIVB] = "Divergence of the variation magnetic field (T/m),\n"
			"divB = ddx(Bx1) + ddy(By1) + ddz(Bz1)";
		vardescr[VAR_DIVB_REL] = "Normalized divergence of the variation magnetic field,\n"
			"divBrel=divB*(dx/B1), where dx is the local grid spacing";
		vardescr[VAR_LARMOR_P] = "Proton Larmor radius (m)";
		vardescr[VAR_LARMOR_P_OVER_DX] = "Proton Larmor radius divided by local grid spacing";
		vardescr[VAR_LARMOR_E] = "Electron Larmor radius (m), computed assuming Ti=4*Te";
		vardescr[VAR_INERTLEN_P] = "Proton inertial length (m, c/omega_pi)";
		vardescr[VAR_INERTLEN_E] = "Electron inertial length (m, c/omega_pe)";
		vardescr[VAR_GYRO_P] = "Proton gyrofrequency (Hz)";
		vardescr[VAR_GYRO_E] = "Electron gyrofrequency (Hz)";
		vardescr[VAR_LOWERHYBRID] = "Lower hybrid frequency (Hz)";
		vardescr[VAR_PLASMAFREQ] = "(Electron) plasma frequency (Hz)";
		vardescr[VAR_DEBYELEN] = "Debye length (m), computed assuming Ti=4*Te";
		vardescr[VAR_JX] = "X-component of the current density (A/m^2), jx = (ddy(Bz1) - ddz(By1))/mu0";
		vardescr[VAR_JY] = "Y-component of the current density (A/m^2), jy = (ddz(Bx1) - ddx(Bz1))/mu0";
		vardescr[VAR_JZ] = "Z-component of the current density (A/m^2), jz = (ddx(By1) - ddy(Bx1))/mu0";
		vardescr[VAR_J] = "Magnitude of the current density (A/m^2), j = sqrt(jx^2 + jy^2 + jz^2)";
		vardescr[VAR_J_OVER_ENV] = "j/(e*n*v)";
		vardescr[VAR_JPAR] = "Parallel current density (A/m^2), jPar = (Bx*jx + By*jy + Bz*jz)/B";
		vardescr[VAR_EX] = "X-component of the B x v electric field (V/m^2)";
		vardescr[VAR_EY] = "Y-component of the B x v electric field (V/m^2)";
		vardescr[VAR_EZ] = "Z-component of the B x v electric field (V/m^2)";
		vardescr[VAR_E] = "Magnitude of the B x v electric field (V/m^2)";
		vardescr[VAR_SX] = "X-component of the Poynting vector (W/m^2)";
		vardescr[VAR_SY] = "Y-component of the Poynting vector (W/m^2)";
		vardescr[VAR_SZ] = "Z-component of the Poynting vector (W/m^2)";
		vardescr[VAR_S] = "Magnitude of the Poynting vector (W/m^2)";
		vardescr[VAR_KX] = "X-component of energy flux vector (W/m^2)";
		vardescr[VAR_KY] = "Y-component of energy flux vector (W/m^2)";
		vardescr[VAR_KZ] = "Z-component of energy flux vector (W/m^2)";
		vardescr[VAR_K] = "Magnitude of energy flux vector (W/m^2)";
		FirstTime = false;
	}
	mintab = new double [VAR_LAST];
	maxtab = new double [VAR_LAST];
	minmaxexiststab = new bool [VAR_LAST];
	int i;
	for (i=0; i<VAR_LAST; i++) {
		mintab[i] = 0;
		maxtab[i] = -1;
		minmaxexiststab[i] = false;
	}
}

int Tvariable::Nvars() const {return VAR_LAST;}

Tvariable::~Tvariable()
{
	delete [] minmaxexiststab;
	delete [] maxtab;
	delete [] mintab;
}

void Tvariable::clearAllMinMax()
{
	int i;
	for (i=0; i<VAR_LAST; i++) minmaxexiststab[i] = false;
}

void Tvariable::getMinMax(double& m, double& M) const
{
	if (!minmaxexiststab[var]) {
		cerr << "*** Tvariable::getMinMax: min and max not defined for var=" << var << "\n";
		return;
	}
	m = mintab[var];
	M = maxtab[var];
}

void Tvariable::select(int i)
{
	if (i < 0 || i > VAR_LAST) {cerr << "*** Tvariable::select(i=" << i << "), out of range\n"; return;}
	var = i;
}

bool Tvariable::select(const char *varname, double gamma1, double invmu01, double mass1)
{
	int i;
	for (i=0; i<VAR_LAST; i++)
		if (!strcmp(varname,varnames[i])) {
			var = i;
			gamma = gamma1;
			invmu0 = invmu01;
		        mass = mass1;
			if (invmu0 == 0) cerr << "*** Tvariable::select: gamma1=" << gamma1 << ", invmu01=" << invmu01 << ", mass1=" << mass1 << "\n";
			return true;
		}
	cerr << "*** Tvariable::select: Unknown variable name '" << varname << "', unchanged\n";
   return false;
}

bool Tvariable::check(const char *varname)
{
   for (int i=0; i<VAR_LAST; i++)
     {	
	if (!strcmp(varname,varnames[i])) { return true; }
     }
   return false;
}

const char *Tvariable::selected() const
{
	return varnames[var];
}

const char *Tvariable::description() const
{
	return vardescr[var];
}

extern real GridDimensionScaling;	// define in gridcache.C

void ComputeCurl(Tmetagrid& g, const Tdimvec& X, smallnat a0, real& jx, real& jy, real& jz)
// jx = ddy(Bz) - ddz(By)
// jy = ddz(Bx) - ddx(Bz)
// jz = ddx(By) - ddy(Bx)
//
//
//                     6    2
//                     |   /
//                     | /
//              4------0----- 3
//                   / |
//                 /   |
//               1     5
{
	// Save old g.CT() table
	// ---------------------
	int a;
	double CTold[MAX_VARS];
	const int ncd = g.Ncelldata();
	for (a=0; a<ncd; a++) CTold[a] = g.CT(a);
	// ---------------------
	const TGridIndex c = g.find(X);
	real Bx[7],By[7],Bz[7];
	Tdimvec X1;
	const int order = 1;
	if (!g.intpol(X,order,true)) {
		jx = 0;
		jy = 0;
		jz = 0;
		// ----- restore CT table ---------------
		for (a=0; a<ncd; a++) g.CT(a) = CTold[a];
		// --------------------------------------
		return;
	} else {
		Bx[0] = g.CT(a0);
		By[0] = g.CT(a0+1);
		Bz[0] = g.CT(a0+2);
	}
	real Xpts[7][3];
	// Fill in points
	const real dx = g.cellsize(c);
	Xpts[1][0] = X(0)-dx; Xpts[1][1] = X(1);    Xpts[1][2] = X(2);
	Xpts[2][0] = X(0)+dx; Xpts[2][1] = X(1);    Xpts[2][2] = X(2);
	Xpts[3][0] = X(0);    Xpts[3][1] = X(1)-dx; Xpts[3][2] = X(2);
	Xpts[4][0] = X(0);    Xpts[4][1] = X(1)+dx; Xpts[4][2] = X(2);
	Xpts[5][0] = X(0);    Xpts[5][1] = X(1);    Xpts[5][2] = X(2)-dx;
	Xpts[6][0] = X(0);    Xpts[6][1] = X(1);    Xpts[6][2] = X(2)+dx;

	// Compute Bx,By,Bz at the 6+1 stencil points
	int pt;
	for (pt=1; pt<7; pt++) {
		X1[0] = Xpts[pt][0];
		X1[1] = Xpts[pt][1];
		X1[2] = Xpts[pt][2];
		if (g.intpol(X1,order,true)) {
			Bx[pt] = g.CT(a0);
			By[pt] = g.CT(a0+1);
			Bz[pt] = g.CT(a0+2);
		} else {
			// If the point is outside domain, use the central magnetic field.
			// This mimics Neumann boundary condition for B.
			Bx[pt] = Bx[0];
			By[pt] = By[0];
			Bz[pt] = Bz[0];
		}
	}

	// GridDimensionScaling is 1/R_E, or something similar.
	// (if hcvis -scale R_E was used).
	const real inv2dx = GridDimensionScaling/(2*dx);
    // jx = ddy(Bz) - ddz(By)
    // jy = ddz(Bx) - ddx(Bz)
    // jz = ddx(By) - ddy(Bx)
	jx = ((Bz[4] - Bz[3]) - (By[6] - By[5]))*inv2dx;
	jy = ((Bx[6] - Bx[5]) - (Bz[2] - Bz[1]))*inv2dx;
	jz = ((By[2] - By[1]) - (Bx[4] - Bx[3]))*inv2dx;
	// ----- restore CT table ---------------
	for (a=0; a<ncd; a++) g.CT(a) = CTold[a];
	// --------------------------------------
}

real ComputeDiv(Tmetagrid& g, const Tdimvec& X, smallnat a0, bool absflag=false)
// divB = ddx(Bx) + ddy(By) + ddz(Bz)
//
//
//                     6    2
//                     |   /
//                     | /
//              4------0----- 3
//                   / |
//                 /   |
//               1     5
{
	// Save old g.CT() table
	// ---------------------
	int a;
	double CTold[MAX_VARS];
	const int ncd = g.Ncelldata();
	for (a=0; a<ncd; a++) CTold[a] = g.CT(a);
	// ---------------------
	const TGridIndex c = g.find(X);
	real Bx[7],By[7],Bz[7];
	Tdimvec X1;
	const int order = 1;
	if (!g.intpol(X,order,true)) {
		// ----- restore CT table ---------------
		for (a=0; a<ncd; a++) g.CT(a) = CTold[a];
		// --------------------------------------
		return 0;
	} else {
		Bx[0] = g.CT(a0);
		By[0] = g.CT(a0+1);
		Bz[0] = g.CT(a0+2);
	}
	real Xpts[7][3];
	// Fill in points
	const real dx = g.cellsize(c);
	Xpts[1][0] = X(0)-dx; Xpts[1][1] = X(1);    Xpts[1][2] = X(2);
	Xpts[2][0] = X(0)+dx; Xpts[2][1] = X(1);    Xpts[2][2] = X(2);
	Xpts[3][0] = X(0);    Xpts[3][1] = X(1)-dx; Xpts[3][2] = X(2);
	Xpts[4][0] = X(0);    Xpts[4][1] = X(1)+dx; Xpts[4][2] = X(2);
	Xpts[5][0] = X(0);    Xpts[5][1] = X(1);    Xpts[5][2] = X(2)-dx;
	Xpts[6][0] = X(0);    Xpts[6][1] = X(1);    Xpts[6][2] = X(2)+dx;

	// Compute Bx,By,Bz at the 6+1 stencil points
	int pt;
	for (pt=1; pt<7; pt++) {
		X1[0] = Xpts[pt][0];
		X1[1] = Xpts[pt][1];
		X1[2] = Xpts[pt][2];
		if (g.intpol(X1,order,true)) {
			Bx[pt] = g.CT(a0);
			By[pt] = g.CT(a0+1);
			Bz[pt] = g.CT(a0+2);
		} else {
			// If the point is outside domain, use the central magnetic field.
			// This mimics Neumann boundary condition for B.
			Bx[pt] = Bx[0];
			By[pt] = By[0];
			Bz[pt] = Bz[0];
		}
	}

	// GridDimensionScaling is 1/R_E, or something similar.
	// (if hcvis -scale R_E was used).
	const real inv2dx = GridDimensionScaling/(2*dx);
    // jx = ddy(Bz) - ddz(By)
    // jy = ddz(Bx) - ddx(Bz)
    // jz = ddx(By) - ddy(Bx)
	// ----- restore CT table ---------------
	for (a=0; a<ncd; a++) g.CT(a) = CTold[a];
	// --------------------------------------
	if (absflag) {
		return (max(fabs(Bx[2]),fabs(Bx[1])) + max(fabs(By[4]),fabs(By[3])) + max(fabs(Bz[6]),fabs(Bz[5])))*inv2dx;
	} else
		return ((Bx[2] - Bx[1]) + (By[4] - By[3]) + (Bz[6] - Bz[5]))*inv2dx;
}


double Tvariable::get(Tmetagrid& g, const Tdimvec& X) const
{
	double result = 0;
	const int ncd = g.Ncelldata();
	switch (var) {
	case VAR_RHO: result = g.CT(0); break;
	case VAR_N: result = g.CT(0)/mass; break;
	case VAR_RHOVX: result = g.CT(1); break;
	case VAR_RHOVY: result = g.CT(2); break;
	case VAR_RHOVZ: result = g.CT(3); break;
	case VAR_NV: result = sqrt( sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)) )/mass; break;
	case VAR_VX: result = g.CT(1)/g.CT(0); break;
	case VAR_VY: result = g.CT(2)/g.CT(0); break;
	case VAR_VZ: result = g.CT(3)/g.CT(0); break;
	case VAR_VR: result = (g.CT(1)*X(0) + g.CT(2)*X(1) + g.CT(3)*X(2))/(g.CT(0)*sqrt(sqr(X(0))+sqr(X(1))+sqr(X(2)))); break;
	case VAR_V2: result = (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/sqr(g.CT(0)); break;
	case VAR_V: result = sqrt((sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/sqr(g.CT(0))); break;
	case VAR_U:
		if (pseudobackgroundflag) {
			result = g.CT(4);
		} else {
			result = g.CT(4) + 0.5*invmu0*(sqr(g.CT(8)) + sqr(g.CT(9)) + sqr(g.CT(10)))
				+ (g.CT(5)*g.CT(8) + g.CT(6)*g.CT(9) + g.CT(7)*g.CT(10))*invmu0;
		}
		break;
	case VAR_U1:
		if (pseudobackgroundflag) {
			result = g.CT(4) - 0.5*invmu0*(sqr(g.CT(8)) + sqr(g.CT(9)) + sqr(g.CT(10)))
				- (g.CT(5)*g.CT(8) + g.CT(6)*g.CT(9) + g.CT(7)*g.CT(10))*invmu0;
		} else
			result = g.CT(4);
		break;
	case VAR_UMINUSDIP:
		if (pseudobackgroundflag) {
			result = g.CT(4);
			cout << "warning: Uminusdip not implemented correctly when pseudobackgroundflag is true\n";
		} else {
			result = g.CT(4) + (g.CT(5)*g.CT(8) + g.CT(6)*g.CT(9) + g.CT(7)*g.CT(10))*invmu0;
		}
		break;
	case VAR_P:
	case VAR_T:
	case VAR_BETA:
	case VAR_FRACT_THERMAL:
	case VAR_FRACT_KINETIC:
	case VAR_FRACT_MAGNETIC:
	{
		// First, compute P and energy density in U
		real U,P;
		if (ncd == 5) {
			U = g.CT(4);
			P = (gamma-1)*(U - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/(2*g.CT(0)));
		} else {
			U = g.CT(4);
			if (pseudobackgroundflag)
				U -= 0.5*invmu0*(sqr(g.CT(8)) + sqr(g.CT(9)) + sqr(g.CT(10)))
					+ (g.CT(5)*g.CT(8) + g.CT(6)*g.CT(9) + g.CT(7)*g.CT(10))*invmu0;
			P = (gamma-1)*(U - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/(2*g.CT(0))
						   - 0.5*invmu0*(sqr(g.CT(5)) + sqr(g.CT(6)) + sqr(g.CT(7))));
		}
		switch (var) {
		case VAR_P:
			result = P;
			break;
		case VAR_T:
		{
			const double n = g.CT(0)/mass;
		        result = P/(cnst::kB*n);
		}
		break;
		case VAR_BETA:
			if (ncd >= 8) {
				const double Bx1 = g.CT(5);
				const double By1 = g.CT(6);
				const double Bz1 = g.CT(7);
				double Bx0,By0,Bz0;
				if (ncd >= 11) {
					Bx0 = g.CT(8);
					By0 = g.CT(9);
					Bz0 = g.CT(10);
				} else {
					Bx0 = By0 = Bz0 = 0;
				}
				const double Bx = Bx0 + Bx1;
				const double By = By0 + By1;
				const double Bz = Bz0 + Bz1;
				const double B2 = sqr(Bx) + sqr(By) + sqr(Bz);
				result = 2*P/(invmu0*B2);
			} else
				result = 0;
			break;
		case VAR_FRACT_THERMAL:
			result = P/((gamma-1)*U);
			break;
		case VAR_FRACT_KINETIC:
		{
			const double rho = g.CT(0);
			const double p2 = sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3));
			result = 0.5*p2/(rho*U);
		}
		break;
		case VAR_FRACT_MAGNETIC:
		{
			const double B12 = sqr(g.CT(5)) + sqr(g.CT(6)) + sqr(g.CT(7));
			result = 0.5*B12*invmu0/U;
		}
		break;
		}
	}
	break;
	case VAR_MACH:		// Ms = v/a, a=sqrt(gamma*P/rho)
	case VAR_A:
	{
		const double rho = g.CT(0);
		const double v = sqrt((sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/sqr(rho));
		const int ncd = g.Ncelldata();
		double P;
		if (ncd == 5) {
			P = (gamma-1)*(g.CT(4) - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/(2*g.CT(0)));
		} else {
			real U = g.CT(4);
			if (pseudobackgroundflag)
				U -= 0.5*invmu0*(sqr(g.CT(8)) + sqr(g.CT(9)) + sqr(g.CT(10)))
					+ (g.CT(5)*g.CT(8) + g.CT(6)*g.CT(9) + g.CT(7)*g.CT(10))*invmu0;
			P = (gamma-1)*(U - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/(2*g.CT(0))
						   - 0.5*invmu0*(sqr(g.CT(5)) + sqr(g.CT(6)) + sqr(g.CT(7))));
		}
		const double a = sqrt(gamma*max(0.0,P)/max(0.0,rho));
		if (var == VAR_A)
			result = a;
		else
			result = v/a;
	}
	break;
	case VAR_ALFVENMACH:	// MA = v/vA, vA=|B|/sqrt(mu0*rho)
	case VAR_VA:
	case VAR_VF:
	case VAR_WORK:
	{
		const double rho = g.CT(0);
		const double v = sqrt((sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/sqr(rho));
		const int ncd = g.Ncelldata();
		if (ncd < 8) {
			if (var == VAR_WORK) {
				const double P = (gamma-1)*(g.CT(4) - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/(2*g.CT(0)));
				const double a = sqrt(gamma*max(0.0,P)/max(0.0,rho));
				const TGridIndex i = g.find(X);
				if (i == NOINDEX)
					result = 0;
				else
					result = a/pow(g.cellsize(i)/GridDimensionScaling,real(g.dimension()+1));
			} else {
				result = 0;
			}
			break;
		}
		const double Bx1 = g.CT(5);
		const double By1 = g.CT(6);
		const double Bz1 = g.CT(7);
		double Bx0,By0,Bz0;
		if (ncd >= 11) {
			Bx0 = g.CT(8);
			By0 = g.CT(9);
			Bz0 = g.CT(10);
		} else {
			Bx0 = By0 = Bz0 = 0;
		}
		const double Bx = Bx0 + Bx1;
		const double By = By0 + By1;
		const double Bz = Bz0 + Bz1;
		const double B2 = sqr(Bx) + sqr(By) + sqr(Bz);
		const double vA = sqrt(B2*invmu0/rho);
		switch (var) {
		case VAR_ALFVENMACH:
			result = v/vA;
			break;
		case VAR_VA:
			result = vA;
			break;
		case VAR_VF:
		case VAR_WORK:
		default:
		{
			const double rho = g.CT(0);
			const int ncd = g.Ncelldata();
			double P;
			if (ncd == 5) {
				P = (gamma-1)*(g.CT(4) - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/(2*g.CT(0)));
			} else {
				real U = g.CT(4);
				if (pseudobackgroundflag)
					U -= 0.5*invmu0*(sqr(g.CT(8)) + sqr(g.CT(9)) + sqr(g.CT(10)))
						+ (g.CT(5)*g.CT(8) + g.CT(6)*g.CT(9) + g.CT(7)*g.CT(10))*invmu0;
				P = (gamma-1)*(U - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/(2*g.CT(0))
							   - 0.5*invmu0*(sqr(g.CT(5)) + sqr(g.CT(6)) + sqr(g.CT(7))));
			}
			const double a = sqrt(gamma*max(0.0,P)/max(0.0,rho));
			const double vf = sqrt(sqr(a) + sqr(vA));
			if (var == VAR_VF)
				result = vf;
			else {
				const TGridIndex i = g.find(X);
				if (i == NOINDEX)
					result = 0;
				else
					result = vf/pow(g.cellsize(i)/GridDimensionScaling,real(g.dimension()+1));
			}
		}
		break;
		}
	}
	break;
	case VAR_BX:
	case VAR_BY:
	case VAR_BZ:
	case VAR_BX0:
	case VAR_BY0:
	case VAR_BZ0:
	case VAR_BX1:
	case VAR_BY1:
	case VAR_BZ1:
	case VAR_B:
	case VAR_B0:
	case VAR_BMINUSB0:
	case VAR_B1:
	case VAR_B2:
	case VAR_B02:
	case VAR_B12:
	case VAR_LARMOR_P:
	case VAR_LARMOR_P_OVER_DX:
	case VAR_LARMOR_E:
	case VAR_INERTLEN_P:
	case VAR_INERTLEN_E:
	case VAR_GYRO_P:
	case VAR_GYRO_E:
	case VAR_LOWERHYBRID:
	case VAR_PLASMAFREQ:
	case VAR_DEBYELEN:
	case VAR_DIVB:
	case VAR_DIVB_REL:
	case VAR_EX:
	case VAR_EY:
	case VAR_EZ:
	case VAR_E:
	case VAR_SX:
	case VAR_SY:
	case VAR_SZ:
	case VAR_S:
	case VAR_KX:
	case VAR_KY:
	case VAR_KZ:
	case VAR_K:
	case VAR_JX:
	case VAR_JY:
	case VAR_JZ:
	case VAR_J:
	case VAR_J_OVER_ENV:
	case VAR_JPAR:
	{
		const int ncd = g.Ncelldata();
		if (ncd < 8) {result=0; break;}
		const double Bx1 = g.CT(5);
		const double By1 = g.CT(6);
		const double Bz1 = g.CT(7);
		double Bx0,By0,Bz0;
		if (ncd >= 11) {
			Bx0 = g.CT(8);
			By0 = g.CT(9);
			Bz0 = g.CT(10);
		} else {
			Bx0 = By0 = Bz0 = 0;
		}
		const double Bx = Bx0 + Bx1;
		const double By = By0 + By1;
		const double Bz = Bz0 + Bz1;
		const double invrho = 1.0/g.CT(0);
		const double vx = g.CT(1)*invrho;
		const double vy = g.CT(2)*invrho;
		const double vz = g.CT(3)*invrho;
		switch (var) {
		case VAR_BX: result=Bx; break;
		case VAR_BY: result=By; break;
		case VAR_BZ: result=Bz; break;
		case VAR_BX0: result=Bx0; break;
		case VAR_BY0: result=By0; break;
		case VAR_BZ0: result=Bz0; break;
		case VAR_BX1: result=Bx1; break;
		case VAR_BY1: result=By1; break;
		case VAR_BZ1: result=Bz1; break;
		case VAR_EX: result=By*vz-Bz*vy; break;
		case VAR_EY: result=Bz*vx-Bx*vz; break;
		case VAR_EZ: result=Bx*vy-By*vx; break;
		case VAR_E: result=sqrt(sqr(By*vz-Bz*vy) + sqr(Bz*vx-Bx*vz) + sqr(Bx*vy-By*vx)); break;
		case VAR_SX: result=invmu0*( sqr(By)*vx - Bx*By*vy + Bz*(Bz*vx-Bx*vz)); break;
		case VAR_SY: result=invmu0*(-Bx*By*vx + sqr(Bx)*vy + Bz*(Bz*vy-By*vz)); break;
		case VAR_SZ: result=invmu0*(-Bx*Bz*vx + sqr(Bx)*vz + By*(By*vz-Bz*vy)); break;
		case VAR_S: result=invmu0*sqrt(
			(sqr(Bx) + sqr(By) + sqr(Bz))*(
				sqr(Bz)*(sqr(vx)+sqr(vy)) - 2*Bx*Bz*vx*vz -
				2*By*vy*(Bx*vx+Bz*vz) + sqr(By)*(sqr(vx)+sqr(vz)) +
				sqr(Bx)*(sqr(vy)+sqr(vz)))); break;
		case VAR_KX:
		case VAR_KY:
		case VAR_KZ:
		case VAR_K:
		{
			// Energy flux vector K = (u+P-B^2/(2mu0))v + S
			const double Sx = invmu0*( sqr(By)*vx - Bx*By*vy + Bz*(Bz*vx-Bx*vz));
			const double Sy = invmu0*(-Bx*By*vx + sqr(Bx)*vy + Bz*(Bz*vy-By*vz));
			const double Sz = invmu0*(-Bx*Bz*vx + sqr(Bx)*vz + By*(By*vz-Bz*vy));
			real P;
			if (ncd == 5) {
				P = (gamma-1)*(g.CT(4) - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3))*0.5*invrho));
			} else {
				real U = g.CT(4);
				if (pseudobackgroundflag)
					U -= 0.5*invmu0*(sqr(g.CT(8)) + sqr(g.CT(9)) + sqr(g.CT(10)))
						+ (g.CT(5)*g.CT(8) + g.CT(6)*g.CT(9) + g.CT(7)*g.CT(10))*invmu0;
				P = (gamma-1)*(U - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))*0.5*invrho
							   - 0.5*invmu0*(sqr(g.CT(5)) + sqr(g.CT(6)) + sqr(g.CT(7))));
			}
			double Utot=0;
			if (ncd >= 11) {
				Utot = g.CT(4) + 0.5*invmu0*(sqr(g.CT(8)) + sqr(g.CT(9)) + sqr(g.CT(10)))
					+ (g.CT(5)*g.CT(8) + g.CT(6)*g.CT(9) + g.CT(7)*g.CT(10))*invmu0;
			} else if (ncd >= 8) {
				Utot = g.CT(4);
			}
			const double scalar = Utot + P - 0.5*invmu0*(sqr(Bx) + sqr(By) + sqr(Bz));
			const double Kx = Sx + vx*scalar;
			const double Ky = Sy + vy*scalar;
			const double Kz = Sz + vz*scalar;
			if (var == VAR_KX) {
				result = Kx;
			} else if (var == VAR_KY) {
				result = Ky;
			} else if (var == VAR_KZ) {
				result = Kz;
			} else if (var == VAR_K) {
				result = sqrt(sqr(Kx) + sqr(Ky) + sqr(Kz));
			}
			break;
		}
		case VAR_B: result=sqrt(sqr(Bx) + sqr(By) + sqr(Bz)); break;
		case VAR_B0: result=sqrt(sqr(Bx0) + sqr(By0) + sqr(Bz0)); break;
		case VAR_BMINUSB0: result=sqrt(sqr(Bx) + sqr(By) + sqr(Bz))-sqrt(sqr(Bx0) + sqr(By0) + sqr(Bz0)); break;
		case VAR_B1: result=sqrt(sqr(Bx1) + sqr(By1) + sqr(Bz1)); break;
		case VAR_B2: result=sqr(Bx) + sqr(By) + sqr(Bz); break;
		case VAR_B02: result=sqr(Bx0) + sqr(By0) + sqr(Bz0); break;
		case VAR_B12: result=sqr(Bx1) + sqr(By1) + sqr(Bz1); break;
		case VAR_LARMOR_P:
		case VAR_LARMOR_P_OVER_DX:
		case VAR_LARMOR_E:
		case VAR_GYRO_P:
		case VAR_GYRO_E:
		case VAR_DEBYELEN:
		{
			const real rho = g.CT(0);
			real P;
			if (ncd == 5) {
				P = (gamma-1)*(g.CT(4) - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/(2*rho));
			} else {
				real U = g.CT(4);
				if (pseudobackgroundflag)
					U -= 0.5*invmu0*(sqr(g.CT(8)) + sqr(g.CT(9)) + sqr(g.CT(10)))
						+ (g.CT(5)*g.CT(8) + g.CT(6)*g.CT(9) + g.CT(7)*g.CT(10))*invmu0;
				P = (gamma-1)*(U - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))/(2*rho)
							   - 0.5*invmu0*(sqr(g.CT(5)) + sqr(g.CT(6)) + sqr(g.CT(7))));
			}
			const real B = sqrt(sqr(Bx) + sqr(By) + sqr(Bz));
		        const real rL = sqrt(3*sqr(cnst::mp)*max(0.0,P/rho))/(cnst::qe*B);
			switch (var) {
			case VAR_LARMOR_P:
				result = rL;
				break;
			case VAR_LARMOR_E:
			        result = 0.5*rL*(cnst::me/cnst::mp);
				break;
			case VAR_GYRO_P:
			        result = (cnst::qe*B/cnst::mp)/(2*M_PI);
				break;
			case VAR_GYRO_E:
			        result = (cnst::qe*B/cnst::me)/(2*M_PI);
				break;
			case VAR_LARMOR_P_OVER_DX:
			{
				const TGridIndex i = g.find(X);
				result = (i != NOINDEX) ? rL*GridDimensionScaling/g.cellsize(i) : 0.0;
			}
			break;
			case VAR_DEBYELEN:
			{
			        const real n = rho/cnst::mp;
				const real Te = 0.25*(P/n);
			        result = sqrt(Te*cnst::epsilon0/n)/cnst::qe;
			}
			break;
			}
		}
		break;
		case VAR_INERTLEN_P:
		case VAR_INERTLEN_E:
		case VAR_PLASMAFREQ:
		{
			const real rho = g.CT(0);
		        const real ne = rho/cnst::mp;
		        const real wpe = sqrt(ne/(cnst::epsilon0*cnst::me))*cnst::qe;
		        const real wpi = sqrt(cnst::me/cnst::mp)*wpe;
			switch (var) {
			case VAR_INERTLEN_P:
			        result = cnst::c/wpi;
				break;
			case VAR_INERTLEN_E:
			        result = cnst::c/wpe;
				break;
			case VAR_PLASMAFREQ:
				result = wpe/(2*M_PI);
				break;
			}
		}
		break;
		case VAR_LOWERHYBRID:
		{
			const real rho = g.CT(0);
		        const real ne = rho/cnst::mp;
		        const real omega_pe = sqrt(ne/(cnst::epsilon0*cnst::me))*cnst::qe;
		        const real omega_pi = sqrt(cnst::me/cnst::mp)*omega_pe;
			const real B = sqrt(sqr(Bx) + sqr(By) + sqr(Bz));
		        const real Omega_e = cnst::qe*B/cnst::me;
		        const real Omega_i = cnst::qe*B/cnst::mp;
			const real omega_LH2 = (sqr(Omega_e)*sqr(Omega_i) + sqr(omega_pi)*sqr(Omega_e))/(sqr(omega_pe) + sqr(Omega_e));
			result = sqrt(omega_LH2)/(2*M_PI);
		}
		break;
		case VAR_DIVB:
		case VAR_DIVB_REL:
		{
			const real divB = ComputeDiv(g,X,5);
			if (var == VAR_DIVB)
				result = divB;
			else {
				real denom = ComputeDiv(g,X,5,true);
				if (denom <= 0) denom = 1;
				result = divB/denom;
			}
		}
		break;
		case VAR_JX:
		case VAR_JY:
		case VAR_JZ:
		case VAR_J:
		case VAR_JPAR:
		case VAR_J_OVER_ENV:
		{
			real jx,jy,jz;
			ComputeCurl(g,X,5,jx,jy,jz);
			jx*= invmu0;
			jy*= invmu0;
			jz*= invmu0;
			if (var == VAR_JX) {
				result = jx;
			} else if (var == VAR_JY) {
				result = jy;
			} else if (var == VAR_JZ) {
				result = jz;
			} else if (var == VAR_J) {
				result = sqrt(sqr(jx) + sqr(jy) + sqr(jz));
			} else if (var == VAR_J_OVER_ENV) {
				const real j = sqrt(sqr(jx) + sqr(jy) + sqr(jz));
				const real rhov = sqrt(sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)));
			        result = j/((cnst::qe/cnst::mp)*rhov);
			} else if (var == VAR_JPAR) {
				const real invB = 1.0/sqrt(sqr(Bx) + sqr(By) + sqr(Bz));
				const real bx = Bx*invB;
				const real by = By*invB;
				const real bz = Bz*invB;
				result = bx*jx + by*jy + bz*jz;
			}
		}
		break;
		default:
			result = 0;
			break;
		}
	}
	}
	if (logflag) {
		if (result < 1e-30) result = 1e-30;
		result = log10(result);
	}
	if (!(-1e30 <= result && result <= 1e30)) result = 0;
	return result;
}

vector<double> Tvariable::getSpectra(Tmetagrid& g, const Tdimvec& X) const
{
   vector<double> result;
   const int ncd = g.Ncelldata();
   for(int i=0;i<ncd;++i) {
      result.push_back(g.CT(i));
   }
   return result;
}
