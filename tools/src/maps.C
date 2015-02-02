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
#  pragma implementation "maps.H"
#endif

#include "maps.H"
#include <math.h>

void RampMap2D(should_be_const real st[MAXDIM][VECLEN], real xy[MAXDIM][VECLEN], smallnat vlen)
{
	real s,t,x,y;
	const real H = 1.0;
	const real W = 1.0;
	const real alpha = 30.0*(M_PI/180.0);
	const real L1 = 0.4*W;
	const real L2 = W - L1;
	const real slimit = L1/(L1+L2);
	const real tanalpha = tan(alpha);
	smallnat v;
	VDIRS;
	for (v=0; v<vlen; v++) {
		s = st[0][v];
		t = st[1][v];
		x = (L1 + L2)*s;
		if (s <= slimit)
			y = H*t;
		else {
			const real h = (x - L1)*tanalpha;
			y = h + (H-h)*t;
		}
		xy[0][v] = x;
		xy[1][v] = y;
	}
}

void RampMap3D(should_be_const real stw[MAXDIM][VECLEN], real xyz[MAXDIM][VECLEN], smallnat vlen)
{
	real s,t,w, x,y,z;
	const real H = 1.0;
	const real W = 1.0;
	const real alpha = 30.0*(M_PI/180.0);
	const real L1 = 0.4*W;
	const real L2 = W - L1;
	const real slimit = L1/(L1+L2);
	const real tanalpha = tan(alpha);
	smallnat v;
	VDIRS;
	for (v=0; v<vlen; v++) {
		s = stw[0][v];
		w = stw[1][v];
		t = stw[2][v];
		x = (L1 + L2)*s;
		y = w;
		if (s <= slimit)
			z = H*t;
		else {
			const real h = (x - L1)*tanalpha;
			z = h + (H-h)*t;
		}
		xyz[0][v] = x;
		xyz[1][v] = y;
		xyz[2][v] = z;
	}
}

real thetacoeff, phicoeff, theta1, phi1;

void ArcMap2D(should_be_const real stw[MAXDIM][VECLEN], real xyz[MAXDIM][VECLEN], smallnat vlen)
{
	real theta0,phi0,s, x,z;
	real Phi,sintheta0,sintheta02,sintheta04,sintheta08,Phi2,F,AUX,denom,costheta,sintheta,r;
	const real R_E = 6378e3;
	smallnat v;
	VDIRS;
	for (v=0; v<vlen; v++) {

		theta0 = theta1 + stw[0][v]/thetacoeff;
		phi0 = 0*M_PI/2;
		s = 1 + stw[1][v];
		Phi = cos(theta0)/(s*s);
		sintheta0 = sin(theta0);
		sintheta02 = sintheta0*sintheta0;
		sintheta04 = sintheta02*sintheta02;
		sintheta08 = sintheta04*sintheta04;
		Phi2 = Phi*Phi;
		F = 9*sintheta04*sqrt(9*sintheta08 + (256.0/3.0)*Phi2)
			+ 27*sintheta08 + 128*Phi2;
		AUX = 32*Phi2 +
			4*pow(2.0,2.0/3.0)*pow(F,1.0/3.0)*pow(Phi,4.0/3.0) +
			pow(2.0,1.0/3.0)*pow(F*Phi,2.0/3.0);
		denom = 12*pow(F,1.0/6.0)*pow(Phi,2.0/3.0)*pow(AUX,1.0/4.0);
		costheta =
			pow(2.0,2.0/3.0)*sqrt(3.0)*
			( pow(AUX,3.0/4.0)
			  - sqrt(12*Phi*(sqrt(3*F)*sintheta04
							 + pow(2.0,2.0/3.0)*sqrt(AUX)*pow(F*Phi,1.0/3.0))
					 - pow(AUX,3.0/2.0) )
				) /denom;
		sintheta = sqrt(1 - costheta*costheta);
		r = sintheta/sintheta0;
		r = r*r;
		x = r*sintheta*cos(phi0);
		z = r*sqrt(1 - sintheta*sintheta);
		xyz[0][v] = R_E*x;
		xyz[1][v] = R_E*z;
	}
}

void ArcMap3D(should_be_const real stw[MAXDIM][VECLEN], real xyz[MAXDIM][VECLEN], smallnat vlen)
{
	real theta0,phi0,s, x,y,z;
	real phi,Phi,sintheta0,sintheta02,sintheta04,sintheta08,Phi2,F,AUX,denom,costheta,sintheta,r;
	const real R_E = 6378e3;
	smallnat v;
	VDIRS;
	for (v=0; v<vlen; v++) {

		theta0 = theta1 + stw[0][v]/thetacoeff;
		phi0 = phi1 + stw[1][v]/phicoeff;
		s = stw[2][v];
		phi = phi0;
		Phi = cos(theta0)/(s*s);
		sintheta0 = sin(theta0);
		sintheta02 = sintheta0*sintheta0;
		sintheta04 = sintheta02*sintheta02;
		sintheta08 = sintheta04*sintheta04;
		Phi2 = Phi*Phi;
		F = 9*sintheta04*sqrt(9*sintheta08 + (256.0/3.0)*Phi2)
			+ 27*sintheta08 + 128*Phi2;
		AUX = 32*Phi2 +
			4*pow(2,2.0/3.0)*pow(F,1.0/3.0)*pow(Phi,4.0/3.0) +
			pow(2,1.0/3.0)*pow(F*Phi,2.0/3.0);
		denom = 12*pow(F,1.0/6.0)*pow(Phi,2.0/3.0)*pow(AUX,1.0/4.0);
		costheta =
			pow(2,2.0/3.0)*sqrt(3.0)*
			( pow(AUX,3.0/4.0)
			  - sqrt(12*Phi*(sqrt(3*F)*sintheta04
							 + pow(2,2.0/3.0)*sqrt(AUX)*pow(F*Phi,1.0/3.0))
					 - pow(AUX,3.0/2.0) )
				) /denom;
		sintheta = sqrt(1 - costheta*costheta);
		r = sintheta/sintheta0;
		r = r*r;
		x = r*sintheta*cos(phi);
		y = r*sintheta*sin(phi);
		z = r*sqrt(1 - sintheta*sintheta);
		xyz[0][v] = R_E*x;
		xyz[1][v] = R_E*y;
		xyz[2][v] = R_E*z;
	}
}

