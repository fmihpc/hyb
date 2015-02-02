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
#  pragma implementation "GLaxis.H"
#endif

#include "GLaxis.H"

#include <cmath>
#include <iostream>
#include <cstring>
#include <cstdio>
using namespace std;

static int nlabels;

static void DrawAxis1(float x1, float y1, float z1,
					  float x2, float y2, float z2,
					  float dx, float dy, float dz,
					  double u1, double u2,
					  int level,
					  bool labels, GLuint font)
{
	if (u2 < u1) {double tmp=u1; u1=u2; u2=tmp;}
	if (u1 == u2 || (u2-u1 < 1e-10*fabs(u2))) u2 = u1 + (fabs(u1) > 1e-20 ? fabs(u1) : 1e-20);
	const double invlog10 = 1.0/log(10.0);
	glBegin(GL_LINE_STRIP);
	glVertex3f(x1,y1,z1);
	glVertex3f(x2,y2,z2);
	glEnd();
	const double du = pow(10.0,floor(log(u2-u1)*invlog10-0.9999999*(level-1)));
	double u;
	const double C = 1/(u2-u1);
	const double umin = u1 - fmod(u1,du);
	glBegin(GL_LINES);
	for (u=umin; u<=u2; u+= du) if (u1 <= u && u <= u2) {
		const float t = (u-u1)*C;
		const float x = x1 + t*(x2-x1);
		const float y = y1 + t*(y2-y1);
		const float z = z1 + t*(z2-z1);
		glVertex3f(x,y,z);
		glVertex3f(x+dx,y+dy,z+dz);
	}
	glEnd();
	if (labels) {
		const double plusratio = 1.2;
		const double minusratio = 2.2;
		for (u=umin; u<=u2; u+= du) if (u1 <= u && u <= u2) {
			const float t = (u-u1)*C;
			float ratio;
			if (fabs(dx) > fabs(dy)) {
				ratio = (dx > 0) ? plusratio : minusratio;
			} else {
				ratio = (dy > 0) ? plusratio : minusratio;
			}
			const float x = x1 + t*(x2-x1) + ratio*dx;
			const float y = y1 + t*(y2-y1) + ratio*dy;
			const float z = z1 + t*(z2-z1) + ratio*dz;
			static char lab[80];
			if (fabs(u) < 1e-7*(fabs(u1) > fabs(u2) ? fabs(u1) : fabs(u2))) {
				lab[0] = '0';
				lab[1] = '\0';
			} else
				sprintf(lab,"%1.5g",u);
			glRasterPos3f(x,y,z);
			glListBase(font);
			glCallLists(strlen(lab),GL_BYTE,lab);
			nlabels++;
		}
	}
}

void DrawAxis(float x1, float y1, float z1,
			  float x2, float y2, float z2,
			  float dx, float dy, float dz,
			  double u1, double u2,
			  bool include_label, GLuint font)
{
	const double ratio = 0.5;
	nlabels = 0;
	DrawAxis1(x1,y1,z1, x2,y2,z2, dx,dy,dz, u1,u2, 1,include_label,font);
	if (include_label && nlabels < 2) {
		GLfloat rgba[4] = {0,0,0,0};
		glGetFloatv(GL_COLOR_CLEAR_VALUE,&rgba[0]);
		// NOTICE: Mesa3.0 is buggy here: it only fills rgba[0].
		// Assume here that rgba[1]=rgba[2]=rgba[0], This only works if the background is either black or white.
		glPushAttrib(GL_CURRENT_BIT);
		glColor3f(rgba[0],rgba[0/*1*/],rgba[0/*2*/]);
		DrawAxis1(x1,y1,z1, x2,y2,z2, dx,dy,dz, u1,u2, 1,true,font);	// overdraw the previous axis with background color
		glPopAttrib();
		DrawAxis1(x1,y1,z1, x2,y2,z2, dx,dy,dz, u1,u2, 1,false,font);	// redraw it with fg color but without label
		DrawAxis1(x1,y1,z1, x2,y2,z2, ratio*dx,ratio*dy,ratio*dz, u1,u2, 2,true,font);	// draw the finer axis with label
	} else
		DrawAxis1(x1,y1,z1, x2,y2,z2, ratio*dx,ratio*dy,ratio*dz, u1,u2, 2,false,font);
}
