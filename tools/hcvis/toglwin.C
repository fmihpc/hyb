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
#  pragma implementation "toglwin.H"
#endif

#include <iostream>
#include <cstdio>
#include <cstdlib>
using namespace std;
#include <GL/glu.h>
#include "toglwin.H"
#include "gridcache.H"
#include "intpolcache.H"
#include "GLaxis.H"

TToglWindow::TToglWindowEntry *TToglWindow::ToglWindowList = 0;
const real antialiased_FL_width = 1.5;

#define Random() (rand()/(RAND_MAX+1.0))		/* needed in RenderVolumetric */

TToglWindow::TToglWindow()
{
	filename=strdup("NONAME");
	variable="rho";
	gridptr = 0;
	LinearInterpolation = false;
	DrawGrid = false;
	DrawGhosts = false;
	HideGhosts = false;
	DrawDead = false;
	DrawCont = false;
	DrawIsosurf = false;
	DrawVolumetric = false;
	VolumetricAntialiasing = false;
	volumetric_npoints = 400000;
	volumetric_alpha = 0.2;
	antialiased_lines = false;
	smooth_isosurfaces = false;
	shiny_isosurfaces = false;
	view3D = false;
	Logarithmic = false;
	verbose = false;
	PreserveAspect = true;
	needs_redraw = true;
	XYZmode_in_3D_params = false;
	display_list_exists = false;
	can_use_display_list = false;
	Mapping = false;
	ColorBar = true;
	Title = true;
	BoxDefined = false;
	WhiteBackground = false;

	relative_radius = 1.0;
	theta = 50;
	phi = 50;
	cx = 0;
	cy = 0;
	cz = 0;
	field_of_view = 30;
	IsoValue[0] = 1.0;
	Nisovalues = 1;

	slicedim = 1;		// Y=const slices are displayed by default for 3D datasets
	slice_xval = 0;		// some default, could be uninitialized also, but better use something here
	objlist = 0;
	oldheight = oldwidth = -1;
	AlphaValue = 1.0;		// transparency: 1 is opaque, 0 is fully transparent (invisible)
	
	datamin = 0;
	datamax = 2;
	one_over_datamax_minus_datamin = 1.0/(datamax - datamin);
	forced_dataminmax = false;

	wmin = wmax = hmin = hmax = 0;		// when a file is loaded, these are initialized from MappedBox, if they are still zero
	
	dirty = true;

}

TToglWindow::~TToglWindow()
{
	if (dirty) return;
	Togl_UnloadBitmapFont(togl,helvetica_10);
	Togl_UnloadBitmapFont(togl,helvetica_12);
	Togl_UnloadBitmapFont(togl,helvetica_18);
	ClearObjects();
}

void TToglWindow::ClearObjects()
{
	T3DObject *p;
	while (objlist) {
		p = objlist;
		objlist = p->next;
		delete p;
	}
	needs_redraw = true;
}

void TToglWindow::AddSlice(int SlicedDim, double Slice_Xval, bool DontShowIn3D)
{
	T3DObject *p;
	p = new T3DSliceObject(SlicedDim,Slice_Xval);
	((T3DSliceObject *)p)->DontShowIn3D = DontShowIn3D;
	p->next = objlist;
	objlist = p;
	needs_redraw = true;
}

void TToglWindow::AddSphere(double r, double x0, double y0, double z0, double DeltaTheta, double DeltaThetaGrid)
{
	T3DObject *p;
	p = new T3DSphereObject(r,x0,y0,z0,DeltaTheta,DeltaThetaGrid);
	p->next = objlist;
	objlist = p;
	needs_redraw = true;
}

void TToglWindow::AddFieldLineBunch
   (int N, const double r1[3], const double r2[3], TVectorField vecfield, TTraceDirection dir,TFieldLineDistribution distr,
	TLoopThresholdType thresholdtype, double threshold)
{
	T3DObject *p;
	p = new T3DFieldLineBunchObject(N,r1,r2,vecfield,dir,distr,thresholdtype,threshold);
	p->next = objlist;
	objlist = p;
	needs_redraw = true;
}

// Count number of d-directed slices in the object list
int TToglWindow::Nslices(int d) const
{
	T3DObject *p;
	int cnt = 0;
	for (p=objlist; p; p=p->next) if (!strcmp(p->type(),"Slice") && ((T3DSliceObject *)p)->SlicedDim() == d) cnt++;
	return cnt;
}

T3DSliceObject *TToglWindow::NthSlice(int d, int n) const
{
	T3DObject *p;
	int cnt;
	for (p=objlist,cnt=0; cnt<n; cnt++,p=p->next) {
		while (strcmp(p->type(),"Slice") || ((T3DSliceObject *)p)->SlicedDim() != d) p=p->next;
	}
	while (strcmp(p->type(),"Slice") || ((T3DSliceObject *)p)->SlicedDim() != d) p=p->next;
	if (p)
		return (T3DSliceObject *)p;
	else {
		cerr << "*** TToglWindow::NthSlice: no such slice (d=" << d << ", n=" << n << ")\n";
		return 0;
	}
}

void TToglWindow::ClearMinMax(const char *reason)
{
	if (forced_dataminmax) return;
	var.clearAllMinMax();
	if (verbose) cout << "*ClearMinMax because of " << reason << "\n";
	real Xmin[3]={0,0,0},Xmax[3]={0,0,0};
	int XX, YY;
	switch (slicedim) {
	case 0: XX=1; YY=2; break;
	case 1: XX=0; YY=2; break;
	default: case 2: XX=0; YY=1; break;
	}
	if (gridptr->dimension()==2) {XX=0; YY=1;}
	MappedBox(*gridptr,Xmin,Xmax);
	wmin_global = Xmin[XX]; wmax_global = Xmax[XX];
	hmin_global = Xmin[YY]; hmax_global = Xmax[YY];
	if (wmin==0 && wmax==0 && hmin==0 && hmax==0) {
		wmin=wmin_global; wmax=wmax_global;
		hmin=hmin_global; hmax=hmax_global;
	}
}


void TToglWindow::SetFilename(const char *fn)
{
	if (!strcmp(fn,filename)) return;	// do nothing if the name does not change (avoid closing/opening the grid file)
	if (filename) free(filename);
	filename = strdup(fn);
	TGridCache::verbose = verbose;
	static bool FirstTime = true;
	if (FirstTime) {
		theGridCache.init();
		FirstTime = false;
	}
	if (gridptr) theGridCache.close(gridptr);
	gridptr = theGridCache.open(filename,gamma,invmu0,mass,pseudobackground);
	if (!gridptr) return;
	ClearMinMax("SetFilename");
}

void TToglWindow::New(const char *ident)
{
	TToglWindowEntry *p = new TToglWindowEntry;
	p->winptr = new TToglWindow;
	p->ident = strdup(ident);
	p->next = ToglWindowList;
	ToglWindowList = p;
}

TToglWindow* TToglWindow::Find(const char *ident)
{
	TToglWindowEntry *p;
	for (p=ToglWindowList; p; p=p->next)
		if (!strcmp(p->ident, ident)) return p->winptr;
	return 0;
}

void TToglWindow::Delete(const char *ident)
{
	TToglWindowEntry *p;
	for (p=ToglWindowList; p; p=p->next) {
		if (!ident && !p->ident) break;
		if (p->ident && ident && !strcmp(p->ident, ident)) break;
	}
	if (!p) {cerr << "*** DeleteToglWindow(" << (ident ? ident : "NULL") << ") failed\n"; return;}
	if (p == ToglWindowList)
		ToglWindowList = p->next;
	else {
		TToglWindowEntry *prev;
		for (prev=ToglWindowList; prev->next!=p; prev=prev->next);
		prev->next = p->next;
	}
	delete p->winptr;
	delete p;
}

const bool osflag = false;

static real map_Xin[MAXDIM][VECLEN], map_Xout[MAXDIM][VECLEN];
static int map_XX = 0, map_YY = 1;

void TToglWindow::vertex2f(GLfloat x, GLfloat y) const
{
	if (Mapping && MapFunctionPtr) {
		map_Xin[map_XX][0] = x; map_Xin[map_YY][0] = y;
		(*MapFunctionPtr)(map_Xin,map_Xout,1);
		glVertex2f(map_Xout[map_XX][0], map_Xout[map_YY][0]);
	} else {
		glVertex2f(x,y);
	}
}

static int vert3_select = -1;		// 0: take yz, 1: take xz, 2: take xy, -1: take xyz

void TToglWindow::vertex3f(GLfloat x, GLfloat y, GLfloat z) const
{
	GLfloat xx,yy,zz;
	if (Mapping && MapFunctionPtr) {
		map_Xin[0][0] = x; map_Xin[1][0] = y; map_Xin[2][0] = z;
		(*MapFunctionPtr)(map_Xin,map_Xout,1);
		xx = map_Xout[0][0];
		yy = map_Xout[1][0];
		zz = map_Xout[2][0];
	} else {
		xx = x;
		yy = y;
		zz = z;
	}
	switch (vert3_select) {
	case 0: glVertex2f(yy,zz); break;
	case 1: glVertex2f(xx,zz); break;
	case 2: glVertex2f(xx,yy); break;
	default: glVertex3f(xx,yy,zz);
	}
}

inline void TToglWindow::SetColorFromScaledQuantity(GLfloat uscaled) const
{
	GLfloat r,g,b;
	palette.RGB(uscaled, r,g,b);
	glColor4f(r,g,b,AlphaValue);
}

inline void TToglWindow::SetColor(GLfloat u) const
{
	double uscaled = (u-datamin)*one_over_datamax_minus_datamin;
	if (uscaled < 0) uscaled = 0;
	if (uscaled > 1) uscaled = 1;
	SetColorFromScaledQuantity(uscaled);
}

//  u[2]  ---   u[3]
//   |           |
//   |           |
//  u[0]  ---   u[1]
// udens, ufacecenters numbering:
//        3
//        |
//   0  -----  1
//        |
//        2
void TToglWindow::DisplaySquare(const Tdimvec& x, GLfloat dx,
								const Tdimvec& ex, const Tdimvec& ey,
								real ucenter, const real u[4], const bool udens[4], const real ufacecenters[4])
{
	const GLfloat halfdx = 0.500001*dx;
	Tdimvec hx,hy;
	hx[0] = halfdx*ex(0); hy[0] = halfdx*ey(0);
	hx[1] = halfdx*ex(1); hy[1] = halfdx*ey(1);
	hx[2] = halfdx*ex(2); hy[2] = halfdx*ey(2);
	glBegin(GL_TRIANGLE_FAN);
	SetColor(ucenter);
	vertex3f(x(0),x(1),x(2));
	SetColor(u[0]);
	vertex3f(x(0)-hx(0)-hy(0), x(1)-hx(1)-hy(1), x(2)-hx(2)-hy(2));
	SetColor(ufacecenters[2]);
	vertex3f(x(0)-hy(0),       x(1)-hy(1),       x(2)-hy(2));
	SetColor(u[1]);
	vertex3f(x(0)+hx(0)-hy(0), x(1)+hx(1)-hy(1), x(2)+hx(2)-hy(2));
	SetColor(ufacecenters[1]);
	vertex3f(x(0)+hx(0),       x(1)+hx(1),       x(2)+hx(2));
	SetColor(u[3]);
	vertex3f(x(0)+hx(0)+hy(0), x(1)+hx(1)+hy(1), x(2)+hx(2)+hy(2));
	SetColor(ufacecenters[3]);
	vertex3f(x(0)+hy(0),       x(1)+hy(1),       x(2)+hy(2));
	SetColor(u[2]);
	vertex3f(x(0)-hx(0)+hy(0), x(1)-hx(1)+hy(1), x(2)-hx(2)+hy(2));
	SetColor(ufacecenters[0]);
	vertex3f(x(0)-hx(0),       x(1)-hx(1),       x(2)-hx(2));
	SetColor(u[0]);
	vertex3f(x(0)-0.5*dx*ex(0)-hy(0), x(1)-0.5*dx*ex(1)-hy(1), x(2)-0.5*dx*ex(2)-hy(2));
	glEnd();
	if (DrawCont) {
		glColor3f(0,0,0);
		const real ex1[3] = {ex(0),ex(1),ex(2)};
		const real ey1[3] = {ey(0),ey(1),ey(2)};
		GLContourSquare(x(0),x(1),x(2),dx,u,udens,ufacecenters,cs,ex1,ey1,vert3_select);
	}
}

inline void TToglWindow::DisplaySquareFlat(const Tdimvec& x, GLfloat dx, const Tdimvec& ex, const Tdimvec& ey, real ucenter)
{
	const GLfloat halfdx = 0.500001*dx;
	Tdimvec hx,hy;
	hx[0] = halfdx*ex(0); hy[0] = halfdx*ey(0);
	hx[1] = halfdx*ex(1); hy[1] = halfdx*ey(1);
	hx[2] = halfdx*ex(2); hy[2] = halfdx*ey(2);
	glBegin(GL_QUADS);
	SetColor(ucenter);
	vertex3f(x(0)-hx(0)-hy(0), x(1)-hx(1)-hy(1), x(2)-hx(2)-hy(2));
	vertex3f(x(0)+hx(0)-hy(0), x(1)+hx(1)-hy(1), x(2)+hx(2)-hy(2));
	vertex3f(x(0)+hx(0)+hy(0), x(1)+hx(1)+hy(1), x(2)+hx(2)+hy(2));
	vertex3f(x(0)-hx(0)+hy(0), x(1)-hx(1)+hy(1), x(2)-hx(2)+hy(2));
	glEnd();
}

inline void TToglWindow::DisplaySquareGray(const Tdimvec& x, GLfloat dx, const Tdimvec& ex, const Tdimvec& ey)
{
	const GLfloat halfdx = 0.500001*dx;
	Tdimvec hx,hy;
	hx[0] = halfdx*ex(0); hy[0] = halfdx*ey(0);
	hx[1] = halfdx*ex(1); hy[1] = halfdx*ey(1);
	hx[2] = halfdx*ex(2); hy[2] = halfdx*ey(2);
	glBegin(GL_QUADS);
	glColor4f(0.5,0.5,0.5,AlphaValue);
	vertex3f(x(0)-hx(0)-hy(0), x(1)-hx(1)-hy(1), x(2)-hx(2)-hy(2));
	vertex3f(x(0)+hx(0)-hy(0), x(1)+hx(1)-hy(1), x(2)+hx(2)-hy(2));
	vertex3f(x(0)+hx(0)+hy(0), x(1)+hx(1)+hy(1), x(2)+hx(2)+hy(2));
	vertex3f(x(0)-hx(0)+hy(0), x(1)-hx(1)+hy(1), x(2)-hx(2)+hy(2));
	glEnd();
}

void TToglWindow::MappedBox(const Tmetagrid& g, real xmin[3], real xmax[3]) const
{
	static real map_Xin[MAXDIM][VECLEN], map_Xout[MAXDIM][VECLEN];
	g.getbox(xmin,xmax);
	const smallnat dim = g.dimension();
	smallnat d;
	if (!Mapping || !MapFunctionPtr) {
		g.getbox(xmin,xmax);
		if (BoxDefined) {
			for (d=0; d<dim; d++) {
				if (xmin[d] < box_xmin[d]) xmin[d] = box_xmin[d];
				if (xmax[d] > box_xmax[d]) xmax[d] = box_xmax[d];
			}
		}
		return;
	}
	real smin[3],smax[3];
	g.getbox(smin,smax);
#if VECLEN < 8
#error VECLEN < 8
#endif
	int npts=2, i;
	if (dim == 2) {
		map_Xin[0][0] = smin[0];	// x-, y-
		map_Xin[1][0] = smin[1];
		map_Xin[0][1] = smax[0];	// x+, y-
		map_Xin[1][1] = smin[1];
		map_Xin[0][2] = smax[0];	// x+, y+
		map_Xin[1][2] = smax[1];
		map_Xin[0][3] = smin[0];	// x-, y+
		map_Xin[1][3] = smax[1];
		npts = 4;
	} else if (dim == 3) {
		map_Xin[0][0] = smin[0];	// x-, y-, z-
		map_Xin[1][0] = smin[1];
		map_Xin[2][0] = smin[2];
		map_Xin[0][1] = smax[0];	// x+, y-, z-
		map_Xin[1][1] = smin[1];
		map_Xin[2][1] = smin[2];
		map_Xin[0][2] = smax[0];	// x+, y+, z-
		map_Xin[1][2] = smax[1];
		map_Xin[2][2] = smin[2];
		map_Xin[0][3] = smin[0];	// x-, y+, z-
		map_Xin[1][3] = smax[1];
		map_Xin[2][3] = smin[2];
		map_Xin[0][4] = smin[0];	// x-, y-, z+
		map_Xin[1][4] = smin[1];
		map_Xin[2][4] = smax[2];
		map_Xin[0][5] = smax[0];	// x+, y-, z+
		map_Xin[1][5] = smin[1];
		map_Xin[2][5] = smax[2];
		map_Xin[0][6] = smax[0];	// x+, y+, z+
		map_Xin[1][6] = smax[1];
		map_Xin[2][6] = smax[2];
		map_Xin[0][7] = smin[0];	// x-, y+, z+
		map_Xin[1][7] = smax[1];
		map_Xin[2][7] = smax[2];
		npts = 8;
	}
	(*MapFunctionPtr)(map_Xin,map_Xout,npts);
	for (d=0; d<dim; d++) xmin[d] = xmax[d] = map_Xout[d][0];
	for (d=0; d<dim; d++) {
		for (i=0; i<npts; i++) {
			xmin[d] = min(xmin[d], map_Xout[d][i]);
			xmax[d] = max(xmax[d], map_Xout[d][i]);
		}
	}
	if (BoxDefined) {
		for (d=0; d<dim; d++) {
			if (xmin[d] < box_xmin[d]) xmin[d] = box_xmin[d];
			if (xmax[d] > box_xmax[d]) xmax[d] = box_xmax[d];
		}
	}
}

bool TToglWindow::IntersectsPlane(Tmetagrid& g, TGridIndex c, int SlicedDim, real Slice_Xval, Tdimvec& X) const
{
	g.centroid(c,X);
	if (BoxDefined) {
		smallnat d;
		bool isoutside = false;
		for (d=0; d<g.dimension(); d++)
			if (X[d] < box_xmin[d] || X[d] > box_xmax[d]) isoutside = true;
		if (isoutside) return false;
	}
	return fabs(X[SlicedDim] - Slice_Xval) < 0.6*g.cellsize(c);
}

inline bool TToglWindow::InsideBox(const Tdimvec& x) const
{
	if (!BoxDefined) return true;
	return (box_xmin[0] <= x(0) && x(0) <= box_xmax[0] &&
			box_xmin[1] <= x(1) && x(1) <= box_xmax[1] &&
			box_xmin[2] <= x(2) && x(2) <= box_xmax[2]);
}

void TToglWindow::DisplayGridForOneCell
    (const Tdimvec& Xc, const GLfloat halfdx,
	 const Tdimvec& ex, const Tdimvec& ey) const
{
	glColor3f(0,0,0);
	glBegin(GL_LINE_STRIP);
	vertex3f(Xc(0)+halfdx*ex(0)-halfdx*ey(0),
			 Xc(1)+halfdx*ex(1)-halfdx*ey(1),
			 Xc(2)+halfdx*ex(2)-halfdx*ey(2));
	vertex3f(Xc(0)-halfdx*ex(0)-halfdx*ey(0),
			 Xc(1)-halfdx*ex(1)-halfdx*ey(1),
			 Xc(2)-halfdx*ex(2)-halfdx*ey(2));
	vertex3f(Xc(0)-halfdx*ex(0)+halfdx*ey(0),
			 Xc(1)-halfdx*ex(1)+halfdx*ey(1),
			 Xc(2)-halfdx*ex(2)+halfdx*ey(2));
	glEnd();
}

void TToglWindow::DisplayGrid(Tmetagrid& g, int SlicedDim, real Slice_Xval)
{
	TGridIndex c;
	const int dim = g.dimension();
	Tdimvec Xc;
	Tdimvec ex(0,0,0), ey(0,0,0);
	int XX,YY;
	switch (SlicedDim) {
	case 0: XX=1; YY=2; break;
	case 1: XX=0; YY=2; break;
	default: case 2: XX=0; YY=1; break;
	}
	ex[XX] = 1;
	ey[YY] = 1;
	for (c=g.first(); !g.isover(c); c=g.next(c)) {
		if (!g.isleaf(c)) continue;
		const TCellType ct = g.celltype(c);
		if (!DrawDead && ct == DEAD_CELL) continue;
		if (HideGhosts && ct == GHOST_CELL) continue;
		if (ct == REMOVED_CELL) continue;
		if (dim==3 && !IntersectsPlane(g,c,SlicedDim,Slice_Xval,Xc)) continue;
		if (dim==2) g.centroid(c,Xc);
		Xc[SlicedDim] = Slice_Xval;
		const GLfloat halfdx = 0.5*g.cellsize(c);
		DisplayGridForOneCell(Xc,halfdx,ex,ey);
	}
}

void TToglWindow::DataSliceMinMax(Tmetagrid& g, int SlicedDim, real Slice_Xval, bool ToOldExtrema)
{
	if (verbose) cout << "DataSliceMinMax(SlicedDim=" << SlicedDim << ", Slice_Xval=" << Slice_Xval
					  << ", " << ToOldExtrema << ")";
	TGridIndex c;
	Tdimvec Xc;
	const int dim = g.dimension();
	bool FirstTime = ToOldExtrema ? false : true;
	for (c=g.first(); !g.isover(c); c=g.next(c)) {
		if (!g.isleaf(c)) continue;
		if (dim==3 && !IntersectsPlane(g,c,SlicedDim,Slice_Xval,Xc)) continue;
		const TCellType ct = g.celltype(c);
		if (!DrawDead && ct == DEAD_CELL) continue;
		if (ct == REMOVED_CELL) continue;
		if (dim==2) g.centroid(c,Xc);
		g > c;
		Xc[SlicedDim] = Slice_Xval;
		const real dataval = var.get(g,Xc);
		if (dataval < datamin || FirstTime) datamin = dataval;
		if (dataval > datamax || FirstTime) datamax = dataval;
		FirstTime = false;
	}
	if (datamax <= datamin) datamax = datamin + 1;
	if (datamax == datamin) datamax = datamin*1.00000001;
	one_over_datamax_minus_datamin = 1.0/(datamax - datamin);
	if (verbose) cout << "--> datamin=" << datamin << ", datamax=" << datamax << "\n" << flush;
}

void TToglWindow::DisplayData2D(Tmetagrid& g, int SlicedDim, real Slice_Xval)
{
	if (verbose) cout << "DisplayData2D(SlicedDim=" << SlicedDim << ", Slice_Xval=" << Slice_Xval << ")\n" << flush;
	TGridIndex c;
	real ucorners[4];
	Tdimvec X,Xc;
	smallnat d;
	for (d=0; d<MAXDIM; d++) X[d] = 0;
	int XX,YY;
	Tdimvec ex(0,0,0), ey(0,0,0);
	switch (SlicedDim) {
	case 0: XX=1; YY=2; break;
	case 1: XX=0; YY=2; break;
	default: case 2: XX=0; YY=1; break;
	}
	X[SlicedDim] = Slice_Xval;
	ex[XX] = 1;
	ey[YY] = 1;
	real xmin[3]={0,0,0},xmax[3]={0,0,0};
	MappedBox(g,xmin,xmax);
	TIntpolCache cache(g.BasegridCellSize(), xmin[XX], xmin[YY]);
	const int dim = g.dimension();
	const bool draw_grid_at_the_same_time = (DrawGrid && g.dimension() == 3 && view3D);
	GLfloat AlphaValue_sav = 1.0;
	if (antialiased_lines) {
		glPushAttrib(GL_ENABLE_BIT);
		glEnable(GL_LINE_SMOOTH);			// If we draw contours, draw them smoothly
		glEnable(GL_BLEND);					// Must enable blending also in that case
		if (!view3D) {		// But since blending is now on, alpha value must be 1 for 2D drawing
			AlphaValue_sav = AlphaValue;
			AlphaValue = 1.0;
		}
	}
	for (c=g.first(); !g.isover(c); c=g.next(c)) {
		if (!g.isleaf(c)) continue;
		if (dim==3 && !IntersectsPlane(g,c,SlicedDim,Slice_Xval,Xc)) continue;
		if (dim==2) g.centroid(c,Xc);
		const TCellType ct = g.celltype(c);
		if (!DrawGhosts && ct == GHOST_CELL) {
			if (!HideGhosts) {
				Xc[SlicedDim] = Slice_Xval;
				DisplaySquareGray(Xc,g.cellsize(c),ex,ey);
			}
			continue;
		}
		if (!DrawDead && ct == DEAD_CELL) continue;
		if (ct == REMOVED_CELL) continue;
		const real halfdx = 0.5*g.cellsize(c);
		if (LinearInterpolation) {
			Xc[SlicedDim] = Slice_Xval;

			X[XX] = Xc(XX) - halfdx;
			X[YY] = Xc(YY) - halfdx;
			if (!cache.read(ucorners[0], X(XX),X(YY))) {
				g.intpol(X,1,true);
				ucorners[0] = var.get(g,X);
				cache.store(ucorners[0],X(XX),X(YY));
			}

			X[XX] = Xc(XX) + halfdx;
			X[YY] = Xc(YY) - halfdx;
			if (!cache.read(ucorners[1], X(XX),X(YY))) {
				g.intpol(X,1,true);
				ucorners[1] = var.get(g,X);
				cache.store(ucorners[1],X(XX),X(YY));
			}

			X[XX] = Xc(XX) + halfdx;
			X[YY] = Xc(YY) + halfdx;
			if (!cache.read(ucorners[3], X(XX),X(YY))) {
				g.intpol(X,1,true);
				ucorners[3] = var.get(g,X);
				cache.store(ucorners[3],X(XX),X(YY));
			}

			X[XX] = Xc(XX) - halfdx;
			X[YY] = Xc(YY) + halfdx;
			if (!cache.read(ucorners[2], X(XX),X(YY))) {
				g.intpol(X,1,true);
				ucorners[2] = var.get(g,X);
				cache.store(ucorners[2],X(XX),X(YY));
			}

			bool udens[4] = {false,false,false,false};
			real ufacecenters[4];

			const smallnat LEFT=0, RIGHT=1;
			udens[0] = g.Nneighbours(c,0,LEFT) > 1;
			udens[1] = g.Nneighbours(c,0,RIGHT) > 1;
			udens[2] = g.Nneighbours(c,1,LEFT) > 1;
			udens[3] = g.Nneighbours(c,1,RIGHT) > 1;

			const real uaverage = 0.25*(ucorners[0]+ucorners[1]+ucorners[2]+ucorners[3]);

			if (udens[0]) {
				X[XX] = Xc(XX) - halfdx;
				X[YY] = Xc(YY);
				if (!cache.read(ufacecenters[0], X(XX),X(YY))) {
					g.intpol(X,1,true);
					ufacecenters[0] = var.get(g,X);
					cache.store(ufacecenters[0],X(XX),X(YY));
				}
			} else
				ufacecenters[0] = 0.5*(ucorners[0] + ucorners[2]);
			if (udens[1]) {
				X[XX] = Xc(XX) + halfdx;
				X[YY] = Xc(YY);
				if (!cache.read(ufacecenters[1], X(XX),X(YY))) {
					g.intpol(X,1,true);
					ufacecenters[1] = var.get(g,X);
					cache.store(ufacecenters[1],X(XX),X(YY));
				}
			} else
				ufacecenters[1] = 0.5*(ucorners[1] + ucorners[3]);
			if (udens[2]) {
				X[XX] = Xc(XX);
				X[YY] = Xc(YY) - halfdx;
				if (!cache.read(ufacecenters[2], X(XX),X(YY))) {
					g.intpol(X,1,true);
					ufacecenters[2] = var.get(g,X);
					cache.store(ufacecenters[2],X(XX),X(YY));
				}
			} else
				ufacecenters[2] = 0.5*(ucorners[0] + ucorners[1]);
			if (udens[3]) {
				X[XX] = Xc(XX);
				X[YY] = Xc(YY) + halfdx;
				if (!cache.read(ufacecenters[3], X(XX),X(YY))) {
					g.intpol(X,1,true);
					ufacecenters[3] = var.get(g,X);
					cache.store(ufacecenters[3],X(XX),X(YY));
				}
			} else
				ufacecenters[3] = 0.5*(ucorners[2] + ucorners[3]);

			DisplaySquare(Xc,g.cellsize(c),ex,ey,
						  uaverage,
						  ucorners,udens,ufacecenters);
		} else {
			g > c;
			Xc[SlicedDim] = Slice_Xval;
			DisplaySquareFlat(Xc,2*halfdx,ex,ey,var.get(g,Xc));
		}
		if (draw_grid_at_the_same_time) {
			const real epsilon = g.BasegridCellSize()*0.03;
			Xc[SlicedDim] = Slice_Xval + epsilon;
			DisplayGridForOneCell(Xc,halfdx,ex,ey);
			Xc[SlicedDim] = Slice_Xval - epsilon;
			DisplayGridForOneCell(Xc,halfdx,ex,ey);
		}
	}
	if (antialiased_lines) {
		glPopAttrib();
		if (!view3D) {
			AlphaValue = AlphaValue_sav;
		}
	}
	if (DrawGrid && !draw_grid_at_the_same_time) {
		if (view3D) {
			const real epsilon = g.BasegridCellSize()*0.03;
			DisplayGrid(g,SlicedDim,Slice_Xval+epsilon);
			DisplayGrid(g,SlicedDim,Slice_Xval-epsilon);
		} else
			DisplayGrid(g,SlicedDim,Slice_Xval);
	}
	if (LinearInterpolation && verbose) cout << "intpol cache: " << cache.hitratio()*100 << " % hit ratio\n" << flush;
}

void TToglWindow::DisplayDataSphere(Tmetagrid& g, real r, real x0, real y0, real z0, real deltatheta,
									real deltatheta_grid,
									bool DrawIt, bool ToOldExtrema)
{
	if (g.dimension()!=3) return;
	const real pi = M_PI;
	Tdimvec X;
	int n,m;
	const real deltaphi = deltatheta;
	const real deltaphi_grid = deltatheta_grid;
	int N = 1 + int(ceil(180/deltatheta));
	if (N > 1000) N = 1000;
	if (N < 5) N = 5;
	const int M = (N-1)*2 + 1;
	real *u = new real [4*N*M];
	if (!u) {cerr << "*** DisplayDataSphere: out of memory\n"; return;}
	real *x = u + N*M;
	real *y = u + 2*N*M;
	real *z = u + 3*N*M;
	bool FirstTime = true;
	real umin=0,umax=0;
	for (n=0; n<N; n++) for (m=0; m<M; m++) {
		const real theta = deltatheta*n;
		const real phi = deltaphi*m;
		const real thetaRad = (pi/180)*theta;
		const real phiRad = (pi/180)*phi;
		X[0] = x0 + r*sin(thetaRad)*cos(phiRad);
		X[1] = y0 + r*sin(thetaRad)*sin(phiRad);
		X[2] = z0 + r*cos(thetaRad);
		g.intpol(X,1,true);
		const real uu = var.get(g,X);
		u[n*M+m] = uu;
		if (FirstTime) {
			umin = umax = uu;
			FirstTime = false;
		} else {
			if (uu < umin) umin = uu;
			if (uu > umax) umax = uu;
		}
		x[n*M+m] = X(0);
		y[n*M+m] = X(1);
		z[n*M+m] = X(2);
	}
	if (!DrawIt && !ToOldExtrema) {
		datamin = umin;
		datamax = umax;
	} else {
		datamin = min(datamin,umin);
		datamax = max(datamax,umax);
	}
	one_over_datamax_minus_datamin = 1.0/(datamax - datamin);
	if (DrawIt) {
		for (n=0; n<N-1; n++) for (m=0; m<M-1; m++) {
			glBegin(GL_QUADS);
			SetColor(u[n*M+m]);
			vertex3f(x[n*M+m],        y[n*M+m],        z[n*M+m]);
			SetColor(u[(n+1)*M+m]);
			vertex3f(x[(n+1)*M+m],    y[(n+1)*M+m],    z[(n+1)*M+m]);
			SetColor(u[(n+1)*M+(m+1)]);
			vertex3f(x[(n+1)*M+(m+1)],y[(n+1)*M+(m+1)],z[(n+1)*M+(m+1)]);
			SetColor(u[n*M+m+1]);
			vertex3f(x[n*M+m+1],      y[n*M+m+1],      z[n*M+m+1]);
			glEnd();
			if (DrawCont) {
				glColor3f(0,0,0);
				const real C = 1.003;
				GLContourTriangle_3D(C*x[n*M+m],      C*y[n*M+m],      C*z[n*M+m],
									 C*x[(n+1)*M+m],  C*y[(n+1)*M+m],  C*z[(n+1)*M+m],
									 C*x[(n+1)*M+m+1],C*y[(n+1)*M+m+1],C*z[(n+1)*M+m+1],
									 u[n*M+m], u[(n+1)*M+m], u[(n+1)*M+m+1],
									 cs);
				GLContourTriangle_3D(C*x[n*M+m],      C*y[n*M+m],      C*z[n*M+m],
									 C*x[n*M+m+1],    C*y[n*M+m+1],    C*z[n*M+m+1],
									 C*x[(n+1)*M+m+1],C*y[(n+1)*M+m+1],C*z[(n+1)*M+m+1],
									 u[n*M+m], u[n*M+m+1], u[(n+1)*M+m+1],
									 cs);
			}
		}
		if (DrawGrid) {
			glColor3f(0,0,0);
			int N1 = 1 + int(ceil(180/deltatheta_grid));
			if (N1 > 1000) N1 = 1000;
			if (N1 < 5) N1 = 5;
			const int M1 = (N1-1)*2 + 1;
			real *x1 = new real [3*N1*M1];
			if (!x1) {cerr << "*** DisplayDataSphere: out of memory\n"; return;}
			real *y1 = x1 + 1*N1*M1;
			real *z1 = x1 + 2*N1*M1;
			for (n=0; n<N1; n++) for (m=0; m<M1; m++) {
				const real theta = deltatheta_grid*n;
				const real phi = deltaphi_grid*m;
				const real thetaRad = (pi/180)*theta;
				const real phiRad = (pi/180)*phi;
				X[0] = x0 + r*sin(thetaRad)*cos(phiRad);
				X[1] = y0 + r*sin(thetaRad)*sin(phiRad);
				X[2] = z0 + r*cos(thetaRad);
				x1[n*M1+m] = X(0);
				y1[n*M1+m] = X(1);
				z1[n*M1+m] = X(2);
			}
			const real coeff = 1.03;
			for (n=0; n<N1-1; n++) for (m=0; m<M1-1; m++) {
				glBegin(GL_LINE_STRIP);
				vertex3f(coeff*x1[(n+1)*M1+m], coeff*y1[(n+1)*M1+m], coeff*z1[(n+1)*M1+m]);
				vertex3f(coeff*x1[n*M1+m],     coeff*y1[n*M1+m],     coeff*z1[n*M1+m]);
				vertex3f(coeff*x1[n*M1+m+1],   coeff*y1[n*M1+m+1],   coeff*z1[n*M1+m+1]);
				glEnd();
			}
			delete [] x1;
		}
	}
	delete [] u;
}

void TToglWindow::VerticalColorBar(GLfloat x, GLfloat y, GLfloat w, GLfloat h) const
{
	int i;
	const int N = 256;
	const GLfloat dy = h/(N-1);
	const GLfloat C = 1.0/(N-1);
	glBegin(GL_QUAD_STRIP);
	for (i=0; i<N; i++) {
		SetColorFromScaledQuantity(i*C);
		glVertex2f(x,y+i*dy);
		glVertex2f(x+w,y+i*dy);
	}
	glEnd();
}

const int colorbar_area_width = 75 /*65*/;
const int margin_width = 40 /*30*/;
const int title_height = 30;
const int margin_width_3D = 10;

void TToglWindow::VerticalColorBarWithLabels(GLint w0, GLint h0, int tmargin, int bmargin) const
{
	const int ticklen = 12;
	const int inclabel = 1;
	const real bar_x = w0 - colorbar_area_width;
	const real bar_h = 0.8*(h0 - bmargin - tmargin);
	const real bar_y = h0 - tmargin - bar_h;
	const real bar_w = 0.35*colorbar_area_width;
	VerticalColorBar(bar_x, bar_y, bar_w, bar_h);
	glColor3f(foreground_r,foreground_g,foreground_b);
	real datamin_scaled, datamax_scaled;
	const char *const lab = var.selected();
	char s[200];
	const real absdatamax = max(fabs(datamax),fabs(datamin));
	if (absdatamax >= 1000 || absdatamax <= 0.01 && absdatamax != 0) {
		const real scale = pow(10.0,floor(log(absdatamax)/log(10.0)));
		datamin_scaled = datamin/scale;
		datamax_scaled = datamax/scale;
		if (Logarithmic)
			sprintf(s,"log10(%s * %g)",lab,scale);
		else
			sprintf(s,"%s *%g",lab,scale);
	} else {
		datamin_scaled = datamin;
		datamax_scaled = datamax;
		if (Logarithmic)
			sprintf(s,"log10(%s)",lab);
		else
			strcpy(s,lab);
	}
	DrawAxis(bar_x+bar_w,bar_y,0,
			 bar_x+bar_w,bar_y+bar_h,0,
			 ticklen,0,0,
			 datamin_scaled, datamax_scaled, inclabel, helvetica_10);
	glBegin(GL_LINE_STRIP);
	glVertex2d(bar_x+bar_w,bar_y);
	glVertex2d(bar_x,bar_y);
	glVertex2d(bar_x,bar_y+bar_h);
	glVertex2d(bar_x+bar_w,bar_y+bar_h);
	glEnd();

	glListBase(helvetica_12);
	glRasterPos2f(bar_x+2,bar_y-18);
	glCallLists(strlen(s),GL_BYTE,s);
}

static void GLrectframe(GLint x1, GLint y1, GLint x2, GLint y2)
{
	glBegin(GL_LINE_STRIP);
	glVertex2i(x1,y1);
	glVertex2i(x2,y1);
	glVertex2i(x2,y2);
	glVertex2i(x1,y2);
	glVertex2i(x1,y1);
	glEnd();
}

static Tmetagrid *compgradient_gridptr = 0;
static real compgradient_dx = 0;
static Tvariable *compgradient_varptr = 0;

static void compgradient(real x, real y, real z, real& fx, real& fy, real& fz)
{
	smallnat d,dir;
	real val[3][2];
	Tdimvec X;
	for (d=0; d<3; d++) for (dir=0; dir<2; dir++) {
		X[0] = x; X[1] = y; X[2] = z;
		X[d]+= compgradient_dx*(dir==0 ? -1 : +1);
		val[d][dir] = compgradient_gridptr->intpol(X,1,true);
		val[d][dir] = compgradient_varptr->get(*compgradient_gridptr,X);
	}
	const real inv2dx = 0.5/compgradient_dx;
	fx = (val[0][1] - val[0][0])*inv2dx;
	fy = (val[1][1] - val[1][0])*inv2dx;
	fz = (val[2][1] - val[2][0])*inv2dx;
}

static real interpolate(real x, real y, real z)
{
	Tdimvec X(x,y,z);
	compgradient_gridptr->intpol(X,1,true);
	return compgradient_varptr->get(*compgradient_gridptr,X);
}

void TToglWindow::DrawIsoSurface()
{
	if (Nisovalues == 0) return;
	Tmetagrid& g = *gridptr;
	TGridIndex c;
	static bool udens[6] = {false,false,false, false,false,false};
	static const int dir[8][3] = {
		{-1,-1,-1},
		{+1,-1,-1},
		{-1,+1,-1},
		{+1,+1,-1},
		{-1,-1,+1},
		{+1,-1,+1},
		{-1,+1,+1},
		{+1,+1,+1}};
	static const GLfloat red[4] = { 1.0, 0.5, 0.5, 0.0 };
	static const GLfloat blue[4] = { 0.5, 0.5, 1.0, 0.0 };
	static const GLfloat material[4] = {1.0, 1.0, 1.0, 1.0};
	static const GLfloat shininess[1] = {10.0};
	static const GLfloat lightpos0[4] = {+1,0,0,0};
	static const GLfloat lightpos1[4] = {-1,0,0,0};
	static const GLfloat ambient[4] = {0.3,0.3,0.3,1.0};
	real u[8],ufacecenters[6];
	const bool early_ignore = false;
	const real early_ignore_factor = 1.5;
	const real ignore_if_larger = early_ignore_factor*IsoValue[0];
	const real ignore_if_smaller = IsoValue[0]/early_ignore_factor;
	TContourSpec cspec;
	if (Nisovalues == 1) {
		cspec.renew(1.01*IsoValue[0],IsoValue[0],2);
	} else {
		cspec.renew(0,1,Nisovalues);
		int i;
		for (i=0; i<Nisovalues; i++) cspec.setcontour(i,IsoValue[i]);
	}
	int a,d;
	Tdimvec Xc,X;
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_CULL_FACE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambient);
	glLightfv(GL_LIGHT0, GL_POSITION, lightpos0);
	glLightfv(GL_LIGHT1, GL_POSITION, lightpos1);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, red);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, blue);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material);
	if (shiny_isosurfaces) {
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
		const GLfloat reduct = 0.7;
		const GLfloat red_reduced[4] = {red[0]*reduct,red[1]*reduct,red[2]*reduct,red[3]*reduct};
		const GLfloat blue_reduced[4] = {blue[0]*reduct,blue[1]*reduct,blue[2]*reduct,blue[3]*reduct};
		glLightfv(GL_LIGHT0, GL_SPECULAR, red_reduced);
		glLightfv(GL_LIGHT1, GL_SPECULAR, blue_reduced);
	}
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glShadeModel(smooth_isosurfaces ? GL_SMOOTH : GL_FLAT);
	compgradient_gridptr = gridptr;
	compgradient_varptr = &var;
	g.set_intpol_cacheing(true);
	for (c=g.first(); !g.isover(c); c=g.next(c)) {
		if (!g.isleaf(c)) continue;
		const TCellType ct = g.celltype(c);
		if (ct != INTERIOR_CELL) continue;
		if (early_ignore) {
			g > c;
			u[0] = var.get(g,X);
			if (u[0] < ignore_if_smaller || u[0] > ignore_if_larger) continue;
		}
		g.centroid(c,Xc);
		if (!InsideBox(Xc)) continue;
		const real halfdx = 0.5*g.cellsize(c);
		for (a=0; a<8; a++) {
			for (d=0; d<3; d++) X[d] = Xc(d) + dir[a][d]*halfdx;
			g.intpol(X,1,true);
			u[a] = var.get(g,X);
		}
		const smallnat LEFT=0, RIGHT=1;
		udens[0] = g.Nneighbours(c,0,LEFT) > 1;
		udens[1] = g.Nneighbours(c,0,RIGHT) > 1;
		udens[2] = g.Nneighbours(c,1,LEFT) > 1;
		udens[3] = g.Nneighbours(c,1,RIGHT) > 1;
		udens[4] = g.Nneighbours(c,2,LEFT) > 1;
		udens[5] = g.Nneighbours(c,2,RIGHT) > 1;
		if (udens[0]) {
			X[0] = Xc(0) - halfdx;
			X[1] = Xc(1);
			X[2] = Xc(2);
			g.intpol(X,1,true);
			ufacecenters[0] = var.get(g,X);
		} else
			ufacecenters[0] = 0.25*(u[0]+u[2]+u[4]+u[6]);
		if (udens[1]) {
			X[0] = Xc(0) + halfdx;
			X[1] = Xc(1);
			X[2] = Xc(2);
			g.intpol(X,1,true);
			ufacecenters[1] = var.get(g,X);
		} else
			ufacecenters[1] = 0.25*(u[1]+u[3]+u[5]+u[7]);
		if (udens[2]) {
			X[0] = Xc(0);
			X[1] = Xc(1) - halfdx;
			X[2] = Xc(2);
			g.intpol(X,1,true);
			ufacecenters[2] = var.get(g,X);
		} else
			ufacecenters[2] = 0.25*(u[0]+u[1]+u[4]+u[5]);
		if (udens[3]) {
			X[0] = Xc(0);
			X[1] = Xc(1) + halfdx;
			X[2] = Xc(2);
			g.intpol(X,1,true);
			ufacecenters[3] = var.get(g,X);
		} else
			ufacecenters[3] = 0.25*(u[2]+u[3]+u[6]+u[7]);
		if (udens[4]) {
			X[0] = Xc(0);
			X[1] = Xc(1);
			X[2] = Xc(2) - halfdx;
			g.intpol(X,1,true);
			ufacecenters[4] = var.get(g,X);
		} else
			ufacecenters[4] = 0.25*(u[0]+u[1]+u[2]+u[3]);
		if (udens[5]) {
			X[0] = Xc(0);
			X[1] = Xc(1);
			X[2] = Xc(2) + halfdx;
			g.intpol(X,1,true);
			ufacecenters[5] = var.get(g,X);
		} else
			ufacecenters[5] = 0.25*(u[4]+u[5]+u[6]+u[7]);
		compgradient_dx = g.cellsize(c);
		GLIsosurfCube(Xc(0),Xc(1),Xc(2),g.cellsize(c),u,udens,ufacecenters,cspec,
					  &interpolate,(smooth_isosurfaces ? &compgradient : 0));
	}
	g.set_intpol_cacheing(false);
	glPopAttrib();
}

void TToglWindow::RenderVolumetric()
{
	Tmetagrid& g = *gridptr;
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDepthMask(GL_FALSE);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE);
	int Npoints = volumetric_npoints;		// number of points to render
	GLfloat alpha = volumetric_alpha;		// alpha value of points
	if (VolumetricAntialiasing) {
		glEnable(GL_POINT_SMOOTH);
		const GLfloat pointsize = 2;
		glPointSize(pointsize);
		Npoints = int(Npoints/pointsize);
		alpha/= pointsize;
	}
	TGridIndex c;
	Tdimvec Xc;
	int cnt;
	// (1) Integrate the quantity over the grid
	real integral = 0;
	for (cnt=0,c=g.first(); !g.isover(c); c=g.next(c)) {
		if (!g.isleaf(c)) continue;
		const TCellType ct = g.celltype(c);
		if (ct != INTERIOR_CELL) continue;
		g > c;
		g.centroid(c,Xc);
		const real u = var.get(g,Xc);
		const real dx = g.cellsize(c);
		const real vol = dx*dx*dx;
		integral+= vol*max(0.0,u);
	}
	const real vol_per_point = integral/Npoints;
	glColor4f(0.9,0.9,1,alpha);
	int ndone = 0;
	glBegin(GL_POINTS);
	for (cnt=0,c=g.first(); !g.isover(c); c=g.next(c)) {
		if (!g.isleaf(c)) continue;
		const TCellType ct = g.celltype(c);
		if (ct != INTERIOR_CELL) continue;
		g > c;
		g.centroid(c,Xc);
		const real u = var.get(g,Xc);
		const real dx = g.cellsize(c);
		const real vol = dx*dx*dx;
		const int npoints = int(max(0.0,u)*vol/vol_per_point+0.5);
		int i;
		for (i=0; i<npoints; i++) {
			glVertex3f(Xc(0)+dx*(Random()-0.5),Xc(1)+dx*(Random()-0.5),Xc(2)+dx*(Random()-0.5));
			ndone++;
		}
	}
	glEnd();
	glPopAttrib();
}

void TToglWindow::Display(Togl *toglptr)
{
	if (dirty) {
		togl = toglptr;
		helvetica_10 = Togl_LoadBitmapFont(togl,TOGL_BITMAP_HELVETICA_10);
		if (!helvetica_10) {
			cout << "*** Failed to load bitmap font HELVETICA_10, trying 8_BY_13" << endl;
			helvetica_10 = Togl_LoadBitmapFont(togl,TOGL_BITMAP_8_BY_13);
			if (!helvetica_10) cout << "*** Also that failed, small text will not be visible!" << endl;
		}
		helvetica_12 = Togl_LoadBitmapFont(togl,TOGL_BITMAP_HELVETICA_12);
		if (!helvetica_12) {
			cout << "*** Failed to load bitmap font HELVETICA_12, trying 9_BY_15" << endl;
			helvetica_12 = Togl_LoadBitmapFont(togl,TOGL_BITMAP_9_BY_15);
			if (!helvetica_12) cout << "*** Also that failed, medium text will not be visible!" << endl;
		}
		helvetica_18 = Togl_LoadBitmapFont(togl,TOGL_BITMAP_HELVETICA_18);
		if (!helvetica_18) {
			cout << "*** Failed to load bitmap font HELVETICA_18, trying TIMES_ROMAN_24" << endl;
			helvetica_18 = Togl_LoadBitmapFont(togl,TOGL_BITMAP_TIMES_ROMAN_24);
			if (!helvetica_18 && helvetica_12) {
				cout << "*** Also that failed, setting large font = medium font" << endl;
				helvetica_18 = helvetica_12;
			}
			if (!helvetica_18) cout << "*** Also that failed, large text will not be visible!" << endl;
		}
		dirty = false;
	}
	glDisable(GL_LINE_SMOOTH);		// Draw grid lines aliased for speed, and it is actually better because coord-aligned
	glDisable(GL_BLEND);			// Thus can disable blending also
	glDepthMask(GL_FALSE);			// And depth buffer updates
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);	// The default blend function (in case blending is turned on)
	var.select(variable,gamma,invmu0,mass);
	var.logarithmic(Logarithmic);
	var.pseudobackground(pseudobackground);

	glShadeModel( LinearInterpolation ? GL_SMOOTH : GL_FLAT );
	if (WhiteBackground) {
		background_r = background_g = background_b = 1;
		foreground_r = foreground_g = foreground_b = 0;
	} else {
		background_r = background_g = background_b = 0;
		foreground_r = foreground_g = foreground_b = 1;
	}
	glClearColor(background_r,background_g,background_b,0);
//	GLfloat rgba[4];
//	glGetFloatv(GL_COLOR_CLEAR_VALUE,&rgba[0]);
//	cerr << "echo: rgb=" << rgba[0] << "," << rgba[1] << "," << rgba[2] << "\n";
	if (!osflag) glClear(GL_COLOR_BUFFER_BIT);
	if (!gridptr || !gridptr->good()) return;
	Tmetagrid& g = *gridptr;

	const GLint x0 = 0;
	const GLint y0 = 0;
	const GLint w0 = Togl_Width(togl);
	const GLint h0 = Togl_Height(togl);

	glViewport(x0,y0,w0,h0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);  		// Start modifying the projection matrix
	glLoadIdentity();	             	// Reset project matrix
	glOrtho(x0, w0, y0, h0, -1, 1);		// Map abstract coords directly to window coords
	const int mw = (view3D ? margin_width_3D : margin_width);
	lmargin = mw;
	rmargin = mw + (ColorBar ? colorbar_area_width : 0);
	tmargin = mw + (Title ? title_height : 0);
	bmargin = mw;
	real Xmin[3]={0,0,0},Xmax[3]={0,0,0};
	MappedBox(g,Xmin,Xmax);
	const real xmin = Xmin[0], ymin = Xmin[1], zmin = Xmin[2];
	const real xmax = Xmax[0], ymax = Xmax[1], zmax = Xmax[2];
	const int SlicedDim = (g.dimension() == 2) ? 2 : slicedim;

	const int width = w0-lmargin-rmargin;
	const int height = h0-tmargin-bmargin;

	static int cycler = 0;		// Used to cycle through different slices by repeatedly pressing X,Y,Z in 2D viewmode and 3D grid

	slice_xval = 0;		// Used if 2D view and 3D grid
	// Find out datamin,datamax
	bool bypass_min_max_finding = false;
	if (var.hasMinMax() || forced_dataminmax) {
		if (!forced_dataminmax) {
			double m,M;
			var.getMinMax(m,M);
			datamin = m;
			datamax = M;
			one_over_datamax_minus_datamin = 1.0/(datamax - datamin);
		}
		bypass_min_max_finding = true;
		if (verbose) cout << "Bypassing min/max finding\n";
	}
	if (view3D) {
		if ((!can_use_display_list || !display_list_exists) && !bypass_min_max_finding) {
			if (g.dimension() == 3) {
				// 3D view mode, grid is 3D
				// Find slice object min/max
				bool ToOldExtrema = false;
				T3DObject *p;
				for (p=objlist; p; p=p->next)
					if (!strcmp(p->type(),"Slice")) {
						if (!((T3DSliceObject *)p)->DontShowIn3D)
							DataSliceMinMax(g,((T3DSliceObject *)p)->SlicedDim(),
											((T3DSliceObject *)p)->Slice_Xval(), ToOldExtrema);
						ToOldExtrema = true;
					} else if (!strcmp(p->type(),"Sphere")) {
						T3DSphereObject *pp = (T3DSphereObject *)p;
						DisplayDataSphere(g,pp->r,pp->x0,pp->y0,pp->z0, pp->deltatheta,pp->deltatheta_grid,false,ToOldExtrema);
						ToOldExtrema = true;
					}
			} else {
				// 3D view mode, grid is 2D
				DataSliceMinMax(g,2,0);
			}
			var.setMinMax(datamin,datamax);
		}
	} else {
		// 2D view mode
		int XX, YY;
		switch (SlicedDim) {
		case 0: XX=1; YY=2; break;
		case 1: XX=0; YY=2; break;
		default: case 2: XX=0; YY=1; break;
		}
		if (PreserveAspect) {
			const real aspect = real(width)/real(height);
			const real W = max(wmax-wmin, aspect*(hmax-hmin));
			const real H = max(hmax-hmin, (wmax-wmin)/aspect);
			wmin_visible = wmin;
			wmax_visible = wmin+W;
			hmin_visible = hmin;
			hmax_visible = hmin+H;
		} else {
			wmin_visible = wmin;
			wmax_visible = wmax;
			hmin_visible = hmin;
			hmax_visible = hmax;
		}
		if (g.dimension() == 3) {
			// 2D view mode, grid is 3D
			const int N = Nslices(SlicedDim);
			if (N > 0) {
				cycler = (cycler + 1) % N;
				T3DSliceObject *obj = NthSlice(SlicedDim,cycler);
				slice_xval = obj->Slice_Xval();
				if (!forced_dataminmax) {
					DataSliceMinMax(g,SlicedDim,slice_xval);
					var.setMinMax(datamin,datamax);
				}
			}
		} else {
			// 2D view mode, grid is 2D
			if (!bypass_min_max_finding) {
				DataSliceMinMax(g,2,0);
				var.setMinMax(datamin,datamax);
			}
		}
	}
	cs.renew(datamin,datamax,15);
	
	if (!view3D) {
		const int inclabel = 1;
		const int ticklen = 12;
		const double coord_scaling = 1.0;
		glColor3f(foreground_r,foreground_g,foreground_b);
		DrawAxis(lmargin-1,bmargin,0,
				 lmargin-1,h0-tmargin,0,
				 -ticklen,0,0,
				 hmin_visible*coord_scaling,hmax_visible*coord_scaling,inclabel,helvetica_10);
		DrawAxis(w0-rmargin,bmargin,0,
				 w0-rmargin,h0-tmargin,0,
				 ticklen,0,0,
				 hmin_visible*coord_scaling,hmax_visible*coord_scaling,inclabel,helvetica_10);
		DrawAxis(lmargin,bmargin-1,0,
				 w0-rmargin,bmargin-1,0,
				 0,-ticklen,0,
				 wmin_visible*coord_scaling,wmax_visible*coord_scaling,inclabel,helvetica_10);
		DrawAxis(lmargin,h0-tmargin,0,
				 w0-rmargin,h0-tmargin,0,
				 0,ticklen,0,
				 wmin_visible*coord_scaling,wmax_visible*coord_scaling,inclabel,helvetica_10);
	} else {
		glColor3f(foreground_r,foreground_g,0.5*(background_b+foreground_b));
		GLrectframe(lmargin-1,bmargin-1,w0-rmargin,h0-tmargin);
	}

	if (ColorBar) VerticalColorBarWithLabels(w0,h0, tmargin,bmargin);

	if (Title) {
		const char *fntail = filename + strlen(filename) - 1;
		while (fntail >= filename && *fntail != '/') fntail--;
		fntail++;
		char *title = new char [strlen(fntail) + 100];
		if (view3D)
			strcpy(title,fntail);
		else {
			if (g.dimension() == 3) {
				sprintf(title,"%s (%s plane, %s=%g)",
						fntail,slicedim==0 ? "YZ" : (slicedim==1 ? "XZ" : "XY"),
						slicedim==0 ? "X" : (slicedim==1 ? "Y" : "Z"),slice_xval);
			} else {
				sprintf(title,"%s",fntail);
			}
		}
		glColor3f(foreground_r,foreground_g,foreground_b);
		glRasterPos2f(10,h0-22);

		glListBase(helvetica_18);
		glCallLists(strlen(title),GL_BYTE,title);
		delete [] title;
	}
	glViewport(x0+lmargin,y0+bmargin,width,height);
	if (view3D) {
		// In 3D viewing mode; grid can be either 2D or 3D
		glMatrixMode(GL_PROJECTION);  // Start modifying the projection matrix
		glLoadIdentity();
		const real aspect = real(width)/real(height);
		real eye_x, eye_y, eye_z;
		if (XYZmode_in_3D_params) {
			eye_x = eyepointX;
			eye_y = eyepointY;
			eye_z = eyepointZ;
		} else {
			const real pi = M_PI;
			real theta_rad = theta*(pi/180.0);
			if (theta_rad < 0.0001) theta_rad = 0.0001;
			if (theta_rad > pi-0.0001) theta_rad = pi-0.0001;
			const real phi_rad = phi*(pi/180.0);
			const real sintheta = sin(theta_rad);
			const real costheta = cos(theta_rad);
			const real cosphi = cos(phi_rad);
			const real sinphi = sin(phi_rad);
			const real radius_eye = relative_radius*2*max(sqrt(sqr(xmin) + sqr(ymin) + sqr(zmin)),
														  sqrt(sqr(xmax) + sqr(ymax) + sqr(zmax)));
			eye_x = cx + radius_eye * sintheta * cosphi;
			eye_y = cy + radius_eye * sintheta * sinphi;
			eye_z = cz + radius_eye * costheta;
		}
		const real dist1_2 = sqr(eye_x - xmin) + sqr(eye_y - ymin) + sqr(eye_z - zmin);
		const real dist2_2 = sqr(eye_x - xmin) + sqr(eye_y - ymin) + sqr(eye_z - zmax);
		const real dist3_2 = sqr(eye_x - xmin) + sqr(eye_y - ymax) + sqr(eye_z - zmin);
		const real dist4_2 = sqr(eye_x - xmin) + sqr(eye_y - ymax) + sqr(eye_z - zmax);
		const real dist5_2 = sqr(eye_x - xmax) + sqr(eye_y - ymin) + sqr(eye_z - zmin);
		const real dist6_2 = sqr(eye_x - xmax) + sqr(eye_y - ymin) + sqr(eye_z - zmax);
		const real dist7_2 = sqr(eye_x - xmax) + sqr(eye_y - ymax) + sqr(eye_z - zmin);
		const real dist8_2 = sqr(eye_x - xmax) + sqr(eye_y - ymax) + sqr(eye_z - zmax);
		const real zFar = 1.3*sqrt(max(max(max(dist1_2,dist2_2),max(dist3_2,dist4_2)),
									   max(max(dist5_2,dist6_2),max(dist7_2,dist8_2))));
		const real zNear = 0.001*zFar;
		gluPerspective(field_of_view,aspect,zNear,zFar);
		glMatrixMode(GL_MODELVIEW);		// Start modifying the viewing matrix
		glLoadIdentity();				// Reset viewing matrix
		gluLookAt(eye_x,eye_y,eye_z, cx,cy,cz, 0,0,1);
		if (g.dimension()==3) {
			// In 3D viewing mode and grid is 3D
			glPushAttrib(GL_ALL_ATTRIB_BITS);
			if (antialiased_lines) {
				glEnable(GL_LINE_SMOOTH);			// Draw antialised field lines
			}
			glEnable(GL_BLEND);					// Must enable blending also (both for antialiased lines and transparency)
			glEnable(GL_DEPTH_TEST);			// We need to use depth buffer in 3D mode (if grid is 3D)
			glDepthMask(GL_TRUE);				// Depth buffer updates on
			glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);	// Usual transparency setting
			glClear(GL_DEPTH_BUFFER_BIT);		// Zero depth buffer
		}
		if (can_use_display_list && display_list_exists) {
			if (verbose) cout << "Calling display list\n" << flush;
			glCallList(DisplayListID);
		} else {
			if (display_list_exists) glDeleteLists(DisplayListID,1);
			DisplayListID = glGenLists(1);
			glNewList(DisplayListID,GL_COMPILE_AND_EXECUTE);
			if (g.dimension() == 3) {
				// In 3D viewing mode, and grid is 3D
				// Drawing order: 1) Isosurface, 2) Spheres,
				// 3) Non-antialiased (or aliased but in transparent mode) field lines
				// 4) Slices, 5) Antialiased field lines in opaque mode
				// 6) Volumetric data
				T3DObject *p;
				// (1) Isosurfaces
				if (DrawIsosurf) DrawIsoSurface();
				// (2) Spheres
				const real AlphaValue_sav = AlphaValue;
				AlphaValue = 1.0;
				for (p=objlist; p; p=p->next) {
					if (!strcmp(p->type(),"Sphere")) {
						T3DSphereObject *pp = (T3DSphereObject *)p;
						DisplayDataSphere(g,pp->r,pp->x0,pp->y0,pp->z0,pp->deltatheta,pp->deltatheta_grid);
					}
				}
				AlphaValue = AlphaValue_sav;
				// (3) Non-antialiased field lines, or aliased in transparent mode
				if (!antialiased_lines || AlphaValue < 1) {
					if (antialiased_lines) {
						glPushAttrib(GL_LINE_BIT);
						glLineWidth(antialiased_FL_width);
						glBlendFunc(GL_SRC_ALPHA,GL_ONE);
						glDepthMask(GL_FALSE);
					}
					for (p=objlist; p; p=p->next) {
						if (!strcmp(p->type(),"FieldLineBunch")) {
							T3DFieldLineBunchObject *pp = (T3DFieldLineBunchObject *)p;
							DrawFieldLineBunch(g,pp->n, pp->r1, pp->r2, pp->VectorField, pp->TraceDirection,
											   pp->Distribution, pp->LoopThresholdType, pp->LoopThreshold);
						}
					}
					if (antialiased_lines) {
						glDepthMask(GL_TRUE);
						glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);			// normal setting
						glPopAttrib();		// restore old line width
					}
				}
				// (4) Slices
				for (p=objlist; p; p=p->next) {
					if (!strcmp(p->type(),"Slice")) {
						if (!((T3DSliceObject *)p)->DontShowIn3D)
							DisplayData2D(g,((T3DSliceObject *)p)->SlicedDim(),((T3DSliceObject *)p)->Slice_Xval());
					}
				}
				// (5) Antialiased field lines in opaque mode
				if (antialiased_lines && AlphaValue == 1) {
					glPushAttrib(GL_LINE_BIT);
					glLineWidth(antialiased_FL_width);
					if (AlphaValue < 1) {
						glBlendFunc(GL_SRC_ALPHA,GL_ONE);
					}
					glDepthMask(GL_FALSE);
					for (p=objlist; p; p=p->next) {
						if (!strcmp(p->type(),"FieldLineBunch")) {
							T3DFieldLineBunchObject *pp = (T3DFieldLineBunchObject *)p;
							DrawFieldLineBunch(g,pp->n, pp->r1, pp->r2, pp->VectorField, pp->TraceDirection,
											   pp->Distribution, pp->LoopThresholdType, pp->LoopThreshold);
						}
					}
					glDepthMask(GL_TRUE);
					if (AlphaValue < 1) {
						glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);			// normal setting
					}
					glPopAttrib();		// restore old line width
				}
				// (6) Volumetric data
				if (DrawVolumetric) RenderVolumetric();
			} else {
				// In 3D viewing mode, and grid is 2D
				DisplayData2D(g);
			}
			glEndList();
			display_list_exists = true;
		}
		if (g.dimension()==3) {
			// In 3D viewing mode, and grid is 3D
			glPopAttrib();
		}
	} else {
		// In 2D viewing mode; grid can be either 2D or 3D
		glMatrixMode(GL_PROJECTION);  // Start modifying the projection matrix
		glLoadIdentity();             // Reset project matrix
		if (PreserveAspect) {
			glOrtho(wmin_visible,wmax_visible, hmin_visible,hmax_visible, -1,1);
		} else {
			glOrtho(wmin, wmax, hmin, hmax, -1, 1);
		}
		glMatrixMode(GL_MODELVIEW);	// Load identity matrix to modelview matrix in 2D case
		glLoadIdentity();
		if (g.dimension() == 3) {
			// In 2D viewing mode, and grid is 3D
			vert3_select = SlicedDim;
			const int N = Nslices(SlicedDim);
			if (N > 0) {
				T3DSliceObject *obj = NthSlice(SlicedDim,cycler);
				slice_xval = obj->Slice_Xval();
				DisplayData2D(g,SlicedDim,slice_xval  /* 0.5*(Xmin[SlicedDim]+Xmax[SlicedDim]) */  );
			}
			vert3_select = -1;
		} else {
			// In 3D viewing mode, and grid is 2D
			DisplayData2D(g);
		}
	}

	can_use_display_list = false;
}

void TToglWindow::getbox(real xmin[3], real xmax[3], int& dim) const
{
	if (!gridptr || !gridptr->good()) return;
	MappedBox(*gridptr,xmin,xmax);
	dim = gridptr->dimension();
}

bool TToglWindow::EvalVectorField(TVectorField VectorField, Tmetagrid& g, const Tdimvec& r, Tdimvec& V)
{
	if (!g.intpol(r,1,true)) return false;
	const smallnat ncd = g.Ncelldata();
	switch (VectorField) {
	case VectorFieldRhoV:
		V[0] = g.CT(1);
		V[1] = g.CT(2);
		V[2] = g.CT(3);
		break;
	case VectorFieldB:
		if (ncd < 8)
			V[0] = V[1] = V[2] = 0;
		else if (ncd == 8) {
			V[0] = g.CT(5);
			V[1] = g.CT(6);
			V[2] = g.CT(7);
		} else {
			V[0] = g.CT(5) + g.CT(8);
			V[1] = g.CT(6) + g.CT(9);
			V[2] = g.CT(7) + g.CT(10);
		}
		break;
	case VectorFieldB0:
		if (ncd < 11)
			V[0] = V[1] = V[2] = 0;
		else {
			V[0] = g.CT(8);
			V[1] = g.CT(9);
			V[2] = g.CT(10);
		}
		break;
	case VectorFieldB1:
		if (ncd < 8)
			V[0] = V[1] = V[2] = 0;
		else {
			V[0] = g.CT(5);
			V[1] = g.CT(6);
			V[2] = g.CT(7);
		}
		break;
	case VectorFieldJ:
		if (ncd < 8)
			V[0] = V[1] = V[2] = 0;
		else  {
			real jx,jy,jz;
			ComputeCurl(g,r,5,jx,jy,jz);
			V[0] = jx;
			V[1] = jy;
			V[2] = jz;
		}
		break;
	case VectorFieldPoynting:
		if (ncd < 8)
			V[0] = V[1] = V[2] = 0;
		else {
			double Bx,By,Bz;
			if (ncd < 11) {
				Bx = g.CT(5);
				By = g.CT(6);
				Bz = g.CT(7);
			} else {
				Bx = g.CT(5) + g.CT(8);
				By = g.CT(6) + g.CT(9);
				Bz = g.CT(7) + g.CT(10);
			}
			const double invrho = 1.0/g.CT(0);
			const double vx = g.CT(1)*invrho;
			const double vy = g.CT(2)*invrho;
			const double vz = g.CT(3)*invrho;
			const double invmu0 = 795775.0;
			V[0] = invmu0*(sqr(By)*vx - Bx*By*vy + Bz*(Bz*vx - Bx*vz));
			V[1] = invmu0*(-Bx*By*vx + sqr(Bx)*vy + Bz*(Bz*vy - By*vz));
			V[2] = invmu0*(-Bx*Bz*vx + sqr(Bx)*vz + By*(-Bz*vy + By*vz));
		}
		break;
	case VectorFieldEnergyFlux:
	{
		const double invmu0 = 795775.0;
		const double gamma = 5.0/3.0;
		double P,Utot,Bx,By,Bz;
		const double invrho = 1.0/g.CT(0);
		if (ncd == 5) {
			P = (gamma-1)*(g.CT(4) - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3))*0.5*invrho));
			Utot = g.CT(4);
			Bx = By = Bz = 0;
		} else {
			P = (gamma-1)*(g.CT(4) - (sqr(g.CT(1)) + sqr(g.CT(2)) + sqr(g.CT(3)))*0.5*invrho
						   - 0.5*invmu0*(sqr(g.CT(5)) + sqr(g.CT(6)) + sqr(g.CT(7))));
			Utot = g.CT(4);
			Bx = g.CT(5);
			By = g.CT(6);
			Bz = g.CT(7);
			if (ncd >= 11) {
				Utot+= 0.5*invmu0*(sqr(g.CT(8)) + sqr(g.CT(9)) + sqr(g.CT(10)))
					+ (g.CT(5)*g.CT(8) + g.CT(6)*g.CT(9) + g.CT(7)*g.CT(10))*invmu0;
				Bx+= g.CT(8);
				By+= g.CT(9);
				Bz+= g.CT(10);
			}
		}
		const double scalar = Utot + P - 0.5*invmu0*(sqr(Bx) + sqr(By) + sqr(Bz));
		const double vx = g.CT(1)*invrho;
		const double vy = g.CT(2)*invrho;
		const double vz = g.CT(3)*invrho;
		const double Sx = invmu0*(sqr(By)*vx - Bx*By*vy + Bz*(Bz*vx - Bx*vz));
		const double Sy = invmu0*(-Bx*By*vx + sqr(Bx)*vy + Bz*(Bz*vy - By*vz));
		const double Sz = invmu0*(-Bx*Bz*vx + sqr(Bx)*vz + By*(-Bz*vy + By*vz));
		const double Kx = Sx + vx*scalar;
		const double Ky = Sy + vy*scalar;
		const double Kz = Sz + vz*scalar;
		V[0] = Kx;
		V[1] = Ky;
		V[2] = Kz;
		break;
	}
	}
	return true;
}

double FieldLineStopSphereRadius = 0;

void TToglWindow::DrawFieldLine(Tmetagrid& g, const double r0[3], int sgn, TVectorField VectorField,
								TLoopThresholdType LoopThresholdType, double LoopThreshold)
{
	int i;
	smallnat d;
	Tdimvec V,r,dr,newr;
	const int maxcnt = 500;
	const double relstep = 0.3*sgn;
	for (d=0; d<3; d++) r[d] = r0[d];
	switch (VectorField) {
	case VectorFieldRhoV:
		glColor3f(0.7,1,1);
		break;
	case VectorFieldPoynting:
		glColor3f(1,0.5,0.5);
		break;
	case VectorFieldB:
	case VectorFieldB1:
	case VectorFieldB0:
		glColor3f(1,1,0.2);
		break;
	case VectorFieldJ:
		glColor3f(1,0.6,1);
		break;
	case VectorFieldEnergyFlux:
		glColor3f(1,0.5,0.5);
		break;
	}
	double maxdist = 0;
	double mindist = LoopThreshold;
	switch (LoopThresholdType) {
	case LocalSpacingLoopThreshold:
	{
		const TGridIndex i = g.find(r);
		if (i != NOINDEX) mindist*= g.cellsize(i);
	}
	break;
	case BasegridSpacingLoopThreshold:
		mindist*= g.BasegridCellSize();
		break;
	case AbsoluteLoopThreshold:
		break;
	}
	glBegin(GL_LINE_STRIP);
	glVertex3f(r(0),r(1),r(2));
	for (i=0; i<maxcnt; i++) {
		const double dist = sqrt(sqr(r(0)-r0[0]) + sqr(r(1)-r0[1]) + sqr(r(2)-r0[2]));
		if (dist > maxdist) maxdist = dist;
		if (i > 0 && dist < maxdist && dist < mindist) {
			cout << "toglwin.C:DrawFieldLine: Loop detected\n";
			break;
		}
		if (!EvalVectorField(VectorField,g,r,V)) {
			break;
		}
		double Vlen = sqrt(sqr(V[0]) + sqr(V[1]) + sqr(V[2]));
		if (Vlen == 0) {
			break;		// we have no possibility but to stop following the field line if V==(0,0,0) at any point
		}
		double invVlen = 1.0/Vlen;
		V[0] = V(0)*invVlen;
		V[1] = V(1)*invVlen;
		V[2] = V(2)*invVlen;
		const TGridIndex c = g.find(r);
		if (g.celltype(c) != INTERIOR_CELL) {
			break;
		}
		const double step = g.cellsize(c)*relstep;
		dr[0] = (0.5*step)*V(0);
		dr[1] = (0.5*step)*V(1);
		dr[2] = (0.5*step)*V(2);
		newr[0] = r(0) + dr(0);
		newr[1] = r(1) + dr(1);
		newr[2] = r(2) + dr(2);
		if (!EvalVectorField(VectorField,g,newr,V)) break;
		Vlen = sqrt(sqr(V[0]) + sqr(V[1]) + sqr(V[2]));
		if (Vlen == 0) {
			break;		// we have no possibility but to stop following the field line if V==(0,0,0) at any point
		}
		invVlen = 1.0/Vlen;
		V[0] = V(0)*invVlen;
		V[1] = V(1)*invVlen;
		V[2] = V(2)*invVlen;
		dr[0] = step*V(0);
		dr[1] = step*V(1);
		dr[2] = step*V(2);
		newr[0] = r(0) + dr(0);
		newr[1] = r(1) + dr(1);
		newr[2] = r(2) + dr(2);
		glVertex3f(newr(0),newr(1),newr(2));
		r[0] = newr(0);
		r[1] = newr(1);
		r[2] = newr(2);
		if (sqr(r[0]) + sqr(r[1]) + sqr(r[2]) < sqr(FieldLineStopSphereRadius)) {
			break;
		}
	}
	glEnd();
}

void TToglWindow::DrawFieldLineBunch(Tmetagrid& g, int n, const double r1[3], const double r2[3],
									 TVectorField VectorField, TTraceDirection TraceDirection,
									 TFieldLineDistribution Distribution,
									 TLoopThresholdType LoopThresholdType, double LoopThreshold)
{
	if (n < 1) return;
	if (verbose) cout << "DrawFieldLineBunch(n=" << n << ")\n";
	double *const t = new double [n];
	const double dt = 1.0/((n==1 ? 2 : n)-1);
	double t1;
	int i;
	for (i=0,t1=0; i<n; i++,t1+=dt) t[i] = t1;
	if (Distribution != UniformDistribution) {
		const int N = max(200,n);
		double *const lambda = new double [N];
		const double dt1 = 1.0/(N-1);
		double weight;
		Tdimvec r,V;
		smallnat d;
		lambda[0] = 0;
		int validN = 0;
		double aveweight = 0;
		for (i=1,t1=dt1; i<N; i++,t1+=dt1) {
			for (d=0; d<3; d++)
				r[d] = (1-t1)*r1[d] + t1*r2[d];
			const bool isvalid = EvalVectorField(VectorField,g,r,V);
			weight = sqrt(sqr(V[0]) + sqr(V[1]) + sqr(V[2]));
			if (isvalid) {aveweight+= weight; validN++;}
			if (Distribution == InvSqrtDistribution) weight = sqrt(weight);
			lambda[i] = lambda[i-1] + dt1*weight;
		}
		if (validN > 0) {
			aveweight/= validN;
			for (i=1,t1=dt1; i<N; i++,t1+=dt1) {
				for (d=0; d<3; d++)
					r[d] = (1-t1)*r1[d] + t1*r2[d];
				if (EvalVectorField(VectorField,g,r,V))
					weight = sqrt(sqr(V[0]) + sqr(V[1]) + sqr(V[2]));
				else
					weight = aveweight;
				if (Distribution == InvSqrtDistribution) weight = sqrt(weight);
				lambda[i] = lambda[i-1] + dt1*weight;
			}
		} else {
			cerr << "toglwin.C:DrawFieldLineBunch warning: No field lines inside box\n";
			delete [] lambda;
			return;
		}
		const double invC = 1.0/lambda[N-1];
		for (i=0; i<N; i++) lambda[i]*= invC;
		// now lambda[0] == 0, lambda[N-1] == 1
		// "invert" lambda vector into t vector
		double lam;
		for (i=0,lam=0; i<n; i++,lam+=dt) {
			// find largest lambda[j] such that lambda[j] <= lam
			int j;
			for (j=N-1; j>=0; j--) if (lambda[j] <= lam) break;
			if (j < 0) break;		// for safety
			t[i] = double(j)/double(N);
		}
		delete [] lambda;
	}
	smallnat d;
	double r[3];
	g.set_intpol_cacheing(true);
	for (i=0; i<n; i++) {
		for (d=0; d<3; d++) r[d] = (1-t[i])*r1[d] + t[i]*r2[d];
		switch (TraceDirection) {
		case NegativeTraceDirection:
			DrawFieldLine(g,r,-1,VectorField,LoopThresholdType,LoopThreshold);
			break;
		case BothTraceDirections:
			DrawFieldLine(g,r,-1,VectorField,LoopThresholdType,LoopThreshold);
		case PositiveTraceDirection:
			DrawFieldLine(g,r,+1,VectorField,LoopThresholdType,LoopThreshold);
			break;
		}
	}
	delete [] t;
	g.set_intpol_cacheing(false);
}

void TToglWindow::ToModelCoords(int mx, int my, double& w, double& h, double& x, double& y, double& z)
{
	const int current_w = Togl_Width(togl) - lmargin - rmargin;
	const int current_h = Togl_Height(togl) - tmargin - bmargin;
	mx-= lmargin;
	my-= tmargin;
	w = wmin_visible + (mx/double(current_w))*(wmax_visible - wmin_visible);
	h = hmin_visible + ((current_h-my)/double(current_h))*(hmax_visible - hmin_visible);
	if (gridptr->dimension() == 2) {
		x = w;
		y = h;
		z = 0;
	} else {
		switch (slicedim) {
		case 0:
			x = slice_xval;
			y = w;
			z = h;
			break;
		case 1:
			x = w;
			y = slice_xval;
			z = h;
			break;
		case 2:
			x = w;
			y = h;
			z = slice_xval;
			break;
		}
	}
}

void TToglWindow::Zoom(real wmin1, real hmin1, real wmax1, real hmax1)
{
	wmin = wmin1;
	hmin = hmin1;
	wmax = wmax1;
	hmax = hmax1;
	real tmp;
	if (hmin > hmax) {tmp=hmin; hmin=hmax; hmax=tmp;}
	if (wmin > wmax) {tmp=wmin; wmin=wmax; wmax=tmp;}
	needs_redraw = true;
}

void TToglWindow::FullZoom()
{
	real Xmin[3]={0,0,0},Xmax[3]={0,0,0};
	int XX, YY;
	switch (slicedim) {
	case 0: XX=1; YY=2; break;
	case 1: XX=0; YY=2; break;
	default: case 2: XX=0; YY=1; break;
	}
	MappedBox(*gridptr,Xmin,Xmax);
	wmin_global = Xmin[XX]; wmax_global = Xmax[XX];
	hmin_global = Xmin[YY]; hmax_global = Xmax[YY];
	if (wmin != wmin_global || wmax != wmax_global || hmin != hmin_global || hmax != hmax_global) {
		wmin=wmin_global; wmax=wmax_global;
		hmin=hmin_global; hmax=hmax_global;
		needs_redraw = true;
	}
}

void TToglWindow::GetCurrentZoom(double& wmin1, double& hmin1, double& wmax1, double& hmax1) const
{
	wmin1 = wmin;
	hmin1 = hmin;
	wmax1 = wmax;
	hmax1 = hmax;
}

bool TToglWindow::ExportASCII(const char *fn) const
{
	if (view3D) {
		cout << "TToglWindow::ExportASCII: Export ASCII works only in 2D viewing mode\n";
		return false;
	}
	FILE *fp;
	fp = fopen(fn,"w");
	if (!fp) return false;
	fprintf(fp,"# x y z %s\n",var.selected());
	int i,j;
	double x,y,z,w,h;
	const real dw = gridptr->MinimumGridSpacing();
	const real dh = dw;
	const int nx = int((wmax-wmin)/dw + 0.5);
	const int ny = int((hmax-hmin)/dh + 0.5);
	cout << "ExportASCII: writing " << nx << "x" << ny << " matrix to \"" << fn << "\" ..";
	Tdimvec X;
	const bool is2d = gridptr->dimension() == 2;
	for (i=0; i<nx; i++) for (j=0; j<ny; j++) {
		w = wmin + dw*i;
		h = hmin + dh*j;
		if (is2d) {
			x = w;
			y = h;
			z = 0;
		} else {
			switch (slicedim) {
			case 0:
				x = slice_xval;
				y = w;
				z = h;
				break;
			case 1:
				x = w;
				y = slice_xval;
				z = h;
				break;
			default: case 2:
				x = w;
				y = h;
				z = slice_xval;
				break;
			}
		}
		X[0] = x; X[1] = y; X[2] = z;
		gridptr->intpol(X,1,true);
		fprintf(fp,"%g %g %g %g\n",x,y,z,double(var.get(*gridptr,X)));
	}
	fclose(fp);
	cout << " done\n";
	return true;
}


static GLvoid *grabPixels(unsigned int width, unsigned int height)
{
   GLvoid *buffer;
   GLint swapbytes, lsbfirst, rowlength;
   GLint skiprows, skippixels, alignment;
   GLenum format;
   unsigned int size;
   const int inColor = 1;

   if (inColor) {
      format = GL_RGB;
      size = width * height * 3;
   } else {
      format = GL_LUMINANCE;
      size = width * height * 1;
   }

   buffer = (GLvoid *) malloc(size);
   if (buffer == NULL)
      return NULL;

   /* Save current modes. */
   glGetIntegerv(GL_PACK_SWAP_BYTES, &swapbytes);
   glGetIntegerv(GL_PACK_LSB_FIRST, &lsbfirst);
   glGetIntegerv(GL_PACK_ROW_LENGTH, &rowlength);
   glGetIntegerv(GL_PACK_SKIP_ROWS, &skiprows);
   glGetIntegerv(GL_PACK_SKIP_PIXELS, &skippixels);
   glGetIntegerv(GL_PACK_ALIGNMENT, &alignment);
   /* Little endian machines (DEC Alpha for example) could
      benefit from setting GL_PACK_LSB_FIRST to GL_TRUE
      instead of GL_FALSE, but this would require changing the
      generated bitmaps too. */
   glPixelStorei(GL_PACK_SWAP_BYTES, GL_FALSE);
   glPixelStorei(GL_PACK_LSB_FIRST, GL_FALSE);
   glPixelStorei(GL_PACK_ROW_LENGTH, 0);
   glPixelStorei(GL_PACK_SKIP_ROWS, 0);
   glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
   glPixelStorei(GL_PACK_ALIGNMENT, 1);

   /* Actually read the pixels. */
   glReadPixels(0, 0, width, height, format,
                GL_UNSIGNED_BYTE, (GLvoid *) buffer);

   /* Restore saved modes. */
   glPixelStorei(GL_PACK_SWAP_BYTES, swapbytes);
   glPixelStorei(GL_PACK_LSB_FIRST, lsbfirst);
   glPixelStorei(GL_PACK_ROW_LENGTH, rowlength);
   glPixelStorei(GL_PACK_SKIP_ROWS, skiprows);
   glPixelStorei(GL_PACK_SKIP_PIXELS, skippixels);
   glPixelStorei(GL_PACK_ALIGNMENT, alignment);
   return buffer;
}


bool TToglWindow::ExportPPM(const char *fn, unsigned int width, unsigned int height) const
{
	unsigned char *buff = (unsigned char *)grabPixels(width,height);
	if (!buff) return false;
	FILE *fp = fopen(fn,"w");
	fprintf(fp,"P6\n");
	fprintf(fp,"%d %d 255\n",int(width),int(height));
	int i,j;
	for (i=int(height)-1; i>=0; i--) for (j=0; j<int(width); j++) {
		putc(buff[(i*width+j)*3 + 0],fp);
		putc(buff[(i*width+j)*3 + 1],fp);
		putc(buff[(i*width+j)*3 + 2],fp);
	}
	fclose(fp);
	free(buff);
	return true;
}
