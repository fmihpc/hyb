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

#include "toglwin.H"
#include <iostream>
#include <cstdlib>
#include <cstring>
using namespace std;
#include <tcl.h>
#include <tk.h>

// The following variable is a special hack that is needed in order for
// Sun shared libraries to be used for Tcl.
//extern "C" int matherr();
//int *tclDummyMathPtr = (int *) matherr;
extern real GridDimensionScaling;		// defined in gridcache.C

struct THcvisFlags {
	bool threeD;
	bool verbose;
        int intpol_order;
	THcvisFlags() {			// these are the defaults
		threeD = false;
		verbose = false;
	        intpol_order = 1;
	}
	void tclset(Tcl_Interp *interp);
};

void THcvisFlags::tclset(Tcl_Interp *interp)
{
	Tcl_SetVar(interp,const_cast<char*>("threeD"),const_cast<char*>(threeD ? "1" : "0"),0);
	Tcl_SetVar(interp,const_cast<char*>("verbose"),const_cast<char*>(verbose ? "1" : "0"),0);
        Tcl_SetVar(interp,const_cast<char*>("intpol_order"),const_cast<char*>(intpol_order ? "1" : "0"),0);
}

static THcvisFlags HcvisFlags;
static Tcl_DString FileNameList;
static char *InitFile = 0;

static void usage()
{
	clog << "usage: hcvis [-3d] [-2d] [-v] [-z] [-s] [-scale scale] [-init init.tcl] [*.hc]\n";
}

/*
 * Togl widget create callback.  This is called by Tcl/Tk when the widget has
 * been realized.  Here's where one may do some one-time context setup or
 * initializations.
 */
void create_cb(Togl *togl)
{
//	FontBase = Togl_LoadBitmapFont(togl, TOGL_BITMAP_8_BY_13);
//	if (!FontBase) {cerr << "*** Couldn't load font!\n"; exit(1);}
	TToglWindow::New(Togl_Ident(togl));
}

void destroy_cb(Togl *togl)
{
	TToglWindow::Delete(Togl_Ident(togl));
}

/*
 * Togl widget reshape callback.  This is called by Tcl/Tk when the widget
 * has been resized.  Typically, we call glViewport and perhaps setup the
 * projection matrix.
 */
void reshape_cb(Togl *togl)
{
	TToglWindow *const winptr = TToglWindow::Find(Togl_Ident(togl));
	const int width = Togl_Width(togl);
	const int height = Togl_Height(togl);
	if (width != winptr->oldwidth || height != winptr->oldheight) {
		if (HcvisFlags.verbose) cout << "- Slow reshape callback\n" << flush;
		winptr->needs_redraw = true;
		winptr->can_use_display_list = true;
		winptr->oldwidth = width;
		winptr->oldheight = height;
	} else {
		if (HcvisFlags.verbose) cout << "- Fast reshape callback\n" << flush;
	}
	glViewport(0,0,width,height);
}

/*
 * Togl widget display callback.  This is called by Tcl/Tk when the widget's
 * contents have to be redrawn.  Typically, we clear the color and depth
 * buffers, render our objects, then swap the front/back color buffers.
 */
void display_cb(Togl *togl)
{
	TToglWindow *const winptr = TToglWindow::Find(Togl_Ident(togl));
	if (!winptr) {cerr << "*** display_db: TToglWindow::Find() failed, internal error\n"; exit(1);}
	if (winptr->needs_redraw) {
		if (HcvisFlags.verbose) cout << "- Slow display callback\n" << flush;
		winptr->verbose = HcvisFlags.verbose;
		winptr->Display(togl);
		winptr->needs_redraw = false;
	} else {
		if (HcvisFlags.verbose) cout << "- Fast display callback\n" << flush;
	}
	Togl_SwapBuffers(togl);
}

static int a2i(const char *s) {if (!s) return 0; else return atoi(s);}
static double a2f(const char *s) {if (!s) return 0; else return atof(s);}

static int hcsetparams_cb(Togl *togl, int argc, char * /*argv*/ [])
// usage: $win hcsetparams, will copy values from global arrays Flags(win,:) and FileNameIndex(win)
// it also copies Params(win,var), Params(win,alpha), and Params(win,slicedim)
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 2) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcsetparams\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcsetparams: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}
	char varname[256];

	sprintf(varname,"FileNameIndex(%s)",win);
	const int i = a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY));
	char cmd[256];
	sprintf(cmd,"lindex $FileNameList %d",i);
	if (Tcl_GlobalEval(interp,cmd) != TCL_OK) {
		cerr << "*** Tcl_Eval in hcsetparams_cb returned error code\n";
	}

	bool newbool;
	char *newstring;
	int newint;

	winptr->verbose = HcvisFlags.verbose;
	winptr->SetFilename(interp->result);
	
	sprintf(varname,"Flags(%s,LinearInterpolation)",win);
	winptr->LinearInterpolation = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,DrawGrid)",win);
	winptr->DrawGrid = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,DrawGhosts)",win);
	winptr->DrawGhosts = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,HideGhosts)",win);
	winptr->HideGhosts = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,DrawDead)",win);
	winptr->DrawDead = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,DrawCont)",win);
	winptr->DrawCont = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,DrawIsosurf)",win);
	winptr->DrawIsosurf = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,DrawVolumetric)",win);
	winptr->DrawVolumetric = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,VolumetricAntialiasing)",win);
	winptr->VolumetricAntialiasing = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,AntialiasedLines)",win);
	winptr->antialiased_lines = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,SmoothIsosurfaces)",win);
	winptr->smooth_isosurfaces = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,ShinyIsosurfaces)",win);
	winptr->shiny_isosurfaces = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,WhiteBackground)",win);
	winptr->WhiteBackground = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,view3D)",win);
	newbool        = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));
	if (newbool != winptr->view3D) {
		winptr->view3D = newbool;
		winptr->ClearMinMax("view3D");
	}
	
	sprintf(varname,"Flags(%s,Logarithmic)",win);
	newbool             = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));
	if (newbool != winptr->Logarithmic) {
		winptr->Logarithmic = newbool;
		winptr->ClearMinMax("Logarithmic");
	}
	
	sprintf(varname,"Flags(%s,PreserveAspect)",win);
	winptr->PreserveAspect = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,Mapping)",win);
	winptr->Mapping = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,ColorBar)",win);
	winptr->ColorBar = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Flags(%s,Title)",win);
	winptr->Title = bool(a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Params(%s,var)",win);
	newstring = (char *)Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY);
	if (winptr->variable && strcmp(newstring,winptr->variable)) {
		winptr->variable = newstring;
		winptr->ClearMinMax("variable change");
	}

	sprintf(varname,"Params(%s,alpha)",win);
	winptr->AlphaValue = a2f(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY));

	int j;
	winptr->Nisovalues = 0;
	for (j=0; j<MAX_ISO_VALUES; j++) {
		sprintf(varname,"Params(%s,IsoValue%d)",win,j);
		char *const s = (char *)Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY);
		if (!s || !*s) break;
		winptr->IsoValue[winptr->Nisovalues++] = a2f(s);
	}

	sprintf(varname,"Params(%s,VolumetricNpoints)",win);
	winptr->volumetric_npoints = int(a2f(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY)));

	sprintf(varname,"Params(%s,VolumetricAlpha)",win);
	winptr->volumetric_alpha = a2f(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY));

	sprintf(varname,"Params(%s,FieldLineStopSphereRadius)",win);
	FieldLineStopSphereRadius = a2f(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY));

	sprintf(varname,"Params(%s,slicedim)",win);
	newint = a2i(Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY));
	if (newint != winptr->slicedim) {
		winptr->slicedim = newint;
		winptr->FullZoom();
		winptr->ClearMinMax("slicedim change");
	}

	sprintf(varname,"Params(%s,min)",win);
	char *const minstr = (char *)Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY);
	sprintf(varname,"Params(%s,max)",win);
	char *const maxstr = (char *)Tcl_GetVar(interp,varname,TCL_GLOBAL_ONLY);
	if (minstr && maxstr && *minstr && *maxstr) {
		const double datamin = a2f(minstr);
		const double datamax = a2f(maxstr);
		if (datamin != winptr->datamin || datamax != winptr->datamax) {
			winptr->ClearMinMax("data min/max set");
			winptr->datamin = a2f(minstr);
			winptr->datamax = a2f(maxstr);
			winptr->one_over_datamax_minus_datamin = 1.0/(winptr->datamax - winptr->datamin);
			winptr->forced_dataminmax = true;
		}
	} else
		winptr->forced_dataminmax = false;
	
	winptr->needs_redraw = true;
	if (HcvisFlags.verbose) cout << "- hcsetparams_cb\n" << flush;

	return TCL_OK;
}

// usage: $win hcsetpalette palette_filename, will read in the palette from the raw palette file
static int hcsetpalette_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 3) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcsetpalette palettefile\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcsetpalette: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}
	if (strlen(argv[2]) == 0 || !strcmp(argv[2],"Grayscale"))
		winptr->palette.grayscale();
	else
		winptr->palette.load(argv[2]);

	winptr->needs_redraw = true;
	if (HcvisFlags.verbose) cout << "- hcsetpalette_cb\n" << flush;
	
	return TCL_OK;
}

// usage: $win hcpalette [invert | band | dim]
// will apply one of the operations {invert,band,dim} to the current palette.
static int hcpalette_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 3 || !strcmp(argv[2],"invert") && !strcmp(argv[2],"band") && !strcmp(argv[2],"dim")) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcpalette [invert | band | dim]\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcpalette: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}
	if (!strcmp(argv[2],"invert"))
		winptr->palette.invert();
	else if (!strcmp(argv[2],"band"))
		winptr->palette.band();
	else
		winptr->palette.dim();

	winptr->needs_redraw = true;
	if (HcvisFlags.verbose) cout << "- hcpalette_cb\n" << flush;
	
	return TCL_OK;
}

// usage: $win hcset3Dparams r theta phi cx cy cz fov
// theta,phi are polar coordinates of the eye point (with respect to the center point)
// fov is the field of view angle (camera angle, typically 30 degrees)
// cx,cy,cz is the center point, default (0,0,0).
// r is the RELATIVE radius of the eye point, default 1. (Relative to the model size.)
// All angles are in degrees.
static int hcset3Dparams_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 9) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcset3Dparams r theta phi cx cy cz fov\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcset3Dparams: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}

	if (a2f(argv[2]) != 0) winptr->relative_radius = a2f(argv[2]);
	if (a2f(argv[3]) != 0) winptr->theta = a2f(argv[3]);
	winptr->phi = fmod(a2f(argv[4]),360.0);
	winptr->cx = a2f(argv[5]);
	winptr->cy = a2f(argv[6]);
	winptr->cz = a2f(argv[7]);
	winptr->field_of_view = a2f(argv[8]);

	winptr->needs_redraw = true;
	winptr->can_use_display_list = true;
	winptr->XYZmode_in_3D_params = false;
	
	return TCL_OK;
}

// usage: $win hcset3Deyepoint x y z
// x y z are Cartesian coordinates. This is an alternative way of selecting the eye point.
// The method in hcset3Dparams uses the polar coordinates and relative radius.
// Call hcset3Dparams first to set the center point (cx cy cz) and camera view angle (fov)
// first, then overwrite the eye point by calling hcset3Deyepoint.
static int hcset3Deyepoint_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 5) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcset3Dparams x y z\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcset3Dparams: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}
	winptr->eyepointX = a2f(argv[2]);
	winptr->eyepointY = a2f(argv[3]);
	winptr->eyepointZ = a2f(argv[4]);
	winptr->needs_redraw = true;
	winptr->can_use_display_list = true;
	winptr->XYZmode_in_3D_params = true;

	return TCL_OK;
}

// usage: $win hcgetbox --> "xmin XMIN xmax XMAX ymin YMIN ymax YMAX zmin ZMIN zmax ZMAX"
static int hcgetbox_cb(Togl *togl, int argc, char * /*argv*/ [])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 2) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcgetbox\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcgetbox: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}

	real xmin[3]={0,0,0},xmax[3]={0,0,0};
	int dim;
	winptr->getbox(xmin,xmax,dim);
	char s[1024];
	sprintf(s,"xmin %g xmax %g ymin %g ymax %g zmin %g zmax %g",
			double(xmin[0]),double(xmax[0]),double(xmin[1]),double(xmax[1]),double(xmin[2]),double(xmax[2]));
	Tcl_SetResult(interp,s,TCL_VOLATILE);
	return TCL_OK;
}

// usage: $win hcclearobj
static int hcclearobj_cb(Togl *togl, int argc, char * /*argv*/ [])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 2) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcclearobj\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcclearobj: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}

	winptr->ClearObjects();
	// We SHOULD perhaps call ClearMinMax here, but doing so would yield ClearMinMax being
	// called almost all the time, since the object list is rebuilt each time.
	// The problem now is that if a new slice is added, the min/max are not properly updated
	// unless a variable change or other change that forces ClearMinMax to be called is also done.
//	winptr->ClearMinMax("hcclearobj");

	return TCL_OK;
}

// usage: $win hcdefobject
// will read the object definition from Object(type), Object(SlicedDim), etc.
static int hcdefobject_cb(Togl *togl, int argc, char * /*argv*/ [])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 2) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcdefobject\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcdefobject: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}

	char *const objtype = (char *)Tcl_GetVar2(interp,"Object","type",TCL_GLOBAL_ONLY);
	if (!objtype) {Tcl_SetResult(interp,const_cast<char*>("hcdefobject: Object(type) not defined"),TCL_STATIC); return TCL_ERROR;}
	if (!strcmp(objtype,"Slice")) {
		const int SlicedDim = a2i(Tcl_GetVar2(interp,"Object","SlicedDim",TCL_GLOBAL_ONLY));
		char *const Slice_Xval_string = (char*)Tcl_GetVar2(interp,"Object","Slice_Xval",TCL_GLOBAL_ONLY);
		double Slice_Xval;
		if (Slice_Xval_string && *Slice_Xval_string == 'c') {
			real xmin[3]={0,0,0},xmax[3]={0,0,0};
			int dim;
			winptr->getbox(xmin,xmax,dim);
			Slice_Xval = 0.5*(xmin[SlicedDim] + xmax[SlicedDim]);
		} else
			Slice_Xval = a2f(Slice_Xval_string);
		bool notin3 = false;
		char *const dont_string = (char*)Tcl_GetVar2(interp,"Object","DontShowInThreeD",TCL_GLOBAL_ONLY);
		if (dont_string && !strcmp(dont_string,"1")) notin3 = true;
		winptr->AddSlice(SlicedDim,Slice_Xval,notin3);
	} else if (!strcmp(objtype,"Sphere")) {
		const double r = a2f(Tcl_GetVar2(interp,"Object","Radius",TCL_GLOBAL_ONLY));
		const double x0 = a2f(Tcl_GetVar2(interp,"Object","x0",TCL_GLOBAL_ONLY));
		const double y0 = a2f(Tcl_GetVar2(interp,"Object","y0",TCL_GLOBAL_ONLY));
		const double z0 = a2f(Tcl_GetVar2(interp,"Object","z0",TCL_GLOBAL_ONLY));
		const double DeltaTheta = a2f(Tcl_GetVar2(interp,"Object","DeltaTheta",TCL_GLOBAL_ONLY));
		const double DeltaThetaGrid = a2f(Tcl_GetVar2(interp,"Object","DeltaThetaGrid",TCL_GLOBAL_ONLY));
		winptr->AddSphere(r,x0,y0,z0,DeltaTheta,DeltaThetaGrid);
	} else if (!strcmp(objtype,"FieldLineBunch")) {
		const int N = a2i(Tcl_GetVar2(interp,"Object","N",TCL_GLOBAL_ONLY));
		double r1[3],r2[3];
		r1[0] = a2f(Tcl_GetVar2(interp,"Object","r1x",TCL_GLOBAL_ONLY));
		r1[1] = a2f(Tcl_GetVar2(interp,"Object","r1y",TCL_GLOBAL_ONLY));
		r1[2] = a2f(Tcl_GetVar2(interp,"Object","r1z",TCL_GLOBAL_ONLY));
		r2[0] = a2f(Tcl_GetVar2(interp,"Object","r2x",TCL_GLOBAL_ONLY));
		r2[1] = a2f(Tcl_GetVar2(interp,"Object","r2y",TCL_GLOBAL_ONLY));
		r2[2] = a2f(Tcl_GetVar2(interp,"Object","r2z",TCL_GLOBAL_ONLY));
		TVectorField VectorField = VectorFieldRhoV;
		char *const VectorFieldString = (char*)Tcl_GetVar2(interp,"Object","VectorField",TCL_GLOBAL_ONLY);
		if (!strcmp(VectorFieldString,"rhov"))
			VectorField = VectorFieldRhoV;
		else if (!strcmp(VectorFieldString,"B"))
			VectorField = VectorFieldB;
		else if (!strcmp(VectorFieldString,"B0"))
			VectorField = VectorFieldB0;
		else if (!strcmp(VectorFieldString,"B1"))
			VectorField = VectorFieldB1;
		else if (!strcmp(VectorFieldString,"j"))
			VectorField = VectorFieldJ;
		else if (!strcmp(VectorFieldString,"S"))
			VectorField = VectorFieldPoynting;
		else if (!strcmp(VectorFieldString,"K"))
			VectorField = VectorFieldEnergyFlux;
		else
			cerr << "hcdefobject: Unknown VectorField type \"" << VectorFieldString << "\"\n";
		TTraceDirection Direction = BothTraceDirections;
		char *const DirectionString = (char*)Tcl_GetVar2(interp,"Object","Direction",TCL_GLOBAL_ONLY);
		if (!strcmp(DirectionString,"Positive"))
			Direction = PositiveTraceDirection;
		else if (!strcmp(DirectionString,"Negative"))
			Direction = NegativeTraceDirection;
		else if (!strcmp(DirectionString,"Both"))
			Direction = BothTraceDirections;
		else
			cerr << "hcdefobject: Unknown Direction type \"" << DirectionString << "\"\n";
		TFieldLineDistribution Distribution = UniformDistribution;
		char *const DistributionString = (char*)Tcl_GetVar2(interp,"Object","Distribution",TCL_GLOBAL_ONLY);
		if (!strcmp(DistributionString,"Uniform"))
			Distribution = UniformDistribution;
		else if (!strcmp(DistributionString,"InvSqrt"))
			Distribution = InvSqrtDistribution;
		else if (!strcmp(DistributionString,"Inv"))
			Distribution = InvDistribution;
		else
			cerr << "hcdefobject: Unknown Distribution type \"" << DistributionString << "\"\n";
		TLoopThresholdType LoopThresholdType = LocalSpacingLoopThreshold;
		char *const LoopThresholdTypeString = (char*)Tcl_GetVar2(interp,"Object","LoopThresholdType",TCL_GLOBAL_ONLY);
		if (!strcmp(LoopThresholdTypeString,"LocalSpacing"))
			LoopThresholdType = LocalSpacingLoopThreshold;
		else if (!strcmp(LoopThresholdTypeString,"BasegridSpacing"))
			LoopThresholdType = BasegridSpacingLoopThreshold;
		else if (!strcmp(LoopThresholdTypeString,"Absolute"))
			LoopThresholdType = AbsoluteLoopThreshold;
		else
			cerr << "hcdefobject: Unknown LoopThreshold type \"" << LoopThresholdTypeString << "\"\n";
		char *const LoopThresholdString = (char*)Tcl_GetVar2(interp,"Object","LoopThreshold",TCL_GLOBAL_ONLY);
		double LoopThreshold = atof(LoopThresholdString);
		winptr->AddFieldLineBunch(N,r1,r2,VectorField,Direction,Distribution,LoopThresholdType,LoopThreshold);
	} else {
		Tcl_SetResult(interp,const_cast<char*>("hcdefobject: Unknown object type"),TCL_STATIC);
		return TCL_ERROR;
	}
	
	return TCL_OK;
}

// usage: $win hcsetbox xmin xmax ymin ymax zmin zmax
static int hcsetbox_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 8) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcsetbox xmin xmax ymin ymax zmin zmax\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcsetbox: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}

	winptr->box_xmin[0] = a2f(argv[2]);
	winptr->box_xmax[0] = a2f(argv[3]);
	winptr->box_xmin[1] = a2f(argv[4]);
	winptr->box_xmax[1] = a2f(argv[5]);
	winptr->box_xmin[2] = a2f(argv[6]);
	winptr->box_xmax[2] = a2f(argv[7]);
	winptr->BoxDefined = true;
	winptr->needs_redraw = true;
	winptr->FullZoom();
	winptr->ClearMinMax("hcsetbox");

	return TCL_OK;
}

// usage: $win hcunsetbox
// clears any previous box setting done with hcsetbox
static int hcunsetbox_cb(Togl *togl, int argc, char * /*argv*/ [])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 2) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcunsetbox\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcunsetbox: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}

	if (winptr->BoxDefined) {
		winptr->BoxDefined = false;
		winptr->needs_redraw = true;
		winptr->FullZoom();
		winptr->ClearMinMax("hcunsetbox");
	}
	return TCL_OK;
}

// usage: $win hcintpol varname x y z
// returns the variable value (real number).
//        $win hcintpol all x y z
// returns all the variables as a name-value list
static int hcintpol_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 6) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcintpol varname x y z\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcunsetbox: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}

	char *const varname = argv[2];
	const real x = a2f(argv[3]);
	const real y = a2f(argv[4]);
	const real z = a2f(argv[5]);
	Tdimvec X(x,y,z);
	Tvariable var;
	var.pseudobackground(winptr->HasPseudoBackground());
	winptr->GridPointer()->intpol(X,HcvisFlags.intpol_order,true);
	if (!strcmp(varname,"all")) {
		int i;
		const int n = var.Nvars();
		var.select("rho",winptr->Gamma(),winptr->InvMu0(),winptr->Mass());		// Select something, in order to pass gamma,mu0,mass
		for (i=0; i<n; i++) {
			var.select(i);
			Tcl_AppendElement(interp,(char *)var.selected());
			char valuestr[80];
			sprintf(valuestr,"%g",var.get(*winptr->GridPointer(),X));
			Tcl_AppendElement(interp,valuestr);
		}
	} else {
		var.select(varname,winptr->Gamma(),winptr->InvMu0(),winptr->Mass());
		sprintf(interp->result,"%g",var.get(*winptr->GridPointer(),X));
	}
	
	return TCL_OK;
}

// usage: $win hcmodelcoords2D mousex mousey
// finds the physical model coordinates (X,Y,Z) corresponding to window coordinates mousex,mousey.
// If the grid is 2D, Z is returned as zero. If the grid is 3D, one of X,Y,Z is zero, depending on
// SlicedDim.
// The function also returns W,H which are the relevant pair of X,Y,Z. Thus the function returns five numbers.
static int hcmodelcoords2D_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 4) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcmodelcoords2D mousex mousey\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcunsetbox: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}

	const int mousex = a2i(argv[2]);
	const int mousey = a2i(argv[3]);
	double w,h,x,y,z;
	winptr->ToModelCoords(mousex,mousey,w,h, x,y,z);
	sprintf(interp->result,"%g %g %g %g %g",x,y,z,w,h);

	return TCL_OK;
}

// usage: $win hczoom wmin hmin wmax hmax
// alters the zoomed rectangle to that defined by physical coordinates (wmin,hmin), (wmax,hmax)
// To obtain these coordinates from mouse coordinates, use hcmodelcoords2D
static int hczoom_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 6) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hczoom mx1 my1 mx2 my2\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hczoom: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}

	const double wmin = a2f(argv[2]);
	const double hmin = a2f(argv[3]);
	const double wmax = a2f(argv[4]);
	const double hmax = a2f(argv[5]);
	winptr->Zoom(wmin,hmin, wmax,hmax);
	winptr->needs_redraw = true;
	
	return TCL_OK;
}

// usage: $win hcfullzoom
static int hcfullzoom_cb(Togl *togl, int argc, char * /*argv*/ [])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 2) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcfullzoom\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcfullzoom: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}
	winptr->FullZoom();
	
	return TCL_OK;
}

// usage: $win hcgetzoom --> {wmin hmin wmax hmax}
static int hcgetzoom_cb(Togl *togl, int argc, char * /*argv*/ [])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 2) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcgetzoom\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcgetzoom: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}
	double wmin,hmin,wmax,hmax;
	winptr->GetCurrentZoom(wmin,hmin,wmax,hmax);
	sprintf(interp->result,"%g %g %g %g",wmin,hmin,wmax,hmax);
	
	return TCL_OK;
}

// usage: $win hcvardescr varname --> "description string"
static int hcvardescr_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	char *const win = Togl_Ident(togl);
	if (argc != 3) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcvardescr varname\""),TCL_STATIC);
		return TCL_ERROR;
	}
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcvardescr: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}
	char *const varname = argv[2];
	if (!strcmp(varname,"x"))
		Tcl_SetResult(interp,const_cast<char*>("X-coordinate (usually meter, but see hcvis -scale option)"),TCL_STATIC);
	else if (!strcmp(varname,"y"))
		Tcl_SetResult(interp,const_cast<char*>("Y-coordinate"),TCL_STATIC);
	else if (!strcmp(varname,"z"))
		Tcl_SetResult(interp,const_cast<char*>("Z-coordinate"),TCL_STATIC);
	else {
		Tvariable var;
		var.select(varname,2,1,3);	// gamma and invmu0 and mass must be passed some valid values, but values don't matter
		Tcl_SetResult(interp,const_cast<char*>(var.description()),TCL_STATIC);
	}
	return TCL_OK;
}

// usage: $win hcexportascii fn --> ()
static int hcexportascii_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	if (argc != 3) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcexportascii filename\""),TCL_STATIC);
		return TCL_ERROR;
	}
	char *const win = Togl_Ident(togl);
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcunsetbox: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}
	char *const fn = argv[2];
	winptr->ExportASCII(fn);
	strcpy(interp->result,"");
	return TCL_OK;
}

// usage: $win hcwriteeps --> ()
static int hcwriteeps_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	if (argc != 3) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcwriteeps filename\""),TCL_STATIC);
		return TCL_ERROR;
	}
	char *const win = Togl_Ident(togl);
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcwriteeps: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}
	char *const fn = argv[2];
	Togl_DumpToEpsFile(togl, fn, 1, (void (*)(const Togl *))(&display_cb));
	cout << "Wrote EPS image in \"" << fn << "\"\n";
	strcpy(interp->result,"");
	return TCL_OK;
}

// usage: $win hcwriteppm --> ()
static int hcwriteppm_cb(Togl *togl, int argc, char *argv[])
{
	Tcl_Interp *const interp = Togl_Interp(togl);
	if (argc != 3) {
		Tcl_SetResult(interp,const_cast<char*>("wrong # args: should be \"pathName hcwriteppm filename\""),TCL_STATIC);
		return TCL_ERROR;
	}
	char *const win = Togl_Ident(togl);
	TToglWindow *const winptr = TToglWindow::Find(win);
	if (!winptr) {
		Tcl_SetResult(interp,const_cast<char*>("hcwriteppm: FindToglWindow failed"),TCL_STATIC);
		return TCL_ERROR;
	}
	char *const fn = argv[2];
	winptr->ExportPPM(fn,Togl_Width(togl),Togl_Height(togl));
	cout << "Wrote PPM image in \"" << fn << "\"\n";
	strcpy(interp->result,"");
	return TCL_OK;
}

void ParseArgs(int argc, char *argv[])
{
	int a;
	for (a=1; a<argc; a++)
		if (argv[a][0] == '-') {
			if (!strcmp(argv[a],"-3d"))
				HcvisFlags.threeD = true;
			else if (!strcmp(argv[a],"-2d"))
				HcvisFlags.threeD = false;
			else if (!strcmp(argv[a],"-v"))
				HcvisFlags.verbose = true;
		        else if (!strcmp(argv[a],"-z"))
		                HcvisFlags.intpol_order = 0;
			else if (!strcmp(argv[a],"-s"))
				HcvisFlags.verbose = false;
			else if (!strcmp(argv[a],"-scale")) {
				GridDimensionScaling = a2f(argv[++a]);
				if (GridDimensionScaling == 0) {
					cerr << "*** hcvis: Bad -scale argument, scaling reset to 1\n";
					GridDimensionScaling = 1;
				}
				GridDimensionScaling = 1.0/GridDimensionScaling;
			} else if (!strcmp(argv[a],"-init")) {
				InitFile = argv[++a];
			} else if (!strcmp(argv[a],"-help") || !strcmp(argv[a],"--help")) {
				usage();
				exit(0);
			} else
				cerr << "*** hcvis: Unknown option \"" << argv[a] << "\" ignored\n";
		} else
			break;

	// Put all remaining args to the list FileNameList
	Tcl_DStringInit(&FileNameList);
	for (; a<argc; a++)
		Tcl_DStringAppendElement(&FileNameList,argv[a]);
}

static bool isfile(const char *fn)
{
	FILE *const fp = fopen(fn,"r");
	if (fp) {
		fclose(fp);
		return true;
	}
	return false;
}

static char HCVIS_ROOT[1024];		// Initialized in main(), used in hcvis_init

#ifdef HAVE_BLT
extern "C" int Blt_Init(Tcl_Interp *);
#endif

// Called by Tk_Main() to let me initialize the modules (Togl) we need.
int hcvis_init(Tcl_Interp *interp)
{
	// Initialize Tcl, Tk, and the Togl widget module.
	if (Tcl_Init(interp) == TCL_ERROR) return TCL_ERROR;
	if (Tk_Init(interp) == TCL_ERROR) return TCL_ERROR;
	if (Togl_Init(interp) == TCL_ERROR) return TCL_ERROR;
#ifdef HAVE_BLT
	if (Blt_Init(interp) != TCL_OK) return TCL_ERROR;
#endif
	Togl_CreateFunc(create_cb);
	Togl_DisplayFunc(display_cb);
	Togl_ReshapeFunc(reshape_cb);
	Togl_DestroyFunc(destroy_cb);
	Togl_CreateCommand("hcsetparams", hcsetparams_cb);
	Togl_CreateCommand("hcsetpalette", hcsetpalette_cb);
	Togl_CreateCommand("hcpalette", hcpalette_cb);
	Togl_CreateCommand("hcset3Dparams",hcset3Dparams_cb);
	Togl_CreateCommand("hcset3Deyepoint",hcset3Deyepoint_cb);
	Togl_CreateCommand("hcgetbox",hcgetbox_cb);
	Togl_CreateCommand("hcclearobj",hcclearobj_cb);
	Togl_CreateCommand("hcdefobject",hcdefobject_cb);
	Togl_CreateCommand("hcsetbox",hcsetbox_cb);
	Togl_CreateCommand("hcunsetbox",hcunsetbox_cb);
	Togl_CreateCommand("hcintpol",hcintpol_cb);
	Togl_CreateCommand("hcmodelcoords2D",hcmodelcoords2D_cb);
	Togl_CreateCommand("hczoom",hczoom_cb);
	Togl_CreateCommand("hcfullzoom",hcfullzoom_cb);
	Togl_CreateCommand("hcgetzoom",hcgetzoom_cb);
	Togl_CreateCommand("hcvardescr",hcvardescr_cb);
	Togl_CreateCommand("hcexportascii",hcexportascii_cb);
	Togl_CreateCommand("hcwriteeps",hcwriteeps_cb);
	Togl_CreateCommand("hcwriteppm",hcwriteppm_cb);
	// Set the global variables (flags and file name list) before entering the script
	HcvisFlags.tclset(interp);
	Tcl_SetVar(interp,"FileNameList",FileNameList.string,0);
	// Set InitFile, if not null (the -init option to hcvis)
	if (InitFile) Tcl_SetVar(interp,"InitFile",InitFile,0);
	// Also set the HCVIS_ROOT variable, from above static string HCVIS_ROOT which was initialized in main() below
	Tcl_SetVar(interp,"HCVIS_ROOT",HCVIS_ROOT,0);
	return TCL_OK;
}

int main(int argc, char *argv[])
{
	ParseArgs(argc,argv);
	if (getenv("HCVIS_ROOT")) {
		strncpy(HCVIS_ROOT,getenv("HCVIS_ROOT"),1024);
		HCVIS_ROOT[1023] = '\0';
	} else {
		strncpy(HCVIS_ROOT,getenv("HOME"),1000);
		HCVIS_ROOT[999] = '\0';
		strcat(HCVIS_ROOT,"/sim/hc/lib");
	}
	char tclfile[1050];
	strcpy(tclfile,"./hcvis.tcl");
	if (!isfile(tclfile)) {
		strcpy(tclfile,HCVIS_ROOT);
		strcat(tclfile,"/hcvis.tcl");
	}
	argv[1] = tclfile;
	Tk_Main(2, argv, hcvis_init);
	return 0;
}
