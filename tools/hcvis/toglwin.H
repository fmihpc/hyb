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
#  pragma interface
#endif

extern "C" {
#include "togl.h"
}
#include <cstring>
#include "metagrid.H"
#include "palette.H"
#include "contour.H"
#include "variables.H"
#include "3Dobj.H"

#define MAX_ISO_VALUES 10

class TToglWindow {
private:
	struct TToglWindowEntry {
		TToglWindow *winptr;		// pointer to ToglWindow
		char *ident;				// the ident string of this ToglWindow
		TToglWindowEntry *next;		// pointer to next entry in the list
	};
	static TToglWindowEntry *ToglWindowList;
	double gamma,invmu0,mass;		// values read from the file by TGridCache::open
	bool pseudobackground;			// flag read from the file by TGridCache::open
	TContourSpec cs;
	Tvariable var;
	Tmetagrid *gridptr;
	bool display_list_exists;
	GLuint DisplayListID;			// meaningful only if display_list_exists
	bool dirty;						// set to true first, reset in Display()
	Togl *togl;						// not initialized if object is dirty
	GLuint helvetica_10, helvetica_12, helvetica_18;	// font bases, not initialized if object is dirty
	GLclampf background_r,background_g,background_b;
	GLclampf foreground_r,foreground_g,foreground_b;

	real wmin,wmax, hmin,hmax;
	real wmin_visible, wmax_visible, hmin_visible, hmax_visible;
	real wmin_global, wmax_global, hmin_global, hmax_global;

	int lmargin,rmargin,bmargin,tmargin;
	
	void vertex2f(GLfloat,GLfloat) const;
	void vertex3f(GLfloat,GLfloat,GLfloat) const;
	int Nslices(int d) const;		// count number of d-directed slices in objlist
	T3DSliceObject *NthSlice(int d, int n) const;		// return pointer to nth d-directed slice
	void MappedBox(const Tmetagrid& g, real xmin[3], real xmax[3]) const;
	void SetColorFromScaledQuantity(GLfloat uscaled) const;
	void SetColor(GLfloat u) const;
	void DisplaySquare(const Tdimvec& x, GLfloat dx, const Tdimvec& ex, const Tdimvec& ey,
					   real ucenter, const real u[4], const bool udens[4], const real ufacecenters[4]);
	void DisplaySquareFlat(const Tdimvec& x, GLfloat dx, const Tdimvec& ex, const Tdimvec& ey, real ucenter);
	void DisplaySquareGray(const Tdimvec& x, GLfloat dx, const Tdimvec& ex, const Tdimvec& ey);
	void VerticalColorBar(GLfloat x, GLfloat y, GLfloat w, GLfloat h) const;
	void VerticalColorBarWithLabels(GLint w0, GLint h0, int tmargin, int bmargin) const;
	bool IntersectsPlane(Tmetagrid& g, TGridIndex c, int SlicedDim, real Slice_Xval, Tdimvec& X) const;
	bool InsideBox(const Tdimvec& x) const;
	void DataSliceMinMax(Tmetagrid& g, int SlicedDim, real Slice_Xval, bool ToOldExtrema=false);
	void DisplayGridForOneCell(const Tdimvec& Xc, const GLfloat halfdx,
							   const Tdimvec& ex, const Tdimvec& ey) const;
	void DisplayGrid(Tmetagrid& g, int SlicedDim, real Slice_Xval);
	void DisplayData2D(Tmetagrid& g, int SlicedDim=2, real Slice_Xval=0.0);
	void DisplayDataSphere(Tmetagrid& g, real r, real x0, real y0, real z0, real deltatheta, real deltatheta_grid,
						   bool DrawIt=true, bool ToOldExtrema=false);
	bool EvalVectorField(TVectorField VectorField, Tmetagrid& g, const Tdimvec& r, Tdimvec& result);
	void DrawFieldLine(Tmetagrid& g, const double r[3], int sgn, TVectorField VectorField,
					   TLoopThresholdType LoopThresholdType, double LoopThreshold);
	void DrawFieldLineBunch(Tmetagrid& g, int n, const double r1[3], const double r2[3],
							TVectorField VectorField, TTraceDirection TraceDirection,
							TFieldLineDistribution Distribution,
							TLoopThresholdType LoopThresholdType, double LoopThreshold);
	void DrawIsoSurface();
	void RenderVolumetric();
public:
	char *filename;
	const char *variable;
	double datamin,datamax,one_over_datamax_minus_datamin;
	bool forced_dataminmax;
	Tpalette palette;
	// Outside visible flags:
	// ----------------------
	bool LinearInterpolation;
	bool DrawGrid;
	bool DrawGhosts;
	bool HideGhosts;
	bool DrawDead;
	bool DrawCont;
	bool DrawIsosurf;
	bool DrawVolumetric;
	bool VolumetricAntialiasing;
	bool PreserveAspect;
	bool verbose;
	bool view3D;
	bool Logarithmic;
	bool needs_redraw;
	bool can_use_display_list;
	bool XYZmode_in_3D_params;
	bool Mapping;
	bool ColorBar;
	bool Title;
	bool BoxDefined;
	bool WhiteBackground;
	bool antialiased_lines;
	bool smooth_isosurfaces;
	bool shiny_isosurfaces;
	// Outside visible 3D viewing parameters:
	// --------------------------------------
	double eyepointX, eyepointY, eyepointZ;		// used only if XYZmode_in_3D_params is true
	double relative_radius,theta,phi;			// used only if XYZmode_in_3D_params is false
	double cx,cy,cz;
	double field_of_view;
	// Outside visible other parameters:
	// ---------------------------------
	int oldheight, oldwidth;	// window geometry
	int slicedim;
	double slice_xval;			// the "current" slice_xval, copied from the slice object, needed in ToModelCoords
	T3DObject *objlist;
	double box_xmin[3], box_xmax[3];		// meaningful if BoxDefined is true
	GLfloat AlphaValue;
	double IsoValue[MAX_ISO_VALUES];
	int Nisovalues;
	int volumetric_npoints;
	double volumetric_alpha;
	// Outside visible member functions:
	// ---------------------------------
	TToglWindow();
	void SetFilename(const char *);
	static void New(const char *ident);
	static TToglWindow* Find(const char *ident);
	static void Delete(const char *ident);
	void Display(Togl *togl);
	void getbox(real xmin[3], real xmax[3], int& dim) const;
	void ClearObjects();
	void ClearMinMax(const char *reason="no reason");
	void AddSlice(int,double,bool=false);
	void AddSphere(double r,double x0=0,double y0=0,double z0=0, double DeltaTheta=5.0, double DeltaThetaGrid=10.0);
	void AddFieldLineBunch(int,const double[3],const double[3],
						   TVectorField vecfield, TTraceDirection dir,
						   TFieldLineDistribution distr, TLoopThresholdType thresholdtype, double threshold);
	double Gamma() const {return gamma;}
	double InvMu0() const {return invmu0;}
	double Mass() const {return mass;}
	Tmetagrid *GridPointer() const {return gridptr;}
	void ToModelCoords(int mousex,int mousey,double& w,double& h, double& x,double& y,double& z);
	void Zoom(real wmin1, real hmin1, real wmax1, real hmax1);
	void FullZoom();
	void GetCurrentZoom(double& wmin1, double& hmin1, double& wmax1, double& hmax1) const;
	bool HasPseudoBackground() const {return pseudobackground;}
	bool ExportASCII(const char *fn) const;
	bool ExportPPM(const char *fn, unsigned int width, unsigned int height) const;
	~TToglWindow();
};

extern real FieldLineStopSphereRadius;
