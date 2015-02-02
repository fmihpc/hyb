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

/*** X Window System headers ***/
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>  /* for XA_RGB_DEFAULT_MAP atom */
#if defined(__vms)
#include <X11/StdCmap.h>  /* for XmuLookupStandardColormap */
#else
#include <X11/Xmu/StdCmap.h>  /* for XmuLookupStandardColormap */
#endif
#include <GL/glx.h>
/*** Standard C headers ***/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#ifdef USE_LOCAL_TK_H
  #include "tk8.4a3.h"
#else
  #include <tk.h>
#endif

#if TK_MAJOR_VERSION<8
#  error Sorry Togl requires Tcl/Tk ver 8.0 or higher.
#endif

/* workaround for bug #123153 in tcl ver8.4a2 (tcl.h) */
#if defined(Tcl_InitHashTable) && defined(USE_TCL_STUBS)
#undef Tcl_InitHashTable
#define Tcl_InitHashTable (tclStubsPtr->tcl_InitHashTable)
#endif
#if (TK_MAJOR_VERSION>=8 && TK_MINOR_VERSION>=4)
#  define HAVE_TK_SETCLASSPROCS
/* pointer to Tk_SetClassProcs function in the stub table */

static void (*SetClassProcsPtr) _ANSI_ARGS_((Tk_Window, Tk_ClassProcs*,ClientData));
#endif

/*
 * Copy of TkClassProcs declarations form tkInt.h
 * (this is needed for Tcl ver =< 8.4a3)
 */

typedef int (TkBindEvalProc) _ANSI_ARGS_((ClientData clientData,
	Tcl_Interp *interp, XEvent *eventPtr, Tk_Window tkwin,
	KeySym keySym));
typedef void (TkBindFreeProc) _ANSI_ARGS_((ClientData clientData));
typedef Window (TkClassCreateProc) _ANSI_ARGS_((Tk_Window tkwin,
	Window parent, ClientData instanceData));
typedef void (TkClassGeometryProc) _ANSI_ARGS_((ClientData instanceData));
typedef void (TkClassModalProc) _ANSI_ARGS_((Tk_Window tkwin,
	XEvent *eventPtr));
typedef struct TkClassProcs {
  TkClassCreateProc *createProc;
  TkClassGeometryProc *geometryProc;
  TkClassModalProc *modalProc;
} TkClassProcs;

#include "togl.h"

/* Defaults */
#define DEFAULT_WIDTH		"400"
#define DEFAULT_HEIGHT		"400"
#define DEFAULT_IDENT		""
#define DEFAULT_FONTNAME	"fixed"
#define DEFAULT_TIME		"1"

#define MAX(a,b)	(((a)>(b))?(a):(b))

#define TCL_ERR(interp, string)			\
   do {						\
      Tcl_ResetResult(interp);			\
      Tcl_AppendResult(interp, string, NULL);	\
      return TCL_ERROR;				\
   } while (0)

/* The constant DUMMY_WINDOW is used to signal window creation 
   failure from the Togl_CreateWindow() */
#define DUMMY_WINDOW -1

#define ALL_EVENTS_MASK 	\
   (KeyPressMask |		\
    KeyReleaseMask |		\
    ButtonPressMask |		\
    ButtonReleaseMask |		\
    EnterWindowMask |		\
    LeaveWindowMask |		\
    PointerMotionMask |		\
    ExposureMask |		\
    VisibilityChangeMask |	\
    FocusChangeMask |		\
    PropertyChangeMask |	\
    ColormapChangeMask)


struct Togl
{
   struct Togl *Next;           /* next in linked list */
   GLXContext GlCtx;		/* Normal planes GLX context */
   Display *display;		/* X's token for the window's display. */
   Tk_Window  TkWin;		/* Tk window structure */
   Tcl_Interp *Interp;		/* Tcl interpreter */
   Tcl_Command widgetCmd;       /* Token for togl's widget command */
#ifndef NO_TK_CURSOR
   Tk_Cursor Cursor;		/* The widget's cursor */
#endif
   int Width, Height;		/* Dimensions of window */
   int TimerInterval;		/* Time interval for timer in milliseconds */
#if (TCL_MAJOR_VERSION * 100 + TCL_MINOR_VERSION) >= 705
   Tcl_TimerToken timerHandler; /* Token for togl's timer handler */
#else
   Tk_TimerToken timerHandler;  /* Token for togl's timer handler */
#endif
   int RgbaFlag;		/* configuration flags (ala GLX parameters) */
   int RgbaRed;
   int RgbaGreen;
   int RgbaBlue;
   int DoubleFlag;
   int DepthFlag;
   int DepthSize;
   int AccumFlag;
   int AccumRed;
   int AccumGreen;
   int AccumBlue;
   int AccumAlpha;
   int AlphaFlag;
   int AlphaSize;
   int StencilFlag;
   int StencilSize;
   int PrivateCmapFlag;
   int OverlayFlag;
   int StereoFlag;
   int AuxNumber;
   int Indirect;
   char *ShareList;             /* name (ident) of Togl to share dlists with */
   char *ShareContext;          /* name (ident) to share OpenGL context with */

   char *Ident;				/* User's identification string */
   ClientData Client_Data;		/* Pointer to user data */

   GLboolean UpdatePending;		/* Should normal planes be redrawn? */

   Togl_Callback *CreateProc;		/* Callback when widget is created */
   Togl_Callback *DisplayProc;		/* Callback when widget is rendered */
   Togl_Callback *ReshapeProc;		/* Callback when window size changes */
   Togl_Callback *DestroyProc;		/* Callback when widget is destroyed */
   Togl_Callback *TimerProc;		/* Callback when widget is idle */

   /* Overlay stuff */
   GLXContext OverlayCtx;		/* Overlay planes OpenGL context */

   Window OverlayWindow;		/* The overlay window, or 0 */
   Togl_Callback *OverlayDisplayProc;	/* Overlay redraw proc */
   GLboolean OverlayUpdatePending;	/* Should overlay be redrawn? */
   Colormap OverlayCmap;		/* colormap for overlay is created */
   int OverlayTransparentPixel;		/* transparent pixel */
   int OverlayIsMapped;

   /* for DumpToEpsFile: Added by Miguel A. de Riera Pasenau 10.01.1997 */
   XVisualInfo *VisInfo;		/* Visual info of the current */
					/* context needed for DumpToEpsFile */
   GLfloat *EpsRedMap;		/* Index2RGB Maps for Color index modes */
   GLfloat *EpsGreenMap;
   GLfloat *EpsBlueMap;
   GLint EpsMapSize;            	/* = Number of indices in our Togl */
};


/* NTNTNT need to change to handle Windows Data Types */
/*
 * Prototypes for functions local to this file
 */
static int Togl_Cmd(ClientData clientData, Tcl_Interp *interp,
                    int argc, char **argv);
static void Togl_EventProc(ClientData clientData, XEvent *eventPtr);
static Window Togl_CreateWindow(Tk_Window, Window, ClientData);
#ifdef MESA_COLOR_HACK
static int get_free_color_cells( Display *display, int screen,
                                 Colormap colormap);
static void free_default_color_cells( Display *display, Colormap colormap);
#endif
static void ToglCmdDeletedProc( ClientData );


/*
 * Setup Togl widget configuration options:
 */

static Tk_ConfigSpec configSpecs[] = {
    {TK_CONFIG_PIXELS, "-height", "height", "Height",
     DEFAULT_HEIGHT, Tk_Offset(struct Togl, Height), 0, NULL},

    {TK_CONFIG_PIXELS, "-width", "width", "Width",
     DEFAULT_WIDTH, Tk_Offset(struct Togl, Width), 0, NULL},

    {TK_CONFIG_BOOLEAN, "-rgba", "rgba", "Rgba",
     "true", Tk_Offset(struct Togl, RgbaFlag), 0, NULL},

    {TK_CONFIG_INT, "-redsize", "redsize", "RedSize",
     "1", Tk_Offset(struct Togl, RgbaRed), 0, NULL},

    {TK_CONFIG_INT, "-greensize", "greensize", "GreenSize",
     "1", Tk_Offset(struct Togl, RgbaGreen), 0, NULL},

    {TK_CONFIG_INT, "-bluesize", "bluesize", "BlueSize",
     "1", Tk_Offset(struct Togl, RgbaBlue), 0, NULL},

    {TK_CONFIG_BOOLEAN, "-double", "double", "Double",
     "false", Tk_Offset(struct Togl, DoubleFlag), 0, NULL},

    {TK_CONFIG_BOOLEAN, "-depth", "depth", "Depth",
     "false", Tk_Offset(struct Togl, DepthFlag), 0, NULL},

    {TK_CONFIG_INT, "-depthsize", "depthsize", "DepthSize",
     "1", Tk_Offset(struct Togl, DepthSize), 0, NULL},

    {TK_CONFIG_BOOLEAN, "-accum", "accum", "Accum",
     "false", Tk_Offset(struct Togl, AccumFlag), 0, NULL},

    {TK_CONFIG_INT, "-accumredsize", "accumredsize", "AccumRedSize",
     "1", Tk_Offset(struct Togl, AccumRed), 0, NULL},

    {TK_CONFIG_INT, "-accumgreensize", "accumgreensize", "AccumGreenSize",
     "1", Tk_Offset(struct Togl, AccumGreen), 0, NULL},

    {TK_CONFIG_INT, "-accumbluesize", "accumbluesize", "AccumBlueSize",
     "1", Tk_Offset(struct Togl, AccumBlue), 0, NULL},

    {TK_CONFIG_INT, "-accumalphasize", "accumalphasize", "AccumAlphaSize",
     "1", Tk_Offset(struct Togl, AccumAlpha), 0, NULL},

    {TK_CONFIG_BOOLEAN, "-alpha", "alpha", "Alpha",
     "false", Tk_Offset(struct Togl, AlphaFlag), 0, NULL},

    {TK_CONFIG_INT, "-alphasize", "alphasize", "AlphaSize",
     "1", Tk_Offset(struct Togl, AlphaSize), 0, NULL},

    {TK_CONFIG_BOOLEAN, "-stencil", "stencil", "Stencil",
     "false", Tk_Offset(struct Togl, StencilFlag), 0, NULL},

    {TK_CONFIG_INT, "-stencilsize", "stencilsize", "StencilSize",
     "1", Tk_Offset(struct Togl, StencilSize), 0, NULL},

    {TK_CONFIG_INT, "-auxbuffers", "auxbuffers", "AuxBuffers",
     "0", Tk_Offset(struct Togl, AuxNumber), 0, NULL},

    {TK_CONFIG_BOOLEAN, "-privatecmap", "privateCmap", "PrivateCmap",
     "false", Tk_Offset(struct Togl, PrivateCmapFlag), 0, NULL},

    {TK_CONFIG_BOOLEAN, "-overlay", "overlay", "Overlay",
     "false", Tk_Offset(struct Togl, OverlayFlag), 0, NULL},

    {TK_CONFIG_BOOLEAN, "-stereo", "stereo", "Stereo",
     "false", Tk_Offset(struct Togl, StereoFlag), 0, NULL},

#ifndef NO_TK_CURSOR
    { TK_CONFIG_ACTIVE_CURSOR, "-cursor", "cursor", "Cursor",
     "", Tk_Offset(struct Togl, Cursor), TK_CONFIG_NULL_OK },
#endif

    {TK_CONFIG_INT, "-time", "time", "Time",
     DEFAULT_TIME, Tk_Offset(struct Togl, TimerInterval), 0, NULL},

    {TK_CONFIG_STRING, "-sharelist", "sharelist", "ShareList",
     NULL, Tk_Offset(struct Togl, ShareList), 0, NULL},

    {TK_CONFIG_STRING, "-sharecontext", "sharecontext", "ShareContext",
     NULL, Tk_Offset(struct Togl, ShareContext), 0, NULL},

    {TK_CONFIG_STRING, "-ident", "ident", "Ident",
     DEFAULT_IDENT, Tk_Offset(struct Togl, Ident), 0, NULL},

    {TK_CONFIG_BOOLEAN, "-indirect", "indirect", "Indirect",
     "false", Tk_Offset(struct Togl, Indirect), 0, NULL},

    {TK_CONFIG_END, (char *) NULL, (char *) NULL, (char *) NULL,
     (char *) NULL, 0, 0, NULL}
};


/*
 * Default callback pointers.  When a new Togl widget is created it
 * will be assigned these initial callbacks.
 */
static Togl_Callback *DefaultCreateProc = NULL;
static Togl_Callback *DefaultDisplayProc = NULL;
static Togl_Callback *DefaultReshapeProc = NULL;
static Togl_Callback *DefaultDestroyProc = NULL;
static Togl_Callback *DefaultOverlayDisplayProc = NULL;
static Togl_Callback *DefaultTimerProc = NULL;
static ClientData DefaultClientData = NULL;
static Tcl_HashTable CommandTable;

/*
 * Head of linked list of all Togl widgets
 */
static struct Togl *ToglHead = NULL;

/*
 * Add given togl widget to linked list.
 */
static void AddToList(struct Togl *t)
{
   t->Next = ToglHead;
   ToglHead = t;
}

/*
 * Remove given togl widget from linked list.
 */
static void RemoveFromList(struct Togl *t)
{
   struct Togl *prev = NULL;
   struct Togl *pos = ToglHead;
   while (pos) {
      if (pos == t) {
         if (prev) {
            prev->Next = pos->Next;
         }
         else {
            ToglHead = pos->Next;
         }
         return;
      }
      prev = pos;
      pos = pos->Next;
   }
}

/*
 * Return pointer to togl widget given a user identifier string.
 */
static struct Togl *FindTogl(const char *ident)
{
   struct Togl *t = ToglHead;
   while (t) {
      if (strcmp(t->Ident, ident) == 0)
         return t;
      t = t->Next;
   }
   return NULL;
}

/*
 * Return an X colormap to use for OpenGL RGB-mode rendering.
 * Input:  dpy - the X display
 *         scrnum - the X screen number
 *         visinfo - the XVisualInfo as returned by glXChooseVisual()
 * Return:  an X Colormap or 0 if there's a _serious_ error.
 */
static Colormap
get_rgb_colormap( Display *dpy,
                  int scrnum,
                  const XVisualInfo *visinfo,
                  Tk_Window tkwin)
{
   Atom hp_cr_maps;
   Status status;
   int numCmaps;
   int i;
   XStandardColormap *standardCmaps;
   Window root = XRootWindow(dpy,scrnum);
   int using_mesa;

   /*
    * First check if visinfo's visual matches the default/root visual.
    */
   if (visinfo->visual==Tk_Visual(tkwin)) {
      /* use the default/root colormap */
      Colormap cmap;
      cmap = Tk_Colormap(tkwin);
#ifdef MESA_COLOR_HACK
      (void) get_free_color_cells( dpy, scrnum, cmap);
#endif
      return cmap;
   }

   /*
    * Check if we're using Mesa.
    */
   if (strstr(glXQueryServerString( dpy, scrnum, GLX_VERSION ), "Mesa")) {
      using_mesa = 1;
   }
   else {
      using_mesa = 0;
   }

   /*
    * Next, if we're using Mesa and displaying on an HP with the "Color
    * Recovery" feature and the visual is 8-bit TrueColor, search for a
    * special colormap initialized for dithering.  Mesa will know how to
    * dither using this colormap.
    */
   if (using_mesa) {
      hp_cr_maps = XInternAtom( dpy, "_HP_RGB_SMOOTH_MAP_LIST", True );
      if (hp_cr_maps
#ifdef __cplusplus
	  && visinfo->visual->c_class==TrueColor
#else
	  && visinfo->visual->class==TrueColor
#endif
	  && visinfo->depth==8) {
	 status = XGetRGBColormaps( dpy, root, &standardCmaps,
				    &numCmaps, hp_cr_maps );
	 if (status) {
	    for (i=0; i<numCmaps; i++) {
	       if (standardCmaps[i].visualid == visinfo->visual->visualid) {
                  Colormap cmap = standardCmaps[i].colormap;
                  XFree( standardCmaps );
		  return cmap;
	       }
	    }
            XFree(standardCmaps);
	 }
      }
   }

   /*
    * Next, try to find a standard X colormap.
    */
#if !HP && !SUN
#ifndef SOLARIS_BUG
   status = XmuLookupStandardColormap( dpy, visinfo->screen,
				       visinfo->visualid, visinfo->depth,
				       XA_RGB_DEFAULT_MAP,
				       /* replace */ False, /* retain */ True);
   if (status == 1) {
      status = XGetRGBColormaps( dpy, root, &standardCmaps,
				 &numCmaps, XA_RGB_DEFAULT_MAP);
      if (status == 1) {
         for (i = 0; i < numCmaps; i++) {
	    if (standardCmaps[i].visualid == visinfo->visualid) {
               Colormap cmap = standardCmaps[i].colormap;
	       XFree(standardCmaps);
	       return cmap;
	    }
	 }
         XFree(standardCmaps);
      }
   }
#endif
#endif

   /*
    * If we get here, give up and just allocate a new colormap.
    */
   return XCreateColormap( dpy, root, visinfo->visual, AllocNone );
}

/*
 * Togl_Init
 *
 *   Called upon system startup to create Togl command.
 */
int Togl_Init(Tcl_Interp *interp)
{
   int major,minor,patchLevel,releaseType;

#ifdef USE_TCL_STUBS
   if (Tcl_InitStubs(interp, "8.1", 0) == NULL) {return TCL_ERROR;}
#endif
#ifdef USE_TK_STUBS
   if (Tk_InitStubs(interp, "8.1", 0) == NULL) {return TCL_ERROR;}
#endif

   /* Skip all this on Tcl/Tk 8.0 or older.  Seems to work */
#if TCL_MAJOR_VERSION * 100 + TCL_MINOR_VERSION > 800
   Tcl_GetVersion(&major,&minor,&patchLevel,&releaseType);

#ifdef HAVE_TK_SETCLASSPROCS
   if (major >=8 && minor >= 4) {
     SetClassProcsPtr = tkStubsPtr->tk_SetClassProcs;
   } else {
     SetClassProcsPtr = NULL;
   }
#else
   if (major >= 8 && minor >= 4) {
     TCL_ERR(interp,"Sorry, this instance of Togl was not compiled to work with Tcl/Tk 8.4 or higher.");
   }
#endif

#endif

   if (Tcl_PkgProvide(interp, "Togl", TOGL_VERSION) != TCL_OK) {
      return TCL_ERROR;
   }

   Tcl_CreateCommand(interp, "togl", (Tcl_CmdProc *)Togl_Cmd,
                     (ClientData) Tk_MainWindow(interp), NULL);

   Tcl_InitHashTable(&CommandTable, TCL_STRING_KEYS);

   return TCL_OK;
}


/*
 * Register a C function to be called when an Togl widget is realized.
 */
void Togl_CreateFunc( Togl_Callback *proc )
{
   DefaultCreateProc = proc;
}


/*
 * Register a C function to be called when an Togl widget must be redrawn.
 */
void Togl_DisplayFunc( Togl_Callback *proc )
{
   DefaultDisplayProc = proc;
}


/*
 * Register a C function to be called when an Togl widget is resized.
 */
void Togl_ReshapeFunc( Togl_Callback *proc )
{
   DefaultReshapeProc = proc;
}


/*
 * Register a C function to be called when an Togl widget is destroyed.
 */
void Togl_DestroyFunc( Togl_Callback *proc )
{
   DefaultDestroyProc = proc;
}


/*
 * Register a C function to be called from TimerEventHandler.
 */
void Togl_TimerFunc( Togl_Callback *proc )
{
   DefaultTimerProc = proc;
}


/*
 * Reset default callback pointers to NULL.
 */
void Togl_ResetDefaultCallbacks( void )
{
   DefaultCreateProc = NULL;
   DefaultDisplayProc = NULL;
   DefaultReshapeProc = NULL;
   DefaultDestroyProc = NULL;
   DefaultOverlayDisplayProc = NULL;
   DefaultTimerProc = NULL;
   DefaultClientData = NULL;
}


/*
 * Chnage the create callback for a specific Togl widget.
 */
void Togl_SetCreateFunc( struct Togl *togl, Togl_Callback *proc )
{
   togl->CreateProc = proc;
}


/*
 * Change the display/redraw callback for a specific Togl widget.
 */
void Togl_SetDisplayFunc( struct Togl *togl, Togl_Callback *proc )
{
   togl->DisplayProc = proc;
}


/*
 * Change the reshape callback for a specific Togl widget.
 */
void Togl_SetReshapeFunc( struct Togl *togl, Togl_Callback *proc )
{
   togl->ReshapeProc = proc;
}


/*
 * Change the destroy callback for a specific Togl widget.
 */
void Togl_SetDestroyFunc( struct Togl *togl, Togl_Callback *proc )
{
   togl->DestroyProc = proc;
}


/*
 * Togl_Timer
 *
 * Gets called from Tk_CreateTimerHandler.
 */
static void Togl_Timer( ClientData clientData )
{
   struct Togl *togl = (struct Togl *) clientData;
   if (togl->TimerProc) {
      togl->TimerProc(togl);

      /* Re-register this callback since Tcl/Tk timers are "one-shot".
       * That is, after the timer callback is called it not normally
       * called again.  That's not the behavior we want for Togl.
       */
#if (TK_MAJOR_VERSION * 100 + TK_MINOR_VERSION) >= 401
      togl->timerHandler =
         Tcl_CreateTimerHandler( togl->TimerInterval, Togl_Timer, (ClientData)togl );
#else
      togl->timerHandler =
         Tk_CreateTimerHandler( togl->TimeInterval, Togl_Timer, (ClientData)togl );
#endif
   }
}


/*
 * Change the timer callback for a specific Togl widget.
 * Pass NULL to disable the callback.
 */
void Togl_SetTimerFunc( struct Togl *togl, Togl_Callback *proc )
{
   togl->TimerProc = proc;
   if (proc) {
#if (TK_MAJOR_VERSION * 100 + TK_MINOR_VERSION) >= 401
      togl->timerHandler =
         Tcl_CreateTimerHandler( togl->TimerInterval, Togl_Timer, (ClientData)togl );
#else
      togl->timerHandler =
         Tk_CreateTimerHandler( togl->TimeInterval, Togl_Timer, (ClientData)togl );
#endif
   }
}



/*
 * Togl_CreateCommand
 *
 *   Declares a new C sub-command of Togl callable from Tcl.
 *   Every time the sub-command is called from Tcl, the
 *   C routine will be called with all the arguments from Tcl.
 */
void Togl_CreateCommand( const char *cmd_name, Togl_CmdProc *cmd_proc)
{
   int new_item;
   Tcl_HashEntry *entry;
   entry = Tcl_CreateHashEntry(&CommandTable, cmd_name, &new_item);
   Tcl_SetHashValue(entry, cmd_proc);
}


/*
 * Togl_MakeCurrent
 *
 *   Bind the OpenGL rendering context to the specified
 *   Togl widget.
 */
void Togl_MakeCurrent( const struct Togl *togl )
{
   glXMakeCurrent( Tk_Display(togl->TkWin),
                   Tk_WindowId(togl->TkWin),
                   togl->GlCtx );
}

/*
 * Called when the widget's contents must be redrawn.  Basically, we
 * just call the user's render callback function.
 *
 * Note that the parameter type is ClientData so this function can be
 * passed to Tk_DoWhenIdle().
 */
static void Togl_Render( ClientData clientData )
{
   struct Togl *togl = (struct Togl *)clientData;

   if (togl->DisplayProc) {
      Togl_MakeCurrent(togl);
      togl->DisplayProc(togl);
   }
   togl->UpdatePending = GL_FALSE;
}


static void RenderOverlay( ClientData clientData )
{
   struct Togl *togl = (struct Togl *)clientData;

   if (togl->OverlayFlag && togl->OverlayDisplayProc) {

      glXMakeCurrent( Tk_Display(togl->TkWin),
		      togl->OverlayWindow,
		      togl->OverlayCtx );
      togl->OverlayDisplayProc(togl);
   }
   togl->OverlayUpdatePending = GL_FALSE;
}


/*
 * It's possible to change with this function or in a script some
 * options like RGBA - ColorIndex ; Z-buffer and so on
 */
int Togl_Configure(Tcl_Interp *interp, struct Togl *togl,
                   int argc, char *argv[], int flags)
{
   int oldRgbaFlag    = togl->RgbaFlag;
   int oldRgbaRed     = togl->RgbaRed;
   int oldRgbaGreen   = togl->RgbaGreen;
   int oldRgbaBlue    = togl->RgbaBlue;
   int oldDoubleFlag  = togl->DoubleFlag;
   int oldDepthFlag   = togl->DepthFlag;
   int oldDepthSize   = togl->DepthSize;
   int oldAccumFlag   = togl->AccumFlag;
   int oldAccumRed    = togl->AccumRed;
   int oldAccumGreen  = togl->AccumGreen;
   int oldAccumBlue   = togl->AccumBlue;
   int oldAccumAlpha  = togl->AccumAlpha;
   int oldAlphaFlag   = togl->AlphaFlag;
   int oldAlphaSize   = togl->AlphaSize;
   int oldStencilFlag = togl->StencilFlag;
   int oldStencilSize = togl->StencilSize;
   int oldAuxNumber   = togl->AuxNumber;

   if (Tk_ConfigureWidget(interp, togl->TkWin, configSpecs,
                          argc, (const char **)argv, (char *)togl, flags) == TCL_ERROR) {
      return(TCL_ERROR);
   }
#ifndef USE_OVERLAY
   if (togl->OverlayFlag) {
     TCL_ERR(interp,"Sorry, overlay was disabled");
   }
#endif
   

   Tk_GeometryRequest(togl->TkWin, togl->Width, togl->Height);
   /* this added per Lou Arata <arata@enya.picker.com> */
   Tk_ResizeWindow(togl->TkWin, togl->Width, togl->Height);

   if (togl->ReshapeProc &&
       togl->GlCtx
       ) {
      Togl_MakeCurrent(togl);
      togl->ReshapeProc(togl);
   }

   if (togl->RgbaFlag != oldRgbaFlag
       || togl->RgbaRed != oldRgbaRed
       || togl->RgbaGreen != oldRgbaGreen
       || togl->RgbaBlue != oldRgbaBlue
       || togl->DoubleFlag != oldDoubleFlag
       || togl->DepthFlag != oldDepthFlag
       || togl->DepthSize != oldDepthSize
       || togl->AccumFlag != oldAccumFlag
       || togl->AccumRed != oldAccumRed
       || togl->AccumGreen != oldAccumGreen
       || togl->AccumBlue != oldAccumBlue
       || togl->AccumAlpha != oldAccumAlpha
       || togl->AlphaFlag != oldAlphaFlag
       || togl->AlphaSize != oldAlphaSize
       || togl->StencilFlag != oldStencilFlag
       || togl->StencilSize != oldStencilSize
       || togl->AuxNumber != oldAuxNumber) {
#ifdef MESA_COLOR_HACK
      free_default_color_cells( Tk_Display(togl->TkWin),
                                Tk_Colormap(togl->TkWin) );
#endif
   }
   return TCL_OK;
}


int Togl_Widget(ClientData clientData, Tcl_Interp *interp,
	       int argc, char *argv[])
{
   struct Togl *togl = (struct Togl *)clientData;
   int result = TCL_OK;
   Tcl_HashEntry *entry;
   Tcl_HashSearch search;
   Togl_CmdProc *cmd_proc;

   if (argc < 2) {
      Tcl_AppendResult(interp, "wrong # args: should be \"",
                       argv[0], " ?options?\"", NULL);
      return TCL_ERROR;
   }

   Tk_Preserve((ClientData)togl);

   if (!strncmp(argv[1], "configure", MAX(1, strlen(argv[1])))) {
      if (argc == 2) {
         /* Return list of all configuration parameters */
         result = Tk_ConfigureInfo(interp, togl->TkWin, configSpecs,
                                   (char *)togl, (char *)NULL, 0);
      }
      else if (argc == 3) {
         if (strcmp(argv[2],"-extensions")==0) {
            /* Return a list of OpenGL extensions available */
            char *extensions;
            extensions = (char *) glGetString(GL_EXTENSIONS);
            Tcl_SetResult( interp, extensions, TCL_STATIC );
            result = TCL_OK;
         }
         else {
            /* Return a specific configuration parameter */
            result = Tk_ConfigureInfo(interp, togl->TkWin, configSpecs,
                                      (char *)togl, argv[2], 0);
         }
      }
      else {
         /* Execute a configuration change */
         result = Togl_Configure(interp, togl, argc-2, argv+2,
                                TK_CONFIG_ARGV_ONLY);
      }
   }
   else if (!strncmp(argv[1], "render", MAX(1, strlen(argv[1])))) {
      /* force the widget to be redrawn */
      Togl_Render((ClientData) togl);
   }
   else if (!strncmp(argv[1], "swapbuffers", MAX(1, strlen(argv[1])))) {
      /* force the widget to be redrawn */
      Togl_SwapBuffers(togl);
   }
   else if (!strncmp(argv[1], "makecurrent", MAX(1, strlen(argv[1])))) {
      /* force the widget to be redrawn */
      Togl_MakeCurrent(togl);
   }
   else if (!strncmp(argv[1], "loadbitmapfont", MAX(1, strlen(argv[1])))){
      if (argc == 3){
         GLuint fontbase;
         Tcl_Obj * fontbaseAsTclObject;
         fontbase = Togl_LoadBitmapFont(togl,argv[2]);
         if (fontbase) {
            fontbaseAsTclObject = Tcl_NewIntObj(fontbase);
            Tcl_SetObjResult(interp, fontbaseAsTclObject);
            result = TCL_OK;
         }
         else {
            Tcl_AppendResult(interp, "Could not allocate font",NULL);
            result = TCL_ERROR;
         }
      }
      else {
         Tcl_AppendResult(interp, "wrong # args",NULL);
         result=TCL_ERROR;
      }
   }
   else if (!strncmp(argv[1], "unloadbitmapfont", MAX(1, strlen(argv[1])))) {
      if (argc == 3) {
         Togl_UnloadBitmapFont(togl, atoi(argv[2]));
         result = TCL_OK;
      }
      else {
         Tcl_AppendResult(interp, "wrong # args",NULL);
         result = TCL_ERROR;
      }
   }
   else {
      /* Probably a user-defined function */
      entry = Tcl_FindHashEntry(&CommandTable, argv[1]);
      if (entry != NULL) {
         cmd_proc = (Togl_CmdProc *)Tcl_GetHashValue(entry);
         result = cmd_proc(togl, argc, argv);
      }
      else {
         Tcl_AppendResult(interp, "Togl: Unknown option: ", argv[1], "\n",
                          "Try: configure or render\n",
                          "or one of the user-defined commands:\n",
                          NULL);
         entry = Tcl_FirstHashEntry(&CommandTable, &search);
         while (entry) {
            Tcl_AppendResult(interp, "  ",
                             Tcl_GetHashKey(&CommandTable, entry),
                             "\n", NULL);
            entry = Tcl_NextHashEntry(&search);
         }
         result = TCL_ERROR;
      }
   }

   Tk_Release((ClientData)togl);
   return result;
}



/*
 * Togl_Cmd
 *
 *   Called when Togl is executed - creation of a Togl widget.
 *     * Creates a new window
 *     * Creates an 'Togl' data structure
 *     * Creates an event handler for this window
 *     * Creates a command that handles this object
 *     * Configures this Togl for the given arguments
 */
static int Togl_Cmd(ClientData clientData, Tcl_Interp *interp,
                    int argc, char **argv)
{
   char *name;
   Tk_Window main_ = (Tk_Window)clientData;
   Tk_Window tkwin;
   struct Togl *togl;

   if (argc <= 1) {
      TCL_ERR(interp, "wrong # args: should be \"pathName read filename\"");
   }

   /* Create the window. */
   name = argv[1];
   tkwin = Tk_CreateWindowFromPath(interp, main_, name, (char *) NULL);
   if (tkwin == NULL) {
      return TCL_ERROR;
   }

   Tk_SetClass(tkwin, "Togl");

   /* Create Togl data structure */
   togl = (struct Togl *)malloc(sizeof(struct Togl));
   if (!togl) {
      return TCL_ERROR;
   }

   togl->Next = NULL;
   togl->GlCtx = NULL;
   togl->OverlayCtx = NULL;
   togl->display = Tk_Display( tkwin );
   togl->TkWin = tkwin;
   togl->Interp = interp;
#ifndef NO_TK_CURSOR
   togl->Cursor = None;
#endif
   togl->Width = 0;
   togl->Height = 0;
   togl->TimerInterval = 0;
   togl->RgbaFlag = 1;
   togl->RgbaRed = 1;
   togl->RgbaGreen = 1;
   togl->RgbaBlue = 1;
   togl->DoubleFlag = 0;
   togl->DepthFlag = 0;
   togl->DepthSize = 1;
   togl->AccumFlag = 0;
   togl->AccumRed = 1;
   togl->AccumGreen = 1;
   togl->AccumBlue = 1;
   togl->AccumAlpha = 1;
   togl->AlphaFlag = 0;
   togl->AlphaSize = 1;
   togl->StencilFlag = 0;
   togl->StencilSize = 1;
   togl->OverlayFlag = 0;
   togl->StereoFlag = 0;
   togl->AuxNumber = 0;
   togl->Indirect = GL_FALSE;
   togl->UpdatePending = GL_FALSE;
   togl->OverlayUpdatePending = GL_FALSE;
   togl->CreateProc = DefaultCreateProc;
   togl->DisplayProc = DefaultDisplayProc;
   togl->ReshapeProc = DefaultReshapeProc;
   togl->DestroyProc = DefaultDestroyProc;
   togl->TimerProc = DefaultTimerProc;
   togl->OverlayDisplayProc = DefaultOverlayDisplayProc;
   togl->ShareList = NULL;
   togl->ShareContext = NULL;
   togl->Ident = NULL;
   togl->Client_Data = DefaultClientData;

   /* for EPS Output */
   togl->EpsRedMap = togl->EpsGreenMap = togl->EpsBlueMap = NULL;
   togl->EpsMapSize = 0;

   /* Create command event handler */
   togl->widgetCmd = Tcl_CreateCommand(interp, Tk_PathName(tkwin),
				       (Tcl_CmdProc *)      Togl_Widget, 
				       (ClientData)         togl,
				       (Tcl_CmdDeleteProc*) ToglCmdDeletedProc);
   /*
     Setup the Tk_ClassProcs callbacks to point at our 
     own window creation function

     We need to check at runtime if we should use the new 
     Tk_SetClassProcs() API or if we need to modify the window 
     structure directly
   */


#ifdef HAVE_TK_SETCLASSPROCS

   if (SetClassProcsPtr != NULL) {        /* use public API (Tk 8.4+) */
     Tk_ClassProcs *procsPtr;
     procsPtr = (Tk_ClassProcs*) Tcl_Alloc(sizeof(Tk_ClassProcs));
     procsPtr->size             = sizeof(Tk_ClassProcs);
     procsPtr->createProc       = Togl_CreateWindow;
     procsPtr->worldChangedProc = NULL;
     procsPtr->modalProc        = NULL;
     /*      Tk_SetClassProcs(togl->TkWin,procsPtr,(ClientData)togl); */
     (SetClassProcsPtr)(togl->TkWin,procsPtr,(ClientData)togl);
   }
   else 
#endif
     {                                  /* use private API */
       /* 
	  We need to set these fields in the Tk_FakeWin structure:
	  dummy17 = classProcsPtr
	  dummy18 = instanceData
       */
       TkClassProcs *procsPtr;
       Tk_FakeWin *winPtr = (Tk_FakeWin*)(togl->TkWin);
       
       procsPtr = (TkClassProcs*)Tcl_Alloc(sizeof(TkClassProcs));
       procsPtr->createProc     = Togl_CreateWindow;
       procsPtr->geometryProc   = NULL;
       procsPtr->modalProc      = NULL;
       winPtr->dummy17 = (char*)procsPtr;
       winPtr->dummy18 = (ClientData)togl;
     }
   
   Tk_CreateEventHandler(tkwin,
                         ExposureMask | StructureNotifyMask,
                         Togl_EventProc,
                         (ClientData)togl);

   /* Configure Togl widget */
   if (Togl_Configure(interp, togl, argc-2, argv+2, 0) == TCL_ERROR) {
      Tk_DestroyWindow(tkwin);
      goto error;
   }

   /*
    * If OpenGL window wasn't already created by Togl_Configure() we
    * create it now.  We can tell by checking if the GLX context has
    * been initialized.
    */
   if (!
       togl->GlCtx
       ) 
     {
       Tk_MakeWindowExist(togl->TkWin);
       if (Tk_WindowId(togl->TkWin)==DUMMY_WINDOW) {
	 return TCL_ERROR;
       }
       Togl_MakeCurrent(togl);
     }
   
   /* If defined, call create callback */
   if (togl->CreateProc) {
     togl->CreateProc(togl);
   }

   /* If defined, call reshape proc */
   if (togl->ReshapeProc) {
      togl->ReshapeProc(togl);
   }

   /* If defined, setup timer */
   if (togl->TimerProc){
      Tk_CreateTimerHandler( togl->TimerInterval, Togl_Timer, (ClientData)togl );
   }

   Tcl_AppendResult(interp, Tk_PathName(tkwin), NULL);

   /* Add to linked list */
   AddToList(togl);

   return TCL_OK;

error:
   Tcl_DeleteCommand(interp, "togl");
   /*free(togl);   Don't free it, if we do a crash occurs later...*/
   return TCL_ERROR;
}


#ifdef USE_OVERLAY

/*
 * Do all the setup for overlay planes
 * Return:   TCL_OK or TCL_ERROR
 */
static int SetupOverlay( struct Togl *togl )
{
#ifdef GLX_TRANSPARENT_TYPE_EXT
   static int ovAttributeList[] = {
      GLX_BUFFER_SIZE, 2,
      GLX_LEVEL, 1,
      GLX_TRANSPARENT_TYPE_EXT, GLX_TRANSPARENT_INDEX_EXT,
      None
   };
#else
   static int ovAttributeList[] = {
      GLX_BUFFER_SIZE, 2,
      GLX_LEVEL, 1,
      None
   };
#endif

   Display *dpy;
   XVisualInfo *visinfo;
   TkWindow *winPtr = (TkWindow *) togl->TkWin;

   XSetWindowAttributes swa;
   Tcl_HashEntry *hPtr;
   int new_flag;

   dpy = Tk_Display(togl->TkWin);

   visinfo = glXChooseVisual( dpy, Tk_ScreenNumber(winPtr), ovAttributeList );
   if (!visinfo){
      Tcl_AppendResult(togl->Interp,Tk_PathName(winPtr),
                       ": No suitable overlay index visual available",
                       (char *) NULL);
      togl->OverlayCtx = 0;
      togl->OverlayWindow = 0;
      togl->OverlayCmap = 0;
      return TCL_ERROR;
   }

#ifdef GLX_TRANSPARENT_INDEX_EXT
   {
      int fail = glXGetConfig(dpy, visinfo,GLX_TRANSPARENT_INDEX_VALUE_EXT,
                              &togl->OverlayTransparentPixel);
      if (fail)
         togl->OverlayTransparentPixel=0; /* maybe, maybe ... */
   }
#else
   togl->OverlayTransparentPixel=0; /* maybe, maybe ... */
#endif

   /*
   togl->OverlayCtx = glXCreateContext( dpy, visinfo, None, GL_TRUE );
   */
   /* NEW in Togl 1.5 beta 3 */
   /* share display lists with normal layer context */
   togl->OverlayCtx = glXCreateContext( dpy, visinfo,
                                        togl->GlCtx, !togl->Indirect );

   swa.colormap = XCreateColormap( dpy, XRootWindow(dpy, visinfo->screen),
                                   visinfo->visual, AllocNone );
   togl->OverlayCmap = swa.colormap;

   swa.border_pixel = 0;
   swa.event_mask = ALL_EVENTS_MASK;
   togl->OverlayWindow = XCreateWindow( dpy, Tk_WindowId(togl->TkWin), 0, 0,
                                        togl->Width, togl->Height, 0,
                                        visinfo->depth, InputOutput,
                                        visinfo->visual,
                                        CWBorderPixel|CWColormap|CWEventMask,
                                        &swa );

   hPtr = Tcl_CreateHashEntry( &winPtr->dispPtr->winTable,
                               (char *) togl->OverlayWindow, &new_flag );
   Tcl_SetHashValue( hPtr, winPtr );

/*   XMapWindow( dpy, togl->OverlayWindow );*/
   togl->OverlayIsMapped = 0;

   /* Make sure window manager installs our colormap */
   XSetWMColormapWindows( dpy, togl->OverlayWindow, &togl->OverlayWindow, 1 );

   return TCL_OK;
}
#endif /* USE_OVERLAY */


/*
 * Togl_CreateWindow
 *
 *   Window creation function, invoked as a callback from Tk_MakeWindowExist.
 *   Creates an OpenGL window for the Togl widget.
 */
static Window Togl_CreateWindow(Tk_Window tkwin,
				Window parent, 
				ClientData instanceData) {
  
  struct Togl *togl = (struct Togl*) instanceData;
  XVisualInfo *visinfo = NULL;
  Display *dpy;
  Tk_Window *winPtr = (Tk_Window *) togl->TkWin;  Colormap cmap;
  int scrnum;
  int directCtx = GL_TRUE;
  Window window;
  int attrib_list[1000];
  int attrib_count;
  int dummy;
  XSetWindowAttributes swa;
#define MAX_ATTEMPTS 12
   static int ci_depths[MAX_ATTEMPTS] = {
      8, 4, 2, 1, 12, 16, 8, 4, 2, 1, 12, 16
   };
   static int dbl_flags[MAX_ATTEMPTS] = {
      0, 0, 0, 0,  0,  0, 1, 1, 1, 1,  1,  1
   };
   dpy = Tk_Display(togl->TkWin);

   /* Make sure OpenGL's GLX extension supported */
   if (!glXQueryExtension(dpy, &dummy, &dummy)) {
     Tcl_SetResult(togl->Interp, "Togl: X server has no OpenGL GLX extension",TCL_STATIC);
     return DUMMY_WINDOW;
   }

   if (togl->ShareContext && FindTogl(togl->ShareContext)) {
      /* share OpenGL context with existing Togl widget */
      struct Togl *shareWith = FindTogl(togl->ShareContext);
      assert(shareWith);
      assert(shareWith->GlCtx);
      togl->GlCtx = shareWith->GlCtx;
      togl->VisInfo = shareWith->VisInfo;
      visinfo = togl->VisInfo;
      printf("SHARE CTX\n");
   }
   else {
      int attempt;
      /* It may take a few tries to get a visual */
      for (attempt=0; attempt<MAX_ATTEMPTS; attempt++) {
         attrib_count = 0;
         attrib_list[attrib_count++] = GLX_USE_GL;
         if (togl->RgbaFlag) {
            /* RGB[A] mode */
            attrib_list[attrib_count++] = GLX_RGBA;
            attrib_list[attrib_count++] = GLX_RED_SIZE;
            attrib_list[attrib_count++] = togl->RgbaRed;
            attrib_list[attrib_count++] = GLX_GREEN_SIZE;
            attrib_list[attrib_count++] = togl->RgbaGreen;
            attrib_list[attrib_count++] = GLX_BLUE_SIZE;
            attrib_list[attrib_count++] = togl->RgbaBlue;
            if (togl->AlphaFlag) {
               attrib_list[attrib_count++] = GLX_ALPHA_SIZE;
               attrib_list[attrib_count++] = togl->AlphaSize;
            }

            /* for EPS Output */
            if ( togl->EpsRedMap) free( ( char *)togl->EpsRedMap);
            if ( togl->EpsGreenMap) free( ( char *)togl->EpsGreenMap);
            if ( togl->EpsBlueMap) free( ( char *)togl->EpsBlueMap);
            togl->EpsRedMap = togl->EpsGreenMap = togl->EpsBlueMap = NULL;
            togl->EpsMapSize = 0;
         }
         else {
            /* Color index mode */
            int depth;
            attrib_list[attrib_count++] = GLX_BUFFER_SIZE;
            depth = ci_depths[attempt];
            attrib_list[attrib_count++] = depth;
         }
         if (togl->DepthFlag) {
            attrib_list[attrib_count++] = GLX_DEPTH_SIZE;
            attrib_list[attrib_count++] = togl->DepthSize;
         }
         if (togl->DoubleFlag || dbl_flags[attempt]) {
            attrib_list[attrib_count++] = GLX_DOUBLEBUFFER;
         }
         if (togl->StencilFlag) {
            attrib_list[attrib_count++] = GLX_STENCIL_SIZE;
            attrib_list[attrib_count++] = togl->StencilSize;
         }
         if (togl->AccumFlag) {
            attrib_list[attrib_count++] = GLX_ACCUM_RED_SIZE;
            attrib_list[attrib_count++] = togl->AccumRed;
            attrib_list[attrib_count++] = GLX_ACCUM_GREEN_SIZE;
            attrib_list[attrib_count++] = togl->AccumGreen;
            attrib_list[attrib_count++] = GLX_ACCUM_BLUE_SIZE;
            attrib_list[attrib_count++] = togl->AccumBlue;
            if (togl->AlphaFlag) {
               attrib_list[attrib_count++] = GLX_ACCUM_ALPHA_SIZE;
               attrib_list[attrib_count++] = togl->AccumAlpha;
            }
         }
         if (togl->AuxNumber != 0) {
            attrib_list[attrib_count++] = GLX_AUX_BUFFERS;
            attrib_list[attrib_count++] = togl->AuxNumber;
         }
         if (togl->Indirect) {
            directCtx = GL_FALSE;
         }

         /* stereo hack */
         /*
           if (togl->StereoFlag) {
           attrib_list[attrib_count++] = GLX_STEREO;
           }
         */
         attrib_list[attrib_count++] = None;

         visinfo = glXChooseVisual(dpy, Tk_ScreenNumber(togl->TkWin),
                                   attrib_list);
         if (visinfo) {
            /* found a GLX visual! */
            break;
         }
      }

      togl->VisInfo = visinfo;

      if (visinfo==NULL) {
	Tcl_SetResult(togl->Interp,"Togl: couldn't get visual",TCL_STATIC);
	return DUMMY_WINDOW;
      }

      /*
       * Create a new OpenGL rendering context.
       */
      if (togl->ShareList) {
         /* share display lists with existing togl widget */
         struct Togl *shareWith = FindTogl(togl->ShareList);
         GLXContext shareCtx;
         if (shareWith)
            shareCtx = shareWith->GlCtx;
         else
            shareCtx = None;
         togl->GlCtx = glXCreateContext(dpy, visinfo, shareCtx, directCtx);
      }
      else {
         /* don't share display lists */
         togl->GlCtx = glXCreateContext(dpy, visinfo, None, directCtx);
      }

      if (togl->GlCtx == NULL) {
         Tcl_SetResult(togl->Interp, "could not create rendering context",TCL_STATIC);
	 return DUMMY_WINDOW;
      }

   }

   /*
    * find a colormap
    */
   scrnum = Tk_ScreenNumber(togl->TkWin);
   if (togl->RgbaFlag) {
      /* Colormap for RGB mode */
      cmap = get_rgb_colormap( dpy, scrnum, visinfo, togl->TkWin );
   }
   else {
      /* Colormap for CI mode */
      if (togl->PrivateCmapFlag) {
         /* need read/write colormap so user can store own color entries */
         cmap = XCreateColormap(dpy, XRootWindow(dpy, visinfo->screen),
                                visinfo->visual, AllocAll);
      }
      else {
         if (visinfo->visual==DefaultVisual(dpy, scrnum)) {
            /* share default/root colormap */
            cmap = Tk_Colormap(togl->TkWin);
         }
         else {
            /* make a new read-only colormap */
            cmap = XCreateColormap(dpy, XRootWindow(dpy, visinfo->screen),
                                   visinfo->visual, AllocNone);
         }
      }
   }

   /* Make sure Tk knows to switch to the new colormap when the cursor
    * is over this window when running in color index mode.
    */
   Tk_SetWindowVisual(togl->TkWin, visinfo->visual, visinfo->depth, cmap);
   swa.colormap = cmap;
   swa.border_pixel = 0;
   swa.event_mask = ALL_EVENTS_MASK;
   window = XCreateWindow(dpy, parent,
                                  0, 0, togl->Width, togl->Height,
                                  0, visinfo->depth,
                                  InputOutput, visinfo->visual,
                                  CWBorderPixel | CWColormap | CWEventMask,
                                  &swa);
   /* Make sure window manager installs our colormap */
   XSetWMColormapWindows( dpy,window, &window, 1 );

#ifdef USE_OVERLAY
   if (togl->OverlayFlag) {
      if (SetupOverlay( togl )==TCL_ERROR) {
         fprintf(stderr,"Warning: couldn't setup overlay.\n");
         togl->OverlayFlag = 0;
      }
   }
#endif /* USE_OVERLAY */

   /* Request the X window to be displayed */
   XMapWindow(dpy, window);

   /* Check for a single/double buffering snafu */
   {
      int dbl_flag;
      if (glXGetConfig( dpy, visinfo, GLX_DOUBLEBUFFER, &dbl_flag )) {
         if (togl->DoubleFlag==0 && dbl_flag) {
            /* We requested single buffering but had to accept a */
            /* double buffered visual.  Set the GL draw buffer to */
            /* be the front buffer to simulate single buffering. */
            glDrawBuffer( GL_FRONT );
         }
      }
   }

   /* for EPS Output */
   if ( !togl->RgbaFlag) {
      int index_size;
      GLint index_bits;
      glGetIntegerv( GL_INDEX_BITS, &index_bits );
      index_size = 1 << index_bits;
      if ( togl->EpsMapSize != index_size) {
         if ( togl->EpsRedMap) free( ( char *)togl->EpsRedMap);
         if ( togl->EpsGreenMap) free( ( char *)togl->EpsGreenMap);
         if ( togl->EpsBlueMap) free( ( char *)togl->EpsBlueMap);
         togl->EpsMapSize = index_size;
         togl->EpsRedMap = ( GLfloat *)calloc( index_size, sizeof( GLfloat));
         togl->EpsGreenMap = ( GLfloat *)calloc( index_size, sizeof( GLfloat));
         togl->EpsBlueMap = ( GLfloat *)calloc( index_size, sizeof( GLfloat));
      }
   }

   return window;
}

/*
 * ToglCmdDeletedProc
 *
 *      This procedure is invoked when a widget command is deleted.  If
 *      the widget isn't already in the process of being destroyed,
 *      this command destroys it.
 *
 * Results:
 *      None.
 *
 * Side effects:
 *      The widget is destroyed.
 *
 *----------------------------------------------------------------------
 */
static void ToglCmdDeletedProc( ClientData clientData )
{
   struct Togl *togl = (struct Togl *)clientData;
   Tk_Window tkwin = togl->TkWin;

   /*
    * This procedure could be invoked either because the window was
    * destroyed and the command was then deleted (in which case tkwin
    * is NULL) or because the command was deleted, and then this procedure
    * destroys the widget.
    */

   /* NEW in togl 1.5 beta 3 */
   if (togl && tkwin) {
      Tk_DeleteEventHandler(tkwin,
                         ExposureMask | StructureNotifyMask,
                         Togl_EventProc,
                         (ClientData)togl);
   }

   /* NEW in togl 1.5 beta 3 */
   if (togl->GlCtx) {
      /* XXX this might be bad if two or more Togl widgets share a context */
      glXDestroyContext( togl->display, togl->GlCtx );
      togl->GlCtx = NULL;
   }
#ifdef USE_OVERLAY
   if (togl->OverlayCtx) {
      Tcl_HashEntry *entryPtr;
      TkWindow *winPtr = (TkWindow *) togl->TkWin;
      if (winPtr) {
         entryPtr = Tcl_FindHashEntry(&winPtr->dispPtr->winTable,
                                      (char *) togl->OverlayWindow );
         Tcl_DeleteHashEntry(entryPtr);
      }
      glXDestroyContext( togl->display, togl->OverlayCtx );
      togl->OverlayCtx = NULL;
   }
#endif /* USE_OVERLAY */

   if (tkwin != NULL) {
      togl->TkWin = NULL;
      Tk_DestroyWindow(tkwin);
   }
}


/*
 * Togl_Destroy
 *
 * Gets called when an Togl widget is destroyed.
 */
#if (TK_MAJOR_VERSION * 100 + TK_MINOR_VERSION) >= 401
static void Togl_Destroy( char *clientData )
#else
static void Togl_Destroy( ClientData clientData )
#endif
{
   struct Togl *togl = (struct Togl *)clientData;

   Tk_FreeOptions(configSpecs, (char *)togl, togl->display, 0);

#ifndef NO_TK_CURSOR
   if (togl->Cursor != None) {
      Tk_FreeCursor(togl->display, togl->Cursor);
   }
#endif

   /* remove from linked list */
   RemoveFromList(togl);
   free(togl);
}



/*
 * This gets called to handle Togl window configuration events
 */
static void Togl_EventProc(ClientData clientData, XEvent *eventPtr)
{
   struct Togl *togl = (struct Togl *)clientData;

   switch (eventPtr->type) {
      case Expose:
         if (eventPtr->xexpose.count == 0) {
            if (!togl->UpdatePending &&
                eventPtr->xexpose.window==Tk_WindowId(togl->TkWin)) {
               Togl_PostRedisplay(togl);
            }
            if (!togl->OverlayUpdatePending && togl->OverlayFlag
                && togl->OverlayIsMapped
                && eventPtr->xexpose.window==togl->OverlayWindow){
               Togl_PostOverlayRedisplay(togl);
            }
         }
         break;
      case ConfigureNotify:
         if (togl->Width != Tk_Width(togl->TkWin) ||
             togl->Height != Tk_Height(togl->TkWin)) {
            togl->Width = Tk_Width(togl->TkWin);
            togl->Height = Tk_Height(togl->TkWin);
            XResizeWindow(Tk_Display(togl->TkWin), Tk_WindowId(togl->TkWin),
                          togl->Width, togl->Height);
            if (togl->OverlayFlag) {
               XResizeWindow( Tk_Display(togl->TkWin), togl->OverlayWindow,
                              togl->Width, togl->Height );
               XRaiseWindow( Tk_Display(togl->TkWin), togl->OverlayWindow );
            }
            Togl_MakeCurrent(togl);
            if (togl->ReshapeProc) {
               togl->ReshapeProc(togl);
            }
            else {
               glViewport(0, 0, togl->Width, togl->Height);
               if (togl->OverlayFlag) {
                  Togl_UseLayer( togl,TOGL_OVERLAY );
                  glViewport( 0, 0, togl->Width, togl->Height );
                  Togl_UseLayer( togl, TOGL_NORMAL );
               }
            }
            Togl_PostRedisplay(togl);
         }
         break;
      case MapNotify:
         break;
      case DestroyNotify:
		  if (togl->DestroyProc) {
			  togl->DestroyProc(togl);
		  }
		  /* ----------------------------------------------------------------- */
	 if (togl->TkWin != NULL) {
	    togl->TkWin = NULL;
#if (TCL_MAJOR_VERSION * 100 + TCL_MINOR_VERSION) >= 800
            /* This function new in Tcl/Tk 8.0 */
            Tcl_DeleteCommandFromToken( togl->Interp, togl->widgetCmd );
#endif
	 }
	 if (togl->TimerProc != NULL) {
#if (TK_MAJOR_VERSION * 100 + TK_MINOR_VERSION) >= 401
	    Tcl_DeleteTimerHandler(togl->timerHandler);
#else
	    Tk_DeleteTimerHandler(togl->timerHandler);
#endif

	 }
	 if (togl->UpdatePending) {
#if (TCL_MAJOR_VERSION * 100 + TCL_MINOR_VERSION) >= 705
            Tcl_CancelIdleCall(Togl_Render, (ClientData) togl);
#else
            Tk_CancelIdleCall(Togl_Render, (ClientData) togl);
#endif
	 }

#if (TK_MAJOR_VERSION * 100 + TK_MINOR_VERSION) >= 401
         Tcl_EventuallyFree( (ClientData) togl, Togl_Destroy );
#else
         Tk_EventuallyFree((ClientData)togl, Togl_Destroy);
#endif

         break;
      default:
         /*nothing*/
         ;
   }
}



void Togl_PostRedisplay( struct Togl *togl )
{
   if (!togl->UpdatePending) {
      togl->UpdatePending = GL_TRUE;
      Tk_DoWhenIdle( Togl_Render, (ClientData) togl );
   }
}



void Togl_SwapBuffers( const struct Togl *togl )
{
   if (togl->DoubleFlag) {
      glXSwapBuffers( Tk_Display(togl->TkWin), Tk_WindowId(togl->TkWin) );
   }
   else {
      glFlush();
   }
}



char *Togl_Ident( const struct Togl *togl )
{
   return togl->Ident;
}


int Togl_Width( const struct Togl *togl )
{
   return togl->Width;
}


int Togl_Height( const struct Togl *togl )
{
   return togl->Height;
}


Tcl_Interp *Togl_Interp( const struct Togl *togl )
{
   return togl->Interp;
}


Tk_Window Togl_TkWin( const struct Togl *togl )
{
   return togl->TkWin;
}

/*
 * A replacement for XAllocColor.  This function should never
 * fail to allocate a color.  When XAllocColor fails, we return
 * the nearest matching color.  If we have to allocate many colors
 * this function isn't too efficient; the XQueryColors() could be
 * done just once.
 * Written by Michael Pichler, Brian Paul, Mark Kilgard
 * Input:  dpy - X display
 *         cmap - X colormap
 *         cmapSize - size of colormap
 * In/Out: color - the XColor struct
 * Output:  exact - 1=exact color match, 0=closest match
 */
static void
noFaultXAllocColor( Display *dpy, Colormap cmap, int cmapSize,
                    XColor *color, int *exact )
{
   XColor *ctable, subColor;
   int i, bestmatch;
   double mindist;       /* 3*2^16^2 exceeds long int precision.
                          */

   /* First try just using XAllocColor. */
   if (XAllocColor(dpy, cmap, color)) {
      *exact = 1;
      return;
   }

   /* Retrieve color table entries. */
   /* XXX alloca candidate. */
   ctable = (XColor *) malloc(cmapSize * sizeof(XColor));
   for (i = 0; i < cmapSize; i++) {
      ctable[i].pixel = i;
   }
   XQueryColors(dpy, cmap, ctable, cmapSize);

   /* Find best match. */
   bestmatch = -1;
   mindist = 0.0;
   for (i = 0; i < cmapSize; i++) {
      double dr = (double) color->red - (double) ctable[i].red;
      double dg = (double) color->green - (double) ctable[i].green;
      double db = (double) color->blue - (double) ctable[i].blue;
      double dist = dr * dr + dg * dg + db * db;
      if (bestmatch < 0 || dist < mindist) {
         bestmatch = i;
         mindist = dist;
      }
   }

   /* Return result. */
   subColor.red = ctable[bestmatch].red;
   subColor.green = ctable[bestmatch].green;
   subColor.blue = ctable[bestmatch].blue;
   free(ctable);
   /* Try to allocate the closest match color.  This should only
    * fail if the cell is read/write.  Otherwise, we're incrementing
    * the cell's reference count.
    */
   if (!XAllocColor(dpy, cmap, &subColor)) {
      /* do this to work around a problem reported by Frank Ortega */
      subColor.pixel = (unsigned long) bestmatch;
      subColor.red   = ctable[bestmatch].red;
      subColor.green = ctable[bestmatch].green;
      subColor.blue  = ctable[bestmatch].blue;
      subColor.flags = DoRed | DoGreen | DoBlue;
   }
   *color = subColor;
}

unsigned long Togl_AllocColor( const struct Togl *togl,
                               float red, float green, float blue )
{
   if (togl->RgbaFlag) {
      fprintf(stderr,"Error: Togl_AllocColor illegal in RGBA mode.\n");
      return 0;
   }
   if (togl->PrivateCmapFlag) {
      fprintf(stderr,"Error: Togl_FreeColor illegal with private colormap\n");
      return 0;
   }

   {
     XColor xcol;
     int exact;

     xcol.red   = (short) (red   * 65535.0);
     xcol.green = (short) (green * 65535.0);
     xcol.blue  = (short) (blue  * 65535.0);
     
     noFaultXAllocColor( Tk_Display(togl->TkWin), Tk_Colormap(togl->TkWin),
			 Tk_Visual(togl->TkWin)->map_entries, &xcol, &exact );
     /* for EPS output */
     togl->EpsRedMap[ xcol.pixel] = (float) xcol.red / 65535.0;
     togl->EpsGreenMap[ xcol.pixel] = (float) xcol.green / 65535.0;
     togl->EpsBlueMap[ xcol.pixel] = (float) xcol.blue / 65535.0;
     
     return xcol.pixel;
   }
}



void Togl_FreeColor( const struct Togl *togl, unsigned long pixel )
{
   if (togl->RgbaFlag) {
      fprintf(stderr,"Error: Togl_AllocColor illegal in RGBA mode.\n");
      return;
   }
   if (togl->PrivateCmapFlag) {
      fprintf(stderr,"Error: Togl_FreeColor illegal with private colormap\n");
      return;
   }

   XFreeColors( Tk_Display(togl->TkWin), Tk_Colormap(togl->TkWin),
                &pixel, 1, 0 );
}



void Togl_SetColor( const struct Togl *togl,
                    unsigned long index, float red, float green, float blue )
{

   if (togl->RgbaFlag) {
      fprintf(stderr,"Error: Togl_AllocColor illegal in RGBA mode.\n");
      return;
   }
   if (!togl->PrivateCmapFlag) {
      fprintf(stderr,"Error: Togl_SetColor requires a private colormap\n");
      return;
   }

   {
     XColor xcol;
     xcol.pixel = index;
     xcol.red   = (short) (red   * 65535.0);
     xcol.green = (short) (green * 65535.0);
     xcol.blue  = (short) (blue  * 65535.0);
     xcol.flags = DoRed | DoGreen | DoBlue;
     
     XStoreColor( Tk_Display(togl->TkWin), Tk_Colormap(togl->TkWin), &xcol );
     
     /* for EPS output */
     togl->EpsRedMap[ xcol.pixel] = (float) xcol.red / 65535.0;
     togl->EpsGreenMap[ xcol.pixel] = (float) xcol.green / 65535.0;
     togl->EpsBlueMap[ xcol.pixel] = (float) xcol.blue / 65535.0;
   }
}

#define MAX_FONTS 1000
static GLuint ListBase[MAX_FONTS];
static GLuint ListCount[MAX_FONTS];



/*
 * Load the named bitmap font as a sequence of bitmaps in a display list.
 * fontname may be one of the predefined fonts like TOGL_BITMAP_8_BY_13
 * or an X font name, or a Windows font name, etc.
 */
GLuint Togl_LoadBitmapFont( const struct Togl *togl, const char *fontname )
{
   static int FirstTime = 1;
   XFontStruct *fontinfo;
   int first, last, count;
   GLuint fontbase;
   const char *name;

   /* Initialize the ListBase and ListCount arrays */
   if (FirstTime) {
      int i;
      for (i=0;i<MAX_FONTS;i++) {
         ListBase[i] = ListCount[i] = 0;
      }
      FirstTime = 0;
   }

   /*
    * This method of selecting X fonts according to a TOGL_ font name
    * is a kludge.  To be fixed when I find time...
    */
   if (fontname==TOGL_BITMAP_8_BY_13) {
      name = "8x13";
   }
   else if (fontname==TOGL_BITMAP_9_BY_15) {
      name = "9x15";
   }
   else if (fontname==TOGL_BITMAP_TIMES_ROMAN_10) {
      name = "-adobe-times-medium-r-normal--10-100-75-75-p-54-iso8859-1";
   }
   else if (fontname==TOGL_BITMAP_TIMES_ROMAN_24) {
      name = "-adobe-times-medium-r-normal--24-240-75-75-p-124-iso8859-1";
   }
   else if (fontname==TOGL_BITMAP_HELVETICA_10) {
      name = "-adobe-helvetica-medium-r-normal--10-100-75-75-p-57-iso8859-1";
   }
   else if (fontname==TOGL_BITMAP_HELVETICA_12) {
      name = "-adobe-helvetica-medium-r-normal--12-120-75-75-p-67-iso8859-1";
   }
   else if (fontname==TOGL_BITMAP_HELVETICA_18) {
      name = "-adobe-helvetica-medium-r-normal--18-180-75-75-p-98-iso8859-1";
   }
   else if (!fontname) {
      name = DEFAULT_FONTNAME;
   }
   else {
      name = (const char *) fontname;
   }

   assert( name );

   fontinfo = (XFontStruct *) XLoadQueryFont( Tk_Display(togl->TkWin), name );
   if (!fontinfo) {
      return 0;
   }
   first = fontinfo->min_char_or_byte2;
   last = fontinfo->max_char_or_byte2;
   count = last-first+1;
   fontbase = glGenLists( (GLuint) (last+1) );
   if (fontbase==0) {
      return 0;
   }
   glXUseXFont( fontinfo->fid, first, count, (int) fontbase+first );
   /* Record the list base and number of display lists
    * for Togl_UnloadBitmapFont().
    */
   {
      int i;
      for (i=0;i<MAX_FONTS;i++) {
         if (ListBase[i]==0) {
            ListBase[i] = fontbase;
            ListCount[i] = last+1;
            break;
         }
      }
   }

   return fontbase;
}



/*
 * Release the display lists which were generated by Togl_LoadBitmapFont().
 */
void Togl_UnloadBitmapFont( const struct Togl *togl, GLuint fontbase )
{
   int i;
   (void) togl;
   for (i=0;i<MAX_FONTS;i++) {
      if (ListBase[i]==fontbase) {
         glDeleteLists( ListBase[i], ListCount[i] );
         ListBase[i] = ListCount[i] = 0;
         return;
      }
   }
}


/*
 * Overlay functions
 */


void Togl_UseLayer( struct Togl *togl, int layer )
{
   if (togl->OverlayWindow) {
      if (layer==TOGL_OVERLAY) {
	 glXMakeCurrent( Tk_Display(togl->TkWin),
			 togl->OverlayWindow,
			 togl->OverlayCtx );
      }
      else if (layer==TOGL_NORMAL) {
	glXMakeCurrent( Tk_Display(togl->TkWin),
			Tk_WindowId(togl->TkWin),
			togl->GlCtx );
      }
      else {
         /* error */
      }
   }
}

void Togl_ShowOverlay( struct Togl *togl )
{
   if (togl->OverlayWindow) {
      XMapWindow( Tk_Display(togl->TkWin), togl->OverlayWindow );
      XInstallColormap(Tk_Display(togl->TkWin),togl->OverlayCmap);
      togl->OverlayIsMapped = 1;
   }
}

void Togl_HideOverlay( struct Togl *togl )
{
   if (togl->OverlayWindow && togl->OverlayIsMapped) {
      XUnmapWindow( Tk_Display(togl->TkWin), togl->OverlayWindow );
      togl->OverlayIsMapped=0;
   }
}



void Togl_PostOverlayRedisplay( struct Togl *togl )
{
   if (!togl->OverlayUpdatePending
       && togl->OverlayWindow && togl->OverlayDisplayProc) {
      Tk_DoWhenIdle( RenderOverlay, (ClientData) togl );
      togl->OverlayUpdatePending = 1;
   }
}


void Togl_OverlayDisplayFunc( Togl_Callback *proc )
{
   DefaultOverlayDisplayProc = proc;
}


int Togl_ExistsOverlay( const struct Togl *togl )
{
   return togl->OverlayFlag;
}


int Togl_GetOverlayTransparentValue( const struct Togl *togl )
{
   return togl->OverlayTransparentPixel;
}


int Togl_IsMappedOverlay( const struct Togl *togl )
{
   return togl->OverlayFlag && togl->OverlayIsMapped;
}

unsigned long Togl_AllocColorOverlay( const struct Togl *togl,
                                      float red, float green, float blue )
{
   if (togl->OverlayFlag && togl->OverlayCmap) {
      XColor xcol;
      xcol.red   = (short) (red* 65535.0);
      xcol.green = (short) (green* 65535.0);
      xcol.blue  = (short) (blue* 65535.0);
      if (!XAllocColor(Tk_Display(togl->TkWin),togl->OverlayCmap,&xcol))
         return (unsigned long) -1;
      return xcol.pixel;
   }
   else {
      return (unsigned long) -1;
   }
}


void Togl_FreeColorOverlay( const struct Togl *togl, unsigned long pixel )
{

   if (togl->OverlayFlag && togl->OverlayCmap) {
      XFreeColors( Tk_Display(togl->TkWin), togl->OverlayCmap,
                   &pixel, 1, 0 );
   }
}


/*
 * User client data
 */

void Togl_ClientData( ClientData clientData )
{
   DefaultClientData = clientData;
}


ClientData Togl_GetClientData( const struct Togl *togl )
{
   return togl->Client_Data;
}


void Togl_SetClientData( struct Togl *togl, ClientData clientData )
{
   togl->Client_Data = clientData;
}

// X11-only functions

Display* Togl_Display( const struct Togl *togl)
{
   return Tk_Display(togl->TkWin);
}

Screen* Togl_Screen( const struct Togl *togl)
{
   return Tk_Screen(togl->TkWin);
}

int Togl_ScreenNumber( const struct Togl *togl)
{
   return Tk_ScreenNumber(togl->TkWin);
}

Colormap Togl_Colormap( const struct Togl *togl)
{
   return Tk_Colormap(togl->TkWin);
}



#ifdef MESA_COLOR_HACK
/*
 * Let's know how many free colors do we have
 */
#if 0
static unsigned char rojo[] = { 4, 39, 74, 110, 145, 181, 216, 251},
                     verde[] = { 4, 39, 74, 110, 145, 181, 216, 251},
		     azul[] = { 4, 39, 74, 110, 145, 181, 216, 251};

unsigned char rojo[] = { 4, 36, 72, 109, 145, 182, 218, 251},
              verde[] = { 4, 36, 72, 109, 145, 182, 218, 251},
              azul[] = { 4, 36, 72, 109, 145, 182, 218, 251};
              azul[] = { 0, 85, 170, 255};
#endif

#define RLEVELS     5
#define GLEVELS     9
#define BLEVELS     5

/* to free dithered_rgb_colormap pixels allocated by Mesa */
static unsigned long *ToglMesaUsedPixelCells = NULL;
static int ToglMesaUsedFreeCells = 0;

static int get_free_color_cells( Display *display, int screen,
                                 Colormap colormap)
{
   if ( !ToglMesaUsedPixelCells) {
      XColor xcol;
      int i;
      int colorsfailed, ncolors = XDisplayCells( display, screen);

      long r, g, b;

      ToglMesaUsedPixelCells = ( unsigned long *)calloc( ncolors, sizeof( unsigned long));

      /* Allocate X colors and initialize color_table[], red_table[], etc */
      /* de Mesa 2.1: xmesa1.c setup_dithered_(...) */
      i = colorsfailed = 0;
      for (r = 0; r < RLEVELS; r++)
         for (g = 0; g < GLEVELS; g++)
            for (b = 0; b < BLEVELS; b++) {
               int exact;
               xcol.red   = ( r*65535)/(RLEVELS-1);
               xcol.green = ( g*65535)/(GLEVELS-1);
               xcol.blue  = ( b*65535)/(BLEVELS-1);
               noFaultXAllocColor( display, colormap, ncolors,
                                   &xcol, &exact );
               ToglMesaUsedPixelCells[ i++] = xcol.pixel;
               if (!exact) {
                  colorsfailed++;
               }
            }
      ToglMesaUsedFreeCells = i;

      XFreeColors( display, colormap, ToglMesaUsedPixelCells,
                   ToglMesaUsedFreeCells, 0x00000000);
   }
   return ToglMesaUsedFreeCells;
}


static void free_default_color_cells( Display *display, Colormap colormap)
{
   if ( ToglMesaUsedPixelCells) {
      XFreeColors( display, colormap, ToglMesaUsedPixelCells,
                   ToglMesaUsedFreeCells, 0x00000000);
      free( ( char *)ToglMesaUsedPixelCells);
      ToglMesaUsedPixelCells = NULL;
      ToglMesaUsedFreeCells = 0;
   }
}
#endif


// Generate EPS file.

/* Function that creates a EPS File from a created pixmap on the current
 * context.
 * Based on the code from Copyright (c) Mark J. Kilgard, 1996.
 * Parameters: name_file, b&w / Color flag, redraw function.
 * The redraw function is needed in order to draw things into the new
 * created pixmap.
 */
static GLvoid *grabPixels(int inColor, unsigned int width, unsigned int height)
{
   GLvoid *buffer;
   GLint swapbytes, lsbfirst, rowlength;
   GLint skiprows, skippixels, alignment;
   GLenum format;
   unsigned int size;

   if (inColor) {
      format = GL_RGB;
      size = width * height * 3;
   }
   else {
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


static int generateEPS(const char *filename, int inColor,
                       unsigned int width, unsigned int height)
{
   FILE *fp;
   GLvoid *pixels;
   unsigned char *curpix;
   unsigned int components, i;
   int pos;
   unsigned char bitpixel;

   pixels = grabPixels(inColor, width, height);
   if (pixels == NULL)
      return 1;
   if (inColor)
      components = 3;     /* Red, green, blue. */
   else
      components = 1;     /* Luminance. */

   fp = fopen(filename, "w");
   if (fp == NULL) {
      return 2;
   }
   fprintf(fp, "%%!PS-Adobe-2.0 EPSF-1.2\n");
   fprintf(fp, "%%%%Creator: OpenGL pixmap render output\n");
   fprintf(fp, "%%%%BoundingBox: 0 0 %d %d\n", width, height);
   fprintf(fp, "%%%%EndComments\n");

   i = ((( width * height) + 7) / 8 ) / 40; /* # of lines, 40 bytes per line */
   fprintf(fp, "%%%%BeginPreview: %d %d %d %d\n%%", width, height, 1, i);
   pos = 0;
   curpix = ( unsigned char *)pixels;
   for ( i = 0; i < width * height * components; ) {
      bitpixel = 0;
      if ( inColor) {
         double pix = 0.0;
         pix = 0.30 * ( double)curpix[i] + 0.59 * ( double)curpix[i+1] + 0.11 * ( double)curpix[i+2];
         i += 3;
         if ( pix > 127.0) bitpixel |= 0x80;
         pix = 0.30 * ( double)curpix[i] + 0.59 * ( double)curpix[i+1] + 0.11 * ( double)curpix[i+2];
         i += 3;
         if ( pix > 127.0) bitpixel |= 0x40;
         pix = 0.30 * ( double)curpix[i] + 0.59 * ( double)curpix[i+1] + 0.11 * ( double)curpix[i+2];
         i += 3;
         if ( pix > 127.0) bitpixel |= 0x20;
         pix = 0.30 * ( double)curpix[i] + 0.59 * ( double)curpix[i+1] + 0.11 * ( double)curpix[i+2];
         i += 3;
         if ( pix > 127.0) bitpixel |= 0x10;
         pix = 0.30 * ( double)curpix[i] + 0.59 * ( double)curpix[i+1] + 0.11 * ( double)curpix[i+2];
         i += 3;
         if ( pix > 127.0) bitpixel |= 0x08;
         pix = 0.30 * ( double)curpix[i] + 0.59 * ( double)curpix[i+1] + 0.11 * ( double)curpix[i+2];
         i += 3;
         if ( pix > 127.0) bitpixel |= 0x04;
         pix = 0.30 * ( double)curpix[i] + 0.59 * ( double)curpix[i+1] + 0.11 * ( double)curpix[i+2];
         i += 3;
         if ( pix > 127.0) bitpixel |= 0x02;
         pix = 0.30 * ( double)curpix[i] + 0.59 * ( double)curpix[i+1] + 0.11 * ( double)curpix[i+2];
         i += 3;
         if ( pix > 127.0) bitpixel |= 0x01;
      }
      else {
         if ( curpix[ i++] > 0x7f) bitpixel |= 0x80;
         if ( curpix[ i++] > 0x7f) bitpixel |= 0x40;
         if ( curpix[ i++] > 0x7f) bitpixel |= 0x20;
         if ( curpix[ i++] > 0x7f) bitpixel |= 0x10;
         if ( curpix[ i++] > 0x7f) bitpixel |= 0x08;
         if ( curpix[ i++] > 0x7f) bitpixel |= 0x04;
         if ( curpix[ i++] > 0x7f) bitpixel |= 0x02;
         if ( curpix[ i++] > 0x7f) bitpixel |= 0x01;
      }
      fprintf(fp, "%02hx", bitpixel);
      if (++pos >= 40) {
         fprintf(fp, "\n%%");
         pos = 0;
      }
   }
   if (pos)
      fprintf(fp, "\n%%%%EndPreview\n");
   else
      fprintf(fp, "%%EndPreview\n");

   fprintf(fp, "gsave\n");
   fprintf(fp, "/bwproc {\n");
   fprintf(fp, "    rgbproc\n");
   fprintf(fp, "    dup length 3 idiv string 0 3 0\n");
   fprintf(fp, "    5 -1 roll {\n");
   fprintf(fp, "    add 2 1 roll 1 sub dup 0 eq\n");
   fprintf(fp, "    { pop 3 idiv 3 -1 roll dup 4 -1 roll dup\n");
   fprintf(fp, "        3 1 roll 5 -1 roll put 1 add 3 0 }\n");
   fprintf(fp, "    { 2 1 roll } ifelse\n");
   fprintf(fp, "    } forall\n");
   fprintf(fp, "    pop pop pop\n");
   fprintf(fp, "} def\n");
   fprintf(fp, "systemdict /colorimage known not {\n");
   fprintf(fp, "    /colorimage {\n");
   fprintf(fp, "        pop\n");
   fprintf(fp, "        pop\n");
   fprintf(fp, "        /rgbproc exch def\n");
   fprintf(fp, "        { bwproc } image\n");
   fprintf(fp, "    } def\n");
   fprintf(fp, "} if\n");
   fprintf(fp, "/picstr %d string def\n", width * components);
   fprintf(fp, "%d %d scale\n", width, height);
   fprintf(fp, "%d %d %d\n", width, height, 8);
   fprintf(fp, "[%d 0 0 %d 0 0]\n", width, height);
   fprintf(fp, "{currentfile picstr readhexstring pop}\n");
   fprintf(fp, "false %d\n", components);
   fprintf(fp, "colorimage\n");

   curpix = (unsigned char *) pixels;
   pos = 0;
   for (i = width * height * components; i > 0; i--) {
      fprintf(fp, "%02hx", *curpix++);
      if (++pos >= 40) {
	 fprintf(fp, "\n");
	 pos = 0;
      }
   }
   if (pos)
      fprintf(fp, "\n");

   fprintf(fp, "grestore\n");
   free(pixels);
   fclose(fp);
   return 0;
}


/* int Togl_DumpToEpsFile( const struct Togl *togl, const char *filename,
                        int inColor, void (*user_redraw)(void)) */
int Togl_DumpToEpsFile( const struct Togl *togl, const char *filename,
                        int inColor, void (*user_redraw)( const struct Togl *))
{
   int using_mesa = 0;
   Display *dpy = Tk_Display( togl->TkWin);
   int retval;
   int scrnum = Tk_ScreenNumber(togl->TkWin);
   unsigned int width = togl->Width, height = togl->Height;

   if (strstr(glXQueryServerString( dpy, scrnum, GLX_VERSION ), "Mesa"))
      using_mesa = 1;
   else
      using_mesa = 0;
   /* I don't use Pixmap do drawn into, because the code should link
    * with Mesa libraries and OpenGL libraries, and the which library
    * we use at run time should not matter, but the name of the calls
    * differs one from another:
    * MesaGl: glXCreateGLXPixmapMESA( dpy, vi, eps_pixmap, Tk_Colormap(togl->TkWin))
    * OpenGl: glXCreateGLXPixmap( dpy, vi, eps_pixmap);
    *
    * instead of this I read direct from back buffer of the screeen.
    */
   if ( !togl->RgbaFlag) {
      glPixelMapfv( GL_PIXEL_MAP_I_TO_R, togl->EpsMapSize, togl->EpsRedMap);
      glPixelMapfv( GL_PIXEL_MAP_I_TO_G, togl->EpsMapSize, togl->EpsGreenMap);
      glPixelMapfv( GL_PIXEL_MAP_I_TO_B, togl->EpsMapSize, togl->EpsBlueMap);
   }
   /*  user_redraw(); */
   user_redraw(togl);
   /* glReadBuffer( GL_FRONT); */
   /* by default it read GL_BACK in double buffer mode*/
   glFlush();
   retval = generateEPS( filename, inColor, width, height);
   return retval;
}

