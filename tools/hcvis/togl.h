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

#ifndef TOGL_H
#define TOGL_H

#define EXPORT(a,b) a b

#include <tcl.h>
#include <tk.h>
#include <GL/gl.h>
#ifdef TOGL_X11
#include <X11/Xlib.h>
#endif

#ifndef NULL
#define NULL    0
#endif


#ifdef __cplusplus
extern "C" {
#endif



#define TOGL_VERSION "1.6"
#define TOGL_MAJOR_VERSION 1
#define TOGL_MINOR_VERSION 6



/*
 * "Standard" fonts which can be specified to Togl_LoadBitmapFont()
 */
#define TOGL_BITMAP_8_BY_13		((char *) 1)
#define TOGL_BITMAP_9_BY_15		((char *) 2)
#define TOGL_BITMAP_TIMES_ROMAN_10	((char *) 3)
#define TOGL_BITMAP_TIMES_ROMAN_24	((char *) 4)
#define TOGL_BITMAP_HELVETICA_10	((char *) 5)
#define TOGL_BITMAP_HELVETICA_12	((char *) 6)
#define TOGL_BITMAP_HELVETICA_18	((char *) 7)
 

/*
 * Normal and overlay plane constants
 */
#define TOGL_NORMAL	1
#define TOGL_OVERLAY	2



struct Togl;


typedef void (Togl_Callback) (struct Togl *togl);
typedef int  (Togl_CmdProc) (struct Togl *togl, int argc, char *argv[]);
  
  EXPORT(int,Togl_Init)(Tcl_Interp *interp);

/*
 * Default/initial callback setup functions
 */

extern void Togl_CreateFunc( Togl_Callback *proc );

extern void Togl_DisplayFunc( Togl_Callback *proc );

extern void Togl_ReshapeFunc( Togl_Callback *proc );

extern void Togl_DestroyFunc( Togl_Callback *proc );

extern void Togl_TimerFunc( Togl_Callback *proc );

extern void Togl_ResetDefaultCallbacks( void );


/*
 * Change callbacks for existing widget
 */

extern void Togl_SetCreateFunc( struct Togl *togl, Togl_Callback *proc );

extern void Togl_SetDisplayFunc( struct Togl *togl, Togl_Callback *proc );

extern void Togl_SetReshapeFunc( struct Togl *togl, Togl_Callback *proc );

extern void Togl_SetDestroyFunc( struct Togl *togl, Togl_Callback *proc );

extern void Togl_SetTimerFunc( struct Togl *togl, Togl_Callback *proc );


/*
 * Miscellaneous
 */

extern int Togl_Configure( Tcl_Interp *interp, struct Togl *togl, 
                           int argc, char *argv[], int flags );

extern void Togl_MakeCurrent( const struct Togl *togl );

extern void Togl_CreateCommand( const char *cmd_name,
                                Togl_CmdProc *cmd_proc );

extern void Togl_PostRedisplay( struct Togl *togl );

extern void Togl_SwapBuffers( const struct Togl *togl );


/*
 * Query functions
 */

extern char *Togl_Ident( const struct Togl *togl );

extern int Togl_Width( const struct Togl *togl );

extern int Togl_Height( const struct Togl *togl );

extern Tcl_Interp *Togl_Interp( const struct Togl *togl );

extern Tk_Window Togl_TkWin( const struct Togl *togl );


/*
 * Color Index mode
 */

extern unsigned long Togl_AllocColor( const struct Togl *togl,
                                      float red, float green, float blue );

extern void Togl_FreeColor( const struct Togl *togl, unsigned long index );

extern void Togl_SetColor( const struct Togl *togl, unsigned long index,
                           float red, float green, float blue );


/*
 * Bitmap fonts
 */

extern GLuint Togl_LoadBitmapFont( const struct Togl *togl,
                                   const char *fontname );

extern void Togl_UnloadBitmapFont( const struct Togl *togl, GLuint fontbase );


/*
 * Overlay functions
 */

extern void Togl_UseLayer( struct Togl *togl, int layer );

extern void Togl_ShowOverlay( struct Togl *togl );

extern void Togl_HideOverlay( struct Togl *togl );

extern void Togl_PostOverlayRedisplay( struct Togl *togl );

extern void Togl_OverlayDisplayFunc( Togl_Callback *proc );

extern int Togl_ExistsOverlay( const struct Togl *togl );

extern int Togl_GetOverlayTransparentValue( const struct Togl *togl );

extern int Togl_IsMappedOverlay( const struct Togl *togl );

extern unsigned long Togl_AllocColorOverlay( const struct Togl *togl,
                                             float red, float green, 
                                             float blue );

extern void Togl_FreeColorOverlay( const struct Togl *togl, 
                                   unsigned long index );

/*
 * User client data
 */

extern void Togl_ClientData( ClientData clientData );

extern ClientData Togl_GetClientData( const struct Togl *togl );

extern void Togl_SetClientData( struct Togl *togl, ClientData clientData );

#ifdef TOGL_X11
extern Display *Togl_Display( const struct Togl *togl );
extern Screen *Togl_Screen( const struct Togl *togl );
extern int Togl_ScreenNumber( const struct Togl *togl );
extern Colormap Togl_Colormap( const struct Togl *togl );
#endif


extern int Togl_DumpToEpsFile( const struct Togl *togl,
                               const char *filename,
                               int inColor,
                               void (*user_redraw)(const struct Togl *) );

#ifdef __cplusplus
}
#endif


#endif
