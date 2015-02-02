# This file is part of the HYB simulation platform.
#
# Copyright 2014- Finnish Meteorological Institute
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Global variables
# 
# Each open data window is labelled by a string identifier which we call win.
# Associated with each window is a set of Flags, and the FileNameIndex.
# There is also a set of Params for the 3D view point parameters for each win.
# All available filenames are stored in FileNameList.
# This list is initialized in hcvis.C before entering this script.
# The palette file associated with win is Palette(win).
# Slave relationships are stored in matrix Slave(slavewin,masterwin).

# The variables FileNameList and HCVIS_ROOT, and possibly others, come from hcvis.C

if {![file isdirectory $HCVIS_ROOT]} {set HCVIS_ROOT "."}
if {[file exists ./hcvisfilesel.tcl]} {set HCVIS_ROOT "."}
if {$verbose} {puts "HCVIS_ROOT = $HCVIS_ROOT"}

source $HCVIS_ROOT/hcvisfilesel.tcl
source $HCVIS_ROOT/hcvisobj.tcl

# Constants
set FlagList {LinearInterpolation DrawGhosts HideGhosts DrawGrid DrawDead DrawCont \
		      DrawIsosurf SmoothIsosurfaces ShinyIsosurfaces \
			  PreserveAspect view3D Logarithmic Mapping ColorBar Title BoxDefined FullWindow WhiteBackground \
			  AntialiasedLines DrawVolumetric VolumetricAntialiasing}
set ParamList {relative_radius theta phi cx cy cz fov var alpha FieldLineStopSphereRadius slicedim zoomstack obj \
			   IsoValue0 IsoValue1 IsoValue2 IsoValue3 IsoValue4 \
			   VolumetricNpoints VolumetricAlpha}
set PaletteDir $HCVIS_ROOT/palettes
set statefile .hcvisrc
set ZoomInFactor 1.2

# Default flags
set Flags(default,LinearInterpolation) 0
set Flags(default,DrawGhosts) 0
set Flags(default,HideGhosts) 0
set Flags(default,DrawGrid) 1
set Flags(default,DrawDead) 0
set Flags(default,DrawCont) 0
set Flags(default,DrawIsosurf) 0
set Flags(default,SmoothIsosurfaces) 0
set Flags(default,ShinyIsosurfaces) 0
set Flags(default,view3D) 0
set Flags(default,Logarithmic) 0
set Flags(default,Mapping) 0
set Flags(default,ColorBar) 1
set Flags(default,Title) 1
set Flags(default,PreserveAspect) 1
set Flags(default,BoxDefined) 0
set Flags(default,FullWindow) 0
set Flags(default,WhiteBackground) 0
set Flags(default,AntialiasedLines) 0
set Flags(default,DrawVolumetric) 0
set Flags(default,VolumetricAntialiasing) 0

# Default 3D viewpoint parameters
set Params(default,relative_radius) 1
set Params(default,theta) 50
set Params(default,phi) 50
set Params(default,cx) 0
set Params(default,cy) 0
set Params(default,cz) 0
set Params(default,fov) 30
# Default other parameters
set Params(default,var) rho
set Params(default,alpha) 1.0
set Params(default,FieldLineStopSphereRadius) 0.0
set Params(default,slicedim) 1
set Params(default,zoomstack) {}
set Params(default,obj) {\
	{Slice {X-Slice 1} active {SlicedDim 0 Slice_Xval cx}}\
	{Slice {Y-Slice 1} active {SlicedDim 1 Slice_Xval cy}}\
	{Slice {Z-Slice 1} active {SlicedDim 2 Slice_Xval cz}}\
}
set Params(default,IsoValue0) 2.0
set Params(default,IsoValue1) {}
set Params(default,IsoValue2) {}
set Params(default,IsoValue3) {}
set Params(default,IsoValue4) {}
set Params(default,VolumetricNpoints) 400000
set Params(default,VolumetricAlpha) 0.2

set FileNameIndex(default) 0
#set Palette(default) Grayscale
set Palette(default) $PaletteDir/rainbow.raw

# Some label texts for 3D view point selector dialog box
set Edit3DViewpoint(rt) {relative r}
set Edit3DViewpoint(thetat) {theta}
set Edit3DViewpoint(phit) {phi}
set Edit3DViewpoint(cxt) {cx}
set Edit3DViewpoint(cyt) {cy}
set Edit3DViewpoint(czt) {cz}
set Edit3DViewpoint(fovt) {view angle}
set Edit3DViewpoint(g1) {Eye point :}
set Edit3DViewpoint(g2) {Center point :}
set Edit3DViewpoint(g3) {Camera :}

if {[info exists InitFile]} {
	source $HCVIS_ROOT/$InitFile
}

# Bookkeeping information
set winID 0
set NWinsOpen 0
set NFileNames 0
set WindowList {}

# State saving and restoring
proc SaveState {fn} {
	global Flags Params FileNameList FileNameIndex Palette
	global winID NWinsOpen NFileNames WindowList
	global WinGeom Slave
	set fp [open $fn w]
	puts $fp [array get Flags]
	puts $fp [array get Params]
	puts $fp $FileNameList
	puts $fp [array get FileNameIndex]
	puts $fp [array get Palette]
	puts $fp $winID
	puts $fp $NWinsOpen
	puts $fp $NFileNames
	puts $fp $WindowList
	array set WinGeom {}
	foreach i $WindowList {
		if [winfo exists $i] {
			set WinGeom($i) [winfo geometry $i]
		} else {
			set WinGeom($i) {}
		}
	}
	puts $fp [array get WinGeom]
	puts $fp [array get Slave]
	puts $fp [array get DefineBox]
	close $fp
}

proc LoadState {fn} {
	global Flags Params FileNameList FileNameIndex Palette
	global winID NWinsOpen NFileNames WindowList
	global WinGeom Slave
	set fp [open $fn r]
	array set Flags {}
	array set Flags [gets $fp]
	array set Params {}
	array set Params [gets $fp]
#	array set FileNameList {}
	set FileNameList [gets $fp]
	array set FileNameIndex {}
	array set FileNameIndex [gets $fp]
	array set Palette {}
	array set Palette [gets $fp]
	set winID [gets $fp]
	set NWinsOpen [gets $fp]
	set NFileNames [gets $fp]
	# Close all open windows now. Assume that the file reads successfully.
	foreach i $WindowList {destroy $i}
	set WindowList [gets $fp]
	# Reset all zoomstacks. Restoring the zoom state automatically would be messy.
	foreach i $WindowList {set Params($i,zoomstack) {}}
	array set WinGeom {}
	array set WinGeom [gets $fp]
	array set Slave {}
	array set Slave [gets $fp]
	array set DefineBox {}
	array set DefineBox [gets $fp]
	close $fp
	foreach i $WindowList {
		if {[string length $WinGeom($i)] > 0} {
			if ![winfo exists $i] {OpenNewWindow $i 1}
			wm geometry $i $WinGeom($i)
			Edit3DApply $i
			# Update is unncessary here since Edit3DApply calls it
		}
	}
}

proc YesNoQuery {str} {
	global YesNoQuery
	set YesNoQuery(result) 1
	set t [toplevel .yesno]
	grab $t
	label $t.l -text $str
	frame $t.buttons
	set ok [button $t.buttons.ok -text OK -underline 0 -command "set YesNoQuery(result) 1; destroy $t"]
	set can [button $t.buttons.cancel -text Cancel -underline 0 -command "set YesNoQuery(result) 0; destroy $t"]
	pack $t.l -side top
	pack $t.buttons -side bottom
	pack $ok -side left
	pack $can -side right
	SetWindowCloseToMouse $t
	tkwait window $t
	grab release $t
	return $YesNoQuery(result)
}

proc QuerySaveState {} {
	global statefile
	set fn [fileselect "Save hcvis state to" $statefile 0 1 0]
	if {[file exists $fn]} {
		set saveok [YesNoQuery "OK to overwrite $fn ?"]
	} else {
		set saveok 1
	}
	if {$saveok && [string length $fn] > 0} {SaveState $fn}
}

proc QueryLoadState {} {
	global statefile
	set fn [fileselect "Restore hcvis state from" $statefile 1 1 0]
	if {[string length $fn] > 0} {LoadState $fn}
}

# Member x L returns 1 if x is member of list L and 0 otherwise
proc Member {x L} {
	if {[lsearch -exact $L $x] >= 0} {return 1} else {return 0}
}

# MaxLength L returns the maximum (string) length of the elements of list L
proc MaxLength {L} {
	set result 0
	foreach f $L {
		set len [string length $f]
		if {$len > $result} {set result $len}
	}
	return $result
}

# TakeEverySecond L takes the first, third, etc., elements of list L and returns a new list
proc TakeEverySecond {L} {
	set result {}
	set n [llength $L]
	for {set i 0} {$i < $n} {incr i 2} {
		lappend result [lrange $L $i $i]
	}
	return $result
}

# ldelete list value deletes one item from a list and returns it, from Welch' book
proc ldelete {list value} {
	set ix [lsearch -exact $list $value]
	if {$ix >= 0} {
		return [lreplace $list $ix $ix]
	} else {
		return $list
	}
}

# IsSlave win tests if a Slave(win,*) definition exists for win
proc IsSlave {win} {
	global Slave
	if {[string length [array names Slave $win,*]] > 0} {
		return 1
	} else {
		return 0
	}
}

# ----------------------------------------------------------------------------------------

proc FullWindow {win} {
	global FullWindow Flags
	if $Flags($win,FullWindow) {
		# Window is now full, restore the old geometry
		set newgeom $FullWindow($win,geom)
		set Flags($win,FullWindow) 0
		$win.menubar.view.menu entryconfigure "Restore*" -label "Full window        Ctrl+F"
	} else {
		# Make window full
		set FullWindow($win,geom) [winfo geometry $win]
		set w [winfo screenwidth $win]
		set h [winfo screenheight $win]
		set newgeom "${w}x${h}+0+0"
		set Flags($win,FullWindow) 1
		$win.menubar.view.menu entryconfigure "Full*" -label "Restore window    Ctrl+F"
	}
	wm geometry $win $newgeom
}

proc CloseWindow {win} {
	global NWinsOpen WindowList
	destroy $win
	set WindowList [ldelete $WindowList $win]
	incr NWinsOpen -1
	if {$NWinsOpen <= 0} exit
}

proc UpdateObjects {win} {
#	Copies Params($win,obj) to C++ using hcclearobj, hcdefobject
	global Params Object
	$win.win hcclearobj
	foreach obj $Params($win,obj) {
		if {[string compare [lindex $obj 2] active]==0} {
			ToObject $obj
			set Object(DontShowInThreeD) 0
			$win.win hcdefobject
		} else {
			ToObject $obj
			if {[string compare $Object(type) Slice]==0} {
				set Object(DontShowInThreeD) 1
				$win.win hcdefobject
			}
		}
	}
}

proc Tail {fn} {
#	file tail $fn
	concat $fn
}

proc Update {win} {
	global FileNameIndex NFileNames FileNameList Flags Params ZoomInFactor verbose
	if {$verbose} {puts "- Calling Update"}
	set Params($win,VolumetricNpoints) [expr round(1e6 * $Params($win,VolumetricNpointsMillions))]
	# If contours are requested, linear interpolation must be turned on also, otherwise the contours won't be visible
	if $Flags($win,DrawCont) {set Flags($win,LinearInterpolation) 1}
	set i $FileNameIndex($win)
	if {![file isfile [lindex $FileNameList $i]]} {
		return
	}
	# Copy parameters of win from Tcl to C++
	$win.win hcsetparams
	# Update objects
	UpdateObjects $win
	# Post redisplay
	$win.win render
	# Disable the "Next file" menu item if this is the last file
	set fm $win.menubar.file.menu
	if {$i >= $NFileNames - 1} {set state disabled} {set state normal}
	$fm entryconfigure "Next*" -state $state
	# Disable the "Prev file" menu item if this is the first file
	if {$i <= 0} {set state disabled} {set state normal}
	$fm entryconfigure "Prev*" -state $state
	# Set window title
	if {[IsSlave $win]} {
		wm title $win "[Tail [lindex $FileNameList $i]] (Slave)"
		# A slave cannot make a new slave, disable the menu and the keyboard equiv
		$fm entryconfigure "New slave*" -state disabled
		bind $win <Control-d> {}
	} else {
		wm title $win [Tail [lindex $FileNameList $i]]
	}
	# Set or unset some key bindings, depending whether we are in 2D or 3D
	set om $win.menubar.options.menu
	set wm $win.menubar.view.menu
	set fm $win.menubar.file.menu
	if {$Flags($win,view3D)} {
		bind $win <Key-Left> [list Change3DViewpoint $win 1 0 -10]
		bind $win <Key-Right> [list Change3DViewpoint $win 1 0 10]
		bind $win <Key-Up> [list Change3DViewpoint $win 1 -10 0]
		bind $win <Key-Down> [list Change3DViewpoint $win 1 10 0]
		bind $win <Control-Key-Left> [list Change3DViewpoint $win 1 0 -2]
		bind $win <Control-Key-Right> [list Change3DViewpoint $win 1 0 2]
		bind $win <Control-Key-Up> [list Change3DViewpoint $win 1 -2 0]
		bind $win <Control-Key-Down> [list Change3DViewpoint $win 1 2 0]
		bind $win <Shift-Key-Left> [list Change3DViewpoint $win 1 0 -90]
		bind $win <Shift-Key-Right> [list Change3DViewpoint $win 1 0 90]
		bind $win <Shift-Key-Up> [list Change3DViewpoint $win 1 -45 0]
		bind $win <Shift-Key-Down> [list Change3DViewpoint $win 1 45 0]
		bind $win <Key-less> [list Change3DViewpoint $win [expr 1/$ZoomInFactor] 0 0]
		bind $win <Key-greater> [list Change3DViewpoint $win $ZoomInFactor 0 0]
		bind $win <Key-o> [list Origin3DViewpoint $win]
		bind $win <Key-c> [list GridCenter3DViewpoint $win]
		bind $win <Key-X> [list Set3DAngles $win 90 0]
		bind $win <Key-Y> [list Set3DAngles $win 90 90]
		bind $win <Key-Z> [list Set3DAngles $win 179.99 180]
		bind $win <Control-Key-3> [list Edit3DViewpoint $win]
		bind $win <Control-Key-t> [list SetTransparency $win]
		# Enable the 3D menu items
		$wm entryconfigure "Zoom in*" -state normal
		$wm entryconfigure "Zoom out*" -state normal
		$wm entryconfigure "Origin viewpoint*" -state normal
		$wm entryconfigure "Center viewpoint*" -state normal
		$wm entryconfigure "Negative X view*" -state normal
		$wm entryconfigure "Negative Y view*" -state normal
		$wm entryconfigure "Negative Z view*" -state normal
		$wm entryconfigure "3D viewpoint*" -state normal
		$om entryconfigure "Transparency*" -state normal
		# Disable the Preserve aspect, Variable values, Full zoom, Previous zoom, and Export 2D menu items
		$om entryconfigure "Preserve aspect*" -state disabled
		$fm entryconfigure "Variable values*" -state disabled
		$wm entryconfigure "Full zoom*" -state disabled
		$wm entryconfigure "Previous zoom*" -state disabled
		$fm entryconfigure "Export 2D*" -state disabled
		bind $win <Control-Key-v> {}
		bind $win <Key-f> {}
		bind $win <Key-p> {}
		# Enable the Run orbit file menu item
		$fm entryconfigure "Run orbit file*" -state normal
		# Unbind mouse clicks
		bind $win.win <ButtonPress-1> {}
		bind $win.win <ButtonRelease-1> {}
	} else {
		bind $win <Key-Left> {}
		bind $win <Key-Right> {}
		bind $win <Key-Up> {}
		bind $win <Key-Down> {}
		bind $win <Control-Key-Left> {}
		bind $win <Control-Key-Right> {}
		bind $win <Control-Key-Up> {}
		bind $win <Control-Key-Down> {}
		bind $win <Shift-Key-Left> {}
		bind $win <Shift-Key-Right> {}
		bind $win <Shift-Key-Up> {}
		bind $win <Shift-Key-Down> {}
		bind $win <Key-less> {}
		bind $win <Key-greater> {}
		bind $win <Key-o> {}
		bind $win <Key-c> {}
		bind $win <Key-X> {}
		bind $win <Key-Y> {}
		bind $win <Key-Z> {}
		bind $win <Control-Key-3> {}
		bind $win <Control-Key-t> {}
		# Disable the 3D menu items
		$wm entryconfigure "Zoom in*" -state disabled
		$wm entryconfigure "Zoom out*" -state disabled
		$wm entryconfigure "Origin viewpoint*" -state disabled
		$wm entryconfigure "Center viewpoint*" -state disabled
		$wm entryconfigure "Negative X view*" -state disabled
		$wm entryconfigure "Negative Y view*" -state disabled
		$wm entryconfigure "Negative Z view*" -state disabled
		$wm entryconfigure "3D viewpoint*" -state disabled
		$om entryconfigure "Transparency*" -state disabled
		# Enabled the Preserve aspect, Variable values, Full zoom, Previous zoom and Export 2D menu items
		$om entryconfigure "Preserve aspect*" -state normal
		$fm entryconfigure "Variable values*" -state normal
		$wm entryconfigure "Full zoom*" -state normal
		$wm entryconfigure "Previous zoom*" -state normal
		$fm entryconfigure "Export 2D*" -state normal
		bind $win <Key-f> [list FullZoom $win]
		bind $win <Key-p> [list PreviousZoom $win]
		bind $win <Control-Key-v> [list VariableWindow $win]
		# Disable the Run orbit file menu item
		$fm entryconfigure "Run orbit file*" -state disabled
		# Bind mouse clicks
		bind $win.win <ButtonPress-1> [list MouseDown $win %x %y]
		bind $win.win <ButtonRelease-1> [list MouseUp $win %x %y]
	}
	# If a 1D view window exists for this window, update it also
	if {[winfo exists $win.oned]} {
		OneDimensionalViewUpdate $win
	}
}

proc UpdatePossibleSlaves {win} {
#	Call this after modifying any window's FileNameIndex
#	For slave windows it does nothing, but for master windows,
#	it sets the slave windows' FileNameIndex equal to the master one's,
#	and updates them.
	global Slave FileNameIndex
	set slaves [TakeEverySecond [array get Slave *,$win]]
	foreach slavemaster $slaves {
		# slavemaster is a comma-separated pair slave,master
		set commapos [string first "," $slavemaster]
		set slave [string range $slavemaster 0 [expr $commapos - 1]]
		if {[winfo exists $slave]} {
			set FileNameIndex($slave) $FileNameIndex($win)
			Update $slave
		}
	}
}

proc NextFile {win} {
	global NFileNames FileNameList FileNameIndex
	if {$FileNameIndex($win) < $NFileNames - 1} {
		incr FileNameIndex($win)
		Update $win
		UpdatePossibleSlaves $win
	}
}

proc PrevFile {win} {
	global NFileNames FileNameList FileNameIndex
	if {$FileNameIndex($win) > 0} {
		incr FileNameIndex($win) -1
		Update $win
		UpdatePossibleSlaves $win
	}
}

proc LoadFile {{update 0} {win {}}} {
	global FileNameList NFileNames FileNameIndex
	set newfn [fileselect "Load an HC file"]
	if {[string length $newfn]!=0} {
		set pwd [pwd]			
		foreach f $newfn {
			if {[string compare [string range $f 0 [expr [string length $pwd] - 1]] $pwd] == 0} {
				set ff [string range $f [expr [string length $pwd] + 1] end]
				# +1 to skip the "/"
			} else {
				set ff $f
			}
			if {![Member $ff $FileNameList]} {
				lappend FileNameList $ff
				incr NFileNames
			} else {
				puts "Ignored file $ff -- already loaded"
			}
		}
		if {$update} {
			set FileNameIndex($win) [expr $NFileNames-1]
			Update $win
			UpdatePossibleSlaves $win
		}
	}
	return $newfn
}

# ------------------------------------------------------------------
# --- GotoFile
# ------------------------------------------------------------------

# GotoFile dialog box with its helpers GotoFileApply, GotoFileUpdate

proc GotoFile {win} {
	global FileNameList NFileNames FileNameIndex
	set fl $win.filelist
	if {[winfo exists $fl]} {
		raise $fl
		return
	}
	toplevel $fl
	frame $fl.area
	set lb $fl.area.list
	set maxlen [MaxLength $FileNameList]
	if {$maxlen > 80} {set maxlen 80}
	set height [llength $FileNameList]
	if {$height < 7} {set height 7}
	if {$height > 20} {set height 20}
	listbox $lb -yscrollcommand [list $fl.area.sy set] -xscrollcommand [list $fl.area.sx set] -width $maxlen -height $height
	scrollbar $fl.area.sy -orient vertical -command [list $lb yview]
	scrollbar $fl.area.sx -orient horizontal -command [list $lb xview]
	GotoFileUpdate $win
	frame $fl.buttons
	button $fl.buttons.ok -text OK -underline 0 -command "[list GotoFileApply $win]; [list destroy $fl]"
	button $fl.buttons.apply -text Apply -underline 0 -command [list GotoFileApply $win]
	button $fl.buttons.load -text Load -underline 0 -command "LoadFile; [list GotoFileUpdate $win]; $lb selection set end"
	button $fl.buttons.cancel -text Cancel -underline 0 -command [list destroy $fl]
	pack $fl.area -side left -fill both -expand true
	pack $fl.area.sy -side right -fill y
	pack $fl.area.sx -side bottom -fill x
	pack $lb -side left -fill both -expand true
	pack $fl.buttons.ok $fl.buttons.apply $fl.buttons.load -side top
	pack $fl.buttons.cancel -side bottom
	pack $fl.buttons -side left -fill y
	bind $fl <Return> "[list GotoFileApply $win]; [list destroy $fl]"
	bind $fl <Double-1> "[list GotoFileApply $win]; [list destroy $fl]"
	bind $fl <Control-c> [list destroy $fl]
	bind $fl <Control-w> [list destroy $fl]
	bind $fl <Control-a> [list GotoFileApply $win]
	bind $fl <Control-l> "LoadFile; [list GotoFileUpdate $win]"
	SetWindowCloseToMouse $fl
	$lb selection set $FileNameIndex($win)
	$lb see $FileNameIndex($win)
}

proc GotoFileApply {win} {
	global FileNameIndex
	set i [$win.filelist.area.list curselection]
	# It is possible that the selection is empty. In that case, do nothing
	if {[string length $i] > 0} {
		set FileNameIndex($win) $i
		UpdatePossibleSlaves $win
	}
	Update $win
}

proc GotoFileUpdate {win} {
	# Fill in the listbox with entries from FileNameList.
	# Call this when first doing GotoFile dialog, and after modifying FileNameList
	# somehow (through Loadfile for instance).
	global FileNameList
	set fl $win.filelist
	set lb $fl.area.list
	$lb delete 0 end
	foreach file $FileNameList {$lb insert end $file}
}

# ------------------------------------------------------------------
# --- The 3D view point selector
# ------------------------------------------------------------------

set Edit3DViewpoint(realtime) 0

proc Edit3DViewpoint {win} {
	global Params Edit3DViewpoint
	set p $win.params
	if {[winfo exists $p]} {
		raise $p
		return
	}
	toplevel $p
	frame $p.a
	set a [frame $p.a.area]
	set la [frame $p.a.labelarea]
	set sa [frame $p.a.scalearea]
	set w 10
	set w2 13
	set g1t    [entry $la.g1     -textvariable Edit3DViewpoint(g1) -width $w2 -relief flat -state disabled]
	set rt     [entry $la.rt     -textvariable Edit3DViewpoint(rt) -width $w -relief flat -state disabled]
	set thetat [entry $la.thetat -textvariable Edit3DViewpoint(thetat) -width $w -relief flat -state disabled]
	set phit   [entry $la.phit   -textvariable Edit3DViewpoint(phit) -width $w -relief flat -state disabled]
	set g2t   [menubutton $la.g2t -text Center -menu $la.g2t.menu]
	set g2tm  [menu $g2t.menu]
	$g2tm add command -label {Grid center} -command "[list Edit3DGridCenterC $win]; [list Edit3DApply $win]"
	$g2tm add command -label {Origin     } -command "[list Edit3DOriginC $win]; [list Edit3DApply $win]"
	set cxt    [entry $la.cxt    -textvariable Edit3DViewpoint(cxt) -width $w -relief flat -state disabled]
	set cyt    [entry $la.cyt    -textvariable Edit3DViewpoint(cyt) -width $w -relief flat -state disabled]
	set czt    [entry $la.czt    -textvariable Edit3DViewpoint(czt) -width $w -relief flat -state disabled]
	set g3t    [entry $la.g3     -textvariable Edit3DViewpoint(g3) -width $w2 -relief flat -state disabled]
	set fovt   [entry $la.fovt   -textvariable Edit3DViewpoint(fovt) -width $w -relief flat -state disabled]
	set w 6
	set g1    [entry $a.g1 -width 1 -state disabled -relief flat]
	set r     [entry $a.r      -textvariable Params($win,relative_radius) -width $w]
	set theta [entry $a.theta  -textvariable Params($win,theta) -width $w]
	set phi   [entry $a.phi    -textvariable Params($win,phi) -width $w]
	set g2    [entry $a.g2 -width 1 -state disabled -relief flat]
	set cx    [entry $a.cx     -textvariable Params($win,cx) -width $w]
	set cy    [entry $a.cy     -textvariable Params($win,cy) -width $w]
	set cz    [entry $a.cz     -textvariable Params($win,cz) -width $w]
	set g3    [entry $a.g3 -width 1 -state disabled -relief flat]
	set fov   [entry $a.fov    -textvariable Params($win,fov) -width $w]
	set w 200
	set rs [scale $sa.rs -from 0.3 -to 3 -orient horizontal -resolution 0.1 \
			-variable Params($win,relative_radius) -command [list Edit3DRealTimeUpdate $win] \
			-label "Relative eye distance" -length $w -tickinterval 0.5 -font fixed]
	set thetas [scale $sa.thetas -from 0.01 -to 179.99 -orient horizontal \
				-variable Params($win,theta) -command [list Edit3DRealTimeUpdate $win] \
				-label "Theta (eye point colatitude)" -length $w -tickinterval 45 -font fixed]
	set phis   [scale $sa.phis   -from 0 -to 360 -orient horizontal \
				-variable Params($win,phi) -command [list Edit3DRealTimeUpdate $win] \
				-label "Phi (eye point longitude)" -length $w -tickinterval 90 -font fixed]
	set fovs   [scale $sa.fovs   -from 10 -to 80 -orient horizontal \
				-variable Params($win,fov) -command [list Edit3DRealTimeUpdate $win] \
				-label "Camera view angle" -length $w -tickinterval 20 -font fixed]
	set buttons [frame $p.buttons]
	set ok     [button $buttons.ok     -text OK       -command "[list Edit3DApply $win]; [list destroy $p]" -underline 0]
	set apply  [button $buttons.apply  -text Apply    -command [list Edit3DApply $win]]
	set rtmode [checkbutton $buttons.rtmode -text {Real time} \
				-variable Edit3DViewpoint(realtime) -command [list Edit3DRealTimeUpdate $win 0]]
	set def    [button $buttons.def    -text Defaults -command [list Default3DViewpoint $win]]
	set cancel [button $buttons.cancel -text Cancel   -command [list destroy $p] -underline 0]
	bind $p <Key-Return> "[list Edit3DApply $win]; [list destroy $p]"
	bind $p <Control-c> [list destroy $p]
	pack $p.a -side top
	pack $la -side left
	pack $g1t $rt $thetat $phit $g2t $cxt $cyt $czt $g3t $fovt -side top
	pack $a -side left -ipadx 10
	pack $g1 $r $theta $phi $g2 $cx $cy $cz $g3 $fov -side top
	pack $sa -side left
	pack $rs $thetas $phis $fovs -side top
	pack $buttons -side bottom
	pack $ok $apply $rtmode $def -side left
	pack $cancel -side right
	wm title $p {3D viewpoint selection}
}

proc Edit3DGridCenterC {win} {
	global Params
	array set A [$win.win hcgetbox]
	set Params($win,cx) [expr 0.5*($A(xmin) + $A(xmax))]
	set Params($win,cy) [expr 0.5*($A(ymin) + $A(ymax))]
	set Params($win,cz) [expr 0.5*($A(zmin) + $A(zmax))]
}

proc Edit3DOriginC {win} {
	global Params
	set Params($win,cx) 0
	set Params($win,cy) 0
	set Params($win,cz) 0
}

proc Edit3DApply {win} {
	global Params
	$win.win hcset3Dparams $Params($win,relative_radius) \
		$Params($win,theta) $Params($win,phi) \
		$Params($win,cx) $Params($win,cy) $Params($win,cz) \
		$Params($win,fov)
	Update $win
}

proc Edit3DRealTimeUpdate {win dummy} {
	global Edit3DViewpoint
	if {$Edit3DViewpoint(realtime)} {
		Edit3DApply $win
	}
}

proc Change3DViewpoint {win rfactor deltatheta deltaphi} {
	global Params
	set Params($win,relative_radius) [expr $rfactor * $Params($win,relative_radius)]
	set Params($win,theta) [expr round(($Params($win,theta) + $deltatheta))]
	if {$Params($win,theta) < 0.01} {set Params($win,theta) 0.01}
	if {$Params($win,theta) > 179.99} {set Params($win,theta) 179.99}
	set Params($win,phi) [expr ($Params($win,phi) + $deltaphi) % 360]
	Edit3DApply $win
}

proc Set3DAngles {win theta phi} {
	global Params
	set Params($win,theta) $theta
	if {$Params($win,theta) < 0.01} {set Params($win,theta) 0.01}
	if {$Params($win,theta) > 179.99} {set Params($win,theta) 179.99}
	set Params($win,phi) [expr $phi % 360]
	Edit3DApply $win
}

proc Xview {win} {
	global Flags Params
	if {$Flags($win,view3D)} {
		Set3DAngles $win 90 180
	} else {
		set Params($win,slicedim) 0
		Update $win
	}
}

proc Yview {win} {
	global Flags Params
	if {$Flags($win,view3D)} {
		Set3DAngles $win 90 270
	} else {
		set Params($win,slicedim) 1
		Update $win
	}
}

proc Zview {win} {
	global Flags Params
	if {$Flags($win,view3D)} {
		Set3DAngles $win 0.01 180
	} else {
		set Params($win,slicedim) 2
		Update $win
	}
}

proc Origin3DViewpoint {win} {
	global Params
	set Params($win,theta) $Params(default,theta)
	set Params($win,phi) $Params(default,phi)
	set Params($win,relative_radius) $Params(default,relative_radius)
	Edit3DOriginC $win
	Edit3DApply $win
}

proc GridCenter3DViewpoint {win} {
	global Params
	set Params($win,theta) $Params(default,theta)
	set Params($win,phi) $Params(default,phi)
	set Params($win,relative_radius) $Params(default,relative_radius)
	Edit3DGridCenterC $win
	Edit3DApply $win
}

proc Default3DViewpoint {win} {Origin3DViewpoint $win}

# ------------------------------------------------------------------
# --- The help function
# ------------------------------------------------------------------

proc abouthcvis {} {
	global HCVIS_ROOT
	set w .help
	if {[winfo exists $w]} {
		raise $w
		return
	}
	toplevel $w
	wm title $w "About hcvis"
	frame $w.buttons
	pack $w.buttons -side bottom -fill x
	button $w.buttons.ok -text OK -command [list destroy $w] -underline 0
	pack $w.buttons.ok -side left -expand true
	text $w.text -relief sunken -yscrollcommand "$w.scroll set" -setgrid true -height 24
	scrollbar $w.scroll -command [list $w.text yview]
	pack $w.scroll -side right -fill y
	pack $w.text -expand yes -fill both
	bind $w <Key-Return> [list destroy $w]
	$w.text insert 0.0 [exec cat $HCVIS_ROOT/hcvis.help]
	$w.text configure -state disabled
}

proc VariableNameList {win} {
#puts [$win.win hcintpol all 0 0 0]
	concat [list x y z] [TakeEverySecond [$win.win hcintpol all 0 0 0]]
}

proc VariableInfo {win} {
	set w .vardocs
	if {[winfo exists $w]} {
		raise $w
		return
	}
	toplevel $w
	wm title $w "Variable descriptions"
	frame $w.buttons
	pack $w.buttons -side bottom -fill x
	button $w.buttons.ok -text OK -command [list destroy $w] -underline 0
	pack $w.buttons.ok -side left -expand true
	text $w.text -relief sunken -yscrollcommand "$w.scroll set" -setgrid true -height 24
	scrollbar $w.scroll -command [list $w.text yview]
	pack $w.scroll -side right -fill y
	pack $w.text -expand yes -fill both
	bind $w <Key-Return> [list destroy $w]
	foreach var [VariableNameList $win] {
		$w.text insert end "$var: [$win.win hcvardescr $var]\n\n"
	}
#	$w.text insert 0.0 [exec cat $HCVIS_ROOT/hcvis.help]
	$w.text configure -state disabled
}

# ------------------------------------------------------------------
# --- The GetInfo function
# ------------------------------------------------------------------

proc GetInfo {win} {
	global FileNameList FileNameIndex
	set i $FileNameIndex($win)
	set fn [lindex $FileNameList $i]
	set w $win.info
	if {[winfo exists $w]} {
		raise $w
		return
	}
	toplevel $w
	wm title $w "Info about $fn"
	frame $w.buttons
	pack $w.buttons -side bottom -fill x
	button $w.buttons.ok -text OK -command [list destroy $w] -underline 0
	pack $w.buttons.ok -side left -expand true
	text $w.text -relief sunken -yscrollcommand "$w.scroll set" -setgrid true -height 24
	scrollbar $w.scroll -command [list $w.text yview]
	pack $w.scroll -side right -fill y
	pack $w.text -expand yes -fill both
	bind $w <Key-Return> [list destroy $w]
	$w.text insert 0.0 [exec awk "/^eoh$/ \{exit\} \{print\}" <$fn]
	$w.text configure -state disabled
}

# ------------------------------------------------------------------
# --- The DefineBox function
# ------------------------------------------------------------------

proc DefineBox {win} {
	global Flags DefineBox

	# If box seems to be undefined, initialize it by calling hcgetbox
	if {![info exists DefineBox(xmin)]} {
		array set DefineBox {}
		array set DefineBox [$win.win hcgetbox]
		# Also set the flag to True in this case.
		set Flags($win,BoxDefined) 1
	}

	set p $win.defbox
	if {[winfo exists $p]} {
		raise $p
		return
	}
	toplevel $p
	frame $p.a
	set a [frame $p.a.area]
	set w 10
	set w2 13
	set font fixed

	set xmin  [frame $a.xmin]
	set xminl [label $xmin.l -text xmin -font $font]
	set xmine [entry $xmin.e -textvariable DefineBox(xmin) -width $w]
	set xmax  [frame $a.xmax]
	set xmaxl [label $xmax.l -text xmax -font $font]
	set xmaxe [entry $xmax.e -textvariable DefineBox(xmax) -width $w]

	set ymin  [frame $a.ymin]
	set yminl [label $ymin.l -text ymin -font $font]
	set ymine [entry $ymin.e -textvariable DefineBox(ymin) -width $w]
	set ymax  [frame $a.ymax]
	set ymaxl [label $ymax.l -text ymax -font $font]
	set ymaxe [entry $ymax.e -textvariable DefineBox(ymax) -width $w]

	set zmin  [frame $a.zmin]
	set zminl [label $zmin.l -text zmin -font $font]
	set zmine [entry $zmin.e -textvariable DefineBox(zmin) -width $w]
	set zmax  [frame $a.zmax]
	set zmaxl [label $zmax.l -text zmax -font $font]
	set zmaxe [entry $zmax.e -textvariable DefineBox(zmax) -width $w]

	set buttons [frame $p.buttons]
	set onoff  [checkbutton $a.onoff -text {Box On} -variable Flags($win,BoxDefined)]
	set cancel [button $buttons.cancel -text Cancel -command [list destroy $p] -underline 0]
	set apply  [button $buttons.apply  -text Apply  -command [list DefineBoxApply $win]]
	set ok     [button $buttons.ok     -text OK     -command "[list DefineBoxApply $win]; [list destroy $p]" -underline 0]
	bind $p <Key-Return> "[list DefineBoxApply $win]; [list destroy $p]"
	bind $p <Control-c> [list destroy $p]

	pack $p.a -side top
	pack $a -side left
	pack $onoff $xmin $xmax $ymin $ymax $zmin $zmax -side top
	pack $xminl $xmine -side left
	pack $xmaxl $xmaxe -side left
	pack $yminl $ymine -side left
	pack $ymaxl $ymaxe -side left
	pack $zminl $zmine -side left
	pack $zmaxl $zmaxe -side left
	pack $buttons -side bottom
	pack $cancel $apply -side left
	pack $ok -side right
	wm title $p {Define box}
	SetWindowCloseToMouse $p
}

proc DefineBoxApply {win} {
	global DefineBox Flags
	if {$Flags($win,BoxDefined)} {
		$win.win hcsetbox $DefineBox(xmin) $DefineBox(xmax) $DefineBox(ymin) $DefineBox(ymax) $DefineBox(zmin) $DefineBox(zmax)
	} else {
		$win.win hcunsetbox
	}
	Update $win
}

# ------------------------------------------------------------------
# --- The SetTransparency function
# ------------------------------------------------------------------

proc SetTransparency {win} {
	global Params SetTransparency

	set p $win.alpha
	if {[winfo exists $p]} {
		raise $p
		return
	}
	toplevel $p
	set a [frame $p.a]
	set s [scale $a.s -from 1 -to 0 -resolution 0.01 -variable Params($win,alpha) \
			-label {Opacity (Alpha)} -orient horizontal -length 150 -tickinterval 0.5]
	set buttons [frame $p.buttons]
	set cancel [button $buttons.cancel -text Cancel -command [list destroy $p] -underline 0]
	set apply  [button $buttons.apply  -text Apply  -command [list Update $win]]
	set ok     [button $buttons.ok     -text OK     -command "[list Update $win]; [list destroy $p]" -underline 0]
	bind $p <Key-Return> "[list Update $win]; [list destroy $p]"
	bind $p <Control-c> [list destroy $p]
	pack $a -side top
	pack $s -side left
	pack $buttons -side bottom
	pack $cancel $apply -side left
	pack $ok -side right
	wm title $p {Set transparency}
	SetWindowCloseToMouse $p
}

# ------------------------------------------------------------------
# --- The SetFieldLineStopRadius function
# ------------------------------------------------------------------

proc SetFieldLineStopRadius {win} {
	global Params SetFieldLineStopRadius

	set p $win.setr
	if {[winfo exists $p]} {
		raise $p
		return
	}
	toplevel $p
	set a [frame $p.a]
	set rl [label $a.l -text {stop radius:}]
	set re [entry $a.e -textvariable Params($win,FieldLineStopSphereRadius)]
	set buttons [frame $p.buttons]
	set cancel [button $buttons.cancel -text Cancel -command [list destroy $p] -underline 0]
	set apply  [button $buttons.apply  -text Apply  -command [list Update $win]]
	set ok     [button $buttons.ok     -text OK     -command "[list Update $win]; [list destroy $p]" -underline 0]
	bind $p <Key-Return> "[list Update $win]; [list destroy $p]"
	bind $p <Control-c> [list destroy $p]
	pack $a -side top
	pack $rl $re -side left
	pack $buttons -side bottom
	pack $cancel $apply -side left
	pack $ok -side right
	wm title $p {Field line stop sphere radius}
#	SetWindowCloseToMouse $p
}

# ------------------------------------------------------------------
# --- The SetIsoValue function
# ------------------------------------------------------------------

proc SetIsoValue {win} {
	global Params Flags
	set p $win.isoval
	if {[winfo exists $p]} {
		raise $p
		return
	}
	toplevel $p
	set a0 [frame $p.a0]
	set rl0 [label $a0.l -text {isosurface 1 value:}]
	set re0 [entry $a0.e -textvariable Params($win,IsoValue0)]
	set a1 [frame $p.a1]
	set rl1 [label $a1.l -text {isosurface 2 value:}]
	set re1 [entry $a1.e -textvariable Params($win,IsoValue1)]
	set a2 [frame $p.a2]
	set rl2 [label $a2.l -text {isosurface 3 value:}]
	set re2 [entry $a2.e -textvariable Params($win,IsoValue2)]
	set a3 [frame $p.a3]
	set rl3 [label $a3.l -text {isosurface 4 value:}]
	set re3 [entry $a3.e -textvariable Params($win,IsoValue3)]
	set a4 [frame $p.a4]
	set rl4 [label $a4.l -text {isosurface 5 value:}]
	set re4 [entry $a4.e -textvariable Params($win,IsoValue4)]
	set buttons [frame $p.buttons]
	set cancel [button $buttons.cancel -text Cancel -command [list destroy $p] -underline 0]
	if {$Flags($win,DrawIsosurf)} {
		set apply  [button $buttons.apply  -text Apply  -command [list Update $win]]
		set ok     [button $buttons.ok     -text OK     -command "[list Update $win]; [list destroy $p]" -underline 0]
		bind $p <Key-Return> "[list Update $win]; [list destroy $p]"
	} else {
		set apply  [button $buttons.apply  -text Apply]
		set ok     [button $buttons.ok     -text OK     -command [list destroy $p] -underline 0]
		bind $p <Key-Return> [list destroy $p]
	}
	bind $p <Control-c> [list destroy $p]
	pack $a0 $a1 $a2 $a3 $a4 -side top
	pack $rl0 $re0 -side left
	pack $rl1 $re1 -side left
	pack $rl2 $re2 -side left
	pack $rl3 $re3 -side left
	pack $rl4 $re4 -side left
	pack $buttons -side bottom
	pack $cancel $apply -side left
	pack $ok -side right
	wm title $p {Isosurface values}
#	SetWindowCloseToMouse $p
}

# ------------------------------------------------------------------
# --- The VolumetricParameters function
# ------------------------------------------------------------------

proc VolumetricParameters {win} {
	global Params Flags VolumetricParameters

	set p $win.volumetric
	if {[winfo exists $p]} {
		raise $p
		return
	}
	toplevel $p
	set a [frame $p.a]
	set Params($win,VolumetricNpointsMillions) [expr 1e-6 * $Params($win,VolumetricNpoints)]
	set s [scale $a.s -from 1 -to 0 -resolution 0.01 -variable Params($win,VolumetricAlpha) \
			-label {Opacity (Alpha)} -orient horizontal -length 250 -tickinterval 0.5]
	set ss [scale $a.ss -from 0 -to 2 -resolution 0.05 -variable Params($win,VolumetricNpointsMillions) \
			-label {Number of points (Millions)} -orient horizontal -length 250 -tickinterval 0.5]
	set aa [checkbutton $a.aa -text Antialiasing -variable Flags($win,VolumetricAntialiasing)]
	set buttons [frame $p.buttons]
	set cancel [button $buttons.cancel -text Cancel -command [list destroy $p] -underline 0]
	set apply  [button $buttons.apply  -text Apply  -command [list Update $win]]
	set ok     [button $buttons.ok     -text OK     -command "[list Update $win]; [list destroy $p]" -underline 0]
	bind $p <Key-Return> "[list Update $win]; [list destroy $p]"
	bind $p <Control-c> [list destroy $p]
	pack $a -side top
	pack $s $ss $aa -side top
	pack $buttons -side bottom
	pack $cancel $apply -side left
	pack $ok -side right
	wm title $p {Volumetric visualization parameters}
#	SetWindowCloseToMouse $p
}

proc SetMinMax {win} {
	global Params SetMinMax
	set p $win.minmax
	if {[winfo exists $p]} {
		raise $p
		return
	}
	toplevel $p
	set a [frame $p.a]
	set min [frame $a.min]
	set max [frame $a.max]
	set minl [label $min.l -text {min:}]
	set mine [entry $min.e -textvariable Params($win,min)]
	set maxl [label $max.l -text {max:}]
	set maxe [entry $max.e -textvariable Params($win,max)]
	set buttons [frame $p.buttons]
	set cancel [button $buttons.cancel -text Cancel -command [list destroy $p] -underline 0]
	set apply  [button $buttons.apply  -text Apply  -command [list Update $win]]
	set ok     [button $buttons.ok     -text OK     -command "[list Update $win]; [list destroy $p]" -underline 0]
	bind $p <Key-Return> "[list Update $win]; [list destroy $p]"
	bind $p <Control-c> [list destroy $p]
	bind $p <Control-a> [list Update $win]
	pack $a -side top
	pack $min $max -side top
	pack $minl $mine -side left
	pack $maxl $maxe -side left
	pack $buttons -side bottom
	pack $cancel $apply -side left
	pack $ok -side right
	wm title $p {Set min/max}
	SetWindowCloseToMouse $p
}

# ------------------------------------------------------------------
# --- Running (Interball tail) orbit data file
# ------------------------------------------------------------------

proc RunOrbitFile {win} {
	set fn [fileselect {Orbit data file (format: x,y,z items 6,7,8 on each line)} {} 1 1]
	if {[string length $fn] == 0} {return}
	set data [exec awk "\{printf \"\{%g %g %g\}\\n\",\$6,\$7,\$8\}" <$fn]
	$win.win hcset3Dparams 1 45 45 0 0 0 60
	set coeff 1
	set counter 1
	exec rm -f orbit?????.ppm
	foreach line $data {
		set x [expr $coeff*[lindex $line 0]]
		set y [expr $coeff*[lindex $line 1]]
		set z [expr $coeff*[lindex $line 2]]
		puts "x=$x, y=$y, z=$z"
		$win.win hcset3Deyepoint $x $y $z
		Update $win
		set ppmfile [format "orbit%.5d.ppm" $counter]
		$win.win hcwriteppm $ppmfile
		puts "Wrote $ppmfile"
		incr counter
	}
}

# ---------------------------------------------------------------------------
# --- Procesing mouse clicks in 2D view windows when Variable window is open
# ---------------------------------------------------------------------------

proc MouseClick {win mx my} {
	global MouseClickVars
	set xyz [$win.win hcmodelcoords2D $mx $my]
	array set vars {}
	set x [lindex $xyz 0]
	set vars(x) $x
	set y [lindex $xyz 1]
	set vars(y) $y
	set z [lindex $xyz 2]
	set vars(z) $z
	array set vars [$win.win hcintpol all $x $y $z]
	foreach var [array names vars] {
		set MouseClickVars($win,$var) $vars($var)
	}
}

proc VariableWindow {win} {
	global MouseClickVars
	set p $win.vars
	if {[winfo exists $p]} {
		raise $p
		return
	}
	toplevel $p
	wm title $p "Variable view"
	set entries [frame $p.entries]
	set VariableNames [VariableNameList $win]
	set nrows 13
	set ncols [expr int(floor([llength $VariableNames] / $nrows)) + 1]
	for {set i 0} {$i < $ncols} {incr i} {frame $p.entries.col$i}
#	puts "ncols=$ncols, nrows=$nrows, [llength $VariableNames] variables"
	set i 0
	foreach var $VariableNames {
		set col [expr int(floor($i / $nrows))]
		set lcvar [string tolower $var]
		set e $entries.col$col
		frame $e.$lcvar -borderwidth 0
#		entry $e.$lcvar.name  -text "$var = " -relief flat -state disabled
		label $e.$lcvar.name -text "$var = " -font fixed -pady 0
		entry $e.$lcvar.value -textvariable MouseClickVars($win,$var) -relief flat \
			-state disabled -width 12 -font fixed -borderwidth 0
		incr i
	}
	set buttons [frame $p.buttons]
	set ok [button $buttons.ok -text Close -command [list destroy $p] -underline 0]
	bind $p <Key-Return> [list destroy $p]
	pack $entries -side top
	for {set i 0} {$i < $ncols} {incr i} {pack $p.entries.col$i -side left -anchor n}
	set i 0
	foreach var $VariableNames {
		set col [expr int(floor($i / $nrows))]
		set lcvar [string tolower $var]
		set e $entries.col$col
		pack $e.$lcvar -side top
		pack $e.$lcvar.name $e.$lcvar.value -side left
		incr i
	}
	pack $buttons -side bottom
	pack $ok -side left -expand true
}

# ------------------------------------------------------------------
# --- Processing mouse clicks for zooming in 2D view windows
# ------------------------------------------------------------------

proc MouseDown {win mx my} {
	global MouseClickVars
	if {[winfo exists $win.vars]} {
		MouseClick $win $mx $my
	}
	set MouseClickVars($win,mx1) $mx
	set MouseClickVars($win,my1) $my
}

proc MouseUp {win mx my} {
	global MouseClickVars Params
	if {[winfo exists $win.vars]} {return}
	set MouseClickVars($win,mx2) $mx
	set MouseClickVars($win,my2) $my
	set limit 10
	set mx1 $MouseClickVars($win,mx1)
	set my1 $MouseClickVars($win,my1)
	if {abs($mx-$mx1) > $limit && abs($my-$my1) > $limit} {
		set saved [$win.win hcgetzoom]
		set Params($win,zoomstack) [concat [list $saved] $Params($win,zoomstack)]
		set xyzwh1 [$win.win hcmodelcoords2D $mx1 $my1]
		set xyzwh2 [$win.win hcmodelcoords2D $mx $my]
		set wmin [lindex $xyzwh1 3]
		set hmin [lindex $xyzwh1 4]
		set wmax [lindex $xyzwh2 3]
		set hmax [lindex $xyzwh2 4]
		$win.win hczoom $wmin $hmin $wmax $hmax
		Update $win
	}
}

proc PreviousZoom {win} {
	global Params
	if {[llength $Params($win,zoomstack)] > 0} {
		set zoomrect [lindex $Params($win,zoomstack) 0]
		set Params($win,zoomstack) [lrange $Params($win,zoomstack) 1 end]
		$win.win hczoom [lindex $zoomrect 0] [lindex $zoomrect 1] [lindex $zoomrect 2] [lindex $zoomrect 3]
		Update $win
	}
}

proc FullZoom {win} {
	global Params
	set Params($win,zoomstack) {}
	$win.win hcfullzoom
	Update $win
}

# ------------------------------------------------------------------
# --- One-dimensional view (uses BLT graph widget)
# ------------------------------------------------------------------

#proc NRange {a b n} {
#	set xx {}
#	set x $a
#	set dx [expr ($b-$a)/double($n)]
#	set i 0
#	while {$i<$n} {
#		lappend xx $x
#		set x [expr $x+$dx]
#		incr i
#	}
#	return $xx
#}

proc NRange3D {x1 y1 z1 x2 y2 z2 n} {
	set L {}
	set dx [expr ($x2-$x1)/double($n)]
	set dy [expr ($y2-$y1)/double($n)]
	set dz [expr ($z2-$z1)/double($n)]
	for {set i 0} {$i<$n} {incr i} {
		set x [expr $x1+$i*$dx]
		set y [expr $y1+$i*$dy]
		set z [expr $z1+$i*$dz]
		lappend L [list $x $y $z]
	}
	return $L
}

proc OneDimensionalViewUpdate {win} {
	# This is called from Update, and also from OneDimensionalView
	global OneDimensionalView Params FileNameIndex FileNameList
	set p $win.oned
	set i $FileNameIndex($win)
	set fn [Tail [lindex $FileNameList $i]]
	wm title $p $fn
	set v $Params($win,var1)
	set w $win.win
	set x1 $OneDimensionalView($win,x1)
	set y1 $OneDimensionalView($win,y1)
	set z1 $OneDimensionalView($win,z1)
	set x2 $x1
	set y2 $y1
	set z2 $z1
	set d $OneDimensionalView($win,dir)
	set sign $OneDimensionalView($win,sign)
	set len $OneDimensionalView($win,len)
	switch $sign {
		Positive {
			switch $d {
				0 {set x2 [expr $x1+$len]}
				1 {set y2 [expr $y1+$len]}
				2 {set z2 [expr $z1+$len]}
			}
		}
		Negative {
			switch $d {
				0 {set x2 [expr $x1-$len]}
				1 {set y2 [expr $y1-$len]}
				2 {set z2 [expr $z1-$len]}
			}
		}
		Both {
			switch $d {
				0 {set x2 [expr $x1+$len]; set x1 [expr $x1-$len]}
				1 {set y2 [expr $y1+$len]; set y1 [expr $y1-$len]}
				2 {set z2 [expr $z1+$len]; set z1 [expr $z1-$len]}
			}
		}
	}	

	set OneDimensionalView($win,vec) [NRange3D  $x1 $y1 $z1  $x2 $y2 $z2  $OneDimensionalView($win,n)]
	set Yvec {}
	set AxisNames {x y z}
	foreach r $OneDimensionalView($win,vec) {
		set x [lindex $r 0]
		set y [lindex $r 1]
		set z [lindex $r 2]
		set Yval [$w hcintpol $v $x $y $z]
		if {![string match NaN* $Yval]} {
			lappend Yvec $Yval
			lappend Xvec [lindex $r $d]
		}
	}
	$p.a.g element configure "elementti1" -xdata $Xvec -ydata $Yvec -label $v
	$p.a.g axis configure x -title [lindex $AxisNames $d]
	$p.a.g configure -title $fn
}

proc OneDimensionalViewExport {win} {
	global OneDimensionalView
	set w $win.win
	set varnames [TakeEverySecond [$w hcintpol all 0 0 0]]
	set fn {export.dat}
	set fn [fileselect "Export 1D ASCII data to" $fn 0 1 0]
	if {[string length $fn] <= 0} {return}
	set fp [open $fn w]
	puts $fp "# x y z $varnames"
	foreach r $OneDimensionalView($win,vec) {
		set x [lindex $r 0]
		set y [lindex $r 1]
		set z [lindex $r 2]
		set varvalues [TakeEverySecond [lrange [$w hcintpol all $x $y $z] 1 end]]
		puts $fp "$x $y $z $varvalues"
	}
	close $fp
}

proc TwoDimensionalViewExport {win} {
	global Params
	set wi $win.win
	set fn {export.dat}
	set fn [fileselect "Export 2D ASCII $Params($win,var) to" $fn 0 1 0]
	$wi hcexportascii $fn
}

proc FastOneDimensionalViewUpdate {win {dummy 0}} {
	global OneDimensionalView
	if {$OneDimensionalView($win,realtime)} {OneDimensionalViewUpdate $win}
}


proc OneDimensionalViewUpdateRealtimeMode {win} {
	global OneDimensionalView
	array set box [$win.win hcgetbox]
	set Lx [expr $box(xmax)-$box(xmin)]
	set Ly [expr $box(ymax)-$box(ymin)]
	set Lz [expr $box(zmax)-$box(zmin)]
	set Lmax [expr $Lx>$Ly ? $Lx : $Ly]
	set Lmax [expr $Lmax>$Lz ? $Lmax : $Lz]
	if {$OneDimensionalView($win,realtime)} {
		set OneDimensionalView($win,rel_resol) 0.05
	} else {
		set OneDimensionalView($win,rel_resol) 0.001
	}
	set scales $win.oned.scales
	$scales.x1 configure -resolution [expr $OneDimensionalView($win,rel_resol)*$Lx]
	$scales.y1 configure -resolution [expr $OneDimensionalView($win,rel_resol)*$Ly]
	$scales.z1 configure -resolution [expr $OneDimensionalView($win,rel_resol)*$Lz]
	$scales.len configure -resolution [expr $OneDimensionalView($win,rel_resol)*$Lmax]
	OneDimensionalViewUpdate $win
}

proc OneDimensionalView {win} {
	global Params OneDimensionalView
	if {[llength [info commands graph]]==0} {
		puts "BLT library not available"
		return
	}
	set p $win.oned
	if {[winfo exists $p]} {
		raise $p
		return
	}
	toplevel $p
	set a [frame $p.a]
	set g [graph $a.g -title "Kuvio"]
	$g element create "elementti1" -symbol none

	array set box [$win.win hcgetbox]
	# --- default parameters	
	set OneDimensionalView($win,x1) $box(xmin)
	set OneDimensionalView($win,y1) [expr 0.5*($box(ymin)+$box(ymax))]
	set OneDimensionalView($win,z1) $box(zmin)
	set OneDimensionalView($win,dir) 0
	set OneDimensionalView($win,sign) Positive
	set OneDimensionalView($win,n) 100
	set OneDimensionalView($win,len) [expr $box(xmax)-$box(xmin)]
	set OneDimensionalView($win,rel_resol) 0.01
	set OneDimensionalView($win,realtime) 0
	# --- call hook proc if defined to set other defaults
	if {[llength [info commands OneDimensionalViewInitParams]]==1} {
		OneDimensionalViewInitParams $win
	}

	set Params($win,var1) $Params($win,var)
	OneDimensionalViewUpdate $win
	set menus [frame $p.menus]
	set buttons [frame $p.buttons]
	set scales [frame $p.scales]

	set w 260

	set Lx [expr $box(xmax)-$box(xmin)]
	set Ly [expr $box(ymax)-$box(ymin)]
	set Lz [expr $box(zmax)-$box(zmin)]
	set Lmax [expr $Lx>$Ly ? $Lx : $Ly]
	set Lmax [expr $Lmax>$Lz ? $Lmax : $Lz]

	set varselect [menubutton $menus.vars -text Variables -menu $menus.vars.menu]
	set vm [menu $menus.vars.menu]
	VariablesMenu $win $vm Params($win,var1) [list OneDimensionalViewUpdate $win]

	scale $scales.x1 -from $box(xmin) -to $box(xmax) -orient horizontal \
		-resolution [expr $OneDimensionalView($win,rel_resol)*$Lx] \
		-variable OneDimensionalView($win,x1) -command [list FastOneDimensionalViewUpdate $win] \
		-label "x1" -length $w

	scale $scales.y1 -from $box(ymin) -to $box(ymax) -orient horizontal \
		-resolution [expr $OneDimensionalView($win,rel_resol)*$Ly] \
		-variable OneDimensionalView($win,y1) -command [list FastOneDimensionalViewUpdate $win] \
		-label "y1" -length $w

	scale $scales.z1 -from $box(zmin) -to $box(zmax) -orient horizontal \
		-resolution [expr $OneDimensionalView($win,rel_resol)*$Lz] \
		-variable OneDimensionalView($win,z1) -command [list FastOneDimensionalViewUpdate $win] \
		-label "z1" -length $w

	scale $scales.len -from [expr 0.01*$Lmax] -to $Lmax -orient horizontal -resolution [expr 0.01*$Lmax] \
			-variable OneDimensionalView($win,len) -command [list FastOneDimensionalViewUpdate $win] \
			-label "Length of 1D axis" -length $w

	set dselect [menubutton $menus.dselect -text Direction -menu $menus.dselect.menu]
	menu $dselect.menu
	$dselect.menu add radiobutton -label {X} -variable OneDimensionalView($win,dir) \
		-value 0 -command [list OneDimensionalViewUpdate $win]
	$dselect.menu add radiobutton -label {Y} -variable OneDimensionalView($win,dir) \
		-value 1 -command [list OneDimensionalViewUpdate $win]
	$dselect.menu add radiobutton -label {Z} -variable OneDimensionalView($win,dir) \
		-value 2 -command [list OneDimensionalViewUpdate $win]
	$dselect.menu add radiobutton -label {Positive} -variable OneDimensionalView($win,sign) \
		-value Positive -command [list OneDimensionalViewUpdate $win]
	$dselect.menu add radiobutton -label {Negative} -variable OneDimensionalView($win,sign) \
		-value Negative -command [list OneDimensionalViewUpdate $win]
	$dselect.menu add radiobutton -label {Both} -variable OneDimensionalView($win,sign) \
		-value Both -command [list OneDimensionalViewUpdate $win]

	set nselect [menubutton $menus.nselect -text Npoints -menu $menus.nselect.menu]
	menu $nselect.menu
	$nselect.menu add radiobutton -label 10 -variable OneDimensionalView($win,n) \
		-value 10 -command [list OneDimensionalViewUpdate $win]
	$nselect.menu add radiobutton -label 20 -variable OneDimensionalView($win,n) \
		-value 20 -command [list OneDimensionalViewUpdate $win]
	$nselect.menu add radiobutton -label 50 -variable OneDimensionalView($win,n) \
		-value 50 -command [list OneDimensionalViewUpdate $win]
	$nselect.menu add radiobutton -label 100 -variable OneDimensionalView($win,n) \
		-value 100 -command [list OneDimensionalViewUpdate $win]
	$nselect.menu add radiobutton -label 200 -variable OneDimensionalView($win,n) \
		-value 200 -command [list OneDimensionalViewUpdate $win]
	$nselect.menu add radiobutton -label 500 -variable OneDimensionalView($win,n) \
		-value 500 -command [list OneDimensionalViewUpdate $win]
	$nselect.menu add radiobutton -label 1000 -variable OneDimensionalView($win,n) \
		-value 1000 -command [list OneDimensionalViewUpdate $win]
	$nselect.menu add radiobutton -label 2000 -variable OneDimensionalView($win,n) \
		-value 2000 -command [list OneDimensionalViewUpdate $win]
	$nselect.menu add radiobutton -label 5000 -variable OneDimensionalView($win,n) \
		-value 5000 -command [list OneDimensionalViewUpdate $win]
	$nselect.menu add radiobutton -label 10000 -variable OneDimensionalView($win,n) \
		-value 10000 -command [list OneDimensionalViewUpdate $win]

	set rtselect [checkbutton $buttons.rtmode -text {Real time} -variable OneDimensionalView($win,realtime) \
		-command [list OneDimensionalViewUpdateRealtimeMode $win]]

	set export [button $buttons.export -text {Export..} -command [list OneDimensionalViewExport $win]]
	set apply [button $buttons.apply -text Apply -command [list OneDimensionalViewUpdate $win] -underline 0]
	set close [button $buttons.close -text Close -command [list destroy $p] -underline 0]
	bind $p <Control-Key-w> [list destroy $p]
	bind $p <Key-Return> [list OneDimensionalViewUpdate $win]
	bind $p <Control-Key-a> [list OneDimensionalViewUpdate $win]
	bind $p <Control-Key-c> [list destroy $p]
	pack $menus $a $scales -side top
	pack $varselect $dselect $nselect -side left
	pack $g -side top
	pack $scales.x1 $scales.y1 $scales.z1 $scales.len -side top
	pack $buttons -side bottom
	pack $rtselect $export $apply $close -side left
}

proc WriteEPS {win} {
	global FileNameList FileNameIndex
#	puts "WriteEPS $win"
	set fn "[file rootname [file tail [lindex $FileNameList $FileNameIndex($win)]]].eps"
	$win.win hcwriteeps $fn
}

proc WritePPM {win} {
	global FileNameList FileNameIndex
#	puts "WritePPM $win"
	set fn "[file rootname [file tail [lindex $FileNameList $FileNameIndex($win)]]].ppm"
	$win.win hcwriteppm $fn
}

proc WriteAllPPM {win} {
	global FileNameIndex NFileNames
	set saveok [YesNoQuery {Go through all files to the end
and dump them as PPM images ?}]
	if {$saveok} {
		WritePPM $win
		while {$FileNameIndex($win) < $NFileNames - 1} {
			NextFile $win
			WritePPM $win
		}
	}
}

# ------------------------------------------------------------------
# --- The main menu
# ------------------------------------------------------------------

proc ToggleFlag {win flag} {
	global Flags
	set Flags($win,$flag) [expr !$Flags($win,$flag)]
	Update $win
}

proc NewVariable {win} {
	global Params
	set Params($win,min) {}
	set Params($win,max) {}
	set Params($win,var1) $Params($win,var)
	Update $win
}

proc VariablesMenu {f vm params update {commandbuttons 0}} {
	$vm add radiobutton -label "rho           " -variable $params -value rho -command $update
	$vm add radiobutton -label "n (density)   " -variable $params -value n -command $update
	$vm add radiobutton -label "vx            " -variable $params -value vx -command $update
	$vm add radiobutton -label "vy            " -variable $params -value vy -command $update
	$vm add radiobutton -label "vz            " -variable $params -value vz -command $update
	$vm add radiobutton -label "vr            " -variable $params -value vr -command $update
	$vm add radiobutton -label "v^2           " -variable $params -value v2 -command $update
	$vm add radiobutton -label "|v|           " -variable $params -value v -command $update
	$vm add radiobutton -label "P             " -variable $params -value P -command $update
	$vm add radiobutton -label "T             " -variable $params -value T -command $update
	set m $vm.conserv
	$vm add cascade     -label "Conservative  " -menu $m
	menu $m
	$m add radiobutton -label "rhovx      " -variable $params -value rhovx -command $update
	$m add radiobutton -label "rhovy      " -variable $params -value rhovy -command $update
	$m add radiobutton -label "rhovz      " -variable $params -value rhovz -command $update
	$m add radiobutton -label "nv         " -variable $params -value nv -command $update
	$m add radiobutton -label "U (energy) " -variable $params -value U -command $update
	$m add radiobutton -label "U1 (energy)" -variable $params -value U1 -command $update
	$m add radiobutton -label "Uminusdip  " -variable $params -value Uminusdip -command $update
	$m add radiobutton -label "Kx         " -variable $params -value Kx -command $update
	$m add radiobutton -label "Ky         " -variable $params -value Ky -command $update
	$m add radiobutton -label "Kz         " -variable $params -value Kz -command $update
	$m add radiobutton -label "|K|        " -variable $params -value K -command $update
	set m $vm.numbers
	$vm add cascade     -label "Numbers etc   " -menu $m
	menu $m
	$m add radiobutton -label "a (sonic)  "  -variable $params -value a -command $update
	$m add radiobutton -label "vA (Alfven)"  -variable $params -value vA -command $update
	$m add radiobutton -label "vf (fast)  "  -variable $params -value vf -command $update
	$m add radiobutton -label "Ms=v/a     "  -variable $params -value Ms -command $update
	$m add radiobutton -label "MA=v/vA    "  -variable $params -value MA -command $update
	$m add radiobutton -label "beta       "  -variable $params -value beta -command $update
	$m add radiobutton -label "Comp. work "  -variable $params -value work -command $update
	set m $vm.plasma
	$vm add cascade     -label "Plasma params " -menu $m
	menu $m
	$m add radiobutton -label "rLp (Prot Larmor) "  -variable $params -value rLp -command $update
	$m add radiobutton -label "rLp/dx            "  -variable $params -value rLp/dx -command $update
	$m add radiobutton -label "rLe (Ele Larmor)  "  -variable $params -value rLe -command $update
	$m add radiobutton -label "fpe (Plasma freq) "  -variable $params -value fpe -command $update
	$m add radiobutton -label "rDebye            "  -variable $params -value rDebye -command $update
	$m add radiobutton -label "fp (Prot gyrofreq)"  -variable $params -value fgyrop -command $update
	$m add radiobutton -label "fe (Ele gyrofreq) "  -variable $params -value fgyroe -command $update
	$m add radiobutton -label "c/wpi (Ion inert) "  -variable $params -value c/wpi -command $update
	$m add radiobutton -label "c/wpe (Ele inert) "  -variable $params -value c/wpe -command $update
	$m add radiobutton -label "fLH (Lower hybrid)"  -variable $params -value fLH -command $update
	set m $vm.magnetic
	$vm add cascade     -label "Magnetic      " -menu $m
	menu $m
	$m add radiobutton -label "Bx         " -variable $params -value Bx -command $update
	$m add radiobutton -label "By         " -variable $params -value By -command $update
	$m add radiobutton -label "Bz         " -variable $params -value Bz -command $update
	$m add radiobutton -label "Bx0        " -variable $params -value Bx0 -command $update
	$m add radiobutton -label "By0        " -variable $params -value By0 -command $update
	$m add radiobutton -label "Bz0        " -variable $params -value Bz0 -command $update
	$m add radiobutton -label "Bx1        " -variable $params -value Bx1 -command $update
	$m add radiobutton -label "By1        " -variable $params -value By1 -command $update
	$m add radiobutton -label "Bz1        " -variable $params -value Bz1 -command $update
	$m add separator
	$m add radiobutton -label "|B|        " -variable $params -value B -command $update
	$m add radiobutton -label "|B0|       " -variable $params -value B0 -command $update
	$m add radiobutton -label "|B1|       " -variable $params -value B1 -command $update
	$m add radiobutton -label "|B|-|B0|   " -variable $params -value BminusB0 -command $update
	$m add radiobutton -label "B^2        " -variable $params -value B2 -command $update
	$m add radiobutton -label "B0^2       " -variable $params -value B02 -command $update
	$m add radiobutton -label "B1^2       " -variable $params -value B12 -command $update
	$m add separator
	$m add radiobutton -label "divB       " -variable $params -value divB -command $update
	$m add radiobutton -label "divB_rel   " -variable $params -value divBrel -command $update
	set m $vm.current
	$vm add cascade    -label "Current        " -menu $m
	menu $m
	$m add radiobutton -label "jx         " -variable $params -value jx -command $update
	$m add radiobutton -label "jy         " -variable $params -value jy -command $update
	$m add radiobutton -label "jz         " -variable $params -value jz -command $update
	$m add radiobutton -label "j          " -variable $params -value j  -command $update
	$m add radiobutton -label "j/(enV)    " -variable $params -value j/env -command $update
	$m add radiobutton -label "jPar       " -variable $params -value jPar -command $update
	set m $vm.electric
	$vm add cascade    -label "Electric       " -menu $m
	menu $m
	$m add radiobutton -label "Ex         " -variable $params -value Ex -command $update
	$m add radiobutton -label "Ey         " -variable $params -value Ey -command $update
	$m add radiobutton -label "Ez         " -variable $params -value Ez -command $update
	$m add radiobutton -label "E          " -variable $params -value E -command $update
	$m add radiobutton -label "Sx (Poynt.)" -variable $params -value Sx -command $update
	$m add radiobutton -label "Sy         " -variable $params -value Sy -command $update
	$m add radiobutton -label "Sz         " -variable $params -value Sz -command $update
	$m add radiobutton -label "S          " -variable $params -value S -command $update
	if {$commandbuttons} {
		$vm add checkbutton -label "Logarithmic   L " -variable Flags($f,Logarithmic) -command $update
		$vm add command     -label "Min/Max.. Ctrl+M" -command [list SetMinMax $f]
	}
}

proc SetupMenus {f} {
	global Flags env Palette PaletteDir statefile ZoomInFactor
	frame $f.menubar -relief raised -bd 2
	pack $f.menubar -side top -fill x

	set menufont fixed
	# File menu
	menubutton $f.menubar.file -text File -menu $f.menubar.file.menu -font $menufont
	set fm [menu $f.menubar.file.menu -font $menufont]
	$fm add command -label "Next file         N,Space" -command [list NextFile $f]
	$fm add command -label "Previous file   Backspace" -command [list PrevFile $f]
	$fm add command -label "Goto file ...      Ctrl+G" -command [list GotoFile $f]
	$fm add separator
	$fm add command -label "Load ...           Ctrl+L" -command [list LoadFile 1 $f]
	$fm add command -label "Close              Ctrl+W" -command [list CloseWindow $f]
	$fm add command -label "New view window    Ctrl+N" -command [list NewWindow $f]
	$fm add command -label "New slave window   Ctrl+D" -command [list NewSlaveWindow $f]
#	$fm add command -label "Write EPS file           " -command [list WriteEPS $f]
	$fm add command -label "Write PPM file           " -command [list WritePPM $f]
	$fm add command -label "Write all files as PPM .." -command [list WriteAllPPM $f]
	$fm add separator
	$fm add command -label "Save state        Shift+S" -command [list SaveState $statefile]
	$fm add command -label "Restore state     Shift+R" -command [list LoadState $statefile]
	$fm add command -label "Save state to ..   Ctrl+S" -command {QuerySaveState}
	$fm add command -label "Restore from ..    Ctrl+R" -command {QueryLoadState}
	$fm add separator
	$fm add command -label "Edit objects ...   Ctrl+E" -command [list EditObjects $f]
	$fm add command -label "Variable values    Ctrl+V" -command [list VariableWindow $f]
	$fm add command -label "1D view window     Ctrl+1" -command [list OneDimensionalView $f]
	$fm add command -label "Run orbit file ...       " -command [list RunOrbitFile $f]
	$fm add command -label "Export 2D data ...       " -command [list TwoDimensionalViewExport $f]
	$fm add separator
	$fm add command -label "Quit               Ctrl+Q" -command {exit}
	bind $f <Key-n>     [list NextFile $f]
	bind $f <Key-space> [list NextFile $f]
	bind $f <BackSpace> [list PrevFile $f]
	bind $f <Control-p> [list PrevFile $f]
	bind $f <Control-g> [list GotoFile $f]
	bind $f <Control-l> [list LoadFile 1 $f]
	bind $f <Control-w> [list CloseWindow $f]
	bind $f <Control-n> [list NewWindow $f]
	bind $f <Control-d> [list NewSlaveWindow $f]
	bind $f <Key-S>     {SaveState $statefile}
	bind $f <Key-R>     {LoadState $statefile}
	bind $f <Control-s> {QuerySaveState}
	bind $f <Control-r> {QueryLoadState}
	bind $f <Control-e> [list EditObjects $f]
	bind $f <Control-Key-1> [list OneDimensionalView $f]
	bind $f <Control-i> [list GetInfo $f]

	# Options menu
	menubutton $f.menubar.options -text Options -menu $f.menubar.options.menu -font $menufont
	set om [menu $f.menubar.options.menu -font $menufont]
	$om add checkbutton -label "Linear interpolation   I " -variable Flags($f,LinearInterpolation) -command [list Update $f]
	$om add checkbutton -label "Show ghost cell data     " -variable Flags($f,DrawGhosts) -command [list Update $f]
	$om add checkbutton -label "Hide ghost cells       H " -variable Flags($f,HideGhosts) -command [list Update $f]
	$om add checkbutton -label "Draw grid              G " -variable Flags($f,DrawGrid) -command [list Update $f]
	$om add checkbutton -label "Show dead cell data      " -variable Flags($f,DrawDead) -command [list Update $f]
	$om add checkbutton -label "Draw contours     Shift+C" -variable Flags($f,DrawCont) -command [list Update $f]
	$om add checkbutton -label "Antialiased lines Shift+A" -variable Flags($f,AntialiasedLines) -command [list Update $f]
	$om add checkbutton -label "Preserve aspect          " -variable Flags($f,PreserveAspect) -command [list Update $f]
	$om add checkbutton -label "Mapping                M " -variable Flags($f,Mapping) -command [list Update $f]
	$om add checkbutton -label "Draw color bar           " -variable Flags($f,ColorBar) -command [list Update $f]
	$om add checkbutton -label "Show title               " -variable Flags($f,Title) -command [list Update $f]
	$om add checkbutton -label "White background         " -variable Flags($f,WhiteBackground) -command [list Update $f]
	$om add command     -label "Transparency ...   Ctrl+T" -command [list SetTransparency $f]
	$om add command     -label "Field line stop radius.. " -command [list SetFieldLineStopRadius $f]
	$om add separator
	$om add checkbutton -label "Draw isosurface          " -variable Flags($f,DrawIsosurf) -command [list Update $f]
	$om add command     -label "Set iso-values ...       " -command [list SetIsoValue $f]
	$om add checkbutton -label "Smooth isosurface        " -variable Flags($f,SmoothIsosurfaces) -command [list Update $f]
	$om add checkbutton -label "Shiny isosurface         " -variable Flags($f,ShinyIsosurfaces) -command [list Update $f]
	$om add separator
	$om add checkbutton -label "Render volumetric        " -variable Flags($f,DrawVolumetric) -command [list Update $f]
#	$om add checkbutton -label "Antialiased pts in volum." -variable Flags($f,VolumetricAntialiasing) -command [list Update $f]
	$om add command     -label "Volumetric parameters .. " -command [list VolumetricParameters $f]

	# View menu	
	menubutton $f.menubar.view -text View -menu $f.menubar.view.menu -font $menufont
	set wm [menu $f.menubar.view.menu -font $menufont]
	$wm add checkbutton -label "3D view               3 " -variable Flags($f,view3D) -command [list Update $f]
	$wm add command     -label "Full window       Ctrl+F" -command [list FullWindow $f]
	$wm add command     -label "Full zoom             F " -command [list FullZoom $f]
	$wm add command     -label "Previous zoom         P " -command [list PreviousZoom $f]
	$wm add command     -label "Zoom in               < " -command [list Change3DViewpoint $f [expr 1/$ZoomInFactor] 0 0]
	$wm add command     -label "Zoom out         Shift+<" -command [list Change3DViewpoint $f $ZoomInFactor 0 0]
	$wm add separator
	$wm add command     -label "Origin viewpoint      O " -command [list Default3DViewpoint $f]
	$wm add command     -label "Center viewpoint      C " -command [list Default3DViewpoint $f]
	$wm add command     -label "X view                X " -command [list Xview $f]
	$wm add command     -label "Negative X view  Shift+X" -command [list Set3DAngles $f 90 0]
	$wm add command     -label "Y view                Y " -command [list Yview $f]
	$wm add command     -label "Negative Y view  Shift+Y" -command [list Set3DAngles $f 90 90]
	$wm add command     -label "Z view                Z " -command [list Zview $f]
	$wm add command     -label "Negative Z view  Shift+Z" -command [list Set3DAngles $f 179.99 180]
	$wm add separator
	$wm add command     -label "3D viewpoint ...  Ctrl+3" -command [list Edit3DViewpoint $f]
	$wm add command     -label "Define box ...    Ctrl+B" -command [list DefineBox $f]

	bind $f <Key-3> [list ToggleFlag $f view3D]
	bind $f <Key-i> [list ToggleFlag $f LinearInterpolation]
	bind $f <Key-g> [list ToggleFlag $f DrawGrid]
	bind $f <Key-h> [list ToggleFlag $f HideGhosts]
	bind $f <Key-C> [list ToggleFlag $f DrawCont]
	bind $f <Key-m> [list ToggleFlag $f Mapping]
	bind $f <Key-A> [list ToggleFlag $f AntialiasedLines]
	bind $f <Control-f> [list FullWindow $f]
	bind $f <Key-f> [list FullZoom $f]
	bind $f <Key-p> [list PreviousZoom $f]
	bind $f <Key-x> [list Xview $f]
	bind $f <Key-y> [list Yview $f]
	bind $f <Key-z> [list Zview $f]
	bind $f <Control-Key-b> [list DefineBox $f]
	bind $f <Control-Key-t> [list SetTransparency $f]

	# Palettes menu
	menubutton $f.menubar.palettes -text Palettes -menu $f.menubar.palettes.menu -font $menufont
	set pm [menu $f.menubar.palettes.menu -font $menufont]
	$pm add radiobutton -label Grayscale -variable Palette($f) -value Grayscale \
		-command "[list $f.win hcsetpalette {}]; [list Update $f]"
	$pm add command -label {Invert palette} -command "[list $f.win hcpalette invert]; [list Update $f]"
	$pm add command -label {Band palette}   -command "[list $f.win hcpalette band]; [list Update $f]"
	$pm add command -label {Dim palette}    -command "[list $f.win hcpalette dim]; [list Update $f]"
	if [file isdirectory $PaletteDir] {
		foreach palfile [glob -nocomplain $PaletteDir/*.raw] {
			$pm add radiobutton -label [file root [file tail $palfile]] -variable Palette($f) -value $palfile \
				-command "[list $f.win hcsetpalette $palfile]; [list Update $f]"
		}
	}
	bind $f <Key-d> "[list $f.win hcpalette dim]; [list Update $f]"

	# Variables menu
	menubutton $f.menubar.vars -text Variables -menu $f.menubar.vars.menu -font $menufont
	set vm [menu $f.menubar.vars.menu -font $menufont]
	VariablesMenu $f $vm Params($f,var) [list NewVariable $f] 1
	bind $f <Key-l> [list ToggleFlag $f Logarithmic]
	bind $f <Control-Key-m> [list SetMinMax $f]

	# Help menu
	menubutton $f.menubar.help -text Help -menu $f.menubar.help.menu -font $menufont
	menu $f.menubar.help.menu -font $menufont
	$f.menubar.help.menu add command -label "About hcvis   Control-H" -command "abouthcvis"
	$f.menubar.help.menu add command -label "File info        Ctrl+I" -command [list GetInfo $f]
	$f.menubar.help.menu add command -label "Variables descriptions " -command [list VariableInfo $f]
	bind $f <Control-Key-h> {abouthcvis}

	pack $f.menubar.file $f.menubar.options $f.menubar.view $f.menubar.vars $f.menubar.palettes -side left
	pack $f.menubar.help -side right
}

proc OpenNewWindow {win {noupdate 0}} {
	global Palette
	toplevel $win
	SetupMenus $win
#   togl $win.win -width 400 -height 400 -rgba true -alpha true -double true -depth true -ident $win
    togl $win.win -width 400 -height 400 -rgba true -double true -depth true -ident $win
	$win.win hcsetpalette $Palette($win)
	if {!$noupdate} {Update $win}
    pack $win.win -side left -padx 0 -pady 0 -fill both -expand true
}

proc NewWindow {oldwin {isslave 0}} {
	global NWinsOpen WindowList winID FileNameIndex Flags Params FlagList ParamList Palette Slave
	set win .win2D_$winID
	if {$isslave} {set Slave($win,$oldwin) 1}
	# Copy FileNameIndex, Flags, Params and Palette from oldwin to win
	set FileNameIndex($win) $FileNameIndex($oldwin)
	foreach flag $FlagList {
		set Flags($win,$flag) $Flags($oldwin,$flag)
	}
	foreach param $ParamList {
		set Params($win,$param) $Params($oldwin,$param)
	}
	set Params($win,VolumetricNpointsMillions) [expr 1e-6 * $Params($win,VolumetricNpoints)]
	set Palette($win) $Palette($oldwin)
	incr NWinsOpen
	lappend WindowList $win
	OpenNewWindow $win
	incr winID
	return $win
}

proc NewSlaveWindow {masterwin} {
	NewWindow $masterwin 1
}

# Execution starts here!
#puts stdout "FileNameList = $FileNameList"
#puts stdout "threeD = $threeD"
if {[info tclversion] >= 8.0} {
	if {[llength [info commands blt::graph]] > 0} {
		namespace import blt::*
	}
}
set NFileNames [llength $FileNameList]
wm withdraw .
if {[llength $FileNameList] == 0} {
	set newfn [LoadFile]
	if {[string length $newfn] == 0} {
		exit 0
	}
}
bind all <Control-q> exit
set winID 1
if $threeD {
	set Flags(default,view3D) 1
}
NewWindow default
# Avoid displaying a prompt by waiting here until main window is closed
# (which does never happen since main window was 'withdrawn' above)
tkwait window .
