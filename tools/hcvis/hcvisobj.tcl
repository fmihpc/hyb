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

# Object list is of the form: {type name flag {field1 val1 field2 val2 ...}}
# where flag is either active or inactive, type is Slice, Sphere or FieldLineBunch.

proc ToObject {obj} {
	# Transfer contents of obj list into global array Object.
	# This is called before calling hcdefobject, and also before calling
	# the object-specific dialog boxes.
	global Object
	array set Object {}
	array set Object [lindex $obj 3]
	set Object(type) [lindex $obj 0]
}

proc SetWindowCloseToMouse {win} {
	set oldgeom [winfo geometry $win]
	set maxw [winfo screenwidth $win]
	set maxh [winfo screenheight $win]
	set x [winfo pointerx $win]
	set y [winfo pointery $win]
	wm geometry $win +$x+$y
	update
	set w [winfo width $win]
	set h [winfo height $win]
	set remap 0
	if {$x+$w > $maxw} {set x [expr $maxw-$w]; set remap 1}
	if {$y+$h > $maxh} {set y [expr $maxh-$h]; set remap 1}
	if {$remap} {wm geometry $win +$x+$y}
}

proc EditObjMakeNewName {root} {
	# Make a new unique name, trying root1, root2,... until a nonexistent is found.
	# Stop after 100 trials. All object names in EditObj(list) are scanned.
	# Notice that this contains only a temporary copy of the current windows' object list.
	global EditObj
	for {set n 1} {$n < 100} {incr n} {
		set trial "${root}$n"
		set unique 1
		foreach obj $EditObj(list) {
			set name [lindex $obj 1]
			if {[string compare $trial $name]==0} {
				set unique 0
				break
			}
		}
		if {$unique} {return $trial}
	}
	return "Failed to make a new name"
}

proc EditObjects {win} {
	# The main object editor window
	global Params EditObj
	set fl $win.edit
	if {[winfo exists $fl]} {
		raise $fl
		return
	}
	toplevel $fl
	wm title $fl {Edit objects}
	frame $fl.area
	set lb $fl.area.list
	set width 20
	set EditObj(list) $Params($win,obj)
	set height [llength $EditObj(list)]
	if {$height < 5} {set height 5}
	if {$height > 20} {set height 20}
	listbox $lb -yscrollcommand [list $fl.area.sy set] -xscrollcommand [list $fl.area.sx set] -width $width -height $height
	scrollbar $fl.area.sy -orient vertical -command [list $lb yview]
	scrollbar $fl.area.sx -orient horizontal -command [list $lb xview]
	EditObjUpdate $win
	set bs [frame $fl.buttons]
	menubutton $bs.add -text Add -menu $bs.add.menu
	menu $bs.add.menu
	$bs.add.menu add command -label {X slice} \
		-command "DefaultObjSlice; set Object(SlicedDim) 0; [list EditObjSlice $win]; [list EditObjUpdate $win]"
	$bs.add.menu add command -label {Y slice} \
		-command "DefaultObjSlice; set Object(SlicedDim) 1; [list EditObjSlice $win]; [list EditObjUpdate $win]"
	$bs.add.menu add command -label {Z slice} \
		-command "DefaultObjSlice; set Object(SlicedDim) 2; [list EditObjSlice $win]; [list EditObjUpdate $win]"
	$bs.add.menu add command -label {Sphere} \
		-command "DefaultObjSphere; [list EditObjSphere $win]; [list EditObjUpdate $win]"
	$bs.add.menu add command -label {rhov field line bunch} \
		-command "DefaultObjFieldLineBunch; set Object(VectorField) rhov; \
				  [list EditObjFieldLineBunch $win]; [list EditObjUpdate $win]"
	$bs.add.menu add command -label {B field line bunch} \
		-command "DefaultObjFieldLineBunch; set Object(VectorField) B; \
				  [list EditObjFieldLineBunch $win]; [list EditObjUpdate $win]"
	$bs.add.menu add command -label {j field line bunch} \
		-command "DefaultObjFieldLineBunch; set Object(VectorField) j; \
				  [list EditObjFieldLineBunch $win]; [list EditObjUpdate $win]"
	$bs.add.menu add command -label {B0 field line bunch} \
		-command "DefaultObjFieldLineBunch; set Object(VectorField) B0; \
				  [list EditObjFieldLineBunch $win]; [list EditObjUpdate $win]"
	$bs.add.menu add command -label {B1 field line bunch} \
		-command "DefaultObjFieldLineBunch; set Object(VectorField) B1; \
				  [list EditObjFieldLineBunch $win]; [list EditObjUpdate $win]"
	$bs.add.menu add command -label {S (Poynting) field line bunch} \
		-command "DefaultObjFieldLineBunch; set Object(VectorField) S; \
				  [list EditObjFieldLineBunch $win]; [list EditObjUpdate $win]"
	button $bs.addsim -text {Add similar} -command "[list EditObjEdit $win 1]; [list EditObjUpdate $win]"
	button $bs.edit -text {Edit} -command [list EditObjEdit $win]
	button $bs.remove -text Remove -command "[list EditObjRemoveItem $win]; [list EditObjUpdate $win]"
	button $bs.select -text Select -command "[list EditObjSelectItem $win 1]; [list EditObjUpdate $win]"
	button $bs.deselect -text Deselect -command "[list EditObjSelectItem $win 0]; [list EditObjUpdate $win]"
	button $bs.ok -text OK -underline 0 -command "[list EditObjApply $win]; [list destroy $fl]"
	button $bs.apply -text Apply -underline 0 -command [list EditObjApply $win]
	button $bs.cancel -text Cancel -underline 0 -command [list destroy $fl]
	pack $fl.area -side left -fill both -expand true
	pack $fl.area.sy -side right -fill y
	pack $fl.area.sx -side bottom -fill x
	pack $lb -side left -fill both -expand true
	pack $bs.add $bs.addsim $bs.edit $bs.remove $bs.select $bs.deselect
	pack $bs.ok $bs.apply $bs.cancel -side bottom
	pack $bs -side left -fill y
	bind $fl <Return> "[list EditObjApply $win]; [list destroy $fl]"
	bind $fl <Double-1> "[list EditObjApply $win]; [list destroy $fl]"
	bind $fl <Control-c> [list destroy $fl]
	bind $fl <Control-w> [list destroy $fl]
	bind $fl <Control-a> [list EditObjApply $win]
	$lb selection set 0
	$lb see 0
}

proc EditObjApply {win} {
	# Write the object list from EditObj dialog box specific storage into
	# final resting place (Params array)
	global Params EditObj
	set Params($win,obj) $EditObj(list)
#	puts "EditObjApply: objectlist=$EditObj(list)"
	Update $win
}

proc EditObjRemoveItem {win} {
	# Remove the (first of the) currently selected item(s) from the EditObj(list)
	# Just modifies the EditObj object list, does not call EditObjApply.
	global EditObj
	set i [lindex [$win.edit.area.list curselection] 0]
	set EditObj(list) [lreplace $EditObj(list) $i $i]
}

proc EditObjSelectItem {win selectflag} {
	global EditObj
	set i [lindex [$win.edit.area.list curselection] 0]
	set obj [lindex $EditObj(list) $i]
	if {$selectflag} {set newflag active} else {set newflag inactive}
	set obj [lreplace $obj 2 2 $newflag]
	set EditObj(list) [lreplace $EditObj(list) $i $i $obj]
}

proc EditObjUpdate {win} {
	# Updates the listbox display to reflect changes in EditObj(list).
	# Call this when first doing EditObj dialog, and after modifying EditObj(list).
	global EditObj
	set fl $win.edit
	set lb $fl.area.list
	set selection [$lb curselection]
	if {[string length $selection] == 0} {set selection 0}
	$lb delete 0 end
	foreach obj $EditObj(list) {
		if {[string compare [lindex $obj 2] active]==0} {set star "*"} else {set star " "}
		$lb insert end "$star [lindex $obj 1]"
	}
	$lb selection set $selection
	$lb see $selection
}

proc EditObjEdit {win {createnew 0}} {
	# Edits the currently selected object
	# if createnew==1, creates a new object based on old, otherwise modifies old.
	# Branches to the correct object-specific routine through a switch statement.
	global EditObj EditObjSlice EditObjSphere
	set selected [$win.edit.area.list curselection]
	if {[llength $selected] == 0} {return}
#	set i [lindex $selected 0]
	set i $selected
	set obj [lindex $EditObj(list) $i]
	if {$createnew} {set objindex {}} else {set objindex $i}
	# Transfer obj into Object
#puts "EditObjEdit: i=$i, obj=$obj"
	ToObject $obj
	# Call either EditObjSlice, EditObjSphere or EditObjFieldLineBunch
	EditObj[lindex $obj 0] $win $objindex
}

proc EditObjCommonApply {objindex} {
	# Update EditObj(list) from fields of Object(). Callback from object-specific dialog boxes.
	# If objindex is {}, creates a new object, otherwise modifies old.
	# Calls EditObj${type}RootName, which expands to EditObjSliceRootName etc.,
	# this function must exist.
	global EditObj Object
	set type $Object(type)
	if {[string length $objindex] > 0} {
		# Modify old object
		set obj [lindex $EditObj(list) $objindex]
		set oldname [lindex $obj 1]
		set activestatus [lindex $obj 2]
		set EditObj(list) \
			[lreplace $EditObj(list) $objindex $objindex \
				[list $type $oldname $activestatus [array get Object]]]
	} else {
		# Create new object
		lappend EditObj(list) [list $type [EditObjMakeNewName [EditObj${type}RootName]] active [array get Object]]
	}
}

# -----------------------------------------------------------------------------
# The slice object dialog (EditObjSlice, EditObjSliceRootName, DefaultObjSlice)
# -----------------------------------------------------------------------------

proc EditObjSlice {win {objindex {}}} {
	# If objindex is {}, creates a new object, otherwise modifies old.
	global Object
	set w [toplevel $win.edit.slice]
	wm title $w "Edit Slice"
	grab $w
	set a [frame $w.area]
	radiobutton $a.x -text {X slice} -variable Object(SlicedDim) -value 0
	radiobutton $a.y -text {Y slice} -variable Object(SlicedDim) -value 1
	radiobutton $a.z -text {Z slice} -variable Object(SlicedDim) -value 2
	set a2 [frame $w.area2]
	label $a2.label -text {Coordinate value}
	entry $a2.entry -textvariable Object(Slice_Xval) -width 10
	set bs [frame $w.buttons]
	set ok [button $bs.ok -text OK -command "[list EditObjCommonApply $objindex]; [list destroy $w]" -underline 0]
	set cancel [button $bs.cancel -text Cancel -command "destroy $w"]
	pack $a -side top
	pack $a.x $a.y $a.z -side left
	pack $a2 -side top
	pack $a2.label $a2.entry -side left
	pack $bs -side bottom
	pack $cancel $ok -side left
	bind $w <Return> "[list EditObjCommonApply $objindex]; [list destroy $w]"
	SetWindowCloseToMouse $w
	tkwait window $w
	grab release $w
}

proc EditObjSliceRootName {} {
	global Object
	switch -exact -- $Object(SlicedDim) {
	0 {set xyz X}
	1 {set xyz Y}
	2 {set xyz Z}
	}
	return "${xyz}-Slice "
}

proc DefaultObjSlice {} {
	global Object
	array set Object {}
	set Object(SlicedDim) 0
	set Object(Slice_Xval) 0.0
	set Object(type) Slice
}

# ---------------------------------------------------------------------------------
# The sphere object dialog (EditObjSphere, EditObjSphereRootName, DefaultObjSphere)
# ---------------------------------------------------------------------------------

proc EditObjSphere {win {objindex {}}} {
	# If objindex is {}, creates a new object, otherwise modifies old.
	global Object
	set w [toplevel $win.edit.sphere]
	wm title $w "Edit Sphere"
	grab $w
	set a [frame $w.area]
	set width 7
	set r   [frame $a.r]
	set rl  [label $r.l -text Radius]
	set re  [entry $r.e -textvariable Object(Radius) -width $width]
	set x0  [frame $a.x0]
	set x0l [label $x0.l -text x0]
	set x0e [entry $x0.e -textvariable Object(x0) -width $width]
	set y0  [frame $a.y0]
	set y0l [label $y0.l -text y0]
	set y0e [entry $y0.e -textvariable Object(y0) -width $width]
	set z0  [frame $a.z0]
	set z0l [label $z0.l -text z0]
	set z0e [entry $z0.e -textvariable Object(z0) -width $width]
	set dt  [scale $a.dt -from 0.5 -to 20 -length 80 -resolution 0.5 -orient horizontal \
				-label {interp step in deg} -showvalue true -variable Object(DeltaTheta)]
	set dtgrid [scale $a.dtgrid -from 0.5 -to 20 -length 80 -resolution 0.5 -orient horizontal \
				-label {grid step in deg} -showvalue true -variable Object(DeltaThetaGrid)]
	set bs [frame $w.buttons]
	set ok [button $bs.ok -text OK -command "[list EditObjCommonApply $objindex]; [list destroy $w]" -underline 0]
	set cancel [button $bs.cancel -text Cancel -command "destroy $w"]
	pack $a -side top
	pack $r $x0 $y0 $z0 $dt $dtgrid -side top
	pack $rl $re -side left
	pack $x0l $x0e -side left
	pack $y0l $y0e -side left
	pack $z0l $z0e -side left
	pack $bs -side bottom
	pack $cancel $ok -side left
	bind $w <Return> "[list EditObjCommonApply $objindex]; [list destroy $w]"
	SetWindowCloseToMouse $w
	tkwait window $w
	grab release $w
}

proc EditObjSphereRootName {} {
	return "Sphere "
}

proc DefaultObjSphere {} {
	global Object
	array set Object {}
	set Object(Radius) 1
	set Object(x0) 0
	set Object(y0) 0
	set Object(z0) 0
	set Object(DeltaTheta) 5
	set Object(DeltaThetaGrid) 10
	set Object(type) Sphere
}

# -------------------------------------------------------------------------------------------------------------------
# The field line bunch object dialog (EditObjFieldLineBunch, EditObjFieldLineBunchRootName, DefaultObjFieldLineBunch)
# -------------------------------------------------------------------------------------------------------------------

proc EditObjFieldLineBunch {win {objindex {}}} {
	# If objindex is {}, creates a new object, otherwise modifies old.
	global Object
	set w [toplevel $win.edit.flbunch]
	wm title $w {Edit field line bunch}
	grab $w
	set r1r2 [frame $w.r1r2]
	set width 7
	set r1x [frame $r1r2.r1x]
	set r1xl [label $r1x.l -text r1x]
	set r1xe [entry $r1x.e -textvariable Object(r1x) -width $width]
	set r1y [frame $r1r2.r1y]
	set r1yl [label $r1y.l -text r1y]
	set r1ye [entry $r1y.e -textvariable Object(r1y) -width $width]
	set r1z [frame $r1r2.r1z]
	set r1zl [label $r1z.l -text r1z]
	set r1ze [entry $r1z.e -textvariable Object(r1z) -width $width]
	set sep1 [label $r1r2.sep1 -text {}]
	set r2x [frame $r1r2.r2x]
	set r2xl [label $r2x.l -text r2x]
	set r2xe [entry $r2x.e -textvariable Object(r2x) -width $width]
	set r2y [frame $r1r2.r2y]
	set r2yl [label $r2y.l -text r2y]
	set r2ye [entry $r2y.e -textvariable Object(r2y) -width $width]
	set r2z [frame $r1r2.r2z]
	set r2zl [label $r2z.l -text r2z]
	set r2ze [entry $r2z.e -textvariable Object(r2z) -width $width]
	set bs [frame $w.bs]
	# N scale
	set N  [scale $w.n -from 50 -to 1 -length 160 -resolution 1 -orient vertical -label Npts -showvalue true -variable Object(N)]
	# Direction
	set dir     [frame $bs.dir]
	set dirl    [label $dir.l -text {Direction:}]
	set pos     [radiobutton $dir.pos       -text Positive -variable Object(Direction) -value Positive]
	set neg     [radiobutton $dir.neg       -text Negative -variable Object(Direction) -value Negative]
	set both    [radiobutton $dir.both      -text Both     -variable Object(Direction) -value Both]
	# Distribution
	set distr   [frame $bs.distr]
	set distrl  [label $distr.l -text {Distribution:}]
	set uniform [radiobutton $distr.uniform -text Uniform -variable Object(Distribution) -value Uniform]
	set invsqrt [radiobutton $distr.invsqrt -text InvSqrt -variable Object(Distribution) -value InvSqrt]
	set inv     [radiobutton $distr.inv     -text Inv     -variable Object(Distribution) -value Inv]
	# Field
	set field   [frame $bs.field]
	set fieldl  [label $field.l -text {Vector field:}]
	set rhov    [radiobutton $field.rhov    -text rhov    -variable Object(VectorField)  -value rhov]
	set B       [radiobutton $field.b       -text B       -variable Object(VectorField)  -value B]
	set j       [radiobutton $field.j       -text j       -variable Object(VectorField)  -value j]
	set B0      [radiobutton $field.b0      -text B0      -variable Object(VectorField)  -value B0]
	set B1      [radiobutton $field.b1      -text B1      -variable Object(VectorField)  -value B1]
	set S       [radiobutton $field.s       -text S       -variable Object(VectorField)  -value S]
	set K       [radiobutton $field.k       -text K       -variable Object(VectorField)  -value K]
	# Loop threshold
	set lt      [frame $bs.lt]
	set ltl     [label $lt.l -text {Loop threshold:}]
	set lts     [scale $lt.s -from 0.0 -to 0.25 -resolution 0.01 -orient horizontal \
					-showvalue true -variable Object(LoopThreshold)]
	# Loop threshold type
	set ltt     [frame $bs.ltt]
	set lttl    [label $ltt.l -text {Threshold:}]
	set local   [radiobutton $ltt.local   -text {Local dx rel} -variable Object(LoopThresholdType) -value LocalSpacing]
	set basegr  [radiobutton $ltt.basegr  -text {Basegrid dx rel} -variable Object(LoopThresholdType) -value BasegridSpacing]
	set abso    [radiobutton $ltt.abso    -text {Absolute}      -variable Object(LoopThresholdType) -value Absolute]
	# Cancel, OK
	set cancelok [frame $bs.cancelok]
	set cancel [button $cancelok.cancel -text Cancel -command "destroy $w"]
	set ok [button $cancelok.ok -text OK -command "[list EditObjCommonApply $objindex]; [list destroy $w]" -underline 0]
	# Pack
	pack $r1r2 -side left
	pack $N -side left -ipadx 5
	pack $bs -side left
	pack $r1x $r1y $r1z $sep1 $r2x $r2y $r2z -side top
	pack $r1xl $r1xe -side left
	pack $r1yl $r1ye -side left
	pack $r1zl $r1ze -side left
	pack $r2xl $r2xe -side left
	pack $r2yl $r2ye -side left
	pack $r2zl $r2ze -side left
	pack $dir $distr $field $lt $ltt -side top
	pack $dirl $pos $neg $both -side left
	pack $distrl $uniform $invsqrt $inv -side left
	pack $fieldl $rhov $B $j $B0 $B1 $S $K -side left
	pack $ltl $lts -side left
	pack $lttl $local $basegr $abso -side left
	pack $cancelok -side bottom
	pack $cancel $ok -side right
	bind $w <Return> "[list EditObjCommonApply $objindex]; [list destroy $w]"
#	SetWindowCloseToMouse $w
	tkwait window $w
	grab release $w
}

proc EditObjFieldLineBunchRootName {} {
	return "FieldLineBunch "
}

proc DefaultObjFieldLineBunch {} {
	global Object
	array set Object {}
	set Object(r1x) 0.0
	set Object(r1y) 0.0
	set Object(r1z) 0.0
	set Object(r2x) 1.0
	set Object(r2y) 0.0
	set Object(r2z) 0.0
	set Object(N) 10
	set Object(Direction) Both
	set Object(Distribution) Uniform
	set Object(VectorField) rhov
	set Object(LoopThreshold) 0.05
	set Object(LoopThresholdType) LocalSpacing
	set Object(type) FieldLineBunch
}
