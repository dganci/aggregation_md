# ---------------------------------------------------------------------------
# vmd_martini_view.tcl - Martini 3 multi-chain visualization for VMD
#
#   source /data/tools/vmd_martini_view.tcl      (or paste into the TkConsole)
#
# Works on the protein-only trajectory (vis_prot_ref.gro + vis_prot_nojump.xtc)
# and on a full-system one; solvent/ion representations are only created if
# those beads are actually present.
# ---------------------------------------------------------------------------

# --- palette: 5 bright, well-separated hues, readable on a dark background ---
# Change these RGB triples to restyle everything; nothing else refers to colours.
set CG_PALETTE {
    {0.910 0.639 0.239}
    {0.310 0.639 0.890}
    {0.271 0.722 0.435}
    {0.878 0.404 0.608}
    {0.663 0.498 0.910}
}
# ColorIDs redefined by this script (VMD has only 33 slots; these are the
# least-used ones). Chain A..E are then mapped onto them.
set CG_SLOTS {23 24 25 26 27}

# Minimal "-flag value" parser: takes the caller's args and a dict of defaults,
# returns a dict. Unknown flags are an error rather than silently ignored.
proc cg_parse {argv defaults} {
    array set o $defaults
    if {[llength $argv] % 2} { error "cg: option \"[lindex $argv end]\" has no value" }
    foreach {k v} $argv {
        if {![string match "-*" $k]} { error "cg: expected -flag, got \"$k\"" }
        set key [string range $k 1 end]
        if {![info exists o($key)]} {
            error "cg: unknown option $k (valid: -[join [lsort [array names o]] { -}])"
        }
        set o($key) $v
    }
    return [array get o]
}

# cg_setup ?-mol id? ?-style tube|licorice|vdw? ?-smooth N?
#          ?-quality high|fast? ?-bg dark|light? ?-box on|off? ?-zoom Z?
proc cg_setup {args} {
    global CG_PALETTE CG_SLOTS
    # tolerate the old positional form cg_setup <molid>
    if {[llength $args] == 1 && ![string match "-*" [lindex $args 0]]} {
        set args [list -mol [lindex $args 0]]
    }
    array set o [cg_parse $args {mol top style licorice smooth 5 \
                                 quality high bg dark box off zoom 1.7}]
    set mol $o(mol)
    if {[llength [molinfo list]] == 0} {
        puts "cg_setup: no molecule loaded. Load the structure first, e.g."
        puts "   cg_load vis_prot_ref.gro vis_prot_nojump.xtc"
        return
    }
    if {$mol eq "top"} { set mol [molinfo top] }

    # ---- wipe existing representations ------------------------------------
    while {[molinfo $mol get numreps] > 0} { mol delrep 0 $mol }

    # ---- inventory ---------------------------------------------------------
    set selP [atomselect $mol "not resname W ION"]
    set selW [atomselect $mol "resname W"]
    set selI [atomselect $mol "resname ION"]
    set selB [atomselect $mol "name BB"]
    set nP [$selP num] ; set nW [$selW num] ; set nI [$selI num] ; set nB [$selB num]

    # chains: resid restarts at 1 in every protomer, so nchains = nBB / max(resid)
    set nchain 0 ; set perchain 0
    if {$nB > 0} {
        set rmax 0
        foreach r [$selB get resid] { if {$r > $rmax} { set rmax $r } }
        set nchain [expr {$nB / $rmax}]
        set perchain [expr {$nP / $nchain}]
    }
    puts [format "protein %d   water %d   ions %d   chains %d (%d beads each)" \
              $nP $nW $nI $nchain $perchain]

    # ---- label the chains --------------------------------------------------
    # .gro carries no chain column, and resid 1..108 repeats in every protomer,
    # so "resid 5" would select five different residues. Index ranges are the
    # only unambiguous handle; write chain/segname once and use those after.
    set letters {A B C D E F G H I J K L}
    for {set c 0} {$c < $nchain} {incr c} {
        set lo [expr {$c * $perchain}]
        set hi [expr {$lo + $perchain - 1}]
        set s [atomselect $mol "index $lo to $hi"]
        $s set chain [lindex $letters $c]
        $s set segname "P[expr {$c + 1}]"
        $s delete
    }

    # ---- physical bead radii ----------------------------------------------
    # VMD guesses ~1.5 A from the atom name; a Martini 3 regular bead is
    # sigma = 0.47 nm, i.e. r = 2.35 A. Without this every bead is too small
    # and the packing looks far looser than it is.
    set all [atomselect $mol all]
    $all set radius 2.35
    $all delete
    if {$nI > 0} { $selI set radius 1.90 }
    if {$nW > 0} { $selW set radius 2.35 }

    # ---- real bonds --------------------------------------------------------
    # VMD guesses bonds from covalent radii: at 3.5 A a BB-BB pair is far
    # outside the cutoff, so a CG system loads with no bonds at all and every
    # bond-based representation renders empty. Build the connectivity here.
    cg_bonds_from_layout $mol

    # ---- colours -----------------------------------------------------------
    for {set c 0} {$c < $nchain} {incr c} {
        set slot [lindex $CG_SLOTS [expr {$c % [llength $CG_SLOTS]}]]
        set rgb  [lindex $CG_PALETTE [expr {$c % [llength $CG_PALETTE]}]]
        color change rgb $slot [lindex $rgb 0] [lindex $rgb 1] [lindex $rgb 2]
        color Chain [lindex $letters $c] [lindex [colorinfo colors] $slot]
    }

    # ---- materials ---------------------------------------------------------
    catch {material delete CGSide}
    if {[catch {material add CGSide copy AOChalky}]} { material add CGSide }
    material change opacity CGSide 0.55

    # ---- representations ---------------------------------------------------
    # -quality fast drops the tessellation of every rep; with tens of thousands
    # of solvent beads this is what actually moves the frame rate, far more than
    # the sphere sizes.
    set res [expr {$o(quality) eq "fast" ? 10 : 24}]
    set rep 0

    # backbone: sticks along the BB-BB bonds built above
    mol selection "name BB"
    mol representation Licorice 1.1 $res $res
    mol color Chain
    mol material AOChalky
    mol addrep $mol ; incr rep

    # side chains: physically sized, semi-transparent so the trace stays visible
    mol selection "name SC1 SC2 SC3 SC4"
    mol representation VDW 1.0 $res
    mol color Chain
    mol material CGSide
    mol addrep $mol ; incr rep

    if {$nW > 0 || $nI > 0} {
        cg_solvent on -mol $mol
        set rep [molinfo $mol get numreps]
    }

    # ---- trajectory smoothing ---------------------------------------------
    # One frame per 100 ps of Martini time is coarse; a 5-frame running mean
    # removes the visual stutter. Display filter only - the data is untouched.
    for {set r 0} {$r < $rep} {incr r} { mol smoothrep $mol $r $o(smooth) }

    $selP delete ; $selW delete ; $selI delete ; $selB delete
    cg_display
    cg_bg $o(bg)
    if {$o(style) ne "licorice"} { cg_style $o(style) }
    cg_quality $o(quality)
    if {$o(box) eq "on"} { cg_box on }
    cg_view $o(zoom)
    puts "reps: $rep  style=$o(style) smooth=$o(smooth) quality=$o(quality)"
    puts "cg_speed X ?-skip N? | cg_style S | cg_smooth N | cg_quality high|fast | cg_view Z | cg_bg dark|light | cg_box on|off"
}

# Bonds from the Martini bead layout: beads of one residue are contiguous and
# ordered BB, SC1, SC2, ... so BB->BB along the chain, BB->SC1, SCn->SCn+1.
# Ring side chains (PHE/TYR/TRP/HIS) lose their closing bond - cosmetic only.
proc cg_bonds_from_layout {mol} {
    set all [atomselect $mol all]
    set names [$all get name]
    set chains [$all get chain]
    set resids [$all get resid]
    $all delete

    set bonds {}
    set prevBB -1 ; set prevChain "" ; set lastSC -1 ; set curBB -1
    set n [llength $names]
    for {set i 0} {$i < $n} {incr i} {
        set nm [lindex $names $i]
        set ch [lindex $chains $i]
        if {$nm eq "BB"} {
            if {$prevBB >= 0 && $ch eq $prevChain} { lappend bonds [list $prevBB $i] }
            set prevBB $i ; set prevChain $ch ; set curBB $i ; set lastSC -1
        } elseif {[string match "SC*" $nm]} {
            if {$nm eq "SC1"} {
                if {$curBB >= 0} { lappend bonds [list $curBB $i] }
            } elseif {$lastSC >= 0} {
                lappend bonds [list $lastSC $i]
            }
            set lastSC $i
        }
    }
    # Write the table with setbonds (core VMD): topotools is a plugin and is
    # not always loaded, and "topo addbond" rebuilds the whole table per call.
    for {set i 0} {$i < $n} {incr i} { set bl($i) {} }
    foreach b $bonds {
        set a [lindex $b 0] ; set c [lindex $b 1]
        lappend bl($a) $c ; lappend bl($c) $a
    }
    set out {}
    for {set i 0} {$i < $n} {incr i} { lappend out $bl($i) }
    set all [atomselect $mol all]
    $all setbonds $out
    $all delete
    puts "bonds: [llength $bonds]"
}

proc cg_display {} {
    display projection Orthographic
    display depthcue off
    display nearclip set 0.01
    axes location Off
    # ColorID 2 (gray) is redefined as the background; no representation here
    # uses it, so nothing else changes appearance.
    color change rgb gray 0.13 0.14 0.16
    color Display Background gray
    catch {display shadows on}
    catch {display ambientocclusion on}
    catch {display aoambient 0.85}
    catch {display aodirect 0.35}
    catch {display antialias on}
}

# Fit the current frame in the window. resetview alone leaves the assembly
# small once the chains have diffused apart.
proc cg_view {{zoom 1.7}} {
    display resetview
    scale by $zoom
}

proc cg_bg {{mode dark}} {
    if {$mode eq "light"} {
        color change rgb gray 0.96 0.96 0.96
    } else {
        color change rgb gray 0.13 0.14 0.16
    }
    color Display Background gray
}

# Load structure + trajectory into one molecule and configure in a single call.
# "vmd file.gro file.xtc -e this.tcl" does NOT work: VMD runs -e scripts before
# it loads the files named on the command line, so the script finds nothing.
# cg_load <struct> ?<traj>? ?-step N? ?-first N? ?-last N? plus any cg_setup flag.
# -step is the cheapest speed-up there is: loading every Nth frame cuts memory
# and makes playback N times faster in wall-clock terms without touching the
# animation rate.
proc cg_load {struct {traj ""} args} {
    array set o [cg_parse $args {step 1 first 0 last -1 style licorice smooth 5 \
                                 quality high bg dark box off zoom 1.7}]
    set mol [mol new $struct waitfor all]
    if {$traj ne ""} {
        mol addfile $traj molid $mol first $o(first) last $o(last) \
            step $o(step) waitfor all
    }
    puts "loaded: [molinfo $mol get numframes] frames (step $o(step))"
    cg_setup -mol $mol -style $o(style) -smooth $o(smooth) -quality $o(quality) \
             -bg $o(bg) -box $o(box) -zoom $o(zoom)
    return $mol
}

# Solvent as a thin skin on the proteins rather than a box full of fog.
#
#   cg_solvent on ?-shell A? ?-mode dots|points? ?-opacity X? ?-ions all|near|off?
#   cg_solvent off
#
# The shell is a distance selection recomputed every frame (mol selupdate),
# so it follows the proteins instead of freezing at frame 0. Note that "within"
# is plain Cartesian: on a -pbc nojump trajectory a water bead that is
# physically adjacent may sit a box vector away in the file, so the shell comes
# out patchy. Use a wrapped/centred trajectory (branch B) for solvent views.
proc cg_solvent {{state on} args} {
    array set o [cg_parse $args {mol top shell 8 mode dots size 0.40 \
                                 opacity 0.50 ions near \
                                 ionsize 0.65 ionopacity 0.85}]
    if {$o(mol) eq "top"} { set mol [molinfo top] } else { set mol $o(mol) }

    # drop any existing solvent/ion rep, whatever its index
    for {set r [expr {[molinfo $mol get numreps] - 1}]} {$r >= 0} {incr r -1} {
        set sel [lindex [molinfo $mol get "{selection $r}"] 0]
        if {[string match "*resname W*" $sel] || [string match "*resname ION*" $sel]} {
            mol delrep $r $mol
        }
    }
    if {$state eq "off"} { puts "solvent off" ; return }

    catch {material delete CGSolvent}
    material add CGSolvent
    material change opacity   CGSolvent $o(opacity)
    material change ambient   CGSolvent 0.05
    material change diffuse   CGSolvent 0.65
    material change specular  CGSolvent 0.00
    material change shininess CGSolvent 0.00

    catch {material delete CGIon}
    material add CGIon
    material change opacity   CGIon $o(ionopacity)
    material change ambient   CGIon 0.10
    material change diffuse   CGIon 0.75
    material change specular  CGIon 0.15
    material change shininess CGIon 0.30

    # Water gets its own colour slot: stock "silver" is too close to the dark
    # background to read at low opacity. Slot 22 is not used by anything else
    # here (chains take 23-27, the background takes 2).
    color change rgb 22 0.60 0.78 0.90

    set prot "not resname W ION"
    set nW [[atomselect $mol "resname W"] num]
    set nI [[atomselect $mol "resname ION"] num]

    if {$nW > 0} {
        mol selection "resname W and within $o(shell) of ($prot)"
        switch -- $o(mode) {
            points { mol representation Points 1.0 }
            dots   { mol representation VDW $o(size) 8 }
            default { error "cg_solvent: -mode must be dots or points" }
        }
        mol color ColorID 22
        mol material CGSolvent
        mol addrep $mol
        mol selupdate [expr {[molinfo $mol get numreps] - 1}] $mol on
    }
    if {$nI > 0 && $o(ions) ne "off"} {
        # ions carry real information (they pair with charged side chains), so
        # they stay a step more visible than the water - but only a step
        if {$o(ions) eq "near"} {
            mol selection "resname ION and within $o(shell) of ($prot)"
        } else {
            mol selection "resname ION"
        }
        mol representation VDW $o(ionsize) 12
        mol color Name
        mol material CGIon
        mol addrep $mol
        mol selupdate [expr {[molinfo $mol get numreps] - 1}] $mol on
    }
    puts "solvent: shell $o(shell) A, $o(mode) size $o(size) opacity $o(opacity) | ions $o(ions) size $o(ionsize) opacity $o(ionopacity)"
}

# Playback rate. VMD's speed runs 0..1, where 1 is "no delay between frames";
# once you are at 1 the only way to go faster is to show fewer frames, which is
# what -skip does (display every Nth frame without reloading anything).
proc cg_speed {{speed 1.0} args} {
    array set o [cg_parse $args {skip "" style ""}]
    animate speed $speed
    if {$o(skip)  ne ""} { animate skip  $o(skip) }
    if {$o(style) ne ""} { animate style $o(style) }
    puts "speed=[animate speed]  skip=[animate skip]  style=[animate style]"
}

# Render an image sequence for a movie.
#
#   cg_movie <outdir> ?-first N? ?-last N? ?-step N? ?-width W? ?-height H?
#            ?-renderer TachyonInternal|snapshot?
#
# TachyonInternal ray-traces off-screen: it honours shadows and ambient
# occlusion and works even if the VMD window is covered, but costs a couple of
# seconds per frame. "snapshot" grabs the OpenGL window instead - instant, but
# the window must be fully visible and unobscured or you capture whatever is on
# top of it.
#
# VMD writes the frames only; there is no encoder inside VMD, so the last step
# is ffmpeg (or QuickTime: File > Open Image Sequence). The command is printed
# at the end.
proc cg_movie {dir args} {
    array set o [cg_parse $args {first 0 last -1 step 1 width 1280 height 960 \
                                 renderer TachyonInternal fps 25}]
    set mol [molinfo top]
    set n [molinfo $mol get numframes]
    if {$o(last) < 0} { set o(last) [expr {$n - 1}] }
    file mkdir $dir
    display resize $o(width) $o(height)
    set i 0
    for {set f $o(first)} {$f <= $o(last)} {incr f $o(step)} {
        animate goto $f
        display update
        render $o(renderer) [format "%s/frame.%05d.tga" $dir $i]
        if {$i % 10 == 0} { puts "  frame $i (traj frame $f)" }
        incr i
    }
    puts "wrote $i frames to $dir"
    puts "encode with:"
    puts "  ffmpeg -framerate $o(fps) -i $dir/frame.%05d.tga \\"
    puts "         -c:v libx264 -pix_fmt yuv420p -crf 18 $dir/movie.mp4"
}

# Frame rate vs looks. "fast" also turns off ambient occlusion and shadows,
# which are the expensive part of the pretty mode.
proc cg_quality {{mode high}} {
    set mol [molinfo top]
    if {$mode eq "fast"} {
        catch {display shadows off}
        catch {display ambientocclusion off}
        set res 10
    } else {
        catch {display shadows on}
        catch {display ambientocclusion on}
        set res 24
    }
    for {set r 0} {$r < [molinfo $mol get numreps]} {incr r} {
        set st [lindex [molinfo $mol get "{rep $r}"] 0]
        set kind [lindex $st 0]
        switch -- $kind {
            Licorice { mol modstyle $r $mol Licorice [lindex $st 1] $res $res }
            VDW      { mol modstyle $r $mol VDW [lindex $st 1] $res }
        }
    }
    puts "quality = $mode (resolution $res)"
}

# ---- switches --------------------------------------------------------------
proc cg_style {what} {
    set mol [molinfo top]
    switch -- $what {
        tube     { mol modstyle 0 $mol Licorice 1.6 24 24
                   mol modstyle 1 $mol VDW 0.45 20 }
        licorice { mol modstyle 0 $mol Licorice 1.1 24 24
                   mol modstyle 1 $mol VDW 1.0 24 }
        vdw      { mol modstyle 0 $mol VDW 1.0 24
                   mol modstyle 1 $mol VDW 1.0 24 }
        default  { puts "cg_style tube|licorice|vdw" }
    }
}

proc cg_smooth {n} {
    set mol [molinfo top]
    for {set r 0} {$r < [molinfo $mol get numreps]} {incr r} { mol smoothrep $mol $r $n }
    puts "smoothing window = $n"
}

# The cell is drawn from the frame's box vectors. On a -pbc nojump trajectory
# the chains legitimately wander outside it, so the box is informative for a
# wrapped/centred file and misleading for a continuous one. Off by default.
proc cg_box {{on off}} {
    if {[catch {package require pbctools} e]} { puts "pbctools unavailable: $e" ; return }
    if {$on eq "on"} { pbc box -on -center bb -color gray -width 1 -style dashed } else { pbc box -off }
}

# Auto-configure if something is already loaded; otherwise wait for cg_load.
if {[llength [molinfo list]] > 0} {
    cg_setup
} else {
    puts "vmd_martini_view loaded.  Next:  cg_load <struct.gro> ?<traj.xtc>?"
}
