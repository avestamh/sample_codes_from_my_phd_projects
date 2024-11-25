# tcl script for protein unfolding movie

set psfdir     "/home/sadra/anal-traj-gfp-go-bias1-3e8/gfp-bulk/psfgen"
set dscdir     "/home/sadra/anal-traj-gfp-go-bias1-3e8/gfp-bulk/trajs/k1-rate001/pulling-go-gfp-rate0.001-k1-traj6"
set fname	"go-gfp"
set datfile     "smd1.1-229_1"
set psf		$psfdir/$fname.psf
set dcd         $dscdir/$datfile.dcd
#set cor		$psfdir/$fname.cor

#load the molecule
mol new       $psf
mol addfile   $dcd  waitfor all

# change VMD setting
color Display Background white
display depthcue off
display projection orthographic
axes location off

# Change the orientation of thr protein
# del the current rep
mol delrep 0 0
# creat a new rep
mol rep Tube 0.7 20
mol color ColorID 19
mol addrep 0

#creat a new rep for 10helix-helix1
mol selection "resid 9 to 11"
mol rep tube 0.9 20
mol color ColorID 17;#color for 10helix-helix1
mol addrep 0;#creat the rep

#creaaat a new rep for B1sheet
mol selection "resid 12 to 22"
mol rep tube 0.9 20
mol color ColorID 29 ;#color of B1
mol addrep 0 ;# creat the rep

# creat a new rep for B3
mol selection "resid 41 to 48"
mol rep tube 0.9 20
mol color ColorID 23
mol addrep 0

# creat a new rep for B10
mol selection "resid 199 to 208"
mol rep tube 0.9 20
mol color ColorID 27
mol addrep 0

################

#change the resoloution

aasamples TachyonInternal 40


# CHANGE THE ORIENTATION OF PROTEIN
#rotate z by 90 ;# I DID NOT ASK THAT DAMN IRANIAN
scale by 0.5



proc renderme { } { 
 set nr_steps [molinfo 0 get numframes ] 
 puts "number of steps: $nr_steps" 
 set file [molinfo 0 get name] 
 for {set t 0} {$t <= 1600} {incr t 10} { 
 animate goto $t 
 display update
# display resetview
 puts "rendering frame: $t" 
render TachyonInternal "gfp-frame$t.tga"
} 
}    
