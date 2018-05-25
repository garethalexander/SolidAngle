# SolidAngle
-----------------------------------------------------------------------------

Computing the Solid Angle Function About a Link

Based on the article "Maxwell's Theory of Solid and the Construction of Knotted Vortex Field", Binysh J. and Alexander G.P (Ref. 1)

-----------------------------------------------------------------------------

Authors  : Jack Binysh, Gareth Alexander
Email    : j.binysh@warwick.ac.uk, g.p.alexander@warwick.ac.uk
Time     : April 2018

If you find a bug, or have ideas for improvements, please email j.binysh@warwick.ac.uk 

-----------------------------------------------------------------------------

This is a code to compute the solid angle function of a link from a user provided specification of its embedding in space. It implements eq. (3) of Ref. 1 to do so. 

If you use this code, please cite:
Binysh J., Alexander G.P, "Maxwell's Theory of Solid Angle and the Construction of Knotted Vortex Fields"


------------------------------------------------------------------------------
SOME REMARKS ON HOW THE CODE WORKS:

- The code is written to be read! Comments on function purposes etc are found in the header for each file, implementation details are in comments in the .cpp files. Start at SolidAngle.cpp and flick through the function calls in main().

- Bear in mind this code expects "friendly" input - it doesn't do much serious input cleaning/checking. But it’s pretty simple, so use it as a white box! 

- The code takes an ordered list of points in space, and outputs a 3D cubic grid of solid angle values around it. The knot embedding should be thought of as lying in "physical space" with some arbitrary units. It has some discretization on this lengthscale. Then there is another lengthscale, that of the grid discretization. These two scales don't have to be linked, though they should be if you want good output (see the Ref. 1). Broadly it is the user’s responsibility to get the scaling right. The "physical" box size is set by N h  where N is the number of gridpoints, h is the gridspacing (see section INPUT).

- The code does not implement the cylindrical mesh around the knot discussed in the section "Remarks on Numerical Implementation" of Ref. 1 - everything is computed on a simple cubic grid. 

- The code does not assume equispaced points of discretization (though it has been tested mainly when this has been approximately the case). 

- Rather than picking n_infty truly at random as in Ref. 1, the code cycles through a predefined list of choices - this means it gives identical answers between runs and across thread number.

------------------------------------------------------------------------------
COMPILATION:

to compile, type:

g++ -o SolidAngle -O3 SolidAngle.cpp Geometry.cpp InputOutput.cpp Constants.cpp  

The code is parallelised using openmp. To take advantage of this, have openmp installed and type

g++ -o SolidAngle -O3 -fopenmp SolidAngle.cpp Geometry.cpp InputOutput.cpp Constants.cpp  

and then set the number of threads using export OMP_NUM_THREADS=<thread number> before running the executable. See any openmp resource for more details.


------------------------------------------------------------------------------
RUNNING

To run the program, type

./SolidAngle <my_knot_name>

where you give the name of your knot folder (see below for details) as a command line input.

------------------------------------------------------------------------------
INPUT

We have provided two example sets of input files inside the "Knots" folder. The first set is found in the folder "Unknot" which gives example files for a twisted unknot. 
The second set is found in the folder "Whitehead_Link" which gives example files for a Whitehead link.

As input, the code expects a folder, <my_knot_name>, inside another folder called Knots (so full path is Knots/<my_knot_name>). Inside THIS folder the code expects:

---

(1) A file called "input.txt" containing user defined settings for the code. I'll run through the option meanings, with reference to the Unknot example file: 

NumComponents=1
An integer. The number of components of the link. See above.

NumRefinements=0
An integer >= 0. As discussed above, units are arbitrary and it’s the user’s responsibility to ensure their knot discretization is sufficient to give good results over their chosen gridspacing scale. That said, if you turn this option on, your input knot will be refined to have new points placed linearly between existing points the specified number of times. Set this to 0 if you don't want it done. A convenience for quick-and-dirty work.  

h=0.5
A double. The Gridspacing.

Nx=201
Ny=201
Nz=201
Three integers giving the number of gridpoints along x,y,z. The box sides don’t have to all be the same.

Scaling=1
A bool (can be 0 or 1). As discussed above, units are arbitrary and it’s the user’s responsibility to ensure their knot discretization is sufficient to give good results over their chosen gridspacing scale. That said, if you turn this option on, your input knot will simply be scaled to some fraction of the given physical box size. A convenience for quick-and-dirty work.

ScaleProportionally=1
A bool (0 or 1). If turned off, the aspect ratio of the knot’s bounding box can change during scaling (if you made the box fractions different). If it’s on, the aspect ratio will be preserved - the "safest" of the 3 scalings will be used.:w

BoxFractionx=0.6
BoxFractiony=0.6
BoxFractionz=0.6
Three doubles between 0 and 1, specifying what fraction of the box size the bounding box of the input knot will be scaled to, along each dimension. They can differ.

Initialninftyx= 0
Initialninftyy= 0
Initialninftyz= 1
Three doubles. The 3 components of the vector n_infty described first in the section "A curve homotopy formula" of Ref 1. Unless you have a good reason not to, just set n_infty = (0,0,1). Make sure it’s a unit vector!

Threshold_ninf_dot_n=-0.95
A double between 0 and -1. The user defined threshold discussed in the section "Remarks on Numerical Implementation". It specifies when you consider the most negative value of n dot n_infty unacceptable (when your discretization starts missing the Lorentzian peaks near the surface of discontinuity). In practice, this is set somewhere between -0.9 and -0.8 --- just try -0.95, and if the answer has unacceptable artefacts, make it less negative.

---

(2) A file containing an ordered list (the ordering defining an orientation) of points also called <my_knot_name>.txt, with newlines between points and spaces separating the coordinates within points. If inputting a link (see the Whitehead Link example), provide each component in a separate file, the naming convention going <my_knot_name>_0.txt, <my_knot_name>_1.txt etc. Note the numbering is 0-based. 

------------------------------------------------------------------------------
OUTPUT 

Output goes inside the folder that you read in. 

The code outputs files in human (or grad student) readable .vtk format - the "legacy" format detailed at https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
These files may be read into e.g. ParaView https://www.paraview.org/ for visualisation.

It outputs 3 files

(1) A .vtk file containing the solid angle data on a 3D grid, in the range 0 to 4 pi, named <my_name>_0_4pi.vtk
(2) "" in the range -2pi to 2pi, named <my_name>_-2pi_2pi.vtk

The two are useful to load simultaneously into ParaView - it will have a little trouble getting the level set omega=0 in file (1), but no trouble in file (2).

(3) The input knot, scaled if scalings were set, named ParaViewBoundaryCurve_<my_name>.obj in .obj format http://paulbourke.net/dataformats/obj/ so it can be loaded into Paraview and visualised along with the solid angle function.

It also outputs progress to standard output, and error messages to standard error.

----------------------------------------------------------------------------
A BRIEF GUIDE TO POST-PROCESSING IN PARAVIEW

Here we briefly outline how the visuals in Ref.1 were constructed in ParaView.

(1) Load in one of the ParaView files, e.g. in the range 0_4pi. Unfortunately, if one simply visualises this function, ParaView will get confused at the cut, and create an unpleasant discontinuity there. We can get around this by filtering the spurious high gradient ParaView associates to this cut, so ...

(2) Apply the "Gradient" filter to the dataset. 

(3) Apply the "Calculator" filter to "Gradient". Set it to calculate sqrt(omegaGradient.omegaGradient), where "omegaGradient" can be found in the Vectors dropdown menu. 

(4) Apply the "Threshold" filter to "Calculator". The Scalar "Result" contains our values of |grad(omega)|. Reduce the maximum threshold value until the cut vanishes (this may be easier to see after the next step...).

(5) apply the "Contour" filter to "Threshold", and select contour by "omega". Also make sure to tick "Compute Scalars", so your contours are coloured by omega.

One may then fiddle with opacities etc. to get a nice picture. Two points about colour:

- The standard colour map is red to blue. One instead may want a colour wheel - "HSV" in the predefined favourites is such a wheel, or define your own.

- If taking a slice through data (e.g. the Slice through the Whitehead Link of Fig. 4 (d)) you don't need to do the above - however, do select a colour map which is circular, and make sure to turn off "interpolate scalars before mapping" to remove the cut artefacts.

ParaView will have trouble getting the omega = 0 (/4 pi) contour from the 0_4pi input file. To get around this, simply apply the above procedure to the -2pi_2pi file. The "ParaViewBoundaryCurve" output may also be imported, so one can see the curve clearly.

We emphasise briefly: the cut ParaView sees in e.g. the 0_4pi dataset has nothing to do with the "surface of discontinuity" discussed in the Ref 1. It’s just a consequence of what arbitrary range we decide to put our circle value function into. 

Happy visualising!

----------------------------------------------------------------------------
REFERENCES

(1) "Maxwell's Theory of Solid and the Construction of Knotted Vortex Field", Binysh J. and Alexander G.P 
