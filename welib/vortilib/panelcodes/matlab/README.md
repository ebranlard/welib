# 2D panels codes

More recent implementations:

- vortexpoints.m : 2D vortex points satisfying no *tangential* velocity inside the body. Inspired by Lewis's method. Similar to python implementation. The code is not written as a function at this point. 

## Older implementations:
Note, these are written as functions, but can also be run as standalone scripts if no arguments are provided, in which case a VonDeVooren airfoil is used.

In expected order of accuracy or "readiness": 

 - DoubletSource: Constant strength doublet that are solved for, constant sources prescribed based on freestream, Neuman BC (velocity). Probably the most accurate.
 - DoubletSourcePot: Use the velocity potential formulation. Inspired by  Katz Plotkin Program 8.
 - Doublet: Constant strength doublet, Neuman BC (velocity). Inspired from Katz Plotkin Program 3, should reproduce Figure 11.22
 - DoubletSource\_ManyOptions: Same as above with more trial and error options, for debugging.
 - DoubletSource: Constant strength sources that are solved for, constant doublet prescribed based on curvilinear coordinates, Neuman BC (velocity). Seems not as accurate.


Old filenames:

- fVortexCode2DDoublet.m
- fVortexCode2DDoubletSource.m
- etc.
