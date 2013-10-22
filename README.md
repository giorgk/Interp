Interp repo provides a class for interpolation of gridded data.

The code is my first attempt to write c++ in a dimension independed way. Supports 1 2 and 3D interpolation.
Besides gridded format elevations Z can vary in space, which is very convinient in porous media simulations.
This option is used only in 3D.
The code ignores missing data. In particular when all data are available the interpolation is linear, 
bilinear or trilinear. When there are gaps in data then an inverse distance weighting interpolation is used.
A similar scheme is also used for extrapolation. 
NOTE: The class is purposufully designed so that always returns a value. If you get a nan or some other 
unexpected results then please report it because its probably a bug.

Contents:
interpND.h is the main c++ header file.

Examples on how to use the class within a c++ can be found in Interpolation.cpp file

I'm also providing wrappers of this class for Matlab and octave. 
In particular both matlab and octave calls the same matlab function which identifies the system and runs 
the appropriate wrapper. Examples on how to call the interpolation class from Matlab/Octave can be found 
Tutorial_interp.m script. The first part of the script gives examples on how to prepare input files.
The second part gives examples on how to call the class form Matlab/Octave.

