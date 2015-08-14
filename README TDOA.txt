This MATLAB script uses the Time-Difference-of-Arrival (TDOA) technique to determine the location of an isotropic point source. 
Four sensors were used in a square geometric configuration. A cartesian co-ordinate system was used and the origin was centered
on the 1st sensor.

---------------------------------------------------------
|							|
|					SOURCE		|
|							|
|							|
|							|
|	4	1					|
|	X	X					|
|							|
|							|
|							|
|	X	X					|
|	3	2					|
|							|
---------------------------------------------------------

Since all sensors are located in the same plane, an ambiguity arises when determining the x-coordinate of the location of the source.

Imagine that the source is infinitesimally small emanating SINE waves in all directions. If we just consider the crests of the SINE
waves then we can draw concentric circles centered on the source. These circles will increase in radius as time increases. They will
be spaced equally apart by a distance exactly equal to one wavelength. 



