/*
	Author:  Joshua Simmons
	Started: September, 2015
	Status:  INCOMPLETE

	PINGER ASSUMED TO BE INTERMITTENT AT FIXED INTERVALS!!!

	Description: Uses Time-Difference-of-Arrival (TDOA) to determine the azimuth to
				 pulsing signal source underwater.

	Source Characteristics
		1. SINE wave
		2. 20 [kHz] to 30 [kHz]
		3. tOn = 1.3 [ms]
		4. PRT = 2 [s]

	3D Cartesian co-ordinate system with the origin centered in the middle of the sensor
	array. Sensor geometry is square shaped residing all in the XY-plane.

	Sensor layout

		  (Top View)
	---------------------
	|                   |
	|         1         |
	|                   |
	|     4       2     |
	|                   |
	|         3         |
	|                   |
	---------------------

	Coordinates
		chan1  = ( 0, d,0)		D = sqrt(2)*d
		chan2  = ( d, 0,0)
		chan3  = ( 0,-d,0)
		chan4  = (-d, 0,0)
		pinger = (xP,yP,zP)

	Sequence of Events
		1. Initialization of parameters
		2. Synchorization with pinger using a single channel
			A. Chan1 sampled
			B. FFT
			C. Ideal bandpass filter
			D. iFFT
			E. Pulse-Repetitive-Period estimated.
				i. Adjusting trigger delay
			F. Heads centered in frame window.
				i. Adjusting trigger delay
		3. All 4 channels sampled and azimuth to source estimated.
			A. All channels sampled
			B. FFT
			C. Ideal bandpass filter
			D. iFFT
			E. Secondary Pseudo-TOAs
				i. Head triggered left-to-right break threshold
			F. Secondary time delays
				i. Difference in secondary pseudo-TOAs
			G. Primary time delays
				i. Cross-correlation sliver
			H. Comparing primary and secondary time delays
			I. Primary Pseudo-TOA
			J. Sphere radii
			K. Pseudo pinger location
			L. Azimuths
				i.   Checking for complex solutions
				ii.  Saving last 10 azimuths
				iii. Taking median of last 10 azimuths
*/

#include <supportFunctions.h>
#include <stdio.h>
#include <math.h>

// Mathematical constants
static const double PI = 3.141592653589793;

// Pinger parameters
static const double fPinger = 10.0E+3;
static const double tOn = 1.3E-3;
static const double PRT = 2.0;

// Hydrophone parameters

// ADC parameters
static const double fADC = 44.1E+3;
static const double tADC = 0.00002267573696;
static const int N0 = 1024;

// FFT parameters


// Ideal digital bandpass filter parameters
static const double fCenter = 30.0E+3;

int main (void){
	double y[N0];

	// Initializing signal array
	for (int i=0; i<N0; i++) { y[i] = 0; }

	// Filling signal array
	for (int i=0; i<N0; i++) { y[i] = sin(2*PI*fPinger*i*tADC); }

	// Padding tail of signal array with zeros for pulsing
	//for (int i=0; i*tADC > PRT-tOn; i++) { y[i] = 0;}

	// Making multiple pulses
	// 

	// Checking results
	for (int i=0; i<N0; i++) { printf("%i \t %f \n", i, y[i]); }

	return 0;
}


