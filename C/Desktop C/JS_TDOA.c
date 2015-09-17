/***************************************************************************************************
Author:      Joshua Simmons
Start Date:  September 11, 2015
Description: This C program generates an array of a pulsed SINE wave for audio speaker transmission.
***************************************************************************************************/

#include <stdio.h>
#include <math.h>

static const double PI = 3.141592653589793;

static const double fPinger = 10E+3;
static const double tOn = 1.3E-3;
static const double PRT = 2;

static const double tADC = 0.00002267573696; // 44.1 kHz
static const int N0 = 1024;

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


