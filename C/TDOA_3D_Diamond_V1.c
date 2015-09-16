/*
*	Author:  Joshua Simmons
*	Started: September, 2015
*	Status:  INCOMPLETE
*
*	PINGER ASSUMED TO BE INTERMITTENT AT FIXED INTERVALS!!!
*
*	Description: Uses Time-Difference-of-Arrival (TDOA) to determine the azimuth to
*				 pulsing signal source underwater.
*
*	Source Characteristics (Teledyne Benthos)
*		1. Waveform: sinusoidal
*		2. Frequency: 25-40 [kHz] in 0.5 [kHz] increments
*		3. Pulse Length: 4.0 [ms]
*		4. Pulse Repetition: 0.5 OR 1 OR 2 [s]
*       5. Acoustic Output Power: 0.125-5.000 [W]
*       6. Depth Rating: 750 [m]
*
*	3D Cartesian co-ordinate system with the origin centered in the middle of the sensor
*	array. Sensor geometry is square shaped residing all in the XY-plane.
*
*	Sensor layout
*
*		  (Top View)
*	---------------------
*	|                   |
*	|         1         |
*	|                   |
*	|     4       2     |
*	|                   |
*	|         3         |
*	|                   |
*	---------------------
*
*	Coordinates
*		chan1  = ( 0, d,0)		D = sqrt(2)*d
*		chan2  = ( d, 0,0)
*		chan3  = ( 0,-d,0)
*		chan4  = (-d, 0,0)
*		pinger = (xP,yP,zP)
*
*	Sequence of Events
*		1. Initialization of parameters
*		2. Synchorization with pinger using a single channel
*			A. Chan1 sampled
*			B. FFT
*			C. Power
*			    i. Adjust PGA gain if necessary
*			D. Ideal bandpass filter
*			E. iFFT
*			F. Pulse-Repetitive-Period estimated.
*				i. Adjust trigger delay
*			G. Heads centered in frame window.
*				i. Adjust trigger delay
*		3. All 4 channels sampled and azimuth to source estimated
*			A. All channels sampled
*			B. FFT
*			C. Power
*			    i. Adjust PGA gain if necessary
*			D. Ideal bandpass filter
*			E. iFFT
*			F. Secondary Pseudo-TOAs
*				i. Head triggered left-to-right break threshold
*			G. Secondary time delays
*				i. Difference in secondary pseudo-TOAs
*			H. Primary time delays
*				i. Cross-correlation sliver
*			I. Comparing primary and secondary time delays
*			J. Primary Pseudo-TOA
*			K. Sphere radii
*			L. Pseudo pinger location
*			M. Azimuths
*				i.   Checking for complex solutions
*				ii.  Saving last 10 azimuths
*				iii. Taking median of last 10 azimuths
*/

#include "supportFunctions.h"
#include <stdio.h>
#include <math.h>

// Mathematical constants
static const double PI = 3.141592653589793;

// Pinger parameters
static const double fPinger = 10.0E+3;      // Pinger Frequency [Hz]
static const double pulseLength = 1.3E-3;   // Pulse Length [s]
static const double PRT = 2.0;              // Pulse-Repetitive-Period [s]
static const int azimuthArraySize = 10;     // For Taking Medians of Azimuths

// Hydrophone parameters
static const double vP = 1482.0;            // Velocity of Propagation [m/s]
static const double lambda = vP/fPinger;    // Wavelength [m]
static const double D = lambda;             // Sensor spacing [m]
static const double d = D/sqrt(2.0);        // For System of Coordinates [m]

// ADC parameters
static const double fADC = 1800.0E+3;   // ADC Sampling Frequency [Hz]
static const double tADC = 1.0/fADC;    // ADC Samplgin Period [s]
static const int N0 = 1024;             // Frame Size [samples]

// FFT parameters
static const double T0 = N0*tADC;   // Truncation Time Interval [s]
static const double f0 = 1/T0;      // Frequency Resolution [Hz]

// Ideal digital bandpass filter parameters
static const double fCenter = 30.0E+3;  // Center frequency [Hz]
static const double halfChan = 5.0E+3;  // Channel half width [Hz]

// MAIN Variable Declaration
boolean pingerSynced = false;

double chan1t[N0];  // Chan1 Time Data
double chan2t[N0];  // Chan2 Time Data
double chan3t[N0];  // Chan3 Time Data
double chan4t[N0];  // Chan4 Time Data

double chan1f[N0];  // Chan1 Frequency Data
double chan2f[N0];  // Chan2 Frequency Data
double chan3f[N0];  // Chan3 Frequency Data
double chan4f[N0];  // Chan4 Frequency Data

double P1 = 0.0;    // Chan1 Power
double P2 = 0.0;    // Chan2 Power
double P3 = 0.0;    // Chan3 Power
double P4 = 0.0;    // Chan4 Power

double primaryTOAs[16]; // Primary Pseudo TOAs ***THIS IS A FUNCTION***
double secondaryTOAs[4] = {0.0; 0.0; 0.0; 0.0}; // Secondary Pseudo TOAs
    
double primaryTDs[4] = {0.0; 0.0; 0.0; 0.0};    // Primary Time Delays
double secondaryTDs[16];                        // Secondary Time Delays ***THIS IS A FUNCTION***

double sphereRadii[4] = {0.0; 0.0; 0.0; 0.0};
double pingerCoordinates[3] = {0.0; 0.0; 0.0};
double azimuthH  = 0.0;
double azimuthV1 = 0.0;
double azituthV2 = 0.0;
    
double azimuthHArray[10];   // For taking medians
double azimuthV1Array[10];
double azimuthV2Array[10];

int main (void){

    /////////////////////////////////
    // INITIALIZING ARRAYS TO ZERO //
    /////////////////////////////////
    for ( int i = 0; i < N0; i++) {
        chan1t[i] = 0.0;
        chan2t[i] = 0.0;
        chan3t[i] = 0.0;
        chan4t[i] = 0.0;
        chan1f[i] = 0.0;
        chan2f[i] = 0.0;
        chan3f[i] = 0.0;
        chan4f[i] = 0.0;
    }

    for ( int i = 0; i < azimuthArraySize; i++ ) {
        azimuthHArray[i] = 0.0;
        azimuthV1Array[i] = 0.0;
        azimuthV2Array[i] = 0.0;
    }

    ////////////////////////
    // START MAIN ROUTINE //
    ////////////////////////
    while (1) {
        if ( pingerSynced == false )
            //Sample chan1
            FFT();
            Power();
            AdjustPGA();
            BPF();
            iFFT();
            syncPinger();
        else
            //Sample all channels
            FFT();
            Power();
            AdjustPGA();
            BPF();
            iFFT();
            SecondaryTOAs();
            SecondaryTimeDelays();
            PrimaryTimeDelays();
            Compare();
            PrimaryTOAs();
            SphereRadii();
            PingerLocation();
            PingerAzimuth();
    }

	return 0;
}


