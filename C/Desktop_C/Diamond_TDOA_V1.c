/*
*	Author:  Joshua Simmons, Zack Goyetche, Michael Daub, and Shane Stroh
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
*		5. Acoustic Output Power: 0.125-5.000 [W]
*		6. Depth Rating: 750 [m]
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
*		1. Constant and variable declarations and initializations
*		2. Synchorization with pinger using a single channel
*			A. Chan1 sampled
*			B. FFT
*			C. Power
*			    i. Adjust PGA gain if necessary
*			D. Ideal bandpass filter
*			E. iFFT
*			F. Pinger synchronized
*				i.  Pulse-Repetitive-Period estimated
*                       a. PRT adjusted
*				ii. Heads centered in frame window
*						a. PRT adjusted
*		3. All 4 channels sampled and azimuth to source estimated
*			A. All channels sampled
*			B. FFT
*			C. Power
*			    i. Adjust PGA gain if necessary
*			D. Ideal bandpass filter
*			E. iFFT
*			F. Breakwall TOAs
*				i. Head triggered left-to-right break threshold detection
*			G. Break time delays
*				i. Difference in breakwall TOAs
*			H. XC time delays
*				i. Bounded cross-correlation 
*			I. Comparing breakwall and XC time delays to find the superior TOA
*			K. Sphere radii
*			L. Pseudo pinger location (could be complex)
*			M. Azimuths
*				i.   Checking for complex solutions (arg(complex number) == 0)
*				ii.  Saving last 10 azimuths if REAL
*				iii. Taking median of last 10 azimuths
*
*   Programming Conventions:
*       1. All arrays have their size stored in their zeroth element.
*           Example: int myArray[1+10] = {10,1,2,3,4,5,6,7,8,9,10}; 
*/

#include "js_tdoa.h"

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// PIN DECLARATIONS /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// NONE!

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// CONSTANT DECLARATION AND INITIALIZATION ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// ADC
static const double fADC = 1800.0E+3; // ADC Sampling Frequency [Hz]
static const double powerMin =  1.0;  // Minimum Signal Power [W]
static const double powerMax = 10.0;  // Maximum Signal Power [W]

// Azimuths
static const int medianSize = 10; // Azimuths only.

// Bandpass Filter
static const double fCenter = 30.0E+3; // Center Frequency [Hz]
static const double halfChan = 5.0E+3; // Channel Half-Width [Hz]

// FFT
static const int N0 = 1024; // Frame Size [samples]
static const double f0 = fADC/N0; // Frequency Resolution [Hz]

// Pinger
static const double fPinger = 37.0E+3;  // Pinger Frequency [Hz]
//  fPinger may become a variable that can be determined by the FFT.
//  The constant declared here my become fPingerFloor.
// fPingerMin
// PRT_Min
// PRT_Max

// Hydrophones
static const double vP = 1482.0; // Velocity of Propagation [m/s]
static const double lambda = vP/fPinger; // Wavelength [m]
static const double D = lambda; // Hydrophone Spacing [m]
static const double d = D/1.414213562373095; // For System of Coordinates [m]

// Time Delays
static const double threshold = 0.5; // CalcTimeDelay() and CenterWindow()
static const int lagBounds = (int) (D*fADC/vP+1);   // XCorr boundary limits
static const int pkCounterMax = (int) (D/lambda+1); // Max number of peaks for Max(XCorr)

int main (void) {
/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// VARIABLE DECLARATION AND INITIALIZATION ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

    // FFT
    double f[1+N0] = {N0}; // Double-Sided Frequency Vector
    double _Complex chanx_f[1+N0] = {N0}; // Chanx in Frequency Domain
    double _Complex H[1+N0] = {N0};       // Ideal Bandpass Filter

    // Pinger
    int pingerSynced = FALSE;
    double pingerLocs[1+3] = {3,0,0,0};
    double PRT = 0.5;   // Pinger Pulse-Repetitive-Period [s]
                        // This initialized value is a fraction of our best guess of the PRT.
                        // In the algorithm that follows, this value will change for synchronization.
                        // Never make this initialized value at or above the actual PRT or else
                        // the algorithm might delay twice as long as it needs to. 

    // TDOA
    int medianCounter = 1;  // For logging single azimuth estimates into their respective arrays.
    double azimuthH  = 0.0; // Single horizontal azimuth estimate [deg]
    double azimuthV1 = 0.0; // Single vertical azimuth 1 estimate [deg]
    double azimuthV2 = 0.0; // Single vertical azimuth 2 estimate [deg]
    double tD2 = 0.0;  // Channel 2 Time Delay [s]
    double tD3 = 0.0;  // Channel 3 Time Delay [s]
    double tD4 = 0.0;  // Channel 4 Time Delay [s]
    double TOA1 = 0.0; // Time-of-Arrival for Channel 1 [s]
    double azimuthHArray [1+medianSize] = {medianSize}; // Array of previous horizontal azimuth estimates
    double azimuthV1Array[1+medianSize] = {medianSize}; // Array of previous vertical azimuth 1 estimates
    double azimuthV2Array[1+medianSize] = {medianSize}; // Array of previous vertical azimuth 2 estimates

    // Time Sampled Data
    double t[1+N0] = {N0}; // Time Vector
    double chan1_t[1+N0] = {N0}; // Channel 1 Time Sampled Data
    double chan2_t[1+N0] = {N0}; // Channel 2 Time Sampled Data
    double chan3_t[1+N0] = {N0}; // Channel 3 Time Sampled Data
    double chan4_t[1+N0] = {N0}; // Channel 4 Time Sampled Data

    // Initialization and Creation of Arrays of Size N0
    for (int i = 1; i <= N0; i++) {

        // Initializing
        f[i] = 0.0;
        chanx_f[i] = 0.0;
        H[i] = 0.0;

        t[i] = 0.0;
        chan1_t[i] = 0.0;
        chan2_t[i] = 0.0;
        chan3_t[i] = 0.0;
        chan4_t[i] = 0.0;

        // Constructing Time Vector
        t[i] = (i-1)/fADC;

        // Constructing Double-Sided Frequency Vector
        f[i] = -N0/2 + (i-1);
        f[i] = f0 * f[i];

        // Constructing Ideal Digital Bandpass Filter
        if ( abs(f[i]) >= fCenter-halfChan && abs(f[i]) <= fCenter+halfChan ) {
            H[i] = 1.0;
        }
        else;
    }

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// START MAIN ROUTINE ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

    while (TRUE) {

        if ( pingerSynced == FALSE ) {
            /////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////// SYNCHRONIZING WITH PINGER ///////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////

            SampleAllChans(fADC, chan1_t, chan2_t, chan3_t, chan4_t);
            FFT(chan1_t,chanx_f);
            AdjustPGA(f,chanx_f, powerMin, powerMax);

            // Bandpass Filtering Channel 1
            for (int i=1; i <= N0; i++) {
                chanx_f[i] = chanx_f[i] * H[i];
            }

            iFFT(chanx_f,chan1_t);
            SyncPinger(chan1_t, pingerSynced, PRT);
        }

        else if ( pingerSynced == TRUE ) {
            /////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////// CHECK FOR LOSS OF SYNCHRONIZATION ////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////

            SampleAllChans(fADC, chan1_t, chan2_t, chan3_t, chan4_t);
            FFT(chan1_t,chanx_f);
            AdjustPGA(f,chanx_f, powerMin, powerMax);

            // Bandpass Filtering Channel 1
            for (int i=1; i <= N0; i++) {
                chanx_f[i] = chanx_f[i] * H[i];
            }

            iFFT(chanx_f,chan1_t);
            CenterWindow(chan1_t, N0, fADC, threshold, pingerSynced, PRT, TOA1);   // Fine adjusting PRT
            // Resynchronize with the pinger
            if ( TOA1 < 0 ) { pingerSynced == FALSE; }
            else;
            
            if ( pingerSynced == TRUE ) {

                /////////////////////////////////////////////////////////////////////////////////////////
                ///////////////////////////////// CHANNEL 2 TIME DELAY //////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////

                FFT(chan2_t,chanx_f);

                // Bandpass Filtering Channel 2
                for (int i=1; i <= H[0]; i++) {
                    chanx_f[i] = chanx_f[i] * H[i];
                }

                iFFT(chanx_f,chan2_t);

                tD2 = CalcTimeDelay(chan1_t, chan2_t, THD, TOA1, lagBounds, pkCounterMax);

                /////////////////////////////////////////////////////////////////////////////////////////
                ///////////////////////////////// CHANNEL 3 TIME DELAY //////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////

                FFT(chan3_t,chanx_f);

                // Bandpass Filtering Channel 3
                for (int i=1; i <= H[0]; i++) {
                    chanx_f[i] = chanx_f[i] * H[i];
                }

                iFFT(chanx_f,chan3_t);

                tD3 = CalcTimeDelay(chan1_t, chan3_t, THD, lagBounds, pkCounterMax);

                /////////////////////////////////////////////////////////////////////////////////////////
                ///////////////////////////////// CHANNEL 4 TIME DELAY //////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////

                FFT(chan4_t,chanx_f);

                // Bandpass Filtering Channel 4
                for (int i=1; i <= H[0]; i++) {
                    chanx_f[i] = chanx_f[i] * H[i];
                }

                iFFT(chanx_f,chan4_t);

                tD4 = CalcTimeDelay(chan1_t, chan4_t, THD, lagBounds, pkCounterMax);

                /////////////////////////////////////////////////////////////////////////////////////////
                ///////////////////////////////// CALCULATING AZIMUTHS //////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////

                CalcDiamondPingerLocation(d, td2, td3, td4, TOA, vP, pingerLocs);

                // Checking for complex pinger coordinates
                if ( cimag(pingerLoc[1]) == 0 && cimag(pingerLoc[2]) && cimag(pingerLoc[3]) == 0 ) {

                    CalcPingerAzimuths(pingerLocs, azimuthH, azimuthV1, azimuthV2);

                    // Indices for running median
                    if ( medianCounter % (medianSize+1) == 0 ) {
                        medianCounter = 1;  // counter reset
                    }
                    else;

                    azimuthHArray [medianCounter] = azimimuthH;
                    azimuthV1Array[medianCounter] = azimimuthV1;
                    azimuthV2Array[medianCounter] = azimimuthV2;

                    azimuthH  = Median(azimuthHArray);  // fix bug Median() overwrites
                    azimuthV1 = Median(azimuthV1Array);
                    azimuthV2 = Median(azimuthV2Array);

                    printf("\nAzimuthA = %f\tAzimuthV1 = %f\tAzimuthV2 = %f", azimuthH, azimuthV1, azimuthV2);

                    medianCounter++;
                }
                else;
            }
            else;
        }
        else;

        DelaySampleTrigger(PRT);
    }

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// END OF PROGRAM ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
