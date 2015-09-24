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

#include "js_cmath.h"
#include "js_tdoa.h"

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// PIN DECLARATIONS /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// NONE!

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// CONSTANT DECLARATION AND INITIALIZATION ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// ADC
static const double fADC = 1800.0E+3;   // ADC Sampling Frequency [Hz]
static const double powerMax = 10.0;    // Maximum Signal Power [W]
static const double powerMin =  1.0;    // Minimum Signal Power [W]

// Bandpass Filter
static const double fCenter = 30.0E+3;  // Center Frequency [Hz]
static const double halfChan = 5.0E+3;  // Channel Half-Width [Hz]

// FFT
static const int N0 = 1024;         // Frame Size [samples]
static const double f0 = fADC/N0;    // Frequency Resolution [Hz]

// Pinger
static const double fPinger = 37.0E+3;  // Pinger Frequency [Hz] THIS MAY BECOME A VARIABLE

// Hydrophones
static const double vP = 1482.0;                    // Velocity of Propagation [m/s]
static const double lambda = vP/fPinger;            // Wavelength [m]
static const double D = lambda;                     // Hydrophone Spacing [m]
static const double d = D/1.414213562373095;   // For System of Coordinates [m]

// Processor
static const int medianSize = 10;       // For Taking Medians of Azimuths
static const double threshold = 0.5;    // Y-Bar of the top lobe of a Sine?

// XCorr
static const int lagBounds = (int) (D*fADC/vP+1);   // XCorr boundary limits
static const int pkCounterMax = (int) (D/lambda+1); // Max number of peaks for Max(XCorr)

int main (void) {
/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// VARIABLE DECLARATION AND INITIALIZATION ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

    // FFT
    double power1;                          // Channel 1 Power [W] (for adjusting PGA)
    double f[1+N0] = {N0};                  // Double-Sided Frequency Vector
    double complex chanx_f[1+N0] = {N0};    // Chanx in Frequency Domain
    double complex H[1+N0] = {N0};          // Ideal Bandpass Filter

    // Pinger
    int pingerSynced = FALSE;
    double PRT = 0.5;                       // Pinger Pulse-Repetitive-Period [s]
    double complex pingerLoc[1+3] = {3};    // Pinger Location Co-ordinates

    // TDOA
    int medianCounter = 1;  // For logging single azimuth estimates into their respective arrays.
    double azimuthH;        // Single horizontal azimuth estimate [deg]
    double azimuthV1;       // Single vertical azimuth 1 estimate [deg]
    double azimuthV2;       // Single vertical azimuth 2 estimate [deg]
    double tD2;             // Channel 2 Time Delay [s]
    double tD3;             // Channel 3 Time Delay [s]
    double tD4;             // Channel 4 Time Delay [s]
    double TOA;             // Time-of-Arrival [s]
    double azimuthHArray [1+medianSize] = {medianSize}; // Array of previous horizontal azimuth estimates
    double azimuthV1Array[1+medianSize] = {medianSize}; // Array of previous vertical azimuth 1 estimates
    double azimuthV2Array[1+medianSize] = {medianSize}; // Array of previous vertical azimuth 2 estimates
    double breakwall_tDs [1+4] = {4,0}; // Breakwall Time Delays
    double breakwall_TOAs[1+4] = {4};   // Breakwall Time-of-Arrivals
    double sphereRadii[1+4] = {4};      // Sphere Radii (multilateration)

    // Time Sampled Data
    double t[1+N0] = {N0};          // Time Vector
    double chan1_t[1+N0] = {N0};    // Channel 1 Time Sampled Data
    double chan2_t[1+N0] = {N0};    // Channel 2 Time Sampled Data
    double chan3_t[1+N0] = {N0};    // Channel 3 Time Sampled Data
    double chan4_t[1+N0] = {N0};    // Channel 4 Time Sampled Data

    // XCorr
    int XCorr1x_Lags[1+(2*lagBounds+1)] = {2*lagBounds+1};  // Cross-correlation lags
    double XCorr1x  [1+(2*lagBounds+1)] = {2*lagBounds+1};  // Cross-correlations
    double XCorr_tDs[1+pkCounterMax] = {pkCounterMax};      // Cross-correlation Time Delays

    // Constructing Double-Sided Frequency Vector
    for (int i = 1; i <= N0; i++) {
        f[i] = -N0/2 + (i-1);
        f[i] = f0 * f[i];
    }

    // Constructing Ideal Digital Bandpass Filter
    for (int i = 1; i <= N0; i++) {
        if ( abs(f[i]) >= fCenter-halfChan && abs(f[i]) <= fCenter+halfChan ) {
            H[i] = 1.0;
        }
        else {
            H[i] = 0.0;
        }
    }

    // Constructing Time Vector
    for (int i = 1; i <= N0; i++) {
        t[i] = (i-1)/fADC;
    }

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// START MAIN ROUTINE ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*
    while (TRUE) {

        if ( pingerSynced == FALSE ) {
            SampleChan1(chan1_t);
            FFT(t,chan1_t,chanx_f);
            power1 = SignalPower(f,chanx_f);

            // Adjust PGA
            if      (power1 >= powerMax)    { adjustPGA(0); }
            else if (power1 <= powerMin)    { adjustPGA(1); }
            else;

            // Bandpass Filtering
            for (int i=1; i <= N0; i++) { chanx_f[i] = chanx_f[i] * H[i];}

            iFFT(f,chanx_f,chan1_t);

            printf("\nAttempting to synchonize with the pinger...");

            PRT = syncPinger();
        }

        else if ( pingerSynced == TRUE ) {
            SampleAllChans(chan1_t, chan2_t, chan3_t, chan4_t);
            FFT(t,chan1_t,chanx_f);
            power1 = SignalPower(f,chanx_f);

            // Adjust PGA
            if      (power1 >= powerMax)    { adjustPGA(0); }
            else if (power1 <= powerMin)    { adjustPGA(1); }
            else;

            // Bandpass Filtering
            for (int i=1; i <= N0; i++) { chanx_f[i] = chanx_f[i] * H[i]; }

            iFFT(f,chanx_f,chan1_t);
            breakwall_TOAs[1] = Breakwall(chan1_t,threshold,0) / fADC;
            
            // Resynchronize with the pinger
            if ( breakwall_TOAs[1] == -1 ) { pingerSynced == FALSE; }
            else;
            
            if ( pingerSynced == TRUE ) {
                FFT(t,chan2_t,chanx_f);
                for (int i=1; i <= H[0]; i++) { chanx_f[i] = chanx_f[i] * H[i];}
                iFFT(f,chanx_f,chan2_t);
                breakwall_TOAs[2] = Breakwall(chan2_t,THD,0) / fADC;
                breakwall_tDs[2] = breakwall_TOAs[1] - breakwall_TOAs[2];
                XCorrBounded(chan1_t,chan2_t,lagBounds, XCorr1x_Lags, XCorr1x);
                localMaxima(XC1x);// Where is the return argument going?
                td2 = Compare(breakwall_tDs[2], XC_tDs);

                FFT(t,chan3_t,chanx_f);
                for (int i=1; i <= H[0]; i++) { chanx_f[i] = chanx_f[i] * H[i]; }
                iFFT(f,chanx_f,chan3_t);
                breakwall_TOAs[3] = Breakwall(chan3_t,THD,0) / fADC;
                breakwall_tDs[3] = breakwall_TOAs[1] - breakwall_TOAs[3];
                XCorrBounded(chan1_t,chan3_t,lagBounds, XCorr1x_Lags, XCorr1x);
                localMaxima(XC1x);// Where is the return argument going?
                td3 = Compare(breakwall_tDs[3], XC_tDs);

                FFT(t,chan4_t,chanx_f);
                for (int i=1; i <= H[0]; i++) { chanx_f[i] = chanx_f[i] * H[i]; }
                iFFT(f,chanx_f,chan4_t);
                breakwall_TOAs[4] = Breakwall(chan4_t,threshold,0) / fADC;
                breakwall_tDs[4] = breakwall_TOAs[1] - breakwall_TOAs[4];
                XCorrBounded(chan1_t,chan4_t,lagBounds, XCorr1x_Lags, XCorr1x);
                localMaxima(XC1x, );// Where is the return argument going?
                td4 = Compare(breakwall_tDs[4], XC_tDs);

                CalcPingerLoc(td2, td3, td4, TOA, vP, pingerLoc);

                // Checking for pingerLoc complex solutions
                if ( cimag(pingerLoc[1]) == 0 && cimag(pingerLoc[2]) && cimag(pingerLoc[3]) == 0 ) {

                    CalcPingerAzimuths(pingerLoc, azimuthH, azimuthV1, azimuthV2);

                    if ( medianCounter % (medianSize+1) == 0 ) { medianCounter = 1; }
                    else;

                    azimuthHArray [medianCounter] = azimimuthH;
                    azimuthV1Array[medianCounter] = azimimuthV1;
                    azimuthV2Array[medianCounter] = azimimuthV2;

                    azimuthH  = Median(azimuthHArray);
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

        Delay(PRT);
    }
*/
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// END OF PROGRAM ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
