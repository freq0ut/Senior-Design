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
static const double T0 = N0/fADC;   // Truncation Time Period [s]
static const double f0 = 1.0/T0;    // Frequency Resolution [Hz]

// Pinger
static const double fPinger = 37.0E+3;  // Pinger Frequency [Hz] THIS MAY BECOME A VARIABLE

// Hydrophones
static const double vP = 1482.0;                    // Velocity of Propagation [m/s]
static const double lambda = vP/fPinger;            // Wavelength [m]
static const double D = lambda;                     // Hydrophone Spacing [m]
static const double d = lambda/1.414213562373095;   // For System of Coordinates [m]

// Processor
static const int medianSize = 10;   // For Taking Medians of Azimuths
static const int N0 = 1024;         // Frame Size [samples]

// XCorr
static const int lagBounds = (int) (D*fADC/vP+1);   // XCorr boundary limits
static const int pkCounterMax = (int) (D/lambda+1); // Max number of peaks for Max(XCorr)

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// VARIABLE DECLARATION AND INITIALIZATION ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// FFT
double powerx_f;                        // Power [W] (for adjusting PGA)
double f[1+N0] = {N0};                  // Double-Sided Frequency Vector
double complex chanx_f[1+N0] = {N0};    // Chanx in Frequency Domain
double complex H[1+N0] = {N0};          // Ideal Bandpass Filter

// Pinger
int pingerSynced = FALSE;               // Processor Initially Unsynchonized
double PRT = 0.5;                       // Pinger Pulse-Repetitive-Period [s]
double complex pingerLoc[1+3] = {3};    // Pinger Location Co-ordinates

// TDOA
double azimuthH;    // Single horizontal azimuth estimate [deg]
double azimuthV1;   // Single vertical azimuth 1 estimate [deg]
double azimuthV2;   // Single vertical azimuth 2 estimate [deg]
double tD2;         // Channel 2 Time Delay [s]
double tD3;         // Channel 3 Time Delay [s]
double tD4;         // Channel 4 Time Delay [s]
double TOA;         // Time-of-Arrival [s]
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

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// START MAIN ROUTINE ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int main (void){

    // Constructing Double-Sided Frequency Vector
    for (int i = 1; i <= N0; i++) {
        f[i] = -N0/2 + (i-1);
        f[i] = f0 * f[i];
    }

    // Constructing Ideal Digital Bandpass Filter
    for (int i = 1; i <= N0; i++) {
        if ( abs(f[i]) >= fCent-halfChan && abs(f[i]) <= fCent+halfChan ) {
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
/*
    while (1) {

        if ( pingerSynced == 0 ) {
            SampleADC();// Output goes to chan1_t, chan2_t, chan3_t, chan4_t
            TDOA_FFT(t,chan1_t,chanx_f);
            powerx_f = Power_f(f,chanx_f);

            if      (powerx_f >= powerMax) { decreasePGA(); }
            else if (powerx_f <= powerMin) { increasePGA(); }
            else;

            // Bandpass Filtering
            for (int i=1; i <= H[0]; i++) { chanx_f[i] = chanx_f[i] * H[i];}

            TDOA_iFFT(f,chanx_f,chan1_t);
            syncPinger();
            Delay(PRT);
        }

        else if ( pingerSynced == 1 ) {
            SampleADC();// Output goes to chan1_t, chan2_t, chan3_t, chan4_t
            TDOA_FFT(t,chan1_t,chanx_f);
            powerx_f = Power_f(f,chanx_f);

            if      (powerx_f >= powerMax) { decreasePGA(); }
            else if (powerx_f <= powerMin) { increasePGA(); }
            else;

            for (int i=1; i <= H[0]; i++) { chanx_f[i] = chanx_f[i] * H[i]; }

            TDOA_iFFT(f,chanx_f,chan1_t);
            breakwall_TOAs[1] = Breakwall(chan1_t,THD,0) * tADC;
            
            // Resynchronize with the pinger
            if ( breakwall_TOAs[1] == -1 ) { pingerSynced == 0; }
            else;
            
            if ( pingerSynced == 1 ) {
                TDOA_FFT(t,chan2_t,chanx_f);
                for (int i=1; i <= H[0]; i++) { chanx_f[i] = chanx_f[i] * H[i]; }
                TDOA_iFFT(f,chanx_f,chan2_t);
                breakwall_TOAs[2] = Breakwall(chan2_t,THD,0) * tADC;
                breakwall_tDs[2] = breakwall_TOAs[1] - breakwall_TOAs[2];
                XC_Bounded(chan1_t,chan2_t,iBound);// Where is the return argument going?
                localMaxima(XC1x);// Where is the return argument going?
                td2 = Compare(breakwall_tDs[2], XC_tDs);

                TDOA_FFT(t,chan3_t,chanx_f);
                for (int i=1; i <= H[0]; i++) { chanx_f[i] = chanx_f[i] * H[i]; }
                TDOA_iFFT(f,chanx_f,chan3_t);
                breakwall_TOAs[3] = Breakwall(chan3_t,THD,0) * tADC;
                breakwall_tDs[3] = breakwall_TOAs[1] - breakwall_TOAs[3];
                XC_Bounded(chan1_t,chan3_t,iBound);// Where is the return argument going?
                localMaxima(XC1x);// Where is the return argument going?
                td3 = Compare(breakwall_tDs[3], XC_tDs);

                TDOA_FFT(t,chan4_t,chanx_f);
                for (int i=1; i <= H[0]; i++) { chanx_f[i] = chanx_f[i] * H[i]; }
                TDOA_iFFT(f,chanx_f,chan4_t);
                breakwall_TOAs[4] = Breakwall(chan4_t,THD,0) * tADC;
                breakwall_tDs[4] = breakwall_TOAs[1] - breakwall_TOAs[4];
                XC_Bounded(chan1_t,chan4_t,iBound);// Where is the return argument going?
                localMaxima(XC1x);// Where is the return argument going?
                td4 = Compare(breakwall_tDs[4], XC_tDs);

                TOA = CalcTOA();
                CalcSphereRadii();
                CalcPingerLoc();
                CalcPingerAzimuths();
                Delay(PRT);
            }
        }

        else;
    }
*/
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// END OF PROGRAM ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

