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
*			F. Break TOAs
*				i. Head triggered left-to-right break threshold
*			G. Break time delays
*				i. Difference in secondary pseudo-TOAs
*			H. XC time delays
*				i. Cross-correlation sliver
*			I. Comparing break and XC time delays
*			J. Primary TOA
*			K. Sphere radii
*			L. Pseudo pinger location
*			M. Azimuths
*				i.   Checking for complex solutions
*				ii.  Saving last 10 azimuths
*				iii. Taking median of last 10 azimuths
*/

//#include "JS_TDOA.h"
//#include "JS_cmath.h"

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// CONSTANT DECLARATIONS //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// Mathematical Constants
static const double PI = acos(-1.0);

// Hydrophone Parameters
static const double vP = 1482.0;        // Velocity of Propagation [m/s]
static const double lambda = vP/fPinger;// Wavelength [m]
static const double D = lambda;         // Sensor spacing [m]
static const double d = D/sqrt(2.0);    // For System of Coordinates [m]

// ADC Parameters
static const double fADC = 1800.0E+3;// ADC Sampling Frequency [Hz]
static const double tADC = 1.0/fADC; // ADC Sampling Period [s]
static const int N0 = 1024;          // Frame Size [samples]

// FFT Parameters
static const double T0 = N0*tADC;// Truncation Time Interval [s]
static const double f0 = 1.0/T0; // Frequency Resolution [Hz]

// Ideal Digital Bandpass Filter Parameters
static const double fCenter = 30.0E+3;// Center frequency [Hz]
static const double halfChan = 5.0E+3;// Channel half width [Hz]

// Processor Parameters
static const int maxPeaks = (int) ceiling(D/lambda) + 1;// Maximum Number of Peaks for XC
static const int medianSize = 10;                       // For Taking Medians of Azimuths

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// VARIABLE DECLARATIONS //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// Estimated Pinger Parameters
double fPinger = 37.0E+3;   // Pinger Frequency [Hz]
double pulseLength = 4.0E-3;// Pulse Length [s]
double PRT = 2.0;           // Pulse-Repetitive-Period [s]

// Are we synchronized with the pinger?
int pingerSynced = 0;

// Time Data
double chan1t[N0];
double chan2t[N0];
double chan3t[N0];
double chan4t[N0];

// Frequency Data
double complex chan1f[N0];
double complex chan2f[N0];
double complex chan3f[N0];
double complex chan4f[N0];

// Power
double power1, power2, power3, power4;

// Time-of-Arrival and Time Delays
double primaryTOAs[maxPeaks];
double secondaryTOAs[4];
double primaryTDs[4];
double secondaryTDs[maxPeaks];

// Pinger Locating
double sphereRadii[4];
double pingerCoordinates[3];
double azimuthH, azimithV1, azimuthV2;
double azimuthHArray[medianSize];
double azimuthV1Array[medianSize];
double azimuthV2Array[medianSize];

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// START MAIN ROUTINE ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int main (void){
/*
    while (1) {

        if ( pingerSynced == 0 ) {
            //Sample chan1
            FFT();
            Power();
            AdjustPGA();
            BPF();
            iFFT();
            syncPinger();
            Delay();
        }

        else if ( pingerSynced == 1 ) {
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
            Delay();
        }

        else;
    }
*/
	return 0;
}

