#include "js_tdoa.h"

void AdjustPGA (double* f, double _Complex* chanx_f, double powerMin, double powerMax) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Calculates the signal power and adjusts the gain of the ADC PGA accordingly.
    *
    *   Status: incomplete
    *
    *   Notes:
    *       1. Will have to wait until familarization with ADC to be implemented.
    */

    double power = 0.0;

    power = SignalPower(f,chanx_f);

    if (power >= powerMax) {
        printf("\nDecreased PGA Gain.");
    }
    else if (power <= powerMin) {
        printf("\nIncreased PGA Gain.");
    }
    else;

    return;
}

void CalcPingerAzimuths(double* pingerLocs, double* azimuthH, double* azimuthV1, double* azimuthV2) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Calculates the pinger azimuths from its Cartesian coordinates.
    *
    *   Status: untested
    *
    *   Notes: NONE! 
    */

    // Determining pinger azimuths in radians
    *azimuthH  = atan2(   pingerLocs[2],pingerLocs[1]);
    *azimuthV1 = atan2(   pingerLocs[3],pingerLocs[1]);
    *azimuthV2 = atan2(-1*pingerLocs[3],pingerLocs[1]);

    // Wrapping angle [0,2pi]
    *azimuthH  = WrapTo2Pi(*azimuthH);
    *azimuthV1 = WrapTo2Pi(*azimuthV1);
    *azimuthV2 = WrapTo2Pi(*azimuthV2);
    
    // Converting to degrees
    *azimuthH  = (*azimuthH)  * (180.0/pi);
    *azimuthV1 = (*azimuthV1) * (180.0/pi);
    *azimuthV2 = (*azimuthV2) * (180.0/pi);

    return;
}

void CalcDiamondPingerLocation(double d, double td2, double td3, double td4, double TOA1, double vP, double* pingerLocs){
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Calculates the pinger coordinates for a diamond sensor configuration using TDOA.
    *
    *   Status: untested
    *
    *   Notes:
    *       1. As long as TOA1 is large then we do not have to worry about complex solutions. 
    */

    double R1 = 0.0;
    double R2 = 0.0;
    double R3 = 0.0;
    double R4 = 0.0;

    R1 = sqrt(vp*(TOA1));
    R1 = sqrt(vp*(TOA1+td2));
    R1 = sqrt(vp*(TOA1+td3));
    R1 = sqrt(vp*(TOA1+td4));

    *pingerLocs[1] = (R4*R4 - R2*R2) / (4.0*d);
    *pingerLocs[2] = (R3*R3 - R1*R1) / (4.0*d);
    *pingerLocs[3] = sqrt( R1*R1 - pingerLocs[1]*pingerLocs[1] - (pingerLocs[2]-d)*(pingerLocs[2]-d) );

    return;
}

double CalcTimeDelay (double* chan1_t, double* chanx_t, double threshold, double TOA1, int lagBounds, int pkCounterMax) {
     /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Computes the time delay for TDOA directional finding via cross-correlation. 
    *
    *   Status: untested
    *
    *   Notes: NONE! 
    */

    double tD = 0.0;
    double tD_BW = 0.0;
    double TOA234 = 0.0;

    double tD_XCs[1+pkCounterMax];
    double XC[2*(1+lagBounds)];
    double XC_Lags[2*(1+lagBounds)];

    tD_XCs[0] = pkCounterMax;
    XC[0] = 2*lagBounds+1;
    XC_Lags[0] = 2*lagBounds+1;

    char* dir = "LR";   // Head triggered. Switch to "RL" for tail triggered. 
    TOA234 = BreakWall(dir,chanx_t,threshold) / fADC;

    tD_BW = TOA1 - TOA234; // Head triggered. Switch to TOA234 - TOA1 for tail triggered. 

    XCorr(chan1_t,chanx_t,lagBounds,XC_Lags,XC);

    LocalMaxima(XC,2,threshold,pkLocs);  // Adjust the second value (iMPD) when fPinger is a variable

    for ( int i=1; i <= pkCounterMax; i++) {
        tD_XCs[i] = XC_Lags[pkLocs[i]] / fADC;
    }

    tD = Compare(tD_BW, tD_XCs);

    return tD;
}

void CenterWindow (double* chan1_t, int N0, double fADC, double threshold, int* pingerSynced, double* PRT, double *TOA1) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Fine adjusts the PRT to maintain pinger synchronization. 
    *
    *   Status: untested
    *
    *   Notes: NONE! 
    */

    double tError = 0.0;
    double tCenter = (N0/2+1)/fADC;

    char* dir = "LR";
    *TOA1 = BreakWall(dir,chan1_t,threshold) / fADC;

    if ( TOA1 > 0 ) {
        tError = tCenter - TOA1;

        if ( TOA1 > tCenter ) {
            *PRT -= tError;
        }
        else {
            *PRT += tError;
        }
    }
    else {
        *pingerSynced = FALSE;
        printf("\nWarning in CenterWindow(). Lost sync with pinger.");
    }

    return;
}

void DelaySampleTrigger(double PRT) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description:
    *
    *   Status: incomplete
    *
    *   Notes: NONE! 
    */

    printf("\nDelayed by %f[ms]", PRT*1E+3);

    return;
}

void SampleAllChans (double fADC, double* chan1_t, double* chan2_t, double* chan3_t, double* chan4_t) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Communicates with the ADC to sample all channels asynchrously. 
    *
    *   Status: incomplete
    *
    *   Notes: NONE! 
    */

    // Simulating ADC sampling
    ReadCSV(fileName,chan1_t);
    ReadCSV(fileName,chan2_t);
    ReadCSV(fileName,chan3_t);
    ReadCSV(fileName,chan4_t);

    printf("\nSampled All Channels.");

    return;
}

void SyncPinger (double* chan1_t, int* pingerSynced, double* PRT) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Synchronizes processor with pinger by changing PRT. 
    *
    *   Status: incomplete
    *
    *   Notes: NONE! 
    */

    if ( pingerSynced == FALSE ) {
        printf("\nAttemping to synchronize with pinger...");
    }
    else {
        printf("\nSuccessfully synchronized with the pinger.");
    }

    return;
}
