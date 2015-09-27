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
        // Decrease PGA gain
    }
    else if (power <= powerMin) {
        // Increase PGA gain
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

    double x = pingerLocs[1];
    double y = pingerLocs[2];
    double z = pingerLocs[3];

    // Determining pinger azimuths in radians
    *azimuthH  = atan2(   y,x);
    *azimuthV1 = atan2(   z,x);
    *azimuthV2 = atan2(-1*z,x);

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

void CalcDiamondPingerLocation(double d, double td2, double td3, double td4, double TOA, double vP, double* pingerLocs){
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
    *       1. As long as TOA is large then we do not have to worry about complex solutions. 
    */

    double x = pingerLocs[1];
    double y = pingerLocs[2];

    double R1 = 0.0;
    double R2 = 0.0;
    double R3 = 0.0;
    double R4 = 0.0;

    R1 = sqrt(vp*(TOA));
    R1 = sqrt(vp*(TOA+td2));
    R1 = sqrt(vp*(TOA+td3));
    R1 = sqrt(vp*(TOA+td4));

    *pingerLocs[1] = (R4*R4 - R2*R2) / (4.0*d);
    *pingerLocs[2] = (R3*R3 - R1*R1) / (4.0*d);
    *pingerLocs[3] = sqrt( R1*R1 - x*x* - (y-d)*(y-d) );

    return;
}

double CalcTimeDelay (double* chanx_t, double threshold, double TOA1, int lagBounds, int pkCounterMax) {
     /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Computes the time delay for TDOA directional finding via cross-correlation. 
    *
    *   Status: INCOMPLETE
    *
    *   Notes: NONE! 
    */

    double tD = 0.0;
    double tD_BW = 0.0;
    double TOA234 = 0.0;

    double tD_XCs[1+pkCounterMax];
    tD_XCs[0] = pkCounterMax;

    char* dir = "LR";
    TOA234 = BreakWall(dir,chanx_t,threshold) / fADC;

    tD_BW = TOA1 - TOA234;

    /*
    breakwall_TOAs[2] = Breakwall(chan2_t,THD,0) / fADC;
    breakwall_tDs[2] = breakwall_TOAs[1] - breakwall_TOAs[2];
    XCorrBounded(chan1_t,chan2_t,lagBounds, XCorr1x_Lags, XCorr1x);
    localMaxima(XC1x);// Where is the return argument going?
    td2 = Compare(breakwall_tDs[2], XC_tDs);
    int XCorr1x_Lags[1+(2*lagBounds+1)] = {2*lagBounds+1}; // Cross-correlation lags
    double XCorr1x  [1+(2*lagBounds+1)] = {2*lagBounds+1}; // Cross-correlations
    */

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

    double tCenter = (N0/2+1) / fADC;

    char* dir = "LR";
    *TOA1 = BreakWall(dir,chan1_t,threshold) / fADC;

    if ( *TOA1 > tCenter ) {
        *PRT -= tError;
    }
    else {
        *PRT += tError;
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
    *   Status: INCOMPLETE
    *
    *   Notes: NONE! 
    */

    return;
}

void SampleAllChans (double fADC, double* chan1_t, double* chan2_t, double* chan3_t, double* chan4_t) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description:
    *
    *   Status: INCOMPLETE
    *
    *   Notes: NONE! 
    */

    return;
}

void SyncPinger (double* chan1_t, int* pingerSynced, double* PRT) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description:
    *
    *   Status: INCOMPLETE
    *
    *   Notes: NONE! 
    */

    return;
}
