#include "js_tdoa.h"

void AdjustPGA (double* f, double _Complex* chanx_f, double powerMin, double powerMax) {
    /*
    * Author: Joshua Simmons
    *
    * Date: September 2015
    *
    * Description: Calculates the signal power and adjusts the gain of the ADC PGA accordingly.
    *
    * Status: INCOMPLETE
    *
    * Notes:
    *   1. Will have to wait until familarization with ADC to be implemented.
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
    * Author: Joshua Simmons
    *
    * Date: September 2015
    *
    * Description:
    *
    * Status: INCOMPLETE
    *
    * Notes:
    *   1. 
    */

    /*
    // Determining pinger azimuths in radians
    *azimuth1 = atan2(   pingerLocation[2],pingerLocation[1]);
    *azimuth2 = atan2(   pingerLocation[3],pingerLocation[1]);
    *azimuth3 = atan2(-1*pingerLocation[3],pingerLocation[1]);

    // Wrapping angle [0,2pi]
    *azimuth1 = WrapTo2Pi(*azimuth1);
    *azimuth2 = WrapTo2Pi(*azimuth2);
    *azimuth3 = WrapTo2Pi(*azimuth3);
    
    // Converting to degrees
    *azimuth1 = (*azimuth1) * (180.0/pi);
    *azimuth2 = (*azimuth2) * (180.0/pi);
    *azimuth3 = (*azimuth3) * (180.0/pi);
    */

    return;
}

void CalcDiamondPingerLocation(double td2, double td3, double td4, double TOA, double vP, double* pingerLocs){
    /*
    * Author: Joshua Simmons
    *
    * Date: September 2015
    *
    * Description:
    *
    * Status: INCOMPLETE
    *
    * Notes:
    *   1. 
    */

    /*
    sphericalRadii[1] = vp*(toa);
    sphericalRadii[2] = vp*(toa+td2);
    sphericalRadii[3] = vp*(toa+td3);
    sphericalRadii[4] = vp*(toa+td4);

    pingerLocation[1] = (     pow(sphericalRadii[4],2) - pow(sphericalRadii[2],2) ) / (4.0*d);
    pingerLocation[2] = (     pow(sphericalRadii[3],2) - pow(sphericalRadii[1],2) ) / (4.0*d);
    pingerLocation[3] = sqrt( pow(sphericalRadii[2],2) - pow(sphericalRadii[1],2) - pow(sphericalRadii[2]-d,2) );
    */

    return;
}

double CalcTimeDelay (double* chan2_t, double THD, int lagBounds, int pkCounterMax) {
    /*
    * Author: Joshua Simmons
    *
    * Date: September 2015
    *
    * Description:
    *
    * Status: INCOMPLETE
    *
    * Notes:
    *   1. 
    */

    /*
    breakwall_TOAs[2] = Breakwall(chan2_t,THD,0) / fADC;
    breakwall_tDs[2] = breakwall_TOAs[1] - breakwall_TOAs[2];
    XCorrBounded(chan1_t,chan2_t,lagBounds, XCorr1x_Lags, XCorr1x);
    localMaxima(XC1x);// Where is the return argument going?
    td2 = Compare(breakwall_tDs[2], XC_tDs);
    double complex pingerLoc[1+3] = {3,0,0,0}; // Pinger Location Co-ordinates
    double breakwall_tDs [1+4] = {4,0,0,0,0}; // Breakwall Time Delays
    double breakwall_TOAs[1+4] = {4,0,0,0,0}; // Breakwall Time-of-Arrivals
    double sphereRadii[1+4] = {4,0,0,0};      // Sphere Radii (multilateration)
    int XCorr1x_Lags[1+(2*lagBounds+1)] = {2*lagBounds+1}; // Cross-correlation lags
    double XCorr1x  [1+(2*lagBounds+1)] = {2*lagBounds+1}; // Cross-correlations
    double XCorr_tDs[1+pkCounterMax] = {pkCounterMax};     // Cross-correlation Time Delays
    */

    return -1.0;
}

void CenterWindow (double* chan1_t, int* pingerSynced, double* PRT) {
    /*
    * Author: Joshua Simmons
    *
    * Date: September 2015
    *
    * Description:
    *
    * Status: INCOMPLETE
    *
    * Notes:
    *   1. 
    */

    // breakwall_TOAs[1] = Breakwall(chan1_t,threshold,0) / fADC;

    return;
}

void DelaySampleTrigger(double PRT) {
    /*
    * Author: Joshua Simmons
    *
    * Date: September 2015
    *
    * Description:
    *
    * Status: INCOMPLETE
    *
    * Notes:
    *   1. 
    */

    return;
}

void SampleAllChans (double fADC, double* chan1_t, double* chan2_t, double* chan3_t, double* chan4_t) {
    /*
    * Author: Joshua Simmons
    *
    * Date: September 2015
    *
    * Description:
    *
    * Status: INCOMPLETE
    *
    * Notes:
    *   1. 
    */

    return;
}

void SyncPinger (double* chan1_t, int* pingerSynced, double* PRT) {
    /*
    * Author: Joshua Simmons
    *
    * Date: September 2015
    *
    * Description:
    *
    * Status: INCOMPLETE
    *
    * Notes:
    *   1. 
    */

    return;
}
