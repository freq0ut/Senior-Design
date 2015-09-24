#include "js_tdoa.h"

void PingerAzimuths (double* pingerLocation, double* azimuth1, double* azimuth2, double* azimuth3) {
    // Author: Joshua Simmons
    //
    // Date: September 2015
    //
    // Description: Computes the horizontal and vertical azimuths to the pinger.
    //
    // Status: UNTESTED
    //
    // Notes: NONE!

    // Determining pinger azimuths in radians
    *azimuth1 = atan2(   pingerLocation[2],pingerLocation[1]);
    *azimuth2 = atan2(   pingerLocation[3],pingerLocation[1]);
    *azimuth3 = atan2(-1*pingerLocation[3],pingerLocation[1]);

    // Wrapping angle [0,2pi]
    *azimuth1 = WrapTo2Pi(azimuth1);
    *azimuth2 = WrapTo2Pi(azimuth2);
    *azimuth3 = WrapTo2Pi(azimuth3);
    
    // Converting to degrees
    *azimuth1 = (*azimuth1) * (180.0/pi);
    *azimuth2 = (*azimuth2) * (180.0/pi);
    *azimuth3 = (*azimuth3) * (180.0/pi);

    return;
}

void PingerLocation (double d, double* sphericalRadii, double* pingerLocation) {
    // Author: Joshua Simmons
    //
    // Date: September 2015
    //
    // Description: Computes the estimated pseudo-location of the pinger.
    //
    // Status: UNTESTED
    //
    // Notes: I might have messed up the formula for the 3rd argument.

    pingerLocation[1] = (     pow(sphericalRadii[4],2) - pow(sphericalRadii[2],2) ) / (4.0*d);
    pingerLocation[2] = (     pow(sphericalRadii[3],2) - pow(sphericalRadii[1],2) ) / (4.0*d);
    pingerLocation[3] = sqrt( pow(sphericalRadii[2],2) - pow(sphericalRadii[1],2) - pow(sphericalRadii[2]-d,2) );

    return;
}

void SphereRadii (double td2, double td3, double td4, double toa, double vp, double* sphericalRadii) {
    // Author: Joshua Simmons
    //
    // Date: September 2015
    //
    // Description: Computes the spherical radii (multilateration).
    //
    // Status: UNTESTED
    //
    // Notes: NONE!

    sphericalRadii[1] = vp*(toa);
    sphericalRadii[2] = vp*(toa+td2);
    sphericalRadii[3] = vp*(toa+td3);
    sphericalRadii[4] = vp*(toa+td4);

    return;
}

