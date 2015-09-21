#include "js_tdoa.h"

void PingerAzimuths (void) {
    // Status: UNTESTED
    //
    // Description: Computes the horizontal and vertical azimuths to the pinger.

    azimuthH  = atan(   pingerLoc[2],pingerLoc[1]);
    azimuthV1 = atan(   pingerLoc[3],pingerLoc[1]);
    azimuthV2 = atan(-1*pingerLoc[3],pingerLoc[1]);

    // Wrapping angle [0,2pi]
    azimuthH  = WrapTo2Pi(azimuthH);
    azimuthV1 = WrapTo2Pi(azimuthV1);
    azimuthV2 = WrapTo2Pi(azimuthV2);
    
    // Converting to Degrees
    azimuthH  = azimuthH  * (180.0/pi);
    azimuthV1 = azimuthV1 * (180.0/pi);
    azimuthV2 = azimuthV2 * (180.0/pi);

    return;
}

void PingerLocation (void) {
    // Status: UNTESTED
    //
    // Description: Computes the estimated pseudo-location of the pinger.

    PingerLoc[1] = (R[4]*R[4]-R[2]*R[2])/(4.0*d);
    PingerLoc[2] = (R[3]*R[3]-R[1]*R[1])/(4.0*d);
    PingerLoc[3] = sqrt(R[1]*R[1]-PingerLoc[1]*PingerLoc[1]-(PingerLoc[2]-d)*(PingerLoc[2]-d));

    return;
}

void SphereRadii (void) {
    // Status: UNTESTED
    //
    // Description: Computes the spherical radii.

    R[1] = vP*(TOA);
    R[2] = vP*(TOA+tD2);
    R[3] = vP*(TOA+tD3);
    R[4] = vP*(TOA+tD4);

    return;
}

