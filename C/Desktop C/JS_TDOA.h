#ifndef JS_TDOA_H
#define JS_TDOA_H

#include "js_cmath.h"

//BreakWall_tDs();
//BreakWall_TOAs();
//Compare_tDs(
void PingerLocation (double d, double* sphericalRadii, double* pingerLocation);
void PingerAzimuths (double* pingerLocation, double* azimuth1, double* azimuth2, double* azimuth3);
void SphereRadii (double td2, double td3, double td4, double toa, double vp, double* sphericalRadii);
//Superior_TOA();
//syncPinger();
//XC_tDs();

#endif
