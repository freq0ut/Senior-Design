#ifndef JS_TDOA_H
#define JS_TDOA_H

#include "js_cmath.h"
#include "js_cfile.h"

void AdjustPGA (double* f, double _Complex* chanx_f, double powerMin, double powerMax);
void CalcDiamondPingerLocation(double d, double td2, double td3, double td4, double TOA1, double vP, double* pingerLocs);
void CalcPingerAzimuths(double* pingerLocs, double* azimuthH, double* azimuthV1, double* azimuthV2);
double CalcTimeDelay (double* chan1_t, double* chanx_t, double threshold, double TOA1, int lagBounds, int pkCounterMax);
void CenterWindow (double* chan1_t, int N0, double fADC, double threshold, int* pingerSynced, double* PRT, double* TOA1);
void DelaySampleTrigger(double PRT);
void SampleAllChans (double fADC, double* chan1_t, double* chan2_t, double* chan3_t, double* chan4_t);
void SyncPinger (double* chan1_t, int* pingerSynced, double* PRT);

#endif
