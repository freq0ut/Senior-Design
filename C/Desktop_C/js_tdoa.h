#ifndef JS_TDOA_H
#define JS_TDOA_H

#include "js_cmath.h"
#include "js_cfile.h"

void AdjustPGA (double* f, double _Complex* chanx_f, double powerMin, double powerMax);
void CalcDiamondPingerLocation(double D, double td2, double td3, double td4, double TOA1, double* pingerLocs);
void CalcPingerAzimuths(double* pingerLocs, double* azimuthH, double* azimuthV1, double* azimuthV2);
double CalcTimeDelay (double* chan1_t, double* chanx_t, double TOA1);
void CenterWindow (double* chan1_t, double* PRT, double *TOA1);
void DelaySampleTrigger(double PRT);
void SampleAllChans (int simFrame, double* chan1_t, double* chan2_t, double* chan3_t, double* chan4_t);
void SyncPinger (int simFrame, double* chan1_t, double* chan2_t, double* chan3_t, double* chan4t, double _Complex H, double* PRT, 
	int* pingerSynced);

#endif
