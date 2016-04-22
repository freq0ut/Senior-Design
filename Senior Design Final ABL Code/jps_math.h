#ifndef JPS_CMATH_H
#define JPS_CMATH_H

	///////////////////////////////////
	////////// PREPROCESSORS //////////
	///////////////////////////////////

    #include <xc.h>
	#include <stdint.h>
    #include <dsp.h>

	////////////////////////////////
	////////// PROTOTYPES //////////
	////////////////////////////////

	float JPS_Average (int16_t N, fractional srcV[]);
	void JPS_BubbleSort (int16_t N, fractional srcV[]);
	void JPS_CircularShift (int16_t, fractional srcV[]);
	void JPS_Convolve (int16_t N, fractional srcV1[], fractional srcV2[], fractional dstV[]);
	int16_t JPS_Max (fractional srcV[], int16_t iLower, int16_t iUpper);
	fractional JPS_Median (int16_t N, fractional srcV[]);
	void JPS_PeaksFinder (int16_t N, fractional srcV[], int16_t Npks, int16_t dstV[], fractional THR, int16_t MPD);
	float JPS_SignalPower (int16_t N, fractional srcV[]);
	float JPS_Simpson (int16_t N, fractional srcV[]);
	float JPS_SquareRoot (float xSquared);
	float JPS_Trapz (int16_t N, fractional srcV[]);
	void JPS_XCorr (int16_t N, fractional srcV1[], fractional srcV2[], fractional dstV[]); // Needs changing

	void JPS_CmpxAdd (int16_t N, fractcomplex srcV[], fractcomplex dstV[]);
    
	void JPS_CmpxMultiply (int16_t N, fractcomplex srcV[], fractcomplex dstV[]);
    
	void JPS_CmpxPeaksFinder (int16_t N, fractional srcV[], int16_t Npks, int16_t dstV[], fractional THR, int16_t MPD);

#endif
