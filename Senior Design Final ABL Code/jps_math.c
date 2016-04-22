#include "jps_math.h"

float JPS_Average (int16_t N, fractional srcV[]) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Computes the arithmetic mean.
    *
    *   Status: untested
    *
    *   Notes: NONE!
    */

    int16_t i;
    float average = 0;

    for (i=0; i < N; i++) {
        average += srcV[i];
    }

    average /= N;

    return (average);
}


void JPS_BubbleSort (int16_t N, fractional srcV[]) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Sorts an array using the Bubble Sort algorithm.
    *   
    *   Status: untested
    *               
    *   Notes: NONE!
    */

    int16_t i;

    for (i=0; i < N; i++) {
        if ( srcV[i] > srcV[i+1] ) {
            srcV[i]   ^= srcV[i+1];
            srcV[i+1] ^= srcV[i];
            srcV[i]   ^= srcV[i+1];
        }
        else;
    }

    return;
}


void JPS_CircularShift (int16_t N, fractional srcV[]) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: April 2016
    *
    *   Description: Swaps the left half with the right half of an array.
    *   
    *   Status: untested
    *               
    *   Notes: NONE!
    */

    int16_t i;

    for(i=0; i<N/2; i++) {
        srcV[i]     ^= srcV[N/2+i];
        srcV[N/2+i] ^= srcV[i];
        srcV[i]     ^= srcV[N/2+i];
    }

    return;
}


void JPS_Convolve (int16_t N, fractional srcV1[], fractional srcV2[], fractional dstV[]) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: April 2016
    *
    *   Description: Convolution of two vectors srcV1 and srcV2. They must both be the same length.
    *
    *   Status: Good to go!
    *
    *   Notes: NONE!
    */

    int16_t k;
    int16_t n;

    for (n=0; n < 2*N-1; n++) {
        for (k=-N+1+n; k <= n; k++) {
            if ( (k >= 0) && (k <= N-1) ) {
                dstV[n] += srcV1[k] * srcV2[n-k];
            }
            else;
        }
    }

    return;
}


int16_t JPS_Max (fractional srcV[], int16_t iLower, int16_t iUpper) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Returns the index corresponding to the maximum value of an array.
    *                   If no max is found then -1 is returned. 
    *
    *   Status: untested
    *
    *   Notes: NONE!
    */

    int16_t i;
    int16_t iMax = -1;
    fractional max = 0;

    for (i=iLower; i < iUpper; i++) {
        if (srcV[i]>max) {
            iMax = i;
            max = srcV[i];
        }
        else;
    }

    return (iMax);
}


fractional JPS_Median (int16_t N, fractional srcV[]) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Sorts from smallest-to-largest an array Y then returns the median.
    *
    *   Status: untested
    *
    *   Notes:
    *       1. The array Y is OVERWRITTEN.
    *       2. The array Y is assumed to be small since Bubble Sort is used.
    */

    fractional median;

    JPS_BubbleSort(N, &srcV[0]);

    if ( N%2 == 0 ) {
        median = srcV[N/2]*0x4000 + srcV[N/2+1]*0x4000;
    }
    else {
        median = srcV[N/2+1];
    }

    return (median);
}


void JPS_PeaksFinder (int16_t N, fractional* srcV, int16_t Npks, int16_t* dstV, fractional THR, int16_t MPD) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: April 9, 2016
    *
    *   Description: Finds the peaks of a signal.
    *                   THR = Threshold
    *                   MPD = Minimum Peak Distance
    *
    *   Status: untested
    *
    *   Notes: NONE!
    */

    int16_t i;
    int16_t iPeak = 0;
    int16_t iMax;
    fractional max = 0x0000;
    fractional yTHR;

    // Initialize peakLocs array to zero
    for (i=0; i < Npks; i++) {
        dstV[i] = 0;
    }
    
    // Find the maximum value of the array
    max = VectorMax (N, &srcV[0], &iMax);
    
    // Store the max as the zeroth peak
    dstV[iPeak] = iMax;
    iPeak++;

    VectorScale(1, &yTHR, &max, THR);
    
    // Search for peaks to the left of iMax
    for (i=iMax-MPD; i > 0; i--) { 
        if ( (srcV[i] > yTHR) && (srcV[i] > srcV[i-1]) && (srcV[i] > srcV[i+1]) && (iPeak < Npks)) {
                dstV[iPeak] = i;
                iPeak++;
                i = i - MPD + 1;
        }
        else;
    }
    
    // Search for peaks to the right of iMax
    for (i=iMax+MPD; i < N-1; i++) { 
        if ( (srcV[i] > yTHR) && (srcV[i] > srcV[i-1]) && (srcV[i] > srcV[i+1]) && (iPeak < Npks) ) {
                dstV[iPeak] = i;
                iPeak++;
                i = i + MPD - 1;
                
        }
        else;
    }

    return;
}


float JPS_SignalPower (int16_t N, fractional srcV[]) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Computes the real power of signal Y.
    *
    *   Status: untested
    *
    *   Notes: NONE!
    */

    int16_t i;
    float power = 0;

    for (i=0; i < N; i++) {
        power += srcV[i]*srcV[i];
    }

    power /= N;

    return (power);
}


float JPS_Simpson (int16_t N, fractional srcV[]) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Computes the area under the curve numerically using Simpsons Rule.
    *
    *   Status: untested
    *
    *   Notes:
    *       1. dx is assumed to be constant.
    *       2. The lengths of X and Y must be even.
    *       3. You will need to multiply the result by deltaX.
    */

    int i;
    float A;
    float Seven = 0;
    float Sodd = 0;

    for (i=1; i < N/2-1; i+=2) {
        Seven = Seven + srcV[i];
        Sodd  = Sodd  + srcV[i+1];
    }

    A = ( srcV[0] + 4*Seven + 2*Sodd + srcV[N-1] ) / 3;

    return (A);
}


float JPS_SquareRoot (float xSquared) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Computes the square root via Newton-Raphson Method.
    *
    *   Status: untested
    *
    *   Notes: NONE!
    */

    int16_t i;
    float x0 = xSquared; // Seed
    float x1;

    for (i=0; i < 10; i++) {
        x1 = x0 - ((x0*x0-xSquared)/x0)/2;
        x0 = x1;
    }

    return (x1);
}


float JPS_Trapz (int16_t N, fractional srcV[]) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: September 2015
    *
    *   Description: Computes the area under the curve numerically using Trapezoidal Rule.
    *
    *   Status: untested
    *
    *   Notes:
    *       1. deltaX is assumed to be constant.
    *       2. You will need to multiply the result by deltaX.
    */

    int16_t  i;
    float A;
    float Smid = 0;

    for (i=1; i < N-1; i++) {
        Smid += srcV[i];
    }

    A = 0.5 * ( srcV[0] + 2*Smid + srcV[N-1] );

    return (A);
}


void JPS_XCorr (int16_t N, fractional srcV1[], fractional srcV2[], fractional dstV[]) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: November 2015
    *
    *   Description: Computes the non-normalized cross-correlation of two arrays.
    *
    *   Status: untested
    *
    *   Notes: NONE!
    */

    int16_t k;
    int16_t n;

    // Vector Reverse then call convolve

    for (n=0; n < 2*N-1; n++) {
        for (k=-N+1+n; k <= n; k++) {
            if ( (k >= 0) && (k <= N-1) ) {
                dstV[n] += srcV1[k] * srcV2[-n+k+N-1];
            }
            else;
        }
    }

    return;
}


void JPS_CmpxAdd (int16_t N, fractcomplex srcV[], fractcomplex dstV[]) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: April 2016
    *
    *   Description: Adds two complex arrays.
    *   
    *   Status: untested
    *               
    *   Notes: NONE!
    */

    int16_t i;

    for (i=0; i < N; i++) {
        dstV[i].real += srcV[i].real;
        dstV[i].imag += srcV[i].imag;
    }

    return;
}



void JPS_CmpxMultiply(int16_t N, fractcomplex* srcV, fractcomplex* dstV) 
{
    int16_t i;
    //int16_t oldCORCON;
    fractcomplex temp;
    //oldCORCON = CORCON;
    //CORCON = 0x0024;
    for (i=0; i<N; i++) 
    {
        temp.imag = srcV[i].real*dstV[i].imag + srcV[i].imag*dstV[i].real;
        
        dstV[i].real = srcV[i].real*dstV[i].real - srcV[i].imag*dstV[i].imag;
        dstV[i].imag = temp.imag;
    }
    //CORCON = oldCORCON;
    return;
}


void JPS_CmpxPeaksFinder (int16_t N, fractional srcV[], int16_t Npks, int16_t dstV[], fractional THR, int16_t MPD) {
    /*
    *   Author: Joshua Simmons
    *
    *   Date: April 9, 2016
    *
    *   Description: Finds the peaks of a signal.
    *                   THR = Threshold
    *                   MPD = Minimum Peak Distance
    *
    *   Status: untested
    *
    *   Notes: NONE!
    */

    int16_t i;
    int16_t iPeak = 0;
    int16_t iMax;
    fractional max = 0x0000;
    fractional yTHR;

    // Initialize peakLocs array to zero
    for (i=0; i < Npks; i++) {
        dstV[i] = 0;
    }

    // Find the maximum y
    for (i=0; i < N; i++) {
        if (srcV[i]>max) {
            iMax = i;
            max = srcV[i];
        }
        else;
    }

    dstV[iPeak] = iMax;
    iPeak++; 

    VectorScale(1, &yTHR, &max, THR);
    
    // Search for peaks to the left of iMax
    for (i=iMax-MPD; i > 0; i--) { 
        if ( (srcV[i] > yTHR) && (srcV[i] > srcV[i-1]) && (srcV[i] > srcV[i+1]) && (iPeak < Npks)) {
                dstV[iPeak] = i;
                iPeak++;
                i = i - MPD + 1;
        }
        else;
    }

    // Search for peaks to the right of iMax
    for (i=iMax+MPD; i < N-1; i++) { 
        if ( (srcV[i] > yTHR) && (srcV[i] > srcV[i-1]) && (srcV[i] > srcV[i+1]) && (iPeak < Npks) ) {
                dstV[iPeak] = i;
                iPeak++;
                i = i + MPD - 1;
                
        }
        else;
    }
    return;
}
