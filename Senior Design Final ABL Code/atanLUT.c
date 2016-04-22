#include "atanLUT.h"

float atan_LUT[90] __attribute__ ((space(auto_psv))) = {
	0.000000,//  0 degrees
	0.017455,//  1 degrees
	0.034921,//  2 degrees
	0.052408,//  3 degrees
	0.069927,//  4 degrees
	0.087489,//  5 degrees
	0.105104,//  6 degrees
	0.122785,//  7 degrees
	0.140541,//  8 degrees
	0.158384,//  9 degrees
	0.176327,// 10 degrees
	0.194380,// 11 degrees
	0.212557,// 12 degrees
	0.230868,// 13 degrees
	0.249328,// 14 degrees
	0.267949,// 15 degrees
	0.286745,// 16 degrees
	0.305731,// 17 degrees
	0.324920,// 18 degrees
	0.344328,// 19 degrees
	0.363970,// 20 degrees
	0.383864,// 21 degrees
	0.404026,// 22 degrees
	0.424475,// 23 degrees
	0.445229,// 24 degrees
	0.466308,// 25 degrees
	0.487733,// 26 degrees
	0.509525,// 27 degrees
	0.531709,// 28 degrees
	0.554309,// 29 degrees
	0.577350,// 30 degrees
	0.600861,// 31 degrees
	0.624869,// 32 degrees
	0.649408,// 33 degrees
	0.674509,// 34 degrees
	0.700208,// 35 degrees
	0.726543,// 36 degrees
	0.753554,// 37 degrees
	0.781286,// 38 degrees
	0.809784,// 39 degrees
	0.839100,// 40 degrees
	0.869287,// 41 degrees
	0.900404,// 42 degrees
	0.932515,// 43 degrees
	0.965689,// 44 degrees
	1.000000,// 45 degrees
	1.035530,// 46 degrees
	1.072369,// 47 degrees
	1.110613,// 48 degrees
	1.150368,// 49 degrees
	1.191754,// 50 degrees
	1.234897,// 51 degrees
	1.279942,// 52 degrees
	1.327045,// 53 degrees
	1.376382,// 54 degrees
	1.428148,// 55 degrees
	1.482561,// 56 degrees
	1.539865,// 57 degrees
	1.600335,// 58 degrees
	1.664279,// 59 degrees
	1.732051,// 60 degrees
	1.804048,// 61 degrees
	1.880726,// 62 degrees
	1.962610,// 63 degrees
	2.050304,// 64 degrees
	2.144507,// 65 degrees
	2.246037,// 66 degrees
	2.355852,// 67 degrees
	2.475087,// 68 degrees
	2.605089,// 69 degrees
	2.747478,// 70 degrees
	2.904211,// 71 degrees
	3.077683,// 72 degrees
	3.270853,// 73 degrees
	3.487414,// 74 degrees
	3.732051,// 75 degrees
	4.010781,// 76 degrees
	4.331476,// 77 degrees
	4.704630,// 78 degrees
	5.144554,// 79 degrees
	5.671282,// 80 degrees
	6.313752,// 81 degrees
	7.115370,// 82 degrees
	8.144346,// 83 degrees
	9.514364,// 84 degrees
	11.430053,// 85 degrees
	14.300666,// 86 degrees
	19.081137,// 87 degrees
	28.636253,// 88 degrees
	57.289963// 89 degrees
};


int16_t JPS_ArcTangent2 (float y, float x) {

    int16_t i;
    int16_t angle = -1;
    float tangent;

    if 		(x >  0 && y == 0) { angle =   0; }
    else if (x == 0 && y >  0) { angle =  90; }
    else if (x <  0 && y == 0) { angle = 180; }
    else if (x == 0 && y <  0) { angle = 270; }
    else {

    	tangent = y/x;

	    // This value must be positive because of the LUT
	    if (tangent < 0) { tangent = -tangent; }
	    else;

	    // Round up to the nearest whole number index. Max error is half a degree.
	    for (i=1; i<90; i++) {
	        if (tangent < atan_LUT[i]) {
	        	angle = i;
	        	break;
	        }
	        else;
	    }

	    // Adjusting Angle
	    if      (x < 0 && y >= 0) { angle = 360 - angle; } // 2nd Quadrant
	    else if (x > 0 && y <  0) { angle = 360 - angle; } // 4th Quadrant
	    else;
	}

    return (angle);
}
