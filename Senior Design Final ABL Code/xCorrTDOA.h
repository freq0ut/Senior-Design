#include <xc.h>
#include "globaldefs.h"
#include "jps_math.h"
#include "atanLUT.h"
#include "twiddleFactors.h"
#include "math.h"

void TDOA (void) 
{
    /////////////////////////////////////////////////////////////
    ////////// SCALING TO PREVENT ACCUMULATOR OVERFLOW //////////
    /////////////////////////////////////////////////////////////

    VectorScale(NUMBER_OF_SAMPLES, &ADC1CH0[0], &ADC1CH0[0], FFT_PRESCALAR);
    VectorScale(NUMBER_OF_SAMPLES, &ADC1CH1[0], &ADC1CH1[0], FFT_PRESCALAR);
    VectorScale(NUMBER_OF_SAMPLES, &ADC1CH2[0], &ADC1CH2[0], FFT_PRESCALAR);
    // VectorScale(NUMBER_OF_SAMPLES, &ADC1CH3[0], &ADC1CH3[0], FFT_PRESCALAR);
    
    ////////////////////////////////////////
    ////////// CROSS-CORRELATIONS //////////
    ////////////////////////////////////////

    VectorConvolve(NUMBER_OF_SAMPLES,NUMBER_OF_SAMPLES,&XCORR_01[0],&ADC1CH0[0],&ADC1CH1[0]); // O to X
    VectorConvolve(NUMBER_OF_SAMPLES,NUMBER_OF_SAMPLES,&XCORR_02[0],&ADC1CH0[0],&ADC1CH2[0]); // O to Y
    // VectorConvolve(NUMBER_OF_SAMPLES,NUMBER_OF_SAMPLES,&XCORR_03[0],&ADC1CH0[0],&ADC1CH3[0]);
    
    //////////////////////////////////////////////
    ////////// FINDING THE CORRECT LAGS //////////
    //////////////////////////////////////////////
    
    VectorMax(2*NUMBER_OF_SAMPLES-1, &XCORR_01[0], &XCORR_01_Index); // X Time Delay
    VectorMax(2*NUMBER_OF_SAMPLES-1, &XCORR_02[0], &XCORR_02_Index); // Y Time Delay
    
    //Calculate Time delay X and Time delay Y
    
    O_X_Td = (XCORR_01_Index - 512)*SAMPLING_PERIOD_us;
    O_Y_Td = (XCORR_02_Index - 512)*SAMPLING_PERIOD_us;

    //Calculate Angle X and Angle Y
    

    //Calculate ArcTan
    
    
    // VectorMax(2*NUMBER_OF_SAMPLES-1, &XCORR_03[0], &XCORR_03_Index);
    
    asm("nop");
    asm("nop");
    asm("nop");
    PGA_gain = PGA_gain;
    asm("nop");
    asm("nop");
    asm("nop"); 
}

void CalcDiamondTDOA_Angles (int16_t iDelay1, int16_t iDelay2, int16_t iDelay3, float d, float vP_Ts) {

    int16_t i;

    float iTOA0;
    float R0, R1, R2, R3;
    float Sx, Sy, Sz, Sz_Squared;
    float x0;

    iTOA0 = iDelay2-iDelay1-iDelay3;

    // 1 over zero protection and fudge factoring
    if ( iTOA0 == 0 ) { iTOA0++; }
    else;

    //R0_Actual = vP_Ts * (iDelay1*iDelay3)/iTOA0; // Necessary if servos not co-located with sensor array

    iTOA0  = (iDelay1*iDelay1+iDelay3*iDelay3-iDelay2*iDelay2) / (2*iTOA0);

    R0 = vP_Ts *  iTOA0;
    R1 = vP_Ts * (iTOA0+iDelay1);
    R2 = vP_Ts * (iTOA0+iDelay2);
    R3 = vP_Ts * (iTOA0+iDelay3);

    Sx = ((R3*R3 - R1*R1)/d)/4;
    Sy = ((R2*R2 - R0*R0)/d)/4;

    Sz_Squared = ( (R0*R0) - (Sx*Sx) - ((Sy-d)*(Sy-d)) );

    // Imaginary Number Protection
    if (Sz_Squared>0) {

        x0 = Sz_Squared; // Seed

        // Square Root Finding Algorithm (Newton-Raphson Method)
        for (i=0; i < 5; i++) {
            Sz = x0 - ((x0*x0-Sz_Squared)/x0)/2;
            x0 = Sz;
        }

        // Inverse Tangent Look-Up-Table
        theta = JPS_ArcTangent2(Sy,Sx); // Angle will be in degrees of data type int16_t
        phi   = JPS_ArcTangent2(Sz,Sx); // Angle will be in degrees of data type int16_t

        //WriteServo(theta,phi);

    }
    else{};
}