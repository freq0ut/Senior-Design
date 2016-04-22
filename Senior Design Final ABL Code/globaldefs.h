#ifndef GLOBALDEFS_H
#define GLOBALDEFS_H

    ///////////////////////////////////
    ////////// PREPROCESSORS //////////
    ///////////////////////////////////

    #include <xc.h>
    #include <dsp.h>
    #include <stdint.h>

	/////////////////////////////////
	////////// DEFINITIONS //////////
	/////////////////////////////////

	#define SAMPLING_FREQUENCY_kHz 267
	#define SAMPLING_PERIOD_us 3.75

    #define DMA_COPY_DEPTH 1 // 1,2,4,8,16,32,64,128 ONLY.
    #define NUMBER_OF_SAMPLES 512 // Must be divisible by 2*DMA_COPY_DEPTH

	#define FFT_PRESCALAR 0x5A2 //haxor // 1/2

	#define PEAK_MAX 0x7CCC   // Automatic Gain Control. 0x7CCC =  0.975 in signed fractional format.
	#define TROUGH_MIN 0x8333 // Automatic Gain Control. 0x8333 = -0.975 in signed fractional format.

    #define MAX_PEAKS 50
	#define THRESHOLD 0.9//0x7333 // 0.928
	#define MINIMUM_PEAK_DISTANCE 1

	#define VELOCITY_OF_PROPAGATION_TIMES_SAMPLING_PERIOD 0.001275 // For TDOA.
	#define SENSOR_SPACING_MINOR 0.070710678 // For TDOA. D/sqrt(2)

    #define PWM_PERIOD 1302 // PERIOD = F_OSC / (F_PWM * Pre)
                             // 41667 = (400/3)*1E6 / (200*16)
    #define DUTY_CYCLE_MIN  4167 // Servo max clockwise position (0 degrees).
    #define DUTY_CYCLE_MAX 20833 // Servo max counterclockwise position (180 degrees).

    #define tick2usConv 0.26046666667

	////////////////////////////
	////////// MACROS //////////
	////////////////////////////

	// Don't forget to configure TRISx in main()

	#define SPI_CS1_HIGH LATEbits.LATE14 = 1
	#define SPI_CS1_LOW  LATEbits.LATE14 = 0

	#define SPI_CS2_HIGH LATEbits.LATE15 = 1
	#define SPI_CS2_LOW  LATEbits.LATE15 = 0

	#define SPI_CLK_HIGH LATBbits.LATB4 = 1
	#define SPI_CLK_LOW  LATBbits.LATB4 = 0

	#define SPI_DIN_HIGH LATAbits.LATA8 = 1
	#define SPI_DIN_LOW  LATAbits.LATA8 = 0

    /////////////////////////////////////////////
    ////////// PUT INTO PROGRAM MEMORY //////////
    /////////////////////////////////////////////

    extern const float atan_LUT[90] __attribute__ ((space(auto_psv)));

    //////////////////////////////////////
    ////////// GLOBAL VARIABLES //////////
    //////////////////////////////////////

    volatile uint16_t triggerdelay_us = 1; // When to start ADC sampling after a "Trigger Event"
    volatile int16_t iTrig = 0;

    char strL1[32];
    volatile int16_t peakLocs[MAX_PEAKS];

    volatile int16_t PGA_gain = 1; // Shared gain setting between both PGAs (1,2,3,4,5,6,7 only)

    volatile int16_t pingPongToggle = 0; // For toggling between each DMA buffer
    volatile int16_t iADC1, iADC2, iADC3; // For placing samples from DMA to RAM

    struct 
    {
        fractional ADC1CH0[DMA_COPY_DEPTH];
        fractional ADC1CH1[DMA_COPY_DEPTH];
        fractional ADC1CH2[DMA_COPY_DEPTH];
        fractional ADC1CH3[DMA_COPY_DEPTH];
    } ADC1_Ping_Buffer __attribute__((space(data), aligned(4*DMA_COPY_DEPTH*2)));

    struct 
    {
        fractional ADC1CH0[DMA_COPY_DEPTH];
        fractional ADC1CH1[DMA_COPY_DEPTH];
        fractional ADC1CH2[DMA_COPY_DEPTH];
        fractional ADC1CH3[DMA_COPY_DEPTH];
    } ADC1_Pong_Buffer __attribute__((space(data), aligned(4*DMA_COPY_DEPTH*2)));
    
    fractional ADC1CH0 [NUMBER_OF_SAMPLES];
    fractional ADC1CH1 [NUMBER_OF_SAMPLES];
    fractional ADC1CH2 [NUMBER_OF_SAMPLES];
    fractional ADC1CH3 [NUMBER_OF_SAMPLES];
    
    fractional XCORR_01 [2*NUMBER_OF_SAMPLES-1];
    fractional XCORR_02 [2*NUMBER_OF_SAMPLES-1];
    fractional XCORR_03 [2*NUMBER_OF_SAMPLES-1];

    volatile double ADC0_MAX;
    volatile double ADC1_MAX;
    volatile double ADC2_MAX;
    volatile double ADC3_MAX;
    
    int MaxHasBroken0;
    int MaxHasBroken1;
    int MaxHasBroken2;
    int MaxHasBroken3;
    
    int breakPoint0;
    int breakPoint1;
    int breakPoint2;
    int breakPoint3;
    
    int TD0, TD1, TD2, TD3;
    
    double float_AVG_TRIG_P = 0;
    double ticks2us = 0;
    double initialPulseSetting;
    double O_X_Td;
    double O_Y_Td;
    
    unsigned long MSW_PR7_32;
    unsigned int MSW_PR7_16;
    unsigned long LSW_PR6_32;
    unsigned int LSW_PR6_16;
    
    int triggerEvent = 1; // Counter for triggers
    int sampleEvent = 0; // Counter for samples
    int vecMax_o, vecMin_o, vecMax_y, vecMin_y, vecMax_x, vecMin_x;
    int triggerActive = 0;
    int captureFronts = 0;
    int ADC_count = 0;
    int portionHigh = 0;
    
    int16_t timerTwo = 0;
    int16_t timerThree = 0;
    uint16_t PWMperiod = 1302;
    uint16_t PWMduty = 651;
    int16_t iMin, iMax;
    
    int16_t XCORR_01_Index;
    int16_t XCORR_02_Index;
    int16_t XCORR_03_Index;
    
    int16_t theta, phi;
    int16_t iDXC01, iDXC02, iDXC03;
    int16_t iLag1, iLag2, iLag3;
    
    double timeDelay1;
    double timeDelay2;
    double timeDelay3;
    
    uint32_t AVG_TRIG_PERIOD = 0;


    /////////////////////////////////////////
    ////////// FUNCTION PROTOTYPES //////////
    /////////////////////////////////////////

    void InitSystemClock140(void);
    void InitSystemClock133(void);
    void InitPWM(void);
    void InitADC1(void);
    void InitNoTriggerWatchdogTimer(void);
    void InitTrigger(void);
    void AutoGainControlAndCentering(void);
    
    void TDOA (void);
    int16_t Find_XC_Lag (fractional XCORR_0X[0]);
    void CalcDiamondTDOA_Angles (int16_t iDelay1, int16_t iDelay2, int16_t iDelay3, float d, float vP_Ts);
    
    void WritePGA(int16_t gain);
    void WriteServo(int16_t theta, int16_t phi);
    void GenPurpTimer32(void); 
    void InitC4TriggerTimer(void);
    void CaptureFrontTimer(void); 
    void ActivateADC(void);
    void UpdateLCD_AcqSig(void);
    void UpdateLCD_Sampling(void);

#endif
