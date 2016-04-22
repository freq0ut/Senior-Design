#include <stdio.h>
#include <libpic30.h>
#include "config.h"
#include "globaldefs.h"
#include "jps_math.h"
#include "atanLUT.h"
#include "LCD.h"
#include "ADC.h"
#include "ISRs.h"
#include "PGAControl.h"
#include "PWM.h"
#include "timingTriggers.h"
#include "xCorrTDOA.h"

///////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// MAIN FUNCTION ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
int main(void) 
{
    ///////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////// SETUP ////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
    // Analog/Digital Pin Configuration
    ANSELA = 0x0000;
    ANSELB = 0x0000;
    ANSELC = 0x0000;
    ANSELD = 0x0000;
    ANSELG = 0x0000;
    ANSELCbits.ANSC1 = 0; // Pin 22 (RC1/AN7) is analog (Trigger)
    ANSELAbits.ANSA0 = 1; // Pin 13 (RA0/AN0) is analog (ADC1CH0)
    ANSELAbits.ANSA1 = 1; // Pin 14 (RA1/AN1) is analog (ADC1CH1)
    ANSELBbits.ANSB0 = 1; // Pin 15 (RB0/AN2) is analog (ADC1CH2)
    ANSELBbits.ANSB1 = 1; // Pin 16 (RB1/AN3) is analog (ADC1CH3)
    
    __builtin_write_OSCCONL(OSCCON & 0xBF); // unlock registers to configure PPS
    RPOR7bits.RP57R = 0b10000;   // configure RP57 as OC1 (PWM Horizontal)
    __builtin_write_OSCCONL(OSCCON | 0x40); // lock registers
    
    // Digital I/O Pin Configuration
    //TRISCbits.TRISC9 = 0;
    
    TRISEbits.TRISE14 = 0; // Pin 29 (RE14) is a digital output (Chip Select 1)
    SPI_CS1_HIGH;

    TRISEbits.TRISE15 = 0; // Pin 30 (RE15) is a digital output (Chip Select 2)
    SPI_CS2_HIGH;

    TRISBbits.TRISB4 = 0; // Pin 32 (RB4) is a digital output (SPI Clock)
    SPI_CLK_LOW;

    TRISAbits.TRISA8 = 0; // Pin 31 (RA8) is a digital output (SPI Data In)
    SPI_DIN_LOW;
    
    TRISGbits.TRISG6 = 0; // Register Select Output for LCD.

    // Initializing Main Loop
    InitSystemClock133(); // Make sure this matches in the config.h (FOSC 133333332LL)
    //InitSystemClock140(); // Make sure this matches in the config.h (FOSC 140000000LL)
    _NSTDIS = 0;
    GenPurpTimer32();
    InitNoTriggerWatchdogTimer();
    InitC4TriggerTimer();
    InitTrigger();
    CaptureFrontTimer();
    
    InitPWM();
    initPMP();
    initLCD();
    InitADC1();
    clearLCD();
    
    WritePGA(PGA_gain); // Initializes both PGAs
    //homeScreen();
    WriteServo(90,90); // Initialize servo positions (whole number degrees from 0 to 180))
    
    while(1) 
    {
    	// Infinite Loop:
        WriteServo(90,90);
        delay_ms(1000);
        WriteServo(90+35,90);
        delay_ms(1000);
        WriteServo(90,90);
        delay_ms(1000);
        WriteServo(90-35,90);
        delay_ms(1000);
    }
}

