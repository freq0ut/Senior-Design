// DSPIC33EP512GM306 Configuration Bit Settings

// 'C' source line config statements



#include <xc.h>

// FICD
#pragma config ICS = PGD2               // ICD Communication Channel Select bits (Communicate on PGEC2 and PGED2)
#pragma config JTAGEN = OFF             // JTAG Enable bit (JTAG is disabled)

// FPOR
#pragma config BOREN = ON               // Brown Out Reset enabled
#pragma config ALTI2C1 = OFF            // Alternate I2C1 pins (I2C1 mapped to SDA1/SCL1 pins)
#pragma config ALTI2C2 = OFF            // Alternate I2C2 pins (I2C2 mapped to SDA2/SCL2 pins)
#pragma config WDTWIN = WIN25           // Watchdog Window Select bits (WDT Window is 25% of WDT period)

// FWDT
#pragma config WDTPOST = PS32768        // Watchdog Timer Postscaler bits (1:32,768)
#pragma config WDTPRE = PR128           // Watchdog Timer Prescaler bit (1:128)
#pragma config PLLKEN = ON              // PLL Lock Enable bit (Clock switch to PLL source will wait until the PLL lock signal is valid.)
#pragma config WINDIS = OFF             // Watchdog Timer Window Enable bit (Watchdog Timer in Non-Window mode)
#pragma config FWDTEN = OFF              // Watchdog Timer Enable bit (Watchdog timer always enabled)

// FOSC
#pragma config POSCMD = HS              // Primary Oscillator Mode Select bits (HS Crystal Oscillator Mode)
#pragma config OSCIOFNC = OFF           // OSC2 Pin Function bit (OSC2 is clock output)
#pragma config IOL1WAY = ON             // Peripheral pin select configuration (Allow only one reconfiguration)
#pragma config FCKSM = CSECMD           // Clock Switching Mode bits (Clock switching is enabled,Fail-safe Clock Monitor is disabled)

// FOSCSEL
#pragma config FNOSC = FRCDIVN          // Oscillator Source Selection (Internal Fast RC (FRC) Oscillator with postscaler)
#pragma config PWMLOCK = ON             // PWM Lock Enabled
#pragma config IESO = ON                // Two-speed Oscillator Start-up Enable bit (Start up device with FRC, then switch to user-selected oscillator source)

// FGS
#pragma config GWRP = OFF               // General Segment Write-Protect bit (General Segment may be written)
#pragma config GCP = OFF                // General Segment Code-Protect bit (General Segment Code protect is Disabled)

#ifndef __DELAY_H

	//#define FOSC        140000000LL
    //#define TCYC        14.2857e-9
    #define FOSC		133333332LL		// clock-frequecy in Hz with suffix LL (64-bit-long), eg. 32000000LL for 32MHz
    #define TCYC        15e-9
    #define FCY      	(FOSC/2)		// MCU is running at FCY MIPS

    #define delay_us(x)	__delay32(((x*FCY)/1000000L))	// delays x us
    #define delay_ms(x)	__delay32(((x*FCY)/1000L))		// delays x ms

	#define __DELAY_H	1

#endif

/*
 * Notes:
 *
 *      Device
 *          dsPIC33EP512GM306
 *          16 MHz oscillator
 *          200/3 MIPs ~ 67 MIPs
 *
 *      Programmable Gain Amplifiers
 *          Linear Technology LTC6912-2
 *              1 => G=1
 *              2 => G=2
 *              3 => G=4
 *              4 => G=8
 *              5 => G=16
 *              6 => G=32
 *              7 => G=64
 *
 *      Comparator 5 Setup
 *          Comparator reference voltage set to 75% of AVDD 
 *          Pin 12 (C5IN3-/RA4) = Trigger signal 
 *
 *      Sample Countdown Timer (TMR1)
 *          16-bit Timer
 *          Prescalar = 8
 *          Maximum Delay = 7.8642 ms
 *
 *      ADC1 Setup
 *          10-bit, simultaneous, quad channel, TAD = 75 ns
 *          Pin 2 (AN0/RA0) = ADC1CH0
 *          Pin 3 (AN1/RA1) = ADC1CH1
 *          Pin 4 (AN2/RB0) = ADC1CH2
 *          Pin 5 (AN3/RB1) = ADC1CH3
 *
 *      Header Files
 *          dsp.h : fractional, fractcomplex
 *          stdint.h : int16_t, uint16_t
 * 
 *      Library Files
 *          /opt/microchip/xc16/v1.25/lib/libdsp-elf-a
 */