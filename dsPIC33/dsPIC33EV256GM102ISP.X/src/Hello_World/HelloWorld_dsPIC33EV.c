/*
*   Joshua Simmons
*
*   January 21, 2016
*
*   Turns an LED on and off using a delay.
*/

#include "config.h"

void ms_delay(int N);

int main(void) {

    /*******************************************************
    ********** OSCILLATOR CHANGE OVER FOR 70 MIPS **********
    *******************************************************/

    // For use with 20 MHz external oscillator only!
    PLLFBD = 54; // M = 56, PLL Feedback Divisor
    CLKDIVbits.PLLPOST = 2; // N1 = 4, PLL Pre-scalar
    CLKDIVbits.PLLPRE = 0; // N2 = 2, PLL Post-scalar

    // Initiate Clock Switch to Primary Oscillator with PLL (NOSC = 0b011)
    __builtin_write_OSCCONH(0x03);
    __builtin_write_OSCCONL(0x01);

    // Wait for Clock switch to occur
    while (OSCCONbits.COSC != 0b011);

    // Wait for PLL to lock
    while(OSCCONbits.LOCK != 1) {};

    /**************************************
    ********** I/O CONFIGURATION **********
    **************************************/

    TRISB = 0x0000;// Output configured
    
    /******************************************
    ********** START OF MAIN PROGRAM **********
    ******************************************/

    while (1) {
        PORTBbits.RB8 = 0;
        ms_delay(1);
        PORTBbits.RB8 = 1;
        ms_delay(1);
    }
    
    return 0;
}

/*******************************************
********** SUPPLEMENTAL FUNCTIONS **********
*******************************************/

void ms_delay(int N) {
/*
*   This function implements a 16-bit timer to delay for N milliseconds.
*
*   TIMER = fOSC * tD / (2*P)
*       fOSC = Oscillator Frequency [Hz]
*       tD = Time Delay [s]
*       P = Pre-scalar (1,8,64,256)
*
*   TIMER = fOSC * N * 1E-3 / (2*P)
*       N = integer
*
*   Let fOSC = 140 [MHz], P = 256
*
*   TIMER = 140E+6 * N * 1E-3 / (2*256)
*   TIMER = 273.4375 * N
*
*   TIMER_MAX = 2^16-1 = 273.4375 * N_MAX
*   N_MAX = 239
*/

    int TIMER;
    T1CON = 0x8030;

    if ( N > 0 && N <= 239 ) {
        TIMER = (int) 273.4 * N;
        TMR1 = 0;//Reset timer

        while ( TMR1 < TIMER ) {}
    }
    else;

    return;
}
