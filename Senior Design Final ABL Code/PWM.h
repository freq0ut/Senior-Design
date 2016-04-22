#include <xc.h>
#include "globaldefs.h"
#include "jps_math.h"
#include "atanLUT.h"

void InitPWM (void) 
{
    T3CON = 0x8010;     // enable timer 3, PreScaler of 8

                        // want low end to be 5 ms
                        // 65535 is full scale of 16 bit counter
    PR3 = PWMperiod-1;     // set the period register
    _T3IF = 0;          // set T3 interrupt flag low
    _T3IE = 1;          // enable the T3 interrupt
    
    OC1CON1 = 0x040E;    // activate PWM module (T3 is clock source, PWM mode on;
    OC1CON2 = 0x0000;
    OC1R = OC1RS = PWMduty;  // set initial duty cycle to 50%
    
    
                        // Fault pin, OCFx disabled.)
    //T3CONbits.TON = 1;
    TMR3 = 0;
}

void WriteServo(int16_t theta, int16_t phi) 
{
    // Duty Cycle for PWM6H (Horizontal Servo)
    if((theta >= 22) && (theta <= 158)) 
    {
        OC1R = theta*((DUTY_CYCLE_MAX - DUTY_CYCLE_MIN)/180) + DUTY_CYCLE_MIN;
    }
    else{};
    // Duty Cycle for PWM6H (Vertical Servo)
    
//    if((phi >= 75) &&   (phi <= 105)) 
//    {
//        PDC6 = phi*((DUTY_CYCLE_MIN - DUTY_CYCLE_MAX)/180) + DUTY_CYCLE_MAX;
//    }
//    else{};
}
