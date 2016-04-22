#include <xc.h>
#include "globaldefs.h"
#include "jps_math.h"
#include "atanLUT.h"

////////////////////////////////////////////////
////////// INTERRUPT SERVICE ROUTINES //////////
////////////////////////////////////////////////
#include "globaldefs.h"

// Interrupt for ADC Enabled (Timer7)
void __attribute__((__interrupt__, auto_psv)) _T7Interrupt(void)
{
    TMR6 = 0; // Reset Timer for Pulse Offset Capture.
    TMR7 = 0; // Reset Timer for Pulse Offset Capture.
    if(captureFronts == 1)
    {
        ActivateADC(); // Turn on ADC and update LCD with info.
        captureFronts = 0;
    }
    _T7IF = 0; // Clear TMR7 interrupt flag
}

// Interrupt for PWM Output Compare
void __attribute__((__interrupt__, auto_psv)) _T3Interrupt(void)
{
    OC1TMR = 0;
    TMR3 = 0;
    _T3IF = 0;  //set T3 flag back to zero 
}

// Comparator Interrupt Service Routine (Trigger)
void __attribute__ ((__interrupt__, auto_psv)) _CM1Interrupt(void) 
{
    TMR6 = 0; // Reset Timer for Pulse Offset Capture.
    TMR7 = 0; // Reset Timer for Pulse Offset Capture.
    
    TMR4 = 0; // Reset 32-bit TMR45 counter value.(used for loss of signal detection).  
    TMR5 = 0; // Reset 32-bit TMR45 counter value (used for loss of signal detection). 
    TMR1 = 0; // Reset timer to prevent numerous triggering within the pulse.
    
    triggerActive = 1; // Enable toggle bit to re-enable comparator in ISR.
    captureFronts = 1; // this allows ADC interrupt (timer 7) to capture leading edge

    /*CLEAR COMPARATOR INTERRUPT FLAG*/
    _CMIF = 0; 
}

// Interrupt for Comparator Event Clear (Timer1)
void __attribute__((__interrupt__, auto_psv)) _T1Interrupt(void)
{
    if(triggerActive == 1)
    {
        triggerActive = 0;
        CM4CONbits.CEVT = 0; // Clear CMP event status
        TMR1 = 0;
    }   
    _T1IF = 0;  //set T3 flag back to zero 
}

// ADC1 DMA0 Interrupt Service Routine (Sampling)
void __attribute__((__interrupt__, auto_psv)) _DMA0Interrupt(void) 
{

    if(pingPongToggle == 0) 
    {
        for(iADC2 = 0; iADC2 < DMA_COPY_DEPTH; iADC2++) 
        {
            ADC1CH0[iADC1+iADC2] = ADC1_Ping_Buffer.ADC1CH3[iADC2];
//            ADC1CH1[NUMBER_OF_SAMPLES-1-(iADC1+iADC2)] = ADC1_Ping_Buffer.ADC1CH2[iADC2];
//            ADC1CH2[NUMBER_OF_SAMPLES-1-(iADC1+iADC2)] = ADC1_Ping_Buffer.ADC1CH1[iADC2];
//            ADC1CH3[NUMBER_OF_SAMPLES-1-(iADC1+iADC2)] = ADC1_Ping_Buffer.ADC1CH0[iADC2];
            
            ADC1CH1[iADC1+iADC2] = ADC1_Ping_Buffer.ADC1CH2[iADC2];
            ADC1CH2[iADC1+iADC2] = ADC1_Ping_Buffer.ADC1CH1[iADC2];
            ADC1CH3[iADC1+iADC2] = ADC1_Ping_Buffer.ADC1CH0[iADC2];
        }
    }
    else 
    {
        for(iADC2 = 0; iADC2 < DMA_COPY_DEPTH; iADC2++) 
        {
            ADC1CH0[iADC1+iADC2] = ADC1_Pong_Buffer.ADC1CH3[iADC2];
//            ADC1CH1[NUMBER_OF_SAMPLES-1-(iADC1+iADC2)] = ADC1_Pong_Buffer.ADC1CH2[iADC2];
//            ADC1CH2[NUMBER_OF_SAMPLES-1-(iADC1+iADC2)] = ADC1_Pong_Buffer.ADC1CH1[iADC2];
//            ADC1CH3[NUMBER_OF_SAMPLES-1-(iADC1+iADC2)] = ADC1_Pong_Buffer.ADC1CH0[iADC2];
            
            ADC1CH1[iADC1+iADC2] = ADC1_Ping_Buffer.ADC1CH2[iADC2];
            ADC1CH2[iADC1+iADC2] = ADC1_Ping_Buffer.ADC1CH1[iADC2];
            ADC1CH3[iADC1+iADC2] = ADC1_Ping_Buffer.ADC1CH0[iADC2];
        }
    }
    
    ADC1CH0[0] = 0;
    ADC1CH1[0] = 0;
    ADC1CH2[0] = 0;
    ADC1CH3[0] = 0;
    
    ADC1CH0[1] = 0;
    ADC1CH1[1] = 0;
    ADC1CH2[1] = 0;
    ADC1CH3[1] = 0;
    
//    ADC1CH0[1023] = 0;
//    ADC1CH1[1023] = 0;
//    ADC1CH2[1023] = 0;
//    ADC1CH3[1023] = 0;
//    
//    ADC1CH0[1022] = 0;
//    ADC1CH1[1022] = 0;
//    ADC1CH2[1022] = 0;
//    ADC1CH3[1022] = 0;
    
    iADC1 += DMA_COPY_DEPTH;
    
    if(iADC1 == NUMBER_OF_SAMPLES) 
    {
        AD1CON1bits.ADON = 0; // Turn off ADC1
        AutoGainControlAndCentering();
        
    }
    else;

    pingPongToggle ^= 1; // Toggle between each DMA buffer
    _DMA0IF = 0; // Clear DMA0 interrupt flag
}

// No Trigger Watchdog Timer Interrupt Service Routine (Timer 4/5)
void __attribute__((interrupt, auto_psv)) _T5Interrupt(void) 
{ 
    captureFronts = 0; //Prevent ADCs from sampling and reacquire signal.
    if(PGA_gain < 7)
    {
        PGA_gain++; // Increase gain one notch 
        WritePGA(PGA_gain);
    }
    UpdateLCD_AcqSig();

    TMR5 = 0; // Reset 32-bit TMR45 counter value. Time since last "Trigger Event"
    TMR4 = 0; // Reset 32-bit TMR45 counter value. Time since last "Trigger Event"
    TMR8 = 0; // Reset 32-bit TMR3 (MSW) (period between pulse counter)
    TMR9 = 0; // Reset 32-bit TMR3 (LSW) (period between pulse counter)
    _T5IF = 0; // Clear TMR5 interrupt flag
}

void ActivateADC(void)
{
    iADC1 = 0; // Initialize ADC sample transfer index
    AD1CON1bits.ADON = 1; // Enable ADC1
    
    ADC_count = ADC_count + 1;
}

void UpdateLCD_Sampling(void)
{
    clearLCD();

    setCursor(1,2);
    putsLCD("Max");
    
    setCursor(1,8);
    putsLCD("Min");
    
    setCursor(1,8);
    putsLCD("Min");
    
    setCursor(1,14);
    putsLCD("Samp#:");
    setCursor(2,16);
    sprintf(strL1,"%d", ADC_count);
    putsLCD(strL1);

    setCursor(2,0);
    putsLCD("O:");
    setCursor(3,0);
    putsLCD("X:");
    setCursor(4,0);
    putsLCD("Y:");
    
    setCursor(2,2);
    sprintf(strL1,"%d", vecMax_o);
    putsLCD(strL1);
    setCursor(2,6);
    putsLCD(" ");
    
    setCursor(3,2);
    sprintf(strL1,"%d", vecMax_x);
    putsLCD(strL1);
    setCursor(3,6);
    putsLCD(" ");
    
    setCursor(4,2);
    sprintf(strL1,"%d", vecMax_y);
    putsLCD(strL1);
    setCursor(4,6);
    putsLCD(" ");
    
    setCursor(2,7);
    sprintf(strL1,"%d", vecMin_o);
    putsLCD(strL1);
    setCursor(3,7);
    sprintf(strL1,"%d", vecMin_x);
    putsLCD(strL1);
    setCursor(4,7);
    sprintf(strL1,"%d", vecMin_y);
    putsLCD(strL1);
    
    setCursor(2,12);
    putsLCD("    ");
    setCursor(3,12);
    putsLCD("       ");
    setCursor(4,12);
    putsLCD("  ");

    setCursor(4,14);
    putsLCD("Gain: ");
    setCursor(4,19);
    sprintf(strL1,"%d", PGA_gain);
    putsLCD(strL1);
    
//    setCursor(2,0);
//    putsLCD("pHigh: ");
//    setCursor(2,7);
//    sprintf(strL1,"%d", portionHigh);
//    putsLCD(strL1);
//
//    setCursor(3,0);
//    putsLCD("PR6: ");
//    setCursor(3,5);
//    sprintf(strL1,"%u", PR6);
//    putsLCD(strL1);
    
//    setCursor(2,0);
//    putsLCD("vMax: ");
//    setCursor(2,6);
//    sprintf(strL1,"%d", portionHigh);
//    putsLCD(strL1);
//
//    setCursor(3,0);
//    putsLCD("vMin: ");
//    setCursor(3,6);
//    sprintf(strL1,"%d", PR6);
//    putsLCD(strL1);

//    setCursor(4,0);
//    putsLCD("Gain: ");
//    setCursor(4,6);
//    sprintf(strL1,"%d", PGA_gain);
//    putsLCD(strL1);
}

void UpdateLCD_AcqSig(void)
{
    clearLCD();
    
    setCursor(1,0);
    putsLCD("Acquiring signal...");
    
    setCursor(3,0);
    putsLCD("Gain set to: ");
    setCursor(3,13);
    sprintf(strL1,"%d", PGA_gain);
    putsLCD(strL1); 
    
    if(PGA_gain == 7)
    {
        setCursor(3,16);
        putsLCD("MAX!");
        PGA_gain = 1;
    }
    else
    {
        setCursor(3,16);
        putsLCD("    ");
    }
}