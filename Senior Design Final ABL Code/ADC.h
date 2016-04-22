#include <xc.h>
#include "globaldefs.h"
#include "jps_math.h"
#include "atanLUT.h"

void InitADC1(void) 
{
    //////////////////////////////////////
    ////////// CONFIGURING ADC1 //////////
    //////////////////////////////////////
    AD1CON1 = 0x0000; // Turn off the ADC until it is configured
    AD1CON1 = 0x03EC; // Continue in idle mode
                      // DMA scatter/gather mode
                      // 10-bit signed fractional
                      // Auto-convert
                      // Simultaneous sample 4-analog input channels
                      // Sampling begins immediately after last conversion
    AD1CON2 = 0x0200; // Use AVDD and AVSS as voltage references
                      // Disable scan mode
                      // Convert SAH's: CH0,CH1,CH2, and CH3
    AD1CON3 = 0x0204; // ADC conversion clock derived from the system clock
                      // Sample for 2 TAD
                      // TAD = 75 ns
    AD1CON4 = 0x0100; // DMA enabled
                      // All results written to ADC1BUF0              
    // DMA Copy Depth (Number of samples per analog channel)

    if(DMA_COPY_DEPTH == 128) 
    { 
    	AD1CON4bits.DMABL = 7; 
    }
    else if(DMA_COPY_DEPTH == 64) 
    { 
    	AD1CON4bits.DMABL = 6; 
    }
    else if(DMA_COPY_DEPTH == 32) 
    { 
    	AD1CON4bits.DMABL = 5; 
    }
    else if(DMA_COPY_DEPTH == 16) 
    { 
    	AD1CON4bits.DMABL = 4; 
    }
    else if(DMA_COPY_DEPTH == 8)
	{ 
		AD1CON4bits.DMABL = 3; 
	}
    else if(DMA_COPY_DEPTH == 4)
    { 
    	AD1CON4bits.DMABL = 2; 
    }
    else if(DMA_COPY_DEPTH == 2)
    { 
    	AD1CON4bits.DMABL = 1; 
    }
    else
    { 
    	AD1CON4bits.DMABL = 0; 
    }

    AD1CHS0 = 0x0003; // CH0 positive input is AN3
                      // CH0 negative input is AVSS
    AD1CHS123 = 0x0000; // Sample A: CH1 positive input is AN0, CH1 negative input is AVSS
                        // Sample A: CH2 positive input is AN1, CH2 negative input is AVSS
                        // Sample A: CH3 positive input is AN2, CH3 negative input is AVSS
    AD1CSSL = 0x0000; // Skip AN0 -AN15 (for scan mode)
    AD1CSSH = 0x0000; // Skip AN16-AN31 (for scan mode)
    _AD1IF = 0; // Clear ADC1 interrupt flag
    _AD1IE = 0; // Disable ADC1 interrupt

    ////////////////////////////////////////////////////////
    ////////// CONFIGURING DMA0 FOR USE WITH ADC1 //////////
    ////////////////////////////////////////////////////////

    DMA0CON = 0x0000; // Disable DMA0 until it is configured
    DMA0CON = 0x0022; // Word size transfer
                      // Initiate interrupt after all data has been MOVed
                      // Peripheral Indirect Addressing mode
                      // Continuous Ping Pong mode 

    DMA0REQ = 0x000D; // Automatic DMA transfer initiation by DMA request
                      // Select ADC1 as DMA request source
    DMA0CNT = (4*DMA_COPY_DEPTH) - 1; // Number of DMA requests
    DMA0PAD = (volatile uint16_t)&ADC1BUF0; // Point DMA to ADC1BUF0
    // Passing addresses of ADC1 buffers
    DMA0STAL = (uint16_t)&ADC1_Ping_Buffer;
    DMA0STAH = 0x0000;
    DMA0STBL = (uint16_t)&ADC1_Pong_Buffer;
    DMA0STBH = 0x0000;
    _DMA0IP = 7;
    _DMA0IF = 0; // Clear DMA0 interrupt flag
    _DMA0IE = 1; // Enable DMA interrupt
    DMA0CONbits.CHEN = 1; // Enable DMA Channel 0
}