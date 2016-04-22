#include <xc.h>
#include "globaldefs.h"
#include "jps_math.h"
#include "atanLUT.h"

void AutoGainControlAndCentering(void) 
{
    portionHigh = 0;
    // Finding the maximum and minimum values of ADC1CH0 Ping Buffer
    vecMax_o = VectorMax(NUMBER_OF_SAMPLES, &ADC1CH0[0], &iMax);
    vecMin_o = VectorMin(NUMBER_OF_SAMPLES, &ADC1CH0[0], &iMin);
    
    vecMax_x = VectorMax(NUMBER_OF_SAMPLES, &ADC1CH1[0], &iMax);
    vecMin_x = VectorMin(NUMBER_OF_SAMPLES, &ADC1CH1[0], &iMin);
    
    vecMax_y = VectorMax(NUMBER_OF_SAMPLES, &ADC1CH2[0], &iMax);
    vecMin_y = VectorMin(NUMBER_OF_SAMPLES, &ADC1CH2[0], &iMin);
    
    UpdateLCD_Sampling();
    //Make sure signal is still present... if not, don't clear watch dog.
    if(vecMax_o >= 1000)
    {
        TMR5 = 0; // Clear watchdog for loss of signal.
        TMR4 = 0; // Clear watchdog for loss of signal.
    }
    
    // If signal is present and too strong, decrease gain.
    if(vecMax_o >= 27000)
    {
        if(PGA_gain > 1)
        {
            PGA_gain = PGA_gain - 1; // Decrease gain one notch
            WritePGA(PGA_gain);
        }
    }
    TDOA();
}

void WritePGA(int16_t gain) 
{
    // Linear Technology LT69122
    int16_t i;
    int16_t PGA_DataByte;
    int16_t PGA_DataBit;
    int16_t PGA_SPI_QuarterClockPeriodMicro = 3; // Originally 3
    // The PGA gain must be 1,2,3,4,5,6, or 7 and nothing else.
    if(gain < 1) 
    { 
    	gain = 1; 
    }
    else if(gain > 7)
    { 
    	gain = 7; 
    }
    else{};
    PGA_DataByte = (gain << 4) + gain;
    SPI_CLK_LOW;
    delay_us(PGA_SPI_QuarterClockPeriodMicro); // Delay for a quarter SPI clock period
    SPI_CS1_LOW;
    SPI_CS2_LOW;
    for (i = 0; i < 8; i++) 
    {
        PGA_DataBit = (PGA_DataByte >> (7-i)) & 0x0001; // MSB First
        if(PGA_DataBit == 0x0001) 
        { 
        	SPI_DIN_HIGH; 
        }
        else
        { 
        	SPI_DIN_LOW;  
        }
        delay_us(PGA_SPI_QuarterClockPeriodMicro); // Delay for a quarter SPI clock period
        SPI_CLK_HIGH;
        delay_us(2*PGA_SPI_QuarterClockPeriodMicro); // Delay for a half SPI clock period
        SPI_CLK_LOW;
        delay_us(PGA_SPI_QuarterClockPeriodMicro); // Delay for a quarter SPI clock period
    }
    SPI_DIN_LOW;
    SPI_CS1_HIGH;
    SPI_CS2_HIGH;
    delay_us(PGA_SPI_QuarterClockPeriodMicro); // Delay for a quarter SPI clock period
    SPI_CLK_HIGH;
    delay_us(2*PGA_SPI_QuarterClockPeriodMicro); // Delay for a half SPI clock period
    SPI_CLK_LOW;
}