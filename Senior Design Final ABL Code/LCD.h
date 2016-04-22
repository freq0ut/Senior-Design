/* 
 * File:   LCD.h
 * Author: Zachary
 *
 * Created on October 10, 2015, 12:31 PM
 */

#ifndef LCD_H
#define	LCD_H

#include "globaldefs.h"

#define LCD_DATA    1
#define LCD_CMD     0
#define PM_DATA     PMDIN1
#define PM_ADDR     LATGbits.LATG6

void initLCD(void);
char readLCD(int addr);
void writeLCD(int addr, char c);
void putsLCD(char *s);
void initPMP(void);
void setCursor(int row, int col);
void clearLCD(void);
void homeScreen(void);

#define busyLCD()   readLCD(LCD_CMD) & 0x80
#define addrLCD()   readLCD(LCD_CMD) & 0x7F
#define getLCD()    readLCD(LCD_DATA)
#define putLCD(d)   writeLCD(LCD_DATA,d)
#define cmdLCD(c)   writeLCD(LCD_CMD, (c))
#define homeLCD()   writeLCD(LCD_CMD, 2)
#define clrLCD()    writeLCD(LCD_CMD, 1)
#define setLCDG(a)  writeLCD(LCD_CMD, (a & 0x3F) | 0x40)
#define setLCDC(a)  writeLCD(LCD_CMD, (a & 0x7F) | 0x80)

void initPMP(void)
{
    PMCON = 0x0000;
    PMCON = 0x033F;     // Enable the PMP, long waits
    PMMODE = 0x03FF;    // Master Mode 1
    PMAEN = 0x0000; 
    PMCON |= 0x8000;
}

void initLCD(void)
{
    while(PMMODEbits.BUSY);     // wait for PMP to be available
    delay_ms(50);
    PM_ADDR = LCD_CMD;
    PM_DATA = 0b00111000;       // 8-bit, 2 lines, 5x7
    while(busyLCD()){};
    delay_us(50);
    PM_DATA = 0b00001000;       // display off
    while(busyLCD()){};
    delay_us(50);
    PM_DATA = 0b00000001;       // clear display
    while(busyLCD()){};
    delay_ms(1.8);
    PM_DATA = 0b00000110;       // increment cursor, no shift
    while(busyLCD()){};
    delay_us(50);
    PM_DATA = 0b00000010;       // home command
    while(busyLCD()){};
    delay_us(50);
    PM_DATA = 0b00001100;       // display on
    while(busyLCD()){};
    delay_us(50);
    PMMODEbits.BUSY = 0;
}

void putsLCD(char *s)
{
    while(*s)
    {
        putLCD(*s++);
    }  
}

char readLCD(int addr)
{
    int i;
    while(PMMODEbits.BUSY);     // wait for PMP to be available
    PM_ADDR = addr;              // select the command address
    i = PM_DATA;            // initiate a read cycle, dummy
    while(PMMODEbits.BUSY);     // wait for PMP to be available
    return(PM_DATA);            // read the status register
}

void writeLCD(int addr, char c)
{
    while(PMMODEbits.BUSY);      // wait for PMP to be available
    PM_ADDR = addr;              // select the command address
    PM_DATA = c;
    while(busyLCD());
    delay_us(50);
    PMMODEbits.BUSY = 0;
}

void clearLCD(void)
{
    while(PMMODEbits.BUSY);
    PM_ADDR = LCD_CMD;
    PM_DATA = 0b00000001; // clear display
    delay_ms(1.8);
    PMMODEbits.BUSY = 0;
}

void setCursor(int row, int col)
{
    while(PMMODEbits.BUSY);
    if(row == 4)
    {
        setLCDC((0x54 + col));
    }   
    else if(row == 3)
    {
        setLCDC((0x14 + col));
    }
    else if(row == 2)
    {
        setLCDC((0x40 + col));
    }  
    else
    {
        setLCDC((0x00 + col));
    }      
    PMMODEbits.BUSY = 0;    
}

void homeScreen(void)
{
    setCursor(1,0);
    putsLCD("Acoustic Beacon Loc.");
    setCursor(2,0);
    putsLCD("   Design Team F");
    setCursor(3,0);
    putsLCD("    U. of Akron");
    setCursor(4,0);
    putsLCD("    Dr. Toonen");
    delay_ms(1000);
}

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef	__cplusplus
}
#endif

#endif	/* LCD_H */
