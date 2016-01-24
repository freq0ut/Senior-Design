/*
*   Joshua Simmons
*
*   January 21, 2016
*
*   Turns an LED on and off using a delay.
*/

#include "config.h"
#include <dsp.h>
#include "fft.h"

// Extern Definitions
extern fractcomplex sigCmpx1[FFT_BLOCK_LENGTH]// Real input signal, stored in a complex array in Y-Space
__attribute__ ((section (".ydata, data, ymemory"),
aligned (FFT_BLOCK_LENGTH*2*2)));

extern fractcomplex sigCmpx2[FFT_BLOCK_LENGTH]// Real input signal, stored in a complex array in Y-Space
__attribute__ ((section (".ydata, data, ymemory"),
aligned (FFT_BLOCK_LENGTH*2*2)));

// Global Definitions
#ifndef FFTTWIDCOEFFS_IN_PROGMEM
	fractcomplex twiddleFactors[FFT_BLOCK_LENGTH/2]// Declare Twiddle Factor array in X-Space
	__attribute__ ((section (".xbss, bss, xmemory"), aligned (FFT_BLOCK_LENGTH*2)));
#else
	extern const fractcomplex twiddleFactors[FFT_BLOCK_LENGTH/2]// Twiddle Factor array in Program memory
	__attribute__ ((space(auto_psv), aligned (FFT_BLOCK_LENGTH*2)));
#endif

int	peakFrequencyBin = 0;
unsigned long peakFrequency = 0;

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
    
    /************************************************
    ********** COMPUTING CROSS-CORRELATION **********
    ************************************************/
    
	int i = 0;
	fractional *p_real = &sigCmpx1[0].real;
	fractcomplex *p_cmpx = &sigCmpx1[0];

	// Generate TwiddleFactor Coefficients (Only once at stary-up)
	#ifndef FFTTWIDCOEFFS_IN_PROGMEM
		TwidFactorInit (LOG2_BLOCK_LENGTH, &twiddleFactors[0], 0);
	#endif

	// Scaling the input data to be in the range of [-0.5, +0.5] by right bit shifting once.
	// For optimization: 
	//		1. Scale when first obtaining the time samples. 
	//		2. Or, within the BitReverseComplex function source code.
	for ( i=0; i < FFT_BLOCK_LENGTH; i++ ) {
		*p_real = *p_real >>1;
		*p_real++;
	}

	// Set up pointers to convert real array to a complex array.
	// The input array initially has all the real input samples followed by a series of zeros. 
	p_real = &sigCmpx1[(FFT_BLOCK_LENGTH/2)-1].real;
	p_cmpx = &sigCmpx1[FFT_BLOCK_LENGTH-1];

	// Convert the Real input sample array to a Complex input sample array.
	// We do this by zeroing out the imaginary part of each data sample.
	for ( i = FFT_BLOCK_LENGTH; i > 0; i-- ) {
		(*p_cmpx).real = (*p_real--);
		(*p_cmpx--).imag = 0x0000;
	}
    
	// Perform FFT operation
	#ifndef FFTTWIDCOEFFS_IN_PROGMEM
		FFTComplexIP (LOG2_BLOCK_LENGTH, &sigCmpx1[0], &twiddleFactors[0], COEFFS_IN_DATA);
	#else
		FFTComplexIP (LOG2_BLOCK_LENGTH, &sigCmpx1[0], (fractcomplex *) __builtin_psvoffset(&twiddleFactors[0]), (int) __builtin_psvpage(&twiddleFactors[0]));
	#endif

	// Store output samples in bit-reversed order of their addresses
	BitReverseComplex (LOG2_BLOCK_LENGTH, &sigCmpx1[0]);

	// Compute the square magnitude of the complex FFT output array so we have a Real output vector
	SquareMagnitudeCplx(FFT_BLOCK_LENGTH, &sigCmpx1[0], &sigCmpx1[0].real);

	// Find the frequency Bin ( = index into the SigCmpx[] array) that has the largest energy
	VectorMax(FFT_BLOCK_LENGTH/2, &sigCmpx1[0].real, &peakFrequencyBin);

	// Compute the frequency (in Hz) of the largest spectral component
	peakFrequency = peakFrequencyBin*(SAMPLING_RATE/FFT_BLOCK_LENGTH);

	// Infinite Loop
	while (1);
}
