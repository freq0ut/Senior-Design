#include "config.h"
#include <p33Exxxx.h>
#include <dsp.h>
#include "fft.h"

/*******************************************************************
** Device configuration register macros for building the hex file **
*******************************************************************/
//_FOSC(CSW_FSCM_OFF & XT_PLL8);// XT with 8xPLL oscillator, Failsafe clock off
//_FWDT(WDT_OFF);               // Watchdog timer disabled
//_FBORPOR(PBOR_OFF & MCLR_EN); // Brown-out reset disabled, MCLR reset enabled
//_FGS(CODE_PROT_OFF);          // Code protect disabled


// Extern Definitions
extern fractcomplex sigCmpx[FFT_BLOCK_LENGTH]// Real input signal, stored in a complex array in Y-Space
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
	int i = 0;
	fractional *p_real = &sigCmpx[0].real;
	fractcomplex *p_cmpx = &sigCmpx[0];

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
	p_real = &sigCmpx[(FFT_BLOCK_LENGTH/2)-1].real;
	p_cmpx = &sigCmpx[FFT_BLOCK_LENGTH-1];

	// Convert the Real input sample array to a Complex input sample array.
	// We do this by zeroing out the imaginary part of each data sample.
	for ( i = FFT_BLOCK_LENGTH; i > 0; i-- ) {
		(*p_cmpx).real = (*p_real--);
		(*p_cmpx--).imag = 0x0000;
	}
    
	// Perform FFT operation
	#ifndef FFTTWIDCOEFFS_IN_PROGMEM
		FFTComplexIP (LOG2_BLOCK_LENGTH, &sigCmpx[0], &twiddleFactors[0], COEFFS_IN_DATA);
	#else
		FFTComplexIP (LOG2_BLOCK_LENGTH, &sigCmpx[0], (fractcomplex *) __builtin_psvoffset(&twiddleFactors[0]), (int) __builtin_psvpage(&twiddleFactors[0]));
	#endif

	// Store output samples in bit-reversed order of their addresses
	BitReverseComplex (LOG2_BLOCK_LENGTH, &sigCmpx[0]);

	// Compute the square magnitude of the complex FFT output array so we have a Real output vector
	SquareMagnitudeCplx(FFT_BLOCK_LENGTH, &sigCmpx[0], &sigCmpx[0].real);

	// Find the frequency Bin ( = index into the SigCmpx[] array) that has the largest energy
	VectorMax(FFT_BLOCK_LENGTH/2, &sigCmpx[0].real, &peakFrequencyBin);

	// Compute the frequency (in Hz) of the largest spectral component
	peakFrequency = peakFrequencyBin*(SAMPLING_RATE/FFT_BLOCK_LENGTH);

	while (1);
}
