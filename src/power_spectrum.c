/*****************************************************************************************
Title: power_spectrum.c
Purpose: Contains functions related to calculating the power spectrum
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "allvars.h"
#include "proto.h"



/********************************************
               power_spectrum
********************************************/
void power_spectrum(POWER_SPECTRUM *ps)
{
    // This function calculates the power spectrum of the quantity F_p = F/mean_flux - 1.0
    // because this quantity is not as sensitive to the mean flux, which is uncertain.

    float *Fp;
    fftwf_complex *Fp_k;
    int i;
    fftwf_plan plan;

    // Allocate memory for signal. This is the array that holds F_p defined above
    if(!(Fp = calloc(ps->npix, sizeof(float))))
    {
        printf("Error, could not allocate memory for Fp on task %d\n", thistask);
        exit_code_local = 0;
        return;
    }

    // Allocate memory for Fp_k. From the fftw3 docs: For a real transform of rank d, 
    // the complex data is an n0 x n1 x ... x (n_(d-1)/2 + 1) array of fftw_complex 
    // values in row-major order
    if(!(Fp_k = (fftwf_complex *)fftwf_malloc(((ps->npix/2) + 1) * \
        sizeof(fftwf_complex))))
    {
        printf("Error, could not allocate memory for Fp_k on task %d\n", thistask);
        exit_code_local = 0;
        free(Fp);
        return;
    }

    // Allocate memory for ps_amp
    if(!(ps->ps_amp = calloc((ps->npix / 2) + 1, sizeof(float))))
    {
        printf("Error, could not allocate memory for ps_amp!\n");
        exit_code_local = 0;
        free(Fp);
        fftwf_free(Fp_k);
        return;
    }

    // Calculate Fp
    for(i = 0; i < ps->npix; i++)
    {
        Fp[i] = (ps->flux[i] / mean_flux) - 1.0;
    }

    // Now take the fft of the data. The data is all real and one dimensional
    plan = fftwf_plan_dft_r2c_1d(ps->npix, (float *)Fp, (fftwf_complex *)Fp_k, \
        FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();

    // Now use the complex data to get the power spectrum. This is because
    // P = A(Fp_k \cdot Fp_k*)^2. Check whether I need to divide by N here!!
    for(i = 0; i < (ps->npix / 2) + 1; i++)
    {
        ps->ps_amp[i] = pow(cabsf(Fp_k[i]), 2.0) / (ps->npix * ps->npix);
    }

    // Clean up
    free(Fp);
    fftwf_free(Fp_k);
}



/********************************************
           bin_power_spectrum
********************************************/
void bin_power_spectrum(POWER_SPECTRUM ps)
{
    // This function calculates the k values for the power spectrum and then adds
    // the corresponding P values to the appropriate bin

    int i;
    int j;
    float T;
    float *freqs;

    // In order to get the k values, we need the period of the data, T. For a DFT, the
    // data is assumed to be periodic over the range in which there is data. This is
    // just the difference between the last and first velocity element.
    T = fabs(ps.vel[ps.npix -1] - ps.vel[0]);

    // Allocate memory for frequencies (k values)
    if(!(freqs = calloc((ps.npix / 2) + 1, sizeof(float))))
    {
        printf("Error, could not allocate memory for freqs on task %d!\n", thistask);
        exit_code_local = 0;
        return;
    }

    // Now get the k values
    for(i = 0; i < (ps.npix / 2) + 1; i++)
    {
        freqs[i] = (2.0 * M_PI * i) / T;

        // Set the bin k here, if using that option
        if(use_spec_k == 1)
        {
            bin_freqs[i] = freqs[i];
        }

        // We also multiply the ps by T.
        // NOTE: I do not understand this. S. Bertone did this in her code, but I don't
        // know why it's done. The only thing I can think of is for units. The PS needs
        // to be normalized and is dimensionless (because the argument we took the FFT of
        // was dimensionless, I think this makes the PS dimensionless) so multiplying by T
        // puts it in the right units
        ps.ps_amp[i] *= T;

        // Now bin the power spectrum
        // This is wrong. It skips the last bin. Also, the bins aren't evenly spaced,
        // which makes getting the width annoying, as each width will be different.

        // This loop does every bin EXCEPT the last one
        for(j = 0; j < ncounts - 1; j++)
        {
            if((freqs[i] >= bin_freqs[j]) && (freqs[i] < bin_freqs[j + 1]))
            {
                // Update the counts for i'th bin and add value of PS to total
                bin_ps[j] += ps.ps_amp[i];
                counts[j]++;
        
                // It can only be in one bin, so there's not point in finishing the loop
                break;
            }
        }

        // Now do the last bin
        if((freqs[i] >= bin_freqs[ncounts - 1]) && (freqs[i] < upper_k_limit))
        {
            bin_ps[ncounts - 1] += ps.ps_amp[i];
            counts[ncounts - 1]++;
        }
    }

    // Clean up
    free(freqs);
}



/********************************************
        normalize_power_spectrum
********************************************/
void normalize_power_spectrum(void)
{
    // This function simply divides out the number of counts in each bin so that the
    // final power spectrum is independent of the counts

    int i;

    for(i = 0; i < ncounts; i++)
    {
        if(counts[i] > 0)
        {
            bin_ps[i] /= counts[i];
        }
    }
}
