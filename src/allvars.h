/*****************************************************************************************
Title: allvars.h
Purpose: Contains global variable and struct definitions
Notes:
*****************************************************************************************/
#ifndef ALLVARS_H
    #define ALLVARS_H

    #include <gsl/gsl_rng.h>
    #include <fftw3.h>

    // MPI parameters
    extern int ntasks;   // Total number of processors
    extern int thistask; // Current processor id

    // Parameter file...parameters
    extern char file_base[256]; // Input spectrum base name (everything but the id number)
    extern int nspec_files;     // Total number of input spectra (gives rng range)
    extern int nfiles_to_use;   // Total number of input spectra to use per ps
    extern int npspec;          // Total number of ps realizations to create
    extern int use_sdss;        // Flag for using SDSS data. Results in using sdss k vals
    extern int use_spec_k;      // Flag for using k found from spectra as binning k
    extern int use_mcdonald;    // Flag for using min k found in McDonald+ 2000 for bins
    extern int use_all;         // If true, use every Ly alpha spectrum to get a ps
    extern int use_noise;       // Noise specs are laid out differently. Sets right format
    extern float redshift;      // z of snapshot. Used for getting the mean flux

    // General globals
    extern int max_iters;       // Max number of iters to use in while loops
    extern int exit_code_local; // Set to 0 when there's a problem
    extern int exit_code_global;// Product of each processor's local exit code
    extern int npspec_local;    // # of ps realizations for each task (= npspec/ntasks)
    extern int snapnum;         // The id # of snapshot. Used for labeling output files
    extern gsl_rng *rng;        // Random number generator. Used for choosing input files

    // Power spectrum globals
    extern float mean_flux;     // Mean flux as determined by observations
    extern float *bin_freqs;    // The x-axis for the final power spectrum
    extern float *bin_ps;       // The normalized, final power spectrum
    extern float upper_k_limit; // The values in bin_freqs are lower lims. This is the
                                // upper limit for the last bin, used in binning the ps
    extern int *counts;         // Number of counts in each bin. Used to norm ps
    extern int ncounts;         // Size of bin_freqs, bin_ps, and counts arrays (# bins)

    // Debugging params
    #ifdef DEBUGGING
        extern int *py_file_ids;
    #endif

    // Power Spectrum struct
    typedef struct POWER_SPECTRUM
    {
        int npix;
        float *vel;
        float *flux;
        float *ps_amp;
    } POWER_SPECTRUM;
#endif
