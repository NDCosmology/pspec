/*****************************************************************************************
Title:   pspec
Author:  Jared Coughlin
Date:    2/27/18
Purpose: c version of powspec_py. This uses the spectra files created by expsec to
        calcualte the flux auto power spectrum.
Notes:   1.) The reason this is written in c is so that I can parallelize it. mpi4py doesn't
            seem to be providing any speed up, for some reason.
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <fftw3.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"



int main(int argc, char **argv)
{
    // Set up MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &thistask);

    int i;
    int j;
    double start;
    double end;
    double tot_time_local;
    double tot_time_max = 0;
    POWER_SPECTRUM ps;

    // Get start time
    start = clock();

    // Initialize run
    if(thistask == 0)
    {
        printf("Initializing run...\n");
    }
    initialize_run(argc, argv, &ps);

    // Loop over the number of power spectra each processor needs to make
    if(thistask == 0)
    {
        printf("Generating power spectra...\n");
    }
    for(i = 0; i < npspec_local; i++)
    {
        // Initalize the counts and bin_ps arrays
        for(j = 0; j < ncounts; j++)
        {
            counts[j] = 0;
            bin_ps[j] = 0.0;
        }

        // Loop over the number of files required to generate a power spectrum
        for(j = 0; j < nfiles_to_use; j++)
        {
            // Read data from file
            if(use_all == 0)
            {
                read_spectrum_file(file_base, &ps);
            }

            else
            {
                read_all_spectrum_files(file_base, &ps, j);
            }

            // Check for errors that may have occurred while reading spectrum file
            check_exit_code(exit_code_local, &ps);

            // Calculate power spectrum for the current spectrum
            power_spectrum(&ps);

            // Check for errors that may have occurred while calculating P
            check_exit_code(exit_code_local, &ps);

            // Bin the power spectrum
            bin_power_spectrum(ps);

            // Check for errors that may have occurred while binning
            check_exit_code(exit_code_local, &ps);

            // Reset ps for next file
            free(ps.vel);
            free(ps.flux);
            free(ps.ps_amp);

            ps.vel = NULL;
            ps.flux = NULL;
            ps.ps_amp = NULL;
        }

        // Normalize the power spectrum
        normalize_power_spectrum();

        // Write power spectrum to a file
        write_power_spectrum(i);

        // Check for errors that may have occurred while writing data
        check_exit_code(exit_code_local, &ps);
    }
    //fclose(debug_file);

    // Clean up
    free(counts);
    free(bin_ps);
    free(bin_freqs);

    #ifdef DEBUGGING
        free(py_file_ids);
    #endif

    gsl_rng_free(rng);

    // Get end time
    end = clock();

    // Get total time
    tot_time_local = (end - start) / (double)CLOCKS_PER_SEC;

    // Get total run time
    MPI_Reduce(&tot_time_local, &tot_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(thistask == 0)
    {
        printf("Total run time: %lf sec\n", tot_time_max);
        printf("Have a nice day!\n");
    }

    MPI_Finalize();

    return 0;
}
