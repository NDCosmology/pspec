/*****************************************************************************************
Title: initialize.c
Purpose: Contains functions for setting up the run and reading the parameter file
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/********************************************
               initialize_run
********************************************/
void initialize_run(int nargs, char **arglist, POWER_SPECTRUM *ps)
{
    // This function does an args check, reads in the parameter file, and sets up
    // some of the global variables needed for the run

    FILE *fd;
    char buffer[256];
    char buffer2[256];
    char *error;
    int elements;
    int iters = 0;
    int i;
    size_t seed;
    float tau_avg;

    // Initialize the arrays in the power spectrum struct (this is because, if
    // uninitialized, if there is an args error, when free_globals gets called in
    // check_exit_code then there is an attempt to free these pointers when they haven't
    // been allocated, which results in a seg fault.
    ps->vel = NULL;
    ps->flux = NULL;
    ps->ps_amp = NULL;
     

    // Do an args check
    if(nargs != 2)
    {
        printf("Usage: ./pspec ./param_file.txt\n");
        exit_code_local = 0;
    }

    check_exit_code(exit_code_local, ps);

    if(thistask == 0)
    {
        // Read in the parameter file
        if(!(fd = fopen(arglist[1], "r")))
        {
            printf("Error, could not open parameter file for reading!\n");
            exit_code_local = 0;
        }

        // The file contains spectrum file base, total number of spectrum files, number of
        // spectrum files to use per power spectrum, redshift, and the number of power 
        // spectra to make

        // File base
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "Input Spectra File Base:%s", file_base);
        exit_code_local = checkforerror(error, elements, &buffer[0]);

        // total number of spectrum files
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "Total Number of Input Spectra:%s", buffer2);
        exit_code_local = checkforerror(error, elements, &buffer[0]);
        nspec_files = atoi(buffer2);

        // number of spectrum files to use per power spectrum
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "Number of Input Spectra To Use:%s", buffer2);
        exit_code_local = checkforerror(error, elements, &buffer[0]);
        nfiles_to_use = atoi(buffer2);

        // redshift
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "Redshift:%s", buffer2);
        exit_code_local = checkforerror(error, elements, &buffer[0]);
        redshift = atof(buffer2);

        // number of power spectra to make (total, across all processors)
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "Total Number of Power Spectra:%s", buffer2);
        exit_code_local = checkforerror(error, elements, &buffer[0]);
        npspec = atoi(buffer2);

        // Snapshot number (because it's easier than parsing it from basename)
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "Snapshot Number:%s", buffer2);
        exit_code_local = checkforerror(error, elements, &buffer[0]);
        snapnum = atoi(buffer2);

        // Use sdss data set (this results in different k values being used)
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "use_sdss:%s", buffer2);
        exit_code_local = checkforerror(error, elements, &buffer[0]);
        use_sdss = atoi(buffer2);

        // Use spectral freqs (this results in the array of k to bin to being the same
        // as the k found in power_spectrum.c, thereby eliminating the empty bins issue
        // without having to have spectra of varying length. The idea is that this does
        // not affect small scale power at all, but since the spectra will be longer,
        // there will be more lines separated by larger scales, so it should hopefully
        // increase the large-scale power)
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "use_spec_k:%s", buffer2);
        exit_code_local = checkforerror(error, elements, &buffer[0]);
        use_spec_k = atoi(buffer2);

        // Use McDonald. Uses The k_min from McDonald+ 2000 for the binning
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "use_mcdonald:%s", buffer2);
        exit_code_local = checkforerror(error, elements, &buffer[0]);
        use_mcdonald = atoi(buffer2);

        // Use all Ly alpha spectra to get a power spectrum
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "use_all:%s", buffer2);
        exit_code_local = checkforerror(error, elements, &buffer[0]);
        use_all = atoi(buffer2);

        // Use the noisy spectra instead of the perfect ones
        error = fgets(buffer, sizeof(buffer), fd);
        elements = sscanf(buffer, "use_noise:%s", buffer2);
        exit_code_local = checkforerror(error, elements, &buffer[0]);
        use_noise = atoi(buffer2);

        // Close file
        fclose(fd);
    }
    check_exit_code(exit_code_local, ps);

    // Send number of power spectra to make, number of spectra to use for each
    // power spectrum, the spectrum file base, the total number of spectrum files
    // that exist, and the redshift to the other processors. And snapsnum
    MPI_Bcast(&npspec, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nfiles_to_use, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nspec_files, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&redshift, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(file_base, 245, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&snapnum, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&use_sdss, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&use_spec_k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&use_mcdonald, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&use_all, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&use_noise, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Now figure out how many power spectra each processor needs to make. In the interest
    // of simplicity, I'm just going to increase nspec until nspec % ntasks == 0 so that
    // every processor has an equal number to make and I don't have to deal with assigning
    // stragglers.
    if(ntasks > 1)
    {
        while((npspec % ntasks) != 0)
        {
            npspec++;
            iters++;

            if(iters > max_iters)
            {
                printf("Error, max iters exceeded when finding npspec_local!\n");
                exit_code_local = 0;
            }
            check_exit_code(exit_code_local, ps);
        }
        npspec_local = npspec / ntasks;
    }

    // Handle serial case
    else if((ntasks == 1) || (ntasks == 0))
    {
        npspec_local = npspec;
    }

    // Adjust nfiles_to_use if use_all is true
    if(use_all == 1)
    {
        nfiles_to_use = nspec_files;
    }

    // Set up the random number generator on each processor. Use the current time +
    // task id for the seed
    seed = time(0) + thistask;
    //seed = 42 + thistask;
    #ifdef DEBUGGING
        seed = 42;
    #endif
    rng = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(rng, seed);

    // Check number of input spectra. I only support up to 10000
    if(nspec_files >= 10000)
    {
        printf("Error, too many input spectra!\n");
        exit_code_local = 0;
    }
    check_exit_code(exit_code_local, ps);

    // Set up the mean observed flux at the given redshift. This is used when getting
    // the power spectrum. The formula for tau_avg is from Kim 2002 sec. 3.3
    tau_avg = 0.0032 * pow(1.0 + redshift, 3.37);
    mean_flux = exp(tau_avg * -1.0);

    // Set ncounts (it's already set to the default of 30, so we only need to change it
    // if using sdss data)
    if(use_sdss == 1)
    {
        ncounts = 12;
    }

    // Set ncounts if using use_spec_k. I'm going to hard code it. Even though I know
    // that's bad, I just want to finish, at this point. The number of pixels is 3000,
    // and the number of k is (N / 2) + 1 because of the real fft symmetry. 
    if(use_spec_k == 1)
    {
        ncounts = 1501;
    }

    // Number of k in McDonald+ 2000
    if(use_mcdonald == 1)
    {
        ncounts = 18;
    }

    // Allocate memory for bin_freqs, bin_counts, and counts
    if(!(bin_freqs = calloc(ncounts, sizeof(float))))
    {
        printf("Error, could not allocate memory for bin_freqs on task %d!\n", thistask);
        exit_code_local = 0;
    }

    if(!(bin_ps = calloc(ncounts, sizeof(float))))
    {
        printf("Error, could not allocate memory for bin_ps on task %d!\n", thistask);
        exit_code_local = 0;
    }

    if(!(counts = calloc(ncounts, sizeof(float))))
    {
        printf("Error, could not allocate memory for counts on task %d!\n", thistask);
        exit_code_local = 0;
    }

    check_exit_code(exit_code_local, ps);

    // Initialize the bin_freqs array (same as Kim 2004 observations)
    // These are the LOWER LIMITS for each bin
    if(use_sdss == 0)
    {
        for(i = 0; i < ncounts; i++)
        {
            bin_freqs[i] = 0.001 * pow(10.0, (float)i/8.0);
        }
    }

    if(use_sdss == 1)
    {
        // The parameters used here were determined using scipy.optimize.curve_fit in
        // Python with the call: curve_fit(relation, x, k). x is simply the ints [0,11],
        // and k are the k values given in sdss.dat (see McDonald04). relation was a
        // function defined as: def relation(x, A, B): val = A * np.exp(B * x). The
        // params below are the best fit. When rounded to five decimal places, they fit
        // the sdss k values exactly.
        for(i = 0; i < ncounts; i++)
        {
            bin_freqs[i] = 0.00141289 * exp(0.23023187 * i);

            // Now round to five decimal places to match sdss
            bin_freqs[i] = roundf(bin_freqs[i] * 1e5) / 1e5;
        }
    }

    if(use_mcdonald == 1)
    {
        // The parameters used here were determined using scipy.optimize.curve_fit in
        // Python with the call: curve_fit(relation, x, k). x is simply the ints [0,11],
        // and k are the k values given in sdss.dat (see McDonald04). relation was a
        // function defined as: def relation(x, A, B): val = A * np.exp(B * x). The
        // params below are the best fit.
        for(i = 0; i < ncounts; i++)
        {
            bin_freqs[i] = 0.0025096 * exp(0.23034602 * i);
        }

    }

    // Get value of upper limit for last bin, used only for binning purposes.
    // This is just the next value in the above sequence (i.e. if the loop were one
    // element longer)
    if(use_sdss == 0)
    {
        upper_k_limit = 0.001 * pow(10.0, (float)ncounts/8.0);
    }

    if(use_sdss == 1)
    {
        upper_k_limit = 0.00141289 * exp(0.23023187 * ncounts);
    }

    if(use_mcdonald == 1)
    {
        // This value comes from the caption to table 4
        upper_k_limit = 0.159;
    }

    #ifdef DEBUGGING
        // I want to compare this code with the python version, so I have the python
        // version generate a file that holds the names of all of the spectra files
        // that it used to generate its power spectrum. So with this option I want to
        // read in that file and use those same spectra files.

        // Read the list of spectrum files
        if(!(fd = fopen("./list_of_used_specs.txt", "r")))
        {
            printf("Error, could not open file for reading python spec list!\n");
            exit(EXIT_FAILURE);
        }

        // Allocate memory for file ids
        if(!(py_file_ids = calloc(nfiles_to_use, sizeof(int))))
        {
            printf("Error, could not allocate memory for py_file_ids!\n");
            exit(EXIT_FAILURE);
        }

        // Read the file ids
        for(i = 0; i < nfiles_to_use; i++)
        {
            fscanf(fd, "%d", &py_file_ids[i]);
        }
        fclose(fd);
    #endif
}



/********************************************
               checkforerror
********************************************/
int checkforerror(char *error, int elements, char *buf)
{
    // This function makes sure the proper number of elements,
    // as well as the proper element, was read from the param file

    char buffer2[256];

    if(error == NULL)
    {
        fprintf(stderr, "Error, could not read from parameter file on task %d.\n", \
            thistask);
        exit_code_local = 0;
    }

    if(elements == 0)
    {
        strncpy(buffer2, buf, strlen(buf));
  
        // Add null terminating
        buffer2[strlen(buf)-1]='\0';
        fprintf(stderr, "Could not read '%s' on task %d.\n", buffer2, thistask);
        exit_code_local = 0;
    }

    return exit_code_local;
}
