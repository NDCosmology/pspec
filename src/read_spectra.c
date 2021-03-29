/*****************************************************************************************
Title: read_spectra.c
Purpose: Contains functions related to reading in the input spectra files
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"



/********************************************
             read_spectrum_file
********************************************/
void read_spectrum_file(char *fbase, POWER_SPECTRUM *ps)
{
    // This function randomly chooses a spectrum file to read

    FILE *fd = NULL;
    int d = 0;
    int nlines = 0;
    int i = 0;
    float dummy = 0.0;
    char f_in[256];

    #ifdef DEBUGGING
        // Choose from the next entry in the py_file_ids list
        static int j = 0;
        d = py_file_ids[j];
        j++;
    #endif

    #ifndef DEBUGGING
        // Choose a file randomly
        d = gsl_rng_uniform_int(rng, nspec_files);
    #endif

    // Build file name
    if(nspec_files < 100)
    {
        sprintf(f_in, "%s%02d.txt", fbase, d); 
    }

    else if(nspec_files < 1000)
    {
        sprintf(f_in, "%s%03d.txt", fbase, d); 
    }

    else if(nspec_files < 10000)
    {
        sprintf(f_in, "%s%04d.txt", fbase, d); 
    }

    // Open file for reading
    if(!(fd = fopen(f_in, "r")))
    {
        printf("Error, could not open file for reading input spectrum %d on task %d!\n",\
            d, thistask);
        exit_code_local = 0;
        return;
    }

    #ifdef DEBUGGING
        FILE *flog;

        // Write the file name to a log file
        if(!(flog = fopen("log.txt", "a")))
        {
            printf("Error, could not open file for writing file name!\n");
            exit(EXIT_FAILURE);
        }

        fprintf(flog, "%s\n", f_in);
        fclose(flog);
    #endif


    // Read file. There are five lines in the header and they all start with '#'. The
    // first line contains the number of pixels, which is how many lines of data there
    // are (and so gives the size of the arrays). Then, in the actual data, there are
    // 8 columns: vel, lambda, nHI, T, vHI, tau, rho_tot, and F. 
    
    // Get number of pixels
    fscanf(fd, "# Npixels: %d\n", &ps->npix);

    // Skip the rest of the header
    while(nlines < 4)
    {
        while(getc(fd) != '\n')
        {
            continue;
        }

        nlines++;
    }

    // Now allocate memory for tau and vel
    if(!(ps->flux = calloc(ps->npix, sizeof(float))))
    {
        printf("Error, could not allocate memory for F on task %d, file %d\n", \
            thistask, d);
        exit_code_local = 0;
        return;
    }

    if(!(ps->vel = calloc(ps->npix, sizeof(float))))
    {
        printf("Error, could not allocate memory for vel on task %d, file %d\n", \
            thistask, d);
        exit_code_local = 0;
        return;
    }

    // Read the rest of the data
    if(use_noise == 0)
    {
        for(i = 0; i < ps->npix; i++)
        {
            fscanf(fd, "%f %f %f %f %f %f %f %f\n", &ps->vel[i], &dummy, &dummy, &dummy, \
                &dummy, &dummy, &dummy, &ps->flux[i]);
        }
    }

    else
    {
        for(i = 0; i < ps->npix; i++)
        {
            fscanf(fd, "%f %f\n", &ps->vel[i], &ps->flux[i]);
        }
    }

    fclose(fd);
}



/********************************************
           read_all_spectrum_files
********************************************/
void read_all_spectrum_files(char *fbase, POWER_SPECTRUM *ps, int file_num)
{
    // This function reads in the file given by file_num

    FILE *fd = NULL;
    int nlines = 0;
    int i = 0;
    float dummy = 0.0;
    char f_in[256];

    // Build file name
    if(nspec_files < 100)
    {
        sprintf(f_in, "%s%02d.txt", fbase, file_num); 
    }

    else if(nspec_files < 1000)
    {
        sprintf(f_in, "%s%03d.txt", fbase, file_num); 
    }

    else if(nspec_files < 10000)
    {
        sprintf(f_in, "%s%04d.txt", fbase, file_num); 
    }

    // Open file for reading
    if(!(fd = fopen(f_in, "r")))
    {
        printf("Error, could not open file for reading input spectrum %d on task %d!\n",\
            file_num, thistask);
        exit_code_local = 0;
        return;
    }

    // Read file. There are five lines in the header and they all start with '#'. The
    // first line contains the number of pixels, which is how many lines of data there
    // are (and so gives the size of the arrays). Then, in the actual data, there are
    // 8 columns: vel, lambda, nHI, T, vHI, tau, rho_tot, and F. 
    
    // Get number of pixels
    fscanf(fd, "# Npixels: %d\n", &ps->npix);

    // Skip the rest of the header
    while(nlines < 4)
    {
        while(getc(fd) != '\n')
        {
            continue;
        }

        nlines++;
    }

    // Now allocate memory for tau and vel
    if(!(ps->flux = calloc(ps->npix, sizeof(float))))
    {
        printf("Error, could not allocate memory for F on task %d, file %d\n", \
            thistask, file_num);
        exit_code_local = 0;
        return;
    }

    if(!(ps->vel = calloc(ps->npix, sizeof(float))))
    {
        printf("Error, could not allocate memory for vel on task %d, file %d\n", \
            thistask, file_num);
        exit_code_local = 0;
        return;
    }

    // Read the rest of the data
    if(use_noise == 0)
    {
        for(i = 0; i < ps->npix; i++)
        {
            fscanf(fd, "%f %f %f %f %f %f %f %f\n", &ps->vel[i], &dummy, &dummy, &dummy, \
                &dummy, &dummy, &dummy, &ps->flux[i]);
        }
    }

    else
    {
        for(i = 0; i < ps->npix; i++)
        {
            fscanf(fd, "%f %f\n", &ps->vel[i], &ps->flux[i]);
        }
    }

    fclose(fd);
}
