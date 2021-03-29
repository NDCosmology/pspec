/*****************************************************************************************
Title: write.c
Purpose: Contains functions related to writing the normalized and binned power spectrum
            to a file.
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"



/********************************************
            write_power_spectrum
********************************************/
void write_power_spectrum(int iter)
{
    // This function does as the name says and writes the binned, normalized ps to a file.

    int i;
    FILE *fd;
    char f_out[256];

    // Build the file name for writing
    sprintf(f_out, "ps_data_task-%d_iter-%d_snap-%03d.txt", thistask, iter, snapnum);

    // Open file for writing
    if(!(fd = fopen(f_out, "w")))
    {
        printf("Error, could not open file for writing ps on task %d iter %d!\n",\
            thistask, iter);
        exit_code_local = 0;
        return;
    }

    // Write header
    fprintf(fd, "k (s/km) \t P(k) (km/s)\n");

    // Write the data
    for(i = 0; i < ncounts; i++)
    {
        fprintf(fd, "%f \t %f\n", bin_freqs[i], bin_ps[i]);
    }

    // Clean up
    fclose(fd);
}
