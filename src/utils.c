/*****************************************************************************************
Title:   utils.c
Purpose: Contains miscellaneous functions
Notes:
*****************************************************************************************/
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <fftw3.h>
#include <mpi.h> 
#include "allvars.h"
#include "proto.h"



/********************************************
             check_exit_code
********************************************/
void check_exit_code(int exit_code_local, POWER_SPECTRUM *ps)
{
    // This function performs an Allreduce product operation on each processor's exit
    // code. If this product is zero then it means that one of the processor's failed
    // and we need each to call MPI_Finalize and then quit
    MPI_Allreduce(&exit_code_local, &exit_code_global, 1, MPI_INT, MPI_PROD, \
        MPI_COMM_WORLD);

    if(exit_code_global == 0)
    {
        free_globals(ps);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}



/********************************************
               free_globals
********************************************/
void free_globals(POWER_SPECTRUM *ps)
{
    // This function checks to see if the applicable globals have been allocated (by
    // comparing them with NULL) and then freeing them if necessary.

    if(bin_freqs != NULL)
    {
        free(bin_freqs);
    }

    if(bin_ps != NULL)
    {
        free(bin_ps);
    }

    if(counts != NULL)
    {
        free(counts);
    }

    if(rng != NULL)
    {
        gsl_rng_free(rng);
    }

    if(ps->vel != NULL)
    {
        free(ps->vel);
    }

    if(ps->flux != NULL)
    {
        free(ps->flux);
    }

    if(ps->ps_amp != NULL)
    {
        free(ps->ps_amp);
    }
}
