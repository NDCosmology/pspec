/*****************************************************************************************
Title: allvars.c
Purpose: Contains actual declarations and initializations of global variables
Notes:
*****************************************************************************************/
#include <gsl/gsl_rng.h>
#include "allvars.h"



// MPI Parameters
int ntasks = 0;
int thistask = 0;

// Parameter file parameters
char file_base[256];
int nspec_files =0;
int nfiles_to_use = 0;
int npspec = 0;
int use_sdss = 0;
int use_spec_k = 0;
int use_mcdonald = 0;
int use_all = 0;
int use_noise = 0;
float redshift = 0.0;

// General globals
int max_iters = 10000;
int exit_code_local = 1;
int exit_code_global = 1;
int npspec_local = 0;
int snapnum = 0;
gsl_rng *rng = NULL;

// Debugging params
#ifdef DEBUGGING
    int *py_file_ids;
#endif

// Power spectrum globals
float mean_flux = 0.0;
float *bin_freqs = NULL;
float *bin_ps = NULL;
float upper_k_limit = 0.0;
int *counts = NULL;
int ncounts = 30;
