/*****************************************************************************************
Title: proto.h
Purpose: Contains function prototypes
Notes:
*****************************************************************************************/
#ifndef ALLVARS_H
    #include "allvars.h"
#endif



/********************************************
                 initialize.c
********************************************/
int checkforerror(char *, int, char *);
void initialize_run(int, char **, POWER_SPECTRUM *);



/********************************************
              power_spectrum.c
********************************************/
void bin_power_spectrum(POWER_SPECTRUM);
void normalize_power_spectrum(void);
void power_spectrum(POWER_SPECTRUM *);



/********************************************
               read_spectra.c
********************************************/
void read_spectrum_file(char *, POWER_SPECTRUM *);
void read_all_spectrum_files(char *, POWER_SPECTRUM *, int);



/********************************************
                   utils.c
********************************************/
void check_exit_code(int, POWER_SPECTRUM *);
void free_globals(POWER_SPECTRUM *);



/********************************************
                   write.c
********************************************/
void write_power_spectrum(int);
