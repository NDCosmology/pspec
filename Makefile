# Makefile for pspec








#--------------------------------------- Compile Time Options
#OPT += -DDEBUGGING

#--------------------------------------- Select Target Computer

#SYSTYPE="CRC-Opteron-long"
#SYSTYPE="phillips"
SYSTYPE="jared_home"

#--------------------------------------- System Specifics

ifeq ($(SYSTYPE),"jared_home")
CC = mpicc
#OPTIMIZE  = -Wall -g3
OPTIMIZE = -Wall -O3
GSL_INCL  = -I/usr/local/include/gsl
GSL_LIBS  = -L/usr/local/lib
FFTW_INCL = -I/usr/local/include
FFTW_LIBS = -L/usr/local/lib
endif



ifeq ($(SYSTYPE),"phillips")
CC = gcc
OPTIMIZE = -Wall -g3
GSL_INCL = -I/opt/local/include
GSL_LIBS = -L/opt/local/lib -Wl
endif



ifeq ($(SYSTYPE),"CRC-Opteron-long")
CC       = /afs/crc.nd.edu/user/j/jcoughl2/.local/bin/mpicc
OPTIMIZE =  -O3 -Wall
#OPTIMIZE = -g3 -Wall
GSL_INCL = -I/afs/crc.nd.edu/user/j/jcoughl2/.local/include
GSL_LIBS = -L/afs/crc.nd.edu/user/j/jcoughl2/.local/lib
FFTW_INCL=  -I/afs/crc.nd.edu/user/j/jcoughl2/.local/include
FFTW_LIBS=  -L/afs/crc.nd.edu/user/j/jcoughl2/.local/lib
MPICHLIB =	-L/afs/crc.nd.edu/user/j/jcoughl2/.local/lib
endif



#-------------------------------------- Bookkeeping
PREFIX = ./src
OBJ_DIR = $(PREFIX)/obj

OPTIONS =  $(OPTIMIZE) $(OPT)

EXEC   = pspec_noise

OBJS   = $(OBJ_DIR)/allvars.o $(OBJ_DIR)/main.o $(OBJ_DIR)/initialize.o \
         $(OBJ_DIR)/power_spectrum.o $(OBJ_DIR)/read_spectra.o $(OBJ_DIR)/utils.o \
         $(OBJ_DIR)/write.o
   
INCL   = $(PREFIX)/allvars.h  $(PREFIX)/proto.h Makefile

INCLUDE = $(GSL_INCL) $(FFTW_INCL)

CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL)

all: $(EXEC)  

LIBS = $(GSL_LIBS) $(FFTW_LIBS) -lgsl -lgslcblas -lfftw3f -lm

#ALI Added 2/6/13  see http://www.apl.jhu.edu/Misc/Unix-info/make/make_10.html#SEC90
#                  for logic 
$(OBJ_DIR)/%.o : $(PREFIX)/%.c $(INCL)
	$(CC) $(OPTIONS) $(INCLUDE) -c $< -o $@

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) -o  $(EXEC)

clean:
	rm -f $(OBJS) *.gch 
