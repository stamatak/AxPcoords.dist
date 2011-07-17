# Axelerated Pcoords Makefile 
# Copyright January 2007 by Alexandros Stamatakis

CC = gcc 

#FOR USE WITH GNU SCIENTIFIC LIBRARY

#REMOVED -static flag in GSL seemd to cause trouble with compilation
# on some systems

GSL  = -lm -lgsl -lgslcblas
GSL_CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -D_USE_GSL

######################################################################################################
# AMD ACML-specific FLAGS
# 
# ATTENTION ! ATTENTION ! adapt ACML and ACML_CFLAGS to your local installation/path

#ACML = -L/home/stamatak/ACML/gnu64/lib/ -lacml -lm -lg2c
ACML = -L../acml3.6.0/gnu32/lib/ -lacml -lm -lg2c
#ACML = -Lc:/msys/1.0/mingw/acml3.6.0/gnu32/lib -lacml -lm -lg2c

#ACML_CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -D_USE_LAPACK -D_USE_ACML -I/home/stamatak/ACML/gnu64/include/
ACML_CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -D_USE_LAPACK -D_USE_ACML -I../acml3.6.0/gnu32/include/
#ACML_CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -D_USE_LAPACK -D_USE_ACML -Ic:/msys/1.0/mingw/acml3.6.0/gnu32/include


######################################################################################################
# Intel MKL-specific FLAGS
# 
# ATTENTION ! ATTENTION ! adapt MKL and MKL_CFLAGS to your local installation/path

MKL = -L../intel-mkl/lib/32 -lmkl -lmkl_lapack -lm -lguide
MKL_CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -D_USE_LAPACK -D_USE_MKL -I../intel-mkl/include/

RM = rm -f

all : AxPcoordsGSL AxPcoordsACML AxPcoordsMKL

#FOR USE WITH GNU SCIENTIFIC LIBRARY

AxPcoordsGSL : AxPcoordsGSL.o
	$(CC) -o AxPcoordsGSL AxPcoordsGSL.o $(GSL)

AxPcoordsGSL.o : AxPcoords.c
	$(CC) $(GSL_CFLAGS) -c AxPcoords.c -o AxPcoordsGSL.o

#FOR USE WITH ACML library

AxPcoordsACML : AxPcoordsACML.o
	$(CC) -o AxPcoordsACML AxPcoordsACML.o $(ACML)

AxPcoordsACML.o : AxPcoords.c
	$(CC) $(ACML_CFLAGS) -c  AxPcoords.c -o AxPcoordsACML.o

#FOR USE WITH MKL LIBRARY

AxPcoordsMKL : AxPcoordsMKL.o
	$(CC) -o AxPcoordsMKL AxPcoordsMKL.o $(MKL)

AxPcoordsMKL.o : AxPcoords.c
	$(CC) $(MKL_CFLAGS) -c  AxPcoords.c -o AxPcoordsMKL.o


clean : 
	$(RM) *.o
	$(RM) AxPcoordsACML AxPcoordsGSL AxPcoordsMKL
	$(RM) AxPcoordsACML.exe AxPcoordsGSL.exe AxPcoordsMKL.exe
