# Makefile created by mkmf $Id: mkmf,v 1.1.1.1 2008-08-28 01:48:18 carbrevi Exp $ 


# choose fortran compiler
#FC = gfortran
#LD = gfortran

FC = ifort
LD = ifort


# Fortran compiler flags
#intel fortran compiler debug mode: 
FFLAGS  = -g -O0 -fpe:0 -warn declarations -warn unused -warn ignore_loc -warn truncated_source -traceback -check all -implicitnone -openmp
LDFLAGS = -mkl
#intel fortran compiler optimized mode: 
#FFLAGS  = -O3 -fpe:0 -implicitnone -fast -ipo -xHost -parallel -openmp
#LDFLAGS = -mkl

#gfortran compiler debug mode: 
#FFLAGS  = -Wall -g -fbounds-check
#LDFLAGS = 
#gfortran compiler optimized mode: 
#FFLAGS  = -O2 -s -fomit-frame-pointer -fexpensive-optimizations -ffast-math
#LDFLAGS = 


.DEFAULT:
	-touch $@
all: a.out
ausm.o: ./ausm.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ausm.f90
eulerflux.o: ./eulerflux.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./eulerflux.f90
explicit_euler.o: ./explicit_euler.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./explicit_euler.f90
harten.o: ./harten.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./harten.f90
implicitbw.o: ./implicitbw.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./implicitbw.f90
jacobian.o: ./jacobian.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./jacobian.f90
laxwendroff.o: ./laxwendroff.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./laxwendroff.f90
maccormack.o: ./maccormack.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./maccormack.f90
main.o: ./main.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./main.f90
output.o: ./output.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./output.f90
preproc.o: ./preproc.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./preproc.f90
roe.o: ./roe.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./roe.f90
shared_vars.o: ./shared_vars.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./shared_vars.f90
steger_warming.o: ./steger_warming.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./steger_warming.f90
vanLeer.o: ./vanLeer.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./vanLeer.f90
SRC = ./main.f90 ./steger_warming.f90 ./harten.f90 ./preproc.f90 ./roe.f90 ./output.f90 ./maccormack.f90 ./vanLeer.f90 ./shared_vars.f90 ./explicit_euler.f90 ./laxwendroff.f90 ./implicitbw.f90 ./jacobian.f90 ./ausm.f90 ./eulerflux.f90
OBJ = main.o steger_warming.o harten.o preproc.o roe.o output.o maccormack.o vanLeer.o shared_vars.o explicit_euler.o laxwendroff.o implicitbw.o jacobian.o ausm.o eulerflux.o
clean: neat
	-rm -f .cppdefs $(OBJ) *.mod a.out *.pdf *.out
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
a.out: $(OBJ) 
	$(LD) $(OBJ) -o a.out  $(LDFLAGS)
