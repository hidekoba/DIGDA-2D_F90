##############################
# Targets
##############################

PROGS=a.out

##############################
# SRC
##############################

SRC=io_netcdf.f90 netcdf_utilities.f90 sdmain_ionc.F90 sdclay.F90 sdpom.F90 sdcaco3.F90 sdopal.F90 sdchemeq.F90 sdsolid.F90 sdliquid.F90
#SRC=netcdf_utilities.f90 sdmain_ionc.f90 io_netcdf.f90

OBJS=sdclay.o sdpom.o sdcaco3.o sdopal.o sdchemeq.o sdsolid.o sdliquid.o sdmain_ionc.o netcdf_utilities.o io_netcdf.o

#OBJS=$(SRC:.f90=.o)
#OBJS=$(SRC:.F=.o) $(SRC:.f90=.o)
#OBJS=$(SRC:.f90=.o) $(SRC:.F=.o)

##############################
# Compiler
##############################

#FC=gfortran
#FC=dclfrt
#FC=ifort
#FC=f77
#FC=f90
#FC=mpif90
FC=mpiifort

##############################
# Option
##############################

#FCFLAGS=-Wall -O2
#FCFLAGS=-O3 -u -convert big_endian -fbounds-check
#FCFLAGS=-O3 -u -convert big_endian -traceback -pg -shared-intel -mcmodel=large -warn declarations -check uninit -warn unused
#FCFLAGS=-convert big_endian -traceback -mcmodel=large

FFLAGS=-O3 -u -convert big_endian -traceback -I/home/hidekoba/local/netcdf4-intel-fortran/4.6.0/include
LDFLAGS=-L/home/hidekoba/local/netcdf4-intel-fortran/4.6.0/lib -lnetcdff -L/home/hidekoba/local/netcdf4-intel/4.9.2/lib -lnetcdf -lnetcdf -lm

##############################
# Suffixes
##############################

.SUFFIXES: .f90 .F90

all:$(PROGS)

# Build executable file
$(PROGS): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for source files to object files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
%.o: %.F90
	$(FC) $(FFLAGS) -c $<

# .f90 file dependencies
io_netcdf.o: netcdf_utilities.o
#sdmain_ionc.o: sdclay.o sdpom.o sdcaco3.o sdopal.o sdchemeq.o sdsolid.o sdliquid.o

# output
OUTPUTs1=scla spom scar sopa sage pcla ppom pcar popa w
OUTPUTs2=scla spom scar sopa sage pcla ppom pcar popa w
OUTPUTs3=ficla fipom ficar fiopa
OUTPUTs4=focla fopom focar foopa
OUTPUTs5=bfcla bfpom bfcar bfopa
OUTPUTs6=w
OUTPUTs7=s13pom s14pom s13car s14car
OUTPUTs8=fo13pom fo14pom fo13car fo14car
OUTPUTs9=fi13pom fi14pom fi13car fi14car
OUTPUTs10=bf13pom bf14pom bf13car bf14car
OUTPUTl1=cnit cdic calk cdo2 csil cco32 satsil satco32
OUTPUTl2=finit fidic fialk fido2 fisil
OUTPUTl3=fonit fodic foalk fodo2 fosil
OUTPUTl4=c13dic c14dic cdc13 cdc14
OUTPUTl5=fi13dic fi14dic
OUTPUTl6=fo13dic fo14dic

# cleanup
clean:
#rm -f $(OBJS) $(OBJS:.o=.mod)
	rm -f $(PROGS) $(OUTPUTs1) $(OUTPUTs2) $(OUTPUTs3) $(OUTPUTs4) $(OUTPUTs5) $(OUTPUTs6) $(OUTPUTs7) $(OUTPUTs8) $(OUTPUTs9) $(OUTPUTs10) $(OUTPUTl1) $(OUTPUTl2) $(OUTPUTil3) $(OUTPUTl4) $(OUTPUTl5) $(OUTPUTl6) *.o *.mod fort.* dcl*
