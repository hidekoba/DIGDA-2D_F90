##############################
# Targets
##############################

PROGS=a.out

##############################
# SRC
##############################

SRC=config.f90 sddim.f90 zcfile.f90 zchem.f90 zcphys.f90 sdclay.f90 sdpom.f90 sdcaco3.f90 sdopal.f90 sdchemeq.f90 sdsolid.f90 sdliquid.f90 sdmain_ionc.f90 netcdf_utilities.f90 io_netcdf.f90
OBJS=config.o sddim.o zcfile.o zchem.o zcphys.o sdclay.o sdpom.o sdcaco3.o sdopal.o sdchemeq.o sdsolid.o sdliquid.o sdmain_ionc.o netcdf_utilities.o io_netcdf.o


#SRC=config.f90 sddim.f90 zcfile.f90 zchem.f90 zcphys.f90 sdclay.f90 sdpom.f90 sdcaco3.f90 sdopal.f90 sdchemeq.f90 sdsolid.f90 sdliquid.f90 sdmain_ionc.f90 netcdf_utilities.f90 io_gtool.f90
#OBJS=config.o sddim.o zcfile.o zchem.o zcphys.o sdclay.o sdpom.o sdcaco3.o sdopal.o sdchemeq.o sdsolid.o sdliquid.o sdmain_ionc.o gtool_utilities.o io_gtool.o


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

FFLAGS=-O3 -u -convert big_endian -traceback -mcmodel=large -I/home/hidekoba/local/netcdf4-intel-fortran/4.6.0/include
LDFLAGS=-L/home/hidekoba/local/netcdf4-intel-fortran/4.6.0/lib -lnetcdff -L/home/hidekoba/local/netcdf4-intel/4.9.2/lib -lnetcdf -lnetcdf -lm

##############################
# Suffixes
##############################

.SUFFIXES: .f90

all:$(PROGS)

# Build executable file
$(PROGS): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for source files to object files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# .f90 file dependencies
io_netcdf.o: netcdf_utilities.o
#sdmain_ionc.o: sdclay.o sdpom.o sdcaco3.o sdopal.o sdchemeq.o sdsolid.o sdliquid.o

# output
OUTPUTs1=scla spom scar sopa sage pcla ppom pcar popa w
OUTPUTs2=scla.nc spom.nc scar.nc sopa.nc sage.nc pcla.nc ppom.nc pcar.nc popa.nc w.nc
OUTPUTs3=ficla fipom ficar fiopa
OUTPUTs4=focla fopom focar foopa
OUTPUTs5=bfcla bfpom bfcar bfopa
OUTPUTs6=w
OUTPUTs7=s13pom s14pom s13car s14car
OUTPUTs8=fo13pom fo14pom fo13car fo14car
OUTPUTs9=fi13pom fi14pom fi13car fi14car
OUTPUTs10=bf13pom bf14pom bf13car bf14car
OUTPUTl1=cnit cdic calk cdo2 csil cco32 satsil satco32
OUTPUTl2=cnit.nc cdic.nc calk.nc cdo2.nc csil.nc cco32.nc satsil.nc satco32.nc
OUTPUTl3=finit fidic fialk fido2 fisil
OUTPUTl4=fonit fodic foalk fodo2 fosil
OUTPUTl5=c13dic c14dic cdc13 cdc14
OUTPUTl6=fi13dic fi14dic
OUTPUTl7=fo13dic fo14dic

# cleanup
clean:
#rm -f $(OBJS) $(OBJS:.o=.mod)
	rm -f $(PROGS) $(OUTPUTs1) $(OUTPUTs2) $(OUTPUTs3) $(OUTPUTs4) $(OUTPUTs5) $(OUTPUTs6) $(OUTPUTs7) $(OUTPUTs8) $(OUTPUTs9) $(OUTPUTs10) $(OUTPUTl1) $(OUTPUTl2) $(OUTPUTil3) $(OUTPUTl4) $(OUTPUTl5) $(OUTPUTl6) $(OUTPUTl7) *.o *.mod fort.* dcl*
