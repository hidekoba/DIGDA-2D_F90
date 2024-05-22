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
sdmain_ionc.o: netcdf_utilities.o
#sdmain_ionc.o: sdclay.o sdpom.o sdcaco3.o sdopal.o sdchemeq.o sdsolid.o sdliquid.o

# output
OUTPUTsg1=scla spom scar sopa sage pcla ppom pcar popa w
OUTPUTsg2=ficla fipom ficar fiopa
OUTPUTsg3=focla fopom focar foopa
OUTPUTsg4=bfcla bfpom bfcar bfopa
OUTPUTsg5=w
OUTPUTsg6=s13pom s14pom s13car s14car
OUTPUTsg7=fo13pom fo14pom fo13car fo14car
OUTPUTsg8=fi13pom fi14pom fi13car fi14car
OUTPUTsg9=bf13pom bf14pom bf13car bf14car

OUTPUTsn1=scla.nc spom.nc scar.nc sopa.nc sage.nc pcla.nc ppom.nc pcar.nc popa.nc w.nc
OUTPUTsn2=ficla.nc fipom.nc ficar.nc fiopa.nc

OUTPUTlg1=cnit cdic calk cdo2 csil cco32 satsil satco32
OUTPUTlg2=finit fidic fialk fido2 fisil
OUTPUTlg3=fonit fodic foalk fodo2 fosil
OUTPUTlg4=c13dic c14dic cdc13 cdc14
OUTPUTlg5=fi13dic fi14dic
OUTPUTlg6=fo13dic fo14dic

OUTPUTln1=cnit.nc cdic.nc calk.nc cdo2.nc csil.nc

# cleanup
clean:
#rm -f $(OBJS) $(OBJS:.o=.mod)
	rm -f $(OUTPUTsg1) $(OUTPUTsg2) $(OUTPUTsg3) $(OUTPUTsg4) $(OUTPUTsg5) $(OUTPUTsg6) $(OUTPUTsg7) $(OUTPUTsg8) $(OUTPUTsg9) $(OUTPUTlg1) $(OUTPUTlg2) $(OUTPUTlg3) $(OUTPUTlg4) $(OUTPUTlg5) $(OUTPUTlg6) fort.* dcl*
#	rm -f $(OUTPUTsn1) $(OUTPUTsn2) $(OUTPUTln1)
	rm -rf tmpdir ; mkdir tmpdir ; mv $(OUTPUTsn1) $(OUTPUTsn2) $(OUTPUTsn3) $(OUTPUTsn4) $(OUTPUTsn5) $(OUTPUTsn6) $(OUTPUTsn7) $(OUTPUTsn8) $(OUTPUTsn9) $(OUTPUTln1) $(OUTPUTln2) $(OUTPUTln3) $(OUTPUTln4) $(OUTPUTln5) $(OUTPUTln6) tmpdir

cleanall:
#rm -f $(OBJS) $(OBJS:.o=.mod)
	rm -f $(PROGS) $(OUTPUTsg1) $(OUTPUTsg2) $(OUTPUTsg3) $(OUTPUTsg4) $(OUTPUTsg5) $(OUTPUTsg6) $(OUTPUTsg7) $(OUTPUTsg8) $(OUTPUTsg9) $(OUTPUTlg1) $(OUTPUTlg2) $(OUTPUTlg3) $(OUTPUTlg4) $(OUTPUTlg5) $(OUTPUTlg6) *.o *.mod fort.* dcl*
	rm -f $(OUTPUTsn1) $(OUTPUTsn2) $(OUTPUTln1)
#	rm -rf tmpdir ; mkdir tmpdir ; mv $(OUTPUTsn1) $(OUTPUTsn2) $(OUTPUTsn3) $(OUTPUTsn4) $(OUTPUTsn5) $(OUTPUTsn6) $(OUTPUTsn7) $(OUTPUTsn8) $(OUTPUTsn9) $(OUTPUTln1) $(OUTPUTln2) $(OUTPUTln3) $(OUTPUTln4) $(OUTPUTln5) $(OUTPUTln6) tmpdir
