include Mkinclude
include Sysdep.$(SYSTEM)$(PARALLEL)
OBJS=digda_config.o digda_mpi.o digda_grid.o digda_parameter.o digda_timestep.o zcfile.o zchem.o zcphys.o digda_clay.o digda_pom.o digda_carbonate.o digda_silicate.o digda_chemeq.o digda_step_solid.o digda_step_liquid.o netcdf_utilities.o gtool_utilities.o digda_cprdc_ionc.o digda_dataio_netcdf.o digda_ased.o

#.f90.o:
#       $(FC) $(FOPTION) -c $<
#
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

#all:lib $(TARGET)
all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FOPTION) -o $@ $(OBJS)
	$(RM) $(OBJS)

clean:
	$(RM) *~ *.o *.mod *.lst

reset:clean
	$(RM) ../lib/*
	$(RM) ../bin/*
