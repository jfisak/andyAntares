#Compiler settings
F90 = gfortran
FCFLAGS =  -g -O2 -cpp -fbounds-check #-mp

#Variables
progname = main
PROJECT = $(progname).exe

linkdate=$(shell date)
linkuser=$(shell whoami)
linkhost=$(shell hostname)

#Objects
OBJECTS = main.o sargc.o sargp.o sargv.o idx.o read_input.o modul.o	  \
          random_unitvector.o ran2.o boundary3.o find_dist.o \
          change_cell.o move_package.o emit_rpackage.o			  \
          init_photsphere.o random_unitvector2.o 		  \
          random_unitvector1.o read_1D_model.o doppler_factor.o		  \
          vec_length.o velo.o angle_aberration.o freq_from_planck.o	  \
          do_rpackage.o event_dist.o get_package_model_index.o		  \
          do_rpackage_event.o update_packages.o do_spectrum.o             \
          read_composition.o read_atomic_data2.o read_transitions2.o      \
          sorting-new.o setup_model_grid.o update_grid.o                  \
          saha_boltzmann_factor.o ionization_fraction.o f_edens.o         \
          find_e_nd.o part_fun.o update_estimators.o freq_from_file.o     \
          acc_rej_montecarlo.o virtual_particles.o setup_grid2.o          \
          create_dynamical_grid_cells.o find_dyn_cell1.o \
          connection_prop_model_grid.o 		\
	  next_cell_down.o next_cell_up.o read_2D_model.o

#Rules
all : $(PROJECT)

modul.o: modul.f90 
	$(F90) $(FCFLAGS) -o $@ -c $<

%.for:
	@printf "      PROGRAM MAIN$(progname)\n"  > compiletime.for
	@printf "C***  Provide Link data for possible use in the programm\n" >> compiletime.for
	@printf "      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60\n" >> compiletime.for 
	@printf "      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST\n">>compiletime.for  
	@printf "      LINK_DATE = '$(linkdate)'\n" >>compiletime.for  
	@printf "      LINK_USER = '$(linkuser)'\n" >>compiletime.for 
	@printf "      LINK_HOST = '$(linkhost)'\n" >>compiletime.for 
	@printf "             \n" >> compiletime.for
	@printf "      CALL $*\n" >> compiletime.for
	@printf "      END\n" >> compiletime.for

%.o: %.f90
	$(F90) $(FCFLAGS)  -c $<

%.exe : $(OBJECTS) %.for
	$(F90) -o $(PROJECT) compiletime.for $(OBJECTS) $(FCFLAGS)



# Utility targets
.PHONY: clean veryclean

clean :
	rm -f $(OBJECTS) *mod compiletime.for

veryclean : clean
	rm -f $(PROJECT) 
