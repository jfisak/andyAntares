#Compiler settings
#F90 = gfortran
# F90=mpif90
# F90=ifort
F90=mpifort.mpich
FCFLAGS = -g -cpp -O0 -Dmpi=1 -check all -llapack 
FCFLAGS = -g -cpp -O0 -fbounds-check -Dmpi=1 -fcheck=all -Wall -pg -llapack # -fprofile-arcs# -ffpe-trap=zero,overflow,invalid,underflow

#Variables
progname = main
PROJECT = $(progname).run

linkdate=$(shell date)
linkuser=$(shell whoami)
linkhost=$(shell hostname)

#Objects
OBJECTS = main.o sargc.o sargp.o sargv.o idx.o read_input.o types.o	  \
          random_unitvector.o ran2.o boundary3.o bound_dist.o rates_i.o\
          change_cell.o move_package.o emit_rpackage.o	rates_k.o rates_r.o \
          init_photsphere.o random_unitvector2.o ran2_class.o save_rates.o\
          random_unitvector1.o read_1D_model.o doppler_factor.o doppler_factor2.o oct_neighbors.o \
          vec_length.o velo.o angle_aberration.o freq_from_planck.o lin_interpolation.o \
          do_rpackage.o do_ipackage.o event_dist.o get_package_model_index.o  \
          do_rpackage_event.o update_packages.o do_kpackage.o r_choose_line.o \
          read_composition.o read_levels.o read_transitions.o read_populations.o \
          sorting-new.o setup_model_grid.o update_grid.o i_ion_recomb.o cmf_freq.o \
          saha_boltzmann_factor.o ionization_fraction.o f_edens.o round_number.o\
          find_e_nd.o part_fun.o update_estimators.o freq_from_file.o rob.o \
          acc_rej_montecarlo.o virtual_points.o setup_propgrid.o velo_vector.o \
          create_dynamical_grid_cells.o find_dyn_cell1.o do_spectrum.o read_3D_nico.o\
          connection_prop_model_grid.o divide_cell_8.o divide_cell_ijk.o do_dpackage.o \
	  next_cell_down.o next_cell_up.o read_2D_model.o populations.o  \
	  i_coltrans.o gamma_function.o exp_int_func.o cool_excit.o read_e_nd.o \
	  read_photcs.o i_radtrans.o i_radion.o i_colion.o find_populations.o \
	  flux_function.o find_element_index.o analyse_input.o read_3D_model.o \
	  r_kappa_cont.o i_freq_recomb.o cool_ff.o warning.o k_freq_ff.o \
	  cool_ionization.o cool_fb.o k_freq_fb.o check_pop.o r_kappa_line.o \
	  save_output.o gauntff.o counters.o saha_factor.o find_photion_elindex.o \
	  mpi_distribute_estimators.o save_temp_packs.o find_unfinished_run.o \
	  lte_pops.o photosphere_interaction.o roverw.o vel_discrete_points.o \
          next_line_bluered.o resonance_distance2.o read_3D_pseudo3D.o save_propmod_grid.o \
	  read_propmod_grid.o oct_virtcube.o diffusion_approximation.o d_choosenextcell.o \
	  read_1D_araya.o calc_tau.o
# end of procedures

#Rules
all : $(PROJECT)

types.o: types.f90 
	$(F90) $(FCFLAGS) -o $@ -c $<
rates_i.o: rates_i.f90
	$(F90) $(FCFLAGS) -o $@ -c $<
rates_k.o: rates_k.f90
	$(F90) $(FCFLAGS) -o $@ -c $<
rates_r.o: rates_r.f90
	$(F90) $(FCFLAGS) -o $@ -c $<
ran2_class.o: ran2_class.f90 
	$(F90) $(FCFLAGS) -o $@ -c $<
counters.o: counters.f90
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

%.run : $(OBJECTS) %.for
	$(F90) -o $(PROJECT) compiletime.for $(OBJECTS) $(FCFLAGS)



# Utility targets
.PHONY: clean veryclean

clean :
	rm -f $(OBJECTS) *mod compiletime.for

veryclean : clean
	rm -f $(PROJECT) 
