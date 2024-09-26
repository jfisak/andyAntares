#Compiler settings
#F90 = gfortran
# F90=mpif90
# F90=ifort
F90=mpifort.mpich
FCFLAGS = -g -cpp -O0 -Dmpi=1 -check all -llapack 
FCFLAGS = -g -cpp -O0 -fbounds-check -Dmpi=1 -fcheck=all -Wall -llapack # -pg -fprofile-arcs# -ffpe-trap=zero,overflow,invalid,underflow

#Variables
progname = main
PROJECT = $(progname).run

linkdate=$(shell date)
linkuser=$(shell whoami)
linkhost=$(shell hostname)

#Objects
MODULES=types.o rates_i.o rates_k.o rates_r.o constants.o dummypacket.o virt_gridAB.o
MODELS= read_1D_model.o read_2D_model.o read_3D_model.o read_1D_araya.o read_3D_nico.o\
        read_3D_pseudo3D.o read_2D_peku.o read_2D_basic.o
RATES= r_kappa_cont.o r_kappa_line.o i_ion_recomb.o i_radtrans.o i_radion.o i_colion.o \
       i_coltrans.o i_freq_recomb.o  k_freq_ff.o k_freq_fb.o read_photcs.o cool_excit.o \
       cool_ff.o cool_ionization.o cool_fb.o 
ATOMIC=read_levels.o read_transitions.o read_e_nd.o 
CONNECTION=connection_prop_model_grid.o connect_2D_peku.o connect_2D_basic.o
VELOCITY=velo.o velo_vector.o vel_discrete_points.o vel_interpolation.o vel_vector_interpolation.o \
         lin_interpolation_3D.o lin_interpolation.o vel_pseudo3D_model.o
PACKETS=do_rpackage.o do_ipackage.o do_rpackage_event.o do_kpackage.o do_dpackage.o \
 do_vpackage.o packet_dynamics.o
SPECTRA=do_brtm_spectrum.o do_spectrum.o
GRIDS=propmodgrid_diagnostics.o setup_model_grid.o update_grid.o setup_propgrid.o create_dynamical_grid_cells.o \
      propgrid_dist.o save_propmod_grid.o read_propmod_grid.o
OBJECTS = main.o sargc.o sargp.o sargv.o idx.o read_input.o 	  \
          random_unitvector.o ran2.o boundary3.o bound_dist.o \
          change_cell.o move_package.o emit_rpackage.o	\
          init_photsphere.o random_unitvector2.o ran2_class.o save_rates.o\
          random_unitvector1.o doppler_factor.o oct_neighbors.o \
          vec_length.o angle_aberration.o freq_from_planck.o \
           event_dist.o get_package_model_index.o cross_product.o \
          update_packages.o r_choose_line.o brtm.o \
          read_composition.o read_populations.o \
          sorting-new.o cmf_freq.o n_closest_points_2D.o\
          saha_boltzmann_factor.o ionization_fraction.o f_edens.o round_number.o\
          find_e_nd.o part_fun.o update_estimators.o freq_from_file.o rob.o \
          acc_rej_montecarlo.o virtual_points.o \
          find_dyn_cell1.o seek_nclosest_points.o\
          divide_cell_8.o divide_cell_ijk.o \
	  next_cell_down.o next_cell_up.o populations.o  \
	  gamma_function.o exp_int_func.o find_populations.o\
	  flux_function.o find_element_index.o analyse_input.o \
	  warning.o calc_tau.o check_pop.o \
	  save_output.o gauntff.o counters.o saha_factor.o find_photion_elindex.o \
	  mpi_distribute_estimators.o save_temp_packs.o find_unfinished_run.o \
	  lte_pops.o photosphere_interaction.o roverw.o \
          next_line_bluered.o resonance_distance2.o \
	  oct_virtcube.o diffusion_approximation.o d_choosenextcell.o \
	  $(MODELS) $(RATES) $(ATOMIC) $(MODULES) $(CONNECTION) $(VELOCITY) $(PACKETS) $(SPECTRA) \
	  $(GRIDS)
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
constants.o: constants.f90
	$(F90) $(FCFLAGS) -o $@ -c $<
dummypacket.o: dummypacket.f90 
	$(F90) $(FCFLAGS) -o $@ -c $<
virt_gridAB.o: virt_gridAB.f90 
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
