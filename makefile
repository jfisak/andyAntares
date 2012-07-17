F90 = f95
FFLAGS = -fdump-core -g -O2 #-m64


#-i8 -r8 -traceback -O3 #-mp

# F90 = gfortran
# FFLAGS = -static -fdefault-real-8 -fdefault-integer-8 -m32

##(cd src/)

progname = main
PROJECT = $(progname).exe

# SOURCES = mainopdat.f
OBJECTS = main.o sargc.o sargp.o sargv.o idx.o read_input.o modul.o random_unitvector.o      \
          ran2.o init_photonpack.o propagation.o boundary.o do_event.o change_cell.o         \
          move_package.o emit_rpackage.o init_photsphere.o random_unitvector2.o setup_grid.o \
          random_unitvector1.o 

#OBJECTS = mainopdat.o


linkdate=$(shell date)
linkuser=$(shell whoami)
linkhost=$(shell hostname)

# Syntax:
# ZIEL : ABHAENGIGKEIT
# <TAB> BEFEHL

all : $(PROJECT)

.SECONDARY: $(OBJECTS)

%.exe : $(OBJECTS) %.for
#	echo null > \(
#	echo null > AS_NEEDED
	$(F90) -o $(PROJECT)   main$*.for $(OBJECTS) $(FFLAGS)

%.o: %.f90
	$(F90) $(FFLAGS)  -c $? 

%.for:
	@printf "      PROGRAM MAIN$(progname)\n"  > main$*.for
	@printf "C***  Provide Link data for possible use in the programm\n" >> main$*.for
	@printf "      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60\n" >> main$*.for 
	@printf "      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST\n">>main$*.for  
	@printf "      LINK_DATE = '$(linkdate)'\n" >>main$*.for  
	@printf "      LINK_USER = '$(linkuser)'\n" >>main$*.for 
	@printf "      LINK_HOST = '$(linkhost)'\n" >> main$*.for 
	@printf "             \n" >> main$*.for
	@printf "      CALL $*\n" >> main$*.for
	@printf "      END\n" >> main$*.for

clean :
	rm -f $(OBJECTS)
	rm -f $(PROJECT) 
