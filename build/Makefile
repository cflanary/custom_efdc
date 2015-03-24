# Makefile for Linux Gfortran EFDC
#
#
# 6/5/2014 Chris Chartrand Created makefile
# 6/25/2014 Chris Flanary Modified for EFDC Gfortran build
# 6/30/2014 Chris Flanary Added support for NetCDF

# Name of executable file to make
EXEC=/usr/local/efdc/alc_bcsaS

# Set compiler flags
FC= gfortran
FFLAGS=-O3 -w -ffree-line-length-none -fbounds-check
NETCDF_INCLDIR=-I/usr/local/netcdf/include
NETCDF_LIBDIR=-L/usr/local/netcdf/lib -lnetcdff

# Module files needed for linking, compiled first
MODOBJS := Var_Global_Mod.o DRIFTER-SCJ.o WINDWAVE.o

# Define object file build directory
OBJSDIR=./Build

# Add module files to MODOBJS for build directory 
MODOBJS := $(addprefix $(OBJSDIR)/,$(MODOBJS))

# Make modules simultaneously
MOD := -J $(OBJSDIR)

# Create ALL_SRCS array with all source code from any Fortran format
ALL_SRCS := $(wildcard *.f)
ALL_SRCS := $(ALL_SRCS) $(wildcard *.f90)
ALL_SRCS := $(ALL_SRCS) $(wildcard *.F90)
ALL_SRCS := $(ALL_SRCS) $(wildcard *.for)
ALL_SRCS := $(ALL_SRCS) $(wildcard *.FOR)

# Set command to run Makefile
MAKE=make

# Convert all Fortran source code into object files
OBJS := $(ALL_SRCS:%.f=$(OBJSDIR)/%.o)
OBJS := $(OBJS:%.f90=$(OBJSDIR)/%.o)
OBJS := $(OBJS:%.F90=$(OBJSDIR)/%.o)
OBJS := $(OBJS:%.for=$(OBJSDIR)/%.o)
OBJS := $(OBJS:%.FOR=$(OBJSDIR)/%.o)

# Link executable
$(EXEC): dircheck $(MODOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS) $(NETCDF_INCLDIR) $(NETCDF_LIBDIR)

# Compile object files from source code
$(OBJSDIR)/%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJSDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@ $(NETCDF_INCLDIR) $(NETCDF_LIBDIR)

$(OBJSDIR)/%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJSDIR)/%.o: %.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJSDIR)/%.o: %.FOR
	$(FC) $(FFLAGS) -c $< -o $@

# Create build directory
dircheck:
	@mkdir -p $(OBJSDIR)

# Clean after make is complete
clean:
	-rm $(OBJSDIR)/*.o
	-rm *.mod
