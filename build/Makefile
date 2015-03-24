# Makefile for EFDC
# Chris Chartrand & Chris Flanary
# 6/5/2014

# Name of executable file to make
EXEC=efdcS_ifort

# Set compiler flags
FC= ifort
FFLAGS=-fast
NETCDF_INCLDIR=-I/usr/local/netcdf/include
NETCDF_LIBDIR=-L/usr/local/netcdf/lib -lnetcdff

# Module files needed for linking
MODOBJS := Var_Global_Mod.o DRIFTER-SCJ.o WINDWAVE.o

# Define object file directory
OBJSDIR=./Build

# ?
MODOBJS := $(addprefix $(OBJSDIR)/,$(MODOBJS))

# ?
MOD := -J $(OBJSDIR)

# ?
ALL_SRCS := $(wildcard *.f)
ALL_SRCS := $(ALL_SRCS) $(wildcard *.f90)
ALL_SRCS := $(ALL_SRCS) $(wildcard *.F90)
ALL_SRCS := $(ALL_SRCS) $(wildcard *.for)
ALL_SRCS := $(ALL_SRCS) $(wildcard *.FOR)

# ?
MAKE=make

# ?
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
