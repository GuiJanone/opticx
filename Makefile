# Compiler and flags
FC = gfortran
FFLAGS = -O2

# Directories
MAINDIR = main
SRCDIR = src
BINDIR = bin
BUILDDIR = build

# Source and object files
SRC_MODULES = $(wildcard $(SRCDIR)/*.f90)
OBJ_MODULES = $(SRC_MODULES:$(SRCDIR)/%.f90=$(BUILDDIR)/%.o)
SRC_MAIN = $(MAINDIR)/opticx.f90
OBJ_MAIN = $(BINDIR)/opticx.o

# Executable
TARGET = $(BINDIR)/opticx

# Libraries
LIBS = -lopenblas -fopenmp -lgfortran -O2

# Default target to build the executable
all: $(TARGET)
	rm -f $(OBJ_MAIN)

# Rule for creating the executable
$(TARGET): $(OBJ_MODULES) $(OBJ_MAIN)
	$(FC) $(FFLAGS) $(OBJ_MODULES) $(OBJ_MAIN) -o $(TARGET) $(LIBS)

# Rule for compiling module files (*.f90) into object files (*.o)
# Use the -J flag to specify the directory for .mod files
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) -J$(BUILDDIR) -c $< $(FFLAGS) -o $@ $(LIBS)

# Dependencies between modules:
# Assume `sigma_first.f90` uses `constants_math.mod` and `parser_input_file.mod`
# These lines ensure that sigma_first.o is compiled after constants_math.o and parser_input_file.o

$(BUILDDIR)/sigma_first.o: $(BUILDDIR)/constants_math.o 
$(BUILDDIR)/sigma_first.o: $(BUILDDIR)/parser_input_file.o $(BUILDDIR)/parser_optics_xatu_dim.o $(BUILDDIR)/ome.o 
# Add more dependencies as needed for other modules
# For example:
# $(BUILDDIR)/some_other_module.o: $(BUILDDIR)/dependency_module.o

# Rule for compiling the main program file (opticx.f90) into opticx.o in the bin folder
$(BINDIR)/opticx.o: $(SRC_MAIN)
	$(FC) -J$(BUILDDIR) -c $(SRC_MAIN) $(FFLAGS) -o $(BINDIR)/opticx.o

# Clean up the object files and the executable
clean:
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.mod $(BINDIR)/opticx.o $(TARGET)
