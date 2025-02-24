# Makefile for compiling the Fortran program

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
LIBS = -lopenblas -fopenmp -lgfortran

# Default target to build the executable
all: $(TARGET)
# Remove the opticx.o after linking
	rm -f $(OBJ_MAIN) 

# Rule for creating the executable
$(TARGET): $(OBJ_MODULES) $(OBJ_MAIN)
	$(FC) $(OBJ_MODULES) $(OBJ_MAIN) -o $(TARGET) $(LIBS)

# Rule for compiling module files (*.f90) into object files (*.o)
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)

# Rule for compiling the main program file (opticx.f90) into opticx.o in bin folder
$(BINDIR)/opticx.o: $(SRC_MAIN)
	$(FC) $(FFLAGS) -c $(SRC_MAIN) -o $(BINDIR)/opticx.o

# Clean up the object files and the executable
clean:
	rm -f $(BUILDDIR)/*.o $(BINDIR)/opticx.o $(TARGET)
