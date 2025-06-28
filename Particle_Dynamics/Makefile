# Compiler and flags
FC = gfortran
FFLAGS = -cpp -O2 -Wall -g -fbacktrace -fcheck=all -fopenmp \
         -ffpe-trap=invalid,overflow,zero -fno-range-check
DEPFLAGS = -MMD -MP

# Directories
SRCDIR = src
MODDIR = $(SRCDIR)/modulos
OBJDIR = obj

# Ensure obj directory exists
$(shell mkdir -p $(OBJDIR))

# List source files in correct dependency order
SRC = 	precisionMod.f90 \
        mzranMod.f90 \
        mzran_threadsafeMod.f90 \
        randomMod.f90 \
        formatsMod.f90 \
        constantsMod.f90 \
        parametersMod.f90 \
        subroutinesMod.f90 \
        potentialsMod.f90 \
        observablesMod.f90 \
        thermostatsMod.f90 \
        updatePositionsMod.f90 \
        parsingMod.f90 \
        forcesMod.f90 \
        writing2filesMod.f90 \
        initializationsMod.f90 \
        integratorsMod.f90 \
        main.f90

# Prefix source files with full path
SRCFULL = $(addprefix $(MODDIR)/, $(filter-out main.f90, $(SRC))) $(SRCDIR)/main.f90
OBJ = $(addprefix $(OBJDIR)/, $(SRC:.f90=.o))
DEP = $(OBJ:.o=.d)

# Target name
TARGET = run.exe

# Default rule
all: $(TARGET)

# Linking step
$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $^ -o $@

# Compile each object file
$(OBJDIR)/%.o: $(MODDIR)/%.f90
	$(FC) $(FFLAGS) $(DEPFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/main.o: $(SRCDIR)/main.f90
	$(FC) $(FFLAGS) $(DEPFLAGS) -J$(OBJDIR) -c $< -o $@

# Include dependency files (if they exist)
-include $(DEP)

# Clean rule
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*.d $(TARGET)

.PHONY: all clean
