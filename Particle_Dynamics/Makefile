# Compiler and flags
FC = gfortran
FFLAGS = -cpp -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=2 \
         -ffpe-trap=invalid,overflow,zero -Wall -fcheck=all -g -fbacktrace \
         -fopenmp -O2
DEPFLAGS = -MMD -MP   # Generates .d files for dependencies

# Directories
MODDIR = ./modulos
OBJDIR = obj

# Ensure obj directory exists
$(shell mkdir -p $(OBJDIR))

# Source files
SRC = 	./modulos/precisionMod.f90\
		./modulos/mzranMod.f90\
		./modulos/mzran_threadsafeMod.f90\
		./modulos/randomMod.f90\
		./modulos/formatsMod.f90\
		./modulos/constantsMod.f90\
		./modulos/parametersMod.f90\
		./modulos/subroutinesMod.f90\
		./modulos/potentialsMod.f90\
		./modulos/observablesMod.f90\
		./modulos/thermostatsMod.f90\
		./modulos/updatePositionsMod.f90\
		./modulos/parsingMod.f90\
		./modulos/forcesMod.f90\
		./modulos/writing2filesMod.f90\
		./modulos/initializationsMod.f90\
		./modulos/integratorsMod.f90\
      	main.f90

# Object and dependency files
OBJ = $(patsubst %.f90, $(OBJDIR)/%.o, $(notdir $(SRC)))
DEP = $(OBJ:.o=.d)

# Target
TARGET = run.exe

# Default rule
all: $(TARGET)

# Linking step
$(TARGET): $(OBJ)
	$(FC) $^ -o $@ $(FFLAGS)

# Compilation with dependency generation
$(OBJDIR)/%.o: $(MODDIR)/%.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(DEPFLAGS) -J$(OBJDIR)

$(OBJDIR)/main.o: main.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(DEPFLAGS) -J$(OBJDIR)

# Include generated dependency files, if present
-include $(DEP)

# Clean rule
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*.d $(TARGET)

.PHONY: all clean
