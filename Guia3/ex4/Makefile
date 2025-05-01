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
SRC = $(MODDIR)/precision.f90 \
      $(MODDIR)/formats.f90 \
      $(MODDIR)/constantes.f90 \
      $(MODDIR)/mzranmod.f90 \
      $(MODDIR)/mzranmod_threadsafe.f90 \
      $(MODDIR)/mtmod.f90 \
      $(MODDIR)/funciones.f90 \
      $(MODDIR)/subrutinas.f90 \
      ex4.f90

# Object and dependency files
OBJ = $(patsubst %.f90, $(OBJDIR)/%.o, $(notdir $(SRC)))
DEP = $(OBJ:.o=.d)

# Target
TARGET = ex4.exe

# Default rule
all: $(TARGET)

# Linking step
$(TARGET): $(OBJ)
	$(FC) $^ -o $@ $(FFLAGS)

# Compilation with dependency generation
$(OBJDIR)/%.o: $(MODDIR)/%.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(DEPFLAGS) -J$(OBJDIR)

$(OBJDIR)/ex4.o: ex4.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(DEPFLAGS) -J$(OBJDIR)

# Include generated dependency files, if present
-include $(DEP)

# Clean rule
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*.d $(TARGET)

.PHONY: all clean
