# Compiler and flags
FC = gfortran

FOXCONFIG = external/FoX/FoX-config

FFLAGS = -cpp -O2 -Wall -g -fbacktrace -fcheck=all -fopenmp \
         -ffpe-trap=invalid,overflow,zero -fno-range-check \
         $(shell $(FOXCONFIG) --fcflags)

LFLAGS = $(shell $(FOXCONFIG) --libs --dom --wxml) #external/FoX/objs/lib/libFoX_dom_parse.a

DEPFLAGS = -MMD -MP

# Directories
SRCDIR = src
MODDIR = $(SRCDIR)/modules
OBJDIR = build

# Ensure obj directory exists
$(shell mkdir -p $(OBJDIR) $(OBJDIR)/RNGs $(OBJDIR)/variables $(OBJDIR)/parsing $(OBJDIR)/phases)

# List source files in correct dependency order
SRC = 	precisionMod.f90 \
        RNGs/mzranMod.f90 \
        RNGs/mzran_threadsafeMod.f90 \
        RNGs/randomMod.f90 \
        formatsMod.f90 \
        variables/constantsMod.f90 \
        variables/parametersMod.f90 \
        variables/propertiesMod.f90 \
        subroutinesMod.f90 \
        potentialsMod.f90 \
        observablesMod.f90 \
        thermostatsMod.f90 \
        updatePositionsMod.f90 \
		parsing/checkParsingMod.f90 \
		parsing/get_parsed_valueMod.f90 \
        parsing/parsingMod.f90 \
        forcesMod.f90 \
        writing2filesMod.f90 \
        integratorsMod.f90 \
        initializationsMod.f90 \
        measurementsMod.f90 \
		phases/phase_initializeMod.f90 \
		phases/phase_runMod.f90 \
		phases/phase_endMod.f90 \
        main.f90

# Prefix source files with full path
MODSFULL = $(addprefix $(MODDIR)/, $(filter-out main.f90, $(SRC)))
SRCFULL = $(MODSFULL) $(SRCDIR)/main.f90
OBJ = $(addprefix $(OBJDIR)/, $(SRC:.f90=.o))
DEP = $(OBJ:.o=.d)

# Target name
TARGET = run.exe

# Default rule
all: $(TARGET)

# Linking step
$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $^ $(LFLAGS) -o $@

# ========== Define Per-Directory Compile Rules ==========
define make_compile_rule
$(OBJDIR)/$(1)/%.o: $(MODDIR)/$(1)/%.f90
	$(FC) $(FFLAGS) $(DEPFLAGS) -J$(OBJDIR) -c $$< -o $$@
endef

$(eval $(call make_compile_rule,RNGs))
$(eval $(call make_compile_rule,variables))
$(eval $(call make_compile_rule,parsing))

# Generic rule for modules at the root of src/modules/
$(OBJDIR)/%.o: $(MODDIR)/%.f90
	$(FC) $(FFLAGS) $(DEPFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/main.o: $(SRCDIR)/main.f90
	$(FC) $(FFLAGS) $(DEPFLAGS) -J$(OBJDIR) -c $< -o $@

# Include dependency files (if they exist)
-include $(DEP)

# Clean rule
clean:
	rm -rf $(OBJDIR)/*

.PHONY: all clean