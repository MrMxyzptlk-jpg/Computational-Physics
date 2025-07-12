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
$(shell mkdir -p $(OBJDIR) \
				$(OBJDIR)/parameters \
				$(OBJDIR)/RNGs \
				$(OBJDIR)/properties \
				$(OBJDIR)/potentials \
				$(OBJDIR)/measurements \
				$(OBJDIR)/thermostats \
				$(OBJDIR)/parsing \
				$(OBJDIR)/forces \
				$(OBJDIR)/integrators \
				$(OBJDIR)/phases)

# List source files in correct dependency order
SRC =  	parameters/precisionMod.f90 \
        parameters/formatsMod.f90 \
        parameters/constantsMod.f90 \
        RNGs/mzranMod.f90 \
        RNGs/mzran_threadsafeMod.f90 \
        RNGs/randomMod.f90 \
        variablesMod.f90 \
        properties/propertiesDeclarationMod.f90 \
        subroutinesMod.f90 \
        potentials/potentialContributionMod.f90 \
        potentials/Lennard_JonesMod.f90 \
        potentials/Coulomb_EwaldMod.f90 \
        potentials/Coulomb_ReactionFieldMod.f90 \
        potentials/potentialsMod.f90 \
        measurements/meanSquareDisplacementMod.f90 \
        measurements/pairCorrelationFunctionMod.f90 \
        measurements/structureFactorMod.f90 \
        measurements/observablesMod.f90 \
        thermostats/MD_rescaleMod.f90 \
        thermostats/MD_BerendsenMod.f90 \
        thermostats/MC_RndStepMod.f90 \
        thermostats/thermostatsMod.f90 \
        properties/updatePositions_velVerletMod.f90 \
        properties/updatePositions_MCMod.f90 \
        properties/updatePositions_BDMod.f90 \
        properties/updateVelocities_velVerletMod.f90 \
        properties/updatePositionsMod.f90 \
		parsing/checkParsingMod.f90 \
		parsing/get_parsed_valueMod.f90 \
        parsing/parsingMod.f90 \
        forces/forces_contributionMod.f90 \
        forces/forces_EwaldMod.f90 \
        forces/forces_LinkedListsMod.f90 \
        forces/forces_AllvsAllMod.f90 \
        forces/forcesMod.f90 \
        writing2filesMod.f90 \
        integrators/velVerlet_stepMod.f90 \
        integrators/MC_stepMod.f90 \
        integrators/BD_stepMod.f90 \
        integrators/integratorsMod.f90 \
        initializationsMod.f90 \
        measurements/measurementsMod.f90 \
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
$(eval $(call make_compile_rule,parameters))
$(eval $(call make_compile_rule,properties))
$(eval $(call make_compile_rule,potentials))
$(eval $(call make_compile_rule,measurements))
$(eval $(call make_compile_rule,thermostats))
$(eval $(call make_compile_rule,parsing))
$(eval $(call make_compile_rule,forces))
$(eval $(call make_compile_rule,integrators))
$(eval $(call make_compile_rule,phases))

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