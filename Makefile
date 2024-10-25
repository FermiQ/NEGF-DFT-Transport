# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $

# Include architecture-specific Makefile
include Makefile.arch

# Default rule for targets that don't exist
.DEFAULT:
	-echo $@ does not exist.

# Main target: build the executable transomat
all: $(BIN_DIR)/transomat

# Rule to compile analysis_mod.f90
# Dependencies: src/analysis_mod.f90, object files for PETSc module, globals, and k_on_demand module
# Prerequisites: object directory and module directory
$(OBJ_DIR)/analysis_mod.o: src/analysis_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/k_on_demand_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/analysis_mod.f90 -o $(OBJ_DIR)/analysis_mod.o

# Rule to compile check_sys.f90
# Dependencies: src/check_sys.f90, object files for kinds, PETSc module, globals, and error_handler module
# Prerequisites: object directory and module directory
$(OBJ_DIR)/check_sys.o: src/check_sys.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/check_sys.f90 -o $(OBJ_DIR)/check_sys.o

# Rule to compile conquest_mod.f90
# Dependencies: src/conquest_mod.f90, object files for kinds, error_handler module, and globals
# Prerequisites: object directory and module directory
$(OBJ_DIR)/conquest_mod.o: src/conquest_mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/conquest_mod.f90 -o $(OBJ_DIR)/conquest_mod.o

# Rule to compile control.f90
# Dependencies: src/control.f90, object files for kinds, misc, and globals
# Prerequisites: object directory and module directory
$(OBJ_DIR)/control.o: src/control.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/control.f90 -o $(OBJ_DIR)/control.o

# Rule to compile dft_sigma_mod.f90
# Prerequisites: object directory and module directory
$(OBJ_DIR)/dft_sigma_mod.o: src/dft_sigma_mod.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/dft_sigma_mod.f90 -o $(OBJ_DIR)/dft_sigma_mod.o

# Rule to compile error_handler_mod.f90
# Dependencies: src/error_handler_mod.f90, object file for kinds
# Prerequisites: object directory and module directory
$(OBJ_DIR)/error_handler_mod.o: src/error_handler_mod.f90 $(OBJ_DIR)/kinds.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/error_handler_mod.f90 -o $(OBJ_DIR)/error_handler_mod.o

# Rule to compile ft_mod.f90
# Dependencies: src/ft_mod.f90, object files for PETSc module, kinds, and globals
# Prerequisites: object directory and module directory
$(OBJ_DIR)/ft_mod.o: src/ft_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/ft_mod.f90 -o $(OBJ_DIR)/ft_mod.o

# Rule to compile globals.f90
# Dependencies: src/globals.f90, object file for kinds
# Prerequisites: object directory and module directory
$(OBJ_DIR)/globals.o: src/globals.f90 $(OBJ_DIR)/kinds.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/globals.f90 -o $(OBJ_DIR)/globals.o

# Rule to compile init_control.f90
# Dependencies: src/init_control.f90, object files for globals and control
# Prerequisites: object directory and module directory
$(OBJ_DIR)/init_control.o: src/init_control.f90 $(OBJ_DIR)/globals.o $(OBJ_DIR)/control.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/init_control.f90 -o $(OBJ_DIR)/init_control.o

# Rule to compile init_ep_coupling_mod.f90
# Dependencies: src/init_ep_coupling_mod.f90, object files for kinds, misc, globals, PETSc module, petsc_wrapper module, and error_handler module
# Prerequisites: object directory and module directory
$(OBJ_DIR)/init_ep_coupling_mod.o: src/init_ep_coupling_mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/init_ep_coupling_mod.f90 -o $(OBJ_DIR)/init_ep_coupling_mod.o

# Rule to compile init_petsc.f90
# Dependencies: src/init_petsc.f90, object file for globals
# Prerequisites: object directory and module directory
$(OBJ_DIR)/init_petsc.o: src/init_petsc.f90 $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/init_petsc.f90 -o $(OBJ_DIR)/init_petsc.o

# Rule to compile init_sys.f90
# Dependencies: src/init_sys.f90, many object files for various modules
# Prerequisites: object directory and module directory
$(OBJ_DIR)/init_sys.o: src/init_sys.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/read_sysinfo.o $(OBJ_DIR)/lattice_mod.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/conquest_mod.o $(OBJ_DIR)/kmat_mod.o $(OBJ_DIR)/slepc_mod.o $(OBJ_DIR)/init_ep_coupling_mod.o $(OBJ_DIR)/phonon_mod.o $(OBJ_DIR)/loe_mod.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/k_on_demand_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/init_sys.f90 -o $(OBJ_DIR)/init_sys.o

# Rule to compile integration_weights.f90
# Prerequisites: object directory and module directory
$(OBJ_DIR)/integration_weights.o: src/integration_weights.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/integration_weights.f90 -o $(OBJ_DIR)/integration_weights.o

# Rule to compile integrator.f90
# Dependencies: src/integrator.f90, many object files for various modules
# Prerequisites: object directory and module directory
$(OBJ_DIR)/integrator.o: src/integrator.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/math_stuff.o $(OBJ_DIR)/surf_gf.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/integration_weights.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/pexsi_wrapper.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/integrator.f90 -o $(OBJ_DIR)/integrator.o

# Rule to compile integrator_scalar.f90
# Dependencies: src/integrator_scalar.f90, object files for kinds, PETSc module, math_stuff, globals, and integrator
# Prerequisites: object directory and module directory
$(OBJ_DIR)/integrator_scalar.o: src/integrator_scalar.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/math_stuff.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/integrator.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/integrator_scalar.f90 -o $(OBJ_DIR)/integrator_scalar.o

# Rule to compile j_dump_mod.f90
# Dependencies: src/j_dump_mod.f90, object files for PETSc module, petsc_wrapper module, kinds, globals, and write_conquest_dump module
# Prerequisites: object directory and module directory
$(OBJ_DIR)/j_dump_mod.o: src/j_dump_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/write_conquest_dump.mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/j_dump_mod.f90 -o $(OBJ_DIR)/j_dump_mod.o

# Rule to compile k_on_demand_mod.f90
# Dependencies: src/k_on_demand_mod.f90, object files for globals, ft_mod, and PETSc module
# Prerequisites: object directory and module directory
$(OBJ_DIR)/k_on_demand_mod.o: src/k_on_demand_mod.f90 $(OBJ_DIR)/globals.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/petsc_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/k_on_demand_mod.f90 -o $(OBJ_DIR)/k_on_demand_mod.o

# Rule to compile kinds.f90
# Prerequisites: object directory and module directory
$(OBJ_DIR)/kinds.o: src/kinds.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/kinds.f90 -o $(OBJ_DIR)/kinds.o

# Rule to compile kmat_from_diag.f90
# Dependencies: src/kmat_from_diag.f90, many object files for various modules
# Prerequisites: object directory and module directory
$(OBJ_DIR)/kmat_from_diag.o: src/kmat_from_diag.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/k_on_demand_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/slepc_mod.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/kmat_from_diag.f90 -o $(OBJ_DIR)/kmat_from_diag.o

# Rule to compile kmat_mod.f90
# Dependencies: src/kmat_mod.f90, object files for PETSc module, petsc_wrapper module, kinds, globals, and write_conquest_dump module
# Prerequisites: object directory and module directory
$(OBJ_DIR)/kmat_mod.o: src/kmat_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/write_conquest_dump.mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/kmat_mod.f90 -o $(OBJ_DIR)/kmat_mod.o

# Rule to compile lattice_mod.f90
# Dependencies: src/lattice_mod.f90, object files for kinds, globals, error_handler module, and math_stuff
# Prerequisites: object directory and module directory
$(OBJ_DIR)/lattice_mod.o: src/lattice_mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/math_stuff.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/lattice_mod.f90 -o $(OBJ_DIR)/lattice_mod.o

# Rule to compile loe_mod.f90
# Dependencies: src/loe_mod.f90, many object files for various modules
# Prerequisites: object directory and module directory
$(OBJ_DIR)/loe_mod.o: src/loe_mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/phonon_mod.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/integrator_scalar.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/loe_mod.f90 -o $(OBJ_DIR)/loe_mod.o

# Rule to compile math_stuff.f90
# Dependencies: src/math_stuff.f90, object file for kinds
# Prerequisites: object directory and module directory
$(OBJ_DIR)/math_stuff.o: src/math_stuff.f90 $(OBJ_DIR)/kinds.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/math_stuff.f90 -o $(OBJ_DIR)/math_stuff.o

# Rule to compile misc.f90
# Dependencies: src/misc.f90, object file for kinds
# Prerequisites: object directory and module directory
$(OBJ_DIR)/misc.o: src/misc.f90 $(OBJ_DIR)/kinds.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/misc.f90 -o $(OBJ_DIR)/misc.o

# Rule to compile petsc_mod.f90
# Dependencies: src/petsc_mod.f90, object files for kinds, globals, error_handler module, and math_stuff
# Prerequisites: object directory and module directory
$(OBJ_DIR)/petsc_mod.o: src/petsc_mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/math_stuff.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/petsc_mod.f90 -o $(OBJ_DIR)/petsc_mod.o

# Rule to compile petsc_test.f90
# Dependencies: src/petsc_test.f90, object files for slepc_mod, PETSc module, petsc_wrapper module, and globals
# Prerequisites: object directory and module directory
$(OBJ_DIR)/petsc_test.o: src/petsc_test.f90 $(OBJ_DIR)/slepc_mod.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/petsc_test.f90 -o $(OBJ_DIR)/petsc_test.o

# Rule to compile petsc_wrapper_mod.f90
# Dependencies: src/petsc_wrapper_mod.f90, many object files for various modules
# Prerequisites: object directory and module directory
$(OBJ_DIR)/petsc_wrapper_mod.o: src/petsc_wrapper_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/timing_mod.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/slepc_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/petsc_wrapper_mod.f90 -o $(OBJ_DIR)/petsc_wrapper_mod.o

# Rule to compile pexsi_wrapper.f90
# Dependencies: src/pexsi_wrapper.f90, object files for kinds, PETSc module, error_handler module, and globals
# Prerequisites: object directory and module directory
$(OBJ_DIR)/pexsi_wrapper.o: src/pexsi_wrapper.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/pexsi_wrapper.f90 -o $(OBJ_DIR)/pexsi_wrapper.o

# Rule to compile phonon_mod.f90
# Dependencies: src/phonon_mod.f90, many object files for various modules
# Prerequisites: object directory and module directory
$(OBJ_DIR)/phonon_mod.o: src/phonon_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/slepc_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/phonon_mod.f90 -o $(OBJ_DIR)/phonon_mod.o

# Rule to compile read_sysinfo.f90
# Dependencies: src/read_sysinfo.f90, object files for kinds, lattice_mod, globals, and PETSc module
# Prerequisites: object directory and module directory
$(OBJ_DIR)/read_sysinfo.o: src/read_sysinfo.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/lattice_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/read_sysinfo.f90 -o $(OBJ_DIR)/read_sysinfo.o

# Rule to compile scf_mod.f90
# Dependencies: src/scf_mod.f90, many object files for various modules
# Prerequisites: object directory and module directory
$(OBJ_DIR)/scf_mod.o: src/scf_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/kmat_mod.o $(OBJ_DIR)/timing_mod.o $(OBJ_DIR)/k_on_demand_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/scf_mod.f90 -o $(OBJ_DIR)/scf_mod.o

# Rule to compile slepc_mod.f90
# Dependencies: src/slepc_mod.f90, object files for PETSc module, globals, kinds, and error_handler module
# Prerequisites: object directory and module directory
$(OBJ_DIR)/slepc_mod.o: src/slepc_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/slepc_mod.f90 -o $(OBJ_DIR)/slepc_mod.o

# Rule to compile surf_gf.f90
# Dependencies: src/surf_gf.f90, many object files for various modules
# Prerequisites: object directory and module directory
$(OBJ_DIR)/surf_gf.o: src/surf_gf.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/math_stuff.o $(OBJ_DIR)/timing_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/surf_gf.f90 -o $(OBJ_DIR)/surf_gf.o

# Rule to compile timing_mod.f90
# Dependencies: src/timing_mod.f90, object file for kinds
# Prerequisites: object directory and module directory
$(OBJ_DIR)/timing_mod.o: src/timing_mod.f90 $(OBJ_DIR)/kinds.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/timing_mod.f90 -o $(OBJ_DIR)/timing_mod.o

# Rule to compile transmission_mod.f90
# Dependencies: src/transmission_mod.f90, many object files for various modules
# Prerequisites: object directory and module directory
$(OBJ_DIR)/transmission_mod.o: src/transmission_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/slepc_mod.o $(OBJ_DIR)/k_on_demand_mod.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/transmission_mod.f90 -o $(OBJ_DIR)/transmission_mod.o

# Rule to compile transomat.f90
# Dependencies: src/transomat.f90, many object files for various modules
# Prerequisites: object directory and module directory
$(OBJ_DIR)/transomat.o: src/transomat.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/init_sys.o $(OBJ_DIR)/init_petsc.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/scf_mod.o $(OBJ_DIR)/loe_mod.o $(OBJ_DIR)/phonon_mod.o $(OBJ_DIR)/transmission_mod.o $(OBJ_DIR)/kmat_from_diag.o $(OBJ_DIR)/kmat_mod.o $(OBJ_DIR)/analysis_mod.o $(OBJ_DIR)/timing_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/transomat.f90 -o $(OBJ_DIR)/transomat.o

# Rule to compile write_conquest_dump.mod.f90
# Dependencies: src/write_conquest_dump.mod.f90, object files for kinds, globals, PETSc module, and error_handler module
# Prerequisites: object directory and module directory
$(OBJ_DIR)/write_conquest_dump.mod.o: src/write_conquest_dump.mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/write_conquest_dump.mod.f90 -o $(OBJ_DIR)/write_conquest_dump.mod.o

# List of source files
SRC = src/check_sys.f90 src/lattice_mod.f90 src/k_on_demand_mod.f90 src/j_dump_mod.f90 src/transmission_mod.f90 src/ft_mod.f90 src/control.f90 src/dft_sigma_mod.f90 src/kinds.f90 src/surf_gf.f90 src/kmat_mod.f90 src/analysis_mod.f90 src/integrator_scalar.f90 src/init_petsc.f90 src/kmat_from_diag.f90 src/pexsi_wrapper.f90 src/loe_mod.f90 src/petsc_mod.f90 src/petsc_test.f90 src/transomat.f90 src/init_control.f90 src/scf_mod.f90 src/conquest_mod.f90 src/init_ep_coupling_mod.f90 src/phonon_mod.f90 src/write_conquest_dump.mod.f90 src/petsc_wrapper_mod.f90 src/integrator.f90 src/read_sysinfo.f90 src/integration_weights.f90 src/globals.f90 src/math_stuff.f90 src/error_handler_mod.f90 src/misc.f90 src/init_sys.f90 src/timing_mod.f90 src/slepc_mod.f90

# List of object files
OBJ = $(OBJ_DIR)/check_sys.o $(OBJ_DIR)/lattice_mod.o $(OBJ_DIR)/k_on_demand_mod.o $(OBJ_DIR)/j_dump_mod.o $(OBJ_DIR)/transmission_mod.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/control.o $(OBJ_DIR)/dft_sigma_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/surf_gf.o $(OBJ_DIR)/kmat_mod.o $(OBJ_DIR)/analysis_mod.o $(OBJ_DIR)/integrator_scalar.o $(OBJ_DIR)/init_petsc.o $(OBJ_DIR)/kmat_from_diag.o $(OBJ_DIR)/pexsi_wrapper.o $(OBJ_DIR)/loe_mod.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_test.o $(OBJ_DIR)/transomat.o $(OBJ_DIR)/init_control.o $(OBJ_DIR)/scf_mod.o $(OBJ_DIR)/conquest_mod.o $(OBJ_DIR)/init_ep_coupling_mod.o $(OBJ_DIR)/phonon_mod.o $(OBJ_DIR)/write_conquest_dump.mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/read_sysinfo.o $(OBJ_DIR)/integration_weights.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/math_stuff.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/init_sys.o $(OBJ_DIR)/timing_mod.o $(OBJ_DIR)/slepc_mod.o

# List of Fortran source files (for localization)
OFF = src/pexsi_wrapper.f90 src/misc.f90 src/check_sys.f90 src/control.f90 src/init_control.f90 src/integration_weights.f90 src/integrator.f90 src/surf_gf.f90 src/petsc_wrapper_mod.f90 src/error_handler_mod.f90 src/integrator_scalar.f90 src/j_dump_mod.f90 src/transomat.f90 src/petsc_test.f90 src/petsc_mod.f90 src/analysis_mod.f90 src/init_ep_coupling_mod.f90 src/timing_mod.f90 src/kmat_mod.f90 src/globals.f90 src/dft_sigma_mod.f90 src/k_on_demand_mod.f90 src/loe_mod.f90 src/init_sys.f90 src/read_sysinfo.f90 src/lattice_mod.f90 src/scf_mod.f90 src/kmat_from_diag.f90 src/phonon_mod.f90 src/conquest_mod.f90 src/kinds.f90 src/slepc_mod.f90 src/math_stuff.f90 src/init_petsc.f90 src/ft_mod.f90 src/write_conquest_dump.mod.f90 src/transmission_mod.f90

# Clean targets
clean: neat
	-rm -f .transomat.cppdefs $(OBJ) $(MOD_DIR)/*.mod $(BIN_DIR)/transomat
neat:
	-rm -f $(TMPFILES)

# Localize source files
localize: $(OFF)
	cp $(OFF) .

# Generate TAGS file
TAGS: $(SRC)
	etags $(SRC)

# Generate tags file
tags: $(SRC)
	ctags $(SRC)

# Create object directory if it doesn't exist
$(OBJ_DIR):
	mkdir $(OBJ_DIR)

# Create module directory if it doesn't exist
$(MOD_DIR):
	mkdir $(MOD_DIR)

# Create binary directory if it doesn't exist
$(BIN_DIR):
	mkdir $(BIN_DIR)

# Link object files to create the executable
# Prerequisites: object files, binary directory, and object directory
$(BIN_DIR)/transomat: $(OBJ) | $(BIN_DIR) $(OBJ_DIR) 
	$(LD) $(OBJ) -o $(BIN_DIR)/transomat  $(LDFLAGS)