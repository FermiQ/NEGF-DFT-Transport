# .SUFFIXES: Defines the default suffixes for dependency rules.  In this case, it's specifying that .f90 files (Fortran 90 source code) are relevant.

# FC: Specifies the Fortran compiler to use.  Here, it's set to ifort (Intel Fortran Compiler).

# FFLAGS: Compiler flags for the Fortran compiler.
# -g: Enables debugging information in the compiled code.
# -traceback: Enables stack trace information for debugging.
# -c: Compiles the code but does not link it into an executable.  This creates object files (.o).

# LINK: Linker flags.  -L. specifies that the linker should search the current directory for libraries.


# OBJS: Lists the object files that will be linked to create the final executable.

# %.o %.mod : %.f90: This is a pattern rule.  It defines how to create object files (.o) and module files (.mod) from Fortran source files (.f90).
# The `$@` represents the target (the .o and .mod files), and `$*.f90` represents the source file (the .f90 file).
# The rule states that to create `$@` (the .o and .mod files), the compiler `$(FC)` should be invoked with the flags `$(FFLAGS)` on the source file `$*.f90`.

# transomat_make_conquest_input: $(OBJS): This is a rule to create the executable `transomat_make_conquest_input`.
# It specifies that the executable depends on the object files listed in `$(OBJS)`.
# The `$@` represents the target (the executable).
# The rule states that to create `$@` (the executable), the compiler `$(FC)` should be invoked with the object files `$(OBJS)` and the linker flags `$(LINK)`.


# clean: This is a rule to clean up the build directory.
# It removes all object files (.o), module files (.mod), and the executable `transomat_make_conquest_input`.



.SUFFIXES: .f90

FC = ifort

FFLAGS = -g -traceback -c

LINK = -L.




OBJS = make_conquest_mod.o transomat_make_conquest_input.o
%.o %.mod : %.f90
	$(FC) $(FFLAGS)  $*.f90

transomat_make_conquest_input: $(OBJS)
	$(FC)  -o $@ $(OBJS) $(LINK)


clean:
	rm -f *.o *.mod transomat_make_conquest_input

