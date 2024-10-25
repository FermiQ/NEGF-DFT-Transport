# Compile the Fortran code fc_to_petsc.f90 using the MPIF90 compiler with various options.

# mpif90: The MPI Fortran compiler.
# -mkl=cluster: Use the Intel Math Kernel Library (MKL) for cluster computing.  This optimizes math routines for parallel processing.
# -fpp: Enable preprocessor directives in the Fortran code (e.g., #ifdef, #include).
# -O2: Optimize the compiled code for speed (level 2 optimization).
# -traceback: Enable stack trace generation for debugging purposes.  Provides more information when errors occur.
# -g: Generate debugging information.  Allows using debuggers like gdb to step through the code.
# -xHost: Optimize the code for the specific host architecture.  Improves performance by tailoring to the CPU.
# -I${PETSC_DIR}/include/petsc: Include directory for PETSc header files.  PETSc (Portable, Extensible Toolkit for Scientific Computation) is a suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations.
# -I${PETSC_DIR}/include/petsc/finclude: Include directory for PETSc Fortran header files.
# -I${PETSC_DIR}/include: Include directory for general PETSc header files.
# fc_to_petsc.f90: The Fortran source code file to be compiled.
# -o fc_to_petsc: Specify the output executable file name.
# -L$PETSC_DIR/lib: Link against libraries located in the PETSc library directory.
# -lpetsc: Link against the PETSc library.
# -lslepc: Link against the SLEPc library (Scalable Library for Eigenvalue Problem Computations), which is often used with PETSc.

mpif90 -mkl=cluster -fpp -O2 -traceback -g -xHost -I${PETSC_DIR}/include/petsc  -I${PETSC_DIR}/include/petsc/finclude -I${PETSC_DIR}/include  fc_to_petsc.f90 -o fc_to_petsc -L$PETSC_DIR/lib -lpetsc -lslepc
