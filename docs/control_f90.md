# `control.f90` - Control Parameter Management

## Overview

This module is central to TransOMat's configuration. It defines the structures for holding all control parameters that steer the simulation, provides mechanisms to initialize their expected key names (as found in the `transp.ini` input file), and implements the logic to parse these parameters from the input file. It essentially acts as the bridge between the user's input specifications and the program's internal variables.

## Key Components

*   **`module control`**: Encapsulates all definitions and subroutines related to control parameters.
    *   **Derived Types for Parameters:**
        *   `type ikey`: Defines a structure for integer parameters, holding its name (`kname`) and default value (`def`).
        *   `type rkey`: Defines a structure for real (double precision) parameters.
        *   `type lkey`: Defines a structure for logical parameters.
        *   `type skey`: Defines a structure for string parameters.
    *   **`type cntrl_`**: A large derived type that aggregates all known control parameters for a TransOMat run. Each field in `cntrl_` is one of the `ikey`, `rkey`, `lkey`, or `skey` types. Examples include `ecc_dir` (string), `estart` (real), `n_energy_steps` (integer), `oneshot` (logical).
    *   **`interface get_cntrl_key`**: Defines a generic interface for type-specific getter subroutines.
        *   `get_skey(keyname, keystr, l_iscritical)`: Reads a string value.
        *   `get_rkey(keyname, kvalue, l_iscritical)`: Reads a real value.
        *   `get_ikey(keyname, kvalue, l_iscritical)`: Reads an integer value.
        *   `get_lkey(keyname, kvalue, l_iscritical)`: Reads a logical value.
    *   **`subroutine init_cntrl(cntrl)`**:
        *   Populates the `kname` field of each parameter within the `cntrl_` structure. This sets up the expected string (e.g., `$ecc_dir=`) that the parser will look for in the input file.
    *   **`subroutine get_key(keyname, keystr, keyint, keyreal, keylogic, l_iscritical)`**:
        *   The core parsing routine. It reads the control input file (identified by `iunit_control` from `globals.f90`), line by line.
        *   Searches for `keyname` in each line.
        *   If found, extracts the value part of the string, converts it to the appropriate type (string, integer, real, or logical), and stores it.
        *   Prints the found key-value pair using `PetscPrintf`.
        *   If a `l_iscritical` key is not found, the program stops. Otherwise, it prints the default or current value.

*   **Standalone Getter Subroutines (outside the module but using it):**
    *   `subroutine get_skey(keyname, kvalue, l_iscritical)`: Wrapper around `get_key` for string parameters.
    *   `subroutine get_rkey(keyname, kvalue, l_iscritical)`: Wrapper for real parameters.
    *   `subroutine get_ikey(keyname, kvalue, l_iscritical)`: Wrapper for integer parameters.
    *   `subroutine get_lkey(keyname, kvalue, l_iscritical)`: Wrapper for logical parameters.

## Important Variables/Constants

*   **`cntrl_` type instance (typically named `cntrl` in `init_control.f90`):** The main structure holding all parameters. The actual global instance of these parameters is usually stored in `globals.f90` and populated by `init_control.f90`.
*   **`keyname` (Character(strln)) in `get_key` and wrappers:** The string identifier for a parameter as expected in the input file (e.g., `$ecc_dir=`).
*   **`l_iscritical` (Logical) in `get_key` and wrappers:** A flag indicating whether a parameter is mandatory. If true and the parameter is not found, the program will terminate.

## Usage Examples

This module is primarily used internally by `init_control.f90`. The workflow is:
1.  An instance of `cntrl_` is declared (usually in `globals.f90`).
2.  `init_cntrl()` (from this module, but typically called by the top-level `init_control()` subroutine in `init_control.f90`) is called to set the `kname` for each parameter.
3.  The type-specific getter subroutines (`get_skey`, `get_rkey`, etc.) are then called (e.g., in `init_control.f90`) for each parameter to parse its value from the input file.

```fortran
! Conceptual usage within init_control.f90
! Assuming 'gcntrl' is a global variable of type 'cntrl_' from globals.f90
! and 'l_optional' is a logical .false.

! First, define the key names (done by this module's init_cntrl)
! call control_module_init_cntrl(gcntrl) ! (Illustrative name)

! Then, read each key
call get_skey(gcntrl%ecc_dir%kname, ecc_dir_val) ! ecc_dir_val is now populated
call get_rkey(gcntrl%estart%kname, estart_val)
call get_ikey(gcntrl%n_energy_steps%kname, n_steps_val)
call get_lkey(gcntrl%oneshot%kname, oneshot_val, l_optional)
! ... and so on for all parameters
```

The input file (`transp.ini`) would look like:
```ini
$ecc_dir=./central_region/
$e_start=-1.0
$n_steps=100
$oneshot=.FALSE.
# ... other parameters
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `USE kinds`: For `dp` (double precision) and `strln` (string length) definitions.
    *   `USE misc`: For `str2` (string to type conversion) and `New_line` utilities.
    *   `USE globals`: For `inode` (MPI rank) and `iunit_control` (file unit for `transp.ini`).
*   **External Library Dependencies:**
    *   `USE petsc, only : PETSC_COMM_WORLD`: Used by `get_key` for `PetscPrintf` to ensure controlled output in parallel environments.
*   **Interactions:**
    *   This module is primarily called by `init_control.f90` to define and parse all runtime parameters.
    *   It reads from the input file specified by the unit `iunit_control` (which is opened in `init_control.f90`).
    *   The parsed values are then stored in global variables (typically in `globals.f90`) for use throughout the TransOMat application.

---
*This documentation was auto-generated by Jules.*
