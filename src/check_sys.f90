! Check if the left/right electrode can be translated onto the l and r ends of the ecc.
! At the time being, no reordering, rotations, and different dimensions between l/r_ecc and electrodes are allowed.

subroutine check_sys(xyz_ecc, xyz_elec_l, xyz_elec_r, nat_ecc, nat_elec_l, nat_elec_r, lcr_info)
  ! Use statements for modules
  use kinds ! Module defining kind type parameters (e.g., dp for double precision)
  use petsc_mod, only: pstr_out, petsc_print_master ! Module for PETSc output functions
  use globals, only: eps_geo, inode ! Module containing global variables, eps_geo (geometric tolerance), inode (likely node index)
  use error_handler ! Module for error handling

  ! Implicit none statement to enforce explicit declaration of variables
  implicit none

  ! Variable declarations
  integer :: nat_ecc, nat_elec_l, nat_elec_r ! Number of atoms in the central region (ecc), left electrode, and right electrode, respectively.
  real(dp) :: xyz_ecc(3, nat_ecc), xyz_elec_l(3, nat_elec_l), xyz_elec_r(3, nat_elec_r) ! Cartesian coordinates of atoms in the ecc, left electrode, and right electrode.  (3: x,y,z; nat_*: number of atoms)
  character(1) :: lcr_info(nat_ecc) ! Character array to store the information of each atom in the ecc (l: left, r: right, c: central).

  ! Local variables
  real(dp) :: transvec(3) ! Translation vector
  integer :: i ! Loop counter

  ! Initialize lcr_info array to 'c' (central)
  lcr_info = "c"

  ! Calculate translation vector from the first atom of the ecc and left electrode
  transvec = xyz_ecc(1:3, 1) - xyz_elec_l(1:3, 1)

  ! Loop through each atom in the left electrode
  do i = 1, nat_elec_l
    ! Check if the coordinates match within the tolerance eps_geo after translation
    if (.not. all(abs(xyz_ecc(1:3, i) - xyz_elec_l(1:3, i) - transvec) .le. eps_geo)) then
      ! Print error messages if coordinates do not match
      write (pstr_out, fmt='(A)') "missmatch between ecc_l and left electrode"; call petsc_print_master()
      write (pstr_out, fmt='(A,3e24.12)') "transvec ", transvec; call petsc_print_master()
      write (pstr_out, fmt='(A,i8)') "ecc_l,elec_l,abs(diff)", i; call petsc_print_master()
      write (pstr_out, fmt='(3e24.12)') xyz_ecc(1:3,  i); call petsc_print_master()
      write (pstr_out, fmt='(3e24.12)') xyz_elec_l(1:3, i); call petsc_print_master()
      write (pstr_out, fmt='(3e24.12)') transvec; call petsc_print_master()
      write (pstr_out, fmt='(3e24.12)') xyz_ecc(1:3,  i) - xyz_elec_l(1:3, i) - transvec; call petsc_print_master()
      errormsg = "lattice error" ! Set error message
      call error() ! Call error handling routine
    end if
    ! Mark the corresponding atom in the ecc as belonging to the left electrode
    lcr_info(i) = "l"
  end do

  ! Calculate translation vector from the last atom of the ecc and right electrode
  transvec = xyz_ecc(1:3, nat_ecc) - xyz_elec_r(1:3, nat_elec_r)

  ! Loop through each atom in the right electrode
  do i = 1, nat_elec_r
    ! Check if the coordinates match within the tolerance eps_geo after translation
    if (.not. all(abs(xyz_ecc(1:3, nat_ecc - nat_elec_r + i) - xyz_elec_r(1:3, i) - transvec) .le. eps_geo)) then
      ! Print error messages if coordinates do not match
      write (pstr_out, fmt='(A)') "missmatch between ecc_r and right electrode"; call petsc_print_master()
      write (pstr_out, fmt='(A,3e24.12)') "transvec ", transvec; call petsc_print_master()
      write (pstr_out, fmt='(A,i8)') "ecc_r,elec_r,abs(diff)", i ; call petsc_print_master()
      write (pstr_out, fmt='(3e24.12)') xyz_ecc(1:3, nat_ecc - nat_elec_r + i); call petsc_print_master()
      write (pstr_out, fmt='(3e24.12)') xyz_elec_r(1:3, i); call petsc_print_master()
      write (pstr_out, fmt='(3e24.12)') transvec; call petsc_print_master()
      write (pstr_out, fmt='(3e24.12)') xyz_ecc(1:3, nat_ecc - nat_elec_r + i) - xyz_elec_r(1:3, i) - transvec; call petsc_print_master()
      errormsg = "lattice error" ! Set error message
      call error() ! Call error handling routine
    end if
    ! Mark the corresponding atom in the ecc as belonging to the right electrode
    lcr_info(nat_ecc - nat_elec_r + i) = "r"
  end do

end subroutine check_sys