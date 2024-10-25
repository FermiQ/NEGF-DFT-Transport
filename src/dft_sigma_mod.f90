      !! DFT+sigma scf
!if (dftsigma) then
!  ndfts=jmu_dftsigma-imu_dftsigma+1  ! Calculate the number of DFT+sigma orbitals.
!  if (inode.eq.0) write(6,*) "imu_dftsigma,jmu_dftsigma",imu_dftsigma,jmu_dftsigma  ! Print indices if on root node.
!  call petsc_split_matrix(p_h00_cc(0,0),p_tmp1,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)  ! Split the Hamiltonian matrix.
!  call petsc_split_matrix(p_s00_cc(0,0),p_tmp2,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)  ! Split the overlap matrix.
!  call petsc_split_matrix(p_k00_cc(0,0),p_tmp3,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)  ! Split the Coulomb matrix.
!  call MatGetType(p_tmp1,sdummy,ierr)  ! Get the type of the matrix p_tmp1.
!  call MatGetInfo(p_tmp1, MAT_LOCAL, info,p_ierr)  ! Get information about the matrix p_tmp1.
!  call MatGetOwnershipRange(p_tmp1,nl1,nl2,ierr)  ! Get the ownership range of the matrix p_tmp1.
!  write(sdummy,fmt='(A,5i8)') "p_tmp1 "//trim(sdummy),inode,int(info(mat_info_nz_allocated)),int(info(mat_info_nz_used)),nl1,nl2  ! Write matrix information.
!  call petsc_print(sdummy)  ! Print the string sdummy.
!
!  call MatGetOwnershipRange(p_tmp1,nl1,nl2,ierr)  ! Get the ownership range of the matrix p_tmp1.
!  call MatCreateDense(PETSC_COMM_WORLD,nl2-nl1,PETSC_DECIDE,ndfts,ndfts,PETSC_NULL_SCALAR ,p_ev,ierr)  ! Create a dense matrix for eigenvalues.
!
!
!
!  call MatCreateVecs(p_ev,p_ew1,PETSC_NULL_vec,ierr)  ! Create vectors associated with the matrix p_ev.
!  call MatCreateVecs(p_ev,p_ew2,PETSC_NULL_vec,ierr)  ! Create vectors associated with the matrix p_ev.
!
!~    !     call diag_mat(p_tmp1,p_tmp2,p_ew1,p_ev,ndfts,ierr)  ! Commented out diagonalization routine.
!
!
!  p_scal=d_virt  ! Assign value to p_scal.
!  call MatAXPY(p_tmp1,p_scal,p_tmp2, DIFFERENT_NONZERO_PATTERN,ierr)  ! h=h+d_virt*s  Matrix addition.
!
!
!
!  call MatMatMatMult(p_tmp2,p_tmp3,p_tmp2,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,p_tmp4,ierr)  ! s*k*s Matrix multiplication.
!
!
!
!  p_scal=d_occ-d_virt  ! Assign value to p_scal.
!
!  call MatAXPY(p_tmp1,p_scal,p_tmp4,DIFFERENT_NONZERO_PATTERN,ierr)  ! h=h+(d_occ-d_virt)*s*k*s Matrix addition.
!
!
!
!~    !     call diag_mat(p_tmp1,p_tmp2,p_ew2,p_ev,ndfts,ierr)  ! Commented out diagonalization routine.
!
!  call MatDestroy(p_tmp4,ierr)  ! Destroy the matrix p_tmp4.
!
!  call petsc_split_matrix(p_h00_cc(0,0),p_tmp4,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)  ! Split the Hamiltonian matrix.
!
!
!
!
!  allocate(cols(ndfts),p_vals(ndfts),stat=ierr)  ! Allocate memory for cols and p_vals arrays.
!  if (ierr.ne.0) then
!    write(0,*) "allocation error init_sys ",ierr  ! Print allocation error message.
!    stop  ! Stop execution.
!  end if
!
!  call MatGetOwnershipRange(p_tmp4,nl1,nl2,ierr)  ! Get the ownership range of the matrix p_tmp4.
!
!  do imu1=nl1,nl2-1  ! Loop over rows.
!    call MatGetRow(p_tmp4,imu1,ncols,cols,p_vals,p_ierr)  ! Get a row from the matrix p_tmp4.
!    do imu2=1,ncols  ! Loop over columns in the row.
!      idxm=imu1  ! Assign row index.
!      idym=cols(imu2)  ! Assign column index.
!      call MatGetValues(p_tmp1,1,idxm,1,idym,p_vals(1),ierr)  ! Get matrix element from p_tmp1.
!      p_scal=p_vals(1)  ! Assign matrix element to p_scal.
!      idxm=imu1+imu_dftsigma-1  ! Adjust row index.
!      idym=cols(imu2)+imu_dftsigma-1  ! Adjust column index.
!      call MatSetValue(p_h00_cc(0,0),idxm,idym,p_scal,INSERT_VALUES,ierr)  ! Set matrix element in p_h00_cc.
!    end do
!    call MatRestoreRow(p_tmp4,imu1,ncols,cols,p_vals,p_ierr)  ! Restore the row.
!  end do
!
!  call MatAssemblyBegin(p_h00_cc(0,0),MAT_FINAL_ASSEMBLY,ierr)  ! Begin matrix assembly.
!  call MatAssemblyEnd(p_h00_cc(0,0),MAT_FINAL_ASSEMBLY,ierr)  ! End matrix assembly.
!
!
!  call MatDestroy(p_tmp1,ierr)  ! Destroy the matrix p_tmp1.
!  call MatDestroy(p_tmp2,ierr)  ! Destroy the matrix p_tmp2.
!
!  call petsc_split_matrix(p_h00_cc(0,0),p_tmp1,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)  ! Split the Hamiltonian matrix.
!  call petsc_split_matrix(p_s00_cc(0,0),p_tmp2,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)  ! Split the overlap matrix.
!~    !     call diag_mat(p_tmp1,p_tmp2,p_ew1,p_ev,ndfts,ierr)  ! Commented out diagonalization routine.
!
!  call MatDestroy(p_tmp3,ierr)  ! Destroy the matrix p_tmp3.
!  call MatDestroy(p_tmp4,ierr)  ! Destroy the matrix p_tmp4.
!  call MatDestroy(p_ev,ierr)  ! Destroy the matrix p_ev.
!  call VecDestroy(p_ew1,ierr)  ! Destroy the vector p_ew1.
!  call VecDestroy(p_ew2,ierr)  ! Destroy the vector p_ew2.
!
!  deallocate(cols)  ! Deallocate memory for cols array.
!
!end if