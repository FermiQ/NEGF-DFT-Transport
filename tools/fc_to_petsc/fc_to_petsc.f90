module fctopetsc_mod
  implicit none
  
  contains
  
      subroutine diag_mat(a,b,ew,ev,nev,ierr)
#include <slepc/finclude/slepceps.h>
      use slepceps    
      
      implicit none
      
      Mat :: a,b,ev,aa,bb,cc,dd
      Vec :: ew
      PetscInt :: nev,i,nl1,nl2,ml1,ml2,j,maxits
      
      PetscInt :: nconv,ncv,mpd,ni(1),mi(1)
      integer :: ierr
      
      PetscInt, allocatable :: idxr(:),idxc(:),idxvec(:)
      EPSProblemType :: epsproblem
      EPS :: eps
      EPSType :: tname
      PetscScalar kr(1), ki(1)
      PetscReal :: tol,error,norm
      Vec :: xr, xi,xd,xdd
      PetscScalar, pointer :: xx_v(:)
    
      
      character(256) :: sdummy

      real(kind=8) :: info(MAT_INFO_SIZE)    
    
    

      
      call VecSetOption(ew, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE,ierr)
    
      if (ev.ne.PETSC_NULL_MAT) then
        call MatCreateVecs(ev,xr,PETSC_NULL_VEC,ierr)      
        call MatCreateVecs(ev,xi,PETSC_NULL_VEC,ierr)
      end if
            
      if (b.eq.PETSC_NULL_MAT) then
        epsproblem=EPS_HEP
      else 
        epsproblem=EPS_GHEP
      end if
      
   !     ** Create eigensolver context      
      call EPSCreate(PETSC_COMM_WORLD,eps,ierr)

      
      call EPSSetConvergenceTest(eps,EPS_CONV_ABS,ierr)

      call EPSSetDimensions(eps,nev,PETSC_DEFAULT_INTEGER,&
     &  PETSC_DEFAULT_INTEGER,ierr)


!     ** Set operators. In this case, it is a standard eigenvalue problem
      call EPSSetOperators(eps,a,b,ierr)        
      call EPSSetProblemType(eps,epsproblem,ierr)
!~       call EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE,ierr)
      call EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL,ierr)

!     ** Set solver parameters at runtime
      call EPSSetFromOptions(eps,ierr)        

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ** Optional: Get some information from the solver and display it
      call EPSGetType(eps,tname,ierr)
 
 
      call EPSSolve(eps,ierr)
      call EPSGetConverged(eps,nconv,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      allocate(idxc(1),idxvec(1),stat=ierr)
      if (ierr.ne.0) then
        write(0,*) "allocation error ",ierr
        stop
      end if
      

      do i=0,nev-1
        if (ev.ne.PETSC_NULL_MAT) then
          call EPSGetEigenpair(eps,i,kr,PETSC_NULL_SCALAR,xr,PETSC_NULL_VEC,ierr)        
        else
          call EPSGetEigenpair(eps,i,kr,PETSC_NULL_SCALAR,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr)                  
        end if
        call EPSComputeError(eps,i,EPS_ERROR_RELATIVE,error,ierr)        
        
        call VecGetOwnershipRange(ew,ml1,ml2,ierr)       
!~         if ((i.ge.ml1).and.(i.le.ml2-1)) then
!~           write(6,*) "ev ",i,"=",kr
!~         end if

        if (ev.ne.PETSC_NULL_MAT) then
          call VecGetOwnershipRange(xr,nl1,nl2,ierr)       
          call VecGetArrayReadF90(xr,xx_v,ierr)        
          allocate(idxr(0:nl2-nl1-1),stat=ierr)        
          if (ierr.ne.0) then
            write(0,*) "allocation error ",ierr
            stop
          end if
          do j=0,nl2-nl1-1
            idxr(j)=j+nl1
          end do
            idxc(1)=i
        
          call MatSetValuesBlocked(ev,nl2-nl1,idxr,1,idxc,xx_v,INSERT_VALUES,ierr)

          deallocate(idxr)
          
        end if

        idxvec(1)=i*sign(1,ml2-i)*sign(1,i-ml1-1)
                     
        call VecSetValues(ew,1,idxvec,kr,INSERT_VALUES,ierr)
        

      end do
      call VecAssemblyBegin(ew,ierr)
      call VecAssemblyEnd(ew,ierr)  
      if (ev.ne.PETSC_NULL_MAT) then
        call MatAssemblyBegin(ev, MAT_FINAL_ASSEMBLY,ierr)    
        call MatAssemblyEnd(ev, MAT_FINAL_ASSEMBLY,ierr)      
        call VecRestoreArrayReadF90(xr,xx_v,ierr)
      end if


        

      if (ev.ne.PETSC_NULL_MAT) then
        call VecDestroy(xr,ierr)
        call VecDestroy(xi,ierr)
      end if


      call EPSDestroy(eps,ierr)

          
          
    end subroutine diag_mat

end module fctopetsc_mod

program fc_to_petsc
  use mpi
  use ISO_C_BINDING
  use petscmat
  use fctopetsc_mod
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscsys.h>
#include <petscsys.h>
  implicit none
  
  PetscScalar :: p_one=1d0, p_mone=-1d0
  real(kind=8) :: lattice(3,3),cutoff,vec(3),fcin,dh
  real(kind=8), allocatable :: atoms(:,:)
  integer :: iat1,iat2,nat,nz,ixyz1,ixyz2,ierr,iunit,n,m,np,inode,ii(1),jj(1),nl1,nl2,irun
  integer :: iroot
  character(256) :: file_coord,file_fc,file_petsc,instr,outstr
  Mat :: p_fc,p_fct,p_invsqrt_mass,p_dynmat,p_ev,p_mass_matrix
  logical :: l_ionode
  PetscScalar :: zfc(1),fcsum,fcdiag,p_scal,fcij,fcji
  PetscReal :: norm
  Vec :: ew
  PetscScalar, pointer :: xx_v(:)
  PetscViewer :: p_out_view
  
  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)

  if (command_argument_count().lt.4) then
    call PetscPrintf(PETSC_COMM_WORLD,"--- \n",ierr)
    call PetscPrintf(PETSC_COMM_WORLD,"\n",ierr)
    call PetscPrintf(PETSC_COMM_WORLD,"fc_to_petsc coord.in force_constant_file dh cutoff \n",ierr)
    call PetscPrintf(PETSC_COMM_WORLD,"\n",ierr)
    call PetscPrintf(PETSC_COMM_WORLD,"convert force_constant file from Conquest to PETSc format \n",ierr)
    call PetscPrintf(PETSC_COMM_WORLD,"\n",ierr)
    call PetscPrintf(PETSC_COMM_WORLD,"dh displacement use for finite difference (f(R+dh)-f(R))/dh \n",ierr)
    call PetscPrintf(PETSC_COMM_WORLD,"cutoff for force constant matrix elements \n",ierr)
    call SlepcFinalize(ierr)
    stop
  end if
  call get_command_argument(1,file_coord)
  call get_command_argument(2,file_fc)
  call get_command_argument(3,instr)
  read(instr,fmt=*) dh
  call get_command_argument(4,instr)
  read(instr,fmt=*) cutoff

  open(newunit=iunit,file=trim(file_coord),action="read",status="old")
  read(iunit,*) lattice(1:3,1)
  read(iunit,*) lattice(1:3,2)
  read(iunit,*) lattice(1:3,3)
  read(iunit,*) nat
  allocate(atoms(3,nat),stat=ierr)
  if (ierr.ne.0) then
    write(0,*) "allocation error ",ierr
    stop
  end if
  do iat1=1,nat
    read(iunit,fmt=*) atoms(1:3,iat1)
    atoms(1:3,iat1)=atoms(1,iat1)*lattice(1:3,1)+atoms(2,iat1)*lattice(1:3,2)+atoms(3,iat1)*lattice(1:3,3)
  end do
  close(iunit)
  
  call MPI_COMM_SIZE(PETSC_COMM_WORLD, np,ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD, inode, ierr )  
  l_ionode=inode.eq.0
  if (l_ionode) then
    write(6,*) "nat ",nat
    write(6,*) "cutoff ",cutoff
  end if
  
  n=3*nat
  m=n
  call MatCreate(PETSC_COMM_WORLD,p_fc,ierr)
  call MatSetType(p_fc, MATMPIAIJ,ierr)
  call MatSetSizes(p_fc,PETSC_DECIDE,PETSC_DECIDE,n,m,ierr)  
  call MatMPIAIJSetPreallocation(p_fc,0,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,ierr)
  call MatSetOption(p_fc,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,ierr)
  
  call MatCreate(PETSC_COMM_WORLD,p_invsqrt_mass,ierr)
  call MatSetType(p_invsqrt_mass, MATMPIAIJ,ierr)
  call MatSetSizes(p_invsqrt_mass,PETSC_DECIDE,PETSC_DECIDE,n,m,ierr)  
  call MatMPIAIJSetPreallocation(p_invsqrt_mass,1,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,ierr)
  call MatSetOption(p_invsqrt_mass,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,ierr)
  
  call MatCreate(PETSC_COMM_WORLD,p_mass_matrix,ierr)
  call MatSetType(p_mass_matrix, MATMPIAIJ,ierr)
  call MatSetSizes(p_mass_matrix,PETSC_DECIDE,PETSC_DECIDE,n,m,ierr)  
  call MatMPIAIJSetPreallocation(p_mass_matrix,1,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,ierr)
  call MatSetOption(p_mass_matrix,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,ierr)


  
  if (l_ionode) then
    open(newunit=iunit,file=trim(file_fc),action="read",status="old")
    nz=0
    do 
      read(unit=iunit,fmt=*,iostat=ierr) iat1,iat2,ixyz1,ixyz2,fcin
      if (ierr.ne.0) exit
      zfc=fcin
      vec=atoms(1:3,iat1)-atoms(1:3,iat2)
      if (dot_product(vec,vec).gt.cutoff*cutoff) cycle
      nz=nz+1
      ii=(iat1-1)*3+ixyz1-1
      jj=(iat2-1)*3+ixyz2-1
!~       write(6,fmt='(6i8,e24.12)') iat1,iat2,ixyz1,ixyz2,ii,jj,fcin
      call MatSetValues(p_fc,1,ii,1,jj,zfc,INSERT_VALUES,ierr)
    end do
    close(iunit)
    write(6,fmt='(A,2i8,e24.12)') "ndim, nz, nz/ndim ",n*n,nz,real(nz,8)/real(n*n,8)
    open(newunit=iunit,file="mass_au.dat",action="read",status="old")
    ii=-1
    do 
      read(unit=iunit,fmt=*,iostat=ierr) fcin
      if (ierr.ne.0) exit
      ii=ii+1
      zfc=fcin
      call MatSetValues(p_mass_matrix,1,ii,1,ii,zfc,INSERT_VALUES,ierr)
      zfc=1d0/dsqrt(fcin)      
      call MatSetValues(p_invsqrt_mass,1,ii,1,ii,zfc,INSERT_VALUES,ierr)
    end do    
    close(iunit)
  end if
  
  call MatAssemblyBegin(p_fc,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(p_fc,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyBegin(p_invsqrt_mass,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(p_invsqrt_mass,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyBegin(p_mass_matrix,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(p_mass_matrix,MAT_FINAL_ASSEMBLY,ierr)


  call MatGetOwnershipRange(p_fc,nl1,nl2,ierr)

  do iat2=1,nat    
    do ixyz2=1,3
      jj=(iat2-1)*3+ixyz2-1        
      do ixyz1=1,3
        fcsum=0d0
        fcdiag=0d0
        do iat1=1,nat        
          ii=(iat1-1)*3+ixyz1-1
          if (ii(1).ge.nl1.and.ii(1).le.nl2-1) then
            call MatGetValues(p_fc,1,ii,1,jj,zfc,ierr)
!~             write(6,fmt='(6i8,2e24.12)') iat1,iat2,ixyz1,ixyz2,ii,jj,zfc
            if (iat1.eq.iat2) then
              fcdiag=fcdiag+zfc(1)
            else
              fcsum=fcsum+zfc(1)              
            end if
          end if
        end do
        
        call MPI_Allreduce(MPI_IN_PLACE,fcsum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,PETSC_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,fcdiag,1,MPI_DOUBLE_COMPLEX,MPI_SUM,PETSC_COMM_WORLD,ierr)
!~         if (l_ionode) write(6,fmt='(A,3i8,4e24.12)') "fcsum ",iat2,ixyz1,ixyz2,fcsum,fcdiag
        ii=(iat2-1)*3+ixyz1-1
        if (ii(1).ge.nl1.and.ii(1).le.nl2-1) then
          zfc=-fcsum
          call MatSetValues(p_fc,1,ii,1,jj,zfc,INSERT_VALUES,ierr)
        end if
        call MatAssemblyBegin(p_fc,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(p_fc,MAT_FINAL_ASSEMBLY,ierr)
      end do
    end do
  end do
  


 
  
  do iat2=1,nat    
    do ixyz2=1,3
      jj=(iat2-1)*3+ixyz2-1        
      do iat1=1,nat
        do ixyz1=1,3
          ii=(iat1-1)*3+ixyz1-1
          if (ii(1).gt.jj(1)) cycle
          fcsum=0d0
          iroot=-1
          if (ii(1).ge.nl1.and.ii(1).le.nl2-1) then
            call MatGetValues(p_fc,1,ii,1,jj,zfc,ierr)
            fcij=zfc(1)
            iroot=inode
!~             write(6,fmt='(6i8,2e24.12)') iat1,iat2,ixyz1,ixyz2,ii,jj,zfc
          end if
          call MPI_Allreduce(MPI_IN_PLACE,iroot,1,MPI_INTEGER,MPI_MAX,PETSC_COMM_WORLD,ierr)
          call MPI_BCAST(fcij,1,MPI_DOUBLE_COMPLEX,iroot,PETSC_COMM_WORLD,ierr)
          iroot=-1
          if (jj(1).ge.nl1.and.jj(1).le.nl2-1) then
            call MatGetValues(p_fc,1,jj,1,ii,zfc,ierr)
            fcji=zfc(1)
            iroot=inode
!~             write(6,fmt='(6i8,2e24.12)') iat1,iat2,ixyz1,ixyz2,ii,jj,zfc
          end if
          call MPI_Allreduce(MPI_IN_PLACE,iroot,1,MPI_INTEGER,MPI_MAX,PETSC_COMM_WORLD,ierr)
          call MPI_BCAST(fcji,1,MPI_DOUBLE_COMPLEX,iroot,PETSC_COMM_WORLD,ierr)          
          zfc=(fcij+fcji)*0.5d0
!~           zfc=666d0
          if (ii(1).ge.nl1.and.ii(1).le.nl2-1) then
!~             write(6,fmt='(A,6i8,3e24.12)') "ii ",iat1,iat2,ixyz1,ixyz2,ii,jj,real(zfc),real(fcij),real(fcji)
            call MatSetValues(p_fc,1,ii,1,jj,zfc,INSERT_VALUES,ierr)
          end if
          if (jj(1).ge.nl1.and.jj(1).le.nl2-1) then
            call MatSetValues(p_fc,1,jj,1,ii,zfc,INSERT_VALUES,ierr)
!~             write(6,fmt='(A,6i8,3e24.12)') "jj ",iat2,iat1,ixyz2,ixyz1,jj,ii,real(zfc),real(fcij),real(fcji)
          end if
          call MatAssemblyBegin(p_fc,MAT_FINAL_ASSEMBLY,ierr)
          call MatAssemblyEnd(p_fc,MAT_FINAL_ASSEMBLY,ierr)            
        end do
        
      end do
    end do
  end do
  

  
  p_scal=1d0/dh
  call MatScale(p_fc,p_scal,ierr)
  
  do iat2=1,nat    
    do ixyz2=1,3
      jj=(iat2-1)*3+ixyz2-1        
      do ixyz1=1,3
        do iat1=1,nat        
          ii=(iat1-1)*3+ixyz1-1
          if (ii(1).ge.nl1.and.ii(1).le.nl2-1) then
            call MatGetValues(p_fc,1,ii,1,jj,zfc,ierr)
!~             write(6,fmt='(A,6i8,2e24.12)') "symmetrized ",iat1,iat2,ixyz1,ixyz2,ii,jj,zfc
          end if
        end do
      end do
    end do
  end do
  
  call MatPtAP(p_fc,p_invsqrt_mass,MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL,p_dynmat,ierr)
  
  call MatCreateDense(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,m,PETSC_NULL_SCALAR,p_ev,ierr)
  call MatCreateVecs(p_ev,ew,PETSC_NULL_VEC,ierr)
  
  call diag_mat(p_dynmat,PETSC_NULL_MAT,ew,p_ev,n,ierr)
  call VecGetArrayReadF90(ew,xx_v,ierr)        
  call VecGetOwnershipRange(ew,nl1,nl2,ierr)       
  do iat1=1,n
    if (iat1.ge.nl1+1.and.iat1.le.nl2) then      
      write(outstr,fmt='(i8,4e24.12,"\n")') iat1,xx_v(iat1-nl1),zsqrt(xx_v(iat1-nl1))*27.2114d0
      call PetscSynchronizedPrintf(PETSC_COMM_WORLD,outstr,ierr)
    end if
  end do
  call PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT,ierr)
  call VecRestoreArrayReadF90(ew,xx_v,ierr)        
  
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"sqrtMinv_petsc.dat", FILE_MODE_WRITE,p_out_view,ierr)
  call Matview(p_invsqrt_mass,p_out_view,ierr)
  call PetscViewerDestroy(p_out_view,ierr)
  
  
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Massmat_petsc.dat", FILE_MODE_WRITE,p_out_view,ierr)
  call Matview(p_mass_matrix,p_out_view,ierr)
  call PetscViewerDestroy(p_out_view,ierr)
  
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"FCmat_0_0_0_petsc.dat", FILE_MODE_WRITE,p_out_view,ierr)
  call Matview(p_fc,p_out_view,ierr)
  call PetscViewerDestroy(p_out_view,ierr)
  
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Dynmat_0_0_0_petsc.dat", FILE_MODE_WRITE,p_out_view,ierr)
  call Matview(p_dynmat,p_out_view,ierr)
  call PetscViewerDestroy(p_out_view,ierr)
  call SlepcFinalize(ierr)
    
end program fc_to_petsc
