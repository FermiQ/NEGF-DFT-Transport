program transomat_make_conquest_input
  use make_conquest_mod
  implicit none
  
  real(kind=8) :: z_left,z_right,basis_radius_max,gridcutoff,grid(3),&
 &bravais(3,3),left_cell_z,right_cell_z,offset_left_elec,offset_right_elec,&
 &bravais_right(3,3),bravais_left(3,3),dgrid(3)
  integer :: nat,nat_left,nat_right,ngrid(3),nblock(3),nspecies,nspecies_left,&
 &nspecies_right,ngrid_left(3),ngrid_right(3),nblock_left(3),nblock_right(3)
  real(kind=8), allocatable :: xyz_coords(:,:),xyz_left(:,:),xyz_right(:,:)
  integer, allocatable :: atms(:,:),atms_left(:,:),atms_right(:,:),&
 &species_left(:),species_right(:), types_check(:)
  
  integer :: n,m,l,i,j,k
  integer :: iunit,ierr
  character(256) :: coordfile,instr,change_str(1),change_str_left(3),&
 &change_str_right(3)
  character(1) :: move123(3)  
  real(kind=8) :: cell_z,rout
  logical :: lfound
  
  
  if (command_argument_count().lt.3) then
    write(0,*) "transomat_make_conquest_input coord.in left_cut right_cut"
    stop
  end if
  
  call get_command_argument(1,coordfile)
  call get_command_argument(2,instr)
  read(instr,*) z_left
  call get_command_argument(3,instr)
  read(instr,*) z_right
  
  write(6,*) "z_left  ",z_left
  write(6,*) "z_right ",z_right
  
  call read_coord(coordfile,nat,bravais,xyz_coords,atms)


  
  call sort_ab(xyz_coords,atms,nat,3,4,1)
  call sort_ab(xyz_coords,atms,nat,3,4,2)
  call sort_ab(xyz_coords,atms,nat,3,4,3)


  call get_layers_z(xyz_coords,z_left,z_right,nat_left,nat_right)  
  call check_unit_cell(xyz_coords,atms,nat,nat_left,nat_right,&
 &left_cell_z,right_cell_z)
 
! adjust unit cell, i.e. center atoms and assume left and right d/2 from cell edge

  offset_left_elec=(xyz_coords(nat_left+1,3)-xyz_coords(nat_left,3))*0.5d0
  xyz_coords(1:nat,3)=xyz_coords(1:nat,3)+offset_left_elec
 
  offset_right_elec=(xyz_coords(nat-nat_right+1,3)-&
 &xyz_coords(nat-nat_right,3))*0.5d0  
  cell_z=xyz_coords(nat,3)+offset_right_elec
  bravais(1:3,3)=(/ 0d0, 0d0, cell_z /)
  write(6,*) "cell_x ",bravais(1:3,1)
  write(6,*) "cell_y ",bravais(1:3,2)
  write(6,*) "cell_z ",bravais(1:3,3)
  
  open(newunit=iunit,file="Conquest_input",action="read",status="old")
  
  basis_radius_max=0d0
  do 
    call get_conquest_real_option("Atom.SupportFunctionRange",rout,iunit,lfound,.false.)    
    basis_radius_max=max(basis_radius_max,rout)
    if (.not.(lfound)) exit
  end do
  write(change_str(1),fmt='(f12.6)') basis_radius_max
  change_str(1)=adjustl(change_str(1))
  write(6,*) "basis_radius_max ",trim(change_str(1))
  gridcutoff=50d0
  call get_conquest_real_option("Grid.GridCutoff",gridcutoff,iunit,lfound,.true.)      
  write(6,*) "gridcutoff ",gridcutoff
  write(6,*) "cell x grid:"
  call get_grid(gridcutoff,bravais(1,1),dgrid(1),ngrid(1),nblock(1),0) 
  write(6,*) "cell y grid:"
  call get_grid(gridcutoff,bravais(2,2),dgrid(2),ngrid(2),nblock(2),0)          
  write(6,*) "cell z grid:"
  call get_grid(gridcutoff,bravais(3,3),dgrid(3),ngrid(3),nblock(3),1)
  call get_conquest_integer_option("General.NumberOfSpecies",nspecies,iunit,&
 &lfound,.true.)      
  write(6,*) "nspecies ",nspecies
  close(iunit)
  
  call system("mkdir ecc")
  call system("mkdir left_electrode")
  call system("mkdir right_electrode")
  call dump_conquest(change_str,ngrid,nblock,"ecc/Conquest_input")


  allocate(xyz_left(nat_left,3),xyz_right(nat_right,3),atms_left(nat_left,4),&
 &atms_right(nat_right,4), types_check(nspecies))

  call make_coord(nat_left,xyz_coords,atms,1,nat,&
   &bravais,cell_z,offset_left_elec,"ecc/coord.in")

  
  xyz_left(1:nat_left,1:3)=xyz_coords(1:nat_left,1:3)
  atms_left(1:nat_left,1:4)=atms(1:nat_left,1:4)
  bravais_left=bravais
  bravais_left(1:3,3)=(/ 0.0d0, 0.0d0, left_cell_z /)

  nspecies_left=1
  types_check(1) = atms_left(1,1)
  do i=2,nat_left    
    do j = 1, nspecies_left      
      if (atms_left(i,1) .eq. types_check(j)) exit      
      if (j .eq. nspecies_left) then
        nspecies_left = nspecies_left + 1
        types_check(nspecies_left) = atms_left(i,1)
        exit
      end if
    end do
  end do
  write(6,*) "nspecies left",nspecies_left
  allocate(species_left(nspecies_left))
  species_left(1)=atms_left(1,1)  
  k = 1
  do i=2,nat_left    
    do j = 1, k
      if (atms_left(i,1) .eq. species_left(j)) exit
      if ( j .eq. k) then
        k = k +1
        species_left(k)=atms_left(i,1)      
        exit
      end if    
    end do
    if (k .eq. nspecies_left) exit
  end do  
  call sort(species_left)
  write(6,fmt='(i3,A,99i3)') nspecies_left," species left ",species_left
  
  do i = 1, nat_left
    do j = 1, nspecies_left
      if (atms_left(i,1) .eq. species_left(j)) then
        atms_left(i,1) = j
        exit
      end if
    end do
  end do
  
  call make_coord(nat_left,xyz_left,atms_left,1,nat_left,&
   &bravais_left,left_cell_z,offset_left_elec,"left_electrode/coord.in")
  ngrid_left=ngrid
  ngrid_left(3)=nint(bravais_left(3,3)/dgrid(3))
  ngrid_left(3)=ngrid_left(3) + mod(ngrid_left(3),2)
  nblock_left=nblock
  nblock_left(3)=get_block_size(ngrid_left(3))
  write(6,*) "ngrid_left ",ngrid_left
  write(6,*) "nblock_left ",nblock_left
  change_str_left(1)=change_str(1)
  change_str_left(2)="101"
  write(change_str_left(3),fmt='(i8)') nspecies_left
  change_str_left(3)=adjustl(change_str_left(3))
  call dump_conquest_electrode(change_str_left,nspecies_left,species_left,&
  nblock_left,ngrid_left,dgrid,"left_electrode/Conquest_input")  

  xyz_right(1:nat_right,1:3)=xyz_coords(nat-nat_right+1:nat,1:3)
  atms_right(1:nat_right,1:4)=atms(nat-nat_right+1:nat,1:4)
  bravais_right=bravais
  bravais_right(1:3,3)=(/ 0.0d0, 0.0d0, right_cell_z /)   
  
  nspecies_right=1
  types_check(1) = atms_right(1,1)
  do i=2,nat_right
    do j = 1, nspecies_right
      if (atms_right(i,1) .eq. types_check(j)) exit      
      if (j .eq. nspecies_right) then
        nspecies_right = nspecies_right + 1
        types_check(nspecies_right) = atms_right(i,1)
        exit
      end if
    end do
  end do
  write(6,*) "nspecies right",nspecies_right 
  
  allocate(species_right(nspecies_right))  
  species_right(1)=atms_right(1,1)  
  k = 1
  do i=2,nat_right    
    do j = 1, k
      if (atms_right(i,1) .eq. species_right(j)) exit
      if ( j .eq. k) then
        k = k +1
        species_right(k)=atms_right(i,1)      
        exit
      end if    
    end do
    if (k .eq. nspecies_right) exit
  end do
  call sort(species_right)
  write(6,fmt='(i3,A,99i3)') nspecies_right," species right ",species_right 
  
    
  do i = 1, nat_right
    do j = 1, nspecies_right
      if (atms_right(i,1) .eq. species_right(j)) then
        atms_right(i,1) = j
        exit
      end if
    end do
  end do 
    
  call make_coord(nat_right,xyz_right,atms_right,1,nat_right,&
   &bravais_right,right_cell_z,offset_right_elec,"right_electrode/coord.in")
  ngrid_right=ngrid
  ngrid_right(3)=nint(bravais_right(3,3)/dgrid(3))
  ngrid_right(3)=ngrid_right(3) + mod(ngrid_right(3), 2)
  nblock_right=nblock
  nblock_right(3)=get_block_size(ngrid_right(3))  
  write(6,*) "ngrid_right ",ngrid_right
  write(6,*) "nblock_right ",nblock_right
  change_str_right(1)=change_str(1)
  change_str_right(2)="101"
  write(change_str_right(3),fmt='(i8)') nspecies_right
  change_str_right(3)=adjustl(change_str_right(3))  
  call dump_conquest_electrode(change_str_left,nspecies_left,species_left,&
  nblock_left,ngrid_left,dgrid,"right_electrode/Conquest_input")   
!~   do i=1,nat
!~     write(6,fmt='(3e24.12,i8,X,3(l,X))') xyz_coords(i,1:3),atms(i,1),atms(i,2:4).eq.1
!~   end do
  
end program transomat_make_conquest_input

