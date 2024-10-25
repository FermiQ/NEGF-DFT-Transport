module make_conquest_mod
  use ISO_FORTRAN_ENV
  implicit none
  
  ! Define parameters for machine precision and pi.
  real(kind=8), parameter :: eps=1d-6, pi=datan(1d0)*4d0
  
contains
  
  ! Subroutine to write coordinates to a file.
  !
  ! Args:
  !   nat: Total number of atoms.
  !   xyz_coords: Array of atomic coordinates (x, y, z).
  !   atms: Array of atom types and movement flags.
  !   nat_elec1: Index of the first atom in the first electrode.
  !   nat_elec2: Index of the last atom in the second electrode.
  !   bravais: Bravais lattice vectors.
  !   cell_z: z-component of the unit cell vector.
  !   offset_z: Offset in the z-direction.
  !   coordfile: Name of the output coordinate file.
  subroutine make_coord(nat,xyz_coords,atms,nat_elec1,nat_elec2,&
  &bravais,cell_z,offset_z,coordfile)
    implicit none
    
    integer ::nat,nat_elec1,nat_elec2
    real(kind=8) :: xyz_coords(:,:),bravais(3,3),cell_z,offset_z
    integer :: atms(:,:)
    character(*) :: coordfile
    
    integer :: iunit,i
    
    ! Open the output file for writing.
    open(newunit=iunit,file=trim(coordfile),action="write",status="replace")
    ! Write the Bravais lattice vectors to the file.
    write(iunit,fmt='(3e24.12)') bravais(1:3,1)
    write(iunit,fmt='(3e24.12)') bravais(1:3,2)
    write(iunit,fmt='(3e24.12)') bravais(1:3,3)
    ! Write the number of atoms in the electrode to the file.
    write(iunit,fmt='(i8)') nat_elec2-nat_elec1+1
    
    ! Adjust the z-coordinates of the atoms.
    xyz_coords(1:nat,3)=xyz_coords(1:nat,3)-xyz_coords(nat_elec1,3)+offset_z
    
    ! Write the atomic coordinates and types to the file.
    do i=nat_elec1,nat_elec2
      write(iunit,fmt='(3e24.12,i4,3(X,A))') xyz_coords(i,1)/bravais(1,1),&
     &xyz_coords(i,2)/bravais(2,2),xyz_coords(i,3)/cell_z,atms(i,1),&
     &merge("T","F",atms(i,2).eq.1),merge("T","F",atms(i,3).eq.1),&
     &merge("T","F",atms(i,4).eq.1)
    end do
    
    ! Close the output file.
    close(iunit)
     
  end subroutine make_coord
  
  ! Function to compare two strings, ignoring case.
  !
  ! Args:
  !   str1: First string.
  !   str2: Second string.
  !
  ! Returns:
  !   .true. if the strings are equal (ignoring case), .false. otherwise.
  function compare_2strings(str1,str2)
    implicit none
    
    logical :: compare_2strings
    character(*) :: str1,str2
  
    character(len=len(str1)) :: str1_low
    character(len=len(str2)) :: str2_low
          
    ! Convert strings to lowercase for case-insensitive comparison.
    str1_low=to_lower(str1)
    str2_low=to_lower(str2)
    
    ! Compare the trimmed and left-adjusted lowercase strings.
    compare_2strings=trim(adjustl(str1_low)).eq.trim(adjustl(str2_low))
                
  
  end function compare_2strings
  
  ! Subroutine to modify and write a Conquest input file for an electrode.
  !
  ! Args:
  !   change_str: Array of strings to replace in the input file.
  !   nspecies: Number of species.
  !   species: Array of species indices.
  !   nblock: Number of grid points in each direction.
  !   ngrid: Number of grid points in each direction.
  !   dgrid: Grid spacing in each direction.
  !   outfile: Name of the output Conquest input file.
  subroutine dump_conquest_electrode(change_str,nspecies,species,nblock,&
  &ngrid,dgrid,outfile)
    implicit none
  
    integer :: nspecies,species(nspecies),nblock(3),ngrid(3)
    real(kind=8) :: dgrid(3)
    character(*) :: outfile
    character(256) :: change_str(*)

    integer :: iunit,junit,ierr,i,nchange,j,k,ispecies
    real(kind=8) :: ddum
    character(256) :: sstr,sstr1,sspecies
    logical :: lfound

    ! Array of strings to skip during input file processing.
    character(30),parameter :: skip(14)=[ character(len=30) :: &
   &"negf.LoadKnegf","negf.CurrentDensity","negf.dx","negf.dy","negf.dz",&
   &"negf.SaveHSK","Grid.InBlockX","Grid.InBlockY","Grid.InBlockZ",&
   &"Grid.PointsAlongX","Grid.PointsAlongY","Grid.PointsAlongZ",&
   &"Grid.Gridsolver_bcz","Grid.Gridsolver" ]
    
    ! Array of strings to change in the input file.
    character(30),parameter :: change(3)=[ character(len=30) :: &
   &"Atom.SupportFunctionRange", "Diag.MPMeshZ", "General.NumberOfSpecies" ]
   
    nchange=size(change)
    
    ! Open input and output files.
    open(newunit=iunit,file="Conquest_input",action="read",status="old")
    open(newunit=junit,file=trim(adjustl(outfile)),action="write",&
   &status="replace")
    
    ! Write header and NEGF parameters to the output file.
    write(junit,fmt='(A)') "# negf begin ---------------"
    write(junit,fmt='(A)') "negf.SaveHSK t # save H, S, and K mat to disc"
    write(junit,fmt='(A,e24.12)') "negf.dx ",dgrid(1)
    write(junit,fmt='(A,e24.12)') "negf.dy ",dgrid(2)
    write(junit,fmt='(A,e24.12)') "negf.dz ",dgrid(3)
    write(junit,fmt='(A)') "Grid.Gridsolver f"      
    write(junit,fmt='(A,i8)') "Grid.PointsAlongX ",ngrid(1)
    write(junit,fmt='(A,i8)') "Grid.PointsAlongY ",ngrid(2)
    write(junit,fmt='(A,i8)') "Grid.PointsAlongZ ",ngrid(3)
    write(junit,fmt='(A,i8)') "Grid.InBlockX ",nblock(1)
    write(junit,fmt='(A,i8)') "Grid.InBlockY ",nblock(2)
    write(junit,fmt='(A,i8)') "Grid.InBlockZ ",nblock(3)      
    write(junit,fmt='(A)') "# negf end -----------------"
    
    ! Read and process the input file line by line.
    read_input: do
      read(unit=iunit,iostat=ierr,fmt='(A)') sstr        
      if (ierr.ne.0) exit
      sstr=adjustl(sstr)
      sstr1=sstr(1:index(sstr," "))
      lfound=.false.
      ! Skip lines containing specified keywords.
      do i=1,size(skip)
        lfound=compare_2strings(sstr1,skip(i))
        if (lfound) cycle read_input
      end do   
      ! Process ChemicalSpeciesLabel block.
      lfound=compare_2strings(sstr,"%block ChemicalSpeciesLabel")
      if (lfound) then
        write(unit=junit,fmt='(A)') "%block ChemicalSpeciesLabel"
        ispecies = 0
        do 
          read(unit=iunit,iostat=ierr,fmt='(A)') sstr        
          if (ierr.ne.0) exit
          sstr=adjustl(sstr)
          lfound=compare_2strings(sstr,"%endblock ChemicalSpeciesLabel")
          if (lfound) then
            write(junit,fmt='(A)') "%endblock ChemicalSpeciesLabel"
            exit
          end if                        
          write(6,*) trim(sstr)
          read(sstr,*) k,ddum,sspecies
          do j=1,nspecies
            if (k.eq.species(j)) then
              ispecies = ispecies + 1
              write(unit=junit,fmt='(i4,f12.6,X,A)') ispecies,ddum,&
             &trim(adjustl(sspecies))
            end if
          end do
          
        end do         
      end if
           
      ! Replace specified keywords with new values.
      if (.not.lfound) then
        do i=1,size(change)            
          lfound=compare_2strings(sstr1,change(i))            
          if (lfound) then
            write(unit=junit,fmt='(A,X,A)') trim(sstr1),trim(change_str(i))              
            exit
          end if             
        end do
        
      end if
      ! Write unchanged lines to the output file.
      if (.not.lfound) write(unit=junit,fmt='(A)') trim(sstr)
    end do read_input
    
    ! Close input and output files.
    close(iunit)
    close(junit)
    
  
  end subroutine dump_conquest_electrode
  
  ! Subroutine to modify and write a Conquest input file.
  !
  ! Args:
  !   change_str: Array of strings to replace in the input file.
  !   ngrid: Number of grid points in each direction.
  !   nblock: Number of blocks in each direction.
  !   outfile: Name of the output Conquest input file.
  subroutine dump_conquest(change_str,ngrid,nblock,outfile)
    implicit none
    
    ! Array of strings to skip during input file processing.
    character(20),parameter :: skip(14)=[ character(len=20) :: &
   &"negf.LoadKnegf","negf.CurrentDensity","negf.dx","negf.dy","negf.dz",&
   &"negf.SaveHSK","Grid.InBlockX","Grid.InBlockY","Grid.InBlockZ",&
   &"Grid.PointsAlongX","Grid.PointsAlongY","Grid.PointsAlongZ",&
   &"Grid.Gridsolver_bcz","Grid.Gridsolver"]
    
    ! Array of strings to change in the input file.
    character(30),parameter :: change(1)=[ character(len=30) :: &
   &"Atom.SupportFunctionRange" ]
   
    character(*) :: outfile
    character(256) :: change_str(*)
    integer :: nchange,ngrid(3),nblock(3)

    integer :: iunit,junit,ierr,i
    character(256) :: sstr,sstr1
    logical :: lfound
    
    nchange=size(change)
    ! Open input and output files.
    open(newunit=iunit,file="Conquest_input",action="read",status="old")
    open(newunit=junit,file=trim(outfile),action="write",status="replace")
    
    ! Write header and NEGF parameters to the output file.
    write(junit,fmt='(A)') "# negf begin ---------------"
    write(junit,fmt='(A)') "negf.SaveHSK t # save H, S, and K mat to disc"
    write(junit,fmt='(A)') "negf.LoadKnegf f # set to t for NEGF-SCF"
    write(junit,fmt='(A)') "negf.CurrentDensity f"
    write(junit,fmt='(A)') "negf.left_electrode.dir ../left_electrode/"
    write(junit,fmt='(A)') "negf.right_electrode.dir ../right_electrode/"
    write(junit,fmt='(A)') "negf.mul          0.000000000"
    write(junit,fmt='(A)') "negf.mur          0.000000000"            
    write(junit,fmt='(A)') "Grid.Gridsolver_bcz 2"
    write(junit,fmt='(A)') "Grid.Gridsolver t"      
    write(junit,fmt='(A,i8)') "Grid.PointsAlongX ",ngrid(1)
    write(junit,fmt='(A,i8)') "Grid.PointsAlongY ",ngrid(2)
    write(junit,fmt='(A,i8)') "Grid.PointsAlongZ ",ngrid(3)
    write(junit,fmt='(A,i8)') "Grid.InBlockX ",nblock(1)
    write(junit,fmt='(A,i8)') "Grid.InBlockY ",nblock(2)
    write(junit,fmt='(A,i8)') "Grid.InBlockZ ",nblock(3)
    write(junit,fmt='(A)') "# negf end -----------------"
    
    ! Read and process the input file line by line.
    do
      read(unit=iunit,iostat=ierr,fmt='(A)') sstr        
      if (ierr.ne.0) exit
      sstr=adjustl(sstr)
      sstr1=sstr(1:index(sstr," "))
      lfound=.false.
      ! Skip lines containing specified keywords.
      do i=1,size(skip)
        lfound=compare_2strings(sstr1,skip(i))
        if (lfound) exit
      end do        
      ! Replace specified keywords with new values.
      if (.not.lfound) then
        do i=1,nchange          
          lfound=compare_2strings(sstr1,change(i))          
          if (lfound) then
            write(junit,fmt='(A,X,A)') trim(sstr1),trim(change_str(i))       
            exit       
          end if 
        end do
        
      end if
      ! Write unchanged lines to the output file.
      if (.not.lfound) write(junit,fmt='(A)') trim(sstr)
    end do
    
    ! Close input and output files.
    close(iunit)
    close(junit)
         
  end subroutine dump_conquest

  
  ! Function to find the largest integer factor of ngrid that is a power of 2.
  !
  ! Args:
  !   ngrid: Integer value.
  !
  ! Returns:
  !   The largest integer factor of ngrid that is a power of 2.
  function get_block_size(ngrid)
    implicit none
    
    integer :: ngrid,get_block_size
    
    ! Iterate downwards from 5 to find the largest power-of-2 factor.
    do get_block_size=5,1,-1
      if (mod(ngrid,get_block_size).eq.0) exit
    end do
    
    
  end function get_block_size
  
  ! Subroutine to determine grid parameters based on cutoff energy and cell dimensions.
  !
  ! Args:
  !   cutoff: Cutoff energy.
  !   cell_d: Cell dimension.
  !   dgrid: Grid spacing.
  !   ngrid: Number of grid points.
  !   nblock: Number of blocks.
  !   evenodd: Flag to ensure ngrid is either even or odd.
  subroutine get_grid(cutoff,cell_d,dgrid,ngrid,nblock,evenodd)      
    implicit none
    
    real(kind=8) :: cutoff,cell_d,dgrid
    integer :: ngrid,nblock,evenodd
    
    real(kind=8) :: sqk
    integer :: ngrid_tmp,i
    
    ! Calculate the number of grid points based on cutoff energy.
    sqk=dsqrt(2d0*cutoff)      
    ngrid_tmp=nint(sqk*cell_d/pi)
    ! Adjust ngrid to be a multiple of 8 for efficiency.
    ngrid=nint(real(ngrid_tmp,8)/real(8,8))*8+evenodd
    ! Determine the number of blocks.
    nblock=get_block_size(ngrid)
    ! Calculate the grid spacing.
    dgrid=cell_d/real(ngrid,8)
    ! Print grid parameters.
    write(6,*) "ngrid ",ngrid
    write(6,*) "nblock ",nblock
    write(6,*) "dgrid ",dgrid
    write(6,*) "cutoff adjusted to ",(ngrid*pi/cell_d)**2*0.5d0
    
  end subroutine
  
  ! Function to convert a string to lowercase.
  !
  ! Args:
  !   strIn: Input string.
  !
  ! Returns:
  !   strOut: Lowercase version of the input string.
  function to_lower(strIn) result(strOut)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
    ! Original author: Clive Page
    
    implicit none
  
    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
    integer :: i,j
  
    ! Iterate through the string and convert uppercase characters to lowercase.
    do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j.ge. iachar("A") .and. j.le.iachar("Z") ) then
              strOut(i:i) = achar(iachar(strIn(i:i))+32)
          else
              strOut(i:i) = strIn(i:i)
          end if
    end do
  
  end function to_lower
  
  ! Subroutine to read a specific line from a Conquest input file.
  !
  ! Args:
  !   optstr: Option string to search for.
  !   outstr: Output string containing the found line.
  !   iunit: Input file unit number.
  !   lfound: Logical flag indicating whether the line was found.
  !   lrewind: Logical flag indicating whether to rewind the file.
  subroutine get_conquest_input_line(optstr,outstr,iunit,lfound,lrewind)    
    implicit none
    
    character(*) :: optstr,outstr
    integer :: iunit
    logical :: lfound,lrewind
    
    character(len=len(optstr)) :: optstr_low
    character(256) :: instr,instr_low
    integer :: ierr
    
    ! Convert option string to lowercase.
    optstr_low=adjustl(to_lower(optstr))
    ! Rewind the file if requested.
    if (lrewind) rewind(iunit)
    lfound=.false.
    ! Read the file line by line until the option string is found.
    do 
      read(unit=iunit,iostat=ierr,fmt='(A)') instr        
      if (ierr.ne.0) exit
      instr_low=adjustl(to_lower(instr))
      ! Check if the option string is found and the line is not a comment.
      if ((index(instr_low,optstr_low).ge.1).and.&
     &(index(instr_low,"#").eq.0)) then
        lfound=.true.
        outstr=instr
        return        
      end if
    end do
    ! Return if the option string is not found.
    if (.not.lfound) return    

    
  end subroutine get_conquest_input_line
  
  ! Subroutine to read a real-valued option from a Conquest input file.
  !
  ! Args:
  !   optstr: Option string to search for.
  !   optval: Output real value of the option.
  !   iunit: Input file unit number.
  !   lfound: Logical flag indicating whether the option was found.
  !   lrewind: Logical flag indicating whether to rewind the file.
  subroutine get_conquest_real_option(optstr,optval,iunit,lfound,lrewind)    
    implicit none
    
    character(*) :: optstr
    real(kind=8) :: optval      
    integer :: iunit
    logical :: lfound,lrewind
    
    character(len=len(optstr)) :: optstr_low
    character(256) :: instr
    integer :: ierr
    
  
    ! Read the line containing the option.
    call get_conquest_input_line(optstr,instr,iunit,lfound,lrewind)    
    ! Read the real value from the line if the option was found.
    if (lfound) read(instr(index(instr,optstr)+len(trim(optstr)):),&
    &fmt=*) optval
    
  end subroutine get_conquest_real_option
  
  ! Subroutine to read an integer-valued option from a Conquest input file.
  !
  ! Args:
  !   optstr: Option string to search for.
  !   optval: Output integer value of the option.
  !   iunit: Input file unit number.
  !   lfound: Logical flag indicating whether the option was found.
  !   lrewind: Logical flag indicating whether to rewind the file.
  subroutine get_conquest_integer_option(optstr,optval,iunit,lfound,lrewind)    
    implicit none
    
    character(*) :: optstr
    integer :: optval      
    integer :: iunit
    logical :: lfound,lrewind
    
    character(len=len(optstr)) :: optstr_low
    character(256) :: instr
    integer :: ierr
    
  
    ! Read the line containing the option.
    call get_conquest_input_line(optstr,instr,iunit,lfound,lrewind)    
    ! Read the integer value from the line if the option was found.
    if (lfound) read(instr(index(instr,optstr)+len(trim(optstr)):),&
    &fmt=*) optval
    
  end subroutine get_conquest_integer_option
  
  ! Subroutine to check the unit cell of the electrodes.
  !
  ! Args:
  !   xyz_coords: Array of atomic coordinates.
  !   atms: Array of atom types.
  !   nat: Total number of atoms.
  !   nat_left: Number of atoms in the left electrode.
  !   nat_right: Number of atoms in the right electrode.
  !   left_cell_z: z-dimension of the left unit cell.
  !   right_cell_z: z-dimension of the right unit cell.
  subroutine check_unit_cell(xyz_coords,atms,nat,nat_left,nat_right,&
  &left_cell_z,right_cell_z)
    implicit none
    
    real(kind=8) :: xyz_coords(:,:),left_cell_z,right_cell_z
    integer :: nat,nat_left,nat_right,atms(:,:)
    
    ! Adjust z-coordinates to start from 0.
    xyz_coords(1:nat,3)=xyz_coords(1:nat,3)-xyz_coords(1,3)
    
    ! Calculate the z-dimension of the unit cells.
    left_cell_z=xyz_coords(nat_left+1,3)-xyz_coords(1,3)
    right_cell_z=xyz_coords(nat_right+1,3)-xyz_coords(1,3)
    
    ! Check if the left electrode unit cell is consistent.
    if (.not.(all(abs(xyz_coords(1:nat_left,1:2)-&
   &xyz_coords(nat_left+1:nat_left*2,1:2)).le.eps).and.&
   &(all(atms(1:nat_left,1).eq.&
   &atms(nat_left+1:nat_left*2,1))))) then
      write(0,*) "cannot determine electrode left unit cell"
      stop
    end if
    write(6,*) "left unit cell_z  ",left_cell_z
    
    ! Check if the right electrode unit cell is consistent.
    if (.not.(all(abs(xyz_coords(nat-nat_right+1:nat,1:2)-&
   &xyz_coords(nat-nat_right*2+1:nat-nat_right,1:2)).le.eps).and.&
   &(all(atms(nat-nat_right+1:nat,1).eq.&
   &atms(nat-nat_right*2+1:nat-nat_right,1))))) then
      write(0,*) "cannot determine electrode right unit cell"
      stop
    end if
    write(6,*) "right unit cell_z ",right_cell_z
    
  
  end subroutine check_unit_cell
  
  ! Subroutine to determine the number of layers in the left and right electrodes.
  !
  ! Args:
  !   xyz_coords: Array of atomic coordinates.
  !   z_left: z-coordinate of the left electrode boundary.
  !   z_right: z-coordinate of the right electrode boundary.
  !   nat_left: Number of atoms in the left electrode.
  !   nat_right: Number of atoms in the right electrode.
  subroutine get_layers_z(xyz_coords,z_left,z_right,nat_left,nat_right)
    implicit none
    
    real(kind=8) :: xyz_coords(:,:)      
    real(kind=8) :: z_left,z_right
    integer :: nat_left,nat_right
    integer :: nat,i,j,k,ierr
    real(kind=8) :: current_layer_z
    
    nat=size(xyz_coords,1)
    
    ! Determine the number of layers in the left electrode.
    current_layer_z=xyz_coords(1,3)
    nat_left=0
    do i=1,nat        
      if ((.not.(abs(xyz_coords(i,3)-current_layer_z).le.eps)).and.&
     &(xyz_coords(i,3).le.z_left)) then
        current_layer_z=xyz_coords(i,3)      
      else if (xyz_coords(i,3).gt.z_left) then          
        exit
      end if
      nat_left=nat_left+1
    end do
          
    ! Determine the number of layers in the right electrode.
    current_layer_z=xyz_coords(nat,3)      
    nat_right=0
    do i=nat,1,-1
      if ((.not.(abs(xyz_coords(i,3)-current_layer_z).le.eps)).and.&
     &(xyz_coords(i,3).ge.z_right)) then
        current_layer_z=xyz_coords(i,3)
      else if (xyz_coords(i,3).lt.z_right) then
        k=i+1
        exit
      end if
      nat_right=nat_right+1
    end do
  
  end subroutine get_layers_z
  
  ! Subroutine to read atomic coordinates and types from a file.
  !
  ! Args:
  !   coordfile: Name of the coordinate file.
  !   nat: Number of atoms.
  !   bravais: Bravais lattice vectors.
  !   xyz_coords: Array of atomic coordinates.
  !   atms: Array of atom types and movement flags.
  subroutine read_coord(coordfile,nat,bravais,xyz_coords,atms)
    implicit none
    
    character(256) :: coordfile
    integer :: nat
    real(kind=8) ::  bravais(3,3)
    real(kind=8), allocatable :: xyz_coords(:,:),xyz_left(:,:),xyz_right(:,:)
    integer, allocatable :: atms(:,:)
    
    integer :: iunit,ierr,i,j
    character(1) :: move123(3)
    
    ! Open the coordinate file for reading.
    open(newunit=iunit,file=trim(coordfile),action="read",status="old")
    ! Read Bravais lattice vectors and number of atoms.
    read(iunit,*) bravais(1:3,1)
    read(iunit,*) bravais(1:3,2)
    read(iunit,*) bravais(1:3,3)
    read(iunit,*) nat
    
    
    ! Allocate memory for atomic coordinates and types.
    allocate(xyz_coords(nat,3),atms(nat,4),stat=ierr)
    if (ierr.ne.0) then
      write(0,*) "allocation error read_coord: xyz_coords, atms ",ierr
      stop
    end if
    atms=0
    ! Read atomic coordinates and types from the file.
    do i=1,nat
      read(iunit,*) xyz_coords(i,1:3),atms(i,1),move123(1:3)
      xyz_coords(i,1:3)=xyz_coords(i,1)*bravais(1:3,1)+&
     &xyz_coords(i,2)*bravais(1:3,2)+xyz_coords(i,3)*bravais(1:3,3)
      ! Set movement flags based on characters read from file.
      do j=1,3
        atms(i,1+j)=merge(1,0,"t".eq.move123(j).or."T".eq.move123(j))
      end do
    end do
    
  end subroutine read_coord
  
  ! Subroutine to sort two arrays simultaneously based on the values in the first array.
  !
  ! Args:
  !   a: Array of real numbers to sort.
  !   b: Array of integers to sort along with a.
  !   n: Length of the arrays.
  !   m1: Number of rows in a.
  !   m2: Number of rows in b.
  !   k: Index of the column to sort by.
  subroutine sort_ab(a,b,n,m1,m2,k)
    implicit none
    
    integer :: n,m1,m2,k
    real(kind=8) :: a(*)
    integer :: b(*)
    real(kind=8) :: x(1:m1)
    integer :: y(1:m2)
    integer :: i,j,off,l
  
    
    off=(k-1)*n
    ! Sort the arrays using insertion sort.
    do i=2,n
    
        do l=1,m1
          x(l) = a((l-1)*n+i)
        end do        
        
        do l=1,m2
          y(l) = b((l-1)*n+i)
        end do
        
        j=i-1
        do while (j.ge.1)
            if (a(off+j).le.x(k)) exit
            do l=1,m1
              a((l-1)*n+j+1)=a((l-1)*n+j)
            end do
            do l=1,m2
              b((l-1)*n+j+1)=b((l-1)*n+j)
            end do
            j=j-1
        end do
        do l=1,m1
          a((l-1)*n+j+1)=x(l)
        end do
        do l=1,m2
          b((l-1)*n+j+1)=y(l)
        end do
    end do      
  end subroutine sort_ab
  
  ! Subroutine to sort an integer array using insertion sort.
  !
  ! Args:
  !   a: Integer array to sort.
  subroutine sort(a)
    implicit none
    integer :: i, j, n
    integer :: a(:), x
    n = size(a,1)
    ! Sort the array using insertion sort.
    do i = 2, n
        x = a(i)
        j = i - 1
        do while (j >= 1)
            if (a(j) <= x) exit
            a(j + 1) = a(j)
            j = j - 1
        end do
        a(j + 1) = x
    end do
  end subroutine sort

end module make_conquest_mod