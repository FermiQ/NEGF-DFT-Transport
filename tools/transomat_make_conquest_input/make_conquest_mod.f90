module make_conquest_mod
  use ISO_FORTRAN_ENV
  implicit none
  
  real(kind=8), parameter :: eps=1d-6, pi=datan(1d0)*4d0
  
  contains
  
    subroutine make_coord(nat,xyz_coords,atms,nat_elec1,nat_elec2,&
   &bravais,cell_z,offset_z,coordfile)
      implicit none
      
      integer ::nat,nat_elec1,nat_elec2
      real(kind=8) :: xyz_coords(:,:),bravais(3,3),cell_z,offset_z
      integer :: atms(:,:)
      character(*) :: coordfile
      
      integer :: iunit,i
      
      open(newunit=iunit,file=trim(coordfile),action="write",status="replace")
      write(iunit,fmt='(3e24.12)') bravais(1:3,1)
      write(iunit,fmt='(3e24.12)') bravais(1:3,2)
      write(iunit,fmt='(3e24.12)') bravais(1:3,3)
      write(iunit,fmt='(i8)') nat_elec2-nat_elec1+1
      
      xyz_coords(1:nat,3)=xyz_coords(1:nat,3)-xyz_coords(nat_elec1,3)+offset_z
      
      do i=nat_elec1,nat_elec2
        write(iunit,fmt='(3e24.12,i4,3(X,A))') xyz_coords(i,1)/bravais(1,1),&
       &xyz_coords(i,2)/bravais(2,2),xyz_coords(i,3)/cell_z,atms(i,1),&
       &merge("T","F",atms(i,2).eq.1),merge("T","F",atms(i,3).eq.1),&
       &merge("T","F",atms(i,4).eq.1)
      end do
      
      close(iunit)
       
    end subroutine make_coord
  
    function compare_2strings(str1,str2)
      implicit none
      
      logical :: compare_2strings
      character(*) :: str1,str2
    
      character(len=len(str1)) :: str1_low
      character(len=len(str2)) :: str2_low
            
      str1_low=to_lower(str1)
      str2_low=to_lower(str2)
      
      compare_2strings=trim(adjustl(str1_low)).eq.trim(adjustl(str2_low))
                  
    
    end function compare_2strings
  
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

      character(30),parameter :: skip(14)=[ character(len=30) :: &
     &"negf.LoadKnegf","negf.CurrentDensity","negf.dx","negf.dy","negf.dz",&
     &"negf.SaveHSK","Grid.InBlockX","Grid.InBlockY","Grid.InBlockZ",&
     &"Grid.PointsAlongX","Grid.PointsAlongY","Grid.PointsAlongZ",&
     &"Grid.Gridsolver_bcz","Grid.Gridsolver" ]
      
      character(30),parameter :: change(3)=[ character(len=30) :: &
     &"Atom.SupportFunctionRange", "Diag.MPMeshZ", "General.NumberOfSpecies" ]
     
      nchange=size(change)
      
      open(newunit=iunit,file="Conquest_input",action="read",status="old")
      open(newunit=junit,file=trim(adjustl(outfile)),action="write",&
     &status="replace")
      
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
      
      read_input: do
        read(unit=iunit,iostat=ierr,fmt='(A)') sstr        
        if (ierr.ne.0) exit
        sstr=adjustl(sstr)
        sstr1=sstr(1:index(sstr," "))
        lfound=.false.
        do i=1,size(skip)
          lfound=compare_2strings(sstr1,skip(i))
          if (lfound) cycle read_input
        end do   
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
             
        if (.not.lfound) then
          
          do i=1,size(change)            
            lfound=compare_2strings(sstr1,change(i))            
            if (lfound) then
              write(unit=junit,fmt='(A,X,A)') trim(sstr1),trim(change_str(i))              
              exit
            end if             
          end do
          
        end if
        if (.not.lfound) write(unit=junit,fmt='(A)') trim(sstr)
      end do read_input
      
      close(iunit)
      close(junit)
      
    
    end subroutine dump_conquest_electrode
  
    subroutine dump_conquest(change_str,ngrid,nblock,outfile)
      implicit none
      
      character(20),parameter :: skip(14)=[ character(len=20) :: &
     &"negf.LoadKnegf","negf.CurrentDensity","negf.dx","negf.dy","negf.dz",&
     &"negf.SaveHSK","Grid.InBlockX","Grid.InBlockY","Grid.InBlockZ",&
     &"Grid.PointsAlongX","Grid.PointsAlongY","Grid.PointsAlongZ",&
     &"Grid.Gridsolver_bcz","Grid.Gridsolver"]
      
      character(30),parameter :: change(1)=[ character(len=30) :: &
     &"Atom.SupportFunctionRange" ]
     
      character(*) :: outfile
      character(256) :: change_str(*)
      integer :: nchange,ngrid(3),nblock(3)

      integer :: iunit,junit,ierr,i
      character(256) :: sstr,sstr1
      logical :: lfound
      
      nchange=size(change)
      open(newunit=iunit,file="Conquest_input",action="read",status="old")
      open(newunit=junit,file=trim(outfile),action="write",status="replace")
      
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
      
      do
        read(unit=iunit,iostat=ierr,fmt='(A)') sstr        
        if (ierr.ne.0) exit
        sstr=adjustl(sstr)
        sstr1=sstr(1:index(sstr," "))
        lfound=.false.
        do i=1,size(skip)
          lfound=compare_2strings(sstr1,skip(i))
          if (lfound) exit
        end do        
        if (.not.lfound) then
          do i=1,nchange          
            lfound=compare_2strings(sstr1,change(i))          
            if (lfound) then
              write(junit,fmt='(A,X,A)') trim(sstr1),trim(change_str(i))       
              exit       
            end if 
          end do
          
        end if
        if (.not.lfound) write(junit,fmt='(A)') trim(sstr)
      end do
      
      close(iunit)
      close(junit)
           
    end subroutine dump_conquest

  
    function get_block_size(ngrid)
      implicit none
      
      integer :: ngrid,get_block_size
      
      do get_block_size=5,1,-1
        if (mod(ngrid,get_block_size).eq.0) exit
      end do
      
      
    end function get_block_size
  
    subroutine get_grid(cutoff,cell_d,dgrid,ngrid,nblock,evenodd)      
      implicit none
      
      real(kind=8) :: cutoff,cell_d,dgrid
      integer :: ngrid,nblock,evenodd
      
      real(kind=8) :: sqk
      integer :: ngrid_tmp,i
      
      sqk=dsqrt(2d0*cutoff)      
      ngrid_tmp=nint(sqk*cell_d/pi)
! grid solver needs grid to be 2^n efficient. here adjust to atleast n=3
      ngrid=nint(real(ngrid_tmp,8)/real(8,8))*8+evenodd
      nblock=get_block_size(ngrid)
      dgrid=cell_d/real(ngrid,8)
      write(6,*) "ngrid ",ngrid
      write(6,*) "nblock ",nblock
      write(6,*) "dgrid ",dgrid
      write(6,*) "cutoff adjusted to ",(ngrid*pi/cell_d)**2*0.5d0
      
    end subroutine
  
    function to_lower(strIn) result(strOut)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
    ! Original author: Clive Page
    
        implicit none
    
        character(len=*), intent(in) :: strIn
        character(len=len(strIn)) :: strOut
        integer :: i,j
    
        do i = 1, len(strIn)
              j = iachar(strIn(i:i))
              if (j.ge. iachar("A") .and. j.le.iachar("Z") ) then
                  strOut(i:i) = achar(iachar(strIn(i:i))+32)
              else
                  strOut(i:i) = strIn(i:i)
              end if
        end do
    
    end function to_lower
  
    subroutine get_conquest_input_line(optstr,outstr,iunit,lfound,lrewind)    
      implicit none
      
      character(*) :: optstr,outstr
      integer :: iunit
      logical :: lfound,lrewind
      
      character(len=len(optstr)) :: optstr_low
      character(256) :: instr,instr_low
      integer :: ierr
      
      optstr_low=adjustl(to_lower(optstr))
      if (lrewind) rewind(iunit)
      lfound=.false.
      do 
        read(unit=iunit,iostat=ierr,fmt='(A)') instr        
        if (ierr.ne.0) exit
        instr_low=adjustl(to_lower(instr))
        if ((index(instr_low,optstr_low).ge.1).and.&
       &(index(instr_low,"#").eq.0)) then
          lfound=.true.
          outstr=instr
          return        
        end if
      end do
      if (.not.lfound) return    

      
    end subroutine get_conquest_input_line
    
    subroutine get_conquest_real_option(optstr,optval,iunit,lfound,lrewind)    
      implicit none
      
      character(*) :: optstr
      real(kind=8) :: optval      
      integer :: iunit
      logical :: lfound,lrewind
      
      character(len=len(optstr)) :: optstr_low
      character(256) :: instr
      integer :: ierr
      
  
      call get_conquest_input_line(optstr,instr,iunit,lfound,lrewind)    
      if (lfound) read(instr(index(instr,optstr)+len(trim(optstr)):),&
      &fmt=*) optval
      
    end subroutine get_conquest_real_option
    
    subroutine get_conquest_integer_option(optstr,optval,iunit,lfound,lrewind)    
      implicit none
      
      character(*) :: optstr
      integer :: optval      
      integer :: iunit
      logical :: lfound,lrewind
      
      character(len=len(optstr)) :: optstr_low
      character(256) :: instr
      integer :: ierr
      
  
      call get_conquest_input_line(optstr,instr,iunit,lfound,lrewind)    
      if (lfound) read(instr(index(instr,optstr)+len(trim(optstr)):),&
      &fmt=*) optval
      
    end subroutine get_conquest_integer_option
    
    subroutine check_unit_cell(xyz_coords,atms,nat,nat_left,nat_right,&
   &left_cell_z,right_cell_z)
      implicit none
      
      real(kind=8) :: xyz_coords(:,:),left_cell_z,right_cell_z
      integer :: nat,nat_left,nat_right,atms(:,:)
      
      xyz_coords(1:nat,3)=xyz_coords(1:nat,3)-xyz_coords(1,3)
      
      left_cell_z=xyz_coords(nat_left+1,3)-xyz_coords(1,3)
      right_cell_z=xyz_coords(nat_right+1,3)-xyz_coords(1,3)
      
      if (.not.(all(abs(xyz_coords(1:nat_left,1:2)-&
     &xyz_coords(nat_left+1:nat_left*2,1:2)).le.eps).and.&
     &(all(atms(1:nat_left,1).eq.&
     &atms(nat_left+1:nat_left*2,1))))) then
        write(0,*) "cannot determine electrode left unit cell"
        stop
      end if
      write(6,*) "left unit cell_z  ",left_cell_z
      
      if (.not.(all(abs(xyz_coords(nat-nat_right+1:nat,1:2)-&
     &xyz_coords(nat-nat_right*2+1:nat-nat_right,1:2)).le.eps).and.&
     &(all(atms(nat-nat_right+1:nat,1).eq.&
     &atms(nat-nat_right*2+1:nat-nat_right,1))))) then
        write(0,*) "cannot determine electrode right unit cell"
        stop
      end if
      write(6,*) "right unit cell_z ",right_cell_z
      
    
    end subroutine check_unit_cell
    
    subroutine get_layers_z(xyz_coords,z_left,z_right,nat_left,nat_right)
      implicit none
      
      real(kind=8) :: xyz_coords(:,:)      
      real(kind=8) :: z_left,z_right
      integer :: nat_left,nat_right
      integer :: nat,i,j,k,ierr
      real(kind=8) :: current_layer_z
      
      nat=size(xyz_coords,1)
      
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
  
    subroutine read_coord(coordfile,nat,bravais,xyz_coords,atms)
      implicit none
      
      character(256) :: coordfile
      integer :: nat
      real(kind=8) ::  bravais(3,3)
      real(kind=8), allocatable :: xyz_coords(:,:),xyz_left(:,:),xyz_right(:,:)
      integer, allocatable :: atms(:,:)
      
      integer :: iunit,ierr,i,j
      character(1) :: move123(3)
      
      open(newunit=iunit,file=trim(coordfile),action="read",status="old")
      read(iunit,*) bravais(1:3,1)
      read(iunit,*) bravais(1:3,2)
      read(iunit,*) bravais(1:3,3)
      read(iunit,*) nat
      
      
      allocate(xyz_coords(nat,3),atms(nat,4),stat=ierr)
      if (ierr.ne.0) then
        write(0,*) "allocation error read_coord: xyz_coords, atms ",ierr
        stop
      end if
      atms=0
      do i=1,nat
        read(iunit,*) xyz_coords(i,1:3),atms(i,1),move123(1:3)
        xyz_coords(i,1:3)=xyz_coords(i,1)*bravais(1:3,1)+&
       &xyz_coords(i,2)*bravais(1:3,2)+xyz_coords(i,3)*bravais(1:3,3)
        do j=1,3
          atms(i,1+j)=merge(1,0,"t".eq.move123(j).or."T".eq.move123(j))
        end do
      end do
      
    end subroutine read_coord
  
    subroutine sort_ab(a,b,n,m1,m2,k)
      implicit none
      
      integer :: n,m1,m2,k
      real(kind=8) :: a(*)
      integer :: b(*)
      real(kind=8) :: x(1:m1)
      integer :: y(1:m2)
      integer :: i,j,off,l
  
      
      off=(k-1)*n
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
  
    subroutine sort(a)
      implicit none
      integer :: i, j, n
      integer :: a(:), x
      n = size(a,1)
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
