module parser_wannier90_tb
  use parser_input_file, only:read_line_numbers_int
  implicit none
  private
  public :: norb,nR,n1,n2,n3,nRvec
  public :: R,shop,hhop,rhop_c
  public :: material_name
  public :: wannier90_get

  integer norb
  integer nR
  integer n1,n2,n3
  integer nRvec

  real(8) R,shop
  complex*16 hhop,rhop_c
  character(100) material_name

  dimension R(3,3)
  allocatable nRvec(:,:)
  allocatable hhop(:,:,:)
  allocatable rhop_c(:,:,:,:)
  allocatable shop(:,:,:)

  contains
  
  subroutine wannier90_get(material_name_in)
    implicit none
    character(100), intent(in) :: material_name_in
    integer iskip
    integer nkk1,nkk2
    integer ialpha,ialphap
    integer iR
    integer nRzero
    integer nalloc,nallocate,nskip
    allocatable nallocate(:)
    real(8) a1,a2,a3,a4,a5,a6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) '2. Entering parser_wannier90_tb'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    material_name = material_name_in !quick fix. "_in" variable is deprecated
                                     !and should be removed in the future
    !open the w90 file 
    open(10,file='../wannier90_files_input/'//trim(material_name)//'_tb.dat')
    read(10,*) nskip
    read(10,*) R(1,1),R(1,2),R(1,3)
    read(10,*) R(2,1),R(2,2),R(2,3)
    read(10,*) R(3,1),R(3,2),R(3,3)
    read(10,*) norb
    read(10,*) nR
    close(10)
    !allocate nR, h-,r-, and s-hoppings
    allocate (nRvec(nR,3))
    allocate (hhop(nR,norb,norb))
    allocate (shop(nR,norb,norb))
    allocate (rhop_c(3,nR,norb,norb))

    !get where the hoppings start
    open(10,file='../wannier90_files_input/'//trim(material_name)//'_tb.dat')
    do iskip=1,nskip
      read(10,*)
    end do   
    !get the hopping matrices
    do iR=1,nR
      read(10,*) nRvec(iR,1),nRvec(iR,2),nRvec(iR,3)
      do ialphap=1,norb
	      do ialpha=1,norb
		      read(10,*) nkk1,nkk2,a1,a2
			    hhop(iR,ialpha,ialphap)=complex(a1,a2)
		    end do	  
	    end do
      read(10,*)
    end do
    
    !get rhoppings
   	do iR=1,nR
      read(10,*) !nRvec is already strored
	    do ialphap=1,norb
	      do ialpha=1,norb
		      read(10,*) nkk1,nkk2,a1,a2,a3,a4,a5,a6
          rhop_c(1,iR,ialpha,ialphap)=complex(a1,a2)
			    rhop_c(2,iR,ialpha,ialphap)=complex(a3,a4)
			    rhop_c(3,iR,ialpha,ialphap)=complex(a5,a6)
			  end do
      end do
      if (iR /= nR) read(10,*) !blank line
    end do
    close(10)
    
    !get orthogonal overlap: this variable is a reminiscent
    !of the interface with the original crystal interface.
    !I maintain the overlap matrix though

    !locate the (0,0,0) element of nRvec
    do iR=1,nR
      if (nRvec(iR,1)==0 .and. nRvec(iR,2)==0 .and. nRvec(iR,3)==0) then
        nRzero=iR
        exit
      end if
    end do
    !wannier functions are orthonormal
	  shop=0.0d0  
	  do ialpha=1,norb
	    shop(nRzero,ialpha,ialpha)=1.0d0
	  end do 
    !convert units: to Hartree and bohrs
	  hhop=hhop/27.211385d0
	  rhop_c=rhop_c/0.52917721067121d0
		R=R/0.52917721067121d0
	  write(*,*) '   Wannier hamiltonian has been read' 
  end subroutine wannier90_get
  

  !subroutine get_a(value)
    !real, intent(out) :: value
    !value = a
  !end subroutine get_a
  
end module parser_wannier90_tb