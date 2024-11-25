module parser_wannier90_tb
  implicit none
  
  integer norb
  integer nR
  integer n1,n2,n3
  integer nRvec

  dimension R(3,3)
  allocatable nRvec(:,:)
  allocatable hhop(:,:,:)
  allocatable rhop_c(:,:,:,:)
  allocatable shop(:,:,:)

  real(8) R,shop
  complex(16) hhop,rhop_c
  character(100) material_name
  
 






  contains
  
  subroutine wannier90_get(material_name_in)
    implicit none
    character(100), intent(in) :: material_name_in
    integer iskip
    integer nkk1,nkk2
    integer ialpha,ialphap
    integer iR

    real(8) a1,a2,a3,a4,a5,a6

    material_name = material_name_in

    !open the w90 file 
    open(10,file='./wannier90_files_input/'//trim(material_name)//'_tb.dat')
    read(10,*) !first line is blank
    read(10,*) R(1,1),R(1,2),R(1,3)
    read(10,*) R(2,1),R(2,2),R(2,3)
    read(10,*) R(3,1),R(3,2),R(3,3)
    read(10,*) norb
    read(10,*) nR
    
    !allocate nR, h-,r-, and s-hoppings
    allocate (nRvec(3,nR))
    allocate (hhop(nR,norb,norb))
    allocate (shop(nR,norb,norb))
    allocate (rhop_c(3,nR,norb,norb))
    !get where the hoppings start
    do iskip=1,1000
      read(10,*) n1,n2,n3 
      if (n1<0 .and. n2<0) then !locate where hopping matrices start
        nRvec(1,1)=n1 !save first lattice integers
        nRvec(1,2)=n2
        nRvec(1,3)=n3
	      do ialphap=1,norb
	        do ialpha=1,norb
		        read(10,*) nkk1,nkk2,a1,a2
			      hhop(1,ialpha,ialphap)=complex(a1,a2)
		      end do	  
	      end do
        exit
      end if
    end do
    
    !get the rest matrices
    do iR=2,nR
      read(10,*)
      do ialphap=1,norb
	      do ialpha=1,norb
		      read(10,*) nkk1,nkk2,a1,a2
			    hhop(iR,ialpha,ialphap)=complex(a1,a2)
		    end do	  
	    end do
    end do
    
    !get rhoppings
   	do iR=1,nR
      read(10,*) !blank line
      read(10,*) !nRvec is already strored
	    do ialphap=1,norb
	      do ialpha=1,norb
		      read(10,*) nkk1,nkk2,a1,a2,a3,a4,a5,a6
          rhop_c(1,iR,ialpha,ialphap)=complex(a1,a2)
			    rhop_c(2,iR,ialpha,ialphap)=complex(a3,a4)
			    rhop_c(3,iR,ialpha,ialphap)=complex(a5,a6)
			  end do
      end do
    end do
    close(10)
 
	  shop=0.0d0  
	  do ialpha=1,norb
	    shop(1,ialpha,ialpha)=1.0d0
	  end do

	  write(*,*) 'wannier hamiltonian has been read'  
    !convert units: to Hartree and bohrs
	  !hhop=hhop/27.211385d0
	  !rhop_c=rhop_c/0.52917721067121d0


  end subroutine wannier90_get
  
  !subroutine get_a(value)
    !real, intent(out) :: value
    !value = a
  !end subroutine get_a
  
end module parser_wannier90_tb