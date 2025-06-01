module ome_ex
  use constants_math
  use parser_input_file, &
  only:nf,e1,e2,eta,nw
  use parser_wannier90_tb, &
  only:material_name,norb
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  norb_ex_cut,nv_ex,nc_ex,nband_ex,e_ex,fk_ex, &
  get_ex_index_first,print_exciton_wf, & !routines
  rkxvector,rkyvector,rkzvector !k-vectors only used for testing
  implicit none

  !ex-vme
  allocatable :: vme_ex(:,:) !V_{0N} 
  allocatable :: vme_ex_inter(:,:,:) !V_{NN'}
  
  !ex-rme
  allocatable :: xme_ex(:,:) !R_{0N}  
  allocatable :: xme_ex_inter(:,:,:) !R_{NN'}
  allocatable :: yme_ex_inter(:,:,:) !Y_{NN'}
  allocatable :: qme_ex_inter(:,:,:) !Q_{NN'}
  
  complex*16 :: vme_ex,vme_ex_inter
  complex*16 :: xme_ex,xme_ex_inter,yme_ex_inter,qme_ex_inter  
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_ex(iflag_norder)
    implicit none
    integer iflag_norder
    integer :: ibz
    
    !auxiliary arrays used to evaluate ex-ome
    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex)
    dimension shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex)

    real*8 :: ek
    real*8 :: shift_vector_ex_band
    complex*16 vme_ex_band
    complex*16 berry_eigen_ex_band
    complex*16 gen_der_ex_band
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) '6. Entering ome_ex'
    !allocate to read k-dependent vme and energies between involved band
    vme_ex_band=0.0d0
    ek=0.0d0

    !allocate arrays for ex-ome
    !linear conductivity
    if (iflag_norder.eq.1 .or. iflag_norder.eq.2) then
      allocate (vme_ex(3,norb_ex_cut))
      allocate (xme_ex(3,norb_ex_cut))
      vme_ex=0.0d0
      xme_ex=0.0d0
    end if
    !second order ones
    if (iflag_norder.eq.2) then
      allocate (yme_ex_inter(3,norb_ex_cut,norb_ex_cut))
      allocate (qme_ex_inter(3,norb_ex_cut,norb_ex_cut))
      allocate (xme_ex_inter(3,norb_ex_cut,norb_ex_cut))
      yme_ex_inter=0.0d0
      qme_ex_inter=0.0d0
      xme_ex_inter=0.0d0
    end if
    
    !read SP optical matrix elements from file
    write(*,*) '   Reading optical matrix elements (sp)'
    if (iflag_norder.eq.1) then
      call read_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
    end if
    if (iflag_norder.eq.2) then
      call read_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,berry_eigen_ex_band, &
                                   gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek)
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
    !k-space integration of excitonic optical marix elements
    do ibz=1,npointstotal
      !fill only V_{0N} for linear conductivity 
      if (iflag_norder.eq.1) then
        call get_vme_ex_sum_k(ibz,vme_ex_band) !sum over k points
      end if
      !fill only V_{0N} 
      if (iflag_norder.eq.2) then
        !fill V_{0N} and others.... 
        !call get_vme_ex_sum_k(ibz,vme_ex_band) !sum over k points
      end if
    end do
    write(*,*) '   Optical matrix elements (ex) have been evaluated' 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !writing excitonic optical matrix elements
    if (iflag_norder.eq.1) then
      call write_ome_ex_linear(vme_ex)
    end if
    write(*,*) '   Optical matrix elements (ex) have been written in file'
  end subroutine get_ome_ex




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_vme_ex_sum_k(ibz,vme_ex_band)
    implicit none
    integer :: ibz

    integer :: nn,ic,iv,iright,nj
    integer :: i_ex_nn

    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)

    complex*16 :: vme_ex_band
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do nn=1,norb_ex_cut		  
      do ic=1,nc_ex
  	    do iv=1,nv_ex
          iright=0
          call get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex_nn,ic,iv)			
          do nj=1,3			   
  		      vme_ex(nj,nn)=vme_ex(nj,nn)+fk_ex(i_ex_nn,nn)*vme_ex_band(ibz,nj,iv,nv_ex+ic)							  				
          end do			  
        end do	
      end do	
    end do
  end subroutine get_vme_ex_sum_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_ome_ex_linear(vme_ex)
    implicit none
    integer nn,nj
    dimension vme_ex(3,norb_ex_cut)
    complex*16 vme_ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  open(10,file='ome_linear_ex_'//trim(material_name)//'.omeex') 
    write(10,*) 1    
    do nn=1,norb_ex_cut
      write(10,*) nn,(realpart(vme_ex(nj,nn)),aimag(vme_ex(nj,nn)), nj=1,3)
    end do
    close(10)
  end subroutine write_ome_ex_linear
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine read_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
    implicit none
    integer iflag_norder
    integer npointstotal,nband_ex
    integer ibz
    integer nj,i,j

    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    
    real*8 :: ek
    real*8 :: a1,a2,a3,b1,b2,b3,b4,b5,b6
    complex*16 vme_ex_band
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  open(10,file='ome_linear_sp_'//trim(material_name)//'.omesp')     
    read(10,*) iflag_norder
    do ibz=1,npointstotal
	    read(10,*) a1,a2,a3,(ek(ibz,j),j=1,nband_ex)	
	    do i=1,nband_ex			
	      do j=1,nband_ex
	        read(10,*) a1,a2,a3,b1,b2,b3,b4,b5,b6
	        vme_ex_band(ibz,1,i,j)=complex(b1,b2)
	        vme_ex_band(ibz,2,i,j)=complex(b3,b4)
	        vme_ex_band(ibz,3,i,j)=complex(b5,b6)
        end do
      end do
    end do

    close(10)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  end subroutine read_ome_sp_linear

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine read_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,berry_eigen_ex_band, &
                                   gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek)
    implicit none
    integer iflag_norder
    integer npointstotal,nband_ex
    integer ibz
    integer nj,i,j

    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex)
    dimension shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex)

    real*8 :: ek
    real*8 :: a1,a2,a3,b1,b2,b3,b4,b5,b6
    real*8 :: shift_vector_ex_band
    complex*16 vme_ex_band
    complex*16 berry_eigen_ex_band
    complex*16 gen_der_ex_band
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  open(10,file='ome_nonlinear_sp_'//trim(material_name)//'.omesp')     
    read(10,*) iflag_norder
    do ibz=1,npointstotal
      read(10,*) a1,a2,a3,(ek(ibz,j),j=1,nband_ex)	
      do i=1,nband_ex			
        do j=1,nband_ex
	        read(10,*) a1,a2,a3,b1,b2,b3,b4,b5,b6
	        vme_ex_band(ibz,1,i,j)=complex(b1,b2)
	        vme_ex_band(ibz,2,i,j)=complex(b3,b4)
	        vme_ex_band(ibz,3,i,j)=complex(b5,b6)
	        read(10,*) a1,a2,a3,b1,b2,b3,b4,b5,b6
	        berry_eigen_ex_band(ibz,1,i,j)=complex(b1,b2)
	        berry_eigen_ex_band(ibz,2,i,j)=complex(b3,b4)
	        berry_eigen_ex_band(ibz,3,i,j)=complex(b5,b6)
	        do nj=1,3
	          read(10,*) a1,a2,a3,b1,b2,b3
	          shift_vector_ex_band(ibz,nj,1,i,j)=b1
	          shift_vector_ex_band(ibz,nj,2,i,j)=b2
	          shift_vector_ex_band(ibz,nj,3,i,j)=b3
          end do
	        do nj=1,3
	          read(10,*) a1,a2,a3,b1,b2,b3,b4,b5,b6
	          gen_der_ex_band(ibz,nj,1,i,j)=complex(b1,b2)
	          gen_der_ex_band(ibz,nj,2,i,j)=complex(b3,b4)
	          gen_der_ex_band(ibz,nj,3,i,j)=complex(b5,b6)
          end do
	      end do
      end do
    end do	  
    close(10)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  end subroutine read_ome_sp_nonlinear
end module ome_ex