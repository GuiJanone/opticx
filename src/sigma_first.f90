module sigma_first
  use constants_math
  use parser_input_file, &
  only:iflag_xatu,nf,e1,e2,eta,nw
  use parser_wannier90_tb, &
  only:material_name,norb
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  norb_ex_cut,nv_ex,nc_ex,nband_ex,e_ex,fk_ex, &
  get_ex_index_first,print_exciton_wf, & !routines
  rkxvector,rkyvector,rkzvector !k-vectors only used for testing
  use ome_ex, &
  only:read_ome_sp_linear !routine
  !use ome_sp, &
  !only:ek,vme_ex_band
  !subroutine
  !use parser_input_file, &
  !only:e1,e2,eta,nw 
  implicit none

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sigma_first()
    implicit none 
    integer iflag_norder
    integer :: ibz,j
    
    !energies and vme in k-mesh and auxiliary arrays (sp)
    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension vme_nband(3,nband_ex,nband_ex)
    dimension e_nband(nband_ex)
    
    !energies and VME (ex)
    dimension e_ex(norb_ex_cut)
    dimension vme_ex(3,norb_ex_cut)

    dimension wp(nw)
    dimension sigma_w_sp(3,3,nw),sigma_w_ex(3,3,nw)
    
    real*8 wp
    real*8 ek,e_nband
    real*8 e_ex
    complex*16 vme_ex_band,vme_nband 
    complex*16 vme_ex
    complex*16 sigma_w_sp,sigma_w_ex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) '8. Entering sigma_first'
    !initialize sp arrays
    vme_ex_band=0.0d0
    ek=0.0d0   
    e_nband=0.0d0
    vme_nband=0.0d0
    !initialize ex arrays
    vme_ex=0.0d0
    
    !read optical matrix elements from file
    write(*,*) '   Reading optical matrix elements...'
    call read_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)

    if (iflag_xatu .eqv. .true.) then
      call read_ome_ex_linear(vme_ex)
    end if

    !allocate conductivity arrays
    call fill_allocate_sigma_arrays(nw,wp,sigma_w_sp,sigma_w_ex)
    
    !open(50,file='coefs_new_sigma.dat')
    write(*,*) '   Evaluating linear conductivity...'

    do ibz=1,npointstotal
      !write(*,*) 'sigma_sp point:',ibz,npointstotal       
      !fill auxiliary arrays
      e_nband(:)=ek(ibz,:)
      vme_nband(:,:,:)=vme_ex_band(ibz,:,:,:)
      !write(50,*) rkxvector(ibz),rkyvector(ibz),e_nband(2)
      !fill sigma(w) for a given k point
      call get_kubo_intens(nband_ex,npointstotal,vcell,e_nband,vme_nband,nw,wp,sigma_w_sp)
    end do
    !close(50)

    !get excitonic frequency tensor
    if (iflag_xatu .eqv. .true.) then
      write(*,*) '   Broadening excitons peaks...'
      call get_kubo_intens_ex(vme_ex,nw,wp,sigma_w_ex)
    end if

    !print conductivity tensor
    write(*,*) '   Printing sigma first...'
    call print_sigma_first(nw,wp,sigma_w_sp,sigma_w_ex)

  end subroutine get_sigma_first

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_ome_ex_linear(vme_ex)
    implicit none
    integer nn,nkaka
    dimension vme_ex(3,norb_ex_cut)
    complex*16 vme_ex

    real*8 :: a1,a2,a3,a4,a5,a6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  open(10,file='ome_linear_ex_'//trim(material_name)//'.omeex') 
    read(10,*)     
    do nn=1,norb_ex_cut
      read(10,*) nkaka,a1,a2,a3,a4,a5,a6
      vme_ex(1,nn)=complex(a1,a2)
      vme_ex(2,nn)=complex(a3,a4)
      vme_ex(3,nn)=complex(a5,a6)
    end do
    close(10)
  end subroutine read_ome_ex_linear
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine fill_allocate_sigma_arrays(nw,wp,sigma_w_sp,sigma_w_ex)
  implicit none
  
  integer :: nw
  integer :: i

  dimension :: wp(nw)
  dimension :: sigma_w_sp(3,3,nw)
  dimension :: sigma_w_ex(3,3,nw)

  real(8) :: wrange,wp
  complex(8) :: sigma_w_sp,sigma_w_ex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  wp=0.0d0
  sigma_w_sp=0.0d0
  sigma_w_ex=0.0d0
  wrange=e2-e1
  do i=1,nw
    wp(i)=(e1+wrange/dble(nw)*dble(i-1))/27.211385d0
  end do  
  eta=eta/27.211385d0 !change units to hartree units

  
  end subroutine fill_allocate_sigma_arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_sigma_first(nw,wp,sigma_w_sp,sigma_w_ex)
    implicit none
    integer :: iw
    integer :: nw
    dimension :: wp(nw)
    dimension :: sigma_w_sp(3,3,nw)
    dimension :: sigma_w_ex(3,3,nw)
    
    real*8 :: wp,feps
    complex*16 :: sigma_w_sp,sigma_w_ex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !write frequency dependent conductivity	  
    open(50,file='sigma_first_sp_real_'//trim(material_name)//'.dat')
    open(55,file='sigma_first_sp_imag_'//trim(material_name)//'.dat')
    if (iflag_xatu .eqv. .true.) then
      open(60,file='sigma_first_ex_real_'//trim(material_name)//'.dat')
      open(65,file='sigma_first_ex_imag_'//trim(material_name)//'.dat')
    end if

    !do iw=1,nw
      !write(*,*) wp(iw),sigma_w_sp(1,2,iw)
    !end do

    do iw=1,nw
      feps=1.0d0 !use atomic units
      write(50,*) wp(iw)*27.211385d0,realpart(feps*sigma_w_sp(1,1,iw)), &
        realpart(feps*sigma_w_sp(1,2,iw)), &
        realpart(feps*sigma_w_sp(1,3,iw)), &
        realpart(feps*sigma_w_sp(2,1,iw)), &
        realpart(feps*sigma_w_sp(2,2,iw)), &
        realpart(feps*sigma_w_sp(2,3,iw)), &
        realpart(feps*sigma_w_sp(3,1,iw)), &
        realpart(feps*sigma_w_sp(3,2,iw)), &
        realpart(feps*sigma_w_sp(3,3,iw))
  
      write(55,*) wp(iw)*27.211385d0,aimag(feps*sigma_w_sp(1,1,iw)), &
          aimag(feps*sigma_w_sp(1,2,iw)), &
          aimag(feps*sigma_w_sp(1,3,iw)), &
          aimag(feps*sigma_w_sp(2,1,iw)), &
          aimag(feps*sigma_w_sp(2,2,iw)), &
          aimag(feps*sigma_w_sp(2,3,iw)), &
          aimag(feps*sigma_w_sp(3,1,iw)), &
          aimag(feps*sigma_w_sp(3,2,iw)), &
          aimag(feps*sigma_w_sp(3,3,iw))	
      if (iflag_xatu .eqv. .true.) then  
        write(60,*) wp(iw)*27.211385d0,realpart(feps*sigma_w_ex(1,1,iw)), &
          realpart(feps*sigma_w_ex(1,2,iw)), &
          realpart(feps*sigma_w_ex(1,3,iw)), &
          realpart(feps*sigma_w_ex(2,1,iw)), &
          realpart(feps*sigma_w_ex(2,2,iw)), &
          realpart(feps*sigma_w_ex(2,3,iw)), &
          realpart(feps*sigma_w_ex(3,1,iw)), &
          realpart(feps*sigma_w_ex(3,2,iw)), &	
          realpart(feps*sigma_w_ex(3,3,iw))
        write(65,*) wp(iw)*27.211385d0,aimag(feps*sigma_w_ex(1,1,iw)), &
          aimag(feps*sigma_w_ex(1,2,iw)), &
          aimag(feps*sigma_w_ex(1,3,iw)), &
          aimag(feps*sigma_w_ex(2,1,iw)), &
          aimag(feps*sigma_w_ex(2,2,iw)), &
          aimag(feps*sigma_w_ex(2,3,iw)), &
          aimag(feps*sigma_w_ex(3,1,iw)), &
          aimag(feps*sigma_w_ex(3,2,iw)), &	
          aimag(feps*sigma_w_ex(3,3,iw))
      end if
    end do

    close(50)
    close(55)
    if (iflag_xatu .eqv. .true.) then
      close(60)
      close(65)
    end if
  end subroutine print_sigma_first
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_kubo_intens_ex(vme_ex,nw,wp,sigma_w_ex)
  implicit none
  dimension skubo_ex_int(3,norb_ex_cut,norb_ex_cut)
  dimension wp(nw),sigma_w_ex(3,3,nw)
  dimension vme_ex(3,norb_ex_cut)
  
  integer nw
  integer iw,nn,nj,njp
  real*8 delta_n_ex
  real*8 wp
  
  complex*16 :: vme_ex
  complex*16 :: skubo_ex_int, sigma_w_ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
    skubo_ex_int=0.0d0
    sigma_w_ex=0.0d0

    do iw=1,nw  
      do nn=1,norb_ex_cut
	      do nj=1,3
	        do njp=1,3
            
            !N integrand
            skubo_ex_int(nj,njp,nn)=1.0d0/(dble(npointstotal)*vcell) &
            *conjg(vme_ex(nj,nn))*vme_ex(njp,nn)/e_ex(nn)   !pick the correct order of operators
            
            !at a given frequency
            !delta function
            delta_n_ex=pi*1.0d0/eta*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta**2)*(wp(iw)-e_ex(nn))**2)
            !sigma_w
			      sigma_w_ex(nj,njp,iw)=sigma_w_ex(nj,njp,iw)+skubo_ex_int(nj,njp,nn)*delta_n_ex
          
	        end do
	      end do
      end do
    end do  

  end subroutine get_kubo_intens_ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_kubo_intens(nband_ex,npointstotal,vcell,e,vme,nw,wp,sigma_w_sp)
    implicit none
    integer :: iw,nn,nnp,nj,njp
    integer :: nw
    integer :: nband_ex,npointstotal

    dimension :: vme(3,nband_ex,nband_ex)
    dimension :: e(nband_ex)
    dimension :: sigma_w_sp(3,3,nw)

    dimension :: wp(nw)

    real*8 :: fnn,fnnp,factor1,vcell
    real*8 :: wp
    real*8 :: e
    
    complex*16 vme,delta_nnp
    complex*16 vme_prod,sigma_w_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  

    do iw=1,nw  
      do nn=1,nband_ex
      !fermi disrtibution
        if (nn.le.nv_ex) then
          fnn=1.0d0
        else
          fnn=0.0d0
        end if
        do nnp=1,nband_ex
          !fermi distribution
          if (nnp.le.nv_ex) then
            fnnp=1.0d0
          else
            fnnp=0.0d0
          end if                    
          !DEDICE PREFACTOR WITH OCCUPATION
          if (abs(fnn-fnnp).lt.0.1d0) then 
            factor1=0.0d0         
          else
            factor1=(fnn-fnnp)/(e(nn)-e(nnp))
          end if
	        !gaussian broadening
	        delta_nnp=-1.0d0/eta*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta**2)*(wp(iw)-e(nn)+e(nnp))**2)  

          !save oscillator stregths
          do nj=1,3
            do njp=1,3   
              vme_prod=vme(nj,nn,nnp)*vme(njp,nnp,nn) 
              sigma_w_sp(nj,njp,iw)=sigma_w_sp(nj,njp,iw)+ &
    	        pi/(dble(npointstotal)*vcell)*factor1*vme_prod*delta_nnp 
            end do
          end do
    	  	  
        end do
      end do
    end do  

  end subroutine get_kubo_intens
end module sigma_first

