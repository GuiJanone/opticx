module sigma_first
  use constants_math
  use parser_input_file, &
  only:nf,e1,e2,eta,nw
  use parser_wannier90_tb, &
  only:material_name,norb
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  norb_ex_cut,nv_ex,nc_ex,nband_ex,e_ex,fk_ex
  use ome, &
  only:ek,vme_ex_band
  
   !subroutine
  !use parser_input_file, &
  !only:e1,e2,eta,nw 
  implicit none
  
  allocatable :: vme_ex(:,:)
  complex*16 :: vme_ex
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sigma_first()
    implicit none
    integer :: ibz
    
    dimension vme_nband(3,nband_ex,nband_ex)
    dimension e_nband(norb)
    
    dimension wp(nw)
    dimension sigma_w_sp(3,3,nw),sigma_w_ex(3,3,nw)
    
    real*8 :: wp
    real*8 e_nband
    complex*16 vme_nband
    complex*16 sigma_w_sp,sigma_w_ex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !eta=eta/27.211385d0
    allocate (vme_ex(3,norb_ex_cut))

    write(*,*) 'computing single-particle sigma1'
    do ibz=1,npointstotal
      write(*,*) 'sigma_sp point:',ibz,npointstotal       
      !fill auxiliary arrays
      e_nband(:)=ek(ibz,:)
      vme_nband(:,:,:)=vme_ex_band(ibz,:,:,:)
      !fill sigma(w) for a given k point
      call get_kubo_intens(nband_ex,npointstotal,vcell,e_nband,vme_nband,nw,wp,sigma_w_sp)
      !fill vme_ex^{(N)} for a given k point
      call get_vme_ex_sum_k(ibz)
    end do

    !get excitonic frequency tensor
    write(*,*) 'broadening excitons peaks...'
    call get_kubo_intens_ex(nw,wp,sigma_w_ex)	
    write(*,*) 'Printing sigma first...'
    call print_sigma_first(nw,wp,sigma_w_sp,sigma_w_ex)

  end subroutine get_sigma_first

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
    open(50,file='bin/sigma_first_sp_real'//trim(material_name)//'.dat')
    open(55,file='bin/sigma_first_sp_imag_'//trim(material_name)//'.dat')
    open(60,file='bin/sigma_first_ex_real_'//trim(material_name)//'.dat')
    open(65,file='bin/sigma_first_ex_imag_'//trim(material_name)//'.dat')
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
    end do

    close(50)
    close(55)
    close(60)
    close(65)
  end subroutine print_sigma_first
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_kubo_intens_ex(nw,wp,sigma_w_ex)
  implicit none
  dimension skubo_ex_int(3,norb_ex_cut,norb_ex_cut)
  dimension wp(nw),sigma_w_ex(3,3,nw)
  
  integer nw
  integer iw,nn,nj,njp
  real*8 delta_n_ex
  real*8 wp
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_vme_ex_sum_k(ibz)
    implicit none
    integer :: ibz

    integer :: nn,ic,iv,iright,nj
    integer :: i_ex_nn
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !this subroutine has not been updated in 2025
  subroutine get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex,ic,iv)
    implicit none
    integer nf,nv_ex,nc_ex,ibz,i_ex,ic,iv
    integer iright 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !get band indeces (respect the fermi level) from A_cv index
	  if (iright.eq.1) then
	    iv=(nf-nv_ex)+i_ex-int((i_ex-1)/nv_ex)*nv_ex-nf
	    ic=(nf+1)+int((i_ex-1)/nv_ex)-nf	
	  end if
    !get A_cv index from band indeces
	  if (iright.ne.1) then
      i_ex=nc_ex*nv_ex*(ibz-1)+nv_ex*(ic-1)+iv
	  end if	  	  
  end subroutine get_ex_index_first

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_kubo_intens(nband_ex,npointstotal,vcell,e,vme,nw,wp,sigma_w_sp)
    implicit none
    integer :: iw,nn,nnp,nj,njp
    integer :: nw
    integer :: nband_ex,npointstotal

    dimension :: vme(3,norb,norb)
    dimension :: e(norb)
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

