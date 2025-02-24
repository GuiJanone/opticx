module ome
  use constants_math
  use parser_input_file, &
  only:nv_ex,nc_ex
  use parser_wannier90_tb, &
  only:nR,nRvec,norb,R,shop,hhop,rhop_c
  use parser_optics_xatu_dim, &
  only:npointstotal,rkxvector,rkyvector,rkzvector,nband_ex

  implicit none

  allocatable :: skernel(:,:)
  allocatable :: hkernel(:,:)
  allocatable :: sderkernel(:,:,:)
  allocatable :: hderkernel(:,:,:)
  allocatable :: akernel(:,:,:)

  allocatable :: hk_ev(:,:)
  allocatable :: e(:)
  allocatable :: vme(:,:,:)
  
  real(8) :: e
  complex(16) :: skernel,hkernel,sderkernel,hderkernel,akernel
  complex(16) :: hk_ev,vme
  
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome()
    implicit none
    integer ibz
    integer i,j,ii,jj,nj
    
    dimension vme_ex_band(npointstotal,3,norb,norb)
    dimension nband_index(nband_ex)

    real(8) rkxp,rkyp,rkzp
    complex(16) vme_ex_band
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    allocate (skernel(norb,norb))
    allocate (hkernel(norb,norb)) 
    allocate (sderkernel(3,norb,norb))
    allocate (hderkernel(3,norb,norb)) 
    allocate (akernel(3,norb,norb)) 
    
    allocate (hk_ev(norb,norb))
    allocate (e(norb)) 
    allocate (vme(3,norb,norb)) 

    hk_ev=0.0d0
    e=0.0d0
    vme=0.0d0

    !noffset_v=0	
	  !noffset_v=18  !wannier core + spin
	  !noffset_c=0
    !call get_band_indexes(nf_v,noffset_v,noffset_c,nv_ex,nc_ex,nband_ex,nband_index)
    call get_band_indexes()

    write(*,*) 'mapping the BZ...'

	  do ibz=1,npointstotal     
      write(*,*) 'point', ibz, npointstotal
      rkxp=rkxvector(ibz)
		  rkyp=rkyvector(ibz)
      rkzp=rkzvector(ibz)		
      !get matrices in the \alpha, \alpha' basis (orbitals,k)    		
      call get_vme_kernels_ome()
      call get_vme_eigen_ome()
      !call get_eigen_vme_ome(iwannier,norb,skernel,hkernel,akernel,hderkernel, &
      !pgaugekernel,hk_ev,e,pgauge,vjseudoa,vjseudob,vme) 
      !call get_vme_kernels_ome(rkxp,rkyp,rkzp,nR,nRvec,norb,R,hkernel,skernel,shop, &
      !hhop,rhop,sderhop,hderhop,sderkernel,hderkernel,akernel,pgaugekernel)
		  do i=1,nband_ex 
		    ii=nband_index(i)
	      ek(ibz,i)=e(ii)
		    do nj=1,3			
		      do j=1,nband_ex
			      jj=nband_index(j)		  
		  	    vme_ex_band(ibz,nj,i,j)=vme(nj,ii,jj)
		  	    !berry_eigen_ex_band(ibz,nj,i,j)=berry_eigen(nj,ii,jj)
			      !shift_vector_ex_band(ibz,nj,1,i,j)=shift_vector(nj,1,ii,jj)
			      !shift_vector_ex_band(ibz,nj,2,i,j)=shift_vector(nj,2,ii,jj)
			      !shift_vector_ex_band(ibz,nj,3,i,j)=shift_vector(nj,3,ii,jj)
            !gen_der_ex_band(ibz,nj,1,i,j)=gen_der(nj,1,ii,jj)
			      !gen_der_ex_band(ibz,nj,2,i,j)=gen_der(nj,2,ii,jj)
			      !gen_der_ex_band(ibz,nj,3,i,j)=gen_der(nj,3,ii,jj)
            !vme_der_ex_band(ibz,nj,1,i,j)=vme_der(nj,1,ii,jj)
			      !vme_der_ex_band(ibz,nj,2,i,j)=vme_der(nj,2,ii,jj)
			      !vme_der_ex_band(ibz,nj,3,i,j)=vme_der(nj,3,ii,jj)
		      end do
        end do
      end do
		  
    end do

    
  end subroutine get_ome

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !subroutine get_band_indexes(nf_v,nf_c,nv_ex,nc_ex,nband_ex,nband_index)
  subroutine get_band_indexes(nf_v,nf_c,nv_ex,nc_ex,nband_ex,nband_index)
	  implicit none
    integer nv_top,nc_top
	  dimension nband_index(nband_ex)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
	  !nband_index=0
	  !nv_top=nv-noffset_v
	  !nc_bot=nv+1 
    !kacum=1

	  nband_index=0
	  nv_top=nf_v-noffset_v
	  nc_bot=nf_v+1 
    kacum=1

	  do iv=1,nv_ex
	    nband_index(kacum)=nv_top-nv_ex+iv
		  kacum=kacum+1
	  end do

	  do ic=1,nc_ex
	    nband_index(kacum)=nc_bot-1+ic
		  kacum=kacum+1	  
	  end do

	end

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   	        
  ! This routine evaluates the alpha,alpha' matrices
  subroutine get_vme_kernels_ome()
    implicit none 

    dimension hderhop(3,nR,norb,norb)
    dimension sderhop(3,nR,norb,norb)
    
    integer ialpha
    integer ialphap
    integer iRp
    integer nj
    
    real(8) Rx,Ry,Rz
    real(8) rkx,rky,rkz

    complex(16) phase,factor
    complex(16) hderhop,sderhop,rhop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  	  	  
    hkernel=0.0d0

    do ialpha=1,norb
      do ialphap=1,ialpha
        do iRp=1,nR
          Rx=dble(nRvec(iRp,1))*R(1,1)+dble(nRvec(iRp,2))*R(2,1)
          Ry=dble(nRvec(iRp,1))*R(1,2)+dble(nRvec(iRp,2))*R(2,2)
          Rz=0.0d0
          phase=complex(0.0d0,rkx*Rx+rky*Ry+rkz*Rz)
          factor=exp(phase)     
          
          hkernel(ialpha,ialphap)=hkernel(ialpha,ialphap)+ &
          factor*hhop(iRp,ialpha,ialphap)                
          skernel(ialpha,ialphap)=skernel(ialpha,ialphap)+ &
          factor*shop(iRp,ialpha,ialphap)   
       
          do nj=1,3 
            hderhop(nj,iRp,ialpha,ialphap)=complex(0.0d0,Rx)*hhop(iRp,ialpha,ialphap)
            sderhop(nj,iRp,ialpha,ialphap)=complex(0.0d0,Rx)*shop(iRp,ialpha,ialphap)

            sderkernel(nj,ialpha,ialphap)=sderkernel(nj,ialpha,ialphap)+ &
            factor*sderhop(nj,iRp,ialpha,ialphap)    
            hderkernel(nj,ialpha,ialphap)=hderkernel(nj,ialpha,ialphap)+ &
            factor*hderhop(nj,iRp,ialpha,ialphap) 
            akernel(nj,ialpha,ialphap)=akernel(nj,ialpha,ialphap)+ &
            factor*(rhop_c(nj,iRp,ialpha,ialphap)+ &
            complex(0.0d0,1.0d0)*sderhop(nj,iRp,ialpha,ialphap))            
          end do 
        end do   
      end do
    end do
  end subroutine get_vme_kernels_ome

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   	        
  ! This routine evaluates the alpha,alpha' matrices
  subroutine get_vme_eigen_ome()
    implicit none 

    integer :: ialpha
    integer :: ialphap
    integer :: iRp
    integer :: nj
    integer :: i,j,ii,jj,nn,nnp
    
    dimension hk_ev(norb,norb)
    dimension :: ecomplex(norb)
    dimension :: vjseudoa(3,norb,norb)
    dimension :: vjseudob(3,norb,norb)
    
    real(8) :: Rx,Ry,Rz
    real(8) :: rkx,rky,rkz
  	real(8) :: arg

    complex(16) :: vjseudoa,vjseudob
    complex(16) :: ecomplex
    complex(16) :: hk_ev
    complex(16) :: amu,amup,aux1,factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
    e=0.0d0
    call diagoz(norb,e,hkernel) 
    do ii=1,norb
      do jj=1,norb
        hk_ev(ii,jj)=hkernel(ii,jj)
      end do
    end do

    !phase election: this is done to smooth the gauge
    do j=1,norb
      aux1=0.0d0
      do i=1,norb
        aux1=aux1+hk_ev(i,j)
      end do
      !argument of the sym
      arg=atan2(aimag(aux1),realpart(aux1))
      factor=exp(complex(0.0d0,-arg))
      !write(*,*) 'sum is now:',aux1*factor
      do ii=1,norb
        hk_ev(ii,j)=hk_ev(ii,j)*factor
    	  !hk_ev(ii,j)=hk_ev(ii,j)*1.0d0
      end do      
    end do   
    
    !write(*,*) 'computing velocity matrix elements'      
    vme=0.0d0
    do nn=1,norb
      do nnp=1,nn 
      !momentums and A and B term
        do ialpha=1,norb
          do ialphap=1,norb
            amu=hk_ev(ialpha,nn)
            amup=hk_ev(ialphap,nnp)
            do nj=1,3          
              vjseudoa(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+ &
              conjg(amu)*amup*hderkernel(nj,ialpha,ialphap)
             
              vjseudob(nj,nn,nnp)=vjseudob(nj,nn,nnp)+conjg(amu)*amup* &
              (e(nn)*akernel(nj,ialpha,ialphap)-e(nnp)*conjg(akernel(nj,ialphap,ialpha)))* &
              complex(0.0d0,1.0d0)         
            end do                      
          end do
        end do

      end do
    end do
    
    do nj=1,3
      vme(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+vjseudob(nj,nn,nnp)
      vme(nj,nnp,nn)=conjg(vme(nj,nn,nnp))
      vjseudoa(nj,nnp,nn)=conjg(vjseudoa(nj,nn,nnp))
      vjseudob(nj,nnp,nn)=conjg(vjseudob(nj,nn,nnp))
    end do
              
  end subroutine get_vme_eigen_ome
end module ome