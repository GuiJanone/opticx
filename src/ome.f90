module ome
  use constants_math
  use parser_wannier90_tb, &
  only:nR,nRvec,norb,R,shop,hhop,rhop_c
  use parser_optics_xatu_dim, &
  only:npointstotal,rkxvector,rkyvector,rkzvector, &
       nband_ex,nband_index,nv_ex,nc_ex

  implicit none
  
  allocatable vme_ex_band(:,:,:,:)
  allocatable ek(:,:)

  real*8 ek
  complex*16 vme_ex_band

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome()
    implicit none
    integer ibz
    integer i,j,ii,jj,nj

    dimension skernel(norb,norb)
    dimension hkernel(norb,norb)
    dimension sderkernel(3,norb,norb)
    dimension hderkernel(3,norb,norb)
    dimension akernel(3,norb,norb)

    dimension e(norb)
    dimension ecomplex(norb)
    dimension hk_ev(norb,norb)
    dimension vme(3,norb,norb)

    real(8) rkx,rky,rkz
    real(8) :: e

    complex*16 :: skernel,hkernel,sderkernel,hderkernel,akernel
    complex*16 :: hk_ev
    complex*16 :: ecomplex
    complex*16 :: vjseudoa,vjseudob  
    complex*16 vme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !save k-dependent vme and energies between involved band
    allocate (vme_ex_band(npointstotal,3,nband_ex,nband_ex))
    allocate (ek(npointstotal,nband_ex))
    vme_ex_band=0.0d0
    ek=0.0d0
    
    !write(*,*) (nband_index(j), j=1,nband_ex)
    !pause
    !!$OMP PRIVATE(vme_der,shift_vector), &
    !!$OMP PRIVATE(berry_eigen1,berry_eigen2,berry_eigen)
    write(*,*) 'mapping the BZ...'
    !Brillouin zone sampling - parallelization
    do ibz=1,npointstotal    
      write(*,*) 'point', ibz, npointstotal
      rkx=rkxvector(ibz)
		  rky=rkyvector(ibz)
      rkz=rkzvector(ibz)		

      !get matrices in the \alpha, \alpha' basis (orbitals,k)    		
      call get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel, &
           hkernel,hderkernel,akernel)
      call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
           hk_ev,e,vme)
      
      !saving eigenvalues and optical matrix elements at this k point for the bandlist
		  do i=1,nband_ex 
		    ii=nband_index(i)
	      ek(ibz,i)=e(ii)
		    do nj=1,3			
		      do j=1,nband_ex
			      jj=nband_index(j)		  
		  	    vme_ex_band(ibz,nj,i,j)=vme(nj,ii,jj)
		      end do
        end do
      end do
    end do
    write(*,*) 'Optical matrix elements have been evaluated'
     
  end subroutine get_ome
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   	        
  ! This routine evaluates the alpha,alpha' matrices
  subroutine get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel, &
             hkernel,hderkernel,akernel)
    implicit none 
    
    integer norb
    integer ialpha
    integer ialphap
    integer iRp
    integer nj

    dimension skernel(norb,norb)
    dimension hkernel(norb,norb)
    dimension sderkernel(3,norb,norb)
    dimension hderkernel(3,norb,norb)
    dimension akernel(3,norb,norb)

    dimension hderhop(3,nR,norb,norb)
    dimension sderhop(3,nR,norb,norb)

    real(8) Rx,Ry,Rz
    real(8) rkx,rky,rkz
    

    complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
    complex*16 phase,factor
    complex*16 hderhop,sderhop,rhop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  	  	  
    hkernel=0.0d0
    hderkernel=0.0d0
    skernel=0.0d0
    sderkernel=0.0d0
    akernel=0.0d0

    hderhop=0.0d0
    sderhop=0.0d0

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
  subroutine get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
             hk_ev,e,vme)        
    implicit none 

    integer :: norb
    integer :: ialpha
    integer :: ialphap
    integer :: iRp
    integer :: nj
    integer :: i,j,ii,jj,nn,nnp
    
    dimension skernel(norb,norb)
    dimension hkernel(norb,norb)
    dimension sderkernel(3,norb,norb)
    dimension hderkernel(3,norb,norb)
    dimension akernel(3,norb,norb)

    dimension vjseudoa(3,norb,norb)
    dimension vjseudob(3,norb,norb)
    dimension e(norb)
    dimension hk_ev(norb,norb)
    dimension vme(3,norb,norb)
    
    real*8 e
    complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
    complex*16 hk_ev,vjseudoa,vjseudob,vme
    complex*16 amu,amup,aux1,factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
    !diagonalization
    e=0.0d0
    call diagoz(norb,e,hkernel) 
    hk_ev(:,:)=hkernel(ii,jj)
    !Multiply the eigenvectors (C_{\alpha=1 n}(k_0), C_{\alpha=1 n}(k_0), ...) 
    !by a phase phi_n(k_0), this is new eigenvectors are 
    !exp(-i*phi_n(k_0))*(C_{\alpha=1 n}(k_0), C_{\alpha=1 n}(k_0), ...) 
    !
    !Note: right now this phase is used to give a locally smooth phase in the BZ
    call phase_eigvec_nk(norb,hk_ev)
    vme=0.0d0
    vjseudoa=0.0d0
    vjseudob=0.0d0
    
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
        !call cpu_time(time2)
        !write(*,*) 'k-sampling time',norb,time2-time1,'s'
        !pause   
        do nj=1,3
          vme(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+vjseudob(nj,nn,nnp)
          vme(nj,nnp,nn)=conjg(vme(nj,nn,nnp))
          vjseudoa(nj,nnp,nn)=conjg(vjseudoa(nj,nn,nnp))
          vjseudob(nj,nnp,nn)=conjg(vjseudob(nj,nn,nnp))
        end do      
      end do
    end do
         
  end subroutine get_vme_eigen_ome
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine phase_eigvec_nk(norb,hk_ev)
  implicit none
  
  integer norb
  integer i,j,ii

  dimension hk_ev(norb,norb)
  
  real*8 :: arg
  complex*16 :: aux1,hk_ev,factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

  end subroutine
end module ome