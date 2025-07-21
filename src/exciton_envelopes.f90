module exciton_envelopes
  use constants_math
  use parser_optics_xatu_dim, &
  only:G,npointstotal,rkxvector,rkyvector,rkzvector, &
       norb_ex,norb_ex_cut,norb_ex_band,fk_ex
  implicit none

  allocatable fk_ex_der(:,:,:)
  complex*16 fk_ex_der

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_fk_ex_der_k()
      implicit none
      
      allocatable :: fk_ex_aux(:,:),fk_ex_der_aux(:,:,:)
      
	  !here
	  dimension rk1vector(npointstotal),rk2vector(npointstotal)
      
      integer :: ibz
      integer :: nn,j,j_aux
      integer :: nj

      real*8:: rk1vector,rk2vector
      real*8 :: xc,yc,xp_bz,yp_bz,xc_bz,yc_bz
      real*8 :: rk1,rk2,rkxp,rkyp
      real*8 :: rkxp_bz,rkyp_bz,rk1_bz,rk2_bz
      
      complex*16 :: fk_ex_aux,fk_ex_der_aux
	  complex*16 :: fk_ex_k_interp,fk_ex_k_interp_back,fk_ex_k_interp_for	
	  complex*16 :: f_grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      


      !auxiliary arrays within this subroutine
      allocate(fk_ex_aux(npointstotal,norb_ex_cut))
      allocate(fk_ex_der_aux(3,npointstotal,norb_ex_cut))    
      fk_ex_aux=0.0d0
	  fk_ex_der_aux=0.0d0
           
      !get_k_kc gives the crystal coordinates of the k-point.
      !we first fill a vector with all kc points. Note that is is trivial (-.5 to 0.5 ...)
      !but we can test the routine with this simple tasl
      do ibz=1,npointstotal
	    call get_k_kc(G,rkxvector(ibz),rkyvector(ibz),rk1vector(ibz), &
		              rk2vector(ibz),xp_bz,yp_bz,xc_bz,yc_bz)
	  end do

      !We evaluate the k-derivative of the exciton envelope function by finite differences with interpolation
      write(*,*) '   Evaluating exciton envelope function derivative with respect to k...'
	  !Sum over e-h pairs
	  do j=1,norb_ex_band	  
	    !write(*,*) j    
	    !fill a matrix with vector with A_vc(npointstotal,norb_ex)		  
	    do ibz=1,npointstotal
		  j_aux=(ibz-1)*norb_ex_band+j
		  do nn=1,norb_ex_cut         
            fk_ex_aux(ibz,nn)=fk_ex(j_aux,nn)	  
          end do
	    end do	  
      
		!Sum over exciton states
	    do nn=1,norb_ex_cut 
          do ibz=1,npointstotal		
		    !exciton_wf
            rkxp=rkxvector(ibz)
		    rkyp=rkyvector(ibz)		  
            call get_k_kc(G,rkxp,rkyp,rk1,rk2,rkxp_bz,rkyp_bz,rk1_bz,rk2_bz)	
		    !null derivative at borders
            if (abs((abs(rk1)-0.5d0)).lt.0.00001d0 .or. abs((abs(rk2)-0.5d0)).lt.0.00001d0) then 
              fk_ex_der_aux(1,ibz,nn)=0.0d0
		      fk_ex_der_aux(2,ibz,nn)=0.0d0
		      !write(*,*) nn,ibz
		    else
              !!pause
              rkxp=rkxvector(ibz)+dk
		      rkyp=rkyvector(ibz)
              call get_k_kc(G,rkxp,rkyp,rk1,rk2,rkxp_bz,rkyp_bz,rk1_bz,rk2_bz)	
	          call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector, &
	               norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,nn,fk_ex_k_interp_for)
              rkxp=rkxvector(ibz)-dk
		      rkyp=rkyvector(ibz)
              call get_k_kc(G,rkxp,rkyp,rk1,rk2,rkxp_bz,rkyp_bz,rk1_bz,rk2_bz)	
	          call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector, &
	                 norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,nn,fk_ex_k_interp_back)
		    	 
		      fk_ex_der_aux(1,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)	 
		  	
              rkxp=rkxvector(ibz)
		      rkyp=rkyvector(ibz)+dk
              call get_k_kc(G,rkxp,rkyp,rk1,rk2,rkxp_bz,rkyp_bz,rk1_bz,rk2_bz)
	          call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector, &
	               norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,nn,fk_ex_k_interp_for)
		    	 
              rkxp=rkxvector(ibz)
		      rkyp=rkyvector(ibz)-dk
              call get_k_kc(G,rkxp,rkyp,rk1,rk2,rkxp_bz,rkyp_bz,rk1_bz,rk2_bz)
	          call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector, &
	               norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,nn,fk_ex_k_interp_back)
		           fk_ex_der_aux(2,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)	 		    
		    end if
          end do			
        end do
	
		!Save e-h pair wf
		do ibz=1,npointstotal
		  j_aux=(ibz-1)*norb_ex_band+j
		  do nn=1,norb_ex_cut
		    do nj=1,3			
			  fk_ex_der(nj,j_aux,nn)=fk_ex_der_aux(nj,ibz,nn)
			    !if (nj.eq.1 .and. nn.eq.1) then
			        !write(*,*) ibz,j_aux,fk_ex_der(nj,j_aux,nn)
			        !end if
			end do
		  end do
		end do

      end do

    end subroutine get_fk_ex_der_k




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    subroutine get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector, &
	             norb_ex,norb_ex_cut,fk_ex,rk1,rk2,nn,fk_ex_k_interp)
      implicit none 
      
	  !in/out
	  integer :: npointstotal,norb_ex,norb_ex_cut,nn
	  
	  dimension rk1vector(npointstotal),rk2vector(npointstotal)
	  dimension fk_ex(npointstotal,norb_ex_cut)
	 
	  real*8 :: rk1vector,rk2vector
	  real*8 :: rk1,rk2
	  complex*16 :: fk_ex,fk_ex_k_interp
      
	  !here
	  integer :: nside,nblock
	  integer :: ibz_q11,ibz_q21,ibz_q12,ibz_q22
      real*8 :: slice
	  real*8 :: x1,x2,y1,y2,x,y
	  complex*16 f_q11,f_q12,f_q21,f_q22
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  	   
	  !identify the block in which the point is
	  nside=nint(sqrt(dble(npointstotal)))
	  slice=1.0d0/dble(nside-1)
	  nblock=int((rk1+0.5d0)/slice)+1
 
	  !get the corners
	  ibz_q11=nblock+int((rk2+0.5d0)/slice)*nside
	  ibz_q21=ibz_q11+1
	  ibz_q12=ibz_q11+nside
	  ibz_q22=ibz_q21+nside  
	  !write(*,*) rk1_bz,rk2_bz
	  !write(*,*) ibz_q11,ibz_q21,ibz_q12,ibz_q22
	  !write(*,*) 'heeey'
	  !pause  
	  if (ibz_q11.gt.npointstotal .or. &
	      ibz_q21.gt.npointstotal .or. &
		  ibz_q12.gt.npointstotal .or. &
		  ibz_q22.gt.npointstotal) then
	      write(*,*) 'interpolation went wrong'	
	      fk_ex_k_interp=0.0d0
          pause
	  end if
	  if (ibz_q11.lt.1 .or. &
	      ibz_q21.lt.1 .or. &
		  ibz_q12.lt.1 .or. &
		  ibz_q22.lt.1) then
	      write(*,*) 'interpolation went wrong'
	      fk_ex_k_interp=0.0d0
		  pause
	  end if
	  
	  !interpolate
      x1=rk1vector(ibz_q11)
	  x2=rk1vector(ibz_q21)
	  y1=rk2vector(ibz_q11)
	  y2=rk2vector(ibz_q12)
	  x=rk1
	  y=rk2
	  
	  f_q11=fk_ex(ibz_q11,nn)
	  f_q21=fk_ex(ibz_q21,nn)	
	  f_q12=fk_ex(ibz_q12,nn)
	  f_q22=fk_ex(ibz_q22,nn)

	  fk_ex_k_interp=1.0d0/((x2-x1)*(y2-y1))* &
	                 (f_q11*(x2-x)*(y2-y)+f_q21*(x-x1)*(y2-y) &
	                  +f_q12*(x2-x)*(y-y1)+f_q22*(x-x1)*(y-y1))
  
	end subroutine get_fk_ex_k_interp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
	  subroutine get_k_kc(G,xp,yp,xc,yc,xp_bz,yp_bz,xc_bz,yc_bz)
	    implicit none
	  
	    dimension G(3,3),G_half(3,3)

        real*8 :: G,G_half
        real*8 :: xp,yp,xc,yc,xp_bz,yp_bz,xc_bz,yc_bz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
        G_half=G
	  
	    xc=(xp*G_half(2,2)-G_half(2,1)*yp)/(G_half(1,1)*G_half(2,2)-G_half(2,1)*G_half(1,2))
	    yc=(G_half(1,1)*yp-xp*G_half(1,2))/(G_half(1,1)*G_half(2,2)-G_half(2,1)*G_half(1,2))
	    !write(*,*) xp,yp,xp_bz,yp_bz,xc,yc
	    xc_bz=xc-dble(int(xc/0.5d0))
	    yc_bz=yc-dble(int(yc/0.5d0))
	    !end if		
	    xp_bz=xc_bz*G(1,1)+yc_bz*G(2,1)
	    yp_bz=xc_bz*G(1,2)+yc_bz*G(2,2)
	  end subroutine get_k_kc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
end module exciton_envelopes

