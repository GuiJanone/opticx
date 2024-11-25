module ome
  use parser_wannier90_tb
  use parser_xatu_ex
  implicit none

  allocatable vme(:,:,:)
  
  complex(16) vme
  
  
  contains
  
  subroutine get_ome()
    implicit none

    integer ibz
    real(8) rkxp,rkyp,rkzp
    allocate(vme(3,norb,norb)) 
    vme=0.0d0
    write(*,*) 'mapping the BZ...'

	  do ibz=1,npointstotal     
      write(*,*) 'point', ibz, npointstotal
      rkxp=rkxvector(ibz)
		  rkyp=rkyvector(ibz)
      rkzp=rkzvector(ibz)		
      !get matrices in the \alpha, \alpha' basis (orbitals,k)    		
      call get_vme_kernels_ome()
      !call get_vme_kernels_ome(rkxp,rkyp,rkzp,nR,nRvec,norb,R,hkernel,skernel,shop, &
      !hhop,rhop,sderhop,hderhop,sderkernel,hderkernel,akernel,pgaugekernel)
    end do

    
  end subroutine get_ome


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   	        
  ! This routine evaluates the alpha,alpha' matrices
  subroutine get_vme_kernels_ome()
    implicit none 

    dimension hkernel(norb,norb),skernel(norb,norb)
    dimension hderhop(3,nR,norb,norb)
    dimension sderhop(3,nR,norb,norb)
    
    dimension sderkernel(3,norb,norb)
    dimension hderkernel(3,norb,norb)
    dimension akernel(3,norb,norb)
    
    integer ialpha
    integer ialphap
    integer iRp
    integer nj
    
    real(8) Rx,Ry,Rz
    real(8) rkx,rky,rkz

    complex(16) phase,factor
    complex(16) hderhop,sderhop,rhop
    complex(16) hderkernel,sderkernel    
    complex(16) hkernel,skernel
    complex(16) akernel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  	  	  

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
end module ome