module parser_optics_xatu_dim
  use constants_math
  use parser_wannier90_tb, &
  only:material_name,R !variables
  use parser_input_file, &
  only:ndim,npointstotal_sq,nv_ex,nc_ex, &
       iflag_xatu !variables
  implicit none

  integer npointstotal
  integer norb_ex
  integer norb_ex_cut
  integer norb_ex_band
  integer nband_ex
  integer naux
  
  real(8) G,vcell
  real(8) rkxvector,rkyvector,rkzvector
  real(8) auxr1
  real(8) e_ex
  complex*16 fk_ex
  
  dimension G(3,3)

  allocatable rkxvector(:)
  allocatable rkyvector(:)
  allocatable rkzvector(:)
  allocatable fk_ex(:,:)
  allocatable e_ex(:)
  allocatable auxr1(:)
	  
  contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !TO BE ADAPTED FOR 3D	
  ! Here we define some BZ variables, either by reading the output
  ! of Xatu or not
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_optics_xatu_dim()
	    implicit none  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    !get reciprocal lattice vectors (to be adapted for 3D!)
      call get_reciprocal_vectors()
      
      !!! CUTOFF: avoid reading all excitons. LESS THAN norb_ex
	    norb_ex_cut=100

      !get band and grid dimensions: from XATU-output or opticx-input
      if (iflag_xatu .eqv. .true.) then
        call get_exciton_dim()
        nband_ex=nc_ex+nv_ex
      else
        norb_ex_band=nv_ex*nc_ex
        norb_ex=norb_ex_band*npointstotal
        nband_ex=nc_ex+nv_ex
        !number of total k-points
        if (ndim==1) npointstotal=npointstotal_sq
        if (ndim==2) npointstotal=npointstotal_sq**2
        if (ndim==3) npointstotal=npointstotal_sq**3
      end if

      !allocate grid and exciton arrays
	    allocate (rkxvector(npointstotal))
	    allocate (rkyvector(npointstotal))
      allocate (rkzvector(npointstotal))
	    allocate (e_ex(norb_ex_cut))
	    allocate (fk_ex(norb_ex,norb_ex_cut))

      !fill exciton and other arrays: from XATU-output or opticx-input
      if (iflag_xatu .eqv. .true.) then
        call get_exciton_data()
      else
        call get_grid()
        !exciton variables set to zero if XATU interface is not requested
        fk_ex=0.0d0
        e_ex=0.0d0
      end if
      write(*,*) "grid, bands, and exciton quantities loaded"
	  end subroutine get_optics_xatu_dim   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This routine needs to be adapted for 3D
    subroutine get_reciprocal_vectors()
      implicit none
      real(8) cx,cy,cz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    G=0.0d0
      G(1,1)=2.0d0*pi*(-R(1,1)*R(2,2)+R(1,2)*R(2,1))**(-1.0d0) &
             *(-R(2,2))
      G(1,2)=2.0d0*pi*(-R(1,1)*R(2,2)+R(1,2)*R(2,1))**(-1.0d0) &
             *R(2,1)            
      G(2,1)=2.0d0*pi*(-R(2,1)*R(1,2)+R(2,2)*R(1,1))**(-1.0d0) &
             *(-R(1,2))
      G(2,2)=2.0d0*pi*(-R(2,1)*R(1,2)+R(2,2)*R(1,1))**(-1.0d0) &
             *R(1,1)    
	  		
      call crossproduct(R(1,1),R(1,2),0.0d0,R(2,1), &
      R(2,2),0.0d0,cx,cy,cz)         
      vcell=sqrt(cx**2+cy**2+cz**2)
    
    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
	  subroutine get_exciton_dim()
      implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		  
      open(10,file='./xatu_files_input/'//trim(material_name)//'.eigval')

	    read(10,*) nv_ex,nc_ex
	    read(10,*) npointstotal_sq
	    read(10,*) naux
    
      close(10)	  
	    npointstotal=npointstotal_sq**2
	    norb_ex_band=int(naux/npointstotal)
	    norb_ex=norb_ex_band*npointstotal
	  	  
	    !return
	  end subroutine get_exciton_dim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
    subroutine get_exciton_data()
	    implicit none
      integer j,nkaka
      integer ib,ibz,jind  
	    real(8) auxr1

	    dimension auxr1(2*norb_ex)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		  
	    !get energies
	    open(10,file='./xatu_files_input/'//trim(material_name)//'.eigval')
	    read(10,*) 
	    read(10,*)
	    read(10,*) nkaka,(e_ex(j), j=1,norb_ex_cut)
      close(10)

	    open(10,file='./xatu_files_input/'//trim(material_name)//'.eigvec')	  	  
	    read(10,*) 

	    !reading k-mesh
	    do ibz=1,npointstotal
		    read(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz)
		    do ib=1,norb_ex_band-1
		      read(10,*) 
	      end do
	    end do

	    !reading exciton-wf  
      do ibz=1,norb_ex_cut
		    write(*,*) 'reading exciton wf',ibz,norb_ex    
		    read(10,*) (auxr1(j),j=1,2*norb_ex)
        do j=1,norb_ex    
          jind=2*j-1
          fk_ex(j,ibz)=complex(auxr1(jind),auxr1(jind+1))
        end do
	    end do
	    close(10)	
	    e_ex=e_ex/27.211385d0
      rkxvector=rkxvector*0.52917721067121d0 
      rkyvector=rkyvector*0.52917721067121d0 
	    rkzvector=rkzvector*0.52917721067121d0 
	  end subroutine get_exciton_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !TO BE ADAPTED FOR 3D
    subroutine get_grid()
      implicit none

      integer :: npointsdir,k,i1,i2
      real(8) :: step,r1,r2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      step=1.0d0/dble(npointsdir-1)
  
      k=1
      do i1=1,npointsdir
		    r1=-0.5d0+dble(i1-1)*step
        do i2=1,npointsdir
          r2=-0.5d0+dble(i2-1)*step
          rkxvector(k)=r1*G(1,1)+r2*G(2,1) 
          rkyvector(k)=r1*G(1,2)+r2*G(2,2) 
          rkzvector(k)=0.0d0
          k=k+1       
        end do
      end do  
                           
    end subroutine get_grid

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_ex_index(nv,nv_ex,nc_ex,iright,ibz,i_ex,ic,iv)
      implicit none
      integer nv,nv_ex,nc_ex,ibz,i_ex,ic,iv
      integer iright 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !get band indeces (respect the fermi level) from A_cv index
	    if (iright.eq.1) then
	      iv=(nv-nv_ex)+i_ex-int((i_ex-1)/nv_ex)*nv_ex-nv
	      ic=(nv+1)+int((i_ex-1)/nv_ex)-nv	
	    end if
        !get A_cv index from band indeces
	    if (iright.ne.1) then
        i_ex=nc_ex*nv_ex*(ibz-1)+nv_ex*(ic-1)+iv
	    end if	  	  
    end subroutine get_ex_index

end module parser_optics_xatu_dim

