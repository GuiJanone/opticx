module parser_xatu_ex
  use parser_wannier90_tb 
  implicit none
  
  integer npointstotal_sq
  integer npointstotal
  integer norb_ex
  integer norb_ex_cut
  integer norb_ex_band
  integer nv_ex
  integer nc_ex
  integer naux
  
  real(8) rkxvector,rkyvector,rkzvector
  real(8) auxr1
  real(8) e_ex
  complex*16 fk_ex
 
  allocatable rkxvector(:)
  allocatable rkyvector(:)
  allocatable rkzvector(:)
  allocatable fk_ex(:,:)
  allocatable e_ex(:)
  allocatable auxr1(:)
	  
  contains
    
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
	  	  
	  !!! CUTOFF: avoid reading all excitons. LESS THAN norb_ex
	  norb_ex_cut=100
	  allocate (e_ex(norb_ex_cut))
	  allocate (rkxvector(npointstotal))
	  allocate (rkyvector(npointstotal))
	  allocate (fk_ex(norb_ex,norb_ex_cut))
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
	  !write(*,*) e_ex(1)
	  !pause
           !exciton_wf_file_name='GeS_states_nk80.out'
	  open(10,file='./xatu_files_input/'//trim(material_name)//'.eigvec')	  	  
	  read(10,*) 
           !norb_ex=npointstotal
	       !reading k-mesh
	  do ibz=1,npointstotal
		read(10,*) rkxvector(ibz),rkyvector(ibz) !,rkzvector(ibz)
		do ib=1,norb_ex_band-1
		  read(10,*) 
	    end do
	  end do
	       !do ib=1,norb_ex_band
	       !reading exciton-wf  
      do ibz=1,norb_ex_cut
	       !do ibz=1,norb_ex
		write(*,*) 'reading exciton wf',ibz,norb_ex    
		   !do ib=1,norb_ex_band
		read(10,*) (auxr1(j),j=1,2*norb_ex)
        do j=1,norb_ex    
          jind=2*j-1
          fk_ex(j,ibz)=complex(auxr1(jind),auxr1(jind+1))
        end do
	  end do
           !end do
	  close(10)	
	  e_ex=e_ex/27.211385d0
      rkxvector=rkxvector*0.52917721067121d0 
      rkyvector=rkyvector*0.52917721067121d0 
	  
	end subroutine get_exciton_data

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


end module parser_xatu_ex