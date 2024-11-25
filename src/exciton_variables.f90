module exciton_variables
  implicit none
    !allocatable e_ex(:)
	  !allocatable fk_ex(:,:),fk_ex_der(:,:,:),fk_ex_sq_k(:,:)
 

  contains
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

  
end module exciton_variables