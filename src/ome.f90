module ome
  use parser_input_file, &
  only:iflag_xatu, &
  iflag_ome_sp_text,iflag_ome_ex_text
  use ome_sp
  use ome_ex
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome()
    implicit none
    integer iflag_norder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !write(*,*) iflag_ome_sp
    !write(*,*) iflag_ome_ex
    write(*,*) '4. Entering ome'
    !here we compute and write into files optical matrix elements
    !single-particle matrix elements
    if (iflag_ome_sp_text == 'linear') then
      iflag_norder=1
      call get_ome_sp(iflag_norder)     
    end if
    if (iflag_ome_sp_text == 'nonlinear') then
      iflag_norder=2
      call get_ome_sp(iflag_norder)  
    end if
    if (iflag_ome_sp_text == 'none' ) then
      write(*,*) '   Optical matrix elements (sp) will be read from file'   
    end if 
    !excitonic matrix elements
    if (iflag_ome_ex_text == 'linear') then
      iflag_norder=1
      call get_ome_ex(iflag_norder)       
    end if
    if (iflag_ome_ex_text == 'nonlinear' ) then
      iflag_norder=2
      call get_ome_ex(iflag_norder)  
    end if
    if (iflag_ome_ex_text == 'none' ) then
      write(*,*) '   Optical matrix elements (ex) will be read from file'   
    end if

  end subroutine get_ome
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
end module ome