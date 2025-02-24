module parser_input_file
  implicit none
  private
  public :: material_name_in,filename_input, &
  iflag_xatu_text,response_text
  public :: iflag_xatu
  public :: ndim,nv_ex,nc_ex,nf_v,nf_c,npointstotal_sq
  public :: e1,e2
  public :: get_input_file

  character(len=100) :: material_name_in
  character(len=100) :: filename_input
  character(len=100) :: iflag_xatu_text
  character(len=100) :: response_text

  logical :: iflag_xatu
  
  integer :: ndim,nv_ex,nc_ex
  integer :: nf_v,nf_c
  integer :: npointstotal_sq
  
  real(8) :: e1,e2

  contains 
    subroutine get_input_file()
      implicit none
    
      call get_command_argument(1, filename_input)
      open(10,file='./bin/'//trim(filename_input)//'')
      read(10,*) !# Periodic dimensions
      read(10,*) ndim
      read(10,*) !# Wannier90_filename
      read(10,*) material_name_in
      read(10,*) !# Xatu_interface
      read(10,*) iflag_xatu_text
      if (iflag_xatu_text == 'true') then
        iflag_xatu= .true.
        nv_ex=0
        nc_ex=0
        npointstotal_sq=0
        read(10,*) !# Nfermi
        read(10,*) nf_v
        read(10,*) nf_c
        read(10,*) !# Response
        read(10,*) response_text
        read(10,*) !# Energy_range
        read(10,*) e1,e2
        close(10)
      else if (iflag_xatu_text  == 'false') then
        iflag_xatu= .false.
        read(10,*) !# Bloch_info
        read(10,*) nv_ex
        read(10,*) nc_ex
        read(10,*) npointstotal_sq
        read(10,*) !# Nfermi
        read(10,*) nf_v
        read(10,*) nf_c
        read(10,*) !# Response
        read(10,*) response_text
        read(10,*) !# Energy_range
        read(10,*) e1,e2
        close(10)
      else
        write(*,*) 'Error: Invalid value in the input file. Expected "true" or "false".'
        stop
      end if


    end subroutine get_input_file
end module parser_input_file