module parser_input_file
  implicit none
  private
  public :: material_name_in
  public :: filename_input
  public :: iflag_xatu_text
  public :: iflag_ome_sp_text
  public :: iflag_ome_ex_text
  public :: response_text
  public :: iflag_xatu
  public :: iflag_ome_sp
  public :: iflag_ome_ex
  public :: ndim,nf,npointstotal_sq
  public :: e1,e2,eta,nw
  public :: get_input_file
  public :: nband_index
  public :: norb_ex_cut
  public :: read_line_numbers_int !subroutine

  character(len=100) :: material_name_in
  character(len=100) :: filename_input
  character(len=100) :: iflag_xatu_text
  character(len=100) :: iflag_ome_sp_text
  character(len=100) :: iflag_ome_ex_text
  character(len=100) :: response_text

  logical :: iflag_xatu
  logical :: iflag_ome_sp
  logical :: iflag_ome_ex

  integer :: nk1
  integer :: ndim
  integer :: nf
  integer :: npointstotal_sq
  integer :: nband_index_aux 
  integer :: norb_ex_cut
  integer :: j !to read stuff
  integer :: nband_index
  integer :: nw
  real(8) :: e1,e2,eta

  dimension :: nband_index_aux(100) !auxiliary array to save band indeces
  allocatable :: nband_index(:)

  contains 
    subroutine get_input_file()
      implicit none
      allocatable :: narray(:) 
      integer :: i,narray,num_values
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) '1. Entering parser_input_file'
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call get_command_argument(1, filename_input)
      !open(10,file='./bin/'//trim(filename_input)//'')
      open(10,file=trim(filename_input))
      read(10,*) !# Periodic dimensions
      read(10,*) ndim
      read(10,*) !# Wannier90_filename
      read(10,*) material_name_in
      read(10,*) !# Xatu_interface
      read(10,*) iflag_xatu_text
      if (iflag_xatu_text == 'true') then
        iflag_xatu= .true.
        npointstotal_sq=0
        read(10,*) !# Exciton_cutoff
        read(10,*) norb_ex_cut
        read(10,*) !# Nfermi
        read(10,*) nf
        read(10,*) !# OME_sp
        read(10,*) iflag_ome_sp_text
        read(10,*) !# OME_ex
        read(10,*) iflag_ome_ex_text    
        read(10,*) !# Response
        read(10,*) response_text
        read(10,*) !# Energy_variables
        read(10,*) e1,e2,eta,nw
        close(10)
      else if (iflag_xatu_text  == 'false') then
        iflag_xatu= .false.
        read(10,*) !# Bandlist
        !reads a line of numbers and stores them in an allocatable array 'narray'
        call read_line_numbers_int(narray,num_values) 
        !write(*,*) (narray(j),j=1,num_values)
        read(10,*) !# Ncells
        read(10,*) npointstotal_sq
        read(10,*) !# Nfermi
        read(10,*) nf
        read(10,*) !# OME_sp
        read(10,*) iflag_ome_sp_text
        read(10,*) !# Response
        read(10,*) response_text
        read(10,*) !# Energy_variables
        read(10,*) e1,e2,eta,nw
        close(10)
        !fill the nband_index array after knowing the value of num_values
        !it will be used in parser_optics_xatu_dim.f90
        allocate (nband_index(num_values)) 
        nband_index(:)=narray(:)
      else
        write(*,*) 'Error: Invalid value in the input file. Expected "true" or "false".'
        stop
      end if
      
      !declare flags from text strings
      if (iflag_ome_sp_text == 'true') then
        iflag_ome_sp= .true.
      else
        iflag_ome_sp= .false.
      end if
      if (iflag_ome_ex_text == 'true') then
        iflag_ome_ex= .true.
      else
        iflag_ome_ex= .false.
      end if
      
      write(*,*) '   Input file has been read'
    end subroutine get_input_file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This routine reads a line of numbers into an array
    subroutine read_line_numbers_int(narray,num_values)
    implicit none
    allocatable :: narray(:)

    integer :: ios,ncount,i
    integer :: narray
    integer :: num_values
    integer :: temp_num 
    integer :: istart,iposition
    character(len=1000) :: line
 
    !Read the line of text
    read(10,'(A)', iostat=ios) line
    if (ios /= 0) then
      print *, 'Error reading file'
    stop
    end if

    ncount=0
    do i=1,len_trim(line)
      if (line(i:i) == ' ') ncount=ncount+1
    end do
    ncount=ncount+1  ! One more than the number of spaces

    !Allocate the array based on the number of values
    allocate(narray(ncount))

    !Reset the number of values counter
    num_values=0
    istart=1
    !Now, sequentially extract numbers from the line
    do i=1,ncount
      !Find the next number in the line
      read(line(istart:),*,iostat=ios) temp_num
      if (ios == 0) then
        num_values=num_values+1
        narray(num_values)=temp_num
        !Move the starting  position to the next number
        iposition=scan(line(istart:),' ')
        if (iposition>0) then
          istart=istart+iposition
        end if
      end if
    end do

    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module parser_input_file