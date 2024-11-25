program main_program
  use constants
  use parser_wannier90_tb
  use parser_xatu_ex
  use ome
  implicit none

  integer :: j
  character(len=100) :: material_name_in
 
  

  !Reading input file
  open(10,file='./bin/input.in')
  read(10,*)
  ! Read the first word (or line) from the file into the variable material_name_in
  !read(10,*) j
  !material_name_in=trim(material_name_in)
  ! Close the file
  read(10,*) material_name_in
  close(10)
  
  !get wannier90 info
  call wannier90_get(material_name_in)
  !get xatu files
  call get_exciton_dim()
  call get_exciton_data()
  !evaluate optical matrix elements
  call get_ome()

  write(*,*) 'It seems that we are ok'
  pause
  !write(*,*) R(1,1)
  pause


end program main_program