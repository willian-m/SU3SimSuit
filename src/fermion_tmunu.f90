
program fermions_prop
use types_params
use lattice, only : lattice_file, init_lattice
!use FoX_sax
use xml_parser, only : read_xml
use cuda_fermion_prop
implicit none
character(len=1024) :: xml_input_path
type(SU3),allocatable :: prop(:,:)

!Get xml input file path, parse it and load the lattice.

call GET_COMMAND_ARGUMENT(1,xml_input_path)
call GET_COMMAND_ARGUMENT(2,lattice_file)

call read_xml(xml_input_path)
!For test, we turn off the gauge field. Thus, the prop should be  prop to 1/(p^2 + m^2) + lat. artifacts
hot_start = .false.
call init_lattice('')

allocate(prop(0:12,0:12*nx*ny*nz*nt))

open(unit=1,file="fermion_prop")
call compute_prop(prop)
do i=0,11
   do j=0,12*nx*ny*nz*nt-1
      print *,
   end do
end do

end program
