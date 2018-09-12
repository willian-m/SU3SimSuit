
program fermions_prop
use types_params
use lattice, only : lattice_file, init_lattice
!use FoX_sax
use xml_parser, only : read_xml
use cuda_fermion_prop
implicit none
character(len=1024) :: xml_input_path

!Get xml input file path, parse it and load the lattice.

call GET_COMMAND_ARGUMENT(1,xml_input_path)
call GET_COMMAND_ARGUMENT(2,lattice_file)

call read_xml(xml_input_path)
call init_lattice(lattice_file)

call compute_prop 

end program
