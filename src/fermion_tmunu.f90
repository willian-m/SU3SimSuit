program fermion_tmunu
    use types_params
    use lattice, only : lattice_file, init_lattice, hot_start
    use fermion_prop

    implicit none
    character(len=1024) :: xml_input_path
    complex(dp), allocatable :: propagator(:,:)
    integer :: i,j


    call GET_COMMAND_ARGUMENT(1,xml_input_path)
    call GET_COMMAND_ARGUMENT(2,lattice_file)
    
    hot_start = .false.
    nx=4
    ny=4
    nz=4
    nt=4

    allocate(propagator(0:12*nx*ny*nz*nt-1,0:11))
    call init_lattice(lattice_file)
    call compute_prop(propagator)
    
    open(unit=1,file="fermion_prop.dat")
    do i=0,12*nx*ny*nz*nt-1
        do j=0,11
            write(1,*) i,j,propagator(i,j)
        end do
    end do
    close(1)
end program