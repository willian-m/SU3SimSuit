program wilson_loop_correlation
    use types_params
    use lattice, only : lattice_file,init_lattice,pos
    use objects, only : wilson_loop
    use math, only : SU3_Tr


    implicit none
    integer :: dim1,dir1,dim2,dir2,x,y,z,t,block_index,t_index,s
    integer :: n_slices,block_size,nw,b1,b2,t1,t2,i
    complex(dp), allocatable :: average(:,:), corr_func(:)
    integer, allocatable :: counter(:)
    type(SU3) :: loop
    character(1024) :: filename

    !Read arguments and initiate lattice
    call read_args
    call init_lattice(lattice_file)

    print *, "Read ok."
    !Computes the average of the Wilson loop on each time-slice inside the block
    block_size=nt/n_slices
    allocate(average(block_size-1,n_slices))
    average = CMPLX(0.0_dp,0.0_dp)
    do t=1,nt
        t_index = mod(t,block_size)
        print *, "t_index is valid:", t_index .le. block_size-1 
        if (t_index .ne. 0) then
            block_index = t/block_size + 1
            print *, "block_index is valid:", block_index .le. n_slices
            do z=1,nz
                do y=1,ny
                    do x=1,nx
                        i = pos(x,y,z,t)
                        call wilson_loop(dir1,dim1,dir2,dim2,i,loop)
                        print *, "Wilson loop computed"
                        average(t_index,block_index) = average(t_index,block_index) + SU3_Tr(loop)
                        print *, average(t_index,block_index)
                    end do
                end do
            end do
        end if
    end do
    print*, "Ok"
    average = average/(nx*ny*nz)

    print *, "Average ok."
    !We compute the correlation function between these averages.
    allocate(corr_func(nt-1),counter(nt-1))
    corr_func = CMPLX(0.0_dp,0.0_dp)
    counter = 0
    do b1=1,n_slices
        do b2=1,n_slices
            if (b2 .ne. b1) then
                do t1=1,block_size-1
                    do t2=1,block_size-1
                        if (b2 .gt.  b1) then
                            s = t2 - t1 + block_size*(b2-b1)
                            corr_func(s) = corr_func(s) + average(t2,b2)*average(t1,b1)
                            counter(s) = counter(s) + 1
                        else if (b1 .gt. b2) then
                            s = t1 - t2 + block_size*(b1-b2)
                            corr_func(s) = corr_func(s) + average(t2,b2)*average(t1,b1)
                            counter(s) = counter(s) + 1
                        end if
                    end do
                end do
            end if
        end do
    end do

    print *, "Correlation ok."
    do s=1,nt-1
        if (counter(s) .gt. 0) then
            corr_func(s) = corr_func(s)/counter(s)
        end if
    end do

    print *, "Correlation normalized."
    !Now we save the data to disk, so it can be used on a MC averager

    write(filename,"('WilsonCorr_nt',I3.3,'.dat')"),nt 
    open(newunit=nw, status='replace', file=filename, form='unformatted')
    write(nw) corr_func

    print *, "Data saved to disk."
    contains

    subroutine read_args
        integer, parameter :: min_num_par = 10 !Change this if input changes
        character(len=50) :: arg_nx,arg_ny,arg_nz,arg_nt,arg_n_slices,arg_dim1,arg_dir1,arg_dim2,arg_dir2
     
        if (COMMAND_ARGUMENT_COUNT() .lt. min_num_par) then
           print*, "It is mandatory to pass at least 10 parameters as input in the format:"
           print*, "nx ny nz nt N_Slices dimension_1, direction_1, dimension_2, direction_2, lattice/Data/FilePath.dat"
           print*,"   -nx,ny,nz,nt: integers representing the lattice dimensions nx, ny, nz and nt"
           print*,"   -N_Slices: Number of time slices that will be frozen during one-level execution"
           print*,"   -dimension_1,2: The size of the sides of the Wilson loop"
           print*,"   -directions_1,2: Determine the plane on which the Wilson loop lives"
           print*,"   -lattice/Data/FilePath.dat: Optional path to a lattice configuration that will be used as initial condition, overridingf the Hot/Cold choice."
           print*, "Exiting now"
           call EXIT(1)
        else
           call GET_COMMAND_ARGUMENT(1,arg_nx)
           call GET_COMMAND_ARGUMENT(2,arg_ny)
           call GET_COMMAND_ARGUMENT(3,arg_nz)
           call GET_COMMAND_ARGUMENT(4,arg_nt)
           call GET_COMMAND_ARGUMENT(5,arg_n_slices)
           call GET_COMMAND_ARGUMENT(6,arg_dim1)
           call GET_COMMAND_ARGUMENT(7,arg_dir1)
           call GET_COMMAND_ARGUMENT(8,arg_dim2)
           call GET_COMMAND_ARGUMENT(9,arg_dir2)
           call GET_COMMAND_ARGUMENT(10,lattice_file)
           read(arg_nx,*) nx
           read(arg_ny,*) ny
           read(arg_nz,*) nz
           read(arg_nt,*) nt
           read(arg_n_slices,*) n_slices
           read(arg_dim1,*) dim1
           read(arg_dim2,*) dim2
           read(arg_dir1,*) dir1
           read(arg_dir2,*) dir2         

        end if
        end subroutine read_args



end program wilson_loop_correlation
