program wilson_loop_correlation
    use types_params
    use lattice, only : lattice_file,init_lattice,pos
    use objects, only : wilson_loop
    use math, only : SU3_Tr
    use IO, only: read_lattice


    implicit none
    integer :: dim1,dir1,dim2,dir2,x,y,z,t,block_index,t_index,s,n,m,Step_Size
    integer :: n_slices,block_size,nw,b1,b2,t1,t2,i,n_sweeps,n_subsweeps
    complex(dp), allocatable :: average_space(:,:), corr_func(:,:), WL(:,:), corr_func_avrg(:), corr_func_error(:)
    integer, allocatable :: counter(:)
    type(SU3) :: loop
    character(1024) :: FilePath, aux

    !Read arguments and initiate lattice
    call read_args
    lattice_file=''
    call init_lattice(lattice_file) !Just to allocate the lattice

    !allocate vector that will store the average_space of Wilson loops
    block_size=nt/n_slices
    allocate(average_space(block_size-1,n_slices),WL(block_size-1,n_slices))
    allocate(corr_func(nt-1,n_sweeps),counter(nt-1))
    allocate(corr_func_avrg(nt-1),corr_func_error(nt-1))
    !We will need to compute the wilson_loops over each sweep and subsweep
    WL=CMPLX(0.0_dp,0.0_dp)
    corr_func=CMPLX(0.0_dp,0.0_dp)
    do n=1,n_sweeps
        do m=Step_Size,n_subsweeps,Step_Size
            !First, we build the filename
            lattice_file= TRIM(FilePath)
            write(aux,"(I7.7,'Subsweep',I7.7,'.dat')") n,m
            lattice_file = TRIM(lattice_file) // trim(aux)
            !Now, we load the lattice configuration to memory
            call read_lattice(lattice_file)
            !And we compute the average_space of Wilson loops
            call wilson_loop_computer
            !Now that we have the Wilson_Loop (averaged over each time-lice), we
            !average over the subsweep
            WL = WL + average_space
        end do
        WL = WL/(n_subsweeps/Step_Size)
        !We are ready now to compute the correlation function associated to this sweep
        call corr_func_computer(n)
    end do

    !Now we can compute the expectation value of the correlation and estimate its error`
    corr_func_avrg=CMPLX(0.0_dp,0.0_dp)
    corr_func_error=CMPLX(0.0_dp,0.0_dp)
    do n=1,n_sweeps
        do m=1,nt-1
            corr_func_avrg(m)= corr_func_avrg(m) + corr_func(m,n)
            corr_func_error(m)= corr_func_error(m) + corr_func(m,n)**2
        end do
    end do
    corr_func_avrg = corr_func_avrg/n_sweeps
    corr_func_error = corr_func_error/n_sweeps

    do m=1,nt-1
        corr_func_error(m) = cmplx(sqrt((real(corr_func_error(m)) - real(corr_func_error(m))**2)/n_sweeps),&
        sqrt((imag(corr_func_error(m)) - imag(corr_func_error(m))**2)/n_sweeps))

    end do

    !Now we write to the file the output
    open(newunit=nw, status='replace', file='Wilson_Lopp_Correlation.out')
    write(nw,'(3A,4X,23A,2X,23A,4X,23A,2X,23A)') '#t ','Re[C(t)]','Im[C(t)]','Re[Delta C(t)]','Im[Delta C(t)]'
    do m=1,nt-1
        write(nw,"(I3.3,4X,ES23.16,2X,ES23.16,4X,ES23.16,2X,ES23.16)") m,corr_func_avrg(m),corr_func_error(m)
    end do

    contains

    subroutine corr_func_computer(sweep)
        integer,intent(in) :: sweep
        counter = 0 !counter(s) stores how many sums was performed when computing corr_function(s)
    
        !Iterate over blocks
        do b1=1,n_slices
            do b2=1,n_slices
                if (b2 .ne. b1) then !If we are not in the same block, we iterate over each combination of time-slice of the blocks
                    do t1=1,block_size-1
                        do t2=1,block_size-1
                            if (b2 .gt.  b1) then !s stores the distance between time-slices. This depends if b1 > b2 or otherwise
                                s = t2 - t1 + block_size*(b2-b1)
                                corr_func(s,sweep) = corr_func(s,sweep) + WL(t2,b2)*WL(t1,b1)
                                counter(s) = counter(s) + 1
                            else if (b1 .gt. b2) then
                                s = t1 - t2 + block_size*(b1-b2)
                                corr_func(s,sweep) = corr_func(s,sweep) + WL(t2,b2)*WL(t1,b1)
                                counter(s) = counter(s) + 1
                            end if
                        end do
                    end do
                end if
            end do
        end do
    
        !Normalize with counter
        do s=1,nt-1
            if (counter(s) .gt. 0) then
                corr_func(s,sweep) = corr_func(s,sweep)/counter(s)
            end if
        end do
    end subroutine corr_func_computer

    subroutine wilson_loop_computer
        !Computes the average_space of the Wilson loop on each time-slice inside the block
        average_space = CMPLX(0.0_dp,0.0_dp)
        do t=1,nt
            t_index = mod(t,block_size)
            if (t_index .ne. 0) then
                block_index = t/block_size + 1
                do z=1,nz
                    do y=1,ny
                        do x=1,nx
                            i = pos(x,y,z,t)
                            call wilson_loop(dir1,dim1,dir2,dim2,i,loop)
                            average_space(t_index,block_index) = average_space(t_index,block_index) + SU3_Tr(loop)
                        end do
                    end do
                end do
            end if
        end do
        average_space = average_space/(nx*ny*nz)
    end subroutine wilson_loop_computer

    subroutine read_args
        integer, parameter :: min_num_par = 13 !Change this if input changes
        character(len=50) :: arg_nx,arg_ny,arg_nz,arg_nt,arg_nslices
        character(len=50) :: arg_nsweeps,arg_nsubsweeps,arg_dim1,arg_dir1,arg_dim2,arg_dir2,arg_StepSize
     
        if (COMMAND_ARGUMENT_COUNT() .lt. min_num_par) then
           print*,"It is mandatory to pass at least 10 parameters as input in the format:"
           print*,"nx ny nz nt N_Slices N_sweeps N_subsweeps Step_Size dimension_1, direction_1, dimension_2, direction_2, lattice/Data/linksxxxyyyzzztttbeta9.99nslices999Sweep"
           print*,"   -nx,ny,nz,nt: integers representing the lattice dimensions nx, ny, nz and nt"
           print*,"   -N_Slices: Number of time slices that will be frozen during one-level execution"
           print*,"   -dimension_1,2: The size of the sides of the Wilson loop"
           print*,"   -directions_1,2: Determine the plane on which the Wilson loop lives"
           print*,"   -lattice/Data/FilePath.dat: Optional path to a lattice configuration that will be used as initial condition, overridingf the Hot/Cold choice."
           print*,"Exiting now"
           call EXIT(1)
        else
           call GET_COMMAND_ARGUMENT(1,arg_nx)
           call GET_COMMAND_ARGUMENT(2,arg_ny)
           call GET_COMMAND_ARGUMENT(3,arg_nz)
           call GET_COMMAND_ARGUMENT(4,arg_nt)
           call GET_COMMAND_ARGUMENT(5,arg_nslices)
           call GET_COMMAND_ARGUMENT(6,arg_nsweeps)
           call GET_COMMAND_ARGUMENT(7,arg_nsubsweeps)
           call GET_COMMAND_ARGUMENT(8,arg_StepSize)
           call GET_COMMAND_ARGUMENT(9,arg_dim1)
           call GET_COMMAND_ARGUMENT(10,arg_dir1)
           call GET_COMMAND_ARGUMENT(11,arg_dim2)
           call GET_COMMAND_ARGUMENT(12,arg_dir2)
           call GET_COMMAND_ARGUMENT(13,FilePath)
           read(arg_nx,*) nx
           read(arg_ny,*) ny
           read(arg_nz,*) nz
           read(arg_nt,*) nt
           read(arg_nslices,*) n_slices
           read(arg_nsweeps,*) n_sweeps
           read(arg_nsubsweeps,*) n_subsweeps
           read(arg_StepSize,*) Step_Size
           read(arg_dim1,*) dim1
           read(arg_dim2,*) dim2
           read(arg_dir1,*) dir1
           read(arg_dir2,*) dir2         

        end if
        end subroutine read_args



end program wilson_loop_correlation
