!this program is devoted to the computation of the correlation function
!of tensors. it is not expected for it to be heavy and will receive all
!configurations as input, outputting the average and estimated error.

include "mkl_dfti.f90"

program tmunu_corr
use types_params
use MKL_DFTI
use statistic, only : average_real, connected_correlation_real, corr_time_int
use xml_parser, only : read_xml,mu,nu,rho,sigma

implicit none
character(len=1024) :: filename,list_of_files
complex(dp), allocatable :: Tmunu(:,:,:), Tmunu_FFT(:), Trhosigma_FFT(:)
real(dp), allocatable :: MC_corr(:,:), avrg(:), observable(:,:)
integer :: record_len,x,nw,number_of_files,unit_list_files,k,file_num,mu_prime,nu_prime,k1,k2,kx,ky,kz,omega
integer, allocatable :: integrated_corr_time(:)
!load parameters (lattice size and lattice file name)
call readargs() 

!allocate vector that will hold tmunu
allocate(Tmunu(4,4,0:nx*ny*nz*nt-1))
allocate(Tmunu_FFT(0:nx*ny*nz*nt-1),Trhosigma_FFT(0:nx*ny*nz*nt-1))
allocate(observable(0:nx*ny*nz*nt-1,number_of_files))
allocate(MC_corr(0:nx*ny*nz*nt-1,number_of_files))
allocate(integrated_corr_time(0:nx*ny*nz*nt-1))
allocate(avrg(0:nx*ny*nz*nt-1))

observable = CMPLX(0.0_dp,0.0_dp)
open(newunit=unit_list_files,file=list_of_files)

!open file with list of data filenames
inquire(iolength=record_len) Tmunu(1,1,1)

!for each file, computes the temporal correlation
do file_num=1,number_of_files
    !read filename
    filename=''
    read(unit_list_files,'(A)') filename
    !Open file
    open(newunit=nw,file=filename,form="unformatted",access='direct',recl=record_len)
    !We load the file to memory
    do x=0,nx*ny*nz*nt-1
        do mu_prime=1,4
            do nu_prime=1,4
                read(nw,rec=nu_prime + 4*(mu_prime-1) + 16*x) Tmunu(nu_prime,mu_prime,x)
            end do
        end do
    end do
    close(nw)

    !Computes the correlation. Typically, it is displayed as a function of time only.
    !Thus, we average over space.
    !Also, since this is a first run, I am not worried about errors. Later I need to implement
    !a more robust method, eg binned jacknife
    call space_time_FFT(Tmunu(nu,mu,:),Tmunu_FFT)
    call space_time_FFT(Tmunu(sigma,rho,:),Trhosigma_FFT)

    !Compute the FFT of the correlation function
    do omega=0,nt-1
        do kz=0,nz-1
            do ky=0,ny-1
                do kx=0,nx-1
                    k1 = 1 + kx + ky*nx + kz*nx*ny + omega*nx*ny*nz
                    k2 = 1 + (nx+1-kx) + (ny+1-ky)*nx + (nz+1-kz)*nx*ny + (nt+1-omega)*nx*ny*nz
                    observable(k1,file_num)=Tmunu_FFT(k1)*Trhosigma_FFT(k2)
                end do
            end do
        end do
    end do

    print *, "Observable for ",trim(filename)," completed."
end do
close(unit_list_files)

!Normalizes data
observable = observable/real(nx*ny*nz*nt)

!Now starts the analysis.
!1) For each point, compute the average and correlation
do k1=0,nx*ny*nz*nt-1
    avrg(k1) = average_real(observable(k1,:))
    call connected_correlation_real(observable(k1,:),MC_corr(k1,:))
    integrated_corr_time(k1) = corr_time_int(MC_corr(k1,:),5)
    !To do: computation of error estimation using binning.
        !bin_size automatically defined by integrated_corr_time 
end do

!Save essential data
write(filename,'("Corr",4I1.1,".dat")') mu,nu,rho,sigma
open(newunit=nw,file=filename)
write(nw, "(A,A,A,A)") "#t, ", "C(t), ","sigma^2/N, ","tau_int, "
do k1=0,nx*ny*nz*nt-1
    
    write(nw,*) k1, avrg(k1), MC_corr(k1,1)/number_of_files, integrated_corr_time(k1)
end do
close(nw)

!do k1=0,nx*ny*nznt-1
!    write(filename,'("MC_Corr",4I1.1,"t=",I7.7,".dat")') mu,nu,rho,sigma,k1
!    open(newunit=nw,file=filename)
!    write(nw,"(A,A)") "#MC_Step, ","C(MC_Step)"
!    do k=1,number_of_files
!       write(nw, *) k, MC_corr(k1,k)
!    end do
!    close(nw)
!end do


contains
    subroutine readargs()      
        integer, parameter :: min_number_parameters = 3
        character(len=1024) :: input_xml
        character(len=50) :: arg_number_of_files
        if(command_argument_count() .ne. min_number_parameters) then
            print*, "it is mandatory to pass 3 arguments in the format"
            print*, "xml/input/file list/of/files/to/average number_of_lines_of_previous files"
            print*, "exiting now"
            call exit(1)
        else
            call get_command_argument(1,input_xml)
            call get_command_argument(2,list_of_files)
            call get_command_argument(3,arg_number_of_files)
            call read_xml(input_xml)
            read(arg_number_of_files,*) number_of_files
        end if
    end subroutine

    subroutine space_time_FFT(input_data,out_DFT)
        complex(dp), intent(in) :: input_data(:)
        complex(dp), intent(out) :: out_DFT(:)
        integer :: stat
        type(DFTI_DESCRIPTOR), pointer :: descHandler
   
        !We need to create a DESCRIPTOR to describe the parameters of the transform
        stat = DftiCreateDescriptor( descHandler, DFTI_DOUBLE, DFTI_COMPLEX, 4, [nx, ny, nz, nt] ) 
        stat = DftiSetValue( descHandler, DFTI_FORWARD_SCALE, 1.0_dp)
        stat = DftiSetValue( descHandler, DFTI_BACKWARD_SCALE, 1.0_dp/two_pi)
        !We do not desire to overwrite the input data
        stat = DftiSetValue( descHandler, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        !Set up the layout of the output
        stat = DftiSetValue( descHandler, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX)
        !Set up data layout
        stat = DftiSetValue( descHandler, DFTI_INPUT_STRIDES, [0,1,nx,nx*ny,nx*ny*nz])
        stat = DftiSetValue( descHandler, DFTI_OUTPUT_STRIDES, [0,1,nx,nx*ny,nx*ny*nz])
        !Commit descriptor!
        stat = DftiCommitDescriptor( descHandler )

        !Compute FFT
        !Need to use slice, in case we are thowing away the last point
        stat = DftiComputeForward( descHandler, input_data, out_DFT)
        stat = DftiFreeDescriptor( descHandler ) !Descriptor, be free!!!
   end subroutine

end program