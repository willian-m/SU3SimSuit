!this program is devoted to the computation of the correlation function
!of tensors. it is not expected for it to be heavy and will receive all
!configurations as input, outputting the average and estimated error.

program tmunu_corr
use types_params
use xml_parser, only : read_xml,mu,nu,rho,sigma
implicit none
character(len=1024) :: filename,list_of_files
complex(dp), allocatable :: Tmunu(:,:,:), corr(:,:,:,:,:)
integer :: record_len,x,nw,number_of_files,unit_list_files,x1,x2,t,k,s1,s2,file_num,t2
!load parameters (lattice size and lattice file name)
call readargs() 

!allocate vector that will hold tmunu
allocate(tmunu(4,4,0:nx*ny*nz*nt-1))
allocate(corr(4,4,4,4,0:nt-1))

corr = CMPLX(0.0_dp,0.0_dp)
open(newunit=unit_list_files,file=list_of_files)

!open file with list of data filenames
inquire(iolength=record_len) Tmunu(1,1,1)

!for each file
do file_num=1,number_of_files
    !read filename
    filename=''
    read(unit_list_files,'(A)') filename
    !Open file
    open(newunit=nw,file=filename,form="unformatted",access='direct',recl=record_len)
    !We load the file to memory
    do x=0,nx*ny*nz*nt-1
        do mu=1,4
            do nu=1,4
                read(nw,rec=nu + 4*(mu-1) + 16*x) Tmunu(nu,mu,x)
            end do
        end do
    end do
    close(nw)

    !Computes the correlation. Typically, it is displayed as a function of time only.
    !Thus, we average over space.
    !Also, since this is a first run, I am not worried about errors. Later I need to implement
    !a more robust method, eg binned jacknife
    do t=0,nt-1
        do k=0,nt-1
            t2 = t + k - nt*((t+k)/nt)
            do s1=0,nx*ny*nz-1
                do s2=0,nx*ny*nz-1
                    x1 = s1 + t*nx*ny*nz
                    x2 = s2 + (t+k -nt*((k+t)/nt))*nx*ny*nz
                    corr(rho,sigma,mu,nu,k) = corr(rho,sigma,mu,nu,k) + Tmunu(mu,nu,x1)*Tmunu(rho,sigma,x2)
                end do
            end do
        end do
    end do
    print *, "Correlation for file ",trim(filename)," completed."
end do
close(unit_list_files)

print *, "Corr computed"
!Normalization: 
corr = corr/(real(nt*number_of_files)*real(nx*ny*nz)**2)

mu=4
nu=1
rho=mu
sigma=nu
corr(rho,sigma,mu,nu,k) = corr(rho,sigma,mu,nu,k) + Tmunu(mu,nu,x1)*Tmunu(rho,sigma,x2)

write(filename,'("Corr",4I1.1,".dat")') mu,nu,rho,sigma
open(newunit=nw,file=filename)
do t=0,nt-1
    write(nw,*) t, real(corr(mu,nu,rho,sigma,t)), imag(corr(mu,nu,rho,sigma,t))
end do
close(nw)

deallocate(tmunu,corr)

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
end program