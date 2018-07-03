!this program is devoted to the computation of the correlation function
!of tensors. it is not expected for it to be heavy and will receive all
!configurations as input, outputting the average and estimated error.

program tmunu_corr
use types_params
implicit none
character(len=1024) :: filename,list_of_files
complex(dp), allocatable :: Tmunu(:,:,:), corr(:,:,:,:,:)
integer :: record_len,x,nw,number_of_files,unit_list_files,x1,x2,mu,nu,rho,sigma,t,k,s1,s2,file_num,t1,t2
!load parameters (lattice size and lattice file name)
call readargs() 

!allocate vector that will hold tmunu
allocate(tmunu(4,4,nx*ny*nz*nt))
allocate(corr(4,4,4,4,0:nt-1))

corr = CMPLX(0.0_dp,0.0_dp)
open(newunit=unit_list_files,file=list_of_files)

!open file with list of data filenames
inquire(iolength=record_len) Tmunu(1,1,1)

!for each file
do file_num=1,number_of_files
    !read filename
    read(unit_list_files,*) filename
    !Open file
    open(newunit=nw,file=filename,form="unformatted",access='direct',recl=record_len)
    !We load the file to memory
    do x=0,nx*ny*nz*nt-1
        do mu=1,4
            do nu=1,4
                print *, nu + 4*(mu-1) + 16*x, filename, list_of_files
                read(nw,rec=nu + 4*(mu-1) + 16*x) Tmunu(nu,mu,x)
            end do
        end do
    end do
    close(nw)

    !Computes the correlation. Typically, it is displayed as a function of time only.
    !Thus, we average over space.
    !Also, since this is a first run, I am not worried about errors. Later I need to implement
    !a more robust method, eg binned jacknife
    do t=1,nt
        do k=0,nt-1
            do s1=1,nx*ny*nz
                do s2=1,nx*ny*nz
                    x1 = s1 + t1*nx*ny*nz
                    x2 = s2 + (t+k -nt*((k+t)/nt))*nx*ny*nz
                    do nu=1,4
                        do mu=1,4
                            do sigma=1,4
                                do rho=1,4
                                    t1=x1/(nx*ny*nz)
                                    t2=x2/(nx*ny*nz)
                                    corr(rho,sigma,mu,nu,t2-t1) = corr(rho,sigma,mu,nu,abs(t1-t2)) + Tmunu(mu,nu,x1)*Tmunu(rho,sigma,x2)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    close(unit_list_files)
end do

!Normalization: 
corr = corr/(real(nt*number_of_files)*real(nx*ny*nz)**2)

!Finally, we save the file
do nu=1,4
    do mu=1,4
        do sigma=1,4
            do rho=1,4
                write(filename,'("Corr",4I1.1,".dat")') mu,nu,rho,sigma
                open(newunit=nw,file=filename)
                do t=0,nt-1
                    write(nw,*) t, corr(mu,nu,rho,sigma,t)
                end do
                close(nw)
            end do
        end do
    end do
end do

deallocate(tmunu,corr)

contains
   subroutine readargs()
      integer, parameter :: minnumberparameters = 5
      character(len=50) :: argnx,argny,argnz,argnt,arg_number_of_files
      if(command_argument_count() .lt. minnumberparameters) then
         print*, "it is mandatory to pass 5 arguments in the format"
         print*, "nx ny nz nt path/to/lattice/file"
         print*, "exiting now"
         call exit(1)
      else
      
         call get_command_argument(1,argnx)
         call get_command_argument(2,argny)
         call get_command_argument(3,argnz)
         call get_command_argument(4,argnt)
         call get_command_argument(5,list_of_files)
         call get_command_argument(6,arg_number_of_files)

         read(argnx,*) nx
         read(argny,*) ny
         read(argnz,*) nz
         read(argnt,*) nt
         read(arg_number_of_files,*) number_of_files
      end if
   end subroutine
end program