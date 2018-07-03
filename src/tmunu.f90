!This program is devoted to the computation of the energy-momentum tensor
!on a single Monte-Carlo configuration
!Tipically, you will run it once for each lattice configuration

program tmunu_corr
use types_params
use lattice, only : lattice_file, init_lattice
use objects, only: CalcTmunu => Tmunu

implicit none
character(len=1024) :: filename = ''
complex(dp), allocatable, dimension(:,:) :: Tmunu
integer :: record_len,x,nw,mu,nu
!Load parameters (lattice size and lattice file name)
call readArgs() 

!Allocates lattice and loads file
call init_lattice(lattice_file)

!Now we compute all components of the tensor

!Computes T0i over the entire lattice
allocate(Tmunu(4,4))

!We compute and output the data
inquire(iolength=record_len) Tmunu(1,1)
open(newunit=nw,file="Tmunu.dat",form="unformatted",access='direct',recl=record_len)
do x=0,nx*ny*nz*nt-1
   call CalcTmunu(x,Tmunu)
   do mu=1,4
      do nu=1,4
         write(nw,rec=nu + 4*(mu-1) + 16*x) Tmunu(nu,mu)
      end do
   end do
end do
close(nw)

contains
   subroutine readArgs()
      integer, parameter :: minNumberParameters = 4
      character(len=50) :: argNx,argNy,argNz,argNt
      if(COMMAND_ARGUMENT_COUNT() .lt. minNumberParameters) then
         print*, "It is mandatory to pass 4 arguments in the format"
         print*, "nx ny nz nt path/to/Lattice/File"
         print*, "Exiting now"
         call EXIT(1)
      else
      
         call GET_COMMAND_ARGUMENT(1,argNx)
         call GET_COMMAND_ARGUMENT(2,argNy)
         call GET_COMMAND_ARGUMENT(3,argNz)
         call GET_COMMAND_ARGUMENT(4,argNt)
         call GET_COMMAND_ARGUMENT(5,lattice_file)

         read(argNx,*) nx
         read(argNy,*) ny
         read(argNz,*) nz
         read(argNt,*) nt

      end if
   end subroutine
end program