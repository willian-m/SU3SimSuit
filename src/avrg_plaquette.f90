!Willian Matioli Serenone
!Institution: Universidade de Sao Paulo
!                  Instituto de Fisica de Sao Carlos
!e-mail: willian.serenone@usp.br
!######################################################################

!==============================
!PROGRAM: avrg_plaquette
!Read a lattice configuration and computes the value of the average plaquette on it
!==============================

!==============================
!Input list:
!  -4 integers representing the lattice dimensions nx, ny, nz and nt
!  -A string with the path to the file
!==============================

program avrg_plaquette
use types_params
use IO, only : read_lattice
use lattice, only : init_lattice, lattice_file, hot_start
use objects, only : plaquette
use math, only : SU3_ReTr
implicit none
real(dp) :: avrg_plaq
type(SU3) :: plaq
integer :: x,y,z,t,d1,d2

!Load parameters (lattice size and lattice file name)
call read_args() 

!Allocates lattice and loads file
hot_start = .false. !Not relevant, but set it anyway so the variable is initialized
call init_lattice(lattice_file)

avrg_plaq = 0.0_dp
!Computes each one of the six plaquettes for each point, and sum their trace to compute the action
do x=0,nx*ny*nz*nt-1
   do d1=1,4
      do d2=d1+1,4
         call plaquette(2*d1,2*d2,x,plaq)
         avrg_plaq = avrg_plaq + SU3_ReTr(plaq)
      end do
   end do
end do

!Normalizes
avrg_plaq = avrg_plaq/(6*nx*ny*nz*nt)

!prints to stdin
write(6,"(ES23.15E3)") avrg_plaq

contains
   subroutine read_args()
      integer, parameter :: minNumberParameters = 5
      character(len=50) :: argNx,argNy,argNz,argNt
      
      if(COMMAND_ARGUMENT_COUNT() .lt. minNumberParameters) then
         print*, "It is mandatory to pass 6 arguments in the format"
         print*, "nx ny nz nt beta path/to/Lattice/File"
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
   end subroutine read_args
end program
