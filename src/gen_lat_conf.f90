!Willian Matioli Serenone
!Institution: Universidade de Sao Paulo
!                  Instituto de Fisica de Sao Carlos
!e-mail: willian.serenone@usp.br
!######################################################################

!==============================
!PROGRAM: gen_lat_conf
!Generates pure-gauge SU(3) configurations
!==============================

!==============================
!Input list:
!  -4 integers representing the lattice dimensions nx, ny, nz and nt
!  -A float corresponding to the value of beta
!  -H or C to choose between a "Hot" start or a "Cold" start
!  -A (large) integer to tell how many MC steps will be performed
!  -A number telling the intervall for which each lattice configuration will be stored. Must be different of 0
!  -Optionally, one may give the name of a lattice configuration that will be used as initial condition, overridingf the Hot/Cold choice.
!==============================

program gen_lat_conf

!==============================
!List of modules
use types_params
use ziggurat, only : zigset
use lattice, only : hot_start,lattice_file,init_lattice
use heat_bath, only : heat_bath_method
use IO, only : write_lattice
!==============================

implicit none
!==============================
!List of subroutines
!  -read_args: read the program arguments
!==============================

!==============================
!List of variables
integer :: clock !Used as seed
integer :: seed(2),n,nmc
integer :: rec_step
!==============================

!read input parameters
call read_args(nx,ny,nz,nt,hot_start,nmc,beta,rec_step,lattice_file)

!Initialize ziggurat random number generator
clock = 123
!call system_clock(count=clock) !Uncomment to different seeds in each call
call zigset(clock)
!This is for the use of the native fortran random number generator
!data seed /123456789, 987654321/
!call random_seed(size=2,put=seed)
call random_seed()

print*, "Lattice size: ",nx,"x",ny,"x",nz,"x",nt
print*, "Beta:", beta

!Initialize the lattice according to user request
call init_lattice(lattice_file)
do n=1,nmc
   call heat_bath_method
   print *, dble(n)*100.0_dp/nmc, "% completed."
   if (mod(nmc,n) .eq. 0) then
      call write_lattice(n)
   end if
end do

!==================================================================================
contains !REMEMBER: Make sure var names declared on the functions and 
         !subroutines are local, even if they match the name of var of
         !the main program. However, we can use the variables of the
         !main program without the need to declare them.

!==Subroutine to read input arguments
   subroutine read_args(nx,ny,nz,nt,hot_start,nmc,beta,rec_step,lattice_file)
   integer, intent(out) :: nx,ny,nz,nt,nmc,rec_step
   real(dp),intent(out) :: beta
   logical, intent(out) :: hot_start
   character(len=1024), intent(out) :: lattice_file
   
   integer, parameter :: min_num_par = 8 !Change this if input changes
   character(len=50) :: arg_nx,arg_ny,arg_nz,arg_nt,arg_hot_start,arg_nmc,arg_beta,arg_rec_step

   if (COMMAND_ARGUMENT_COUNT() .lt. min_num_par) then
      print*, "It is mandatory to pass at least 8 parameters as input in the format:"
      print*, "nx ny nz nt beta H/C NMC recordStep [lattice/Data/FilePath.dat"
      print*,"  -nx,ny,nz,nt: integers representing the lattice dimensions nx, ny, nz and nt"
      print*,"  -beta: a float corresponding to the value of beta"
      print*,'  -H/C: choose between a "Hot"(H) start or a "Cold"(C) start'
      print*,"  -NMC: a (large) integer to tell how many MC steps will be performed"
      print*,"  -recordStep: a number telling the intervall for which each lattice configuration will be stored. Must be different of 0"
      print*,"  -lattice/Data/FilePath.dat: Optional path to a lattice configuration that will be used as initial condition, overridingf the Hot/Cold choice."

      print*, "Exiting now"
      call EXIT(1)
   else
      call GET_COMMAND_ARGUMENT(1,arg_nx)
      call GET_COMMAND_ARGUMENT(2,arg_ny)
      call GET_COMMAND_ARGUMENT(3,arg_nz)
      call GET_COMMAND_ARGUMENT(4,arg_nt)
      call GET_COMMAND_ARGUMENT(6,arg_hot_start)
      call GET_COMMAND_ARGUMENT(7,arg_nmc)
      call GET_COMMAND_ARGUMENT(5,arg_beta)
      call GET_COMMAND_ARGUMENT(8,arg_rec_step)
      call GET_COMMAND_ARGUMENT(9,lattice_file)
      read(arg_nx,*) nx
      read(arg_ny,*) ny
      read(arg_nz,*) nz
      read(arg_nt,*) nt
      if (arg_hot_start(1:1) .eq. "H") then
         hot_start = .true.
      else if (arg_hot_start(1:1) .eq. "C") then
         hot_start = .false.
      else
         print*,'WARNING: 6th argument must be "H" or "C"'
         print*,'=======I will proceed assuming a cold lattice start'
         hot_start = .false.
      end if

      read(arg_nmc,*) nmc
      read(arg_beta,*) beta
      read(arg_rec_step,*) rec_step
   end if
   end subroutine read_args

end program


