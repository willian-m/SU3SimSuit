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
use heat_bath, only : heat_bath_method,accepted,total
use IO, only : write_lattice
use xml_parser, only : read_xml, nmc, therm, rec_step
use objects, only : S, wilson_action
!==============================

implicit none
!==============================
!List of subroutines
!  -read_args: read the program arguments
!==============================

!==============================
!List of variables
integer :: clock !Used to set ziggurat seed
integer :: n
character(len=1024) :: xml_input_path
!==============================

!read input parameters
call GET_COMMAND_ARGUMENT(1,xml_input_path)
call read_xml(xml_input_path)
!call read_args(nx,ny,nz,nt,hot_start,nmc,beta,rec_step,lattice_file)

!Initialize ziggurat random number generator
clock = 123
!call system_clock(count=clock) !Uncomment to different seeds in each call
call zigset(clock) !Set seed for ziggurat random number generator
call init_random_seed !Set seed for fortran random number generator 

!Initialize the lattice according to user request
call init_lattice(lattice_file)
!S=wilson_action()
print *,0!,",", S
call write_lattice(0)

!Start counting the efficiency of Metropolis algorithm in the overrelaxation method
accepted=0
total=0
do n=1,nmc
   call heat_bath_method
   if (mod(n,rec_step) .eq. 0 .and. n .gt. therm) then
     call write_lattice(n)
   end if
   print *, n
end do

print *, "Finished. Overrelaxation method accepted ",dble(accepted)*100/dble(total),"percent."
print *, "My work is finished. You may want to run the avrg_plaquette or the tmunu script now."
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
      print*, "nx ny nz nt beta H/C NMC recordStep [lattice/Data/FilePath.dat]"
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
      call GET_COMMAND_ARGUMENT(5,arg_beta)
      call GET_COMMAND_ARGUMENT(6,arg_hot_start)
      call GET_COMMAND_ARGUMENT(7,arg_nmc)
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

!==init_random_seed and lcg are adapted from gnu fortran manual
!==The modification allows to always generate the same sequence
!Source link
!==https://gcc.gnu.org/onlinedocs/gcc-6.4.0/gfortran/RANDOM_005fSEED.html
!==Accessed in 7th, July 2018.
   subroutine init_random_seed()
      use iso_fortran_env, only: int64
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid
      integer(int64) :: t

      call random_seed(size = n)
      allocate(seed(n))
      t=123456
      do i = 1, n
        seed(i) = lcg(t)
        t=t+123
      end do
      call random_seed(put=seed)
  end subroutine init_random_seed

  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    use iso_fortran_env, only: int64
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg

end program


