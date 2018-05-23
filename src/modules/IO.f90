!Willian Matioli Serenone
!Institution: Universidade de Sao Paulo
!                  Instituto de Fisica de Sao Carlos
!e-mail: willian.serenone@usp.br
!######################################################################

!=====THE LATTICE MODULE======
!This module stores global variables and declares the lattice

!-----------------------------------------------------
!=====VARIABLES AND TYPES=====
!-----------------------------------------------------

!     -lastWrite: keeps track of when was the last time the lattice
!     configuration was recorded in a file. We index the configuration
!     by number of times that the algoritm has sweeped over the lattice.

!-----------------------------------------------------
!=====SUBROUTINES AND FUNCTIONS======
!-----------------------------------------------------
!
!==============================
!MODULE: IO
!Contains read and write operations from disk
!==============================
module IO
use types_params
implicit none

private
public write_lattice,read_lattice
contains
      
subroutine write_lattice(sweepNum)
integer,intent(in) :: sweepNum
character(len=50) :: fileName
   write(filename,"('links',4I3.3,'beta',F4.2,'Sweep',I9.9,'.dat')") nx,ny,nz,nt,beta,sweepNum
   open(unit=1,status='replace',file=filename,form='unformatted')
   write(1) U 
   close(1)
end subroutine

subroutine read_lattice(filename)
character(len=1024), intent(in) :: filename
   open(unit=10,status='old',file=filename,form='unformatted')
   read(10) U
   close(10)
end subroutine read_lattice

end module IO
