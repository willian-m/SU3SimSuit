!Willian Matioli Serenone
!Institution: Universidade de Sao Paulo
!             Instituto de Fisica de Sao Carlos
!e-mail: willian.serenone@usp.br
!######################################################################


!==============================
!MODULE: objects
!Contains lattice related objects, such as plaquettes and staples
!==============================
module objects
use types_params
use lattice, only : increment_table
use math, only : SU3mult
implicit none

!==============================
!List of subroutines
!  -plaquette: Computes a single plaquette
!  -compute_staples: Computes the staples surrounding a link
!  -mean_field: Uses a staple set to calculate a 'mean field'.
!   Usefull for heat-bath simulations
!==============================

!==============================
!List of public variables

!List private variables
!==============================

private
public plaquette,compute_staple

contains

!==============================
!Returns plaquette U_{mu,nu}(x)
subroutine plaquette(dir1,dir2,x,plaq)
   integer, intent(in) :: dir1,dir2,x
   type(SU3),intent(out) :: plaq
   integer :: x1,x2,x3
   type(SU3) :: U_temp1,U_temp2

   x1 = x + increment_table(x,dir1)
   x2 = x1 + increment_table(x1,dir2)
   x3 = x + increment_table(x,dir2)
      

!
!     dir1
!      ^   x1            x2
!      |    ------->------
!      |    |            |
!      |    |            |
!      |    |            |
!      |    ^            v 
!      |    |            |
!      |    |            |
!      |    |            |
!      |    -------<------
!      |   x            x3
!      |
!      ----------------------------> dir2

   call SU3mult( U(dir1,x), U(dir2,x1) , U_temp1 )
   call SU3mult( U_temp1, U( dir1+2*mod(dir1,2)-1, x2 ) ,U_temp2 )
   call SU3mult( U_temp2, U( dir2+2*mod(dir2,2)-1, x3 ) , plaq )
      
end subroutine plaquette
!==============================

!==============================
!Returns staples surrounding link U_{mu}(x)
!NOTE: It only makes sense dir even.
!WARNING: For performance reasons, we do not check if dir is even
!         You are free to pass an odd value for it, but the result will 
!         probably garbage. Know what you are doing!
subroutine compute_staple(d,x,staple)
   type(SU3), intent(out) :: staple
   integer, intent(in) :: d,x
   type(SU3) :: aux1,aux2
   integer :: d2,x1,x2,x3,x4,x5

   staple%re=0.0_dp
   staple%im=0.0_dp
   do d2=2,8,2
      if (d2 .ne. d) then

         x1 = x + increment_table(x,d)
         x2 = x1 + increment_table(x1,d2)
         x3 = x + increment_table(x,d2)
         x4 = x1 + increment_table(x1,d2-1)
         x5 = x + increment_table(x,d2-1)

         call SU3mult(U(d2,x1),U(d-1,x2),aux1)
         call SU3mult(aux1,U(d2-1,x3),aux2)

         staple%re = staple%re + aux2%re
         staple%im = staple%im + aux2%im

         call SU3mult(U(d2-1,x1),U(d-1,x4),aux1)
         call SU3mult(aux1,U(d2,x5),aux2)

         staple%re = staple%re + aux2%re
         staple%im = staple%im + aux2%im

      end if
   end do
end subroutine compute_staple
!==============================


end module objects
