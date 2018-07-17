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
use math, only : SU3mult,SU3_dagger,SU3_Tr
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
public plaquette,compute_staple,Tmunu

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

   staple%a= cmplx(0.0_dp, 0.0_dp)
   do d2=2,8,2
      if (d2 .ne. d) then

         x1 = x + increment_table(x,d)
         x2 = x1 + increment_table(x1,d2)
         x3 = x + increment_table(x,d2)
         x4 = x1 + increment_table(x1,d2-1)
         x5 = x + increment_table(x,d2-1)

         call SU3mult(U(d2,x1),U(d-1,x2),aux1)
         call SU3mult(aux1,U(d2-1,x3),aux2)

         staple%a = staple%a + aux2%a

         call SU3mult(U(d2-1,x1),U(d-1,x4),aux1)
         call SU3mult(aux1,U(d2,x5),aux2)

         staple%a = staple%a + aux2%a

      end if
   end do
end subroutine compute_staple
!==============================

!==============================
!Compute a clover = Q_{mu nu}(x)
subroutine clover(mu, nu, x, Q)
    type(SU3), intent(out) :: Q
    integer, intent(in) :: mu,nu,x
    integer :: negative_mu,negative_nu
    type(SU3) :: aux
    
    negative_mu = mu+2*mod(mu,2)-1
    negative_nu = nu+2*mod(nu,2)-1

    call plaquette(mu,nu,x,Q)
    call plaquette(nu,negative_mu,x,aux)
    Q%a = Q%a + aux%a
    call plaquette(negative_mu,negative_nu,x,aux)
    Q%a = Q%a + aux%a
    call plaquette(negative_nu,mu,x,aux)
    Q%a = Q%a + aux%a

end subroutine clover
!==============================

!==============================
!Compute tensor Fmunu
subroutine Fmunu(mu, nu, x, F)
    type(SU3), intent(out) :: F
    integer, intent(in) :: mu,nu,x
    type(SU3) :: aux
    !From here on, it does not make sense to consider negative directions, 
    !thus we check all directions are positive

    if(mod(mu,2) .ne. 0 .or. mod(nu,2) .ne. 0) then
        print *, "ERROR! Computing F_{mu nu} for negative directions. Aborting!"
        call EXIT(-2)
    else
        call clover(mu,nu,x,F)
    end if

    !Now we take the dagger
    call SU3_dagger(F,aux)

    !Compute a^2F_{mu nu}
    F%a = (F%a - aux%a)/cmplx(0.0_dp,8.0_dp)
    
end subroutine Fmunu
!==============================

!==============================
!Compute tensor Tmunu
subroutine Tmunu(x, T)
    complex(dp), intent(out) :: T(4,4)
    integer, intent(in) :: x
    type(SU3) :: F(4,4), aux
    integer :: rho,sigma,mu,nu
    complex(dp) :: diag_term

    !Compute the upper triangle of F_{\mu \nu} (without diagonals, since it is 0 and we do not use them)
    !And fill the lower triangle with the negative of the result
    do rho=2,8,2
        do sigma=rho+2,8,2
            call Fmunu(rho,sigma,x,F(sigma/2,rho/2))
            F(rho/2,sigma/2)%a = -F(sigma/2,rho/2)%a
        end do
    end do

    !Compute the upper triangle tensor components (including diagonals)
    T=CMPLX(0.0_dp,0.0_dp)
    do mu=1,4
        do nu=mu,4
            do sigma=1,4
                if (mu .ne. sigma .and. nu .ne. sigma) then !We skip diagonals terms of F_{mu nu}, since it is zero
                    call SU3mult(F(mu,sigma),F(nu,sigma),aux)
                    T(nu,mu) = T(nu,mu) + 2.0_dp*SU3_Tr(aux)
                end if
            end do
            T(mu,nu) = T(nu,mu) !Populate the lower triangle components, using its symmetry
        end do
    end do

    !Compute the global term that will be subtracted from diagonal
    diag_term = CMPLX(0.0_dp,0.0_dp)
    do rho=1,4
        do sigma=rho+1,4
            call SU3mult(F(sigma,rho),F(sigma,rho),aux)
            diag_term = diag_term + 2.0_dp*SU3_Tr(aux)
        end do
    end do

    !Subtract it from the diagonal terms
    do mu=1,4
        T(mu,mu) = T(mu,mu) - diag_term/4.0_dp
    end do
    
end subroutine Tmunu
!==============================

end module objects