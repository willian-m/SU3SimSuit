!Willian Matioli Serenone
!Institution: Universidade de Sao Paulo
!             Instituto de Fisica de Sao Carlos
!e-mail: willian.serenone@usp.br
!######################################################################


!==============================
!MODULE: heat_bath
!Contains everything needed to execute one heat_bath hit in each lattice site
!==============================
module heat_bath
use types_params
use lattice, only : increment_table
use math, only : SU3mult,SU2mult,detSU2_like,subgroup,embed_in_SU3,SU3_dagger,SU3projector
use objects, only : compute_staple
implicit none

!==============================
!List of subroutines
!  -heat_bath_hit: Executes a heat bath hit on the site x and direction d
!==============================

!==============================
!List of public variables

!List private variables
type(SU3) :: staple !Staple of the links to be updated
type(SU3) :: W !Mean field variable
!==============================

private
public heat_bath_method

contains

!==============================
!Executes the heat bath hit on the entire lattice
subroutine heat_bath_method
integer :: d,i

do i=0,nx*ny*nz*nt-1
   do d=2,8,2
      call heat_bath_hit(d,i)
   end do
end do
end subroutine heat_bath_method
!==============================

!==============================
!Executes a heat bath hit on the site x and direction d
subroutine heat_bath_hit(d,y)
integer :: d,y,a,b
type(SU2) :: X,V,Rw
type(SU3) :: R,Aux
real(dp) :: det


!1) Compute Staple
call compute_staple(d,y,staple)

!For each one the 3 SU(2) subgroups of SU(3)
do a=1,2
   do b=a+1,3
      !2) Close staple with link U
      call SU3mult(U(d,y),staple,W)

      !3) Factor out only the relevant terms in W
      call subgroup(W,a,b,Rw)

      !4)Projects Rw in the SU(2)
      det = sqrt(detSU2_like(Rw))
      V%a = Rw%a/det

      !5)Select a new random SU(2) element X following Boltzmann distribution
         !This boils down to choosing a(4) following an exponential distribution
         !And then randomly choosing a direction for the next 3 vectors.
         !I implemented two methods of choosing a(4). Creutz method is found in
         !Creutz book 'Quarks, gluons and lattices'. GL method is found in
         !Gattringer & Lang book 'Quantum Chromodynamics on the Lattice'.
         !Use one method or another.
      call random_a4_GL(det,X%a(4))
         !Now that we have a(4), we pick the other three by randomly picking a vector
         !On the surface of an sphere of radius sqrt(1-a(4)**2)
      call rand_pnt_sphere_marsaglia(sqrt(1.0_dp - X%a(4)**2),X%a(1:3))

      !6)Takes the dagger of V
      V%a(1:3) = -V%a(1:3)
 
      !7)The new Rw' is Rw' = X*V^\dagger
      call SU2mult(X,V,Rw)

      !8)Embeds Rw' in SU(3)
      call embed_in_SU3(Rw,a,b,R)
      !9)Update the link
      call SU3mult(R,U(d,y),Aux)
      !10)Look for possible rounding error in the updated link and fix it
      call SU3projector(Aux)
     
      !11)Saves the link for next iteration
      U(d,y) = Aux
      call SU3_dagger(U(d,y),U(d-1,y+increment_table(y,d)))
   end do
end do

end subroutine heat_bath_hit
!==============================

!==============================
!Randomly pick a(4) following the right distribution - Gattringer method
subroutine random_a4_Creutz(k,a4)
real(dp), intent(in) :: k
real(dp), intent(inout) :: a4
real(dp) :: z, r

r = exp(beta*k*2.0_dp)

call random_number(z) !Fortran's random number generator
z = z*(r-1.0_dp/r) + 1.0_dp/r
if ( z .ge. sqrt(1.0_dp - (log(z)/beta/k)**2) ) then
      a4 = log(z)/beta/k
end if

!random_a4_GL = dlog(z)/beta/k
end subroutine random_a4_Creutz
!==============================

!==============================
!Randomly pick a(4) following the right distribution - GL method
subroutine random_a4_GL(det,a4)
use ziggurat
real(dp), intent(in) :: det
real(dp), intent(inout) :: a4
real(dp) :: r(4), lambda2
integer :: i
logical :: not_accepted
not_accepted = .true.

do while ( not_accepted )
call random_number(r)
do i=1,4
   r(i) = 1.0_dp - r(i)
end do
lambda2 = -(log(r(1)) + log(r(3))*(cos(two_pi*r(2)))**2)/(det*beta*2.0_dp)
if ( r(4)**2 .le. 1.0_dp - lambda2) then
   a4 = 1.0_dp - lambda2*2.0_dp
   not_accepted = .false.
end if
end do
end subroutine random_a4_GL
!==============================

!==============================
!Randomly pick a vector in the surface of an sphere of radius R
subroutine rand_pnt_sphere_marsaglia(R,A)
real(dp),intent(in) :: R
real(dp),intent(out) :: A(3)
logical :: accepted
real(dp) :: z(2) 
real*8,dimension(3) :: temp

accepted = .false.
do while (.not. accepted)
   call random_number(z)
   z(1) = 2.0_dp*z(1)-1.0_dp
   z(2) = 2.0_dp*z(2)-1.0_dp
   if ( z(1)**2 + z(2)**2 .lt. 1.0_dp) then
      temp(1) = 2.0_dp*z(1)*sqrt(1.0_dp-z(1)**2-z(2)**2)
      temp(2) = 2.0_dp*z(2)*sqrt(1.0_dp-z(1)**2-z(2)**2)
      temp(3) = 1.0_dp - 2.0_dp*(z(1)**2+z(2)**2)
      accepted = .true.
   end if
end do

A(1) = temp(1)*R
A(2) = temp(2)*R
A(3) = temp(3)*R

end subroutine rand_pnt_sphere_marsaglia

end module heat_bath
