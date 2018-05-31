!Willian Matioli Serenone
!Institution: Universidade de Sao Paulo
!                  Instituto de Fisica de Sao Carlos
!e-mail: willian.serenone@usp.br
!######################################################################


!==============================
!MODULE: types_params
!Contains definitions of the types used as well as global parameters.
!Must be used on all modules (except ziggurat, which is an external lib).
!==============================
module types_params
implicit none

integer,parameter :: dp=kind(0.d0) !We use double precision by default.
                                   !change to 0.0 to use simple precision
real(dp), parameter :: pi = acos(-1.0_dp)
real(dp), parameter :: two_pi = pi*2.0_dp
real(dp), parameter :: tol = 1.0e-15_dp !Numerical error tolerance used on our checks
integer :: nx,ny,nz,nt ! Lattice dimensions
real(dp) :: beta ! Beta parameter of the simulation

type SU2
   real(dp) :: a(4)
end type SU2

type SU3
   real(dp) :: re(3,3)
   real(dp) :: im(3,3)
end type


type(SU3),allocatable :: U(:,:) !Gauge links 

end module types_params

!==============================
!About the parametrizations of the groups.
!For SU(2), we require a1**2 + a2**2 + a3**2 + a4**2 = 1
!With U = dot(a,sigma) + a4 and sigma are the Pauli matrices

!An alternative is to interpret
!x = a4 + i*a3
!y = a2 + i*a1
!And the group element
!    (  x  y  )
!    ( -y* x* )
!
!For the SU(3) element, we consider
!    (    u   )
!U = (    v   )
!    ( u* ^ v*)
!u and v are a 3-vector of complex numbers
!I store them as 6 real numbers
!Thus, u_i = dcmplx( u(2*i-1), u(2*i) )
!We require that
!u u* = 1
!v v* = 1
!dot(u,v) = 0
!Also, we store the result of u* ^ v* in w
!==============================

!==============================
!About the indexing of the gauge arrays on the U variable
!First index is direction. It ranges from 1 to 8.
!  -Even indexes is pointing forward in direction 2j, where
!   j is the direction index
!  -Odd indexes is pointing backwards in direction 2j-1

!Second index is the linearized lattice position i.
!It relates to the four coordinates x,y,z,t as
!
! i=
