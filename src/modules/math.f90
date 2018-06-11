!Willian Matioli Serenone
!Institution: Universidade de Sao Paulo
!             Instituto de Fisica de Sao Carlos
!e-mail: willian.serenone@usp.br
!######################################################################

!==============================
!MODULE: math
!Contains implementation of some math operations to be executed
!==============================
module math
use types_params
use ziggurat, only : rnor
implicit none

!==============================
!List of subroutines
!  -gen_su3_element: generates a random SU(3) group element
!  -gen_su2_element: generates a random SU(2) group element
!  -SU3mult: multiplies two SU3 type matrix
!  -SU3_dagger: computes the dagger of an SU3 element
!==============================

!==============================
!List of public variables


!Private variables
integer :: levi_civita(3,3,3) = reshape((/0,0,0, & !i=1,j=1
                                          0,0,1, & !i=1,j=2
                                          0,-1,0,& !i=1,j=3
                                          0,0,-1,& !i=2,j=1
                                          0,0,0, & !i=2,j=2
                                          1,0,0, & !i=2,j=3
                                          0,1,0, & !i=3,j=1
                                          -1,0,0, &!i=3,j=2
                                          0,0,0/),&!i=3,j=3
                                         (/3,3,3/),order=(/1,3,2/))

!==============================

private
public gen_su3_element,SU3mult,SU2mult,detSU2_like,SU3_dagger,SU3_ReTr,embed_in_SU3,subgroup
public SU3projector,is_not_SU3,invert_3x3_complex
contains

!==============================
!Generates a random SU(3) group element
subroutine gen_su3_element(V)
type(SU3),intent(out) :: V
type(SU2) :: R,S,T
type(SU3) :: R3,S3,T3,aux

!Generates SU(2) matrices R,S,T
call gen_su2_element(R)
call gen_su2_element(S)
call gen_su2_element(T)

!Now we embed these matrices in 3x3 matrices
!parametrized as su3 matrices

call embed_in_SU3(R,1,2,R3)
call embed_in_SU3(S,1,3,S3)
call embed_in_SU3(T,2,3,T3)

!Finally, we set V=R3*S3*T3
call SU3mult(R3,S3,aux)
call SU3mult(aux,T3,V)

end subroutine gen_su3_element
!==============================

!==============================
!Generates a random SU(2) group element
subroutine gen_su2_element(V)
type(SU2),intent(out) :: V
real(dp) :: norm

V%a(1) = cmplx( rnor( ), rnor( ) )
V%a(2) = cmplx( rnor( ), rnor( ) )
Norm = sqrt( real(V%a(1)*conjg(V%a(1)) +  V%a(2)*conjg(V%a(2))) )
V%a = V%a/Norm
end subroutine gen_su2_element
!==============================

!==============================
!Multiplies SU(3) group elements
subroutine SU3mult(A,B,C)
type(SU3), intent(in) :: A,B
type(SU3), intent(out) :: C
real(dp) :: aux(3,3)
integer :: i,j,l

C%a = matmul(A%a,B%a)

end subroutine SU3mult
!==============================

!==============================
!Computes the dagger of an SU3 group element
subroutine SU3_dagger(A,B)
type(SU3), intent(in) :: A
type(SU3), intent(out) :: B
integer :: i,j

do j=1,3
   do i=1,3
      B%a(i,j) = conjg(A%a(j,i))
   end do
end do

end subroutine SU3_dagger
!==============================

!==============================
subroutine subgroup(V,i,j,R)
type(SU3), intent(in) :: V
integer, intent(in) :: i, j
type(SU2), intent(out) :: R

!I assume i < j but the routine does not explicitly check for this. Beware!
R%a(1) = cmplx( real(V%a(i,i) + V%a(j,j) ), aimag(V%a(i,i) - V%a(j,j)) )/2.0_dp
R%a(2) = cmplx( real(V%a(i,j) - V%a(j,i) ), aimag(V%a(i,j) + V%a(j,i)) )/2.0_dp

end subroutine subgroup
!==============================

!==============================
subroutine embed_in_SU3(R,i,j,V)
type(SU3), intent(out) :: V
integer, intent(in) :: i, j
type(SU2), intent(in) :: R
integer :: k

V%a = cmplx(0.0_dp,0.0_dp)

do k=1,3
   V%a(k,k) = cmplx(1.0_dp,0.0_dp)
end do

V%a(i,i) = R%a(1)
V%a(j,i) = -conjg(R%a(2))
V%a(i,j) = R%a(2)
V%a(j,j) = conjg(R%a(1))
 
end subroutine
!==============================

!==============================
!Computes the determinant of an 2x2 matrix
!which can be parametrized as an SU(2) matrix
!(whithout the contraint a_i*a_i = 1)
real(dp) function detSU2_like(U)
type(SU2), intent(in) :: U

detSU2_like = U%a(1)*conjg(U%a(1)) + U%a(2)*conjg(U%a(2)) 
end function detSU2_like 
!==============================


!==============================
!Computes the determinant of an 3x3 complex matrix
!encapsulated in a SU3 struct
complex(dp) function detSU3_like(U)
type(SU3), intent(in) :: U
integer :: i,j,k

detSU3_like = cmplx(0.0_dp, 0.0_dp)
do i=1,3
   do j=1,3
      do k=1,3
         detSU3_like = detSU3_like + levi_civita(i,j,k)*U%a(i,1)*U%a(j,2)*U%a(k,3)
      end do
   end do
end do 
end function detSU3_like 
!==============================

!==============================
!Multiply two SU(2)-like matrices
subroutine SU2mult(U1,U2,U3)
type(SU2), intent(in)  :: U1,U2
type(SU2), intent(out) :: U3

U3%a(1) = U1%a(1)*U2%a(1) - U1%a(2)*conjg(U2%a(2))
U3%a(2) = U1%a(1)*U2%a(2) + U1%a(2)*conjg(U2%a(1))

end subroutine SU2mult
!==============================

!==============================
!Takes the real part of the trace of a SU(3) matrix
real(dp) function SU3_ReTr(V)
type(SU3), intent(in) :: V
integer :: i
SU3_ReTr = real(V%a(1,1)) + real(V%a(2,2)) + real(V%a(3,3))
end function SU3_ReTr
!==============================

!==============================
!Checks if the input has deviated from SU(3) and in case
!so, projects it again on the group
subroutine SU3projector(A) !Needs finish implementing
type(SU3), intent(inout) :: A
real(dp), parameter :: tol = 1.0e-15
real(dp) :: norm
complex(dp) :: angle

integer :: i

If ( is_not_SU3(A) ) then
  !Computes the norm of the first column
   norm = sqrt(dot_product(real(A%a(:,1)),real(A%a(:,1))) + dot_product(imag(A%a(:,1)),imag(A%a(:,1))))
   A%a(:,1) = A%a(:,1)/norm
   !Computes the dot product between first and second columns
   angle = cmplx(dot_product(real(A%a(:,1)),real(A%a(:,2))) + dot_product(imag(A%a(:,1)),imag(A%a(:,2))),&
                 dot_product(imag(A%a(:,2)),real(A%a(:,1))) - dot_product(real(A%a(:,2)),imag(A%a(:,1))))
   A%a(:,2) = A%a(:,2) - A%a(:,1)*angle
   !Normalize the second column
   norm = sqrt(dot_product(real(A%a(:,2)),real(A%a(:,2))) + dot_product(imag(A%a(:,2)),imag(A%a(:,2))))
   A%a(:,2) = A%a(:,2)/norm
   !Third column is obtained as the cross product between first and second columns   
   do i=1,3
      A%a(i,3) =  cross_product(conjg(A%a(:,1)),conjg(A%a(:,2)),i) 
   end do
   !For debug purposes, I do one more check. In case this fails. I print an error
!   if (is_not_SU3(A)) then
!      print *, "Reunitarization of group element failed. Exitting now."
!      stop -1
!   end if
end if

end subroutine SU3projector
!==============================

!==============================
!Checks if the input has the properties of a SU(3) matrix
logical function is_not_SU3(A)
type(SU3), intent(in) :: A

if( abs(sqrt(dot_product(real(A%a(:,1)),real(A%a(:,1))) + dot_product(imag(A%a(:,1)),imag(A%a(:,1)))) - 1.0_dp) .lt. tol .and. & 
    abs(sqrt(dot_product(real(A%a(:,2)),real(A%a(:,2))) + dot_product(imag(A%a(:,2)),imag(A%a(:,2)))) - 1.0_dp) .lt. tol .and. & 
    abs(sqrt(dot_product(real(A%a(:,3)),real(A%a(:,3))) + dot_product(imag(A%a(:,3)),imag(A%a(:,3)))) - 1.0_dp) .lt. tol .and. &
    abs(dot_product(real(A%a(:,1)),real(A%a(:,2))) + dot_product(imag(A%a(:,1)),imag(A%a(:,2)))) .lt. tol .and. & 
    abs(dot_product(imag(A%a(:,1)),real(A%a(:,2))) - dot_product(real(A%a(:,1)),imag(A%a(:,2)))) .lt. tol .and. & 
    abs(dot_product(real(A%a(:,1)),real(A%a(:,3))) + dot_product(imag(A%a(:,1)),imag(A%a(:,3)))) .lt. tol .and. & 
    abs(dot_product(imag(A%a(:,1)),real(A%a(:,3))) - dot_product(real(A%a(:,1)),imag(A%a(:,3)))) .lt. tol .and. & 
    abs(dot_product(real(A%a(:,2)),real(A%a(:,3))) + dot_product(imag(A%a(:,2)),imag(A%a(:,3)))) .lt. tol .and. & 
    abs(dot_product(imag(A%a(:,2)),real(A%a(:,3))) - dot_product(real(A%a(:,2)),imag(A%a(:,3)))) .lt. tol) then
   is_not_SU3 = .false.
else
   is_not_SU3 = .true.
   !print *, "List of failed tests"
   !write(6,*) abs(sqrt(dot_product(real(A%a(:,1)),real(A%a(:,1))) + dot_product(imag(A%a(:,1)),imag(A%a(:,1)))) - 1.0_dp)   
   !write(6,*) abs(sqrt(dot_product(real(A%a(:,2)),real(A%a(:,2))) + dot_product(imag(A%a(:,2)),imag(A%a(:,2)))) - 1.0_dp)   
   !write(6,*) abs(sqrt(dot_product(real(A%a(:,3)),real(A%a(:,3))) + dot_product(imag(A%a(:,3)),imag(A%a(:,3)))) - 1.0_dp)  
   !write(6,*) abs(dot_product(real(A%a(:,1)),real(A%a(:,2))) + dot_product(imag(A%a(:,1)),imag(A%a(:,2))))   
   !write(6,*) abs(dot_product(imag(A%a(:,1)),real(A%a(:,2))) - dot_product(real(A%a(:,1)),imag(A%a(:,2)))) 
   !write(6,*) abs(dot_product(real(A%a(:,1)),real(A%a(:,3))) + dot_product(imag(A%a(:,1)),imag(A%a(:,3))))   
   !write(6,*) abs(dot_product(imag(A%a(:,1)),real(A%a(:,3))) - dot_product(real(A%a(:,1)),imag(A%a(:,3))))   
   !write(6,*) abs(dot_product(real(A%a(:,2)),real(A%a(:,3))) + dot_product(imag(A%a(:,2)),imag(A%a(:,3))))   
   !write(6,*) abs(dot_product(imag(A%a(:,2)),real(A%a(:,3))) - dot_product(real(A%a(:,2)),imag(A%a(:,3))))
end if
   
end function is_not_SU3
!==============================
   
!==============================
complex(dp) function cross_product(A,B,i)
complex(dp), intent(in) :: A(3), B(3) !A cross product is usually defined only in dimension 3
integer, intent(in) :: i !We compute only the i-th component 
integer :: j,k
   
cross_product = 0.0_dp
do j=1,3
   do k=1,3
      cross_product = cross_product + levi_civita(i,j,k)*A(j)*B(k)
   end do
end do

end function cross_product
!==============================

!==============================
!Inverts an 3x3 complex matrix
subroutine invert_3x3_complex(U,V)
type(SU3),intent(in) :: U
type(SU3),intent(out) :: V
complex(dp) :: detU

detU = detSU3_like(U)
if ( abs(real(detU)) .gt. tol .or. abs(imag(detU)) .gt. tol) then !If matrix not singular
   V%a(1,1) = U%a(2,2)*U%a(3,3) - U%a(3,2)*U%a(2,3)
   V%a(1,2) = U%a(1,3)*U%a(3,2) - U%a(3,3)*U%a(1,2)
   V%a(1,3) = U%a(1,2)*U%a(2,3) - U%a(2,2)*U%a(1,3)
   V%a(2,1) = U%a(2,3)*U%a(3,1) - U%a(3,3)*U%a(2,1)
   V%a(2,2) = U%a(1,1)*U%a(3,3) - U%a(3,1)*U%a(1,3)
   V%a(2,3) = U%a(1,3)*U%a(2,1) - U%a(2,3)*U%a(1,1)
   V%a(3,1) = U%a(2,1)*U%a(3,2) - U%a(3,1)*U%a(2,2)
   V%a(3,2) = U%a(1,2)*U%a(3,1) - U%a(3,2)*U%a(1,1)
   V%a(3,3) = U%a(1,1)*U%a(2,2) - U%a(2,1)*U%a(1,2)
   V%a = V%a/detU
end if
end subroutine invert_3x3_complex
!==============================

end module math
