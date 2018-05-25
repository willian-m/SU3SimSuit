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
!==============================

private
public gen_su3_element,SU3mult,SU2mult,detSU2_like,SU3_dagger,SU3_ReTr
public embedR,embedS,embedT,subgroupR,subgroupS,subgroupT
!public SU3check !Used for debug. If not debugging, comment this line
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

call embedR(R,R3)
call embedS(S,S3)
call embedT(T,T3)

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

V%a(1) = rnor( )
V%a(2) = rnor( )
V%a(3) = rnor( )
V%a(4) = rnor( )
Norm = dsqrt(V%a(1)**2+V%a(2)**2+V%a(3)**2+V%a(4)**2)
V%a(1) = V%a(1)/Norm
V%a(2) = V%a(2)/Norm
V%a(3) = V%a(3)/Norm
V%a(4) = V%a(4)/Norm
end subroutine gen_su2_element
!==============================

!==============================
!Multiplies SU(3) group elements
subroutine SU3mult(A,B,C)
type(SU3), intent(in) :: A,B
type(SU3), intent(out) :: C
real(dp) :: aux(3,3)
integer :: i,j,l

C%re = matmul(A%re,B%re)
aux = matmul(A%im,B%im)
C%re = C%re - aux

C%im = matmul(A%im,B%re)
aux  = matmul(A%re,B%im)
C%im = C%im + aux

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
      B%re(i,j) = A%re(j,i)
      B%im(i,j) = -A%im(j,i)
   end do
end do

end subroutine SU3_dagger
!==============================

!==============================
!Takes the 2x2 upper-left block matrix of U
!and write it in the shape of an SU(2) matrix.
subroutine subgroupR(U,R)
type(SU3), intent(in) :: U
type(SU2), intent(out) :: R

R%a(4) = U%re(1,1) + U%re(2,2)
R%a(1) = -(U%im(2,1) + U%im(1,2)
R%a(3) = -(U%im(1,1) - U%im(2,2)
R%a(2) = -(U%re(1,2) - U%re(2,1)

end subroutine subgroupR
!==============================

!==============================
!Takes the corner elements of U
!and write it in the shape of an SU(2) matrix.
subroutine subgroupS(U,S)
type(SU3), intent(in) :: U
type(SU2), intent(out) :: S

S%a(4) = U%re(1,1) + U%re(3,3)
S%a(1) = -(U%im(3,1) + U%im(1,3))
S%a(3) = -(U%im(1,1) - U%im(3,3))
S%a(2) = -(U%re(1,3) - U%re(3,1))

end subroutine subgroupS
!==============================

!==============================
!Takes the 2x2 lower-right block matrix of U
!and write it in the shape of an SU(2) matrix.
subroutine subgroupT(U,T)
type(SU3), intent(in) :: U
type(SU2), intent(out) :: T

T%a(4) = U%re(2,2) + U%re(3,3)
T%a(1) = -(U%im(3,2) + U%im(2,3)
T%a(3) = -(U%im(2,2) - U%im(3,3)
T%a(2) = -(U%re(2,3) - U%re(3,2)

end subroutine subgroupT
!==============================

!==============================
!Embeds an SU(2) group element as an SU(3)
!element belonging on the R type SU(2) subgroup
subroutine embedR(R,R3)
type(SU2), intent(in) :: R
type(SU3), intent(out) :: R3

!First line of R3
R3%re(1,1) = R%a(4)
R3%im(1,1) = R%a(3)
R3%re(1,2) = R%a(2)
R3%im(1,2) = R%a(1)
R3%re(1,3) = 0.0_dp
R3%im(1,3) = 0.0_dp
!Second line of R3
R3%re(2,1) = -R%a(2)
R3%im(2,1) = R%a(1)
R3%re(2,2) = R%a(4)
R3%im(2,2) = -R%a(3)
R3%re(2,3) = 0.0_dp
R3%im(2,3) = 0.0_dp
!Third line of R3
R3%re(3,1) = 0.0_dp
R3%im(3,1) = 0.0_dp
R3%re(3,2) = 0.0_dp
R3%im(3,2) = 0.0_dp
R3%re(3,3) = 1.0_dp
R3%im(3,3) = 0.0_dp

end subroutine embedR
!==============================

!==============================
!Embeds an SU(2) group element as an SU(3)
!element belonging on the S type SU(2) subgroup
subroutine embedS(S,S3)
type(SU2), intent(in) :: S
type(SU3), intent(out) :: S3

!First line of S3
S3%re(1,1) = S%a(4)
S3%im(1,1) = S%a(3)
S3%re(1,2) = 0.0_dp
S3%im(1,2) = 0.0_dp
S3%re(1,3) = S%a(2)
S3%im(1,3) = S%a(1)
!Second line of S3
S3%re(2,1) = 0.0_dp
S3%im(2,1) = 0.0_dp
S3%re(2,2) = 1.0_dp
S3%im(2,2) = 0.0_dp
S3%re(2,3) = 0.0_dp
S3%im(2,3) = 0.0_dp
!Third line of S3
S3%re(3,1) = -S%a(2)
S3%im(3,1) = S%a(1)
S3%re(3,2) = 0.0_dp
S3%im(3,2) = 0.0_dp
S3%re(3,3) = S%a(4)
S3%im(3,3) = -S%a(3)

end subroutine embedS
!==============================

!==============================
!Embeds an SU(2) group element as an SU(3)
!element belonging on the T type SU(2) subgroup
subroutine embedT(T,T3)
type(SU2), intent(in) :: T
type(SU3), intent(out) :: T3

!First line of T3
T3%re(1,1) = 1.0_dp
T3%im(1,1) = 0.0_dp
T3%re(1,2) = 0.0_dp
T3%im(1,2) = 0.0_dp
T3%re(1,3) = 0.0_dp
T3%im(1,3) = 0.0_dp
!Second line of T3
T3%re(2,1) = 0.0_dp
T3%im(2,1) = 0.0_dp
T3%re(2,2) = T%a(4)
T3%im(2,2) = T%a(3)
T3%re(2,3) = T%a(2)
T3%im(2,3) = T%a(1)

!Third line of T3
T3%re(3,1) = 0.0_dp
T3%im(3,1) = 0.0_dp
T3%re(3,2) = -T%a(2)
T3%im(3,2) = T%a(1)
T3%re(3,3) = T%a(4)
T3%im(3,3) = -T%a(3)

end subroutine embedT
!==============================

!==============================
!Computes the determinant of an 2x2 matrix
!which can be parametrized as an SU(2) matrix
!(whithout the contraint a_i*a_i = 1)
real(dp) function detSU2_like(U)
type(SU2), intent(in) :: U
integer :: i

detSU2_like = 0.0_dp
do i=1,4
   detSU2_like = detSU2_like + U%a(i)**2
end do
end function detSU2_like 
!==============================

!==============================
!Multiply two SU(2)-like matrices
subroutine SU2mult(U1,U2,U3)
type(SU2), intent(in)  :: U1,U2
type(SU2), intent(out) :: U3

U3%a(4) = U1%a(4)*U2%a(4)-U1%a(3)*U2%a(3)-U1%a(2)*U2%a(2)-U1%a(1)*U2%a(1)
!c_k  = a_4   * b_k  +  a_k   * b_4   - a_i*b_j*eps_{i j k}
U3%a(1) = U1%a(4) * U2%a(1) + U1%a(1) * U2%a(4) - U1%a(2) * U2%a(3) + U1%a(3) * U2%a(2)
U3%a(2) = U1%a(4) * U2%a(2) + U1%a(2) * U2%a(4) - U1%a(3) * U2%a(1) + U1%a(1) * U2%a(3)
U3%a(3) = U1%a(4) * U2%a(3) + U1%a(3) * U2%a(4) - U1%a(1) * U2%a(2) + U1%a(2) * U2%a(1)
end subroutine SU2mult
!==============================

!==============================
!Takes the real part of the trace of a SU(3) matrix
real(dp) function SU3_ReTr(V)
type(SU3), intent(in) :: V
integer :: i
SU3_ReTr = 0.0_dp
do i=1,3
   SU3_ReTr = SU3_ReTr + V%re(i,i)
end do
end function SU3_ReTr
!==============================

!==============================
!Checks if the input has the properties of a SU(3) matrix
!Created for debug purposes
subroutine SU3check(A)
type(SU3), intent(in) :: A
type(SU3) :: A_dagger,unity
integer :: i,j,k,l

!Check 1: Multiplying by its dagger should result in unity
call SU3_dagger(A,A_dagger)
call SU3mult(A,A_dagger,unity)

do i=1,3
   do j=1,3
      if ( dabs(unity%im(i,j)) .gt. 1.0e-14_dp ) then
            print *, "Warning: I generated a non-SU(3) element!"
            print *, "Here is the real part of the offending matrix:"
            do k=1,3
               write(6,*) (A%re(l,k),l=1,3)
            end do
            print *, "Here is the imaginary part of the offending matrix:"
            do k=1,3
               write(6,*) (A%im(l,k),l=1,3)
            end do
            print *, "This matrix is not SU(3) because V . V^\dagger is:"
            print *, "Real part:"
            do k=1,3
               write(6,*) (unity%re(l,k),l=1,3)
            end do
            print *, "Imaginary part:"
            do k=1,3
               write(6,*) (unity%im(l,k),l=1,3)
            end do
            print *, "Stopping now"
            stop(1)
      end if
      if ( i .ne. j ) then
         if ( dabs(unity%re(i,j)) .gt. 1.0e-14_dp ) then
            print *, "Warning: I generated a non-SU(3) element!"
            print *, "Here is the real part of the offending matrix:"
            do k=1,3
               write(6,*) (A%re(l,k),l=1,3)
            end do
            print *, "Here is the imaginary part of the offending matrix:"
            do k=1,3
               write(6,*) (A%im(l,k),l=1,3)
            end do
            print *, "This matrix is not SU(3) because V . V^\dagger is:"
            print *, "Real part:"
            do k=1,3
               write(6,*) (unity%re(l,k),l=1,3)
            end do
            print *, "Imaginary part:"
            do k=1,3
               write(6,*) (unity%im(l,k),l=1,3)
            end do
            print *, "Stopping now"
            stop(1)
         end if
      else
         if ( dabs(unity%re(i,j) - 1.0_dp ) .gt. 1.0e-14_dp ) then
            print *, "Warning: I generated a non-SU(3) element!"
            print *, "Here is the real part of the offending matrix:"
            do k=1,3
               write(6,*) (A%re(l,k),l=1,3)
            end do
            print *, "Here is the imaginary part of the offending matrix:"
            do k=1,3
               write(6,*) (A%im(l,k),l=1,3)
            end do
            print *, "This matrix is not SU(3) because V . V^\dagger is:"
            print *, "Real part:"
            do k=1,3
               write(6,*) (unity%re(l,k),l=1,3)
            end do
            print *, "Imaginary part:"
            do k=1,3
               write(6,*) (unity%im(l,k),l=1,3)
            end do
            print *, "Stopping now"
            stop(1)
         end if
      end if
   end do
end do


end subroutine
!==============================
end module math
