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
public SU3projector
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

V%a(1) = rnor( )
V%a(2) = rnor( )
V%a(3) = rnor( )
V%a(4) = rnor( )
Norm = sqrt(V%a(1)**2+V%a(2)**2+V%a(3)**2+V%a(4)**2)
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
subroutine subgroup(V,i,j,R)
type(SU3), intent(in) :: V
real(dp), parameter :: f = 1.0_dp/3.0_dp
integer, intent(in) :: i, j
type(SU2), intent(out) :: R

!I assume i < j but the routine does not explicitly check for this. Beware!
R%a(1) = (V%im(j,i) + V%im(i,j))*f
R%a(2) = (V%re(i,j) - V%re(j,i))*f
R%a(3) = (V%im(i,i) - V%im(j,j))*f
R%a(4) = (V%re(i,i) + V%re(j,j))*f

end subroutine subgroup
!==============================

!==============================
subroutine embed_in_SU3(R,i,j,V)
type(SU3), intent(out) :: V
integer, intent(in) :: i, j
type(SU2), intent(in) :: R
integer :: k
V%re = 0.0_dp
V%im = 0.0_dp
do k=1,3
   if ( k .ne. i .and. k .ne. j) then
      V%re(k,k) = 1.0_dp
   end if
end do
V%re(i,i) = R%a(4)
V%im(i,i) = R%a(3)
V%re(j,i) = -R%a(2)
V%im(j,i) = R%a(1)
V%re(i,j) = R%a(2)
V%im(i,j) = R%a(1)
V%re(j,j) = R%a(4)
V%im(j,j) = -R%a(3)
 
end subroutine
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
!Checks if the input has deviated from SU(3) and in case
!so, projects it again on the group
subroutine SU3projector(A) !Needs finish implementing
type(SU3), intent(inout) :: A
real(dp), parameter :: tol = 1.0e-15
real(dp) :: norm,normIm
integer :: i

If ( is_not_SU3(A) ) then
   norm = sqrt( dot_product(A%re(1,:),A%re(1,:)) + dot_product(A%im(1,:),A%im(1,:)) )
   A%re(1,:) = A%re(1,:)/norm
   A%im(1,:) = A%im(1,:)/norm
   norm   = dot_product(A%re(2,:),A%re(1,:)) + dot_product(A%im(2,:),A%im(1,:))
   normIm = dot_product(A%im(2,:),A%re(1,:)) - dot_product(A%re(2,:),A%im(1,:))
   A%re(2,:) = A%re(2,:) - (A%re(1,:)*norm - A%im(1,:)*normIm)
   A%im(2,:) = A%im(2,:) - (A%re(1,:)*normIm + A%im(1,:)*norm)
   norm = sqrt( dot_product(A%re(2,:),A%re(2,:)) + dot_product(A%im(2,:),A%im(2,:)) )
   A%re(2,:) = A%re(2,:)/norm
   A%im(2,:) = A%im(2,:)/norm
   do i=1,3
      A%re(3,i) =  cross_product(A%re(1,:),A%re(2,:),i) - cross_product(A%im(1,:),A%im(2,:),i)
      A%im(3,i) = -cross_product(A%re(1,:),A%im(2,:),i) - cross_product(A%im(1,:),A%re(2,:),i)
   end do
   !For debug purposes, I do one more check. In case this fails. I print an error
   !if (is_not_SU3(A)) then
   !   print *, "Reunitarization of group element failed. Exitting now."
   !   stop -1
   !end if
end if

end subroutine SU3projector
!==============================

!==============================
!Checks if the input has the properties of a SU(3) matrix
logical function is_not_SU3(A)
type(SU3), intent(in) :: A

!The idea is to check if the lines form an orthonormal basis

!First, we check if all 3 vectors has modulus one
!I will adopt a cascate of ifs, since if one of the conditions is
!violated, I do not need to continue and returns that the element
!does not belongs to SU3


if ( abs( sqrt( dot_product(A%re(1,:),A%re(1,:)) + dot_product(A%im(1,:),A%im(1,:)) ) - 1.0_dp) .lt. tol .and. &   !Checks if first row has modulo 1
     abs( sqrt( dot_product(A%re(2,:),A%re(2,:)) + dot_product(A%im(2,:),A%im(2,:)) ) - 1.0_dp) .lt. tol .and. &   !Checks if second row has modulo 1 
     abs( dot_product(A%re(1,:),A%re(2,:)) + dot_product(A%im(1,:),A%im(2,:)) ) .lt. tol .and. &           !Checks if first line is orthogonal to second line (real part)
     abs( dot_product(A%im(1,:),A%re(2,:)) - dot_product(A%re(1,:),A%im(2,:)) ) .lt. tol .and. &           !Checks if first line is orthogonal to second line (imaginary part)
     abs( A%re(3,1) - cross_product(A%re(1,:),A%re(2,:),1) + cross_product(A%im(1,:),A%im(2,:),1) ) .lt. tol .and. &!Checks if first element of third line can be computed from the other two (real part)
     abs( A%im(3,1) + cross_product(A%re(1,:),A%im(2,:),1) + cross_product(A%im(1,:),A%re(2,:),1) ) .lt. tol .and. &!Checks if first element of third line can be computed from the other two (imaginary part)
     abs( A%re(3,2) - cross_product(A%re(1,:),A%re(2,:),2) + cross_product(A%im(1,:),A%im(2,:),2) ) .lt. tol .and. &!Checks if second element of third line can be computed from the other two (real part)
     abs( A%im(3,2) + cross_product(A%re(1,:),A%im(2,:),2) + cross_product(A%im(1,:),A%re(2,:),2) ) .lt. tol .and. &!Checks if second element of third line can be computed from the other two (imaginary part)  
     abs( A%re(3,3) - cross_product(A%re(1,:),A%re(2,:),3) + cross_product(A%im(1,:),A%im(2,:),3) ) .lt. tol .and. &!Checks if third element of third line can be computed from the other two (real part)
     abs( A%im(3,3) + cross_product(A%re(1,:),A%im(2,:),3) + cross_product(A%im(1,:),A%re(2,:),3) ) .lt. tol ) then !Checks if third element of third line can be computed from the other two (imaginary part)  
   
   is_not_SU3 = .false.
else
   is_not_SU3 = .true.
end if

end function is_not_SU3
!==============================
   
!==============================
real(dp) function cross_product(A,B,i)
real(dp), intent(in) :: A(3), B(3) !A cross product is usually defined only in dimension 3
integer, intent(in) :: i !We compute only the i-th component 
integer :: j,k
   
cross_product = 0.0_dp
do j=1,3
   do k=1,3
      cross_product = cross_product + levi_civita(i,j,k)*A(j)*B(k)
   end do
end do

end function cross_product
   
end module math
