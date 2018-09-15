module cuda_fermion_prop
use types_params, only : nx,ny,nz.nt,U,tol
use lattice, only : increment_table
use cudafor !Needed for CUDA Fortran
use cuBlas_v2 !Library to perform linear algebra operations (scalar product, and matrix-vector multiplication)
implicit none

!Public vars

!Private vars

type(SU3), device, allocatable :: Prop(:,:),U_d(:,:)
type(SU3), texture, pointer :: Utex(:,:) !Use texture, since U does not change during kernel execution and
                                         !needs strided access
complex(dp),device,allocatable :: source(:,:),r(:),r_tilde(:),dirac_vector(:),p(:),v(:),t(:),s(:)
complex(dp),dimension(4,4,8) :: gamma_mu,gamma_mu_star
complex(dp),dimension(4,4,8),constant :: gamma_mu_d,gamma_mu_star_d
integer, device, allocatable :: increment_table_d(:,:)

integer, constant :: nx_d,ny_d,nz_d,nt_d
real(dp), constant :: m_d

private
public :: compute_prop

contains

!=============================
!Compute the propagator - Only this subroutine should be public
subroutine compute_prop(out_prop)
   type(SU3), intent(out) :: out_prop(:,:)
   type(dim3) :: grid_spec,block_spec
   type(cublasHandle) :: handle
   integer, parameter :: max_step = 1000
   integer :: a,alpha, istat,i,array_dim,step
   logical :: stop_condition
   complex(dp) :: rho0,rho1,omega,beta_cg, alpha_cg
   real(dp) :: s_norm


   call set_gamma_mu

   array_dim = 12*nx*ny*nz*nt
   !Allocate GPU memory
   allocate(U_d(8,0:nx*ny*nz*nt-1))
   allocate(increment_table_d(0:nx*ny*nz*nt-1,8))
   allocate(source(0:array_dim-1,0:3*4))
   allocate(r(0:array_dim-1),r_tilde(0:array_dim-1),dirac_vector(0:array_dim-1))
   allocate(v(0:array_dim-1),t(0:array_dim-1),s(0:array_dim-1),p(0:array_dim-1))
   allocate(Prop(0:array_dim,0:3*4))

   !Load data to GPU
   U_d = U
   gamma_mu_d = gamma_mu
   gamma_mu_star_d = gamma_mu_star
   increment_table_d = increment_table
   nx_d = nx
   ny_d = ny
   nz_d = nz
   nt_d = nt
   m_d = 1.0_dp !Placeholder value

   !Bind texture
   Utex => U_d

   !Execute GPU code
   grid_spec=dim3(nx,ny,nz*nt)
   block_spec=dim3(12,12,1)

   !Set source
   !call set_source_checkboard<<<grid_spec,block_spec>>>()
   call set_source_point_source<<<grid_spec,block_spec>>>()
   !Prepare to use cublas
   istat = cublasCreate(handle)

   !For each spin-color index
   do alpha=0,3
      do a=0,2
         istat = cublasZcopy_v2(handle,array_dim,source(:,a+3*alpha),1,Prop(:,a+3*alpha),1)
         !Sets r0
         call dirac_vector_multiply(Prop(:,a+3*alpha),r)
         istat = cublasZaxpy_v2(handle,array_dim,cmplx(-1.0_dp,0.0_dp),source(:,a+3*alpha),1,r,1)
         istat = cublasZdscal_v2(handle,array_dim,cmplx(-1.0_dp,0.0_dp),r,1)
         istat = cublasZcopy_v2(handle,array_dim,r,1,r_dagger,1)

         stop_condition = .true.
         setp = 1
         do while (stop_condition)
            !Sets rho(i-1) = r_dagger x r(i-1)
            istat = cublasZdotc_v2(handle,array_dim,r_tilde,1,r,1,rho1)
            if ( sqrt( real(rho1)**2 + aimag(rho1)**2 ) > tol) then
               if (step .eq. 1) then
                  istat = cublasZcopy_v2(handle,array_dim,r,1,p,1)
               else
                  beta_cg = alpha_cg*rho1/(omega*rho0)
                  !p(i) = r(i-1) + beta(i-1)*( p(i-1) - omega(i-1)*v(i-1) )
                  !istat = cublasZdscal(handle,array_dim,-omega,v,1)
                  istat = cublasZaxpy_v2(handle,array_dim,-omega,v,1,p,1)
                  istat = cublasZdscal(handle,array_dim,beta_cg,p,1)
                  istat = cublasZaxpy(handle,array_dim,cmplx(-1.0_dp,0.0_dp),r,1,p,1,1)
                  !istat = cublasZcopy_v2(handle,array_dim,v,1,p,1)
               end if
               !TO DO: In case of preconditioned algorithm, invert M matrix and find a p_hat
               call dirac_vector_multiply(p,v) !v(i) = A p(i)
               !alpha(i) = rho(i-1)/(r_tilde_dagger x v(i))
               istat = cublasZdotc_v2(handle,array_dim,r_tilde,1,v,1,alpha_cg)
               alpha_cg=rho1/alpha_cg
               !s = r(i-1) - alpha(i)v(i)
               istat = cublasZaxpy_v2(handle,array_dim,-alpha_cg,v,1,r,1) !A new value of r will be set later
               istat = cublasDZnrm2_v2(handle,array_dim,r,1,s_norm)
               if (s_norm .gt. tol) then
                  !TO DO: In case of preconditioned algorithm, invert matrix M and find s_hat
                  call dirac_vector_multiply(r,t)
                  istat = cublasDZnrm2_v2(handle,array_dim,t,1,s_norm) !t_dagger x t stored in s_norm
                  istat = cublasZdotc_v2(handle,array_dim,t,1,r,1,omega) !t_dagger x s stored in omega
                  omega = omega/s_norm ! set omega
                  !r = s - omega t
                  istat = cublasZaxpy_v2(handle,array_dim,-omega,t,1,s,1)
                  istat = cublasZcopy_v2(handle,array_dim,s,1,r,1)
                  istat = cublasZaxpy_v2(handle,array_dim,alpha_cg,p,1,Prop(:,a+3*alpha_cg),1)
                  istat = cublasZaxpy_v2(handle,array_dim,omega,s,1,Prop(:,a+3*alpha_cg),1)
                  step = step + 1
                  if (step .gt. max_step) then
                     print *, "Max step reached without finding a convergence"
                     stop_condition = .false.
                  end if
               else
                  istat = cublas zaxpy_v2(handle,array_dim,alpha_cg,p,1,Prop(:,a+3*alpha),1)
                  stop_condition = .false.
               end if
            else
               print *, "ERROR: Division by zero found. Aborting"
               stop_condition = .false.
            end do
         end do
      end do
   end do

   !Return the result
   out_prop = Prop
   !Deallocate GPU memory
   deallocate(U_d,increment_table_d,source,r,r_tilde,dirac_vector,v,t,s,p,Prop)

   end subroutine compute_prop
!=============================

!=============================
!Multiply a vector by the dirac operator
subroutine dirac_vector_multiply(v,w)
complex(dp), device, intent(in) :: v
complex(dp), device, intent(out) :: w
integer :: i, beta,b, istat
    do i=0,12*nx*ny*nz*nt-1
        do beta=0,3
           do b=0,3
              call make_dirac_vector<<<dim3(nx,ny,nz*nt),dim3(12,1,1)>>>(b,beta,i,dirac_vector)
              istat = cublasZdotu_v2(handle,array_dim,dirac_vector,1,v,w(b + 3*beta + 12*i))
           end do
        end do
     end do
end subroutine

!=============================
!Sets the source on GPU memory
   !Checkboard pattern on x is used since we are
   !interested only on the propagator G(x| x + a mu), mu=1,2,3,4,5,6,7,8
attributes(global) subroutine set_source_checkboard()
   integer :: bx,by,bz,bt

   bx=blockIdx%x-1
   by=blockIdx%y-1
   bt=(blockIdx%z-1)/nt_d
   bz=blockIdx%z-1-bt*nz_d

   if ( (mod(bx,2) .eq. 0) .and. (mod(by,2) .eq. 0) .and. (mod(bz,2) .eq. 0) .and. (mod(bz,2) .eq. 0) ) then
      if (threadIdx%x .eq. threadIdx%y) then
         source((threadIdx%x-1) + bx*12 + by*12*nx_d + bz*12*nx_d*ny_d + bt*12*nx_d*ny_d*nz_d,threadIdx%y-1) = 1
      else
         source((threadIdx%x-1) + bx*12 + by*12*nx_d + bz*12*nx_d*ny_d + bt*12*nx_d*ny_d*nz_d,threadIdx%y-1) = 0
      end if
   else
      source((threadIdx%x-1) + bx*12 + by*12*nx_d + bz*12*nx_d*ny_d + bt*12*nx_d*ny_d*nz_d,threadIdx%y-1) = 0
   end if

end subroutine set_source
!=============================

!=============================
!Sets the source on GPU memory
   !Point source on x = 0
attributes(global) subroutine set_source_point_source()
   integer :: bx,by,bz,bt

   bx=blockIdx%x-1
   by=blockIdx%y-1
   bt=(blockIdx%z-1)/nt_d
   bz=blockIdx%z-1-bt*nz_d

   if ( bx + nx_d*by + nx_d*ny_d*bz + nx_d*ny_d*nz_d*bt .eq. 0) then
      if (threadIdx%x .eq. threadIdx%y) then
         source((threadIdx%x-1) + bx*12 + by*12*nx_d + bz*12*nx_d*ny_d + bt*12*nx_d*ny_d*nz_d,threadIdx%y-1) = 1
      else
         source((threadIdx%x-1) + bx*12 + by*12*nx_d + bz*12*nx_d*ny_d + bt*12*nx_d*ny_d*nz_d,threadIdx%y-1) = 0
      end if
   else
      source((threadIdx%x-1) + bx*12 + by*12*nx_d + bz*12*nx_d*ny_d + bt*12*nx_d*ny_d*nz_d,threadIdx%y-1) = 0
   end if

end subroutine set_source
!=============================

!=============================
!Gets the row a + 3*alpha + 12x of the Dirac operator
attributes(global) subroutine make_dirac_vector(a,alpha,x,w)
   integer, value :: a, alpha, x
   complex(dp), intent(out) :: w(:)
   integer :: bx,by,bz,bt,y,mu,beta,b

   bx=blockIdx%x-1
   by=blockIdx%y-1
   bt=(blockIdx%z-1)/nt_d
   bz=blockIdx%z-1-bt*nz_d
   y = bx + by*nx_d + bz*nx_d*ny_d + bt*nx_d*ny_d*nz_d

   beta = (threadIdx%x-1)/4
   b = threadIdx%x - 1 - alpha*4

   w(b+3*beta + 12*y) = cmplx(0.0_dp,0.0_dp)
   do mu=1,8
      if (y .eq. increment_table_d(mu,x)) then
         w(b+3*beta + 12*y) = w(b+3*beta + 12*y) + gamma_mu_star_d(alpha,beta,mu)*Utex(mu,x)%a(a,b)
      end if
   end do

   if ( (y .eq. x) .and. (a .eq. b) .and. (alpha .eq. beta) ) then
      w(b+3*beta + 12*y) = w(b+3*beta + 12*y) +  dcmplx(4.0_dp + m_d,0.0_dp)
   end if

end subroutine
!=============================

!=============================
!invert sign of array
attributes(global) subroutine invert_sign(v)
complex(dp) :: v(:)
integer :: bx,by,bz,bt,y,beta,b

   bx=blockIdx%x-1
   by=blockIdx%y-1
   bt=(blockIdx%z-1)/nt_d
   bz=blockIdx%z-1-bt*nz_d
   y = bx + by*nx_d + bz*nx_d*ny_d + bt*nx_d*ny_d*nz_d

   beta = (threadIdx%x-1)/4
   b = threadIdx%x - 1 - alpha*4

   v(b + 3*beta + y) = -v(b + 3*beta + y)

end subroutine
!=============================

!=============================
!Set up gamma matrices on color-spin space
subroutine set_gamma_mu()
    complex(dp),dimension(4,4,8) :: g, g_star
    integer :: a,alpha,beta,i,j
       g = cmplx(0.0_dp,0.0_dp)
       g_star = cmplx(0.0_dp,0.0_dp)


       !gamma_1
       g(0,3,2) = dcmplx(0.d0,1.d0)
       g(1,2,2) = dcmplx(0.d0,1.d0)
       g(2,1,2) = dcmplx(0.d0,-1.d0)
       g(3,0,2) = dcmplx(0.d0,-1.d0)

       g(:,:,1) = -g(:,:,2)

       !gamma_2
       g(0,3,4) = dcmplx(1.d0,0.d0)
       g(1,2,4) = dcmplx(-1.d0,0.d0)
       g(2,1,4) = dcmplx(-1.d0,0.d0)
       g(3,0,4) = dcmplx(1.d0,0.d0)

       g(:,:,3) = -g(:,:,4)

       !gamma_3
       g(0,2,6) = dcmplx(0.d0,1.d0)
       g(1,3,6) = dcmplx(0.d0,-1.d0)
       g(2,0,6) = dcmplx(0.d0,-1.d0)
       g(3,1,6) = dcmplx(0.d0,1.d0)

       g(:,:,5) = -g(:,:,6)

       !gamma_4
       g(0,2,8) = dcmplx(1.d0,0.d0)
       g(1,3,8) = dcmplx(1.d0,0.d0)
       g(2,0,8) = dcmplx(1.d0,0.d0)
       g(3,1,8) = dcmplx(1.d0,0.d0)

       g(:,:,7) = -g(:,:,8)

       !g_star is the negative of g
       g_star = -g
       !except on the diagonal components
       do mu=1,8
          do alpha=1,4
          g_star(alpha,alpha,mu) = 1 - g(alpha,alpha,mu)
          end do
       end do

       gamma_mu = g
       gamma_mu_star = g_star

    end subroutine
    !=============================
end module
