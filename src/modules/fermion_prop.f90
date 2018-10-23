module fermion_prop
   !use cudafor
   use cublas
   use types_params, only : dp,nx,ny,nz,nt,SU3,U,tol
   use lattice, only : increment_table
   implicit none
   
   
   integer :: i_one,array_dim
   complex(dp),dimension(4,4,8) :: gamma_mu,gamma_mu_star !GPU variables to store dirac matrices
   complex(dp) :: mass
   
   
   private
   public :: compute_prop

   contains

   !=============================
   !Compute the propagator using the BI-CGStab algorithm
   subroutine compute_prop(prop)
      integer, parameter :: max_step = 100
      complex(dp), intent(out) :: prop(:,:)
      complex(dp), allocatable :: source(:,:)
      type(cublasHandle) :: handle
      integer :: istat,alpha,a,step,spin_color_index
      logical :: stop_condition
      complex(dp) :: rho1, omega, alpha_cg, rho0
      real(dp) :: s_norm
      complex(dp) :: Z_minus_one,Z_one

      !Auxiliary vectors used on the Bi-CGStab
      complex(dp), allocatable :: r(:),r_tilde(:),p(:),v(:),t(:)
      complex(dp) :: beta_cg

      
      call set_gamma_mu !Set up gamma matrices on device memory

      array_dim = 12*nx*ny*nz*nt

      !Allocate resources
      allocate(source(0:array_dim - 1,12))
      allocate(r(array_dim),r_tilde(array_dim),p(array_dim),v(array_dim))
      allocate(t(array_dim))
      
      !Load memory
      
      mass = cmplx(1.0_dp,0.0_dp)
      z_minus_one = cmplx(-1.0_dp,0.0_dp)
      z_one = cmplx(1.0_dp,0.0_dp)
      i_one = 1
      
      call set_source_point_source(source)
      !call set_source_checkboard(source)
      
      !For each spin-color index
      do alpha=1,4
         do a=1,3
            spin_color_index = a + 3*(alpha-1)
            call init_prop(a, alpha, prop(:,spin_color_index))
            
            !Initialization
            call dirac_vector_multiply(prop(:,spin_color_index),r)
            call zaxpy(array_dim,z_minus_one,source(:,spin_color_index),i_one,r,i_one) !constant times a vector plus a vector
            call zscal(array_dim,z_minus_one,r,i_one) !scales a vector by a constant
            call zcopy(array_dim,r,i_one,r_tilde,i_one)
            
            stop_condition = .false.
            step = 1
            do 
               if (stop_condition) exit
               !Begin computes p
               if (step .eq. 1) then
                  rho1 = zdotc(array_dim,r_tilde,i_one,r,i_one)   !Needed for next step
                  call zcopy(array_dim,r,i_one,p,i_one)
               else
                  rho0 = rho1
                  rho1 = zdotc(array_dim,r_tilde,i_one,r,i_one)
                  if ( sqrt( real(rho1)**2 + aimag(rho1)**2 ) > tol) then
                     beta_cg = alpha_cg*rho1/(omega*rho0)
                     call zaxpy(array_dim,-omega,v,i_one,p,i_one)
                     call zscal(array_dim,beta_cg,p,i_one)
                     call zaxpy(array_dim,z_one,r,i_one,p,i_one)
                  else
                     print *, "ERROR: Division by zero found at iteration", step,". Aborting"
                     exit
                  end if
               end if

               !TO DO: In case of preconditioned algorithm, invert M matrix and find a p_hat
               call dirac_vector_multiply(p,v) !v(i) = A p(i)
               
               alpha_cg = zdotc(array_dim,r_tilde,i_one,v,i_one)
               alpha_cg=rho1/alpha_cg
               call zaxpy(array_dim,-alpha_cg,v,i_one,r,i_one) !r = s for now. Later, it will be discarded and a new r computed
               !End compute s
               s_norm = DZnrm2(array_dim,r,i_one)
               if (s_norm .gt. tol) then
                  !TO DO: In case of preconditioned algorithm, invert matrix M and find s_hat
                  call dirac_vector_multiply(r,t)
                  s_norm = DZnrm2(array_dim,t,i_one) !t_dagger x t stored in s_norm
                  omega = zdotc(array_dim,t,i_one,r,i_one) !t_dagger x s stored in omega
                  omega = omega/(s_norm**2) ! set omega
                  !Update propagator
                  call zaxpy(array_dim,alpha_cg,p,i_one,prop(:,spin_color_index),i_one)
                  call zaxpy( array_dim,omega,r,i_one,prop(:,spin_color_index),i_one)
                  !Compute r for next iteration
                  call zaxpy(array_dim,-omega,t,i_one,r,i_one)
                  step = step + 1
                  if (step .gt. max_step) then
                     print *, "Max step reached without finding a convergence"
                     stop_condition = .true.
                  end if
               else
                  call zaxpy(array_dim,alpha_cg,p,i_one,prop(:,spin_color_index),i_one)
                  stop_condition = .true.
               end if
            end do
         end do
      end do

      print *, "Finished. Copying result to output vector"

      deallocate(source,U,r,r_tilde,p,v,t)
      

   end subroutine
   !=============================
   
   !=============================
   !Multiply a vector by the dirac operator
   subroutine dirac_vector_multiply(in_array,out_array)
      complex(dp), intent(in) :: in_array(:)
      complex(dp), intent(out) :: out_array(:)
      complex(dp) :: reduce_var
      integer :: x,y,alpha,beta,a,b,mu,i,j
      
      out_array = cmplx(0.0_dp,0.0_dp)
      
      !$acc parallel copyin(gamma_mu_star(1:4,1:4,1:8), U(8,0:array_dim/12-1),in_array(1:array_dim),i,x,alpha,a) copyout(out_array(1:array_dim))
      !$acc loop
      do x=0,array_dim/12-1
         !$acc loop
         do alpha=1,4 
            !$acc loop
            do a=1,3
               i = a + 3*(alpha-1) + 12*x
               
               !We build the i-th line of the dirac matrix
               reduce_var = cmplx(0.0_dp,0.0_dp)
               !$acc loop reduction(+:reduce_var) private(j,y)
               do beta=1,4
                  !$acc loop 
                  do b=1,3 
                     !Gauge dependent term
                     do mu=1,8
                        y = x + increment_table(x,mu)
                        j = b + 3*(beta-1) + 12*y
                        reduce_var = reduce_var + cmplx(-0.5_dp,0.0_dp)*gamma_mu_star(alpha,beta,mu)*U(mu,y)%a(a,b)*in_array(j)
                        !out_array(i) = out_array(i) - cmplx(0.5_dp,0.0_dp)*gamma_mu_star(alpha,beta,mu)*U(mu,y)%a(a,b)*in_array(j)
                     end do
                  end do
               end do

               !Mass term uses only diagonals, i.e. a = b, beta = alpha and y = x
               j = a + 3*(alpha-1) + 12*x
               reduce_var = reduce_var + (mass + cmplx(4.0_dp,0.0_dp))*in_array(j)
               !out_array(i) = out_array(i) + (mass + cmplx(4.0_dp,0.0_dp))*in_array(j)
            
               out_array(i) = reduce_var
            end do
         end do
      end do
      !$acc end parallel
   end subroutine
   !=============================

   !=============================
   !Gets the row a + 3*(alpha-1) + 12*x of the Dirac operator
   subroutine init_prop(a,alpha, v)
      integer, intent(in) :: a, alpha
      complex(dp), intent(out) :: v(:)
      integer :: i,x,y,z,t,b,beta

      do i=1,array_dim
         t=i/(12*nx*ny*nz)
         z = (i - t*12*nx*ny*nz)/(12*nx*ny)
         y = (i - t*12*nx*ny*nz - z*12*nx*ny)/(12*nx)
         x = (i - t*12*nx*ny*nz - z*12*nx*ny - y*12*nx)/12
         beta = (i - t*12*nx*ny*nz - z*12*nx*ny - y*12*nx - x*12)/3
         b = i - t*12*nx*ny*nz - z*12*nx*ny - y*12*nx - x*12 - 3*beta
         b = b + 1
         beta=beta+1
         if ((b .eq. a) .and. (beta .eq. alpha)) then
            v(i) = cmplx(mass+2.0_dp,0.0_dp)/(sin(real(x)) **2 + sin(real(y))**2 + sin(real(z))**2 + sin(real(t))**2 + mass + cmplx(2.0_dp))
         else
            v(i) = cmplx(0.0_dp,0.0_dp)
         end if
      end do
   end subroutine
   !=============================

   !=============================
   !Sets the source on GPU memory
   !Checkboard pattern on x is used since we are
   !interested only on the propagator G(x| x + a mu), mu=1,2,3,4,5,6,7,8
   subroutine set_source_checkboard(source)
      integer :: i,j,bx,by,bz,bt,x
      complex(dp), intent(out) :: source(:,:)

      do x=0,array_dim/12 - 1
         bt = x/(nx*ny*nz)
         bz = (x - bt*nz*ny*nx)/(ny*nx)
         by = (x - bz*nx*ny - bt*nz*ny*nx)/nx
         bx = x - by*nx - bz*nx*ny - bt*nz*ny*nx
         
         do i=1,12
            do j=1,12
               if ( (mod(bx,2) .eq. 0) .and. (mod(by,2) .eq. 0) .and. (mod(bz,2) .eq. 0) .and. (mod(bz,2) .eq. 0) ) then
                  if (i .eq. j) then
                    source(i + 12*x,j) = cmplx(1.0_dp,0.0_dp)
                  else
                    source(i + 12*x,j) = cmplx(0.0_dp,0.0_dp)
                  end if
               else
                  source(i + 12*x,j) = cmplx(0.0_dp,0.0_dp)
               end if
            end do
         end do
      end do

   end subroutine set_source_checkboard
   !=============================

   !=============================
   !Sets the source on GPU memory
   !Point source on x = 0
   subroutine set_source_point_source(source)
      integer :: i,j,x
      complex(dp), intent(out) :: source(:,:)


      do x=0,array_dim/12 - 1
         do i=1,12
            do j=1,12
               if ( x .eq. 0) then
                  if (i .eq. j) then
                     source(i + x*12,j) = CMPLX(1.0_dp,0.0_dp)
                  else
                     source(i + x*12,j) = CMPLX(0.0_dp,0.0_dp)
                  end if
               else
                  source(i + x*12,j) = CMPLX(0.0_dp,0.0_dp)
               end if
            end do
         end do
      end do
   end subroutine set_source_point_source
   !=============================

   !=============================
   !Set up gamma matrices on color-spin space
   subroutine set_gamma_mu()
      
      integer :: alpha,mu

      gamma_mu = cmplx(0.0_dp,0.0_dp)
      gamma_mu_star = cmplx(0.0_dp,0.0_dp)


      !gamma_mu_1
      gamma_mu(1,4,2) = dcmplx(0.0_dp,1.0_dp)
      gamma_mu(2,3,2) = dcmplx(0.0_dp,1.0_dp)
      gamma_mu(3,2,2) = dcmplx(0.0_dp,-1.0_dp)
      gamma_mu(4,1,2) = dcmplx(0.0_dp,-1.0_dp)

      gamma_mu(:,:,1) = -gamma_mu(:,:,2)!gamma_{-mu} = -gamma_mu.

      !gamma_mu_2
      gamma_mu(1,4,4) = dcmplx(1.0_dp,0.0_dp)
      gamma_mu(2,3,4) = dcmplx(-1.0_dp,0.0_dp)
      gamma_mu(3,2,4) = dcmplx(-1.0_dp,0.0_dp)
      gamma_mu(4,1,4) = dcmplx(1.0_dp,0.0_dp)

      gamma_mu(:,:,3) = -gamma_mu(:,:,4)!gamma_{-mu} = -gamma_mu.

      !gamma_mu_3
      gamma_mu(1,3,6) = dcmplx(0.0_dp,1.0_dp)
      gamma_mu(2,4,6) = dcmplx(0.0_dp,-1.0_dp)
      gamma_mu(3,1,6) = dcmplx(0.0_dp,-1.0_dp)
      gamma_mu(4,2,6) = dcmplx(0.0_dp,1.0_dp)

      gamma_mu(:,:,5) = -gamma_mu(:,:,6)!gamma_{-mu} = -gamma_mu.

      !gamma_mu_4
      gamma_mu(1,3,8) = dcmplx(1.0_dp,0.0_dp)
      gamma_mu(2,4,8) = dcmplx(1.0_dp,0.0_dp)
      gamma_mu(3,1,8) = dcmplx(1.0_dp,0.0_dp)
      gamma_mu(4,2,8) = dcmplx(1.0_dp,0.0_dp)

      gamma_mu(:,:,7) = -gamma_mu(:,:,8)!gamma_{-mu} = -gamma_mu.

      !gamma_mu_star is the negative of gamma_mu
      gamma_mu_star = -gamma_mu
      !except on the diagonal components
      do mu=1,8
         do alpha=1,4
         gamma_mu_star(alpha,alpha,mu) = cmplx(1.0_dp,0.0_dp) - gamma_mu(alpha,alpha,mu)
         end do
      end do


   end subroutine
   !=============================
end module