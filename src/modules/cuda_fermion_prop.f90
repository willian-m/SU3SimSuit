module cuda_fermion_prop
use types_params
use cudafor !Needed for CUDA Fortran
implicit none

!Public vars

!Local vars
type(SU3), device, allocatable :: U_d(:,:), Prop(:,:)
integer,device,allocatable :: source(:,:) !Integer, since it only admits zero and 1

private
public :: compute_prop

contains

!=============================
!Compute the propagator - Only this subroutine should be public
subroutine compute_prop
   type(dim3) :: grid_spec,block_spec

   !Allocate and copy data to GPU
   allocate(U_d(8,0:nx*ny*nz*nt-1))
   allocate(source(0:3*4,0:3*4*nx*ny*nz*nt-1))
   U_d = U
   
   !Execute GPU code
   grid_spec=dim3(nx,ny,nz*nt)
   block_spec=dim3(12,12,1)

   call set_source<<<grid_spec,block_spec>>>(nx,ny,nz,nt)
   
   !Deallocate GPU memory
   deallocate(U_d,source)

   end subroutine compute_prop
!=============================

!=============================
!Sets the source on GPU memory
attributes(global) subroutine set_source(nx_d,ny_d,nz_d,nt_d)
   integer :: bx,by,bz,bt
   integer, value, intent(in) :: nx_d,ny_d,nz_d,nt_d
   
   bx=blockIdx%x-1
   by=blockIdx%y-1
   bt=(blockIdx%z-1)/nt_d
   bz=blockIdx%z-1-bt*nz_d

   if ( (mod(bx,2) .eq. 0) .and. (mod(by,2) .eq. 0) .and. (mod(bz,2) .eq. 0) .and. (mod(bz,2) .eq. 0) ) then
      if (threadIdx%x .eq. threadIdx%y) then
         source(threadIdx%x-1,(threadIdx%y-1) + bx*12 + by*12*nx_d + bz*12*nx_d*ny_d + bt*nx_d*ny_d*nz_d) = 1
      else
         source(threadIdx%x-1,(threadIdx%y-1) + bx*12 + by*12*nx_d + bz*12*nx_d*ny_d + bt*nx_d*ny_d*nz_d) = 0
      end if
   else
      source(threadIdx%x-1,(threadIdx%y-1) + bx*12 + by*12*nx_d + bz*12*nx_d*ny_d + bt*nx_d*ny_d*nz_d) = 0
   end if

end subroutine set_source
!=============================   

!=============================   
!Computes the (GPU) vector resulting from the action of the Dirac op. on a (GPU) vector
attributes(device) subroutine dirac_op_on_vector(v,u)
   complex(dp), intent(in) :: v
   complex(dp), intent(out) :: u

   

end subroutine
!=============================   


end module
