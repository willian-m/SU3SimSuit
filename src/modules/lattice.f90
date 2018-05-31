!==============================
!MODULE: lattice
!Contains variables and subroutines related to the lattice definition
!==============================
module lattice
use types_params
use math, only : gen_su3_element,SU3_dagger
use IO, only : read_lattice
implicit none

!==============================
!List of subroutines
!  -initialize_lattice: starts the lattice configuration
!  -init_hot: sets the lattice in a hot state
!  -init_cold: sets the lattice in a cold state
!  -init_aux_tables: sets the table of neighbours
!==============================

!==============================
!List of public variables
character(len=1024) :: lattice_file = '' !Path to file with lattice configuration for initialization
logical :: hot_start !If true, we will use a 'hot_start. Ortherwise, we use the cold start.'
integer, allocatable :: pos(:,:,:,:) !Array that maps the cartesian coordinates into the lattice index
integer, allocatable :: increment_table(:,:)

!Private variables
integer :: x,y,z,t,d
!==============================

private

!Public variables
public hot_start,increment_table,pos,lattice_file
!Public subroutines
public init_lattice

contains
!==============================
!Load the lattice configuration according to the parameters passed
subroutine init_lattice(lattice_file)
   character(len=1024), intent(in) :: lattice_file
   allocate(U(8,0:nx*ny*nz*nt-1))
   call init_aux_tables()

   if (lattice_file .eq. '') then
      if (hot_start) then
         call init_hot
      else
         call init_cold
      end if
   else
      call read_lattice(lattice_file)
   end if
end subroutine init_lattice
!==============================

!==============================
!Initialize the lattice in a hot state
subroutine init_hot
   integer :: y
   do x=0,nx*ny*nz*nt-1
      do d=2,8,2
         y = x + increment_table(x,d)
         call gen_su3_element(U(d,x))
         call SU3_dagger(U(d,x),U(d-1,y))
      end do
   end do
end subroutine
!==============================

!==============================
!Initialize the lattice in a cold state
subroutine init_cold
   type(SU3) :: V
   !call gen_su3_element(V)
   V%re = 0.0_dp
   V%im = 0.0_dp
   V%re(1,1) = 1.0_dp
   V%re(2,2) = 1.0_dp
   V%re(3,3) = 1.0_dp
   do x=0,nx*ny*nz*nt-1
      do d=1,8
         U(d,x) = V
      end do
   end do
end subroutine
!==============================

!==============================
!Initialize table of neighbours
subroutine init_aux_tables()
   integer :: i,j,k,l,x

   !We allocate a vector that maps the pos vector into an 1D array
   allocate(pos(nx,ny,nz,nt))
   do l=1,nt
      do k=1,nz
         do j=1,ny
            do i=1,nx
              pos(i,j,k,l) = (i-1) + nx*(j-1)+nx*ny*(k-1) + nx*ny*nz*(l-1)
            end do
         end do
      end do
   end do

   !For each one of the n points in direction x, we store the x coordinate
   !of the point that is m sites in front/behind it. Thus, -nx <= m <= nx
      
   allocate(increment_table(0:nx*ny*nz*nt-1,8))
      
   !Backaward at x direction
   do l=1,nt
      do k=1,nz
         do j=1,ny
            x = pos(1,j,k,l)
            increment_table(x,1) = nx-1
            do i=2,nx
               x = pos(i,j,k,l)
               increment_table(x,1) = -1
            end do
         end do
     end do
   end do

   !Forward at x direction
   do l=1,nt
      do k=1,nz
         do j=1,ny
            x = pos(nx,j,k,l)
            increment_table(x,2) = 1-nx
            do i=1,nx-1
               x = pos(i,j,k,l)
               increment_table(x,2) = 1
            end do
         end do
      end do
   end do

   !Backaward at y direction
   do l=1,nt
      do k=1,nz
         do i=1,nx
            x = pos(i,1,k,l)
            increment_table(x,3) = nx*(ny-1)
            do j=2,ny
               x = pos(i,j,k,l)
               increment_table(x,3) = -nx
            end do
         end do
      end do
   end do

   !Forward at y direction
   do l=1,nt
      do k=1,nz
         do i=1,nx
            x = pos(i,ny,k,l)
            increment_table(x,4) = nx*(1-ny)
            do j=1,ny-1
               x = pos(i,j,k,l)
               increment_table(x,4) = nx
            end do
         end do
      end do
   end do

   !Backaward at z direction
   do l=1,nt
      do j=1,ny
         do i=1,nx
            x = pos(i,j,1,l)
            increment_table(x,5) = nx*ny*(nz-1)
            do k=2,nz
               x = pos(i,j,k,l)
               increment_table(x,5) = -nx*ny
            end do
         end do
      end do
   end do

   !Forward at z direction
   do l=1,nt
      do j=1,ny
         do i=1,nx
            x = pos(i,j,nz,l)
            increment_table(x,6) = nx*ny*(1-nz)
            do k=1,nz-1
               x = pos(i,j,k,l)
               increment_table(x,6) = nx*ny
            end do
         end do
      end do
   end do

   !Backaward at t direction
   do k=1,nz
      do j=1,ny
         do i=1,nx
            x = pos(i,j,k,1)
            increment_table(x,7) = nx*ny*nz*(nt-1)
            do l=2,nt
               x = pos(i,j,k,l)
               increment_table(x,7) = -nx*ny*nz
            end do
         end do
      end do
   end do

   !Forward at t direction
   do k=1,nz
      do j=1,ny
         do i=1,nx
            x = pos(i,j,k,nt)
            increment_table(x,8) = nx*ny*nz*(1-nt)
            do l=1,nt-1
               x = pos(i,j,k,l)
               increment_table(x,8) = nx*ny*nz
            end do
         end do
      end do
   end do

end subroutine init_aux_tables
!==============================

end module 
