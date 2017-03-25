module partitioner 

use mpi
use settings
implicit none


! Data structures
type :: mesh_t
    double precision, dimension(:,:), allocatable :: x, y
end type mesh_t

type :: field_t
    double precision, dimension(:,:), allocatable :: u, v, P
    double precision, dimension(:,:), allocatable :: u_last, v_last
    double precision, dimension(:,:), allocatable :: u_star, v_star
    double precision, dimension(:,:), allocatable :: u_star2, v_star2
end type field_t

type :: forces_t
    double precision, dimension(:,:), allocatable :: F				! Pressure driven force
    double precision, dimension(:,:), allocatable :: Hx, Hy			! IBM force
	double precision, dimension(:,:), allocatable :: Bx, By			! Total body force (sum of all other forces)
end type forces_t


! global data variables
type(mesh_t) :: mesh
type(field_t) :: field
type(forces_t) :: forces

! Number of nodes in x-direction on this processor
integer :: m

! global node number of first vertex on this processor in x direction
integer :: iStartx

contains

subroutine setPartitions()

    integer :: Nproc, id, ierr

    call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    

    m = Nx/Nproc
    if (id<(Nx-m*Nproc)) then
        iStartx = id*(m+1)+1
        m = m+1
    else
        iStartx = id*m + (Nx-m*Nproc) +1
    end if

    allocate(mesh%x(Ny,m))
    allocate(mesh%y(Ny,m))

    allocate(field%u(Ny,0:m+1))
    allocate(field%v(Ny,0:m+1))
    allocate(field%P(Ny,0:m+1))

    !allocate(field%u_last(Ny,0:m+1))
    !allocate(field%v_last(Ny,0:m+1))
    
    allocate(field%u_star(Ny,0:m+1))
    allocate(field%v_star(Ny,0:m+1))
    
    allocate(field%u_star2(Ny,0:m+1))
    allocate(field%v_star2(Ny,0:m+1))

    allocate(forces%F(Ny,m))
    allocate(forces%Hx(Ny,m))
    allocate(forces%Hy(Ny,m))
	allocate(forces%Bx(Ny,m))
	allocate(forces%By(Ny,m))

    if (id==0) then
        write(*,'(A,I0,A)') 'Memory allocated on ', Nproc, ' processors'
    end if

end subroutine setPartitions


end module partitioner
