module communicator

use mpi
use settings, only: ny
use partitioner, only: m

implicit none

contains

subroutine communicate(P)
    implicit none
    double precision, dimension(ny,0:m+1) :: P
    integer :: id, Nproc, ierr
    integer :: n        ! size of data to be send/received
    integer :: left_flap_tag, right_flap_tag ! mpi tags for data meant for right and left ghost flaps on each processor

    left_flap_tag = 0
    right_flap_tag = 1
    
    n = ny

    !print*, n

    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)

    call mpi_send(P(:,m),n,MPI_DOUBLE_PRECISION,mod(id+1,Nproc),left_flap_tag,MPI_COMM_WORLD,ierr)
    if (id /= 0) call mpi_send(P(:,1),n,MPI_DOUBLE_PRECISION,id-1,right_flap_tag,MPI_COMM_WORLD,ierr)

    if (id==0) then
        call mpi_recv(P(:,0),n,MPI_DOUBLE_PRECISION,Nproc-1,left_flap_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        call mpi_send(P(:,1),n,MPI_DOUBLE_PRECISION,Nproc-1,right_flap_tag,MPI_COMM_WORLD,ierr)
    end if

    if (id /= 0) call mpi_recv(P(:,0),n,MPI_DOUBLE_PRECISION,id-1,left_flap_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call mpi_recv(P(:,m+1),n,MPI_DOUBLE_PRECISION,mod(id+1,Nproc),right_flap_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)


end subroutine communicate


end module communicator
