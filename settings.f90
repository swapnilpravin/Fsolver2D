module settings

use mpi

    implicit none
    
    double precision, parameter :: pi = 4*atan(1.0)
    
    integer :: nx,ny
    double precision :: dx,dy,dt
    double precision :: rho,mu,nu
    integer :: nt,nit
    double precision :: omega_P, omega_U, omega_V, tol
    double precision :: Lx,Ly
    double precision :: U0
    integer:: nLog, nWrite		! nLog: write log to terminal every nLog timesteps
    							! write: results to file every nwrite timesteps
    double precision :: RADIUS	! for testing: circle radius
    
    
    
    contains
        subroutine setup(filename)

	        character(len=*) :: filename
	        character(len=30) :: temp

            integer :: id, Nproc, ierr

            call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)

            if (id==0) then
	        
                open(unit=10,file=filename,status='old',action='read')
                
                read (10,*) temp, Lx
                read (10,*) temp, Ly
                read (10,*) temp, dx
                read (10,*) temp, dy
                read (10,*) temp, dt
                read (10,*) temp, nt
                read (10,*) temp, nit
                read (10,*) temp, rho
                read (10,*) temp, nu
                read (10,*) temp, omega_P
                read (10,*) temp, omega_U
                read (10,*) temp, omega_V
                read (10,*) temp, tol
                read (10,*) temp, U0
                read (10,*) temp, nLog
                read (10,*) temp, nWrite
                read (10,*) temp, RADIUS
                
			    close(10)

            end if

            ! broadcast to all
            call mpi_bcast(Lx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(Ly,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(dx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(dy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(nt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(nit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(rho,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(nu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(omega_P,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(omega_U,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(omega_V,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(U0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(nLog,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(nWrite,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

            !call mpi_barrier(MPI_COMM_WORLD,ierr)

			nx = floor(Lx/dx)
			ny = floor(Ly/dy)

            !print*, Lx, Ly, Lz
			
			mu = nu*rho
			
            
        end subroutine setup
        
        
        
        subroutine setTimestep(a)
        	double precision :: a
        	dt = a
        end subroutine
        

end module settings
