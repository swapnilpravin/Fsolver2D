module solver

use mpi
use schemes
use boundaryConditions
use communicator

implicit none

contains
	
	
    subroutine PressureSolver_MPI(P,E_rms,NITS_final,u,v)
        use settings, only: nx, ny, dx, dy, dt, rho, nit, omega_P, tol

        double precision, dimension(ny,0:m+1) :: u,v,P
        double precision :: E_rms, E_sqr_loc, E_sqr
        integer :: NITS_final
        
        double precision, dimension(ny-2,1:m) :: DIVU, b
        double precision, dimension(ny,1:m) :: E
        double precision, dimension(ny,0:m+1) :: P_new
        double precision :: ap,ae,aw,an,as,P_max
        integer :: i,n

        integer :: ierr

        call DIV(DIVU,u,v)

		ap = 2/dx**2 + 2/dy**2
		ae = 1/dx**2;
		aw = 1/dx**2;
		an = 1/dy**2;
		as = 1/dy**2;
		b = -rho*(DIVU)/dt; ! dt affects possibility of convergence

		n = nx*ny

		!P_new = P

		do i=1,nit
			!print*, i
			P_new(2:ny-1,1:m) = (an*P(3:ny,1:m)+as*P(1:ny-2,1:m) &
				+ae*P(2:ny-1,2:m+1)+aw*P(2:ny-1,0:m-1) + b)/ap

			call setPressureBC_MPI(P_new)

            ! COMMUNICATE(P)
            call communicate(P_new)

			if (mod(i,50)==0) then

				E = P_new(:,1:m) - P(:,1:m)
				E_sqr_loc = sum(E**2)

                call mpi_allreduce(E_sqr_loc,E_sqr,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
                E_rms = sqrt(E_sqr/n)

				if (E_rms<tol) then
					exit
				end if

			end if

            ! MPI_Barrier ??

			P = P + omega_P * (P_new - P)
			
			call communicate(P)

		end do
		
		NITS_final = i
		
		! R = ap*P(2:nz-1,2:ny-1,2:nx-1)-(an*P(2:nz-1,3:ny,2:nx-1)+as*P(2:nz-1,1:ny-2,2:nx-1) &
		!	+ae*P(2:nz-1,2:ny-1,3:nx)+aw*P(2:nz-1,2:ny-1,1:nx-2) &
		!	+at*P(3:nz,2:ny-1,2:nx-1)+ab*P(1:nz-2,2:ny-1,2:nx-1)+b);
		! R_rms = sqrt(sum(R**2)/n)



    end subroutine PressureSolver_MPI



    !----------------------------------------------
    ! PressureSolver_MPI_BiCGStab
    !----------------------------------------------    
	subroutine PressureSolver_MPI_BiCGStab(X,E_rms,NITS_final,u,v)
        use settings, only: nx, ny, dx, dy, dt, rho, nit, tol
        double precision, dimension(ny,0:m+1) :: u,v,X
        double precision :: E_rms, E_sqr_loc, E_sqr
        integer :: NITS_final
        
        double precision, dimension(ny,0:m+1) :: DIVU, b ! includes boundary nodes and flaps
        double precision, dimension(ny,1:m) :: E
        double precision, dimension(ny,0:m+1) :: X_new
        double precision :: ap,ae,aw,an,as,X_max
        integer :: i,n

        ! BiCGStab Variables
        double precision, dimension(ny,0:m+1) :: r, r0_star, p, s, r_new, p_new,  AX_temp
        double precision :: temp_loc, temp1, temp2
        double precision :: alpha, beta, omega


        integer :: id, Nproc, ierr

        call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
        call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)

		!call setPressureBC_MPI(X)
        !call communicate(X)
        call DIV_All(DIVU,u,v)

		!ap = 2/dx**2 + 2/dy**2
		!ae = 1/dx**2;
		!aw = 1/dx**2;
		!an = 1/dy**2;
		!as = 1/dy**2;
		b = -rho*(DIVU)/dt; ! dt affects possibility of convergence

        ! for BC's
        b(1,:) = 0; b(ny,:) = 0 ! top and bottom
        ! for cavity, left and right boundaries
        !if (id==0) b(:,1) = 0
        !if (id==Nproc-1) b(:,m) = 0


		n = nx*ny

		!P_new = P

        r = b - AX(X)
        r0_star = r
        p = r
		do i=1,nit
			!print*, i
            
            AX_temp = AX(p)
            temp_loc = dot_product(pack(r(:,1:m),.true.),pack(r0_star(:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            temp_loc = dot_product(pack(AX_temp(:,1:m),.true.),pack(r0_star(:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            alpha = temp1/temp2

            !print*, temp_loc, temp1, temp2, alpha
            !if(i==50) stop

            s = r - alpha * AX(p)
            call communicate(s)

            AX_temp = AX(s)
            temp_loc = dot_product(pack(AX_temp(:,1:m),.true.),pack(s(:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            temp_loc = dot_product(pack(AX_temp(:,1:m),.true.),pack(AX_temp(:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            omega = temp1/temp2

            ! update X
            X_new = X + alpha*p + omega*s
            call communicate(X_new)

            ! update r
            r_new = s - omega*AX(s)
            call communicate(r_new)

            temp_loc = dot_product(pack(r_new(:,1:m),.true.),pack(r0_star(:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            temp_loc = dot_product(pack(r(:,1:m),.true.),pack(r0_star(:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            beta = (temp1/temp2)*(alpha/omega)

            ! update p
            p_new = r_new + beta * (p-omega*AX(p))
            call communicate(p_new)




            ! COMMUNICATE(P)
            !call communicate(X_new)

			!call setPressureBC_MPI(X_new)


			if (mod(i,1)==0) then

				E = X_new(:,1:m) - X(:,1:m)
				E_sqr_loc = sum(E**2)

                call mpi_allreduce(E_sqr_loc,E_sqr,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
                E_rms = sqrt(E_sqr/n)

                !print*, 'Pressure Error:',E_rms

				if (E_rms<tol) then
					exit
				end if

			end if

            ! MPI_Barrier ??

			!X = X_new
            !r = r_new
            !p = p_new

			X = X + omega_P*(X_new-X)
			r = r + omega_P*(r_new-r)
			p = p + omega_P*(p_new-p)

		end do
		
		NITS_final = i
		
		! R = ap*P(2:nz-1,2:ny-1,2:nx-1)-(an*P(2:nz-1,3:ny,2:nx-1)+as*P(2:nz-1,1:ny-2,2:nx-1) &
		!	+ae*P(2:nz-1,2:ny-1,3:nx)+aw*P(2:nz-1,2:ny-1,1:nx-2) &
		!	+at*P(3:nz,2:ny-1,2:nx-1)+ab*P(1:nz-2,2:ny-1,2:nx-1)+b);
		! R_rms = sqrt(sum(R**2)/n)



    end subroutine PressureSolver_MPI_BiCGStab



    function AX(X)
        use settings, only: nx, ny, dx, dy, dt

        double precision, dimension(ny,0:m+1) :: X
        double precision, dimension(ny,0:m+1) :: AX
        double precision :: ap, ae, aw, an, as

        integer :: id, Nproc, ierr

        call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
        call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)

    
        ap = 2/dx**2 + 2/dy**2
		ae = 1/dx**2;
		aw = 1/dx**2;
		an = 1/dy**2;
		as = 1/dy**2;

        AX(2:ny-1,1:m) = ap*X(2:ny-1,1:m) - ( an*X(3:ny,1:m) + as*X(1:ny-2,1:m) + aw*X(2:ny-1,0:m-1) + ae*X(2:ny-1,2:m+1) )
         AX(1,1:m) = X(1,1:m) - X(2,1:m)
        AX(ny,1:m) = X(ny,1:m) - X(ny-1,1:m)

		AX(1,:) = X(1,:) - X(2,:)
        AX(ny,:) = X(ny,:) - X(ny-1,:)

        ! for cavity case, bc at left and right edges
        !if (id==0) AX(:,1) = X(:,1) - X(:,2)
        !if (id==Nproc-1) AX(:,m) = X(:,m) - X(:,m-1)

        call communicate(AX)

    end function AX



    subroutine calcUstar_AB(u_star,v_star,u,v,u_last,v_last)
        use settings, only: nx, ny, dx, dy, dt, rho, nit, tol

        double precision, dimension(ny,0:m+1) :: u,v,u_star,v_star,u_last,v_last
        double precision, dimension(2:ny-1,1:m) :: SU,SV,SU_last,SV_last

        double precision, dimension(2:ny-1,1:m) :: D_CONV_U,D_CONV_V,D_CONV_U_last, D_CONV_V_last
        double precision, dimension(2:ny-1,1:m) :: DEL2U,DEL2V,DEL2U_last, DEL2V_last

        call CONV_U_UPWIND(D_CONV_U,u,v)
        call CONV_V_UPWIND(D_CONV_V,u,v)
        call CONV_U_UPWIND(D_CONV_U_last,u_last,v_last)
        call CONV_V_UPWIND(D_CONV_V_last,u_last,v_last)

        call DEL2(DEL2U,u)
        call DEL2(DEL2V,v)
        call DEL2(DEL2U_last,u_last)
        call DEL2(DEL2V_last,v_last)

        SU = -D_CONV_U + nu*DEL2U
        SV = -D_CONV_V + nu*DEL2V
        SU_last = -D_CONV_U_last + nu*DEL2U_last
        SV_last = -D_CONV_V_last + nu*DEL2V_last

        u_star(2:ny-1,1:m) = u(2:ny-1,1:m) + dt * (1.5*SU-0.5*SU_last)
        v_star(2:ny-1,1:m) = v(2:ny-1,1:m) + dt * (1.5*SV-0.5*SV_last)

        call setVelocityBC_MPI(u,v)

    end subroutine
	
	
	
	subroutine calcUstar_Implicit(u_star,v_star,E_rms,NITS_final,u,v,Bx,By)
	
		double precision, dimension(ny,0:m+1) :: u, v, u_star,v_star
		double precision, dimension(ny,0:m+1) :: u_star_new,v_star_new
		
		double precision ,dimension(ny-2,1:m) :: D_CONV_U, D_CONV_V, DEL2U, DEL2V
		
		double precision ,dimension(ny,1:m) :: Bx, By
		
		double precision, dimension(ny,1:m) :: E
		double precision :: E_rms, E_sqr_loc, E_sqr
        integer :: NITS_final
		integer :: i, n, Nproc, ierr, id

		call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
        call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)
		
		
		n = ny*nx
		
		u_star = u
		v_star = v
		
		u_star_new = 0
		v_star_new = 0
		
		do i=1,nit
			call CONV_U_UPWIND(D_CONV_U,u_star,v_star)
			call CONV_V_UPWIND(D_CONV_V,u_star,v_star)
	
			call DEL2(DEL2U,u_star)
			call DEL2(DEL2V,v_star)

			u_star_new(2:ny-1,1:m) = u(2:ny-1,1:m) + dt*(-D_CONV_U + nu*DEL2U + Bx(2:ny-1,:))
			v_star_new(2:ny-1,1:m) = v(2:ny-1,1:m) + dt*(-D_CONV_V + nu*DEL2V + By(2:ny-1,:))

			call setVelocityBC_MPI(u_star_new,v_star_new)
			
			!if (id==0) print*, id, v_star_new(1:5,0)
			!if (id==3) print*, id, v_star_new(1:5,m)

			! COMMUNICATE (u_star)
			call communicate(u_star_new)
			call communicate(v_star_new)
			
			!if (id==0) print*, id, v_star_new(1:5,0)
			!if (id==3) print*, id, v_star_new(1:5,m)

			!if (id==0) read(*,*)

			E = u_star_new(:,1:m) - u_star(:,1:m)
			E_sqr_loc = sum(E**2)

			call mpi_allreduce(E_sqr_loc,E_sqr,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
			E_rms = sqrt(E_sqr/n)

			!if (id==0) print*, 'U_star error:', E_rms
			
			if (E_rms<tol) then
				exit
			end if
			
			u_star = u_star + omega_U * (u_star_new - u_star)
			v_star = v_star + omega_V * (v_star_new - v_star)

			call communicate(u_star)
			call communicate(v_star)

			
		end do
		
		NITS_final = i
	
	end subroutine calcUstar_Implicit
	

end module solver
