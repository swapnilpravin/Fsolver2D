subroutine channel()
    use settings
    use partitioner
	use io
	use schemes
	use solver
	use boundaryConditions
	use IBM

	implicit none

	
    !double precision ,dimension(Np,3) :: Yp, Vp, Vp0		! Particle variables

	double precision ,dimension(ny-2,1:m) :: D_CONV_U, D_CONV_V, DEL2U, DEL2V, DpDx, DpDy
	
    double precision :: R, E_P, E_U
	integer :: NITS_final_P, NITS_final_U

	integer :: i, j

    integer :: id, Nproc, ierr

    ! for io
    double precision :: Umax, Pmax, Umax_loc, Pmax_loc

    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)

    
	! Initialize Flow Variables
	field%u=0; field%v=0; field%P=0
	!field%u_last=0; field%v_last=0;
	field%u_star=0; field%v_star=0
	field%u_star2=0; field%v_star2=0

    !if(id==0) print*, 'flow initialized.'

	! Pressure driven force
	forces%F = 8*nu*U0/Ly**2
    !forces%F = 0

	forces%Hx=0; forces%Hy=0
	forces%Bx=0; forces%By=0;

	! construct x,y,z matrices (mesh)
	do i=1,m; mesh%x(:,i) = (iStartx+i-1)*dx;	end do
	do i=1,ny; mesh%y(i,:) = (i-1)*dy;	end do

	! Velocity IC
	call setIC(field%u,field%v,field%P,mesh%x,mesh%y)

	![u,v] = setVelocityBC1(u,v,dx,dy,dt,nu);

	! write data before starting simulation
	call writeToTecplot2D_MPI(mesh%x,mesh%y,field%u,field%v,field%P,0)
	
	do i=1,nt
			
		!call CONV_U_UPWIND(D_CONV_U,field%u,field%v)
		!call CONV_V_UPWIND(D_CONV_V,field%u,field%v)
	
		!call DEL2(DEL2U,field%u)
		!call DEL2(DEL2V,field%v)

		!field%u_star(2:ny-1,1:m) = field%u(2:ny-1,1:m) + dt*(-D_CONV_U + nu*DEL2U)
		!field%v_star(2:ny-1,1:m) = field%v(2:ny-1,1:m) + dt*(-D_CONV_V + nu*DEL2V)

		! call setVelocityBC_MPI(u,v,w)
		
		!call calcUstar_AB(field%u_star,field%v_star,field%u,field%v,field%u_last,field%v_last)
		
		! Total body force
		forces%Bx = forces%F + forces%Hx
		forces%By = forces%Hy
		
		call calcUstar_Implicit(field%u_star,field%v_star,E_U,NITS_final_U,field%u,field%v,forces%Bx,forces%By)

        !field%u_last = field%u
        !field%v_last = field%v

        ! COMMUNICATE (u_star)
        !call communicate(field%u_star)
        !call communicate(field%v_star)

        !if(id==0) print*, 'communication done'
        
        !call setVelocityBC_MPI(field%u_star,field%v_star)

		call PressureSolver_MPI(field%P,E_P,NITS_final_P,field%u_star,field%v_star)
		!call PressureSolver_MPI_BiCGStab(field%P,E_P,NITS_final_P,field%u_star,field%v_star)

        !if(id==0) print*, 'Pressure solver done: ', NITS_final
		call GRAD(DpDx,DpDy,field%P)

        !if(id==0) print*, 'GRAD calculated'

		field%u_star2(2:ny-1,1:m) = field%u_star(2:ny-1,1:m) - dt/rho*DpDx
		field%v_star2(2:ny-1,1:m) = field%v_star(2:ny-1,1:m) - dt/rho*DpDy

        call setVelocityBC_MPI(field%u_star2,field%v_star2)

		call communicate(field%u_star2)
		call communicate(field%v_star2)


        !if(id==0) print*, 'u_star2 calculated'

		! call setVelocityBC2(u_star2,v_star2,w_star2,u_star,v_star,w_star,P)

		! Calculate force for IBM
		!call IBMForceCircle(forces%Hx,forces%Hy,field%u,field%v,mesh%x,mesh%y)

		! march in time using force
		! (Entire x dimension is used, since periodic BC in that direction)
		field%u = field%u_star2
		field%v = field%v_star2

        !call setVelocityBC_MPI(field%u,field%v)
        
        ! COMMUNICATE (u)
        !call communicate(field%u)
        !call communicate(field%v)


        !if(id==0) print*, 'field updated'

        Umax_loc = maxval(field%u(:,1:m))
        Pmax_loc = maxval(field%P(:,1:m))
        call mpi_reduce(Umax_loc,Umax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call mpi_reduce(Pmax_loc,Pmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        !if(id==0) print*, 'reduction for io done'

        if (id==0) then       
            if (mod(i,nLog)==0) then
                call writeLogToTerminal_MPI(i, &
                'Umax','double',Umax, &
                'Pmax','double',Pmax, &
                'Pressure Error','double',E_P, &
                'Pressure Iterations','integer',dble(NITS_final_P), &
				'Velocity Error','double',E_U, &
				'Velocity Iterations','integer',dble(NITS_final_U) )
            end if
            
            
        end if

        if (mod(i,nWrite)==0) then
            call writeToTecplot2D_MPI(mesh%x,mesh%y,field%u,field%v,field%P,i)
        end if

        call mpi_barrier(MPI_COMM_WORLD,ierr)
		
	end do


end subroutine channel
