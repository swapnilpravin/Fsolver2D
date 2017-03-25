module schemes
use settings, only: nx, ny, dx, dy
use partitioner, only : m
use communicator

implicit none
    
contains
    subroutine CONV_U_UPWIND(D,u,v)
    
		double precision, dimension(ny,0:m+1) :: u,v
		double precision, dimension(ny-2,1:m) :: D 
		
		double precision, dimension(ny-2,1:m) :: u_plus,u_minus,v_plus,v_minus
		double precision, dimension(ny-2,1:m) :: Du2Dx,DuvDy
		
		u_plus=0;u_minus=0;v_plus=0
		v_minus=0;

		where (u(2:ny-1,1:m)>0) u_plus = 1; where (v(2:ny-1,1:m)>=0) v_plus = 1;
		where (u(2:ny-1,1:m)<0) u_minus = 1; where (v(2:ny-1,1:m)<0) v_minus = 1;

		Du2Dx = u_plus*(u(2:ny-1,1:m)*u(2:ny-1,1:m)-u(2:ny-1,0:m-1)*u(2:ny-1,0:m-1))/dx &
		    + u_minus*(u(2:ny-1,2:m+1)*u(2:ny-1,2:m+1)-u(2:ny-1,1:m)*u(2:ny-1,1:m))/dx
		DuvDy = v_plus*(u(2:ny-1,1:m)*v(2:ny-1,1:m)-u(1:ny-2,1:m)*v(1:ny-2,1:m))/dy &
		    + v_minus*(u(3:ny,1:m)*v(3:ny,1:m)-u(2:ny-1,1:m)*v(2:ny-1,1:m))/dy

		D = Du2Dx + DuvDy

    end subroutine CONV_U_UPWIND


    subroutine CONV_V_UPWIND(D,u,v)
    
		double precision, dimension(ny,0:m+1) :: u,v
		double precision, dimension(ny-2,1:m) :: D 
		
		double precision, dimension(ny-2,1:m) :: u_plus,u_minus,v_plus,v_minus
		double precision, dimension(ny-2,1:m) :: DvuDx,Dv2Dy
		
		u_plus=0;u_minus=0;v_plus=0
		v_minus=0;

		where (u(2:ny-1,1:m)>0) u_plus = 1; where (v(2:ny-1,1:m)>=0) v_plus = 1;
		where (u(2:ny-1,1:m)<0) u_minus = 1; where (v(2:ny-1,1:m)<0) v_minus = 1;

		DvuDx = u_plus*(v(2:ny-1,1:m)*u(2:ny-1,1:m)-v(2:ny-1,0:m-1)*u(2:ny-1,0:m-1))/dx &
		    + u_minus*(v(2:ny-1,2:m+1)*u(2:ny-1,2:m+1)-v(2:ny-1,1:m)*u(2:ny-1,1:m))/dx
		Dv2Dy = v_plus*(v(2:ny-1,1:m)*v(2:ny-1,1:m)-v(1:ny-2,1:m)*v(1:ny-2,1:m))/dy &
		    + v_minus*(v(3:ny,1:m)*v(3:ny,1:m)-v(2:ny-1,1:m)*v(2:ny-1,1:m))/dy

		D = DvuDx + Dv2Dy

    end subroutine CONV_V_UPWIND
    

	subroutine DEL2(D,T)

		double precision, dimension(ny,0:m+1) :: T
		double precision, dimension(ny-2,1:m) :: D
		
		double precision, dimension(ny-2,1:m) :: D2TDx2, D2TDy2
		
		D2TDx2 = (T(2:ny-1,2:m+1)-2*T(2:ny-1,1:m)+T(2:ny-1,0:m-1))/dx**2
		D2TDy2 = (T(3:ny,1:m)-2*T(2:ny-1,1:m)+T(1:ny-2,1:m))/dy**2

		D = D2TDx2 + D2TDy2

	end subroutine DEL2
	
	
	subroutine DIV(D,u,v)
	
		double precision, dimension(ny,0:m+1) :: u,v
		double precision, dimension(ny-2,1:m) :: D 
		
		double precision, dimension(ny-2,1:m) :: DuDx,DvDy

        ! Central difference

		! DuDx
		DuDx = (u(2:ny-1,2:m+1)-u(2:ny-1,0:m-1))/(2*dx)

		! DvDy
		DvDy = (v(3:ny,1:m)-v(1:ny-2,1:m))/(2*dy)


		D = DuDx + DvDy

	end subroutine DIV
	
	
    !------------------------------------------------------
    ! Divergence at all points including boundary nodes
    !------------------------------------------------------
    subroutine DIV_All(D,u,v)
	
		double precision, dimension(ny,0:m+1) :: u,v
		double precision, dimension(ny,0:m+1) :: D 
		
		double precision, dimension(ny,0:m+1) :: DuDx,DvDy

        ! Central difference at internal nodes and forward/backward difference
        ! at boundary nodes

		! DuDx
		DuDx(:,1:m) = (u(:,2:m+1)-u(:,0:m-1))/(2*dx)

		! DvDy
		DvDy(2:ny-1,1:m) = (v(3:ny,1:m)-v(1:ny-2,1:m))/(2*dy)
        DvDy(1,1:m) = (v(2,1:m)-v(1,1:m))/dy
        DvDy(ny,1:m) = (v(ny,1:m)-v(ny-1,1:m))/dy


		D = DuDx + DvDy

        call communicate(D)

	end subroutine DIV_All




	subroutine GRAD(DpDx,DpDy,P)
	
		double precision, dimension(ny,0:m+1) :: P
		double precision, dimension(ny-2,1:m) :: DpDx,DpDy

        ! Central difference

		DpDx = (P(2:ny-1,2:m+1)-P(2:ny-1,0:m-1))/(2*dx)
		DpDy = (P(3:ny,1:m)-P(1:ny-2,1:m))/(2*dy)

	end subroutine GRAD
	
	

	


end module schemes
