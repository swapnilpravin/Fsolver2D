module IBM

use settings, only : ny,nx,Ly,Lx,RADIUS,dt,rho
use partitioner, only : m

contains

	subroutine IBMforce(Hx,Hy,u,v,x,y)
	
	end subroutine IBMforce


	subroutine IBMforceCircle(Hx,Hy,u,v,x,y)

		double precision, dimension(ny,m) :: Hx,Hy,x,y
		double precision, dimension(ny,0:m+1) :: u,v

		integer, dimension(ny,m) :: eta
		double precision :: yc,xc

		yc = Ly/2
		xc = Lx/2

		eta = 0
		where ((x-xc)**2 + (y-yc)**2 < RADIUS**2 ) eta=1

		Hx = -eta*u(:,1:m)/(dt)
		Hy = -eta*v(:,1:m)/(dt)

	end subroutine IBMforceCircle
	

end module IBM
