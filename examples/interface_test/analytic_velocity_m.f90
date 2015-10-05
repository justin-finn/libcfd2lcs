module analytic_velocity_m
	use user_data_m
	implicit none
	contains

	subroutine abc_velocity(ni,nj,nk,x,y,z,u,v,w,A,B,C,D,time)
		implicit none
		!-----
		integer,intent(in)::ni,nj,nk
		real(LCSRP),intent(in):: x(1:ni,1:nj,1:nk)
		real(LCSRP),intent(in):: y(1:ni,1:nj,1:nk)
		real(LCSRP),intent(in):: z(1:ni,1:nj,1:nk)
		real(LCSRP),intent(out):: u(1:ni,1:nj,1:nk)
		real(LCSRP),intent(out):: v(1:ni,1:nj,1:nk)
		real(LCSRP),intent(out):: w(1:ni,1:nj,1:nk)
		real(LCSRP),intent(in):: A,B,C,D
		real(LCSRP),intent(in):: time
		!-----
		integer:: i,j,k
		!-----
		! Assumes periodicity of multiple 2pi in X,Y,Z
		!-----

		do k =1,nk
		do j =1,nj
		do i =1,ni
			u(i,j,k) = (A+D*sin(time))*sin(z(i,j,k)) + C*cos(y(i,j,k))
			v(i,j,k) = B*sin(x(i,j,k)) + (A+D*sin(time))*cos(z(i,j,k))
			w(i,j,k) = C*sin(y(i,j,k)) + B*cos(x(i,j,k))
		enddo
		enddo
		enddo

	end subroutine abc_velocity

	subroutine double_gyre(ni,nj,nk,x,y,z,u,v,w,amplitude,eps,omega,time)
		implicit none
		!-----
		integer,intent(in)::ni,nj,nk
		real(LCSRP),intent(in):: x(1:ni,1:nj,1:nk)
		real(LCSRP),intent(in):: y(1:ni,1:nj,1:nk)
		real(LCSRP),intent(in):: z(1:ni,1:nj,1:nk)
		real(LCSRP),intent(out):: u(1:ni,1:nj,1:nk)
		real(LCSRP),intent(out):: v(1:ni,1:nj,1:nk)
		real(LCSRP),intent(out):: w(1:ni,1:nj,1:nk)
		real(LCSRP),intent(in)::  amplitude, eps,omega
		real(LCSRP),intent(in):: time
		!-----
		integer:: i,j,k
		real(LCSRP):: aoft,boft
		real(LCSRP):: fofxt(1:ni,1:nj,1:nk)
		real(LCSRP):: dfdx(1:ni,1:nj,1:nk)
		!-----
		! Assumes periodicity of multiple 2 in X, 1 in y, Z thickness is arbitrary
		!-----

		aoft = eps*sin(omega*time)
		boft = 1.0 -2.0*eps*sin(omega*time)
		do k =1,nk
		do j =1,nj
		do i =1,ni

			fofxt(i,j,k) = aoft*x(i,j,k)**2.0 + boft*x(i,j,k)
			dfdx(i,j,k) = 2.0*aoft*x(i,j,k) + boft

			u(i,j,k) = -PI*amplitude*sin(pi*fofxt(i,j,k))*cos(pi*y(i,j,k))
			v(i,j,k) = PI*amplitude*cos(pi*fofxt(i,j,k))*sin(pi*y(i,j,k))*dfdx(i,j,k)
			w(i,j,k) = 0.0
		enddo
		enddo
		enddo

	end subroutine double_gyre
end module analytic_velocity_m
