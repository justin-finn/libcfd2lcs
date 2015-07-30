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

end module analytic_velocity_m
