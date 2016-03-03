module geometry_m
	use data_m
	implicit none
	contains

	!Return a unit vector cross product of a,b:
	!Optionally, make the result a unit vector.
	function cross_product(a, b, unitvector)
		implicit none
		real (LCSRP) :: cross_product(3)
		real (LCSRP), intent(in) :: a(3), b(3)
		logical,optional:: unitvector
		real (LCSRP):: c(3),mag
		c(1) = a(2)*b(3) - a(3)*b(2)
		c(2) = a(3)*b(1) - a(1)*b(3)
		c(3) = a(1)*b(2) - a(2)*b(1)

		if(present(unitvector) .and. unitvector) then
			mag = sqrt(c(1)*c(1)+c(2)*c(2)+c(3)*c(3))
			cross_product(1) = c(1)/mag
			cross_product(2) = c(2)/mag
			cross_product(3) = c(3)/mag
		else
			cross_product = c
		endif
	end function cross_product

	function grid_refine(v0,v1,v2,v3,v4,v5,v6,v7,s,t,u)
		implicit none
		!-----
		real(LCSRP):: grid_refine(3)
		real(LCSRP):: v0(3),v1(3),v2(3),v3(3),v4(3),v5(3),v6(3),v7(3)
		real(LCSRP):: s,t,u
		!-----
		real(LCSRP):: va(3), vb(3), vc(3) !vertices along edges of hex.
		real(LCSRP):: n1(3), n2(3), n3(3) !normal of plane	
		real(LCSRP):: d1, d2, d3 !point dotted with normal	
		real(LCSRP):: x(3)  !point of intersection
		!-----

		!First plane
		va = v0 + s*(v1-v0)
		vb = v4 + s*(v5-v4)
		vc = v2 + s*(v3-v2)
		n1 = cross_product( vb-va,vc-va, unitvector=.true.)
		d1 = dot_product(n1,va)

		!second plane
		va = v0 + t*(v2-v0)
		vb = v1 + t*(v3-v1)
		vc = v4 + t*(v6-v4)
		n2 = cross_product(vb-va,vc-va, unitvector=.true.)
		d2 = dot_product(n2,va)
		
		!third plane
		va = v0 + u*(v4-v0)
		vb = v1 + u*(v5-v1)
		vc = v2 + u*(v6-v2)
		n3 = cross_product(vb-va,vc-va, unitvector=.true.)
		d3 = dot_product(n3,va)
		
		x = (d1*cross_product(n2,n3) + d2*cross_product(n3,n1) &
			+d3*cross_product(n1,n2)) /	(dot_product(n1,cross_product(n2,n3)))
		

!if(lcsrank==0) write(*,*) '******************************'
!if(lcsrank==0) write(*,*) 'v0',v0
!if(lcsrank==0) write(*,*) 'v1',v1
!if(lcsrank==0) write(*,*) 'v2',v2
!if(lcsrank==0) write(*,*) 'v3',v3
!if(lcsrank==0) write(*,*) 'v4',v4
!if(lcsrank==0) write(*,*) 'v5',v5
!if(lcsrank==0) write(*,*) 'v6',v6
!if(lcsrank==0) write(*,*) 'v7',v7
!if(lcsrank==0) write(*,*) 'x',x
!if(lcsrank==0) write(*,*) 'Delta',x(1) - v0(1), x(2)-v0(2),x(3)-v0(3),s,t,u
!if(lcsrank==0) write(*,*) '******************************'
		grid_refine = x
	end function grid_refine

	


end module geometry_m
