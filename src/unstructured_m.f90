module unstructured_m
	use data_m
	implicit none
	contains

	subroutine init_ur0(r0,n,label)
		implicit none
		!-----
		type(ur0_t):: r0
		integer:: n
		character(len=*):: label
		!-----
		if(lcsrank==0)&
			write(*,*) 'in init_ur0... ', trim(label)
		call destroy_ur0(r0)
		r0%n = n
		allocate(r0%r(1:n))
		r0%r = 0.0_LCSRP
		r0%label = label
	end subroutine init_ur0
	subroutine destroy_ur0(r0)
		implicit none
		!-----
		type(ur0_t):: r0
		!-----
		r0%n = 0
		if(allocated (r0%r)) deallocate(r0%r)
		r0%label = 'unused_ur0'
	end subroutine destroy_ur0

	subroutine init_ur1(r1,n,label)
		implicit none
		!-----
		type(ur1_t):: r1
		integer:: n
		character(len=*):: label
		!-----
		if(lcsrank==0)&
			write(*,*) 'in init_ur1... ', trim(label)
		call destroy_ur1(r1)
		r1%n = n
		allocate(r1%x(1:n))
		allocate(r1%y(1:n))
		allocate(r1%z(1:n))
		r1%x = 0.0_LCSRP
		r1%y = 0.0_LCSRP
		r1%z = 0.0_LCSRP
		r1%label = label
	end subroutine init_ur1
	subroutine destroy_ur1(r1)
		implicit none
		!-----
		type(ur1_t):: r1
		!-----
		r1%n = 0
		if(allocated (r1%x)) deallocate(r1%x)
		if(allocated (r1%y)) deallocate(r1%y)
		if(allocated (r1%z)) deallocate(r1%z)
		r1%label = 'unused_ur1'
	end subroutine destroy_ur1

	subroutine init_ur2(r2,n,label)
		implicit none
		!-----
		type(ur2_t):: r2
		integer:: n
		character(len=*):: label
		!-----
		if(lcsrank==0)&
			write(*,*) 'in init_ur2... ', trim(label)
		call destroy_ur2(r2)
		r2%n = n
		allocate(r2%xx(1:n))
		allocate(r2%xy(1:n))
		allocate(r2%xz(1:n))
		allocate(r2%yx(1:n))
		allocate(r2%yy(1:n))
		allocate(r2%yz(1:n))
		allocate(r2%zx(1:n))
		allocate(r2%zy(1:n))
		allocate(r2%zz(1:n))
		r2%xx = 0.0_LCSRP
		r2%xy = 0.0_LCSRP
		r2%xz = 0.0_LCSRP
		r2%yx = 0.0_LCSRP
		r2%yy = 0.0_LCSRP
		r2%yz = 0.0_LCSRP
		r2%zx = 0.0_LCSRP
		r2%zy = 0.0_LCSRP
		r2%zz = 0.0_LCSRP
		r2%label = label
	end subroutine init_ur2
	subroutine destroy_ur2(r2)
		implicit none
		!-----
		type(ur2_t):: r2
		!-----
		r2%n = 0
		if(allocated (r2%xx)) deallocate(r2%xx)
		if(allocated (r2%xy)) deallocate(r2%xy)
		if(allocated (r2%xz)) deallocate(r2%xz)
		if(allocated (r2%yx)) deallocate(r2%yx)
		if(allocated (r2%yy)) deallocate(r2%yy)
		if(allocated (r2%yz)) deallocate(r2%yz)
		if(allocated (r2%zx)) deallocate(r2%zx)
		if(allocated (r2%zy)) deallocate(r2%zy)
		if(allocated (r2%zz)) deallocate(r2%zz)
		r2%label = 'unused_ur2'
	end subroutine destroy_ur2

end module unstructured_m
