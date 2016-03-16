module unstructured_m
	use data_m
	implicit none
	!JRF:  March 14, 2016:  A few changes here to make sure we preserve the label on unstrucured arrays,
	!when they have n = 0 for size. This can happen for lp tracer sets that do not span all processors.
	integer,parameter:: MEM_INC = 1024

	contains

	subroutine init_ur0(r0,n,label)
		implicit none
		!-----
		type(ur0_t):: r0
		integer:: n
		character(len=*):: label
		!-----
		integer:: s
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in init_ur0... ', trim(label)
		call destroy_ur0(r0)
		r0%n = n
		s = n+MEM_INC
		allocate(r0%r(1:s))
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
	subroutine resize_ur0(r0,n)
		implicit none
		!-----
		type(ur0_t):: r0
		integer:: n
		!-----
		real(LCSRP),allocatable:: tmp(:)
		integer:: s
		character(len=LCS_NAMELEN):: label
		!-----
		!if(LCS_VERBOSE)&
		!	write(*,*) lcsrank,'in resize_ur0...',trim(r0%label)

		label= trim(r0%label)

		!check that this data initialized
		if(r0%n==0) then
			call init_ur0(r0,n,r0%label)
		endif
		!deallocate if n <=0
		if (n<=0) then
			!call destroy_ur0(r0)  !JRF:  Dont allocate
			!return
		endif
		!current size ok
		if( n <= size(r0%r) .AND. size(r0%r) - MEM_INC <= n)then
			!if(LCS_VERBOSE) write(*,*) lcsrank,'size ok'
			r0%n = n
			!return
		endif
		!increase size and set new n
		if (n > size(r0%r)) then
			allocate(tmp(1:r0%n))
			tmp(1:r0%n) = r0%r(1:r0%n)
			deallocate(r0%r)
			s = n+MEM_INC
			allocate(r0%r(1:s))
			r0%r(1:r0%n) = tmp(1:r0%n)
			r0%r(r0%n+1:s) = 0.0_LCSRP
			r0%n = n
			deallocate(tmp)
		endif
		!reduce the size and set new n
		if(size(r0%r)-MEM_INC > n) then
			allocate(tmp(1:n))
			tmp(1:n) = r0%r(1:n)
			deallocate(r0%r)
			s = n+MEM_INC/2
			allocate(r0%r(1:s))
			r0%r(1:n) = tmp(1:n)
			r0%r(n+1:s) = 0.0_LCSRP
			r0%n = n
			deallocate(tmp)
		endif
		r0%label = trim(label)
	end subroutine resize_ur0

	subroutine init_ur1(r1,n,label)
		implicit none
		!-----
		type(ur1_t):: r1
		integer:: n
		character(len=*):: label
		!-----
		integer:: s
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in init_ur1... ', trim(label)
		call destroy_ur1(r1)
		r1%n = n
		s = n+MEM_INC
		allocate(r1%x(1:s))
		allocate(r1%y(1:s))
		allocate(r1%z(1:s))
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
	subroutine resize_ur1(r1,n)
		implicit none
		!-----
		type(ur1_t):: r1
		integer:: n
		!-----
		real(LCSRP),allocatable:: tmp(:)
		integer:: s
		character(len=LCS_NAMELEN):: label
		!-----
		!if(LCS_VERBOSE)&
		!	write(*,*) lcsrank,'in resize_ur1...',trim(r1%label)

		label= trim(r1%label)

		!check that this data initialized
		if(r1%n==0) then
			call init_ur1(r1,n,r1%label)
		endif
		!deallocate if n <=0
		if (n<=0) then
			!call destroy_ur1(r1)  !JRF: Don't destroy
			!return
		endif
		!current size ok
		if( n <= size(r1%x) .AND. size(r1%x) - MEM_INC <= n)then
			!if(LCS_VERBOSE) write(*,*) lcsrank,'size ok'
			r1%n = n
			!return
		endif
		!increase size and set new n
		if (n > size(r1%x)) then
			s = n+MEM_INC
			allocate(tmp(1:r1%n))

			!x
			tmp(1:r1%n) = r1%x(1:r1%n)
			deallocate(r1%x)
			allocate(r1%x(1:s))
			r1%x(1:r1%n) = tmp(1:r1%n)
			r1%x(r1%n+1:s) = 0.0_LCSRP
			!y
			tmp(1:r1%n) = r1%y(1:r1%n)
			deallocate(r1%y)
			allocate(r1%y(1:s))
			r1%y(1:r1%n) = tmp(1:r1%n)
			r1%y(r1%n+1:s) = 0.0_LCSRP
			!z
			tmp(1:r1%n) = r1%z(1:r1%n)
			deallocate(r1%z)
			allocate(r1%z(1:s))
			r1%z(1:r1%n) = tmp(1:r1%n)
			r1%z(r1%n+1:s) = 0.0_LCSRP

			r1%n = n
			deallocate(tmp)
		endif
		!reduce the size and set new n
		if(size(r1%x)-MEM_INC > n) then
			s = n+MEM_INC/2
			allocate(tmp(1:n))

			!x
			tmp(1:n) = r1%x(1:n)
			deallocate(r1%x)
			allocate(r1%x(1:s))
			r1%x(1:n) = tmp(1:n)
			r1%x(n+1:s) = 0.0_LCSRP
			!y
			tmp(1:n) = r1%y(1:n)
			deallocate(r1%y)
			allocate(r1%y(1:s))
			r1%y(1:n) = tmp(1:n)
			r1%y(n+1:s) = 0.0_LCSRP
			!z
			tmp(1:n) = r1%z(1:n)
			deallocate(r1%z)
			allocate(r1%z(1:s))
			r1%z(1:n) = tmp(1:n)
			r1%z(n+1:s) = 0.0_LCSRP

			r1%n = n
			deallocate(tmp)
		endif
		r1%label = trim(label)
	end subroutine resize_ur1


	subroutine init_ur2(r2,n,label)
		implicit none
		!-----
		type(ur2_t):: r2
		integer:: n
		character(len=*):: label
		!-----
		integer:: s
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in init_ur2... ', trim(label)
		call destroy_ur2(r2)
		r2%n = n
		s = n+MEM_INC
		allocate(r2%xx(1:s))
		allocate(r2%xy(1:s))
		allocate(r2%xz(1:s))
		allocate(r2%yx(1:s))
		allocate(r2%yy(1:s))
		allocate(r2%yz(1:s))
		allocate(r2%zx(1:s))
		allocate(r2%zy(1:s))
		allocate(r2%zz(1:s))
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
	subroutine resize_ur2(r2,n)
		implicit none
		!-----
		type(ur2_t):: r2
		integer:: n
		!-----
		real(LCSRP),allocatable:: tmp(:)
		integer:: s
		character(len=LCS_NAMELEN):: label
		!-----
		!if(LCS_VERBOSE)&
		!	write(*,*) 'in resize_ur2...',trim(r2%label)
		label= trim(r2%label)

		!check that this data initialized
		if(r2%n==0) then
			call init_ur2(r2,n,r2%label)
		endif
		!deallocate if n <=0
		if (n<=0) then
			!call destroy_ur2(r2)  !JRF Don't destroy
			!return
		endif
		!current size ok
		if( n <= size(r2%xx) .AND. size(r2%xx) - MEM_INC <= n)then
			!if(LCS_VERBOSE) write(*,*) lcsrank,'size ok'
			r2%n = n
			!return
		endif
		!increase size and set new n
		if (n > size(r2%xx)) then
			s = n+MEM_INC
			allocate(tmp(1:r2%n))

			!xx
			tmp(1:r2%n) = r2%xx(1:r2%n)
			deallocate(r2%xx)
			allocate(r2%xx(1:s))
			r2%xx(1:r2%n) = tmp(1:r2%n)
			r2%xx(r2%n+1:s) = 0.0_LCSRP
			!xy
			tmp(1:r2%n) = r2%xy(1:r2%n)
			deallocate(r2%xy)
			allocate(r2%xy(1:s))
			r2%xy(1:r2%n) = tmp(1:r2%n)
			r2%xy(r2%n+1:s) = 0.0_LCSRP
			!xz
			tmp(1:r2%n) = r2%xz(1:r2%n)
			deallocate(r2%xz)
			allocate(r2%xz(1:s))
			r2%xz(1:r2%n) = tmp(1:r2%n)
			r2%xz(r2%n+1:s) = 0.0_LCSRP
			!yx
			tmp(1:r2%n) = r2%yx(1:r2%n)
			deallocate(r2%yx)
			allocate(r2%yx(1:s))
			r2%yx(1:r2%n) = tmp(1:r2%n)
			r2%yx(r2%n+1:s) = 0.0_LCSRP
			!yy
			tmp(1:r2%n) = r2%yy(1:r2%n)
			deallocate(r2%yy)
			allocate(r2%yy(1:s))
			r2%yy(1:r2%n) = tmp(1:r2%n)
			r2%yy(r2%n+1:s) = 0.0_LCSRP
			!yz
			tmp(1:r2%n) = r2%yz(1:r2%n)
			deallocate(r2%yz)
			allocate(r2%yz(1:s))
			r2%yz(1:r2%n) = tmp(1:r2%n)
			r2%yz(r2%n+1:s) = 0.0_LCSRP
			!zx
			tmp(1:r2%n) = r2%zx(1:r2%n)
			deallocate(r2%zx)
			allocate(r2%zx(1:s))
			r2%zx(1:r2%n) = tmp(1:r2%n)
			r2%zx(r2%n+1:s) = 0.0_LCSRP
			!zy
			tmp(1:r2%n) = r2%zy(1:r2%n)
			deallocate(r2%zy)
			allocate(r2%zy(1:s))
			r2%zy(1:r2%n) = tmp(1:r2%n)
			r2%zy(r2%n+1:s) = 0.0_LCSRP
			!zz
			tmp(1:r2%n) = r2%zz(1:r2%n)
			deallocate(r2%zz)
			allocate(r2%zz(1:s))
			r2%zz(1:r2%n) = tmp(1:r2%n)
			r2%zz(r2%n+1:s) = 0.0_LCSRP

			r2%n = n
			deallocate(tmp)
		endif
		!reduce the size and set new n
		if(size(r2%xx)-MEM_INC > n) then
			s = n+MEM_INC/2
			allocate(tmp(1:n))

			!xx
			tmp(1:n) = r2%xx(1:n)
			deallocate(r2%xx)
			allocate(r2%xx(1:s))
			r2%xx(1:n) = tmp(1:n)
			r2%xx(n+1:s) = 0.0_LCSRP
			!xy
			tmp(1:n) = r2%xy(1:n)
			deallocate(r2%xy)
			allocate(r2%xy(1:s))
			r2%xy(1:n) = tmp(1:n)
			r2%xy(n+1:s) = 0.0_LCSRP
			!xz
			tmp(1:n) = r2%xz(1:n)
			deallocate(r2%xz)
			allocate(r2%xz(1:s))
			r2%xz(1:n) = tmp(1:n)
			r2%xz(n+1:s) = 0.0_LCSRP
			!yx
			tmp(1:n) = r2%yx(1:n)
			deallocate(r2%yx)
			allocate(r2%yx(1:s))
			r2%yx(1:n) = tmp(1:n)
			r2%yx(n+1:s) = 0.0_LCSRP
			!yy
			tmp(1:n) = r2%yy(1:n)
			deallocate(r2%yy)
			allocate(r2%yy(1:s))
			r2%yy(1:n) = tmp(1:n)
			r2%yy(n+1:s) = 0.0_LCSRP
			!yz
			tmp(1:n) = r2%yz(1:n)
			deallocate(r2%yz)
			allocate(r2%yz(1:s))
			r2%yz(1:n) = tmp(1:n)
			r2%yz(n+1:s) = 0.0_LCSRP
			!zx
			tmp(1:n) = r2%zx(1:n)
			deallocate(r2%zx)
			allocate(r2%zx(1:s))
			r2%zx(1:n) = tmp(1:n)
			r2%zx(n+1:s) = 0.0_LCSRP
			!zy
			tmp(1:n) = r2%zy(1:n)
			deallocate(r2%zy)
			allocate(r2%zy(1:s))
			r2%zy(1:n) = tmp(1:n)
			r2%zy(n+1:s) = 0.0_LCSRP
			!zz
			tmp(1:n) = r2%zz(1:n)
			deallocate(r2%zz)
			allocate(r2%zz(1:s))
			r2%zz(1:n) = tmp(1:n)
			r2%zz(n+1:s) = 0.0_LCSRP

			r2%n = n
			deallocate(tmp)
		endif
		r2%label = trim(label)
	end subroutine resize_ur2

	subroutine init_ui0(i0,n,label)
		implicit none
		!-----
		type(ui0_t):: i0
		integer:: n
		character(len=*):: label
		!-----
		integer:: s
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in init_ui0... ', trim(label)
		call destroy_ui0(i0)
		i0%n = n
		s = n+MEM_INC
		allocate(i0%i(1:s))
		i0%i = 0
		i0%label = label
	end subroutine init_ui0
	subroutine destroy_ui0(i0)
		implicit none
		!-----
		type(ui0_t):: i0
		!-----
		i0%n = 0
		if(allocated (i0%i)) deallocate(i0%i)
		i0%label = 'unused_ui0'
	end subroutine destroy_ui0
	subroutine resize_ui0(i0,n)
		implicit none
		!-----
		type(ui0_t):: i0
		integer:: n
		!-----
		integer(LCSIP),allocatable:: tmp(:)
		integer:: s
		character(len=LCS_NAMELEN):: label
		!-----
		!if(LCS_VERBOSE)&
		!	write(*,*) 'in resize_ui0...',trim(i0%label)
		
		label= trim(i0%label)

		!check that this data initialized
		if(i0%n==0) then
			call init_ui0(i0,n,i0%label)
		endif
		!deallocate if n <=0
		if (n<=0) then
			!call destroy_ui0(i0)
			!return
		endif
		!current size ok
		if( n <= size(i0%i) .AND. size(i0%i) - MEM_INC <= n)then
			!if(LCS_VERBOSE) write(*,*) 'size ok'
			i0%n = n
			!return
		endif
		!increase size and set new n
		if (n > size(i0%i)) then
			allocate(tmp(1:i0%n))
			tmp(1:i0%n) = i0%i(1:i0%n)
			deallocate(i0%i)
			s = n+MEM_INC
			allocate(i0%i(1:s))
			i0%i(1:i0%n) = tmp(1:i0%n)
			i0%i(i0%n+1:s) = 0
			i0%n = n
			deallocate(tmp)
		endif
		!reduce the size and set new n
		if(size(i0%i)-MEM_INC > n) then
			allocate(tmp(1:n))
			tmp(1:n) = i0%i(1:n)
			deallocate(i0%i)
			s = n+MEM_INC/2
			allocate(i0%i(1:s))
			i0%i(1:n) = tmp(1:n)
			i0%i(n+1:s) = 0
			i0%n = n
			deallocate(tmp)
		endif
		i0%label = trim(label)
	end subroutine resize_ui0

	subroutine init_ui1(i1,n,label)
		implicit none
		!-----
		type(ui1_t):: i1
		integer:: n
		character(len=*):: label
		!-----
		integer:: s
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in init_ui1... ', trim(label)
		call destroy_ui1(i1)
		i1%n = n
		s = n+MEM_INC
		allocate(i1%x(1:s))
		allocate(i1%y(1:s))
		allocate(i1%z(1:s))
		i1%x = 0
		i1%y = 0
		i1%z = 0
		i1%label = label
	end subroutine init_ui1
	subroutine destroy_ui1(i1)
		implicit none
		!-----
		type(ui1_t):: i1
		!-----
		i1%n = 0
		if(allocated (i1%x)) deallocate(i1%x)
		if(allocated (i1%y)) deallocate(i1%y)
		if(allocated (i1%z)) deallocate(i1%z)
		i1%label = 'unused_ui1'
	end subroutine destroy_ui1
	subroutine resize_ui1(i1,n)
		implicit none
		!-----
		type(ui1_t):: i1
		integer:: n
		!-----
		integer(LCSIP),allocatable:: tmp(:)
		integer:: s
		character(len=LCS_NAMELEN):: label
		!-----
		!if(LCS_VERBOSE)&
		!	write(*,*) 'in resize_ui1...',trim(i1%label)
		label= trim(i1%label)

		!check that this data initialized
		if(i1%n==0) then
			call init_ui1(i1,n,i1%label)
		endif
		!deallocate if n <=0
		if (n<=0) then
			!call destroy_ui1(i1)
			!return
		endif
		!current size ok
		if( n <= size(i1%x) .AND. size(i1%x) - MEM_INC <= n)then
		!	if(LCS_VERBOSE) write(*,*) 'size ok'
			i1%n = n
			!return
		endif
		!increase size and set new n
		if (n > size(i1%x)) then
			s = n+MEM_INC
			allocate(tmp(1:i1%n))

			!x
			tmp(1:i1%n) = i1%x(1:i1%n)
			deallocate(i1%x)
			allocate(i1%x(1:s))
			i1%x(1:i1%n) = tmp(1:i1%n)
			i1%x(i1%n+1:s) = 0
			!y
			tmp(1:i1%n) = i1%y(1:i1%n)
			deallocate(i1%y)
			allocate(i1%y(1:s))
			i1%y(1:i1%n) = tmp(1:i1%n)
			i1%y(i1%n+1:s) = 0
			!z
			tmp(1:i1%n) = i1%z(1:i1%n)
			deallocate(i1%z)
			allocate(i1%z(1:s))
			i1%z(1:i1%n) = tmp(1:i1%n)
			i1%z(i1%n+1:s) = 0

			i1%n = n
			deallocate(tmp)
		endif
		!reduce the size and set new n
		if(size(i1%x)-MEM_INC > n) then
			s = n+MEM_INC/2
			allocate(tmp(1:n))

			!x
			tmp(1:n) = i1%x(1:n)
			deallocate(i1%x)
			allocate(i1%x(1:s))
			i1%x(1:n) = tmp(1:n)
			i1%x(n+1:s) = 0
			!y
			tmp(1:n) = i1%y(1:n)
			deallocate(i1%y)
			allocate(i1%y(1:s))
			i1%y(1:n) = tmp(1:n)
			i1%y(n+1:s) = 0
			!z
			tmp(1:n) = i1%z(1:n)
			deallocate(i1%z)
			allocate(i1%z(1:s))
			i1%z(1:n) = tmp(1:n)
			i1%z(n+1:s) = 0

			i1%n = n
			deallocate(tmp)
		endif
		i1%label = trim(label)
	end subroutine resize_ui1

end module unstructured_m
