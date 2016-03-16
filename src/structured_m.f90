module structured_m
	use data_m
	implicit none

	contains

	subroutine init_sr0(r0,ni,nj,nk,ng,label)
		implicit none
		!-----
		type(sr0_t):: r0
		integer:: ni,nj,nk,ng
		character(len=*):: label
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in init_sr0... ', trim(label)
		call destroy_sr0(r0)
		r0%ni = ni
		r0%nj = nj
		r0%nk = nk
		r0%ng = ng
		allocate(r0%r(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		r0%r = 0.0_LCSRP
		r0%label = label
	end subroutine init_sr0
	subroutine destroy_sr0(r0)
		implicit none
		!-----
		type(sr0_t):: r0
		!-----
		r0%ni = 0
		r0%nj = 0
		r0%nk = 0
		r0%ng = 0
		if(allocated (r0%r)) deallocate(r0%r)
		r0%label = 'unused_sr0'
	end subroutine destroy_sr0

	subroutine init_si0(i0,ni,nj,nk,ng,label)
		implicit none
		!-----
		type(si0_t):: i0
		integer:: ni,nj,nk,ng
		character(len=*):: label
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in init_si0... ', trim(label)
		call destroy_si0(i0)
		i0%ni = ni
		i0%nj = nj
		i0%nk = nk
		i0%ng = ng
		allocate(i0%i(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		i0%i = 0
		i0%label = label
	end subroutine init_si0
	subroutine destroy_si0(i0)
		implicit none
		!-----
		type(si0_t):: i0
		!-----
		i0%ni = 0
		i0%nj = 0
		i0%nk = 0
		i0%ng = 0
		if(allocated (i0%i)) deallocate(i0%i)
		i0%label = 'unused_si0'
	end subroutine destroy_si0

	subroutine init_sr1(r1,ni,nj,nk,ng,label,translate)
		implicit none
		!-----
		type(sr1_t):: r1
		integer:: ni,nj,nk,ng
		character(len=*):: label
		logical:: translate
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in init_sr1... ', trim(label)
		call destroy_sr1(r1)
		r1%ni = ni
		r1%nj = nj
		r1%nk = nk
		r1%ng = ng
		allocate(r1%x(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r1%y(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r1%z(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		r1%x = 0.0_LCSRP
		r1%y = 0.0_LCSRP
		r1%z = 0.0_LCSRP
		r1%label = label
		r1%periodic_translate = translate
	end subroutine init_sr1
	subroutine destroy_sr1(r1)
		implicit none
		!-----
		type(sr1_t):: r1
		!-----
		r1%ni = 0
		r1%nj = 0
		r1%nk = 0
		r1%ng = 0
		if(allocated (r1%x)) deallocate(r1%x)
		if(allocated (r1%y)) deallocate(r1%y)
		if(allocated (r1%z)) deallocate(r1%z)
		r1%periodic_translate = .false.
		r1%label = 'unused_sr1'
	end subroutine destroy_sr1

	subroutine init_sr2(r2,ni,nj,nk,ng,label)
		implicit none
		!-----
		type(sr2_t):: r2
		integer:: ni,nj,nk,ng
		character(len=*):: label
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in init_sr2... ', trim(label)
		call destroy_sr2(r2)
		r2%ni = ni
		r2%nj = nj
		r2%nk = nk
		r2%ng = ng
		allocate(r2%xx(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r2%xy(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r2%xz(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r2%yx(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r2%yy(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r2%yz(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r2%zx(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r2%zy(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r2%zz(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
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
	end subroutine init_sr2
	subroutine destroy_sr2(r2)
		implicit none
		!-----
		type(sr2_t):: r2
		!-----
		r2%ni = 0
		r2%nj = 0
		r2%nk = 0
		r2%ng = 0
		if(allocated (r2%xx)) deallocate(r2%xx)
		if(allocated (r2%xy)) deallocate(r2%xy)
		if(allocated (r2%xz)) deallocate(r2%xz)
		if(allocated (r2%yx)) deallocate(r2%yx)
		if(allocated (r2%yy)) deallocate(r2%yy)
		if(allocated (r2%yz)) deallocate(r2%yz)
		if(allocated (r2%zx)) deallocate(r2%zx)
		if(allocated (r2%zy)) deallocate(r2%zy)
		if(allocated (r2%zz)) deallocate(r2%zz)
		r2%label = 'unused_sr2'
	end subroutine destroy_sr2

end module structured_m










