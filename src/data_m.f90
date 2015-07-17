module data_m
	use precision_m
	use mpi_m
	implicit none
	!-------
	!Data structure definition
	!-------
	integer,parameter:: MAX_NAMELEN = 32
	integer,parameter:: NGHOST_CFD = 1

	!Real, structured Cartesian rank 0 scalar
	type sr0_t
		integer:: ni,nj,nk,ng
		character(len=MAX_NAMELEN):: label
		real(WP), allocatable:: r(:,:,:)
	end type sr0_t

	!Real, structured Cartesian rank 1 vector
	type sr1_t
		integer:: ni,nj,nk,ng
		character(len=MAX_NAMELEN):: label
		real(WP), allocatable:: x(:,:,:)
		real(WP), allocatable:: y(:,:,:)
		real(WP), allocatable:: z(:,:,:)
	end type sr1_t

	!Real, structured Cartesian rank 2 tensor
	type sr2_t
		integer:: ni,nj,nk,ng
		character(len=MAX_NAMELEN):: label
		real(WP), allocatable:: xx(:,:,:)
		real(WP), allocatable:: xy(:,:,:)
		real(WP), allocatable:: xz(:,:,:)
		real(WP), allocatable:: yx(:,:,:)
		real(WP), allocatable:: yy(:,:,:)
		real(WP), allocatable:: yz(:,:,:)
		real(WP), allocatable:: zx(:,:,:)
		real(WP), allocatable:: zy(:,:,:)
		real(WP), allocatable:: zz(:,:,:)
	end type sr2_t

	!Structured, Cartesian grid:
	type cart_t
		integer:: ni, nj, nk, ng
		integer:: gni, gnj, gnk
		integer:: offset_i, offset_j, offset_k
		real(WP),allocatable:: x(:), y(:), z(:)
		!Store the communications buffers here?
	end type cart_t

	!Structured cfd data:
	type scfd_t
		character(len=MAX_NAMELEN):: label
		type(cart_t):: cart
		type(sr1_t):: u

	end type scfd_t

	!The CFD side data:
	type(scfd_t):: cfd

	contains

	subroutine init_cfd(cfd,label,ni,nj,nk,offset_i,offset_j,offset_k,x,y,z)
		implicit none
		!-----
		type(scfd_t):: cfd
		character(len=*):: label
		integer:: ni,nj,nk,offset_i,offset_j,offset_k
		real(WP):: x(1:ni), y(1:nj), z(1:nk)
		!-----
		integer:: ng
		!-----

		if(lcsrank==0) &
			write(*,*) 'in init_cfd... ', trim(label)

		!Set the label
		cfd%label = trim(label)

		!Initialize the grid:
		ng = NGHOST_CFD
		call init_cart(cfd%cart,ni,nj,nk,ng,offset_i,offset_j,offset_k,x,y,z)

		!Initialize velocity field:
		call init_sr1(cfd%u,ni,nj,nk,ng,'U')

	end subroutine init_cfd
	subroutine destroy_cfd(cfd)
		implicit none
		!-----
		type(scfd_t):: cfd
		!-----
		if(lcsrank==0)&
			write(*,*) 'in destroy_cfd...'

		call destroy_cart(cfd%cart)
		call destroy_sr1(cfd%u)
		cfd%label = 'Unused scfd'

	end subroutine destroy_cfd



	subroutine init_cart(cart,ni,nj,nk,ng,offset_i,offset_j,offset_k,x,y,z)
		implicit none
		!----
		type(cart_t):: cart
		integer:: ni,nj,nk,ng
		integer:: offset_i,offset_j,offset_k
		real(WP):: x(1:ni), y(1:nj), z(1:nk)
		!----
		integer:: ierr
		integer:: my_nmax(3),nmax(3)
		!----
		!Initialize the default Cartesian data structure
		!----

		if(lcsrank==0) &
			write(*,*) 'in init_cart...'

		cart%ni = ni
		cart%nj = nj
		cart%nk = nk
		cart%offset_i = offset_i
		cart%offset_j = offset_j
		cart%offset_k = offset_k
		my_nmax(1) = offset_i + ni
		my_nmax(2) = offset_j + nj
		my_nmax(3) = offset_k + nk
		call MPI_ALLREDUCE(my_nmax,nmax,3,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
		cart%gni = nmax(1)
		cart%gnj = nmax(2)
		cart%gnk = nmax(3)
		allocate(cart%x(1-ng:ni+ng))
		allocate(cart%y(1-ng:nj+ng))
		allocate(cart%z(1-ng:nk+ng))
		cart%x(1:ni) = x(1:ni)
		cart%y(1:nj) = y(1:nj)
		cart%z(1:nk) = z(1:nk)

		!initialize communications and exchange ghost coordinates:
		!TODO

	end subroutine init_cart
	subroutine destroy_cart(cart)
		implicit none
		!-----
		type(cart_t):: cart
		!-----
		cart%ni = 0
		cart%nj = 0
		cart%nk = 0
		cart%ng = 0
		cart%gni = 0
		cart%gnj = 0
		cart%gnk = 0
		cart%offset_i = 0
		cart%offset_j = 0
		cart%offset_k = 0
		if(allocated(cart%x))deallocate(cart%x)
		if(allocated(cart%y))deallocate(cart%y)
		if(allocated(cart%z))deallocate(cart%z)
	end subroutine destroy_cart


	subroutine init_sr0(r0,ni,nj,nk,ng,label)
		implicit none
		!-----
		type(sr0_t):: r0
		integer:: ni,nj,nk,ng
		character(len=*):: label
		!-----
		if(lcsrank==0)&
			write(*,*) 'in init_sr0... ', trim(label)

		call destroy_sr0(r0)
		r0%ni = ni
		r0%nj = nj
		r0%nk = nk
		r0%ng = ng
		allocate(r0%r(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		r0%r = 0.0_WP
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


	subroutine init_sr1(r1,ni,nj,nk,ng,label)
		implicit none
		!-----
		type(sr1_t):: r1
		integer:: ni,nj,nk,ng
		character(len=*):: label
		!-----
		if(lcsrank==0)&
			write(*,*) 'in init_sr1... ', trim(label)

		call destroy_sr1(r1)
		r1%ni = ni
		r1%nj = nj
		r1%nk = nk
		r1%ng = ng
		allocate(r1%x(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r1%y(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		allocate(r1%z(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		r1%x = 0.0_WP
		r1%y = 0.0_WP
		r1%z = 0.0_WP
		r1%label = label
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
		r1%label = 'unused_sr1'
	end subroutine destroy_sr1

	subroutine init_sr2(r2,ni,nj,nk,ng,label)
		implicit none
		!-----
		type(sr2_t):: r2
		integer:: ni,nj,nk,ng
		character(len=*):: label
		!-----
		if(lcsrank==0)&
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
		r2%xx = 0.0_WP
		r2%xy = 0.0_WP
		r2%xz = 0.0_WP
		r2%yx = 0.0_WP
		r2%yy = 0.0_WP
		r2%yz = 0.0_WP
		r2%zx = 0.0_WP
		r2%zy = 0.0_WP
		r2%zz = 0.0_WP
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

end module data_m
