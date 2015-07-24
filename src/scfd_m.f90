!
! Routines associated with structured CFD data
!
module scfd_m
	use data_m
	use comms_m

	contains
	subroutine init_scfd(scfd,label,ni,nj,nk,offset_i,offset_j,offset_k,x,y,z,bc_list)
		implicit none
		!-----
		type(scfd_t):: scfd
		character(len=*):: label
		integer:: ni,nj,nk,offset_i,offset_j,offset_k
		real(LCSRP):: x(1:ni,1:nj,1:nk)
		real(LCSRP):: y(1:ni,1:nj,1:nk)
		real(LCSRP):: z(1:ni,1:nj,1:nk)
		integer(LCSIP),dimension(6):: bc_list
		!-----
		integer:: i,j,k
		integer:: my_nmax(3), nmax(3)
		integer:: ierr
		!-----

		if(lcsrank==0) &
			write(*,*) 'in init_scfd... ', trim(label)

		!Set the label
		scfd%label = trim(label)

		!Initialize the structured data:
		scfd%ng = NGHOST_CFD
		scfd%ni = ni
		scfd%nj = nj
		scfd%nk = nk
		scfd%offset_i = offset_i
		scfd%offset_j = offset_j
		scfd%offset_k = offset_k
		my_nmax(1) = offset_i + ni
		my_nmax(2) = offset_j + nj
		my_nmax(3) = offset_k + nk
		call MPI_ALLREDUCE(my_nmax,nmax,3,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
		scfd%gni = nmax(1)
		scfd%gnj = nmax(2)
		scfd%gnk = nmax(3)

		!Initialize the communication pattern
		call init_scomm(scfd%scomm_face,ni,nj,nk,scfd%ng,offset_i,offset_j,offset_k,bc_list,FACE_CONNECT,'Struct. face-nbr comms' )
		call init_scomm(scfd%scomm_max,ni,nj,nk,scfd%ng,offset_i,offset_j,offset_k,bc_list,MAX_CONNECT,'Struct. max-nbr comms' )

		!Initialize the structured grid coordinates:
		call init_sr1(scfd%grid,scfd%ni,scfd%nj,scfd%nk,scfd%ng,'GRID')
		do k=1,nk
		do j=1,nj
		do i=1,ni
			scfd%grid%x(i,j,k) = x(i,j,k)
			scfd%grid%y(i,j,k) = y(i,j,k)
			scfd%grid%z(i,j,k) = z(i,j,k)
		enddo
		enddo
		enddo

		!Initialize velocity field:
		call init_sr1(scfd%u,scfd%ni,scfd%nj,scfd%nk,scfd%ng,'U')

	end subroutine init_scfd
	subroutine destroy_scfd(scfd)
		implicit none
		!-----
		type(scfd_t):: scfd
		!-----
		if(lcsrank==0)&
			write(*,*) 'in destroy_scfd...'

		scfd%ni = 0
		scfd%nj = 0
		scfd%nk = 0
		scfd%ng = 0
		scfd%gni = 0
		scfd%gnj = 0
		scfd%gnk = 0
		scfd%offset_i = 0
		scfd%offset_j = 0
		scfd%offset_k = 0

		call destroy_scomm(scfd%scomm_face)
		call destroy_scomm(scfd%scomm_max)

		call destroy_sr1(scfd%grid)
		call destroy_sr1(scfd%u)
		scfd%label = 'Unused scfd'

	end subroutine destroy_scfd

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
		r1%x = 0.0_LCSRP
		r1%y = 0.0_LCSRP
		r1%z = 0.0_LCSRP
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

end module scfd_m















