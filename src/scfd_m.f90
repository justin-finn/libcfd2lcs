!
! Routines associated with structured CFD data
!
module scfd_m
	use data_m
	use comms_m

	contains
	subroutine init_scfd(scfd,label,n,offset,x,y,z,bc_list,lperiodic)
		implicit none
		!-----
		type(scfd_t):: scfd
		character(len=*):: label

		integer:: n(3),offset(3)
		real(LCSRP):: x(1:n(1),1:n(2),1:n(3))
		real(LCSRP):: y(1:n(1),1:n(2),1:n(3))
		real(LCSRP):: z(1:n(1),1:n(2),1:n(3))
		real(LCSRP):: tmpx(1:n(1),1:n(2),1:n(3))
		real(LCSRP):: tmpy(1:n(1),1:n(2),1:n(3))
		real(LCSRP):: tmpz(1:n(1),1:n(2),1:n(3))
		integer(LCSIP),dimension(6):: bc_list
		real(LCSRP):: lperiodic(3)
		!-----
		integer:: i,j,k,ig
		integer:: my_nmax(3), nmax(3), offset_min(3)
		integer:: ierr
		!-----

		if(lcsrank==0) &
			write(*,*) 'in init_scfd... ', trim(label)

		!Set the label
		scfd%label = trim(label)

		!
		!Initialize the structured data:
		!Regardless of what gets passed in, the data stored by cfd2lcs follows the following rules:
		!	1.  Each processor owns arrays with bounds [ 1:scfd%ni, 1:scfd%nj, 1:scfd%nk ]
		!   2.  The minimum offset in each direction is 0
		!	3.  The global number of grid points is gni  gnj  gnk.
		!   4.  ng must not be greater than ni,nj, or nk
		!
		scfd%ng = NGHOST_CFD
		scfd%ni = n(1)
		scfd%nj = n(2)
		scfd%nk = n(3)
		if(scfd%ng > scfd%ni .OR. scfd%ng > scfd%nj .OR. scfd%ng > scfd%nk) then
			write(*,*) 'ERROR:  lcsrank[',lcsrank,'] has ng > ni,nj or nk', scfd%ni,scfd%nj,scfd%nk,scfd%ng
			CFD2LCS_ERROR = 1
			return
		endif
		call MPI_ALLREDUCE(offset,offset_min,3,MPI_INTEGER,MPI_MIN,lcscomm,ierr)
		scfd%offset_i = offset(1) - offset_min(1)
		scfd%offset_j = offset(2) - offset_min(2)
		scfd%offset_k = offset(3) - offset_min(3)
		my_nmax(1) = scfd%offset_i + scfd%ni
		my_nmax(2) = scfd%offset_j + scfd%nj
		my_nmax(3) = scfd%offset_k + scfd%nk
		call MPI_ALLREDUCE(my_nmax,nmax,3,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
		scfd%gni = nmax(1)
		scfd%gnj = nmax(2)
		scfd%gnk = nmax(3)

		!
		!Check the sign of lperiodic.
		!the convention is that x_i(n) =  x_i(1) + lperiodic(i)
		!
		if(scfd%ni > 1 .AND. x(1,1,1) > x(2,1,1)) then
			lperiodic(1) = -lperiodic(1)
		endif
		if(scfd%nj > 1 .AND. y(1,1,1) > y(1,2,1))then
			lperiodic(2) = -lperiodic(2)
		endif
		if(scfd%nk > 1 .AND. z(1,1,1) > z(1,1,2))then
			lperiodic(3) = -lperiodic(3)
		endif

		!Initialize the communication patterns
		call init_scomm(scfd%scomm_face_r0,scfd%ni,scfd%nj,scfd%nk,scfd%ng,&
			scfd%offset_i,scfd%offset_j,scfd%offset_k,bc_list,lperiodic,FACE_CONNECT,R0_COMM,'R0 face-nbr comms' )
		call init_scomm(scfd%scomm_max_r0,scfd%ni,scfd%nj,scfd%nk,scfd%ng,&
			scfd%offset_i,scfd%offset_j,scfd%offset_k,bc_list,lperiodic,MAX_CONNECT,R0_COMM,'R0 max-nbr comms' )
		call init_scomm(scfd%scomm_face_r1,scfd%ni,scfd%nj,scfd%nk,scfd%ng,&
			scfd%offset_i,scfd%offset_j,scfd%offset_k,bc_list,lperiodic,FACE_CONNECT,R1_COMM,'R1 face-nbr comms' )
		call init_scomm(scfd%scomm_max_r1,scfd%ni,scfd%nj,scfd%nk,scfd%ng,&
			scfd%offset_i,scfd%offset_j,scfd%offset_k,bc_list,lperiodic,MAX_CONNECT,R1_COMM,'R1 max-nbr comms' )
		call init_scomm(scfd%scomm_face_r2,scfd%ni,scfd%nj,scfd%nk,scfd%ng,&
			scfd%offset_i,scfd%offset_j,scfd%offset_k,bc_list,lperiodic,FACE_CONNECT,R2_COMM,'R2 face-nbr comms' )
		call init_scomm(scfd%scomm_max_r2,scfd%ni,scfd%nj,scfd%nk,scfd%ng,&
			scfd%offset_i,scfd%offset_j,scfd%offset_k,bc_list,lperiodic,MAX_CONNECT,R2_COMM,'R2 max-nbr comms' )

		!Initialize the structured grid coordinates:
		call init_sr1(scfd%grid,scfd%ni,scfd%nj,scfd%nk,scfd%ng,'GRID',translate=.true.)

		!Set interior points:
		do k=1,scfd%nk
		do j=1,scfd%nj
		do i=1,scfd%ni
			scfd%grid%x(i,j,k) = x(i,j,k)
			scfd%grid%y(i,j,k) = y(i,j,k)
			scfd%grid%z(i,j,k) = z(i,j,k)
		enddo
		enddo
		enddo

		!TODO:  Set fake boundary coordinates
		!for the case of non-periodic external boundaries

		!Exchange to set ghost coordinates:
		call exchange_sdata(scfd%scomm_max_r1,r1=scfd%grid)

		!Initialize velocity field:
		call init_sr1(scfd%u,scfd%ni,scfd%nj,scfd%nk,scfd%ng,'U',translate=.false.)

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

		call destroy_scomm(scfd%scomm_face_r0)
		call destroy_scomm(scfd%scomm_max_r0)
		call destroy_scomm(scfd%scomm_face_r1)
		call destroy_scomm(scfd%scomm_max_r1)
		call destroy_scomm(scfd%scomm_face_r2)
		call destroy_scomm(scfd%scomm_max_r2)


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


	subroutine init_sr1(r1,ni,nj,nk,ng,label,translate)
		implicit none
		!-----
		type(sr1_t):: r1
		integer:: ni,nj,nk,ng
		character(len=*):: label
		logical:: translate
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


	subroutine grad_sr1(ni,nj,nk,ng,grid,sr1,grad)
		implicit none
		!----
		integer,intent(in):: ni,nj,nk,ng
		type(sr1_t),intent(in):: grid
		type(sr1_t),intent(in):: sr1
		type(sr2_t),intent(inout):: grad
		!----
		integer:: i,j,k
		!----

		if (lcsrank==0) &
			write(*,*) 'in grad_sr1... ',trim(sr1%label),' => ',trim(grad%label)

		!2nd order central scheme:
		do k = 1,nk
		do j = 1,nj
		do i = 1,ni
			grad%xx(i,j,k) = (sr1%x(i+1,j,k)-sr1%x(i-1,j,k)) / (grid%x(i+1,j,k) - grid%x(i-1,j,k))
			grad%xy(i,j,k) = (sr1%x(i,j+1,k)-sr1%x(i,j-1,k)) / (grid%y(i,j+1,k) - grid%y(i,j,k-1))
			grad%xz(i,j,k) = (sr1%x(i,j,k+1)-sr1%x(i,j,k-1)) / (grid%z(i,j,k+1) - grid%z(i,j,k-1))

			grad%yx(i,j,k) = (sr1%y(i+1,j,k)-sr1%y(i-1,j,k)) / (grid%x(i+1,j,k) - grid%x(i-1,j,k))
			grad%yy(i,j,k) = (sr1%y(i,j+1,k)-sr1%y(i,j-1,k)) / (grid%y(i,j+1,k) - grid%y(i,j,k-1))
			grad%yz(i,j,k) = (sr1%y(i,j,k+1)-sr1%y(i,j,k-1)) / (grid%z(i,j,k+1) - grid%z(i,j,k-1))

			grad%zx(i,j,k) = (sr1%z(i+1,j,k)-sr1%z(i-1,j,k)) / (grid%x(i+1,j,k) - grid%x(i-1,j,k))
			grad%zy(i,j,k) = (sr1%z(i,j+1,k)-sr1%z(i,j-1,k)) / (grid%y(i,j+1,k) - grid%y(i,j,k-1))
			grad%zz(i,j,k) = (sr1%z(i,j,k+1)-sr1%z(i,j,k-1)) / (grid%z(i,j,k+1) - grid%z(i,j,k-1))
		enddo
		enddo
		enddo

	end subroutine grad_sr1

end module scfd_m















