!
! Routines associated with structured CFD data
!
module scfd_m
	use data_m
	use structured_m
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
		!Save the boundary conditions for each proc.
		!
		scfd%bc_list(1:6) = LCS_PERIODIC
		if(scfd%offset_i==0) scfd%bc_list(1) = bc_list(1)
		if(scfd%offset_j==0) scfd%bc_list(2) = bc_list(2)
		if(scfd%offset_k==0) scfd%bc_list(3) = bc_list(3)
		if(scfd%offset_i+scfd%ni==scfd%gni) scfd%bc_list(4) = bc_list(4)
		if(scfd%offset_j+scfd%nj==scfd%gnj) scfd%bc_list(5) = bc_list(5)
		if(scfd%offset_k+scfd%nk==scfd%gnk) scfd%bc_list(6) = bc_list(6)

		!
		!Check the sign of lperiodic.
		!the convention is that x_i(n) =  x_i(1) + lperiodic(i)
		!
		if(scfd%ni > 1)then
		if(x(1,1,1) > x(2,1,1)) then
			lperiodic(1) = -lperiodic(1)
		endif
		endif
		if(scfd%nj > 1)then
		if(y(1,1,1) > y(1,2,1))then
			lperiodic(2) = -lperiodic(2)
		endif
		endif
		if(scfd%nk > 1) then
		if(z(1,1,1) > z(1,1,2))then
			lperiodic(3) = -lperiodic(3)
		endif
		endif

		!
		!Initialize the communication patterns
		!
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

		!
		!Initialize the structured grid coordinates:
		!
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

		!Set fake boundary coordinates by default
		!for the case of non-periodic external boundaries
		!These will be over-written in the case of periodicity
		!or internal boundaries below
		do k = 1-scfd%ng,scfd%nk+scfd%ng
		do j = 1-scfd%ng,scfd%nj+scfd%ng
		do i = 1,scfd%ng
			if(scfd%gni==1) then
				scfd%grid%x(1-i,j,k) 		= scfd%grid%x(1,j,k)-1.0
				scfd%grid%x(scfd%ni+i,j,k) 	= scfd%grid%x(scfd%ni,j,k)+1.0
			else
				scfd%grid%x(1-i,j,k) 		= scfd%grid%x(1,j,k) 		- (scfd%grid%x(i+1,j,k) 	- scfd%grid%x(1,j,k))
				scfd%grid%x(scfd%ni+i,j,k) 	= scfd%grid%x(scfd%ni,j,k) 	+ (scfd%grid%x(scfd%ni,j,k) - scfd%grid%x(scfd%ni-i,j,k))
			endif
		enddo
		enddo
		enddo
		do k = 1-scfd%ng,scfd%nk+scfd%ng
		do j = 1,scfd%ng
		do i = 1-scfd%ng,scfd%ni+scfd%ng
			if(scfd%gnj==1) then
				scfd%grid%y(i,1-j,k) 		= scfd%grid%y(i,1,k)-1.0
				scfd%grid%y(i,scfd%nj+i,k) 	= scfd%grid%y(i,scfd%nj,k)+1.0
			else
				scfd%grid%y(i,1-j,k) 		= scfd%grid%y(i,1,k) 		- (scfd%grid%y(i,j+1,k) 	- scfd%grid%y(i,1,k))
				scfd%grid%y(i,scfd%nj+j,k) 	= scfd%grid%y(i,scfd%nj,k) 	+ (scfd%grid%y(i,scfd%nj,k) - scfd%grid%y(i,scfd%nj-j,k))
			endif
		enddo
		enddo
		enddo
		do k = 1,scfd%ng
		do j = 1-scfd%ng,scfd%nj+scfd%ng
		do i = 1-scfd%ng,scfd%ni+scfd%ng
			if(scfd%gnk==1) then
				scfd%grid%z(i,j,1-k) 		= scfd%grid%z(i,j,1)-1.0
				scfd%grid%z(i,j,scfd%nk+k) 	= scfd%grid%z(i,j,scfd%nk) +1.0
			else
				scfd%grid%z(i,j,1-k) 		= scfd%grid%z(i,j,1) 		- (scfd%grid%z(i,j,k+1) 	- scfd%grid%z(i,j,1))
				scfd%grid%z(i,j,scfd%nk+k) 	= scfd%grid%z(i,j,scfd%nk) 	+ (scfd%grid%z(i,j,scfd%nk) - scfd%grid%z(i,j,scfd%nk-k))
			endif
		enddo
		enddo
		enddo

		!Exchange to set ghost coordinates:
		call exchange_sdata(scfd%scomm_max_r1,r1=scfd%grid)

		!
		!Initialize velocity field:
		!
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

	subroutine set_velocity_bc(bc_list,vel)
		implicit none
		!-----
		integer, intent(in):: bc_list(6)
		type(sr1_t):: vel
		!-----
		integer:: i,j,k,ni,nj,nk,ng
		!-----

		if(lcsrank==0) &
			write(*,*) 'In set_velocity_bc... '

		ni = vel%ni
		nj = vel%nj
		nk = vel%nk
		ng = vel%ng


		!
		! X0
		!
		select case(bc_list(1))
			case(LCS_WALL) !No slip...
				do i = 0,1-ng
					vel%x(i,:,:) = 0.0
					vel%y(i,:,:) = 0.0
					vel%z(i,:,:) = 0.0
				enddo
			case(LCS_INFLOW,LCS_OUTFLOW,LCS_SLIP)
				do i = 0,1-ng
					vel%x(i,:,:) = vel%x(1,:,:)
					vel%y(i,:,:) = vel%y(1,:,:)
					vel%z(i,:,:) = vel%z(1,:,:)
				enddo
			case default
				!do nothing...
		end select
		!
		! X1
		!
		select case(bc_list(4))
			case(LCS_WALL) !No slip...
				do i = ni+1,ni+ng
					vel%x(i,:,:) = 0.0
					vel%y(i,:,:) = 0.0
					vel%z(i,:,:) = 0.0
				enddo
			case(LCS_INFLOW,LCS_OUTFLOW,LCS_SLIP) !zero gradient
				do i = ni+1,ni+ng
					vel%x(i,:,:) = vel%x(ni,:,:)
					vel%y(i,:,:) = vel%y(ni,:,:)
					vel%z(i,:,:) = vel%z(ni,:,:)
				enddo
			case default
				!do nothing...
		end select
		!
		! Y0
		!
		select case(bc_list(2))
			case(LCS_WALL) !No slip...
				do j = 0,1-ng
					vel%x(:,j,:) = 0.0
					vel%y(:,j,:) = 0.0
					vel%z(:,j,:) = 0.0
				enddo
			case(LCS_INFLOW,LCS_OUTFLOW,LCS_SLIP)
				do j = 0,1-ng
					vel%x(:,j,:) = vel%x(:,1,:)
					vel%y(:,j,:) = vel%y(:,1,:)
					vel%z(:,j,:) = vel%z(:,1,:)
				enddo
			case default
				!do nothing...
		end select
		!
		! Y1
		!
		select case(bc_list(5))
			case(LCS_WALL) !No slip...
				do j = nj+1,nj+ng
					vel%x(:,j,:) = 0.0
					vel%y(:,j,:) = 0.0
					vel%z(:,j,:) = 0.0
				enddo
			case(LCS_INFLOW,LCS_OUTFLOW,LCS_SLIP) !zero gradient
				do j = nj+1,nj+ng
					vel%x(:,j,:) = vel%x(:,nj,:)
					vel%y(:,j,:) = vel%y(:,nj,:)
					vel%z(:,j,:) = vel%z(:,nj,:)
				enddo
			case default
				!do nothing...
		end select
		!
		! Z0
		!
		select case(bc_list(3))
			case(LCS_WALL) !No slip...
				do k = 0,1-ng
					vel%x(:,:,k) = 0.0
					vel%y(:,:,k) = 0.0
					vel%z(:,:,k) = 0.0
				enddo
			case(LCS_INFLOW,LCS_OUTFLOW,LCS_SLIP)
				do k = 0,1-ng
					vel%x(:,:,k) = vel%x(:,:,1)
					vel%y(:,:,k) = vel%y(:,:,1)
					vel%z(:,:,k) = vel%z(:,:,1)
				enddo
			case default
				!do nothing...
		end select
		!
		! Z1
		!
		select case(bc_list(6))
			case(LCS_WALL) !No slip...
				do k = nk+1,nk+ng
					vel%x(:,:,k) = 0.0
					vel%y(:,:,k) = 0.0
					vel%z(:,:,k) = 0.0
				enddo
			case(LCS_INFLOW,LCS_OUTFLOW,LCS_SLIP) !zero gradient
				do k = nk+1,nk+ng
					vel%x(:,:,k) = vel%x(:,:,nk)
					vel%y(:,:,k) = vel%y(:,:,nk)
					vel%z(:,:,k) = vel%z(:,:,nk)
				enddo
			case default
				!do nothing...
		end select

	end subroutine set_velocity_bc


end module scfd_m















