!
! Routines associated with structured grids
!
module sgrid_m
	use data_m
	use structured_m
	use comms_m
	contains
	subroutine init_sgrid(sgrid,label,n,offset,x,y,z,bc_list,lperiodic)
		implicit none
		!-----
		type(sgrid_t),pointer:: sgrid
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
		integer:: i,j,k,isg,ilcs
		integer:: my_nmax(3), nmax(3), offset_min(3)
		integer:: ierr
		type(sgrid_t),allocatable:: sgrid_c_tmp(:)
		integer:: scfdptr
		integer,allocatable:: lcsptr(:)

		real(LCSRP),allocatable:: xb(:),yb(:),zb(:)
		!-----

		if(lcsrank==0) &
			write(*,*) 'in init_sgrid... ', trim(label)

		!----
		!Add a new item to the sgrid collection (sgrid_c array)
		!Check the association of scfd%sgrid, lcs%sgrid,
		!and any other structures that have an sgrid, so you can
		!preserve the pointer assignment after expansion
		!----
		if(NSGRID == 0 ) then
			NSGRID = NSGRID + 1
			allocate(sgrid_c(NSGRID))
		else
			allocate(sgrid_c_tmp(NSGRID))

			!scfd  ptr
			scfdptr = -1
			do isg = 1, NSGRID
				if (associated(scfd%sgrid,sgrid_c(isg))) then
					scfdptr = isg
				endif
			enddo

			!lcs ptrs
			if(NLCS>0) then
				allocate(lcsptr(1:NLCS))
				lcsptr = -1
				do ilcs = 1,NLCS
					if (associated(lcs_c(ilcs)%sgrid,scfd%sgrid)) then
						lcsptr(ilcs) = 0
					endif
					do isg = 1,NSGRID
						if (associated(lcs_c(ilcs)%sgrid,sgrid_c(isg))) then
							lcsptr(ilcs) = isg
						endif
					enddo
				enddo
			endif

			!expand array of structures
			sgrid_c_tmp = sgrid_c
			deallocate(sgrid_c)
			allocate(sgrid_c(NSGRID+1))
			sgrid_c(1:NSGRID) = sgrid_c_tmp(1:NSGRID)

			!fix the old lcs ptrs
			if(scfdptr>0) then
				scfd%sgrid => sgrid_c(scfdptr)
			endif
			do ilcs = 1,NLCS
				if(lcsptr(ilcs) == 0) then
					lcs_c(ilcs)%sgrid => scfd%sgrid
				endif
				if(lcsptr(ilcs) > 0) then
					lcs_c(ilcs)%sgrid => sgrid_c(lcsptr(ilcs))
				endif
			enddo

			NSGRID = NSGRID + 1
		endif

		!the new scfd ptr
		sgrid => sgrid_c(NSGRID)

		!Set the label
		sgrid%label = trim(label)

		!
		!Initialize the structured data:
		!Regardless of what gets passed in, the data stored by cfd2lcs follows the following rules:
		!	1.  Each processor owns arrays with bounds [ 1:sgrid%ni, 1:sgrid%nj, 1:sgrid%nk ]
		!   2.  The minimum offset in each direction is 0
		!	3.  The global number of grid points is gni  gnj  gnk.
		!   4.  ng must not be greater than ni,nj, or nk
		!
		sgrid%ng = NGHOST_CFD
		sgrid%ni = n(1)
		sgrid%nj = n(2)
		sgrid%nk = n(3)
		if(sgrid%ng > sgrid%ni .OR. sgrid%ng > sgrid%nj .OR. sgrid%ng > sgrid%nk) then
			write(*,*) 'ERROR:  lcsrank[',lcsrank,'] has ng > ni,nj or nk', sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng
			CFD2LCS_ERROR = 1
			return
		endif
		call MPI_ALLREDUCE(offset,offset_min,3,MPI_INTEGER,MPI_MIN,lcscomm,ierr)
		sgrid%offset_i = offset(1) - offset_min(1)
		sgrid%offset_j = offset(2) - offset_min(2)
		sgrid%offset_k = offset(3) - offset_min(3)
		my_nmax(1) = sgrid%offset_i + sgrid%ni
		my_nmax(2) = sgrid%offset_j + sgrid%nj
		my_nmax(3) = sgrid%offset_k + sgrid%nk
		call MPI_ALLREDUCE(my_nmax,nmax,3,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
		sgrid%gni = nmax(1)
		sgrid%gnj = nmax(2)
		sgrid%gnk = nmax(3)

		!
		!Save the boundary conditions for each proc.
		!
		sgrid%global_bc_list = bc_list
		sgrid%bc_list(1:6) = LCS_PERIODIC
		if(sgrid%offset_i==0) sgrid%bc_list(1) = bc_list(1)
		if(sgrid%offset_j==0) sgrid%bc_list(2) = bc_list(2)
		if(sgrid%offset_k==0) sgrid%bc_list(3) = bc_list(3)
		if(sgrid%offset_i+sgrid%ni==sgrid%gni) sgrid%bc_list(4) = bc_list(4)
		if(sgrid%offset_j+sgrid%nj==sgrid%gnj) sgrid%bc_list(5) = bc_list(5)
		if(sgrid%offset_k+sgrid%nk==sgrid%gnk) sgrid%bc_list(6) = bc_list(6)

		!
		!Check the sign of lperiodic.
		!the convention is that x_i(n) =  x_i(1) + lperiodic(i)
		!
		sgrid%lperiodic(1:3) = lperiodic(1:3)
		if(sgrid%ni > 1)then
		if(x(1,1,1) > x(2,1,1)) then
			lperiodic(1) = -lperiodic(1)
		endif
		endif
		if(sgrid%nj > 1)then
		if(y(1,1,1) > y(1,2,1))then
			lperiodic(2) = -lperiodic(2)
		endif
		endif
		if(sgrid%nk > 1) then
		if(z(1,1,1) > z(1,1,2))then
			lperiodic(3) = -lperiodic(3)
		endif
		endif

		!
		!Initialize the communication patterns
		!
		call init_scomm(sgrid%scomm_face_r0,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,bc_list,lperiodic,FACE_CONNECT,R0_COMM,'R0 face-nbr comms' )
		call init_scomm(sgrid%scomm_max_r0,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,bc_list,lperiodic,MAX_CONNECT,R0_COMM,'R0 max-nbr comms' )
		call init_scomm(sgrid%scomm_face_r1,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,bc_list,lperiodic,FACE_CONNECT,R1_COMM,'R1 face-nbr comms' )
		call init_scomm(sgrid%scomm_max_r1,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,bc_list,lperiodic,MAX_CONNECT,R1_COMM,'R1 max-nbr comms' )
		call init_scomm(sgrid%scomm_face_r2,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,bc_list,lperiodic,FACE_CONNECT,R2_COMM,'R2 face-nbr comms' )
		call init_scomm(sgrid%scomm_max_r2,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,bc_list,lperiodic,MAX_CONNECT,R2_COMM,'R2 max-nbr comms' )

		!
		!Initialize the structured grid coordinates:
		!
		call init_sr1(sgrid%grid,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,'GRID',translate=.true.)

		!Set interior points:
		do k=1,sgrid%nk
		do j=1,sgrid%nj
		do i=1,sgrid%ni
			sgrid%grid%x(i,j,k) = x(i,j,k)
			sgrid%grid%y(i,j,k) = y(i,j,k)
			sgrid%grid%z(i,j,k) = z(i,j,k)
		enddo
		enddo
		enddo

		!Set fake boundary coordinates by default
		!for the case of non-periodic external boundaries
		!These will be over-written in the case of periodicity
		!or internal boundaries below
		!JRF:  set all coordinates for handling lp-bc
		do k = 1-sgrid%ng,sgrid%nk+sgrid%ng
		do j = 1-sgrid%ng,sgrid%nj+sgrid%ng
		do i = 1,sgrid%ng
			if(sgrid%gni==1) then
				sgrid%grid%x(1-i,j,k) 		= sgrid%grid%x(1,j,k)-1.0
				sgrid%grid%y(1-i,j,k) 		= sgrid%grid%y(1,j,k)!-1.0
				sgrid%grid%z(1-i,j,k) 		= sgrid%grid%z(1,j,k)!-1.0
				sgrid%grid%x(sgrid%ni+i,j,k) 	= sgrid%grid%x(sgrid%ni,j,k)+1.0
				sgrid%grid%y(sgrid%ni+i,j,k) 	= sgrid%grid%y(sgrid%ni,j,k)!+1.0
				sgrid%grid%z(sgrid%ni+i,j,k) 	= sgrid%grid%z(sgrid%ni,j,k)!+1.0
			else
				sgrid%grid%x(1-i,j,k) 		= sgrid%grid%x(1,j,k) &
				- (sgrid%grid%x(i+1,j,k) 	- sgrid%grid%x(1,j,k))
				sgrid%grid%y(1-i,j,k) 		= sgrid%grid%y(1,j,k) &
				- (sgrid%grid%y(i+1,j,k) 	- sgrid%grid%y(1,j,k))
				sgrid%grid%z(1-i,j,k) 		= sgrid%grid%z(1,j,k) &
				- (sgrid%grid%z(i+1,j,k) 	- sgrid%grid%z(1,j,k))
				sgrid%grid%x(sgrid%ni+i,j,k) 	= sgrid%grid%x(sgrid%ni,j,k) &
				+ (sgrid%grid%x(sgrid%ni,j,k) - sgrid%grid%x(sgrid%ni-i,j,k))
				sgrid%grid%y(sgrid%ni+i,j,k) 	= sgrid%grid%y(sgrid%ni,j,k) &
				+ (sgrid%grid%y(sgrid%ni,j,k) - sgrid%grid%y(sgrid%ni-i,j,k))
				sgrid%grid%z(sgrid%ni+i,j,k) 	= sgrid%grid%z(sgrid%ni,j,k) &
				+ (sgrid%grid%z(sgrid%ni,j,k) - sgrid%grid%z(sgrid%ni-i,j,k))
			endif
		enddo
		enddo
		enddo
		do k = 1-sgrid%ng,sgrid%nk+sgrid%ng
		do j = 1,sgrid%ng
		do i = 1-sgrid%ng,sgrid%ni+sgrid%ng
			if(sgrid%gnj==1) then
				sgrid%grid%x(i,1-j,k) 		= sgrid%grid%x(i,1,k)!-1.0
				sgrid%grid%y(i,1-j,k) 		= sgrid%grid%y(i,1,k)-1.0
				sgrid%grid%z(i,1-j,k) 		= sgrid%grid%z(i,1,k)!-1.0
				sgrid%grid%x(i,sgrid%nj+i,k) 	= sgrid%grid%x(i,sgrid%nj,k)!+1.0
				sgrid%grid%y(i,sgrid%nj+i,k) 	= sgrid%grid%y(i,sgrid%nj,k)+1.0
				sgrid%grid%z(i,sgrid%nj+i,k) 	= sgrid%grid%z(i,sgrid%nj,k)!+1.0
			else
				sgrid%grid%x(i,1-j,k) 		= sgrid%grid%x(i,1,k) &
				- (sgrid%grid%x(i,j+1,k) 	- sgrid%grid%x(i,1,k))
				sgrid%grid%y(i,1-j,k) 		= sgrid%grid%y(i,1,k) &
				- (sgrid%grid%y(i,j+1,k) 	- sgrid%grid%y(i,1,k))
				sgrid%grid%z(i,1-j,k) 		= sgrid%grid%z(i,1,k) &
				- (sgrid%grid%z(i,j+1,k) 	- sgrid%grid%z(i,1,k))
				sgrid%grid%x(i,sgrid%nj+j,k) 	= sgrid%grid%x(i,sgrid%nj,k)&
				+ (sgrid%grid%x(i,sgrid%nj,k) - sgrid%grid%x(i,sgrid%nj-j,k))
				sgrid%grid%y(i,sgrid%nj+j,k) 	= sgrid%grid%y(i,sgrid%nj,k)&
				+ (sgrid%grid%y(i,sgrid%nj,k) - sgrid%grid%y(i,sgrid%nj-j,k))
				sgrid%grid%z(i,sgrid%nj+j,k) 	= sgrid%grid%z(i,sgrid%nj,k)&
				+ (sgrid%grid%z(i,sgrid%nj,k) - sgrid%grid%z(i,sgrid%nj-j,k))
			endif
		enddo
		enddo
		enddo
		do k = 1,sgrid%ng
		do j = 1-sgrid%ng,sgrid%nj+sgrid%ng
		do i = 1-sgrid%ng,sgrid%ni+sgrid%ng
			if(sgrid%gnk==1) then
				sgrid%grid%x(i,j,1-k) 		= sgrid%grid%x(i,j,1)!-1.0
				sgrid%grid%y(i,j,1-k) 		= sgrid%grid%y(i,j,1)!-1.0
				sgrid%grid%z(i,j,1-k) 		= sgrid%grid%z(i,j,1)-1.0
				sgrid%grid%x(i,j,sgrid%nk+k) 	= sgrid%grid%x(i,j,sgrid%nk)! +1.0
				sgrid%grid%y(i,j,sgrid%nk+k) 	= sgrid%grid%y(i,j,sgrid%nk)! +1.0
				sgrid%grid%z(i,j,sgrid%nk+k) 	= sgrid%grid%z(i,j,sgrid%nk) +1.0
			else
				sgrid%grid%x(i,j,1-k) 		= sgrid%grid%x(i,j,1) &
				- (sgrid%grid%x(i,j,k+1) 	- sgrid%grid%x(i,j,1))
				sgrid%grid%y(i,j,1-k) 		= sgrid%grid%y(i,j,1) &
				- (sgrid%grid%y(i,j,k+1) 	- sgrid%grid%y(i,j,1))
				sgrid%grid%z(i,j,1-k) 		= sgrid%grid%z(i,j,1) &
				- (sgrid%grid%z(i,j,k+1) 	- sgrid%grid%z(i,j,1))
				sgrid%grid%x(i,j,sgrid%nk+k) 	= sgrid%grid%x(i,j,sgrid%nk) &
				+ (sgrid%grid%x(i,j,sgrid%nk) - sgrid%grid%x(i,j,sgrid%nk-k))
				sgrid%grid%y(i,j,sgrid%nk+k) 	= sgrid%grid%y(i,j,sgrid%nk) &
				+ (sgrid%grid%y(i,j,sgrid%nk) - sgrid%grid%y(i,j,sgrid%nk-k))
				sgrid%grid%z(i,j,sgrid%nk+k) 	= sgrid%grid%z(i,j,sgrid%nk) &
				+ (sgrid%grid%z(i,j,sgrid%nk) - sgrid%grid%z(i,j,sgrid%nk-k))
			endif
		enddo
		enddo
		enddo

		!Exchange to set ghost coordinates:
		call exchange_sdata(sgrid%scomm_max_r1,r1=sgrid%grid)

		!Finally, if you are using periodicity with just one grid point,
		!make sure the ghosts are offset somewhat.  This gets rid of problems
		!when dx,dy,dz=0 for interp and gradient calcs.
		if (sgrid%gni ==1 .AND. sgrid%bc_list(1) == LCS_PERIODIC) then
			sgrid%grid%x(0,:,:) = sgrid%grid%x(1,:,:) -1.0_LCSRP
			sgrid%grid%x(2,:,:) = sgrid%grid%x(1,:,:) +1.0_LCSRP
		endif
		if (sgrid%gnj ==1 .AND. sgrid%bc_list(2) == LCS_PERIODIC) then
			sgrid%grid%y(:,0,:) = sgrid%grid%y(:,1,:) -1.0_LCSRP
			sgrid%grid%y(:,2,:) = sgrid%grid%y(:,1,:) +1.0_LCSRP
		endif
		if (sgrid%gnk ==1 .AND. sgrid%bc_list(3) == LCS_PERIODIC) then
			sgrid%grid%z(:,:,0) = sgrid%grid%z(:,:,1) -1.0_LCSRP
			sgrid%grid%z(:,:,2) = sgrid%grid%z(:,:,1) +1.0_LCSRP
		endif

		!Check the rectilinear:
		call check_rectilinear(sgrid)

		!Compute the (geometric) least squares gradient construction weights
		if(.NOT. sgrid%rectilinear) then
			call compute_lsg_wts(sgrid)
		endif

	end subroutine init_sgrid

	subroutine new_sgrid_from_sgrid(sgrid_new,sgrid,label,res)
		implicit none
		!-----
		type(sgrid_t),pointer:: sgrid_new
		type(sgrid_t),pointer:: sgrid
		character(len=*):: label
		integer(LCSIP):: res
		!-----
		integer:: n_new(3), offset_new(3), gn_new(3)
		integer:: i,j,k,ii,jj,kk,i_new,j_new,k_new,ir,jr,kr
		real(LCSRP):: dx,dy,dz
		integer:: restest,ierr
		real(LCSRP),allocatable:: x(:,:,:)
		real(LCSRP),allocatable:: y(:,:,:)
		real(LCSRP),allocatable:: z(:,:,:)
		integer(LCSIP),dimension(6):: bc_list
		real(LCSRP):: lperiodic(3)
		!-----
		!Add or remove points to create a new sgrid structure
		!If res > 1, we add points
		!If res < 1, remove points
		!-----

		if(lcsrank==0)&
			write(*,*) 'In new_sgrid_from_sgrid: ',trim(sgrid%label),'=>',trim(label)

		!A couple checks...
		if (sgrid%ni >1 .AND. abs(res) >= sgrid%ni) then
			if (lcsrank==0) write(*,*) 'WARN:  cannot remove enough grid pts in x.  Setting Res=0'
			res = 0
		endif
		if (sgrid%nj >1 .AND. abs(res) >= sgrid%nj) then
			if (lcsrank==0) write(*,*) 'WARN:  cannot remove enough grid pts in y.  Setting Res=0'
			res = 0
		endif
		if (sgrid%nk >1 .AND. abs(res) >= sgrid%nk) then
			if (lcsrank==0) write(*,*) 'WARN:  cannot remove enough grid pts in z.  Setting Res=0'
			res = 0
		endif
		call MPI_ALLREDUCE(abs(res),restest,1,MPI_INTEGER,MPI_MIN,lcscomm,ierr)
		if(restest ==0) res = 0


		if (res == 0) then
			if(lcsrank==0)&
				 write(*,*) 'WARN: Resolution = 0, new_sgrid_from_sgrid does not need to be called'
			sgrid_new => sgrid
			return
		endif

		!-----
		!Determine the number of points in each direction and the global offsets:
		!If n = 1 in any direction, just keep the same grid for that coordinate.
		!In the event of globally non-periodic conditions, dont insert any new points
		!outside of the original grid
		!-----
		n_new(1:3) = 0
		gn_new(1:3) = 0
		offset_new(1:3) = 0
		!i
		if (sgrid%ni==1) then
			n_new(1) = 1
			offset_new(1) = sgrid%offset_i
		else
			!insert points
			if (res > 0) then
				do i = 1,sgrid%gni
					gn_new(1)=gn_new(1)+1
					do ii = 1,res
						if(sgrid%global_bc_list(4)/=LCS_PERIODIC .AND. i == sgrid%gni) cycle
						gn_new(1) =gn_new(1)+1
					enddo
					if(i > sgrid%offset_i .AND. i <= sgrid%offset_i + sgrid%ni) then
						n_new(1) = n_new(1) + 1
						do ii = 1,res
							if(sgrid%global_bc_list(4)/=LCS_PERIODIC .AND. i == sgrid%gni) cycle
							n_new(1) = n_new(1) + 1
						enddo
					elseif(i <= sgrid%offset_i) then
						offset_new(1) = gn_new(1)
					endif
				enddo
			endif
			!remove points
			if (res < 0) then
				do i = 1,sgrid%gni
					if(mod(i,abs(res)+1)==0) then
						if(i > sgrid%offset_i .AND. i <= sgrid%offset_i+sgrid%ni) then
							n_new(1)=n_new(1)+1
						endif
						if(i <=sgrid%offset_i) then
							offset_new(1) = offset_new(1) + 1
						endif
					endif
				enddo
			endif

		endif
		!j
		if (sgrid%nj==1) then
			n_new(2) = 1
			offset_new(2) = sgrid%offset_j
		else
			!insert points
			if (res > 0) then
				do j = 1,sgrid%gnj
					gn_new(2)=gn_new(2)+1
					do jj = 1,res
						if(sgrid%global_bc_list(5)/=LCS_PERIODIC .AND. j == sgrid%gnj) cycle
						gn_new(2) =gn_new(2)+1
					enddo
					if(j > sgrid%offset_j .AND. j <= sgrid%offset_j + sgrid%nj) then
						n_new(2) = n_new(2) + 1
						do jj = 1,res
							if(sgrid%global_bc_list(5)/=LCS_PERIODIC .AND. j == sgrid%gnj) cycle
							n_new(2) = n_new(2) + 1
						enddo
					elseif(j <= sgrid%offset_j) then
						offset_new(2) = gn_new(2)
					endif
				enddo
			endif
			!remove points
			if (res < 0) then
				do j = 1,sgrid%gnj
					if(mod(j,abs(res)+1)==0) then
						if(j > sgrid%offset_j .AND. j <= sgrid%offset_j+sgrid%nj) then
							n_new(2)=n_new(2)+1
						endif
						if(j <=sgrid%offset_j) then
							offset_new(2) = offset_new(2) + 1
						endif
					endif
				enddo
			endif
		endif
		!k
		if (sgrid%nk==1) then
			n_new(3) = 1
			offset_new(3) = sgrid%offset_k
		else
			!insert points
			if (res > 0) then
				do k = 1,sgrid%gnk
					gn_new(3)=gn_new(3)+1
					do kk = 1,res
						if(sgrid%global_bc_list(6)/=LCS_PERIODIC .AND. k == sgrid%gnk) cycle
						gn_new(3) =gn_new(3)+1
					enddo
					if(k> sgrid%offset_k .AND. k <= sgrid%offset_k + sgrid%nk) then
						n_new(3) = n_new(3) + 1
						do kk = 1,res
							if(sgrid%global_bc_list(6)/=LCS_PERIODIC .AND. k == sgrid%gnk) cycle
							n_new(3) = n_new(3) + 1
						enddo
					elseif(k <= sgrid%offset_k) then
						offset_new(3) = gn_new(3)
					endif
				enddo
			endif
			!remove points
			if (res < 0) then
				do k = 1,sgrid%gnk
					if(mod(k,abs(res)+1)==0) then
						if(k > sgrid%offset_k .AND. k <= sgrid%offset_k+sgrid%nk) then
							n_new(3)=n_new(3)+1
						endif
						if(k <=sgrid%offset_k) then
							offset_new(3) = offset_new(3) + 1
						endif
					endif
				enddo
			endif
		endif

		!-----
		!Now populate the New X,Y,Z arrays:
		!-----
		allocate(x(1:n_new(1),1:n_new(2),1:n_new(3)))
		allocate(y(1:n_new(1),1:n_new(2),1:n_new(3)))
		allocate(z(1:n_new(1),1:n_new(2),1:n_new(3)))

		!-----
		!Case of adding points:
		!Strictly speaking, this will only work well for regular cuboid grids.
		!Could generalize to other shapes in the future.
		!-----
		if (res > 0) then
			do k = 1,sgrid%gnk
			do j = 1,sgrid%gnj
			do i = 1,sgrid%gni
				if(i > sgrid%offset_i .AND. i <= sgrid%offset_i + sgrid%ni) then
				if(j > sgrid%offset_j .AND. j <= sgrid%offset_j + sgrid%nj) then
				if(k > sgrid%offset_k .AND. k <= sgrid%offset_k + sgrid%nk) then
					ii = i - sgrid%offset_i !index on this proc
					jj = j - sgrid%offset_j !index on this proc
					kk = k - sgrid%offset_k !index on this proc
					!insert new pts within existing grid:
					do kr = 0,res
					do jr = 0,res
					do ir = 0,res
						!i_new = ii + ir; j_new = jj + jr; k_new = kk + kr
						i_new = ii + (ii-1)*res +ir;
						j_new = jj + (jj-1)*res +jr;
						k_new = kk + (kk-1)*res +kr;
						!Spacing between the new grid points in each dir:  Assumes ng >= 1
						dx = (sgrid%grid%x(ii+1,jj,kk) - sgrid%grid%x(ii,jj,kk)) / real(res +1)
						dy = (sgrid%grid%y(ii,jj+1,kk) - sgrid%grid%y(ii,jj,kk)) / real(res +1)
						dz = (sgrid%grid%z(ii,jj,kk+1) - sgrid%grid%z(ii,jj,kk)) / real(res +1)
						!New grid points
						if (i_new <=n_new(1) .and. j_new <= n_new(2) .and. k_new <= n_new(3)) then
							x(i_new,j_new,k_new) = sgrid%grid%x(ii,jj,kk) + dx*real(ir)
							y(i_new,j_new,k_new) = sgrid%grid%y(ii,jj,kk) + dy*real(jr)
							z(i_new,j_new,k_new) = sgrid%grid%z(ii,jj,kk) + dz*real(kr)
						endif
					enddo
					enddo
					enddo
				endif
				endif
				endif
			enddo
			enddo
			enddo
		endif

		!-----
		!Case of removing points:
		!If the original grid has n=1 in any direction, keep the same
		!-----
		if(res < 0) then
			i_new = 0; j_new=0; k_new = 0
			do k = 1,sgrid%gnk
				i_new = 0; j_new = 0
				if(sgrid%nk==1) then
					k_new = 1
					kk = 1
				elseif(mod(k,abs(res)+1)==0 .AND. k > sgrid%offset_k .AND. k <= sgrid%offset_k+sgrid%nk) then
					k_new = k_new + 1 !new index
					kk = k - sgrid%offset_k !index on old sgrid
				else
					cycle
				endif
				do j = 1,sgrid%gnj
					i_new = 0
					if(sgrid%nj==1) then
						j_new = 1
						jj = 1
					elseif(mod(j,abs(res)+1)==0 .AND. j > sgrid%offset_j .AND. j <= sgrid%offset_j+sgrid%nj) then
						j_new = j_new + 1 !new index
						jj = j - sgrid%offset_j !index on old sgrid
					else
						cycle
					endif
					do i = 1,sgrid%gni
						if(sgrid%ni==1) then
							i_new = 1
							ii = 1
						elseif(mod(i,abs(res)+1)==0 .AND. i > sgrid%offset_i .AND. i <= sgrid%offset_i+sgrid%ni) then
							i_new = i_new + 1 !new index
							ii = i - sgrid%offset_i !index on old sgrid
						else
							cycle
						endif
						!This grid point is saved:
						x(i_new,j_new,k_new) = sgrid%grid%x(ii,jj,kk)
						y(i_new,j_new,k_new) = sgrid%grid%y(ii,jj,kk)
						z(i_new,j_new,k_new) = sgrid%grid%z(ii,jj,kk)
					enddo
				enddo
			enddo
		endif

		!-----
		!Copy the other properties and initialize another sgrid structure
		!-----
		bc_list = sgrid%global_bc_list
		lperiodic = sgrid%lperiodic
		call init_sgrid(sgrid_new,label,n_new,offset_new,x,y,z,bc_list,lperiodic)

		!-----
		!Deallocate/cleanup
		!-----
		deallocate(x)
		deallocate(y)
		deallocate(z)

	end subroutine new_sgrid_from_sgrid



	subroutine destroy_sgrid(sgrid)
		implicit none
		!-----
		type(sgrid_t):: sgrid
		!-----

		if(lcsrank==0)&
			write(*,*) 'in destroy_sgrid...',sgrid%label

		sgrid%ni = 0
		sgrid%nj = 0
		sgrid%nk = 0
		sgrid%ng = 0
		sgrid%gni = 0
		sgrid%gnj = 0
		sgrid%gnk = 0
		sgrid%offset_i = 0
		sgrid%offset_j = 0
		sgrid%offset_k = 0

		call destroy_scomm(sgrid%scomm_face_r0)
		call destroy_scomm(sgrid%scomm_max_r0)
		call destroy_scomm(sgrid%scomm_face_r1)
		call destroy_scomm(sgrid%scomm_max_r1)
		call destroy_scomm(sgrid%scomm_face_r2)
		call destroy_scomm(sgrid%scomm_max_r2)

		call destroy_sr1(sgrid%grid)

		sgrid%label = 'Unused sgrid'

	end subroutine destroy_sgrid

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

	subroutine grad_sr1(sgrid,sr1,grad)
		implicit none
		!----
		type(sgrid_t),intent(in):: sgrid
		type(sr1_t),intent(in):: sr1
		type(sr2_t),intent(inout):: grad
		!----
		integer:: i,j,k
		!----

		if (lcsrank==0 .AND. LCS_VERBOSE) &
			write(*,*) 'in grad_sr1... ',trim(sr1%label),' => ',trim(grad%label)

		!Check that the grid is rectilinear:
		if(.NOT. sgrid%rectilinear) then
			if(lcsrank==0) &
				write(*,*) 'WARNING, in grad_sr1: sgrid is not rectilinear, calling grad_sr1_ls instead'
			call grad_sr1_ls(sgrid,sr1,grad)
			return
		endif

		!2nd order central scheme:
		do k = 1,sgrid%nk
		do j = 1,sgrid%nj
		do i = 1,sgrid%ni
			grad%xx(i,j,k) = (sr1%x(i+1,j,k)-sr1%x(i-1,j,k)) / (sgrid%grid%x(i+1,j,k) - sgrid%grid%x(i-1,j,k))
			grad%xy(i,j,k) = (sr1%x(i,j+1,k)-sr1%x(i,j-1,k)) / (sgrid%grid%y(i,j+1,k) - sgrid%grid%y(i,j-1,k))
			grad%xz(i,j,k) = (sr1%x(i,j,k+1)-sr1%x(i,j,k-1)) / (sgrid%grid%z(i,j,k+1) - sgrid%grid%z(i,j,k-1))

			grad%yx(i,j,k) = (sr1%y(i+1,j,k)-sr1%y(i-1,j,k)) / (sgrid%grid%x(i+1,j,k) - sgrid%grid%x(i-1,j,k))
			grad%yy(i,j,k) = (sr1%y(i,j+1,k)-sr1%y(i,j-1,k)) / (sgrid%grid%y(i,j+1,k) - sgrid%grid%y(i,j-1,k))
			grad%yz(i,j,k) = (sr1%y(i,j,k+1)-sr1%y(i,j,k-1)) / (sgrid%grid%z(i,j,k+1) - sgrid%grid%z(i,j,k-1))

			grad%zx(i,j,k) = (sr1%z(i+1,j,k)-sr1%z(i-1,j,k)) / (sgrid%grid%x(i+1,j,k) - sgrid%grid%x(i-1,j,k))
			grad%zy(i,j,k) = (sr1%z(i,j+1,k)-sr1%z(i,j-1,k)) / (sgrid%grid%y(i,j+1,k) - sgrid%grid%y(i,j-1,k))
			grad%zz(i,j,k) = (sr1%z(i,j,k+1)-sr1%z(i,j,k-1)) / (sgrid%grid%z(i,j,k+1) - sgrid%grid%z(i,j,k-1))
		enddo
		enddo
		enddo

	end subroutine grad_sr1

	subroutine compute_lsg_wts(sgrid)
		implicit none
		!-----
		type(sgrid_t):: sgrid
		!-----
		integer:: i,j,k,ii,jj,kk
		integer:: nbr,nbr_f,nbr_l
		real(lcsrp):: swdx2,swdy2,swdz2,swdxdy,swdxdz,swdydz,weight,dx(3),denom
		character(len=32):: label
		!-----
		if (lcsrank == 0) &
			write(*,*) 'in calc_lsg_weights...', trim(sgrid%label)

		!Determine the nbr range depending on the desired connectivity
		if(FULL_GRADIENT_CONNECTIVITY) then
			nbr_f= 2
			nbr_l= 27
		else
			nbr_f=2
			nbr_l=7
		endif

		!Allocate space for the weights:
		allocate(sgrid%lsg_wts(nbr_f:nbr_l))
		do nbr = nbr_f,nbr_l
			write(label,'(a,i2.2)') 'LSG_WTS_',nbr
			call init_sr1(sgrid%lsg_wts(nbr),sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,trim(label),translate=.false.)
		enddo

		do k= 1,sgrid%nk
		do j= 1,sgrid%nj
		do i= 1,sgrid%ni
			swdx2 = 0.0_LCSRP
			swdy2 = 0.0_LCSRP
			swdz2 = 0.0_LCSRP
			swdxdy = 0.0_LCSRP
			swdxdz = 0.0_LCSRP
			swdydz = 0.0_LCSRP

			do nbr = nbr_f,nbr_l
				ii = i+NBR_OFFSET(1,nbr)
				jj = j+NBR_OFFSET(2,nbr)
				kk = k+NBR_OFFSET(3,nbr)

				dx(1) = sgrid%grid%x(ii,jj,kk) - sgrid%grid%x(i,j,k)
				dx(2) = sgrid%grid%y(ii,jj,kk) - sgrid%grid%y(i,j,k)
				dx(3) = sgrid%grid%z(ii,jj,kk) - sgrid%grid%z(i,j,k)
				!Inverse distance weight:
				weight = 1.0_LCSRP/sum(dx(1:3)**2)

				swdx2 = swdx2 + weight*dx(1)**2
				swdy2 = swdy2 + weight*dx(2)**2
				swdz2 = swdz2 + weight*dx(3)**2
				swdxdy = swdxdy + weight*dx(1)*dx(2)
				swdxdz = swdxdz + weight*dx(1)*dx(3)
				swdydz = swdydz + weight*dx(2)*dx(3)
			enddo

			denom = 2.0_LCSRP*swdxdy*swdxdz*swdydz + &
				swdx2*swdy2*swdz2 - &
				swdx2*swdydz**2 - &
				swdy2*swdxdz**2 - &
				swdz2*swdxdy**2

			do nbr = nbr_f,nbr_l
				ii = i+NBR_OFFSET(1,nbr)
				jj = j+NBR_OFFSET(2,nbr)
				kk = k+NBR_OFFSET(3,nbr)
				dx(1) = sgrid%grid%x(ii,jj,kk) - sgrid%grid%x(i,j,k)
				dx(2) = sgrid%grid%y(ii,jj,kk) - sgrid%grid%y(i,j,k)
				dx(3) = sgrid%grid%z(ii,jj,kk) - sgrid%grid%z(i,j,k)

				!Inverse distance weight:
				weight = 1.0_LCSRP/sum(dx(1:3)**2)
				! x
				sgrid%lsg_wts(nbr)%x(i,j,k) = weight*( &
					(swdy2*swdz2-swdydz**2)*dx(1) + &
					(swdxdz*swdydz-swdxdy*swdz2)*dx(2) + &
					(swdxdy*swdydz-swdxdz*swdy2)*dx(3) )/denom
				! y
				sgrid%lsg_wts(nbr)%y(i,j,k) = weight*( &
					(swdxdz*swdydz-swdxdy*swdz2)*dx(1) + &
					(swdx2*swdz2-swdxdz**2)*dx(2) + &
					(swdxdy*swdxdz-swdydz*swdx2)*dx(3) )/denom
				! z
				sgrid%lsg_wts(nbr)%z(i,j,k) = weight*( &
					(swdxdy*swdydz-swdxdz*swdy2)*dx(1) + &
					(swdxdy*swdxdz-swdydz*swdx2)*dx(2) + &
					(swdx2*swdy2-swdxdy**2)*dx(3) )/denom
			end do
		enddo
		enddo
		enddo

	end subroutine compute_lsg_wts

	subroutine grad_sr1_ls(sgrid,sr1,grad)
		implicit none
		!----
		type(sgrid_t),intent(in):: sgrid
		type(sr1_t),intent(in):: sr1
		type(sr2_t),intent(inout):: grad
		!----
		integer:: i,j,k
		integer:: ni,nj,nk,ng
		integer:: nbr,nbr_f,nbr_l
		type(sr1_t):: tmp
		!----

		if (lcsrank==0 .AND. LCS_VERBOSE) &
			write(*,*) 'in grad_sr1_ls... ',trim(sr1%label),' => ',trim(grad%label),FULL_GRADIENT_CONNECTIVITY

		if(FULL_GRADIENT_CONNECTIVITY) then
			nbr_f= 2
			nbr_l= 27
		else
			nbr_f=2
			nbr_l=7
		endif

		grad%xx = 0.0_LCSRP
		grad%xy = 0.0_LCSRP
		grad%xz = 0.0_LCSRP
		grad%yx = 0.0_LCSRP
		grad%yy = 0.0_LCSRP
		grad%yz = 0.0_LCSRP
		grad%zx = 0.0_LCSRP
		grad%zy = 0.0_LCSRP
		grad%zz = 0.0_LCSRP

		ni = sgrid%ni
		nj = sgrid%nj
		nk = sgrid%nk
		ng = sgrid%ng
		call init_sr1(tmp,ni,nj,nk,ng,'TMP',translate=.false.)
		!These should all vectorize (confirmed with gfortran)
		do nbr = nbr_f,nbr_l
			i = NBR_OFFSET(1,nbr)
			j = NBR_OFFSET(2,nbr)
			k = NBR_OFFSET(3,nbr)
			tmp%x(1:ni,1:nj,1:nk) = sr1%x(1+i:ni+i, 1+j:nj+j, 1+k:nk+k)
			tmp%y(1:ni,1:nj,1:nk) = sr1%y(1+i:ni+i, 1+j:nj+j, 1+k:nk+k)
			tmp%z(1:ni,1:nj,1:nk) = sr1%z(1+i:ni+i, 1+j:nj+j, 1+k:nk+k)
			tmp%x = tmp%x-sr1%x
			grad%xx = grad%xx + tmp%x*sgrid%lsg_wts(nbr)%x
			grad%xy = grad%xy + tmp%x*sgrid%lsg_wts(nbr)%y
			grad%xz = grad%xz + tmp%x*sgrid%lsg_wts(nbr)%z
			tmp%y = tmp%y-sr1%y
			grad%yx = grad%yx + tmp%y*sgrid%lsg_wts(nbr)%x
			grad%yy = grad%yy + tmp%y*sgrid%lsg_wts(nbr)%y
			grad%yz = grad%yz + tmp%y*sgrid%lsg_wts(nbr)%z
			tmp%z = tmp%z-sr1%z
			grad%zx = grad%zx + tmp%z*sgrid%lsg_wts(nbr)%x
			grad%zy = grad%zy + tmp%z*sgrid%lsg_wts(nbr)%y
			grad%zz = grad%zz + tmp%z*sgrid%lsg_wts(nbr)%z
		enddo
		call destroy_sr1(tmp)

	end subroutine grad_sr1_ls

	subroutine check_rectilinear(sgrid)
		implicit none
		!-----
		type(sgrid_t):: sgrid	
		!-----
		integer:: i,j,k,ni,nj,nk
		logical:: ortho_x,ortho_y, ortho_z
		real(LCSRP):: biggest_dim
		real(LCSRP):: xmax,xmin,ymax,ymin,zmax,zmin
		real(LCSRP),parameter:: TOL =1e-4
		!-----

		if(lcsrank==0) &	
			write(*,*) 'In check_rectilinear... ', trim(sgrid%label)

		ni = sgrid%ni
		nj = sgrid%nj
		nk = sgrid%nk
		
		biggest_dim = max( maxval(sgrid%grid%x)-minval(sgrid%grid%x),&
					 maxval(sgrid%grid%y)-minval(sgrid%grid%y),&
					 maxval(sgrid%grid%z)-minval(sgrid%grid%z))

		!x
		ortho_x  = .true.
		xloop:do k = 1,nk
		do j = 1,nj
			xmax = maxval(sgrid%grid%x(1:ni,j,k))
			xmin = minval(sgrid%grid%x(1:ni,j,k))
			ymax = maxval(sgrid%grid%y(1:ni,j,k))
			ymin = minval(sgrid%grid%y(1:ni,j,k))
			zmax = maxval(sgrid%grid%z(1:ni,j,k))
			zmin = minval(sgrid%grid%z(1:ni,j,k))
			if( abs(ymax - ymin) / biggest_dim > TOL .OR. abs(zmax-zmin)/biggest_dim > TOL) then
				ortho_x = .false.
				exit xloop
			endif
		enddo
		enddo xloop
		!y
		ortho_y  = .true.
		yloop:do k = 1,nk
		do i = 1,ni
			xmax = maxval(sgrid%grid%x(i,1:nj,k))
			xmin = minval(sgrid%grid%x(i,1:nj,k))
			ymax = maxval(sgrid%grid%y(i,1:nj,k))
			ymin = minval(sgrid%grid%y(i,1:nj,k))
			zmax = maxval(sgrid%grid%z(i,1:nj,k))
			zmin = minval(sgrid%grid%z(i,1:nj,k))
			if( abs(xmax - xmin) / biggest_dim > TOL .OR. abs(zmax-zmin)/biggest_dim > TOL) then
				ortho_y = .false.
				exit yloop
			endif
		enddo
		enddo yloop
		!z
		ortho_z  = .true.
		zloop:do j = 1,nj
		do i = 1,ni
			xmax = maxval(sgrid%grid%x(i,j,1:k))
			xmin = minval(sgrid%grid%x(i,j,1:k))
			ymax = maxval(sgrid%grid%y(i,j,1:k))
			ymin = minval(sgrid%grid%y(i,j,1:k))
			zmax = maxval(sgrid%grid%z(i,j,1:k))
			zmin = minval(sgrid%grid%z(i,j,1:k))
			if( abs(xmax - xmin) / biggest_dim > TOL .OR. abs(ymax-ymin)/biggest_dim > TOL) then
				ortho_z = .false.
				exit zloop
			endif
		enddo
		enddo zloop

		if(ortho_x .and. ortho_y .and. ortho_z) then
			sgrid%rectilinear = .true.
			if(lcsrank==0) write(*,*) 'Grid ', trim(sgrid%label), ' IS rectilinear'
		else
			sgrid%rectilinear = .false.
			if(lcsrank==0) write(*,*) 'Grid ', trim(sgrid%label), ' IS NOT rectilinear', ortho_x, ortho_y, ortho_z, biggest_dim
		endif

	end subroutine check_rectilinear

end module sgrid_m
