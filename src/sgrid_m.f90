!
! Routines associated with structured grids
!
module sgrid_m
	use data_m
	use structured_m
	use comms_m
	use gradient_m
	use geometry_m
	use io_m
	contains
	subroutine init_sgrid(sgrid,label,n,offset,x,y,z,flag)
		implicit none
		!-----
		type(sgrid_t),pointer:: sgrid
		character(len=*):: label
		integer:: n(3),offset(3)
		real(LCSRP):: x(1:n(1),1:n(2),1:n(3))
		real(LCSRP):: y(1:n(1),1:n(2),1:n(3))
		real(LCSRP):: z(1:n(1),1:n(2),1:n(3))
		integer(LCSIP):: flag(1:n(1),1:n(2),1:n(3))
		!-----
		integer:: i,j,k,ii,jj,kk,isg,ilcs,ierr
		integer:: my_nmax(3), nmax(3), offset_min(3)
		type(sgrid_t),allocatable:: sgrid_c_tmp(:)
		integer:: scfdptr
		integer,allocatable:: lcsptr(:)
		logical:: periodic(6),my_periodic(6)
		real(LCSRP):: v1(3),v2(3),delta(3),mag
		integer:: i_ib, j_ib, k_ib, i_b, j_b, k_b
		integer:: im1,km1,jm1,ip1,jp1,kp1
		type(sr1_t):: fake
		integer:: ni,nj,nk,ng
		real(LCSRP):: periodic_shift(-1:1,-1:1,-1:1,1:3)
		real(LCSRP),allocatable:: mybb(:,:)
		real(LCSRP),allocatable:: myps(:,:,:,:,:)
		type(sr0_t):: bc
		real(LCSRP):: myscaling,scaling
		integer:: gn(3)
		type(sr0_t):: tmp
		character(len=128):: fname
		!-----

		if(lcsrank==0) &
			write(*,*) 'in init_sgrid... ', trim(label)

		!brevity
		ni =n(1)
		nj =n(2)
		nk =n(3)
		ng =NGHOST_CFD

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
		!	5.  If gni, gnj, or gnk = 1, then we cannot have anything but LCS_INTERNAL BC for those nodes
		!
		sgrid%ng = ng
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
		!Check for global periodicity in i,j,k directions:
		!Note, we dont allow periodicity with 1 node in i,j,or k.
		!
		my_periodic(1:6) = .false.
		if (sgrid%offset_i ==0 .and. any(flag(1,:,:)==LCS_INTERNAL)) my_periodic(1) = .true.
		if (sgrid%offset_j ==0 .and. any(flag(:,1,:)==LCS_INTERNAL)) my_periodic(2) = .true.
		if (sgrid%offset_k ==0 .and. any(flag(:,:,1)==LCS_INTERNAL)) my_periodic(3) = .true.
		if (sgrid%offset_i+ni ==sgrid%gni .and. any(flag(ni,:,:)==LCS_INTERNAL)) my_periodic(4) = .true.
		if (sgrid%offset_j+nj ==sgrid%gnj .and. any(flag(:,nj,:)==LCS_INTERNAL)) my_periodic(5) = .true.
		if (sgrid%offset_k+nk ==sgrid%gnk .and. any(flag(:,:,nk)==LCS_INTERNAL)) my_periodic(6) = .true.
		call MPI_ALLREDUCE(my_periodic,periodic,6,MPI_LOGICAL,MPI_LOR,lcscomm,ierr)
		if(periodic(1) .NEQV. periodic(4))then
			write(*,*) 'ERROR:  X Periodicity does not appear to match.  Check flag values passed to cfd2lcs_init.',periodic
			CFD2LCS_ERROR = 1
			return
		endif
		if(periodic(2) .NEQV. periodic(5))then
			write(*,*) 'ERROR:  Y Periodicity does not appear to match.  Check flag values passed to cfd2lcs_init.',periodic
			CFD2LCS_ERROR = 1
			return
		endif
		if(periodic(3) .NEQV. periodic(6))then
			write(*,*) 'ERROR:  Z Periodicity does not appear to match.  Check flag values passed to cfd2lcs_init.',periodic
			CFD2LCS_ERROR = 1
			return
		endif
		sgrid%periodic_i = periodic(1)
		sgrid%periodic_j = periodic(2)
		sgrid%periodic_k = periodic(3)
		!NO comm if only one node in a given direction:
		if(sgrid%gni==1) sgrid%periodic_i = .false.
		if(sgrid%gnj==1) sgrid%periodic_j = .false.
		if(sgrid%gnk==1) sgrid%periodic_k = .false.
		if(lcsrank==0) write(*,*) 'lcsrank[',lcsrank, '] periodic=', sgrid%periodic_i,sgrid%periodic_j,sgrid%periodic_k

		!
		!Initialize the communication patterns
		!
		call init_scomm(sgrid%scomm_face_r0,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,sgrid%periodic_i,&
			sgrid%periodic_j,sgrid%periodic_k,FACE_CONNECT,R0_COMM,'R0 face-nbr comms' )
		call init_scomm(sgrid%scomm_max_r0,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,sgrid%periodic_i,&
			sgrid%periodic_j,sgrid%periodic_k,MAX_CONNECT,R0_COMM,'R0 max-nbr comms' )
		call init_scomm(sgrid%scomm_face_r1,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,sgrid%periodic_i,&
			sgrid%periodic_j,sgrid%periodic_k,FACE_CONNECT,R1_COMM,'R1 face-nbr comms' )
		call init_scomm(sgrid%scomm_max_r1,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,sgrid%periodic_i,&
			sgrid%periodic_j,sgrid%periodic_k,MAX_CONNECT,R1_COMM,'R1 max-nbr comms' )
		call init_scomm(sgrid%scomm_face_r2,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,sgrid%periodic_i,&
			sgrid%periodic_j,sgrid%periodic_k,FACE_CONNECT,R2_COMM,'R2 face-nbr comms' )
		call init_scomm(sgrid%scomm_max_r2,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,&
			sgrid%offset_i,sgrid%offset_j,sgrid%offset_k,sgrid%periodic_i,&
			sgrid%periodic_j,sgrid%periodic_k,MAX_CONNECT,R2_COMM,'R2 max-nbr comms' )

		!
		!Initialize the structured grid coordinates, BC flag, and BC inward normal:
		!
		call init_sr1(sgrid%grid,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,'GRID',translate=.true.)
		call init_si0(sgrid%bcflag,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,'BCFLAG')
		call init_sr1(sgrid%bcnorm,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,'BCNORM',translate=.false.)

		!Set interior points:
		do k=1,sgrid%nk
		do j=1,sgrid%nj
		do i=1,sgrid%ni
			sgrid%grid%x(i,j,k) = x(i,j,k)
			sgrid%grid%y(i,j,k) = y(i,j,k)
			sgrid%grid%z(i,j,k) = z(i,j,k)
			sgrid%bcflag%i(i,j,k) = flag(i,j,k)
		enddo
		enddo
		enddo

		!Set fake boundary coordinates.
		!These will be over-written in the case of periodicity.
		call init_sr1(fake,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,'FAKE',translate=.false.)
		do k = 1-sgrid%ng,sgrid%nk+sgrid%ng
		do j = 1-sgrid%ng,sgrid%nj+sgrid%ng
		do i = 1-sgrid%ng,sgrid%ni+sgrid%ng

			!find the nearest boundary node:
			i_b = max(min(i,sgrid%ni),1)
			j_b = max(min(j,sgrid%nj),1)
			k_b = max(min(k,sgrid%nk),1)

			!If this node is not a fake, cycle:
			if(i==i_b .and. j==j_b .and. k==k_b) cycle

			!The interior: node
			i_ib = max(min(i_b+(i_b-i),sgrid%ni),1)
			j_ib = max(min(j_b+(j_b-j),sgrid%nj),1)
			k_ib = max(min(k_b+(k_b-k),sgrid%nk),1)

			!offset for the fake: x,y,z
			delta(1) = sgrid%grid%x(i_b,j_b,k_b) - sgrid%grid%x(i_ib,j_ib,k_ib)
			delta(2) = sgrid%grid%y(i_b,j_b,k_b) - sgrid%grid%y(i_ib,j_ib,k_ib)
			delta(3) = sgrid%grid%z(i_b,j_b,k_b) - sgrid%grid%z(i_ib,j_ib,k_ib)

			!Fake coord:
			sgrid%grid%x(i,j,k) = sgrid%grid%x(i_b,j_b,k_b) + delta(1)
			sgrid%grid%y(i,j,k) = sgrid%grid%y(i_b,j_b,k_b) + delta(2)
			sgrid%grid%z(i,j,k) = sgrid%grid%z(i_b,j_b,k_b) + delta(3)

			!Fake bcflag
			sgrid%bcflag%i(i,j,k) = sgrid%bcflag%i(i_b,j_b,k_b)
		enddo
		enddo
		enddo

		!Handle 2d conditions with 1 grid in i,j,or k:
		if(sgrid%gni==1)then
			do k = 1-ng,nk+ng
			do j = 1-ng,nj+ng
				if(sgrid%bcflag%i(1,j,k) == LCS_INTERNAL) then
					sgrid%bcflag%i(1-ng:0,j,k) = LCS_2D
					sgrid%bcflag%i(2:ni+ng,j,k) = LCS_2D
				else
					sgrid%bcflag%i(1-ng:0,j,k) = sgrid%bcflag%i(1,j,k)
					sgrid%bcflag%i(2:ni+ng,j,k) = sgrid%bcflag%i(1,j,k)
				endif
			enddo
			enddo
		endif
		if(sgrid%gnj==1)then
			do k = 1-ng,nk+ng
			do i = 1-ng,ni+ng
				if(sgrid%bcflag%i(i,1,k) == LCS_INTERNAL) then
					sgrid%bcflag%i(i,1-ng:0,k) = LCS_2D
					sgrid%bcflag%i(i,2:nj+ng,k) = LCS_2D
				else
					sgrid%bcflag%i(i,1-ng:0,k) = sgrid%bcflag%i(i,1,k)
					sgrid%bcflag%i(i,2:nj+ng,k) = sgrid%bcflag%i(i,1,k)
				endif
			enddo
			enddo
		endif
		if(sgrid%gnk==1)then
			do j = 1-ng,nj+ng
			do i = 1-ng,ni+ng
				if(sgrid%bcflag%i(i,j,1) == LCS_INTERNAL) then
					sgrid%bcflag%i(i,j,1-ng:0) = LCS_2D
					sgrid%bcflag%i(i,j,2:nk+ng) = LCS_2D
				else
					sgrid%bcflag%i(i,j,1-ng:0) = sgrid%bcflag%i(i,j,1)
					sgrid%bcflag%i(i,j,2:nk+ng) = sgrid%bcflag%i(i,j,1)
				endif
			enddo
			enddo
		endif

		!Exchange to set ghost coordinates, but copy the fakes before you do:
		fake%x = sgrid%grid%x
		fake%y = sgrid%grid%y
		fake%z = sgrid%grid%z
		call exchange_sdata(sgrid%scomm_max_r1,r1=sgrid%grid)


		!Figure out the periodic shift for each communication direction:
		!Note, we assume that the shift is constant for ALL nodes along a
		!given periodic direction.  We double check that this is true here:
		fake%x = sgrid%grid%x -fake%x
		fake%y = sgrid%grid%y -fake%y
		fake%z = sgrid%grid%z -fake%z
		do k = -1,1
			if(k==-1) kk = 0
			if(k==0) kk = 1
			if(k==1) kk = nk+1
		do j = -1,1
			if(j==-1)jj = 0
			if(j==0) jj = 1
			if(j==1) jj = nj+1
		do i = -1,1
			if(i==-1)ii = 0
			if(i==0) ii = 1
			if(i==1) ii = ni+1
			periodic_shift(i,j,k,1) = fake%x(ii,jj,kk)
			periodic_shift(i,j,k,2) = fake%y(ii,jj,kk)
			periodic_shift(i,j,k,3) = fake%z(ii,jj,kk)
		enddo
		enddo
		enddo
		sgrid%scomm_face_r0%periodic_shift = periodic_shift
		sgrid%scomm_face_r1%periodic_shift = periodic_shift
		sgrid%scomm_face_r2%periodic_shift = periodic_shift
		sgrid%scomm_max_r0%periodic_shift = periodic_shift
		sgrid%scomm_max_r1%periodic_shift = periodic_shift
		sgrid%scomm_max_r2%periodic_shift = periodic_shift

		!Exchange the grid again, to set the true fake coords in the
		!case of globally periodic bc:
		call exchange_sdata(sgrid%scomm_max_r1,r1=sgrid%grid)


		!Two dimensionality test:  In the case of one direction having only 1 node,
		!offset the fake by unit normal distance for each node.  This makes tracking,
		! and gradient calcs the same as in 3D.  JRF:  This delta is scaled by some
		!measure of the mean grid spacing, to keep the grid geometry nice.
		myscaling = max(&
			 maxval(sgrid%grid%x(1:ni,1:nj,1:nk)) - minval(sgrid%grid%x(1:ni,1:nj,1:nk)),&
			 maxval(sgrid%grid%y(1:ni,1:nj,1:nk)) - minval(sgrid%grid%y(1:ni,1:nj,1:nk)),&
			 maxval(sgrid%grid%z(1:ni,1:nj,1:nk)) - minval(sgrid%grid%z(1:ni,1:nj,1:nk)) &
			 ) / real(max(ni,nj,nk),LCSRP)
		call MPI_ALLREDUCE(myscaling,scaling,1,MPI_LCSRP,MPI_MAX,lcscomm,ierr)
		if(lcsrank==0) write(*,*) '2D scaling factor: ',scaling

		do k = 1-sgrid%ng,sgrid%nk+sgrid%ng
		do j = 1-sgrid%ng,sgrid%nj+sgrid%ng
		do i = 1-sgrid%ng,sgrid%ni+sgrid%ng
			im1 = max(i-1,1-sgrid%ng)
			ip1 = min(i+1,sgrid%ni+sgrid%ng)
			jm1 = max(j-1,1-sgrid%ng)
			jp1 = min(j+1,sgrid%nj+sgrid%ng)
			km1 = max(k-1,1-sgrid%ng)
			kp1 = min(k+1,sgrid%nk+sgrid%ng)

			delta = 0.0_LCSRP
			if(sgrid%gni==1 .and. i/=1 ) then
				v1(1) = sgrid%grid%x(1,jp1,1) - sgrid%grid%x(1,jm1,1)
				v1(2) = sgrid%grid%y(1,jp1,1) - sgrid%grid%y(1,jm1,1)
				v1(3) = sgrid%grid%z(1,jp1,1) - sgrid%grid%z(1,jm1,1)
				v2(1) = sgrid%grid%x(1,1,kp1) - sgrid%grid%x(1,1,km1)
				v2(2) = sgrid%grid%y(1,1,kp1) - sgrid%grid%y(1,1,km1)
				v2(3) = sgrid%grid%z(1,1,kp1) - sgrid%grid%z(1,1,km1)
				delta = real(i-1,LCSRP)*cross_product(v1,v2,unitvector=.true.)
			endif
			if(sgrid%gnj==1 .and. j/=1 ) then
				v1(1) = sgrid%grid%x(ip1,1,1) - sgrid%grid%x(im1,1,1)
				v1(2) = sgrid%grid%y(ip1,1,1) - sgrid%grid%y(im1,1,1)
				v1(3) = sgrid%grid%z(ip1,1,1) - sgrid%grid%z(im1,1,1)
				v2(1) = sgrid%grid%x(1,1,kp1) - sgrid%grid%x(1,1,km1)
				v2(2) = sgrid%grid%y(1,1,kp1) - sgrid%grid%y(1,1,km1)
				v2(3) = sgrid%grid%z(1,1,kp1) - sgrid%grid%z(1,1,km1)
				delta = real(j-1,LCSRP)*cross_product(v1,v2,unitvector=.true.)
			endif
			if(sgrid%gnk==1 .and. k/=1 ) then
				v1(1) = sgrid%grid%x(ip1,1,1) - sgrid%grid%x(im1,1,1)
				v1(2) = sgrid%grid%y(ip1,1,1) - sgrid%grid%y(im1,1,1)
				v1(3) = sgrid%grid%z(ip1,1,1) - sgrid%grid%z(im1,1,1)
				v2(1) = sgrid%grid%x(1,jp1,1) - sgrid%grid%x(1,jm1,1)
				v2(2) = sgrid%grid%y(1,jp1,1) - sgrid%grid%y(1,jm1,1)
				v2(3) = sgrid%grid%z(1,jp1,1) - sgrid%grid%z(1,jm1,1)
				delta = real(k-1,LCSRP)*cross_product(v1,v2,unitvector=.true.)
			endif

			!Scale delta
			delta = delta*scaling

			sgrid%grid%x(i,j,k) = sgrid%grid%x(i,j,k) + delta(1)
			sgrid%grid%y(i,j,k) = sgrid%grid%y(i,j,k) + delta(2)
			sgrid%grid%z(i,j,k) = sgrid%grid%z(i,j,k) + delta(3)
		enddo
		enddo
		enddo

		!Share the bounding box and periodic_shift with everyone:
		!TODO:  modify for twodimensions?
		allocate(sgrid%bb(0:nprocs-1,1:6))
		allocate(sgrid%ps(0:nprocs-1,-1:1,-1:1,-1:1,1:3))
		allocate(mybb(0:nprocs-1,1:6))
		allocate(myps(0:nprocs-1,-1:1,-1:1,-1:1,1:3))
		mybb = 0.0_LCSRP
		mybb(lcsrank,1) = minval(sgrid%grid%x(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		mybb(lcsrank,2) = minval(sgrid%grid%y(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		mybb(lcsrank,3) = minval(sgrid%grid%z(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		mybb(lcsrank,4) = maxval(sgrid%grid%x(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		mybb(lcsrank,5) = maxval(sgrid%grid%y(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		mybb(lcsrank,6) = maxval(sgrid%grid%z(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng))
		call MPI_ALLREDUCE(mybb(0,1),sgrid%bb(0,1),nprocs*6,MPI_LCSRP,MPI_SUM,lcscomm,ierr)
		myps = 0.0_LCSRP
		myps(lcsrank,:,:,:,:) = periodic_shift
		call MPI_ALLREDUCE(myps(0,-1,-1,-1,1),sgrid%ps(0,-1,-1,-1,1),nprocs*81,MPI_LCSRP,MPI_SUM,lcscomm,ierr)

		!Check the rectilinearity:
		call check_rectilinear(sgrid)


		!Set the outward unit normal (into the boundary).  This is constructed
		!by computing the gradient of the real valued bc flag.  Use the least
		!squares gradient with full connectivity here, to handle normal along corners,etc.
		call init_sr0(bc,ni,nj,nk,ng,'BC')
		do k=1-ng,nk+ng
		do j=1-ng,nj+ng
		do i=1-ng,ni+ng
			if(sgrid%bcflag%i(i,j,k) == LCS_2D) then
				bc%r(i,j,k) = real(LCS_INTERNAL,LCSRP)
			else
				bc%r(i,j,k) = real(min(sgrid%bcflag%i(i,j,k),1),LCSRP)
			endif
		enddo
		enddo
		enddo
		call exchange_sdata(sgrid%scomm_max_r0,r0=bc)
		call compute_lsg_wts(sgrid,full_conn=.true.)
		call grad_sr0_ls(sgrid,bc,sgrid%bcnorm)
		call destroy_lsg_wts(sgrid)
		do k=1,nk
		do j=1,nj
		do i=1,ni
			mag = sgrid%bcnorm%x(i,j,k)*sgrid%bcnorm%x(i,j,k)&
				+ sgrid%bcnorm%y(i,j,k)*sgrid%bcnorm%y(i,j,k)&
				+ sgrid%bcnorm%z(i,j,k)*sgrid%bcnorm%z(i,j,k)
			if(mag<=0.0_LCSRP) then
				sgrid%bcnorm%x(i,j,k) = 0.0_LCSRP
				sgrid%bcnorm%y(i,j,k) = 0.0_LCSRP
				sgrid%bcnorm%z(i,j,k) = 0.0_LCSRP
			else
				sgrid%bcnorm%x(i,j,k) = sgrid%bcnorm%x(i,j,k) / sqrt(mag)
				sgrid%bcnorm%y(i,j,k) = sgrid%bcnorm%y(i,j,k) / sqrt(mag)
				sgrid%bcnorm%z(i,j,k) = sgrid%bcnorm%z(i,j,k) / sqrt(mag)
			endif
		enddo
		enddo
		enddo
		call destroy_sr0(bc)

		!Compute the least squares gradient weights, if you need them
		if(.NOT. sgrid%rectilinear) then
			call compute_lsg_wts(sgrid,full_conn=.false.)
		endif

		!For debugging, write the structured grid to a file:
		if(DEBUG_SGRID) then
			gn = (/sgrid%gni,sgrid%gnj,sgrid%gnk/)
			offset = (/sgrid%offset_i,sgrid%offset_j,sgrid%offset_k/)
			write(fname, '(a,a,a,a,a)') './',trim(TEMP_DIR),'/',trim(sgrid%label),trim(FILE_EXT)
			call structured_io(trim(fname),IO_WRITE,gn,offset,r1=sgrid%grid)	!write the grid
			call structured_io(trim(fname),IO_APPEND,gn,offset,r1=sgrid%bcnorm)	!write the gridnorm
			!Append the flag:
			call init_sr0(tmp,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,'FLAG')
			tmp%r  = real(sgrid%bcflag%i)
			call structured_io(trim(fname),IO_APPEND,gn,offset,r0=tmp)
			call destroy_sr0(tmp)
		endif

		!cleanup
		call destroy_sr1(fake)

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
		real(LCSRP):: dx(3)
		integer:: restest,ierr
		real(LCSRP),allocatable:: x(:,:,:)
		real(LCSRP),allocatable:: y(:,:,:)
		real(LCSRP),allocatable:: z(:,:,:)
		integer(LCSIP),allocatable:: flag(:,:,:)
		real(LCSRP):: v0(3),v1(3),v2(3),v3(3),v4(3),v5(3),v6(3),v7(3)
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
						if(.NOT. sgrid%periodic_i .AND. i == sgrid%gni) cycle
						gn_new(1) =gn_new(1)+1
					enddo
					if(i > sgrid%offset_i .AND. i <= sgrid%offset_i + sgrid%ni) then
						n_new(1) = n_new(1) + 1
						do ii = 1,res
							if(.NOT. sgrid%periodic_i .AND. i == sgrid%gni) cycle
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
						if(.NOT. sgrid%periodic_j .AND. j == sgrid%gnj) cycle
						gn_new(2) =gn_new(2)+1
					enddo
					if(j > sgrid%offset_j .AND. j <= sgrid%offset_j + sgrid%nj) then
						n_new(2) = n_new(2) + 1
						do jj = 1,res
							if(.NOT. sgrid%periodic_j .AND. j == sgrid%gnj) cycle
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
						if(.NOT. sgrid%periodic_k .AND. k == sgrid%gnk) cycle
						gn_new(3) =gn_new(3)+1
					enddo
					if(k> sgrid%offset_k .AND. k <= sgrid%offset_k + sgrid%nk) then
						n_new(3) = n_new(3) + 1
						do kk = 1,res
							if(.NOT. sgrid%periodic_k .AND. k == sgrid%gnk) cycle
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
		allocate(flag(1:n_new(1),1:n_new(2),1:n_new(3)))

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
						i_new = ii + (ii-1)*res +ir;
						j_new = jj + (jj-1)*res +jr;
						k_new = kk + (kk-1)*res +kr;
						!New grid points
						if (i_new <=n_new(1) .and. j_new <= n_new(2) .and. k_new <= n_new(3)) then
							!Spacing between the new grid points in each dir:  Assumes ng >= 1
							if(sgrid%rectilinear)then
								dx(1) = (sgrid%grid%x(ii+1,jj,kk) - sgrid%grid%x(ii,jj,kk)) / real(res +1,LCSRP)
								dx(2) = (sgrid%grid%y(ii,jj+1,kk) - sgrid%grid%y(ii,jj,kk)) / real(res +1,LCSRP)
								dx(3) = (sgrid%grid%z(ii,jj,kk+1) - sgrid%grid%z(ii,jj,kk)) / real(res +1,LCSRP)
								x(i_new,j_new,k_new) = sgrid%grid%x(ii,jj,kk) + dx(1)*real(ir,LCSRP)
								y(i_new,j_new,k_new) = sgrid%grid%y(ii,jj,kk) + dx(2)*real(jr,LCSRP)
								z(i_new,j_new,k_new) = sgrid%grid%z(ii,jj,kk) + dx(3)*real(kr,LCSRP)
							else
								v0(1)=sgrid%grid%x(ii,jj,kk);v0(2)=sgrid%grid%y(ii,jj,kk);v0(3)=sgrid%grid%z(ii,jj,kk);
								v1(1)=sgrid%grid%x(ii+1,jj,kk);v1(2)=sgrid%grid%y(ii+1,jj,kk);v1(3)=sgrid%grid%z(ii+1,jj,kk);
								v2(1)=sgrid%grid%x(ii,jj+1,kk);v2(2)=sgrid%grid%y(ii,jj+1,kk);v2(3)=sgrid%grid%z(ii,jj+1,kk);
								v3(1)=sgrid%grid%x(ii+1,jj+1,kk);v3(2)=sgrid%grid%y(ii+1,jj+1,kk);v3(3)=sgrid%grid%z(ii+1,jj+1,kk);
								v4(1)=sgrid%grid%x(ii,jj,kk+1);v4(2)=sgrid%grid%y(ii,jj,kk+1);v4(3)=sgrid%grid%z(ii,jj,kk+1);
								v5(1)=sgrid%grid%x(ii+1,jj,kk+1);v5(2)=sgrid%grid%y(ii+1,jj,kk+1);v5(3)=sgrid%grid%z(ii+1,jj,kk+1);
								v6(1)=sgrid%grid%x(ii,jj+1,kk+1);v6(2)=sgrid%grid%y(ii,jj+1,kk+1);v6(3)=sgrid%grid%z(ii,jj+1,kk+1);
								v7(1)=sgrid%grid%x(ii+1,jj+1,kk+1);v7(2)=sgrid%grid%y(ii+1,jj+1,kk+1);v7(3)=sgrid%grid%z(ii+1,jj+1,kk+1);

								dx = grid_refine(v0,v1,v2,v3,v4,v5,v6,v7,real(ir,LCSRP)/real(res+1,LCSRP)&
									,real(jr,LCSRP)/real(res+1,LCSRP),real(kr,LCSRP)/real(res+1,LCSRP))
								x(i_new,j_new,k_new) =  dx(1)
								y(i_new,j_new,k_new) =  dx(2)
								z(i_new,j_new,k_new) =  dx(3)
							endif
							!inherit flag
							flag(i_new,j_new,k_new) = sgrid%bcflag%i(ii,jj,kk)
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
						!inherit flag
						flag(i_new,j_new,k_new) = sgrid%bcflag%i(ii,jj,kk)
					enddo
				enddo
			enddo
		endif

		!-----
		!Copy the other properties and initialize another sgrid structure
		!-----
		call init_sgrid(sgrid_new,label,n_new,offset_new,x,y,z,flag)

		!-----
		!Deallocate/cleanup
		!-----
		deallocate(x)
		deallocate(y)
		deallocate(z)
		deallocate(flag)

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
		sgrid%periodic_i = .false.
		sgrid%periodic_j = .false.
		sgrid%periodic_k = .false.
		sgrid%rectilinear = .false.

		if(allocated(sgrid%bb))deallocate(sgrid%bb)
		if(allocated(sgrid%ps))deallocate(sgrid%ps)

		call destroy_scomm(sgrid%scomm_face_r0)
		call destroy_scomm(sgrid%scomm_max_r0)
		call destroy_scomm(sgrid%scomm_face_r1)
		call destroy_scomm(sgrid%scomm_max_r1)
		call destroy_scomm(sgrid%scomm_face_r2)
		call destroy_scomm(sgrid%scomm_max_r2)

		call destroy_sr1(sgrid%grid)
		call destroy_si0(sgrid%bcflag)
		call destroy_sr1(sgrid%bcnorm)

		call destroy_lsg_wts(sgrid)

		sgrid%label = 'Unused sgrid'

	end subroutine destroy_sgrid

	subroutine set_velocity_bc(sgrid,vel)
		implicit none
		!-----
		type(sgrid_t):: sgrid
		type(sr1_t):: vel
		!-----
		integer:: i,j,k,ni,nj,nk,ng
		integer:: i_b,j_b,k_b
		!-----
		!Set the velocity BC for the ghost/fake nodes
		!Note, we assume the user's velocity field
		!Respects the desired boundary condition AT the IB nodes
		!-----

		if(lcsrank==0 .and. LCS_VERBOSE) &
			write(*,*) 'In set_velocity_bc... '

		ni = vel%ni
		nj = vel%nj
		nk = vel%nk
		ng = vel%ng

		do k = 1-ng,nk+ng
		do j = 1-ng,nj+ng
		do i = 1-ng,ni+ng
			select case(sgrid%bcflag%i(i,j,k))
			case(LCS_INTERNAL)
				cycle
			case(LCS_MASK)
				!Always zero velocity in a masked node
				vel%x(i,j,k) = 0.0
				vel%y(i,j,k) = 0.0
				vel%z(i,j,k) = 0.0
			case(LCS_WALL)
				!Make sure any velocity is zeroed in any ghost/fake that is a wall:
				if(i>=1 .and. j>=1 .and. k>=1 .and. i<=ni .and. j<=nj .and. k<=nk) cycle
				vel%x(i,j,k) = 0.0
				vel%y(i,j,k) = 0.0
				vel%z(i,j,k) = 0.0
			CASE(LCS_SLIP,LCS_INFLOW,LCS_OUTFLOW,LCS_2D)
				!For these, just set zero gradient.
				if(i>=1 .and. j>=1 .and. k>=1 .and. i<=ni .and. j<=nj .and. k<=nk) cycle
				i_b = max(min(i,ni),1)
				j_b = max(min(j,nj),1)
				k_b = max(min(k,nk),1)
				vel%x(i,j,k) = vel%x(i_b,j_b,k_b)
				vel%y(i,j,k) = vel%y(i_b,j_b,k_b)
				vel%z(i,j,k) = vel%z(i_b,j_b,k_b)
			case default
				write(*,*) 'lcsrank[',lcsrank,'] ERROR: Unknown bcflag:',sgrid%bcflag%i(i,j,k)
				CFD2LCS_ERROR=1
			end select
		enddo
		enddo
		enddo

	end subroutine set_velocity_bc


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
