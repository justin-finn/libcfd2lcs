module comms_m
	use mpi_m
	use data_m
	implicit none
	!-----

	!Checkerboard pattern:
	integer(LCSIP),parameter:: &
		RED = 0, &
		BLACK = 1

	!Communication flag:
	integer(LCSIP),parameter:: &
		NBR_COMM = 0, &
		SELF_COMM = 1,&
		NO_COMM = 2

	!Connectivity type:
	integer(LCSIP),parameter:: &
		FACE_CONNECT = 0 ,&
		MAX_CONNECT = 1

	!Starting ID for point-to-point comms
	!make larger than nproc
	integer(LCSIP), parameter:: TAG_START = 123456

	!Flag for datatype
	integer(LCSIP), parameter:: &
		R0_COMM = 1, &
		R1_COMM = 3, &
		R2_COMM = 9

	!i,j,k Cartesian offsets for each communication vector
	integer(LCSIP),parameter:: COMM_OFFSET(3,26) = reshape((/&
		  -1,   0,   0, &  !face
		   0,  -1,   0, &  !face
		   0,   0,  -1, &  !face
		   1,   0,   0, &  !face
		   0,   1,   0, &  !face
		   0,   0,   1, &  !face
		  -1,  -1,  -1, &
		   0,  -1,  -1, &
		   1,  -1,  -1, &
		  -1,   0,  -1, &
		   1,   0,  -1, &
		  -1,   1,  -1, &
		   0,   1,  -1, &
		   1,   1,  -1, &
		  -1,  -1,   0, &
		   1,  -1,   0, &
		  -1,   1,   0, &
		   1,   1,   0, &
		  -1,  -1,   1, &
		   0,  -1,   1, &
		   1,  -1,   1, &
		  -1,   0,   1, &
		   1,   0,   1, &
		  -1,   1,   1, &
		   0,   1,   1, &
		   1,   1,   1 &
		/),(/3,26/))

	contains

	subroutine init_scomm(scomm,ni,nj,nk,ng,offset_i,offset_j,offset_k,bc_list,lperiodic,connectivity,datatype,label)
		implicit none
		!-----
		type(scomm_t):: scomm
		integer:: ni,nj,nk,ng
		integer:: offset_i,offset_j,offset_k
		integer(LCSIP),dimension(6):: bc_list
		real(LCSIP):: lperiodic(3)
		integer(LCSIP):: connectivity
		integer(LCSIP):: datatype
		character(len=*) label
		!-----
		integer:: imin,imax,jmin,jmax,kmin,kmax
		integer:: rank_i,rank_j,rank_k
		integer:: npi,npj,npk
		!-----

		!
		! Make sure we are starting from clean data sturcture:
		!
		call destroy_scomm(scomm)

		scomm%label = trim(label)
		if(lcsrank==0) &
			write(*,*) 'in init_scomm... ',trim(scomm%label)

		!
		! Set the connectivity,datatype
		!
		if(connectivity == MAX_CONNECT .or. connectivity == FACE_CONNECT ) then
			scomm%connectivity = connectivity !ok...
		else
			write(*,*) 'in init_scomm... ERROR:  connectivity must be, "FACE_CONNECT" or "MAX_CONNECT"'
			CFD2LCS_ERROR = 1
			return
		endif
		if(datatype == R0_COMM  .or. datatype == R1_COMM .or. datatype == R2_COMM ) then
			scomm%datatype = datatype !ok...
		else
			write(*,*) 'in init_scomm... ERROR:  datatype must be "R0_COMM", "R1_COMM" or "R2_COMM"'
			CFD2LCS_ERROR = 1
			return
		endif

		!
		!Check periodicity:
		!
		if (bc_list(1) == LCS_PERIODIC .AND. bc_list(4) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
		if (bc_list(4) == LCS_PERIODIC .AND. bc_list(1) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
		if (bc_list(2) == LCS_PERIODIC .AND. bc_list(5) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
		if (bc_list(5) == LCS_PERIODIC .AND. bc_list(2) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
		if (bc_list(6) == LCS_PERIODIC .AND. bc_list(6) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
		if (bc_list(3) == LCS_PERIODIC .AND. bc_list(3) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
		if(CFD2LCS_ERROR ==1)then
			if(lcsrank==0)&
				write(*,*) 'ERROR: Periodic BCs do not match'
			return
		endif



		!
		! Establish the number of processors in each direction, and each processor's i,j,k rank:
		!
		call establish_nbrs()

		!
		! Set the checkerboard pattern
		!
		call set_checkerboard()

		!
		! Set the pack/unpack ranges
		!
		call set_pack_unpack_list()

		!
		! Test the comm pattern:
		!
		call handshake()


		contains


		subroutine establish_nbrs()
			implicit none
			!-----
			integer,allocatable:: tmp(:),i0(:),j0(:),k0(:),i1(:),j1(:),k1(:)
			integer:: i,j,k,ierr,rank
			integer:: im1, ip1, jm1, jp1, km1, kp1
			integer:: isearch,jsearch,ksearch
			integer:: found
			integer:: low_i,low_j,low_k
			!-----

			!
			!Everybody shares their local indices:
			!
			allocate(i0(0:nprocs-1))
			allocate(j0(0:nprocs-1))
			allocate(k0(0:nprocs-1))
			allocate(i1(0:nprocs-1))
			allocate(j1(0:nprocs-1))
			allocate(k1(0:nprocs-1))
			allocate(tmp(0:nprocs-1))
			tmp = -10000; tmp(lcsrank) = (offset_i+1)
			call MPI_ALLREDUCE(tmp,i0,nprocs,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
			tmp = -10000; tmp(lcsrank) = (offset_j+1)
			call MPI_ALLREDUCE(tmp,j0,nprocs,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
			tmp = -10000; tmp(lcsrank) = (offset_k+1)
			call MPI_ALLREDUCE(tmp,k0,nprocs,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
			tmp = -10000; tmp(lcsrank) = (offset_i+ni)
			call MPI_ALLREDUCE(tmp,i1,nprocs,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
			tmp = -10000; tmp(lcsrank) = (offset_j+nj)
			call MPI_ALLREDUCE(tmp,j1,nprocs,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
			tmp = -10000; tmp(lcsrank) = (offset_k+nk)
			call MPI_ALLREDUCE(tmp,k1,nprocs,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
			!Global max/min of grid index
			call MPI_ALLREDUCE(minval(i0),imin,1,MPI_INTEGER,MPI_MIN,lcscomm,ierr)
			call MPI_ALLREDUCE(minval(j0),jmin,1,MPI_INTEGER,MPI_MIN,lcscomm,ierr)
			call MPI_ALLREDUCE(minval(k0),kmin,1,MPI_INTEGER,MPI_MIN,lcscomm,ierr)
			call MPI_ALLREDUCE(maxval(i1),imax,1,MPI_INTEGER,MPI_MIN,lcscomm,ierr)
			call MPI_ALLREDUCE(maxval(j1),jmax,1,MPI_INTEGER,MPI_MIN,lcscomm,ierr)
			call MPI_ALLREDUCE(maxval(k1),kmax,1,MPI_INTEGER,MPI_MIN,lcscomm,ierr)

			!
			!Establish rank_i, rank_j, rank_k
			!Check that the number of procs in each direction is even (or 1)
			!
			rank_i = 1; rank_j =1; rank_k = 1
			low_i = -1; low_j = -1; low_k = -1
			do rank = 0,nprocs-1
				if(i0(rank) < i0(lcsrank) .and. i0(rank) > low_i) then
					rank_i = rank_i + 1
					low_i = i0(rank)
				endif
				if(j0(rank) < j0(lcsrank) .and. j0(rank) > low_j) then
					rank_j = rank_j + 1
					low_j = j0(rank)
				endif
				if(k0(rank) < k0(lcsrank) .and. k0(rank) > low_k) then
					rank_k = rank_k + 1
					low_k = k0(rank)
				endif
			enddo
			call MPI_ALLREDUCE(rank_i,npi,1,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
			call MPI_ALLREDUCE(rank_j,npj,1,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
			call MPI_ALLREDUCE(rank_k,npk,1,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
			if(lcsrank==0)&
				write(*,*) 'Num Procs in I,J,K=', npi,npj,npk
			if(mod(npi,2)/=0 .AND. npi/=1) then
				if(lcsrank==0)&
					write(*,*) 'ERROR:', 'Number of processors in "X" direction must be even (or 1)'
					CFD2LCS_ERROR = 1
					return
			endif
			if(mod(npj,2)/=0 .AND. npj/=1) then
				if(lcsrank==0)&
					write(*,*) 'ERROR:', 'Number of processors in "Y" direction must be even (or 1)'
					CFD2LCS_ERROR = 1
					return
			endif
			if(mod(npk,2)/=0 .AND. npk/=1) then
				if(lcsrank==0)&
					write(*,*) 'ERROR:', 'Number of processors in "Z" direction must be even (or 1)'
					CFD2LCS_ERROR = 1
					return
			endif

			!
			!Search for the neighbors in each direction, handling the
			!possibility of global periodicity:
			!
			im1 = i0(lcsrank) - 1
			ip1 = i1(lcsrank) + 1
			jm1 = j0(lcsrank) - 1
			jp1 = j1(lcsrank) + 1
			km1 = k0(lcsrank) - 1
			kp1 = k1(lcsrank) + 1
			!Handle potential periodicity...
			if(i0(lcsrank) == imin)im1	= imax
			if(i1(lcsrank) == imax)ip1	= imin
			if(j0(lcsrank) == jmin)jm1	= jmax
			if(j1(lcsrank) == jmax)jp1	= jmin
			if(k0(lcsrank) == kmin)km1	= kmax
			if(k1(lcsrank) == kmax)kp1	= kmin
			do k=-1,1
			do j=-1,1
			do i=-1,1
				if(i==-1)	isearch = im1
				if(i==0)  	isearch = i0(lcsrank)
				if(i==1)	isearch = ip1
				if(j==-1)	jsearch = jm1
				if(j==0)  	jsearch = j0(lcsrank)
				if(j==1)	jsearch = jp1
				if(k==-1)	ksearch = km1
				if(k==0)  	ksearch = k0(lcsrank)
				if(k==1)	ksearch = kp1
				found = 0
				do rank = 0,nprocs-1
					if( isearch >= i0(rank) .and. isearch <= i1(rank) .and. &
						jsearch >= j0(rank) .and. jsearch <= j1(rank) .and. &
						ksearch >= k0(rank) .and. ksearch <= k1(rank) ) then
						found = found + 1
						scomm%nbr_rank(i,j,k) = rank
					endif
				enddo

				!Handle the periodic_shift:
				if(i==-1 .AND. bc_list(1) == LCS_PERIODIC .AND. i0(lcsrank) == imin) &
					scomm%periodic_shift(i,j,k,1) = lperiodic(1)
				if(i==1 .AND. bc_list(4) == LCS_PERIODIC .AND. i1(lcsrank) == imax) &
					scomm%periodic_shift(i,j,k,1) = -lperiodic(1)
				if(j==-1 .AND. bc_list(2) == LCS_PERIODIC .AND. j0(lcsrank) == jmin) &
					scomm%periodic_shift(i,j,k,2) = lperiodic(2)
				if(j==1 .AND. bc_list(5) == LCS_PERIODIC .AND. j1(lcsrank) == jmax) &
					scomm%periodic_shift(i,j,k,2) = -lperiodic(2)
				if(k==-1 .AND. bc_list(3) == LCS_PERIODIC .AND. k0(lcsrank) == kmin) &
					scomm%periodic_shift(i,j,k,3) = lperiodic(3)
				if(k==1 .AND. bc_list(6) == LCS_PERIODIC .AND. k1(lcsrank) == kmax) &
					scomm%periodic_shift(i,j,k,3) = -lperiodic(3)

				if (found /= 1) then
					write(*,*) 'Error:  lcsrank[',lcsrank,&
						'] cant establish unique communication neighbors in direction',i,j,k
					CFD2LCS_ERROR = 1
				endif
			enddo
			enddo
			enddo

			deallocate(i0)
			deallocate(j0)
			deallocate(k0)
			deallocate(i1)
			deallocate(j1)
			deallocate(k1)
			deallocate(tmp)

		end subroutine establish_nbrs

		subroutine set_checkerboard()
			implicit none
			!-----
			integer:: i,j,k
			integer:: checker_i,checker_j,checker_k
			!-----
			!Set the checkerboard pattern for each point-to-point communication
			!There are three possible patterns to use (stripes in each direction),
			!depending on the spatial communication vector.
			!-----

			if(mod(rank_i,2)==0) then
				checker_i = RED
			else
				checker_i = BLACK
			endif
			if(mod(rank_j,2)==0) then
				checker_j = RED
			else
				checker_j = BLACK
			endif
			if(mod(rank_k,2)==0) then
				checker_k = RED
			else
				checker_k = BLACK
			endif


			!Set the checkerboard for this i,j,k vector
			do k = -1,1
			do j = -1,1
			do i = -1,1
				!Choose a checker pattern based on comm vector
				if(i==0 .OR. npi ==1) then
					if(j==0 .OR. npj==1) then
						scomm%checker(i,j,k) = checker_k
					else
						scomm%checker(i,j,k) = checker_j
					endif
				else
					scomm%checker(i,j,k) = checker_i
				endif
			enddo
			enddo
			enddo

			!Make sure periodicity matches:
			if(bc_list(1) == LCS_PERIODIC .AND. bc_list(4) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
			if(bc_list(4) == LCS_PERIODIC .AND. bc_list(1) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
			if(bc_list(2) == LCS_PERIODIC .AND. bc_list(5) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
			if(bc_list(5) == LCS_PERIODIC .AND. bc_list(2) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
			if(bc_list(3) == LCS_PERIODIC .AND. bc_list(6) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
			if(bc_list(6) == LCS_PERIODIC .AND. bc_list(3) /= LCS_PERIODIC) CFD2LCS_ERROR = 1
			if(CFD2LCS_ERROR ==1) then
				write(*,*) 'ERROR: Periodicity does not match in bc_list'
			endif

			!
			! Set the flag to either:
			!  NBR_COMM:  point to point with another proc
			!  SELF_COMM:  exchange buffers with self
			!  NO_COMM:  no comm needed in this direction (because of BC)
			!
			scomm%flag = NBR_COMM !(default)

			do k = -1,1
			do j = -1,1
			do i = -1,1
				if (scomm%nbr_rank(i,j,k) == lcsrank)then
					scomm%flag(i,j,k) = SELF_COMM
				endif
			enddo
			enddo
			enddo

			if (rank_i == 1 .AND. bc_list(1)/=LCS_PERIODIC )then
				scomm%flag(-1,:,:) = NO_COMM
			endif
			if (rank_j == 1 .AND. bc_list(2)/=LCS_PERIODIC )then
				scomm%flag(:,-1,:) = NO_COMM
			endif
			if (rank_k == 1 .AND. bc_list(3)/=LCS_PERIODIC )then
				scomm%flag(:,:,-1) = NO_COMM
			endif
			if (rank_i == npi .AND. bc_list(4)/=LCS_PERIODIC )then
				scomm%flag(1,:,:) = NO_COMM
			endif
			if (rank_j == npj .AND. bc_list(5)/=LCS_PERIODIC )then
				scomm%flag(:,1,:) = NO_COMM
			endif
			if (rank_k == npk .AND. bc_list(6)/=LCS_PERIODIC )then
				scomm%flag(:,:,1) = NO_COMM
			endif

		end subroutine set_checkerboard

		subroutine set_pack_unpack_list()
			implicit none
			!-----
			integer:: i,j,k,icomm,ncomm
			!-----
			! Set the ranges for grid data exchange
			! Allocate buffers
			!-----

			if(ng < 0) then
				write(*,*) 'ERROR:  ng = 0'
				CFD2LCS_ERROR = 1
				return
			endif

			do k = -1,1
			do j = -1,1
			do i = -1,1
				select case(i)
					case(-1)
					scomm%pack_list_min(i,j,k,1) = 1
					scomm%pack_list_max(i,j,k,1) = ng
					scomm%unpack_list_min(i,j,k,1) = 1-ng
					scomm%unpack_list_max(i,j,k,1) = 0
					case(0)
					scomm%pack_list_min(i,j,k,1) = 1
					scomm%pack_list_max(i,j,k,1) = ni
					scomm%unpack_list_min(i,j,k,1) = 1
					scomm%unpack_list_max(i,j,k,1) = ni
					case(1)
					scomm%pack_list_min(i,j,k,1) = ni-ng+1
					scomm%pack_list_max(i,j,k,1) = ni
					scomm%unpack_list_min(i,j,k,1) = ni+1
					scomm%unpack_list_max(i,j,k,1) = ni+ng
				end select

				select case(j)
					case(-1)
					scomm%pack_list_min(i,j,k,2) = 1
					scomm%pack_list_max(i,j,k,2) = ng
					scomm%unpack_list_min(i,j,k,2) = 1-ng
					scomm%unpack_list_max(i,j,k,2) = 0
					case(0)
					scomm%pack_list_min(i,j,k,2) = 1
					scomm%pack_list_max(i,j,k,2) = nj
					scomm%unpack_list_min(i,j,k,2) = 1
					scomm%unpack_list_max(i,j,k,2) = nj
					case(1)
					scomm%pack_list_min(i,j,k,2) = nj-ng+1
					scomm%pack_list_max(i,j,k,2) = nj
					scomm%unpack_list_min(i,j,k,2) = nj+1
					scomm%unpack_list_max(i,j,k,2) = nj+ng
				end select

				select case(k)
					case(-1)
					scomm%pack_list_min(i,j,k,3) = 1
					scomm%pack_list_max(i,j,k,3) = ng
					scomm%unpack_list_min(i,j,k,3) = 1-ng
					scomm%unpack_list_max(i,j,k,3) = 0
					case(0)
					scomm%pack_list_min(i,j,k,3) = 1
					scomm%pack_list_max(i,j,k,3) = nk
					scomm%unpack_list_min(i,j,k,3) = 1
					scomm%unpack_list_max(i,j,k,3) = nk
					case(1)
					scomm%pack_list_min(i,j,k,3) = nk-ng+1
					scomm%pack_list_max(i,j,k,3) = nk
					scomm%unpack_list_min(i,j,k,3) = nk+1
					scomm%unpack_list_max(i,j,k,3) = nk+ng
				end select
			enddo
			enddo
			enddo


			!checks...?

			if(scomm%connectivity == FACE_CONNECT) ncomm = 6
			if(scomm%connectivity == MAX_CONNECT) ncomm = 26
			scomm%n_pack = 0
			scomm%n_unpack = 0
			scomm%pack_bufsize = 0
			scomm%unpack_bufsize = 0
			scomm%pack_start = 0
			scomm%unpack_start = 0
			do icomm = 1,ncomm
				i = COMM_OFFSET(1,icomm)
				j = COMM_OFFSET(2,icomm)
				k = COMM_OFFSET(3,icomm)
				if(scomm%flag(i,j,k) == NO_COMM) cycle

				scomm%n_pack(i,j,k) = (scomm%pack_list_max(i,j,k,1)-scomm%pack_list_min(i,j,k,1)+1) &
									* (scomm%pack_list_max(i,j,k,2)-scomm%pack_list_min(i,j,k,2)+1) &
									* (scomm%pack_list_max(i,j,k,3)-scomm%pack_list_min(i,j,k,3)+1) &
									* scomm%datatype
				scomm%n_unpack(i,j,k) = (scomm%unpack_list_max(i,j,k,1)-scomm%unpack_list_min(i,j,k,1)+1) &
									*   (scomm%unpack_list_max(i,j,k,2)-scomm%unpack_list_min(i,j,k,2)+1) &
									*   (scomm%unpack_list_max(i,j,k,3)-scomm%unpack_list_min(i,j,k,3)+1) &
									* scomm%datatype

				scomm%pack_bufsize = scomm%pack_bufsize + scomm%n_pack(i,j,k)
				scomm%unpack_bufsize = scomm%unpack_bufsize + scomm%n_unpack(i,j,k)

				scomm%pack_start(i,j,k)   = scomm%pack_bufsize   - scomm%n_pack(i,j,k) + 1
				scomm%unpack_start(i,j,k) = scomm%unpack_bufsize - scomm%n_unpack(i,j,k) + 1
			enddo

		if(PREALLOCATE_BUFFERS) then
			if(lcsrank==0) then
				write(*,*) 'Preallocating Pack Buffer, Size (MB):', real(scomm%pack_bufsize * LCSRP)/1000000.0
				write(*,*) 'Preallocating Unpack Buffer, Size (MB):', real(scomm%unpack_bufsize * LCSRP)/1000000.0
			endif
			allocate(scomm%pack_buffer(scomm%pack_bufsize))
			allocate(scomm%unpack_buffer(scomm%unpack_bufsize))
		endif

		end subroutine set_pack_unpack_list


		subroutine handshake()
			implicit none
			!-----
			integer:: i,j,k,icomm,ncomm
			integer:: ibuf,ii,jj,kk
			integer:: ierr, status(MPI_STATUS_SIZE)
			integer:: tag_red,tag_black,comm_id
			integer:: nsend,nrecv,send_count,recv_count
			integer:: ihash,jhash,khash,hashval,bufval
			integer:: ipack,iunpack
			!-----
			! A basic check to make sure that all point-to-point communications
			! are working as expected.  This is also a skeleton for all other
			! new comm routines.
			!-----

			if(lcsrank==0) &
				write(*,*) 'in handshake...'

			if(.NOT. PREALLOCATE_BUFFERS) then
				allocate(scomm%pack_buffer(scomm%pack_bufsize))
				allocate(scomm%unpack_buffer(scomm%unpack_bufsize))
			endif

			!
			! Fill the pack buffer
			! Test the send/recv by passing a global hash for each grid i,j,k.
			!
			scomm%unpack_buffer=-1
			scomm%pack_buffer=-1
			if(scomm%connectivity == FACE_CONNECT) ncomm = 6
			if(scomm%connectivity == MAX_CONNECT) ncomm = 26
			do icomm = 1,ncomm

				i = COMM_OFFSET(1,icomm)
				j = COMM_OFFSET(2,icomm)
				k = COMM_OFFSET(3,icomm)

				if(scomm%flag(i,j,k) == NO_COMM) cycle

				ibuf = scomm%pack_start(i,j,k)
				do kk = scomm%pack_list_min(i,j,k,3),scomm%pack_list_max(i,j,k,3)
				do jj = scomm%pack_list_min(i,j,k,2),scomm%pack_list_max(i,j,k,2)
				do ii = scomm%pack_list_min(i,j,k,1),scomm%pack_list_max(i,j,k,1)

					!*********
					!Here, substitute the hash for real data
					ihash = ii + offset_i
					jhash = jj + offset_j
					khash = kk + offset_k
					hashval = HASH(ihash,jhash,khash,imax,jmax)  !Send the global i,j,k coords.
					!*********

					scomm%pack_buffer(ibuf) = real(hashval,LCSRP)

					ibuf = ibuf + 1
				enddo
				enddo
				enddo
			enddo

			!
			! Send/Recv
			!
			comm_id = 0
			send_count = 0
			recv_count = 0
			if(scomm%connectivity == FACE_CONNECT) ncomm = 6
			if(scomm%connectivity == MAX_CONNECT) ncomm = 26
			do icomm = 1,ncomm
				i = COMM_OFFSET(1,icomm)
				j = COMM_OFFSET(2,icomm)
				k = COMM_OFFSET(3,icomm)

				!Index the tags:
				comm_id = comm_id +1
				tag_red = TAG_START + comm_id
				tag_black = 2*TAG_START + comm_id

				!Red sends first then receives:
				if (scomm%checker(i,j,k) == RED) then
					if(scomm%flag(i,j,k) == NBR_COMM) then
						call MPI_SEND(scomm%pack_buffer(scomm%pack_start(i,j,k)),scomm%n_pack(i,j,k),&
							MPI_LCSRP,scomm%nbr_rank(i,j,k),TAG_RED,lcscomm,ierr)
						send_count = send_count + 1
					endif
					if (scomm%flag(-i,-j,-k) == NBR_COMM) then
						call MPI_RECV(scomm%unpack_buffer(scomm%unpack_start(-i,-j,-k)),scomm%n_unpack(-i,-j,-k),&
							MPI_LCSRP,scomm%nbr_rank(-i,-j,-k),TAG_BLACK,lcscomm,status,ierr)
						recv_count = recv_count + 1
					endif
				endif

				!Black receives first then sends:
				if (scomm%checker(i,j,k) == BLACK) then
					if (scomm%flag(-i,-j,-k) == NBR_COMM) then
						call MPI_RECV(scomm%unpack_buffer(scomm%unpack_start(-i,-j,-k)),scomm%n_unpack(-i,-j,-k),&
							MPI_LCSRP,scomm%nbr_rank(-i,-j,-k),TAG_RED,lcscomm,status,ierr)
						recv_count = recv_count + 1
					endif
					if(scomm%flag(i,j,k) == NBR_COMM) then
						call MPI_SEND(scomm%pack_buffer(scomm%pack_start(i,j,k)),scomm%n_pack(i,j,k),&
							MPI_LCSRP,scomm%nbr_rank(i,j,k),TAG_BLACK,lcscomm,ierr)
						send_count = send_count + 1
					endif
				endif

				!Self communication:
				if (scomm%flag(i,j,k) == SELF_COMM) then
					ipack = scomm%pack_start(i,j,k)
					iunpack = scomm%unpack_start(-i,-j,-k)
					do  ibuf = 1,scomm%n_unpack(-i,-j,-k)
						scomm%unpack_buffer(iunpack) = scomm%pack_buffer(ipack)
						ipack = ipack + 1
						iunpack = iunpack + 1
					enddo
					recv_count = recv_count + 1
					send_count = send_count + 1
				endif
			enddo

			!
			! Unpack/Check the buffers:
			!
			do icomm = 1,ncomm
				i = COMM_OFFSET(1,icomm)
				j = COMM_OFFSET(2,icomm)
				k = COMM_OFFSET(3,icomm)
				if(scomm%flag(i,j,k) == NO_COMM) cycle

				ibuf = scomm%unpack_start(i,j,k)
				do kk = scomm%unpack_list_min(i,j,k,3),scomm%unpack_list_max(i,j,k,3)
				do jj = scomm%unpack_list_min(i,j,k,2),scomm%unpack_list_max(i,j,k,2)
				do ii = scomm%unpack_list_min(i,j,k,1),scomm%unpack_list_max(i,j,k,1)

					!************
					!Check the global hash key for this ii,jj,kk (account for global periodicity)
					ihash = ii + offset_i
					jhash = jj + offset_j
					khash = kk + offset_k
					if(ihash < imin) ihash = imax-(imin-ihash-1)
					if(jhash < jmin) jhash = jmax-(jmin-jhash-1)
					if(khash < kmin) khash = kmax-(kmin-khash-1)
					if(ihash > imax) ihash = imin+(ihash-imax-1)
					if(jhash > jmax) jhash = jmin+(jhash-jmax-1)
					if(khash > kmax) khash = kmin+(khash-kmax-1)
					hashval = HASH(ihash,jhash,khash,imax,jmax)
					bufval = int(scomm%unpack_buffer(ibuf))
					!************

					if(hashval /= bufval) then
						write(*,*) 'ERROR: myrank[',lcsrank,'] Has unexpected value in unpack buffer'
						CFD2LCS_ERROR = 1
					endif

					ibuf = ibuf + 1
				enddo
				enddo
				enddo
			enddo

			!Check the number of comms:
			call MPI_REDUCE(send_count,nsend,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)
			call MPI_REDUCE(recv_count,nrecv,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)
			if(lcsrank==0) then
				if(nsend == nrecv)then
					write(*,*) '     handshake OK for ',nsend, 'communications'
				else
					write(*,*) '    ERROR: nsend /= nrecv :',nsend,nrecv
					CFD2LCS_ERROR = 1
				endif
			endif

			if(.NOT. PREALLOCATE_BUFFERS) then
				deallocate(scomm%pack_buffer)
				deallocate(scomm%unpack_buffer)
			endif

		end subroutine handshake

		integer function HASH(i,j,k,ni,nj)
			implicit none
			integer::i,j,k,ni,nj
			HASH= (k-1)*ni*nj + (j-1)*ni +i + 1
		end function HASH

	end subroutine init_scomm
	subroutine destroy_scomm(scomm)
		implicit none
		type(scomm_t):: scomm

		!if(lcsrank==0) &
		!	write(*,*) 'in destroy_scomm...'

		scomm%connectivity = -1
		scomm%datatype= -1
		scomm%nbr_rank= -1
		scomm%flag= -1
		scomm%checker= -1
		scomm%pack_start= -1
		scomm%unpack_start= -1
		scomm%n_pack= -1
		scomm%n_unpack= -1
		scomm%pack_list_min = -1
		scomm%pack_list_max = -1
		scomm%unpack_list_min = -1
		scomm%unpack_list_max = -1
		scomm%periodic_shift = 0.0
		if(allocated(scomm%pack_buffer))   deallocate(scomm%pack_buffer)
		if(allocated(scomm%unpack_buffer)) deallocate(scomm%unpack_buffer)

	end subroutine destroy_scomm


	subroutine exchange_sdata(scomm,r0,r1,r2)
		implicit none
		!-----
		type(scomm_t):: scomm
		type(sr0_t),optional:: r0
		type(sr1_t),optional:: r1
		type(sr2_t),optional:: r2
		!-----
		integer:: icomm,ncomm,comm_id
		integer:: i,j,k,ii,jj,kk,ibuf
		integer:: ipack,iunpack
		integer:: tag_red,tag_black
		integer:: ierr, status(MPI_STATUS_SIZE)
		!-----
		! Exchange structured grid data or r0, r1,or r2 type
		!-----

		!
		! See what we have:
		!
		if(present(r0)) then
			if(lcsrank==0) &
				write(*,'(a,a)') 'In exchange_sdata... ',trim(r0%label)
		elseif(present(r1)) then
			if(lcsrank==0) &
				write(*,'(a,a)') 'In exchange_sdata... ',trim(r1%label)
		elseif(present(r2)) then
			if(lcsrank==0) &
				write(*,'(a,a)') 'In exchange_sdata... ',trim(r2%label)
		else
			if(lcsrank==0) &
				write(*,'(a,a)') 'In exchange_sdata... ', ' No data present'
			return
		endif

		!Set ncomm...
		if(scomm%connectivity == FACE_CONNECT) ncomm = 6
		if(scomm%connectivity == MAX_CONNECT) ncomm = 26

		!Allocate if we need to
		if(.NOT. PREALLOCATE_BUFFERS) then
			allocate(scomm%pack_buffer(scomm%pack_bufsize))
			allocate(scomm%unpack_buffer(scomm%unpack_bufsize))
		endif

		!
		! Pack the data:
		!
		do icomm = 1,ncomm
			i = COMM_OFFSET(1,icomm)
			j = COMM_OFFSET(2,icomm)
			k = COMM_OFFSET(3,icomm)

			if(scomm%flag(i,j,k) == NO_COMM) cycle

			ibuf = scomm%pack_start(i,j,k)
			do kk = scomm%pack_list_min(i,j,k,3),scomm%pack_list_max(i,j,k,3)
			do jj = scomm%pack_list_min(i,j,k,2),scomm%pack_list_max(i,j,k,2)
			do ii = scomm%pack_list_min(i,j,k,1),scomm%pack_list_max(i,j,k,1)
				select case(scomm%datatype)
				case(R0_COMM)
					scomm%pack_buffer(ibuf + 0) = r0%r(ii,jj,kk)
					ibuf = ibuf + 1
				case(R1_COMM)
					scomm%pack_buffer(ibuf + 0) = r1%x(ii,jj,kk)
					scomm%pack_buffer(ibuf + 1) = r1%y(ii,jj,kk)
					scomm%pack_buffer(ibuf + 2) = r1%z(ii,jj,kk)

					if(r1%periodic_translate) then
						scomm%pack_buffer(ibuf + 0) = r1%x(ii,jj,kk) + scomm%periodic_shift(i,j,k,1)
						scomm%pack_buffer(ibuf + 1) = r1%y(ii,jj,kk) + scomm%periodic_shift(i,j,k,2)
						scomm%pack_buffer(ibuf + 2) = r1%z(ii,jj,kk) + scomm%periodic_shift(i,k,k,3)
					else
						scomm%pack_buffer(ibuf + 0) = r1%x(ii,jj,kk)
						scomm%pack_buffer(ibuf + 1) = r1%y(ii,jj,kk)
						scomm%pack_buffer(ibuf + 2) = r1%z(ii,jj,kk)
					endif

					ibuf = ibuf + 3
				case(R2_COMM)
					scomm%pack_buffer(ibuf + 0) = r2%xx(ii,jj,kk)
					scomm%pack_buffer(ibuf + 1) = r2%xy(ii,jj,kk)
					scomm%pack_buffer(ibuf + 2) = r2%xz(ii,jj,kk)
					scomm%pack_buffer(ibuf + 3) = r2%yx(ii,jj,kk)
					scomm%pack_buffer(ibuf + 4) = r2%yy(ii,jj,kk)
					scomm%pack_buffer(ibuf + 5) = r2%yz(ii,jj,kk)
					scomm%pack_buffer(ibuf + 6) = r2%zx(ii,jj,kk)
					scomm%pack_buffer(ibuf + 7) = r2%zy(ii,jj,kk)
					scomm%pack_buffer(ibuf + 8) = r2%zz(ii,jj,kk)
					ibuf = ibuf + 9
				end select
			enddo
			enddo
			enddo
		enddo

		!
		! Send/Recv
		!
		comm_id = 0
		do icomm = 1,ncomm
			i = COMM_OFFSET(1,icomm)
			j = COMM_OFFSET(2,icomm)
			k = COMM_OFFSET(3,icomm)

			!Index the tags:
			comm_id = comm_id +1
			tag_red = TAG_START + comm_id
			tag_black = 2*TAG_START + comm_id

			!Red sends first then receives:
			if (scomm%checker(i,j,k) == RED) then
				if(scomm%flag(i,j,k) == NBR_COMM) then
					call MPI_SEND(scomm%pack_buffer(scomm%pack_start(i,j,k)),scomm%n_pack(i,j,k),&
						MPI_INTEGER,scomm%nbr_rank(i,j,k),TAG_RED,lcscomm,ierr)
				endif
				if (scomm%flag(-i,-j,-k) == NBR_COMM) then
					call MPI_RECV(scomm%unpack_buffer(scomm%unpack_start(-i,-j,-k)),scomm%n_unpack(-i,-j,-k),&
						MPI_INTEGER,scomm%nbr_rank(-i,-j,-k),TAG_BLACK,lcscomm,status,ierr)
				endif
			endif

			!Black receives first then sends:
			if (scomm%checker(i,j,k) == BLACK) then
				if (scomm%flag(-i,-j,-k) == NBR_COMM) then
					call MPI_RECV(scomm%unpack_buffer(scomm%unpack_start(-i,-j,-k)),scomm%n_unpack(-i,-j,-k),&
						MPI_INTEGER,scomm%nbr_rank(-i,-j,-k),TAG_RED,lcscomm,status,ierr)
				endif
				if(scomm%flag(i,j,k) == NBR_COMM) then
					call MPI_SEND(scomm%pack_buffer(scomm%pack_start(i,j,k)),scomm%n_pack(i,j,k),&
						MPI_INTEGER,scomm%nbr_rank(i,j,k),TAG_BLACK,lcscomm,ierr)
				endif
			endif

			!Self communication:
			if (scomm%flag(i,j,k) == SELF_COMM) then
				ipack = scomm%pack_start(i,j,k)
				iunpack = scomm%unpack_start(-i,-j,-k)
				do  ibuf = 1,scomm%n_unpack(-i,-j,-k)
					scomm%unpack_buffer(iunpack) = scomm%pack_buffer(ipack)
					ipack = ipack + 1
					iunpack = iunpack + 1
				enddo
			endif
		enddo

		!
		! Unpack...
		!
		do icomm = 1,ncomm
			i = COMM_OFFSET(1,icomm)
			j = COMM_OFFSET(2,icomm)
			k = COMM_OFFSET(3,icomm)
			if(scomm%flag(i,j,k) == NO_COMM) cycle

			ibuf = scomm%unpack_start(i,j,k)
			do kk = scomm%unpack_list_min(i,j,k,3),scomm%unpack_list_max(i,j,k,3)
			do jj = scomm%unpack_list_min(i,j,k,2),scomm%unpack_list_max(i,j,k,2)
			do ii = scomm%unpack_list_min(i,j,k,1),scomm%unpack_list_max(i,j,k,1)
				select case(scomm%datatype)
				case(R0_COMM)
					r0%r(ii,jj,kk) = scomm%unpack_buffer(ibuf + 0)
					ibuf = ibuf + 1
				case(R1_COMM)
					r1%x(ii,jj,kk) = scomm%unpack_buffer(ibuf + 0)
					r1%y(ii,jj,kk) = scomm%unpack_buffer(ibuf + 1)
					r1%z(ii,jj,kk) = scomm%unpack_buffer(ibuf + 2)
					ibuf = ibuf + 3
				case(R2_COMM)
					r2%xx(ii,jj,kk) = scomm%unpack_buffer(ibuf + 0)
					r2%xy(ii,jj,kk) = scomm%unpack_buffer(ibuf + 1)
					r2%xz(ii,jj,kk) = scomm%unpack_buffer(ibuf + 2)
					r2%yx(ii,jj,kk) = scomm%unpack_buffer(ibuf + 3)
					r2%yy(ii,jj,kk) = scomm%unpack_buffer(ibuf + 4)
					r2%yz(ii,jj,kk) = scomm%unpack_buffer(ibuf + 5)
					r2%zx(ii,jj,kk) = scomm%unpack_buffer(ibuf + 6)
					r2%zy(ii,jj,kk) = scomm%unpack_buffer(ibuf + 7)
					r2%zz(ii,jj,kk) = scomm%unpack_buffer(ibuf + 8)
					ibuf = ibuf + 9
				end select
			enddo
			enddo
			enddo
		enddo

		if(.NOT. PREALLOCATE_BUFFERS) then
			deallocate(scomm%pack_buffer)
			deallocate(scomm%unpack_buffer)
		endif

	end subroutine exchange_sdata

end module comms_m
!write(*,'(a,i4,a,i4,a,i6)') 'lcsrank[', lcsrank,'] sending to [',scomm%nbr_rank(i,j,k),'] tag=',TAG_RED
!write(*,'(a,i4,a,i4,a,i6)') 'lcsrank[', lcsrank,'] recv from [',scomm%nbr_rank(-i,-j,-k),'] tag=',TAG_BLACK
!write(*,'(a,i4,a,i4,a,i6)') 'lcsrank[', lcsrank,'] recv from [',scomm%nbr_rank(-i,-j,-k),'] tag=',TAG_RED
!write(*,'(a,i4,a,i4,a,i6)') 'lcsrank[', lcsrank,'] send to [',scomm%nbr_rank(i,j,k),'] tag=',TAG_BLACK

