module comms_m
	use data_m
	implicit none

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

	! Set the Preallocate flag for structured communication buffers
	logical,parameter:: PREALLOCATE_BUFFERS = .TRUE.


	!Set the All-to-all strategy
	logical,parameter:: ALLTOALL_REDUCE = .TRUE.

	contains

	subroutine init_lcs_mpi(cfdcomm)
		implicit none
		!----
		integer, intent(in):: cfdcomm
		!----
		integer:: ierr
		!----

		!
		!Duplicate the MPI comm, so that we have our
		!own communications for LCS related stuff.
		!

		call MPI_COMM_DUP(cfdcomm,lcscomm,ierr)
		call MPI_COMM_RANK(lcscomm,lcsrank,ierr)
		call MPI_COMM_SIZE(lcscomm,nprocs,ierr)

		if (lcsrank==0) &
			write(*,'(a,i6,a)') 'in init_lcs_mpi... Using ',nprocs, ' MPI processes.'

		!
		! Set the real precision
		!
		select case (LCSRP)
		case (4)
			MPI_LCSRP = MPI_REAL
		case (8)
			MPI_LCSRP = MPI_DOUBLE_PRECISION
		case default
			if (lcsrank==0) &
				write(*,'(a,i6,a)') 'Cant figure out MPI real precision', LCSRP
			stop
		end select

	end subroutine init_lcs_mpi

	subroutine init_scomm(scomm,ni,nj,nk,ng,offset_i,offset_j,offset_k,per_i,per_j,per_k,connectivity,datatype,label)
		implicit none
		!-----
		type(scomm_t):: scomm
		integer:: ni,nj,nk,ng
		integer:: offset_i,offset_j,offset_k
		logical:: per_i,per_j,per_k
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
			write(*,*) 'in init_scomm... ',trim(scomm%label),per_i,per_j,per_k

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
			!Global max/min of grid index: 
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
			if(i0(lcsrank) == imin .and. per_i)im1	= imax
			if(i1(lcsrank) == imax .and. per_i)ip1	= imin
			if(j0(lcsrank) == jmin .and. per_j)jm1	= jmax
			if(j1(lcsrank) == jmax .and. per_j)jp1	= jmin
			if(k0(lcsrank) == kmin .and. per_k)km1	= kmax
			if(k1(lcsrank) == kmax .and. per_k)kp1	= kmin
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
						if(rank == lcsrank) then
							scomm%flag(i,j,k) = SELF_COMM
						else
							scomm%flag(i,j,k) = NBR_COMM
						endif
					!write(*,*) 'lcsrank[',lcsrank,'] has comm neighbor in direction',i,j,k,per_i,per_j,per_k
					endif
				enddo
				if (found > 1) then
					write(*,*) 'ERROR:  lcsrank[',lcsrank,'] Non-unique comm pattern in direction',i,j,k
				endif
				if (found ==0) then
					scomm%flag(i,j,k) = NO_COMM
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

			if(scomm%connectivity == FACE_CONNECT) ncomm = 7
			if(scomm%connectivity == MAX_CONNECT) ncomm = 27
			scomm%n_pack = 0
			scomm%n_unpack = 0
			scomm%pack_bufsize = 0
			scomm%unpack_bufsize = 0
			scomm%pack_start = 0
			scomm%unpack_start = 0
			do icomm = 2,ncomm
				i = NBR_OFFSET(1,icomm)
				j = NBR_OFFSET(2,icomm)
				k = NBR_OFFSET(3,icomm)
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
			if(lcsrank==0 .AND. LCS_VERBOSE) then
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
			integer:: ihash,jhash,khash
			integer(8):: hashval,bufval
			integer:: ipack,iunpack
			!-----
			! A basic check to make sure that all point-to-point communications
			! are working as expected.  This is also a skeleton for all other
			! new comm routines.
			!-----

			if(lcsrank==0 .AND. LCS_VERBOSE) &
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
			if(scomm%connectivity == FACE_CONNECT) ncomm = 7
			if(scomm%connectivity == MAX_CONNECT) ncomm = 27
			do icomm = 2,ncomm

				i = NBR_OFFSET(1,icomm)
				j = NBR_OFFSET(2,icomm)
				k = NBR_OFFSET(3,icomm)

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
			if(scomm%connectivity == FACE_CONNECT) ncomm = 7
			if(scomm%connectivity == MAX_CONNECT) ncomm = 27
			do icomm = 2,ncomm
				i = NBR_OFFSET(1,icomm)
				j = NBR_OFFSET(2,icomm)
				k = NBR_OFFSET(3,icomm)

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
			do icomm = 2,ncomm
				i = NBR_OFFSET(1,icomm)
				j = NBR_OFFSET(2,icomm)
				k = NBR_OFFSET(3,icomm)
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
					bufval = nint(scomm%unpack_buffer(ibuf),8) !JRF CHANGED TO NINT
					!************

					if(hashval /= bufval) then
						write(*,*) 'ERROR: lcsrank[',lcsrank,'] Has unexpected value in unpack buffer',hashval,bufval
						CFD2LCS_ERROR = 1
						return
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
					if(LCS_VERBOSE) write(*,*) '     handshake OK for ',nsend, 'communications'
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

		integer(8) function HASH(i,j,k,ni,nj)
			implicit none
			integer::i,j,k,ni,nj
			HASH= (k-1)*ni*nj + (j-1)*ni +i !+ 1  JRF GET RID OF +1
		end function HASH

	end subroutine init_scomm
	subroutine destroy_scomm(scomm)
		implicit none
		type(scomm_t):: scomm
		!-----
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
			if(lcsrank==0 .AND. LCS_VERBOSE) &
				write(*,'(a,a)') 'In exchange_sdata... ',trim(r0%label)
		elseif(present(r1)) then
			if(lcsrank==0 .AND. LCS_VERBOSE) &
				write(*,'(a,a)') 'In exchange_sdata... ',trim(r1%label)
		elseif(present(r2)) then
			if(lcsrank==0 .AND. LCS_VERBOSE) &
				write(*,'(a,a)') 'In exchange_sdata... ',trim(r2%label)
		else
			if(lcsrank==0) &
				write(*,'(a,a)') 'In exchange_sdata... ', ' No data present'
			return
		endif

		!Set ncomm...
		if(scomm%connectivity == FACE_CONNECT) ncomm = 7
		if(scomm%connectivity == MAX_CONNECT) ncomm = 27

		!Allocate if we need to
		if(.NOT. PREALLOCATE_BUFFERS) then
			allocate(scomm%pack_buffer(scomm%pack_bufsize))
			allocate(scomm%unpack_buffer(scomm%unpack_bufsize))
		endif

		!
		! Pack the data:
		!
		do icomm = 2,ncomm
			i = NBR_OFFSET(1,icomm)
			j = NBR_OFFSET(2,icomm)
			k = NBR_OFFSET(3,icomm)

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
		do icomm = 2,ncomm
			i = NBR_OFFSET(1,icomm)
			j = NBR_OFFSET(2,icomm)
			k = NBR_OFFSET(3,icomm)

			!Index the tags:
			comm_id = comm_id +1
			tag_red = TAG_START + comm_id
			tag_black = 2*TAG_START + comm_id

			!Red sends first then receives:
			if (scomm%checker(i,j,k) == RED) then
				if(scomm%flag(i,j,k) == NBR_COMM) then
					call MPI_SEND(scomm%pack_buffer(scomm%pack_start(i,j,k)),scomm%n_pack(i,j,k),&
						MPI_LCSRP,scomm%nbr_rank(i,j,k),TAG_RED,lcscomm,ierr)
				endif
				if (scomm%flag(-i,-j,-k) == NBR_COMM) then
					call MPI_RECV(scomm%unpack_buffer(scomm%unpack_start(-i,-j,-k)),scomm%n_unpack(-i,-j,-k),&
						MPI_LCSRP,scomm%nbr_rank(-i,-j,-k),TAG_BLACK,lcscomm,status,ierr)
				endif
			endif

			!Black receives first then sends:
			if (scomm%checker(i,j,k) == BLACK) then
				if (scomm%flag(-i,-j,-k) == NBR_COMM) then
					call MPI_RECV(scomm%unpack_buffer(scomm%unpack_start(-i,-j,-k)),scomm%n_unpack(-i,-j,-k),&
						MPI_LCSRP,scomm%nbr_rank(-i,-j,-k),TAG_RED,lcscomm,status,ierr)
				endif
				if(scomm%flag(i,j,k) == NBR_COMM) then
					call MPI_SEND(scomm%pack_buffer(scomm%pack_start(i,j,k)),scomm%n_pack(i,j,k),&
						MPI_LCSRP,scomm%nbr_rank(i,j,k),TAG_BLACK,lcscomm,ierr)
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
		do icomm = 2,ncomm
			i = NBR_OFFSET(1,icomm)
			j = NBR_OFFSET(2,icomm)
			k = NBR_OFFSET(3,icomm)
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


	subroutine exchange_lpdata(lp,sgrid,no)
		use lp_m
		implicit none
		!-----
		type(lp_t):: lp
		type(sgrid_t):: sgrid
		type(ui1_t):: no !lp node corresponding to sgrid
		!-----
		integer,parameter:: NCOMM_LP = 27
		integer,parameter:: NREAL_LPCOMM = 17  !x,y,z,u,v,w,dx,dy,dz,no0,proc0,no-x,no-y,no-z,no_scfd-x,no_scfd-y,no_scfd-z
		real(LCSRP),parameter::MAGIC_INIT = 12345678.0_LCSRP
		!-----
		integer:: ip
		integer:: commflag(1:3,1:lp%np)
		integer:: pack_buffer_size,unpack_buffer_size
		integer:: np_pack_total, np_unpack_total,pack_global,unpack_global
		integer:: comm_id, i,j,k, icomm,tag_red,tag_black
		integer:: ierr, status(MPI_STATUS_SIZE)
		integer:: pack_end,unpack_end
		integer:: ipack,iunpack,ibuf
		!-----
		integer:: nbr_rank(-1:1,-1:1,-1:1) !Rank of processor in each direction
		integer:: flag(-1:1,-1:1,-1:1) !tells us what to do in each direction
		integer:: checker(-1:1,-1:1,-1:1) !Checkerboard pattern for each comm direction
		integer:: np_pack(-1:1,-1:1,-1:1), np_unpack(-1:1,-1:1,-1:1) !num particles to pack/unpack
		integer:: pack_start(-1:1,-1:1,-1:1), unpack_start(-1:1,-1:1,-1:1), tmp(-1:1,-1:1,-1:1)  !First index into pack/unpack buffers
		real(LCSRP), allocatable :: pack_buffer(:), unpack_buffer(:)  !Exchange buffers
		real(LCSRP):: periodic_shift(-1:1,-1:1,-1:1,1:3)  !For shifting coordinates across periodic boundaries
		!-----
		!Exchange data based on the node of lp
		!-----

		if(lcsrank==0 .AND. LCS_VERBOSE) &
			write(*,*) 'in exchange_lpdata...',trim(lp%label)

		!-----
		!Copy some atributes from the scomm
		!-----
		nbr_rank = sgrid%scomm_max_r1%nbr_rank
		flag = sgrid%scomm_max_r1%flag
		flag(0,0,0) = NO_COMM !explicitly set no comm for 0,0,0
		checker = sgrid%scomm_max_r1%checker
		periodic_shift = sgrid%ps(lcsrank,:,:,:,:)

		!-----
		!Determine what direction each particle goes:
		!and set the pack buffer size.  Careful with NO_COMM bc.
		!-----
		commflag = 0
		np_pack = 0
		np_unpack = 0
		do ip = 1,lp%np
			if(no%x(ip)<1) commflag(1,ip) = -1
			if(no%x(ip)>sgrid%ni) commflag(1,ip) = 1
			if(no%y(ip)<1) commflag(2,ip) = -1
			if(no%y(ip)>sgrid%nj) commflag(2,ip) = 1
			if(no%z(ip)<1) commflag(3,ip) = -1
			if(no%z(ip)>sgrid%nk) commflag(3,ip) = 1
			np_pack(commflag(1,ip),commflag(2,ip),commflag(3,ip)) = &
				np_pack(commflag(1,ip),commflag(2,ip),commflag(3,ip)) + 1
		enddo
		do icomm = 2,NCOMM_LP
			i = NBR_OFFSET(1,icomm)
			j = NBR_OFFSET(2,icomm)
			k = NBR_OFFSET(3,icomm)
			if (flag(i,j,k) == NO_COMM)	np_pack(i,j,k) = 0
		enddo
		np_pack(0,0,0) = 0
		np_pack_total = sum(np_pack)
		pack_buffer_size = NREAL_LPCOMM*np_pack_total
		allocate(pack_buffer(1:pack_buffer_size))
		pack_buffer = MAGIC_INIT

		!-----
		!Exchange pack buffer size
		!-----
		comm_id = 0
		do icomm = 2,NCOMM_LP
			i = NBR_OFFSET(1,icomm)
			j = NBR_OFFSET(2,icomm)
			k = NBR_OFFSET(3,icomm)
			!Index the tags:
			comm_id = comm_id +1
			tag_red = TAG_START + comm_id
			tag_black = 2*TAG_START + comm_id

			!Red sends first then receives:
			if (checker(i,j,k) == RED) then
				if(flag(i,j,k) == NBR_COMM) then
					call MPI_SEND(np_pack(i,j,k),1,MPI_INTEGER,nbr_rank(i,j,k),TAG_RED,lcscomm,ierr)
				endif
				if (flag(-i,-j,-k) == NBR_COMM) then
					call MPI_RECV(np_unpack(-i,-j,-k),1,MPI_INTEGER,nbr_rank(-i,-j,-k),TAG_BLACK,lcscomm,status,ierr)
				endif
			endif
			!Black receives first then sends:
			if (checker(i,j,k) == BLACK) then
				if (flag(-i,-j,-k) == NBR_COMM) then
					call MPI_RECV(np_unpack(-i,-j,-k),1,MPI_INTEGER,nbr_rank(-i,-j,-k),TAG_RED,lcscomm,status,ierr)
				endif
				if(flag(i,j,k) == NBR_COMM) then
					call MPI_SEND(np_pack(i,j,k),1,MPI_INTEGER,nbr_rank(i,j,k),TAG_BLACK,lcscomm,ierr)
				endif
			endif
			!Self communication:
			if (flag(i,j,k) == SELF_COMM) then
				np_unpack(i,j,k) = np_pack(-i,-j,-k)
			endif
		enddo
		np_unpack_total = sum(np_unpack)
		unpack_buffer_size = NREAL_LPCOMM*np_unpack_total
		allocate(unpack_buffer(1:unpack_buffer_size))
		unpack_buffer = MAGIC_INIT

		!-----
		!Check the buffer sizes match:
		!-----
		if(LCS_VERBOSE) then
			call MPI_ALLREDUCE(np_pack_total,pack_global,1,MPI_INTEGER,MPI_SUM,lcscomm,ierr)
			call MPI_ALLREDUCE(np_unpack_total,unpack_global,1,MPI_INTEGER,MPI_SUM,lcscomm,ierr)
			if(pack_global /= unpack_global) then
				if(lcsrank==0)&
					write(*,*) 'ERROR:  global pack/unpack not the same:',pack_global,unpack_global
				CFD2LCS_ERROR = 1
				return
			else
				if(lcsrank==0)&
					write(*,*) 'Exchange of:',pack_global,'particles'
			endif
			if(pack_global==0 ) return
		endif

		!-----
		!Set the buffer start point for each communication
		!-----
		pack_start(:,:,:) = 0
		unpack_start(:,:,:) = 0
		pack_end = 0
		unpack_end = 0
		do icomm = 2,NCOMM_LP
			i = NBR_OFFSET(1,icomm)
			j = NBR_OFFSET(2,icomm)
			k = NBR_OFFSET(3,icomm)
			if (flag(i,j,k) == NO_COMM) cycle
			if(np_pack(i,j,k) > 0)then
				pack_start(i,j,k) = pack_end + 1
				pack_end = pack_end + NREAL_LPCOMM*np_pack(i,j,k)
			endif
			if(np_unpack(i,j,k) > 0)then
				unpack_start(i,j,k) = unpack_end + 1
				unpack_end = unpack_end + NREAL_LPCOMM*np_unpack(i,j,k)
			endif
		enddo
		if (pack_end /= pack_buffer_size) then
			write(*,*) 'lcsrank[',lcsrank,'] ERROR: pack buffer inconsistent',pack_end, pack_buffer_size
			CFD2LCS_ERROR = 1
		endif
		if (unpack_end /= unpack_buffer_size) then
			write(*,*) 'lcsrank[',lcsrank,'] ERROR: unpack buffer inconsistent',unpack_end, unpack_buffer_size
			CFD2LCS_ERROR = 1
		endif

		!-----
		!fill the pack buffer
		!We pass integers as real values, and then convert back using the nint function on unpack
		!Pass the global node index, paying attention to the periodicity.
		!once the data is packed, flag the particle to be recycled.
		!-----
		tmp = pack_start
		do ip = 1,lp%np
			i = commflag(1,ip)
			j = commflag(2,ip)
			k = commflag(3,ip)
			if(flag(i,j,k) == NO_COMM) cycle
			pack_buffer(pack_start(i,j,k)+0) = lp%xp%x(ip) + periodic_shift(i,j,k,1)
			pack_buffer(pack_start(i,j,k)+1) = lp%xp%y(ip) + periodic_shift(i,j,k,2)
			pack_buffer(pack_start(i,j,k)+2) = lp%xp%z(ip) + periodic_shift(i,j,k,3)
			pack_buffer(pack_start(i,j,k)+3) = lp%up%x(ip)
			pack_buffer(pack_start(i,j,k)+4) = lp%up%y(ip)
			pack_buffer(pack_start(i,j,k)+5) = lp%up%z(ip)
			pack_buffer(pack_start(i,j,k)+6) = real(lp%no0%i(ip),LCSRP)
			pack_buffer(pack_start(i,j,k)+7) = real(lp%proc0%i(ip),LCSRP)
			!no
			pack_buffer(pack_start(i,j,k)+8) =  real(sgrid%offset_i+lp%no%x(ip),LCSRP)
			pack_buffer(pack_start(i,j,k)+9) =  real(sgrid%offset_j+lp%no%y(ip),LCSRP)
			pack_buffer(pack_start(i,j,k)+10) = real(sgrid%offset_k+lp%no%z(ip),LCSRP)
			if(pack_buffer(pack_start(i,j,k)+8) > sgrid%gni) &
				pack_buffer(pack_start(i,j,k)+8) =  pack_buffer(pack_start(i,j,k)+8) - real(sgrid%gni)
			if(pack_buffer(pack_start(i,j,k)+9) > sgrid%gnj) &
				pack_buffer(pack_start(i,j,k)+9) =  pack_buffer(pack_start(i,j,k)+9) - real(sgrid%gnj)
			if(pack_buffer(pack_start(i,j,k)+10) > sgrid%gnk) &
				pack_buffer(pack_start(i,j,k)+10) =  pack_buffer(pack_start(i,j,k)+10) - real(sgrid%gnk)
			if(pack_buffer(pack_start(i,j,k)+8) < 1) &
				pack_buffer(pack_start(i,j,k)+8) =  pack_buffer(pack_start(i,j,k)+8) + real(sgrid%gni)
			if(pack_buffer(pack_start(i,j,k)+9) < 1) &
				pack_buffer(pack_start(i,j,k)+9) =  pack_buffer(pack_start(i,j,k)+9) + real(sgrid%gnj)
			if(pack_buffer(pack_start(i,j,k)+10) < 1) &
				pack_buffer(pack_start(i,j,k)+10) =  pack_buffer(pack_start(i,j,k)+10) + real(sgrid%gnk)
			!dx
			pack_buffer(pack_start(i,j,k)+11) = lp%dx%x(ip)
			pack_buffer(pack_start(i,j,k)+12) = lp%dx%y(ip)
			pack_buffer(pack_start(i,j,k)+13) = lp%dx%z(ip)
			!no_scfd
			pack_buffer(pack_start(i,j,k)+14) =  real(scfd%sgrid%offset_i+lp%no_scfd%x(ip),LCSRP)
			pack_buffer(pack_start(i,j,k)+15) =  real(scfd%sgrid%offset_j+lp%no_scfd%y(ip),LCSRP)
			pack_buffer(pack_start(i,j,k)+16) =  real(scfd%sgrid%offset_k+lp%no_scfd%z(ip),LCSRP)
			if(pack_buffer(pack_start(i,j,k)+14) > scfd%sgrid%gni) &
				pack_buffer(pack_start(i,j,k)+14) =  pack_buffer(pack_start(i,j,k)+14) - real(scfd%sgrid%gni)
			if(pack_buffer(pack_start(i,j,k)+15) > scfd%sgrid%gnj) &
				pack_buffer(pack_start(i,j,k)+15) =  pack_buffer(pack_start(i,j,k)+15) - real(scfd%sgrid%gnj)
			if(pack_buffer(pack_start(i,j,k)+16) > scfd%sgrid%gnk) &
				pack_buffer(pack_start(i,j,k)+16) =  pack_buffer(pack_start(i,j,k)+16) - real(scfd%sgrid%gnk)
			if(pack_buffer(pack_start(i,j,k)+14) < 1) &
				pack_buffer(pack_start(i,j,k)+14) =  pack_buffer(pack_start(i,j,k)+14) + real(scfd%sgrid%gni)
			if(pack_buffer(pack_start(i,j,k)+15) < 1) &
				pack_buffer(pack_start(i,j,k)+15) =  pack_buffer(pack_start(i,j,k)+15) + real(scfd%sgrid%gnj)
			if(pack_buffer(pack_start(i,j,k)+16) < 1) &
				pack_buffer(pack_start(i,j,k)+16) =  pack_buffer(pack_start(i,j,k)+16) + real(scfd%sgrid%gnk)

			pack_start(i,j,k) = pack_start(i,j,k)+NREAL_LPCOMM
			lp%flag%i(ip) = LP_RECYCLE
		enddo
		pack_start = tmp

		!-----
		!Exchange buffers
		!-----
		comm_id = 0
		do icomm = 2,NCOMM_LP
			i = NBR_OFFSET(1,icomm)
			j = NBR_OFFSET(2,icomm)
			k = NBR_OFFSET(3,icomm)
			!Index the tags:
			comm_id = comm_id +1
			tag_red = TAG_START + comm_id
			tag_black = 2*TAG_START + comm_id

			!Red sends first then receives:
			if (checker(i,j,k) == RED) then
				if(flag(i,j,k) == NBR_COMM .AND. np_pack(i,j,k) > 0) then
					call MPI_SEND(pack_buffer(pack_start(i,j,k)),np_pack(i,j,k)*NREAL_LPCOMM,&
					MPI_LCSRP,nbr_rank(i,j,k),TAG_RED,lcscomm,ierr)
				endif
				if (flag(-i,-j,-k) == NBR_COMM .AND. np_unpack(-i,-j,-k) > 0) then
					call MPI_RECV(unpack_buffer(unpack_start(-i,-j,-k)),np_unpack(-i,-j,-k)*NREAL_LPCOMM,&
					MPI_LCSRP,nbr_rank(-i,-j,-k),TAG_BLACK,lcscomm,status,ierr)
				endif
			endif
			!Black receives first then sends:
			if (checker(i,j,k) == BLACK) then
				if (flag(-i,-j,-k) == NBR_COMM .AND. np_unpack(-i,-j,-k) > 0) then
					call MPI_RECV(unpack_buffer(unpack_start(-i,-j,-k)),np_unpack(-i,-j,-k)*NREAL_LPCOMM,&
					MPI_LCSRP,nbr_rank(-i,-j,-k),TAG_RED,lcscomm,status,ierr)
				endif
				if(flag(i,j,k) == NBR_COMM .and. np_pack(i,j,k) > 0) then
					call MPI_SEND(pack_buffer(pack_start(i,j,k)),np_pack(i,j,k)*NREAL_LPCOMM,&
					MPI_LCSRP,nbr_rank(i,j,k),TAG_BLACK,lcscomm,ierr)
				endif
			endif
			!Self communication:
			if (flag(i,j,k) == SELF_COMM .AND. np_pack(i,j,k) > 0) then
				ipack = pack_start(i,j,k)
				iunpack = unpack_start(-i,-j,-k)
				do  ibuf = 1,np_unpack(-i,-j,-k)*NREAL_LPCOMM
					unpack_buffer(iunpack) = pack_buffer(ipack)
					ipack = ipack + 1
					iunpack = iunpack + 1
				enddo
			endif
		enddo

		!-----
		!Reorder and remove holes in the lp list
		!-----
		if(np_pack_total > 0) then
			call reorder_lp(lp)
		endif

		!-----
		!Resize lp
		!-----
		ip = lp%np  !starting point for unpacking
		call resize_lp(lp,lp%np+np_unpack_total)

		!-----
		!Unpack data
		!-----
		ibuf = 1
		do iunpack = 1,np_unpack_total
			ip = ip + 1
			lp%xp%x(ip) = unpack_buffer(ibuf+0)
			lp%xp%y(ip) = unpack_buffer(ibuf+1)
			lp%xp%z(ip) = unpack_buffer(ibuf+2)
			lp%up%x(ip) = unpack_buffer(ibuf+3)
			lp%up%y(ip) = unpack_buffer(ibuf+4)
			lp%up%z(ip) = unpack_buffer(ibuf+5)
			lp%no0%i(ip)= nint(unpack_buffer(ibuf+6))
			lp%proc0%i(ip)= nint(unpack_buffer(ibuf+7))
			lp%no%x(ip) = nint(unpack_buffer(ibuf+8)) - sgrid%offset_i  !convert back to local ind.
			lp%no%y(ip) = nint(unpack_buffer(ibuf+9)) - sgrid%offset_j	!convert back to local ind.
			lp%no%z(ip) = nint(unpack_buffer(ibuf+10)) - sgrid%offset_k !convert back to local ind.
			lp%dx%x(ip) = unpack_buffer(ibuf+11)
			lp%dx%y(ip) = unpack_buffer(ibuf+12)
			lp%dx%z(ip) = unpack_buffer(ibuf+13)
			lp%no_scfd%x(ip) = nint(unpack_buffer(ibuf+14)) - scfd%sgrid%offset_i  !convert back to local ind.
			lp%no_scfd%y(ip) = nint(unpack_buffer(ibuf+15)) - scfd%sgrid%offset_j	!convert back to local ind.
			lp%no_scfd%z(ip) = nint(unpack_buffer(ibuf+16)) - scfd%sgrid%offset_k !convert back to local ind.
			!Set the flag
			lp%flag%i(ip) = LP_IB
			ibuf = ibuf + NREAL_LPCOMM

		enddo

		!Cleanup
		deallocate(pack_buffer)
		deallocate(unpack_buffer)

	end subroutine exchange_lpdata

	subroutine exchange_lpmap(lp)
		implicit none
		!-----
		type(lp_t):: lp
		!-----
		integer:: ip,ierr,proc,is,ir,i,j,k,ibuf
		integer,allocatable:: nsend(:), nrecv(:)
		integer,allocatable:: sendstart(:), recvstart(:), tmp(:)
		real(LCSRP),allocatable:: sendbuf(:),recvbuf(:)
		integer:: nsend_total,nrecv_total
		integer:: r_start,r_end, s_start,s_end
		integer:: tag_f  !unique tag for each send/rec pair
		integer :: status(MPI_STATUS_SIZE)
		integer,parameter:: LPMAP_SIZE = 4  !xp,yp,zp,real(node0)
		integer:: no0
		integer:: visited(1:lp%fm%ni,1:lp%fm%nj,1:lp%fm%nk)
		real:: t0,t1
		!-----

		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,'(a)') 'in exchange_lpmap...'


		t0 = cputimer(lcscomm,SYNC_TIMER)

		!-----
		!Count whose particles each proc has and exchange with all
		!-----
		allocate(nsend(0:nprocs-1))
		allocate(nrecv(0:nprocs-1))
		allocate(sendstart(0:nprocs-1))
		allocate(recvstart(0:nprocs-1))
		allocate(tmp(0:nprocs-1))
		nsend = 0; nrecv = 0
		sendstart=-1; recvstart = -1
		do ip = 1,lp%np
			nsend(lp%proc0%i(ip)) = nsend(lp%proc0%i(ip)) + 1
		enddo

		!-----
		!Now communicate to receiving processors
		!-----
		if(ALLTOALL_REDUCE) then
			!Use several reduce ops:
			do proc = 0, nprocs-1
				tmp = 0
				tmp(lcsrank) = nsend(proc)
				call MPI_REDUCE(tmp,nrecv,nprocs,MPI_INTEGER,MPI_SUM,proc,lcscomm,ierr)
			enddo
		else
			!One-to-one comms:
			do proc = 0, nprocs-1
				if(proc == lcsrank) then
					nrecv(lcsrank) = nsend(lcsrank)
				elseif (lcsrank < proc) then
					! lower ranks first send and then receive...
					tag_f = TAG_START + proc
					call MPI_SEND(nsend(proc),1,MPI_INTEGER,proc,tag_f,lcscomm,ierr)
					tag_f = TAG_START + lcsrank
					call MPI_RECV(nrecv(proc),1,MPI_INTEGER,proc,tag_f,lcscomm,status,ierr)
				else
					! higher ranks first receive and then send...
					tag_f = TAG_START + lcsrank
					call MPI_RECV(nrecv(proc),1,MPI_INTEGER,proc,tag_f,lcscomm,status,ierr)
					tag_f = TAG_START + proc
					call MPI_SEND(nsend(proc),1,MPI_INTEGER,proc,tag_f,lcscomm,ierr)
				endif
			end do
		endif
		nsend_total = sum(nsend)
		nrecv_total = sum(nrecv)
		if(nrecv_total /= lp%fm%ni*lp%fm%nj*lp%fm%nk) then
			write(*,*) lcsrank,'ERROR: incomplete map returned:',nrecv_total,lp%fm%ni*lp%fm%nj*lp%fm%nk
		endif

		!-----
		!Determine where the send/recieve will start for each rank
		!This will allow for easy indexing into a single vector when communicating.
		!-----
		sendstart(:) = -1
		recvstart(:) = -1
		is = 0; ir = 0
		do proc = 0, nprocs-1
			if (nsend(proc) > 0) then
				sendstart(proc) = is+1
				is = is + nsend(proc)*LPMAP_SIZE
			endif
			if (nrecv(proc) > 0) then
				recvstart(proc) = ir+1
				ir = ir + nrecv(proc)*LPMAP_SIZE
			endif
		enddo


		!check
		if (is /= nsend_total*LPMAP_SIZE) then
			write(*,*) '[',lcsrank,'] ERROR, expecting to send',nsend_total*LPMAP_SIZE,'data.  Will actually send', is
			CFD2LCS_ERROR = 1
		endif
		if (ir /= nrecv_total*LPMAP_SIZE) then
			write(*,*) '[',lcsrank,'] ERROR, expecting to recv',nrecv_total*LPMAP_SIZE,'data.  Will actually recv', ir
			CFD2LCS_ERROR = 1
		endif

		!-----
		!Pack the buffers
		!Note, we send dx and not xp to handle periodic domains.
		!-----
		allocate(sendbuf(1:nsend_total*LPMAP_SIZE))
		allocate(recvbuf(1:nrecv_total*LPMAP_SIZE))
		tmp = sendstart
		do ip = 1,lp%np
			proc = lp%proc0%i(ip)
			sendbuf(sendstart(proc)+0) = real(lp%no0%i(ip),LCSRP)
			sendbuf(sendstart(proc)+1) = lp%dx%x(ip)
			sendbuf(sendstart(proc)+2) = lp%dx%y(ip)
			sendbuf(sendstart(proc)+3) = lp%dx%z(ip)
			sendstart(proc) = sendstart(proc) + LPMAP_SIZE
		enddo
		sendstart = tmp

		!-----
		!Exchange buffers
		!-----
		do proc = 0, nprocs-1
			s_start = sendstart(proc)
			s_end = s_start + nsend(proc) * LPMAP_SIZE-1
			r_start = recvstart(proc)
			r_end = r_start + nrecv(proc) * LPMAP_SIZE-1
			if(proc == lcsrank .AND. nsend(proc)>0) then
				recvbuf(r_start:r_end) = sendbuf(s_start:s_end)
			elseif (lcsrank < proc) then
				! lower ranks first send and then receive...
				if(nsend(proc) >0)then
					tag_f = TAG_START + proc
					call MPI_SEND(sendbuf(s_start),nsend(proc)*LPMAP_SIZE,MPI_LCSRP,proc,tag_f,lcscomm,ierr)
				endif
				if(nrecv(proc) >0)then
					tag_f = TAG_START + lcsrank
					call MPI_RECV(recvbuf(r_start),nrecv(proc)*LPMAP_SIZE,MPI_LCSRP,proc,tag_f,lcscomm,status,ierr)
				endif
			else
				! higher ranks first receive and then send...
				if(nrecv(proc) >0)then
					tag_f = TAG_START + lcsrank
					call MPI_RECV(recvbuf(r_start),nrecv(proc)*LPMAP_SIZE,MPI_LCSRP,proc,tag_f,lcscomm,status,ierr)
				endif
				if(nsend(proc) >0)then
					tag_f = TAG_START + proc
					call MPI_SEND(sendbuf(s_start),nsend(proc)*LPMAP_SIZE,MPI_LCSRP,proc,tag_f,lcscomm,ierr)
				endif
			endif
		end do

		!-----
		!Unpack the data
		!Check that each IB grid point recieves exactly 1 data point.
		!-----
		ibuf = 1
		visited = 0
		do ip = 1,nrecv_total
			no0  = nint(recvbuf(ibuf+0))
			i = l2i(no0,lp%fm%ni)
			j = l2j(no0,lp%fm%ni,lp%fm%nj)
			k = l2k(no0,lp%fm%ni,lp%fm%nj)
			lp%fm%x(i,j,k) = recvbuf(ibuf+1)
			lp%fm%y(i,j,k) = recvbuf(ibuf+2)
			lp%fm%z(i,j,k) = recvbuf(ibuf+3)
			ibuf = ibuf + LPMAP_SIZE
			visited(i,j,k) = visited(i,j,k)+1
		enddo
		do k =1,lp%fm%nk
		do j =1,lp%fm%nj
		do i =1,lp%fm%ni
			if(visited(i,j,k)/=1) then
				write(*,*) 'lcsrank[',lcsrank,'] ERROR:  i,j,k, visited',visited(i,j,k),'times'
				CFD2LCS_ERROR = 1
			endif
		enddo
		enddo
		enddo

		t1 = cputimer(lcscomm,SYNC_TIMER)
		cpu_lpmap = cpu_lpmap + (t1-t0)

	end subroutine exchange_lpmap

	subroutine exchange_lp_alltoall(lp,sgrid)
		use lp_m
		implicit none
		!-----
		type(lp_t):: lp
		type(sgrid_t):: sgrid
		!-----
		integer:: ni, nj, nk, ng, ip, proc, is, ir
		integer:: i,j,k,ierr
		integer,allocatable:: nsend(:),nrecv(:)
		integer:: nsend_global,nrecv_global, np_pack_total, np_unpack_total
		integer:: tag_f  !unique tag for each send/rec pair
		integer :: status(MPI_STATUS_SIZE)
		integer,allocatable:: sendstart(:), recvstart(:), tmp(:)
		real(LCSRP),allocatable:: sendbuf(:),recvbuf(:)
		integer,parameter:: ALLTOALL_SIZE = 8  !xp,yp,zp,fmx,fmy,fmz,real(node0),real(proc0)
		integer:: r_start,r_end, s_start,s_end
		integer:: ibuf,iunpack
		real(LCSRP):: ijk_xmin(3),ijk_xmax(3),ijk_ymin(3),ijk_ymax(3),ijk_zmin(3),ijk_zmax(3)
		real(LCSRP):: s,t,u
		real(LCSRP),parameter:: onethird = 1.0_LCSRP/3.0_LCSRP
		type(ui0_t):: found
		!-----
		!Used for unstructured exchange or flow map reconstruction
		!-----

		if(lcsrank==0 .AND. LCS_VERBOSE) &
			write(*,*) 'In exchange_lp_alltoall...', trim(lp%label),trim(sgrid%label)

		!brevity...
		ni = sgrid%ni; nj = sgrid%nj; nk = sgrid%nk; ng = sgrid%ng

		!Make sure we do a recursive search during construction
		lp%recursive_tracking = .true.

		!Now create a list of particles to send to each processor
		allocate(nsend(0:nprocs-1))
		allocate(nrecv(0:nprocs-1))
		allocate(tmp(0:nprocs-1))
		call init_ui0(found,lp%np,'found')
		nsend(:) = 0
		nrecv(:) = 0
		found%i(:) = LP_UNKNOWN
		iploop: do ip = 1,lp%np

			!Send to the first proc which has the particle in it's
			!minimum bounding box:
			bbloop: do proc = 0,nprocs-1
				if (&
					lp%xp%x(ip) >= sgrid%bb(proc,1) .and. &
					lp%xp%x(ip) <= sgrid%bb(proc,4) .and. &
					lp%xp%y(ip) >= sgrid%bb(proc,2) .and. &
					lp%xp%y(ip) <= sgrid%bb(proc,5) .and. &
					lp%xp%z(ip) >= sgrid%bb(proc,3) .and. &
					lp%xp%z(ip) <= sgrid%bb(proc,6) &
				) then
					found%i(ip) = proc
					nsend(proc) = nsend(proc)+1
					!exit
					cycle iploop
				endif
			enddo bbloop

			!If you didnt find it, check the (first) periodic image domains
			!for each processor.
			if(found%i(ip)==LP_UNKNOWN) then
				imageloop: do proc = 0,nprocs-1
					do k=-1,1
					do j=-1,1
					do i=-1,1
						if (&
							lp%xp%x(ip) - sgrid%ps(proc,i,j,k,1) >= sgrid%bb(proc,1) .and. &
							lp%xp%x(ip) - sgrid%ps(proc,i,j,k,1) <= sgrid%bb(proc,4) .and. &
							lp%xp%y(ip) - sgrid%ps(proc,i,j,k,2) >= sgrid%bb(proc,2) .and. &
							lp%xp%y(ip) - sgrid%ps(proc,i,j,k,2) <= sgrid%bb(proc,5) .and. &
							lp%xp%z(ip) - sgrid%ps(proc,i,j,k,3) >= sgrid%bb(proc,3) .and. &
							lp%xp%z(ip) - sgrid%ps(proc,i,j,k,3) <= sgrid%bb(proc,6) &
						) then
							lp%xp%x(ip) = lp%xp%x(ip) - sgrid%ps(proc,i,j,k,1)
							lp%xp%y(ip) = lp%xp%y(ip) - sgrid%ps(proc,i,j,k,2)
							lp%xp%z(ip) = lp%xp%z(ip) - sgrid%ps(proc,i,j,k,3)
							found%i(ip) = proc
							nsend(proc) = nsend(proc)+1
							!exit imageloop
							cycle iploop
						endif
					enddo
					enddo
					enddo
				enddo imageloop
			endif

			!Default procedure is to keep the particle on current proc
			!and brute force search if needed.
			if(found%i(ip)==LP_UNKNOWN) then
				found%i(ip) = lcsrank
				nsend(lcsrank) = nsend(lcsrank)+1
			endif
		enddo iploop

		!-----
		!Now communicate to receiving processors the number of send/recv
		!-----
		if(ALLTOALL_REDUCE) then
			!Use several reduce ops:
			do proc = 0, nprocs-1
				tmp = 0
				tmp(lcsrank) = nsend(proc)
				call MPI_REDUCE(tmp,nrecv,nprocs,MPI_INTEGER,MPI_SUM,proc,lcscomm,ierr)
			enddo
		else
			do proc = 0, nprocs-1
				if(proc == lcsrank) then
					nrecv(lcsrank) = nsend(lcsrank)
				elseif (lcsrank < proc) then
					! lower ranks first send and then receive...
					tag_f = TAG_START + proc
					call MPI_SEND(nsend(proc),1,MPI_INTEGER,proc,tag_f,lcscomm,ierr)
					tag_f = TAG_START + lcsrank
					call MPI_RECV(nrecv(proc),1,MPI_INTEGER,proc,tag_f,lcscomm,status,ierr)
				else
					! higher ranks first receive and then send...
					tag_f = TAG_START + lcsrank
					call MPI_RECV(nrecv(proc),1,MPI_INTEGER,proc,tag_f,lcscomm,status,ierr)
					tag_f = TAG_START + proc
					call MPI_SEND(nsend(proc),1,MPI_INTEGER,proc,tag_f,lcscomm,ierr)
				endif
			end do
		endif
		np_pack_total = sum(nsend)
		np_unpack_total = sum(nrecv)
		call MPI_REDUCE(np_pack_total,nsend_global,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)
		call MPI_REDUCE(np_unpack_total,nrecv_global,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)
		if(lcsrank==0) then
			if(nrecv_global /= nsend_global) then
				 write(*,*) 'ERROR: global nsend /= nrecv'
				CFD2LCS_ERROR = 1
			else
				if(LCS_VERBOSE) &
				 write(*,*) 'Exchange of ',nsend_global,' flowmap particles'
			endif
		endif

		!-----
		!Determine where the send/recieve will start for each rank
		!This will allow for easy indexing into a single vector when communicating.
		!-----
		allocate(sendstart(0:nprocs-1))
		allocate(recvstart(0:nprocs-1))
		sendstart(:) = -1
		recvstart(:) = -1
		is = 0; ir = 0
		do proc = 0, nprocs-1
			if (nsend(proc) > 0) then
				sendstart(proc) = is+1
				is = is + nsend(proc)*ALLTOALL_SIZE
			endif
			if (nrecv(proc) > 0) then
				recvstart(proc) = ir+1
				ir = ir + nrecv(proc)*ALLTOALL_SIZE
			endif
		enddo

		!check
		if (is /= np_pack_total*ALLTOALL_SIZE) then
			write(*,*) '[',lcsrank,'] ERROR, bad send indexing'
			CFD2LCS_ERROR = 1
		endif
		if (ir /= np_unpack_total*ALLTOALL_SIZE) then
			write(*,*) '[',lcsrank,'] ERROR, bad recv indexing'
			CFD2LCS_ERROR = 1
		endif

		!-----
		!Pack the buffers
		!Note, we send dx and not xp to handle periodic domains.
		!-----
		allocate(sendbuf(1:np_pack_total*ALLTOALL_SIZE))
		allocate(recvbuf(1:np_unpack_total*ALLTOALL_SIZE))
		tmp = sendstart
		do ip = 1,lp%np
			proc = found%i(ip)

			sendbuf(sendstart(proc)+0) = lp%xp%x(ip)
			sendbuf(sendstart(proc)+1) = lp%xp%y(ip)
			sendbuf(sendstart(proc)+2) = lp%xp%z(ip)
			sendbuf(sendstart(proc)+3) = lp%dx%x(ip)
			sendbuf(sendstart(proc)+4) = lp%dx%y(ip)
			sendbuf(sendstart(proc)+5) = lp%dx%z(ip)
			sendbuf(sendstart(proc)+6) = real(lp%no0%i(ip),LCSRP)
			sendbuf(sendstart(proc)+7) = real(lp%proc0%i(ip),LCSRP)
			sendstart(proc) = sendstart(proc) + ALLTOALL_SIZE
			lp%flag%i(ip) = LP_RECYCLE
		enddo
		sendstart = tmp



		!-----
		!Exchange buffers
		!-----
		do proc = 0, nprocs-1
			s_start = sendstart(proc)
			s_end = s_start + nsend(proc) * ALLTOALL_SIZE-1
			r_start = recvstart(proc)
			r_end = r_start + nrecv(proc) * ALLTOALL_SIZE-1
			if(proc == lcsrank .AND. nsend(proc)>0) then
				recvbuf(r_start:r_end) = sendbuf(s_start:s_end)
			elseif (lcsrank < proc) then
				! lower ranks first send and then receive...
				if(nsend(proc) >0)then
					tag_f = TAG_START + proc
					call MPI_SEND(sendbuf(s_start),nsend(proc)*ALLTOALL_SIZE,MPI_LCSRP,proc,tag_f,lcscomm,ierr)
				endif
				if(nrecv(proc) >0)then
					tag_f = TAG_START + lcsrank
					call MPI_RECV(recvbuf(r_start),nrecv(proc)*ALLTOALL_SIZE,MPI_LCSRP,proc,tag_f,lcscomm,status,ierr)
				endif
			else
				! higher ranks first receive and then send...
				if(nrecv(proc) >0)then
					tag_f = TAG_START + lcsrank
					call MPI_RECV(recvbuf(r_start),nrecv(proc)*ALLTOALL_SIZE,MPI_LCSRP,proc,tag_f,lcscomm,status,ierr)
				endif
				if(nsend(proc) >0)then
					tag_f = TAG_START + proc
					call MPI_SEND(sendbuf(s_start),nsend(proc)*ALLTOALL_SIZE,MPI_LCSRP,proc,tag_f,lcscomm,ierr)
				endif
			endif
		end do

		!-----
		!Reorder and remove holes in the lp list
		!-----
		if(np_pack_total > 0) then
			call reorder_lp(lp)
		endif

		!-----
		!Resize lp
		!-----
		ip = lp%np  !starting point for unpacking
		call resize_lp(lp,lp%np+np_unpack_total)

		!-----
		!Unpack data
		!-----
		if(.NOT. sgrid%rectilinear) then
			ijk_xmin = real(minloc(sgrid%grid%x),LCSRP)
			ijk_ymin = real(minloc(sgrid%grid%y),LCSRP)
			ijk_zmin = real(minloc(sgrid%grid%z),LCSRP)
			ijk_xmax = real(maxloc(sgrid%grid%x),LCSRP)
			ijk_ymax = real(maxloc(sgrid%grid%y),LCSRP)
			ijk_zmax = real(maxloc(sgrid%grid%z),LCSRP)
		endif
		ibuf = 1
		do iunpack = 1,np_unpack_total
			ip = ip + 1
			lp%xp%x(ip) = recvbuf(ibuf+0)
			lp%xp%y(ip) = recvbuf(ibuf+1)
			lp%xp%z(ip) = recvbuf(ibuf+2)
			lp%dx%x(ip) = recvbuf(ibuf+3)
			lp%dx%y(ip) = recvbuf(ibuf+4)
			lp%dx%z(ip) = recvbuf(ibuf+5)
			lp%no0%i(ip)= nint(recvbuf(ibuf+6))
			lp%proc0%i(ip)= nint(recvbuf(ibuf+7))
			lp%up%x(ip) = 0.0_LCSRP
			lp%up%y(ip) = 0.0_LCSRP
			lp%up%z(ip) = 0.0_LCSRP
			!Guess the node that this particle belongs to:
			!This should work well for rectilinear grids.  For non-uniform, you may want something
			!more advanced (like oct-tree search perhaps).
			if(sgrid%rectilinear) then
				s =max(min((lp%xp%x(ip)-sgrid%bb(lcsrank,1))&
					/(sgrid%bb(lcsrank,4)-sgrid%bb(lcsrank,1)),1.0_LCSRP),0.0_LCSRP)
				t =max(min((lp%xp%y(ip)-sgrid%bb(lcsrank,2))&
					/(sgrid%bb(lcsrank,5)-sgrid%bb(lcsrank,2)),1.0_LCSRP),0.0_LCSRP)
				u =max(min((lp%xp%z(ip)-sgrid%bb(lcsrank,3))&
					/(sgrid%bb(lcsrank,6)-sgrid%bb(lcsrank,3)),1.0_LCSRP),0.0_LCSRP)
				lp%no%x(ip) = nint(real(ni)*s)
				lp%no%y(ip) = nint(real(nj)*t)
				lp%no%z(ip) = nint(real(nk)*u)
			else
				s =max(min((lp%xp%x(ip)-sgrid%bb(lcsrank,1))&
					/(sgrid%bb(lcsrank,4)-sgrid%bb(lcsrank,1)),1.0_LCSRP),0.0_LCSRP)
				t =max(min((lp%xp%y(ip)-sgrid%bb(lcsrank,2))&
					/(sgrid%bb(lcsrank,5)-sgrid%bb(lcsrank,2)),1.0_LCSRP),0.0_LCSRP)
				u =max(min((lp%xp%z(ip)-sgrid%bb(lcsrank,3))&
					/(sgrid%bb(lcsrank,6)-sgrid%bb(lcsrank,3)),1.0_LCSRP),0.0_LCSRP)

				lp%no%x(ip) = nint(onethird*(ijk_xmin(1) + ijk_ymin(1) + ijk_zmin(1))) &
					+ nint(onethird*(s*(ijk_xmax(1)-ijk_xmin(1)) &
					+t*(ijk_ymax(1)-ijk_ymin(1)) &
					+u*(ijk_zmax(1)-ijk_zmin(1))))
				lp%no%y(ip) = nint(onethird*(ijk_xmin(2) + ijk_ymin(2) + ijk_zmin(2))) &
					+ nint(onethird*(s*(ijk_xmax(2)-ijk_xmin(2)) &
					+t*(ijk_ymax(2)-ijk_ymin(2)) &
					+u*(ijk_zmax(2)-ijk_zmin(2))))
				lp%no%z(ip) = nint(onethird*(ijk_xmin(3) + ijk_ymin(3) + ijk_zmin(3))) &
					+ nint(onethird*(s*(ijk_xmax(3)-ijk_xmin(3)) &
					+t*(ijk_ymax(3)-ijk_ymin(3)) &
					+u*(ijk_zmax(3)-ijk_zmin(3))))

				!A sort of brute force approach...
				!lp%no%x(ip) = max(ni/2,1)
				!lp%no%y(ip) = max(nj/2,1)
				!lp%no%z(ip) = max(nk/2,1)
			endif

			!Limit guess to ib nodes:
			lp%no%x(ip) = max(min(lp%no%x(ip),ni),1)
			lp%no%y(ip) = max(min(lp%no%y(ip),nj),1)
			lp%no%z(ip) = max(min(lp%no%z(ip),nk),1)

			!Set the flag
			lp%flag%i(ip) = LP_UNKNOWN
			ibuf = ibuf + ALLTOALL_SIZE
		enddo


		!cleanup
		call destroy_ui0(found)
		deallocate(nsend)
		deallocate(nrecv)
		deallocate(sendstart)
		deallocate(recvstart)
		deallocate(tmp)
		deallocate(sendbuf)
		deallocate(recvbuf)


	end subroutine exchange_lp_alltoall

end module comms_m
