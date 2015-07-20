!JRF:  A simple program to test the cfd2lcs interface
!Set all global parameters in user_data_m.
module user_data_m
	implicit none

	!-----
	!Working Precision (4 = single, 8 = double)
	!Make sure this matches the precision you compile the library with
	!-----
	integer,parameter:: WP = 4

	!-----
	!Domain dimensions
	!-----
	real(WP), parameter:: LX = 2.0
	real(WP), parameter:: LY = 1.0
	real(WP), parameter:: LZ = 1.0

	!-----
	!Total number of grid pointimestep in each direction
	!-----
	integer, parameter:: NX = 64
	integer, parameter:: NY = 32
	integer, parameter:: NZ = 32

	!-----
	!Number of processors in each direction
	!-----
	integer, parameter:: NPROC_X = 4
	integer, parameter:: NPROC_Y = 1
	integer, parameter:: NPROC_Z = 1

	!-----
	!Possible boundary conditions for the domain exterior:
	!Will generalize to arbitrary patches later on...
	!-----
	integer,parameter:: &
		PERIODIC = 0, &
		WALL = 1, &
		INFLOW = 2, &
		OUTFLOW = 3

	integer,parameter:: BC_X0 = PERIODIC
	integer,parameter:: BC_X1 = PERIODIC
	integer,parameter:: BC_Y0 = PERIODIC
	integer,parameter:: BC_Y1 = PERIODIC
	integer,parameter:: BC_Z0 = PERIODIC
	integer,parameter:: BC_Z1 = PERIODIC

	!-----
	!"Simulation" parameters
	!-----
	real(WP),parameter:: DT = 1e-2_WP
	real(WP),parameter:: START_TIME = 0.0_WP
	real(WP),parameter:: END_TIME = 0.001_WP



end module user_data_m


program interface_test
	use user_data_m
	implicit none
	!-----
	INCLUDE 'mpif.h'
	integer:: mycomm,nprocs,ierr
	integer:: myrank, myrank_i, myrank_j, myrank_k
	integer:: MASTER = 0
	integer:: ni, nj, nk, offset_i, offset_j, offset_k
	integer:: timestep
	real(WP):: time
	real(WP), allocatable:: x(:,:,:), y(:,:,:), z(:,:,:)
	real(WP), allocatable:: u(:,:,:), v(:,:,:), w(:,:,:)
	!-----

	!-----
	!Initialize MPI:
	!-----
	call MPI_INIT(ierr)
	mycomm = MPI_COMM_WORLD
	call MPI_COMM_RANK(mycomm,myrank,ierr)
	call MPI_COMM_SIZE(mycomm,nprocs,ierr)
	call MPI_BCAST(MASTER,1,MPI_INTEGER,0,mycomm,ierr)

	if (myrank == MASTER) then
	write(*,'(a)') '******************************'
	write(*,'(a)') 'libcfd2lcs interface test'
	write(*,'(a)') '******************************'
	endif

	!-----
	!Setup the parallel partition:
	!-----
	call set_partition()

	!-----
	!Set the grid pointimestep:
	!-----
	call set_grid()

	!-----
	!Now we initialize cfd2lcs
	!-----
	call cfd2lcs_init(mycomm,ni,nj,nk,offset_i,offset_j,offset_k,x,y,z)



	!-----
	!***Start of Main Flow solver loop***
	!-----
	timestep = 0
	time = START_TIME
	call set_velocity(time)
	do while (time <= END_TIME)
		timestep = timestep + 1
		if(myrank == MASTER) then
			write(*,'(a)') '------------------------------------------------------------------'
			write(*,'(a,i10.0,a,ES10.4,a,ES10.4)') 'STARTING TIMESTEP #',timestep,': time = ',time,', DT = ',DT
			write(*,'(a)') '------------------------------------------------------------------'
		endif


		time = time + DT
		call set_velocity(time) !CFD Solve for the velocity field

		call cfd2lcs_update(ni,nj,nk,u,v,w,time)  !Update LCS fields
	enddo


	!-----
	!Cleanup
	!-----
	deallocate(x)
	deallocate(y)
	deallocate(z)
	deallocate(u)
	deallocate(v)
	deallocate(w)
	call cfd2lcs_finalize(ierr)
	call MPI_FINALIZE(ierr)


	contains

	subroutine set_partition()
		implicit none
		!----
		integer:: proc, guess_ni, guess_nj, guess_nk
		integer:: leftover_ni, leftover_nj, leftover_nk
		!----
		!Partition the domain into chunks for each processor
		!based on the number of processors specified in X,Y,Z direction.  Unequal sizes are permitted.
		!Actual grid coordinates are set in set_grid().
		!----

		if (myrank ==MASTER)&
			write(*,'(a)') 'in set_partition...'

		!Check the number of procs:
		if (nprocs /= NPROC_X*NPROC_Y*NPROC_Z) then
			if(myrank==MASTER) write(*,'(a)') 'ERROR: Number of processors incompatible with partition'
			call MPI_FINALIZE()
			stop
		endif

		!Set local proc coordinates:
		myrank_i = rank2i(myrank,NPROC_X)
		myrank_j = rank2j(myrank,NPROC_X,NPROC_Y)
		myrank_k = rank2k(myrank,NPROC_X,NPROC_Y)

		!Allocate the array range on each processor, allowing for unequal sizes:
		guess_ni = NX/NPROC_X
		guess_nj = NY/NPROC_Y
		guess_nk = NZ/NPROC_Z
		leftover_ni = NX - (guess_ni*NPROC_X)
		leftover_nj = NY - (guess_nj*NPROC_Y)
		leftover_nk = NZ - (guess_nk*NPROC_Z)
		if (myrank_i.LE.leftover_ni) then
			ni = guess_ni + 1
			offset_i = (guess_ni+1)*(myrank_i-1)
		else
			ni = guess_ni
			offset_i = (guess_ni+1)*(leftover_ni) + guess_ni*(myrank_i-leftover_ni-1)
		end if

		if (myrank_j.LE.leftover_nj) then
			nj = guess_nj + 1
			offset_j = (guess_nj+1)*(myrank_j-1)
		else
			nj = guess_nj
			offset_j = (guess_nj+1)*(leftover_nj) + guess_nj*(myrank_j-leftover_nj-1)
		end if

		if (myrank_k.LE.leftover_nk) then
			nk = guess_nk + 1
			offset_k = (guess_nk+1)*(myrank_k-1)
		else
			nk = guess_nk
			offset_k = (guess_nk+1)*(leftover_nk) + guess_nk*(myrank_k-leftover_nk-1)
		end if

		!----
		!Check...
		!----
		do proc = 0, nprocs-1
			if(myrank == proc) then
				write(*,'(a,i6,a,i4,i4,i4)') 'Proc [',myrank,'] has subdomain',myrank_i,myrank_j,myrank_k
				write(*,'(a,i6,a,i4,a,i4)') 'Proc [',myrank,'] has "X" nodes', offset_i+1 ,' - ',offset_i+ni
				write(*,'(a,i6,a,i4,a,i4)') 'Proc [',myrank,'] has "Y" nodes', offset_j+1 ,' - ',offset_j+nj
				write(*,'(a,i6,a,i4,a,i4)') 'Proc [',myrank,'] has "Z" nodes', offset_k+1 ,' - ',offset_k+nk
			endif
			call MPI_BARRIER(mycomm,ierr)
		enddo

	end subroutine set_partition


	subroutine set_grid()
		implicit none
		!----
		integer:: i,j,k,ii,jj,kk
		!----
		!Here set the grid coordinates x,y,z.
		!Just use a uniform Cartesian grid for now, arbitrary spacing is possible.
		!----
		if (myrank ==MASTER)&
			write(*,'(a)') 'in set_grid...'

		allocate(x(1:ni,1:nj,1:nk))
		allocate(y(1:ni,1:nj,1:nk))
		allocate(z(1:ni,1:nj,1:nk))

		kk = 0
		do k = offset_k+1,offset_k+nk
			kk = kk + 1
			jj = 0
			do j = offset_j+1,offset_j+nj
				jj = jj + 1
				ii = 0
				do i = offset_i+1,offset_i+ni
					ii = ii + 1
					x(ii,jj,kk) = real(i-1)*real(LX)/real(NX-1)
					y(ii,jj,kk) = real(j-1)*real(LY)/real(NY-1)
					z(ii,jj,kk) = real(k-1)*real(LZ)/real(NZ-1)
				enddo
			enddo
		enddo

	end subroutine set_grid


	subroutine set_velocity(time)
		implicit none
		real(WP),intent(in):: time
		!----
		integer:: i,j,k
		!----

		if (myrank ==MASTER)&
			write(*,'(a)') 'in set_velocity...'

		if (.NOT. allocated(u)) allocate(u(ni,nj,nk))
		if (.NOT. allocated(v)) allocate(v(ni,nj,nk))
		if (.NOT. allocated(w)) allocate(w(ni,nj,nk))


		u = real(myrank,WP)
		v = 2.0_WP* real(myrank,WP)
		w = real(-myrank,WP)

	end subroutine set_velocity

	!----
	!These function set the rules for the relationship
	!between 3D (i,j,k) indices and 1D (l) index. Note, 1D rank, is assumed to start at 0
	!----
	integer function ijk2rank(i,j,k,npi,npj)
		implicit none
		integer::i,j,k,npi,npj
		ijk2rank= (k-1)*npi*npj + (j-1)*npi +i -1
	end function ijk2rank
	integer function rank2i(rank,npi)
		implicit none
		integer::rank,npi
	    rank2i= mod(rank,npi)+1
	end function rank2i
	integer function rank2j(rank,npi,npj)
		implicit none
		integer::rank,npi,npj
		rank2j = mod((rank)/npi,npj)+1
	end function rank2j
	integer function rank2k(rank,npi,npj)
		implicit none
		integer::rank,npi,npj
		rank2k = (rank)/(npi*npj)+1
	end function rank2k


end program interface_test
