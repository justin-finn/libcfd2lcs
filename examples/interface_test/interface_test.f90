!
! A simple program to test the cfd2lcs interface
!
program interface_test
	use user_data_m
	use analytic_velocity_m
	implicit none
	!-----
	INCLUDE 'mpif.h'
	integer narg
	character(len=32):: arg
	integer:: mycomm,nprocs,ierr
	integer:: nproc_x,nproc_y,nproc_z
	integer:: myrank, myrank_i, myrank_j, myrank_k
	integer:: ni, nj, nk, offset_i, offset_j, offset_k
	integer:: timestep
	real(LCSRP):: time
	real(LCSRP), allocatable:: x(:,:,:), y(:,:,:), z(:,:,:)
	real(LCSRP), allocatable:: u(:,:,:), v(:,:,:), w(:,:,:)
	integer:: n(3),offset(3)
	real(LCSRP)::lperiodic(3)
	integer(LCSIP):: BC_LIST(6)
	!-----


	!-----
	!Initialize MPI:
	!-----
	call MPI_INIT(ierr)
	mycomm = MPI_COMM_WORLD
	call MPI_COMM_RANK(mycomm,myrank,ierr)
	call MPI_COMM_SIZE(mycomm,nprocs,ierr)

	if (myrank == 0) then
	write(*,'(a)') '******************************'
	write(*,'(a)') 'libcfd2lcs interface test'
	write(*,'(a)') '******************************'
	endif


	!-----
	!Parse the input
	!-----
	narg=command_argument_count()
	if(narg/=3)then
		if (myrank== 0) then
			write(*,*) 'Error: must supply arguments for NPROCS_X, NPROCS_Y, NPROCS_Z'
			write(*,*) '	Example usage:  mpirun -np 8 CFD2LCS_TEST 4 2 1'
		endif
		call mpi_barrier(mycomm,ierr)
		stop
	endif
	call getarg(1,arg)
	read(arg,*) nproc_x
	call getarg(2,arg)
	read(arg,*) nproc_y
	call getarg(3,arg)
	read(arg,*) nproc_z

	!-----
	!Setup the parallel partition:
	!-----
	call set_partition()

	!-----
	!Set the grid points:
	!-----
	call set_grid()

	!-----
	!Now we initialize cfd2lcs
	!-----
	n = (/ni,nj,nk/)  !number of grid points for THIS partition
	offset = (/offset_i,offset_j,offset_k/)  !Global offset of these grid points
	lperiodic = (/LX,LY,LZ/)  !Periodic length of the domain in x,y,z
	BC_LIST = (/BC_IMIN,BC_JMIN,BC_KMIN,BC_IMAX,BC_JMAX,BC_KMAX/) !List of boundary conditions
	call cfd2lcs_init(mycomm,n,offset,x,y,z,BC_LIST,lperiodic)

	!-----
	!***Start of Main Flow solver loop***
	!-----
	timestep = 0
	time = START_TIME
	call set_velocity(time)
	do while (time <= END_TIME)
		timestep = timestep + 1
		if(myrank == 0) then
			write(*,'(a)') '------------------------------------------------------------------'
			write(*,'(a,i10.0,a,ES10.4,a,ES10.4)') 'STARTING TIMESTEP #',timestep,': time = ',time,', DT = ',DT
			write(*,'(a)') '------------------------------------------------------------------'
		endif

		time = time + DT
		call set_velocity(time) !CFD Solve for the velocity field

		call cfd2lcs_update(n,u,v,w,time)  !Update LCS fields
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

		if (myrank ==0)&
			write(*,'(a)') 'in set_partition...'

		!Check the number of procs:
		if (nprocs /= nproc_x*nproc_y*nproc_z) then
			if(myrank==0) write(*,'(a)') 'ERROR: Number of processors incompatible with partition'
			call MPI_FINALIZE()
			stop
		endif

		if (nproc_x > NX .OR. nproc_y > NY .OR. nproc_z > NZ) then
			if(myrank==0) write(*,'(a)') 'ERROR: More processors than grid points...'
			call MPI_FINALIZE()
			stop
		endif

		!Set local proc coordinates:
		myrank_i = rank2i(myrank,nproc_x)
		myrank_j = rank2j(myrank,nproc_x,nproc_y)
		myrank_k = rank2k(myrank,nproc_x,nproc_y)

		!Allocate the array range on each processor, allowing for unequal sizes:
		guess_ni = NX/nproc_x
		guess_nj = NY/nproc_y
		guess_nk = NZ/nproc_z
		leftover_ni = NX - (guess_ni*nproc_x)
		leftover_nj = NY - (guess_nj*nproc_y)
		leftover_nk = NZ - (guess_nk*nproc_z)
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

	end subroutine set_partition


	subroutine set_grid()
		implicit none
		!----
		integer:: i,j,k,ii,jj,kk
		real(LCSRP):: dx,dy,dz
		!----
		!Here set the grid coordinates x,y,z.
		!Just use a uniform Cartesian grid for now, arbitrary spacing is possible.
		!----
		if (myrank ==0)&
			write(*,'(a)') 'in set_grid...'

		allocate(x(1:ni,1:nj,1:nk))
		allocate(y(1:ni,1:nj,1:nk))
		allocate(z(1:ni,1:nj,1:nk))

		dx = LX / real(NX)
		dy = LY / real(NY)
		dz = LZ / real(NZ)

		kk = 0
		do k = offset_k+1,offset_k+nk
			kk = kk + 1
			jj = 0
			do j = offset_j+1,offset_j+nj
				jj = jj + 1
				ii = 0
				do i = offset_i+1,offset_i+ni
					ii = ii + 1
					x(ii,jj,kk) = 0.5*dx + real(i-1)*dx
					y(ii,jj,kk) = 0.5*dy + real(j-1)*dy
					z(ii,jj,kk) = 0.5*dx + real(k-1)*dz
				enddo
			enddo
		enddo

	end subroutine set_grid


	subroutine set_velocity(time)
		implicit none
		real(LCSRP),intent(in):: time
		!----
		integer:: i,j,k
		real(LCSRP),parameter:: ABC_A = sqrt(3.0_LCSRP)
		real(LCSRP),parameter:: ABC_B = sqrt(2.0_LCSRP)
		real(LCSRP),parameter:: ABC_C = 1.0_LCSRP
		real(LCSRP),parameter:: ABC_D = 0.0_LCSRP
		!----

		if (myrank ==0)&
			write(*,'(a)') 'in set_velocity...'

		if (.NOT. allocated(u)) allocate(u(ni,nj,nk))
		if (.NOT. allocated(v)) allocate(v(ni,nj,nk))
		if (.NOT. allocated(w)) allocate(w(ni,nj,nk))

		call abc_velocity(ni,nj,nk,x,y,z,u,v,w,ABC_A, ABC_B, ABC_C, ABC_D,time)

	end subroutine set_velocity

	!----
	!These function set the rules for the relationship
	!between 3D (i,j,k) indices of each processor and the processor rank.
	!----
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
