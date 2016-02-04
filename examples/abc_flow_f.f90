!
! A Simple working example to show how to call 
! CFD2LCS.  Subroutines starting with "your_" 
! are independent of the cfd2lcs functionality, 
! and are used only to create a simple dataset
! for this example.
!
!
! FLOW FIELD:  3D, time dependent "ABC Flow"
! User parameters included directly below.
!
!
program abc_flow 
	implicit none
	!-----
	include 'cfd2lcs_inc.f90'  !This includes parameter definitions needed for cfd2lcs
	include 'mpif.h'
	!******BEGIN USER INPUT********************
	!-----
	!Domain dimensions
	!-----
	real(LCSRP), parameter:: PI = 4.0*atan(1.0)
	real(LCSRP), parameter:: LX = 2.0*PI
	real(LCSRP), parameter:: LY = 2.0*PI
	real(LCSRP), parameter:: LZ = 2.0*PI
	!-----
	!Total number of grid points in each direction
	!-----
	integer, parameter:: NX =  128 
	integer, parameter:: NY =  128
	integer, parameter:: NZ =  128
	!-----
	!Boundary conditions for the domain exterior:
	!-----
	integer(LCSIP),parameter:: BC_IMIN = LCS_PERIODIC
	integer(LCSIP),parameter:: BC_JMIN = LCS_PERIODIC
	integer(LCSIP),parameter:: BC_KMIN = LCS_PERIODIC
	integer(LCSIP),parameter:: BC_IMAX = LCS_PERIODIC
	integer(LCSIP),parameter:: BC_JMAX = LCS_PERIODIC
	integer(LCSIP),parameter:: BC_KMAX = LCS_PERIODIC
	!-----
	!"Simulation" parameters
	!-----
	real(LCSRP),parameter:: DT = 0.01
	real(LCSRP),parameter:: START_TIME = 0.0
	real(LCSRP),parameter:: END_TIME = 10.1
	real(LCSRP),parameter:: ABC_A = sqrt(3.0)
	real(LCSRP),parameter:: ABC_B = sqrt(2.0)
	real(LCSRP),parameter:: ABC_C = 1.0
	real(LCSRP),parameter:: ABC_D = 0.0
	!----
	!******END USER INPUT************************
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
	integer(LCSIP):: id_fwd,id_bkwd,id_tracer
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
	write(*,'(a)') 'libcfd2lcs double gyre test'
	write(*,'(a)') '******************************'
	endif

	!-----
	!Parse the input on rank 0 and broadcast.
	!Assume that we want the last 3 arguments for NPX, NPY, NPZ,
	!in case there are hidden args picked up by the mpirun or other process.
	!-----
	if(myrank==0) then
		narg=command_argument_count()
		if(narg<3)then
			write(*,*) 'Error: must supply arguments for NPROCS_X, NPROCS_Y, NPROCS_Z'
			write(*,*) '	Example usage:  mpirun -np 8 ./CFD2LCS_TEST 4 2 1'
			stop
		endif
		call getarg(narg-2,arg)
		read(arg,*) nproc_x
		call getarg(narg-1,arg)
		read(arg,*) nproc_y
		call getarg(narg-0,arg)
		read(arg,*) nproc_z
		write(*,*) 'Will partition domain using:',nproc_x,nproc_y,nproc_z,' sub-domains'
	endif
	call MPI_BCAST(nproc_x,1,MPI_INTEGER,0,mycomm,ierr)
	call MPI_BCAST(nproc_y,1,MPI_INTEGER,0,mycomm,ierr)
	call MPI_BCAST(nproc_z,1,MPI_INTEGER,0,mycomm,ierr)

	!-----
	!Setup the parallel partition:
	!-----
	call your_partition_function()

	!-----
	!Set the grid points:
	!-----
	call your_grid_function()

	!-----
	!Initialize cfd2lcs for your data
	!-----
	n = (/ni,nj,nk/)  !number of grid points for THIS partition
	offset = (/offset_i,offset_j,offset_k/)  !Global offset of these grid points
	lperiodic = (/LX,LY,LZ/)  !Periodic length of the domain in x,y,z
	BC_LIST = (/BC_IMIN,BC_JMIN,BC_KMIN,BC_IMAX,BC_JMAX,BC_KMAX/) !List of boundary conditions
	call cfd2lcs_init(mycomm,n,offset,x,y,z,BC_LIST,lperiodic)

	!-----
	!Initialize LCS diagnostics
	!-----
	call cfd2lcs_diagnostic_init(id_fwd,FTLE_FWD,0,15.0,1.0,0.0,0.0,'fwdFTLE')
	call cfd2lcs_diagnostic_init(id_bkwd,FTLE_BKWD,0,15.0,1.0,0.0,0.0,'bkwdFTLE')
	!call cfd2lcs_diagnostic_init(id_tracer,LP_TRACER,0,15.0,1.0,0.0,0.0,'Tracers')

	!-----
	!***Start of your flow solver timestepping loop***
	!-----
	timestep = 0
	time = START_TIME
	call your_flow_solver(time)
	do while (time <= END_TIME)
		if(myrank == 0) then
			write(*,'(a)') '------------------------------------------------------------------'
			write(*,'(a,i10.0,a,ES11.4,a,ES11.4)') 'STARTING TIMESTEP #',timestep,': time = ',time,', DT = ',DT
			write(*,'(a)') '------------------------------------------------------------------'
		endif
		
		!Produce the new velocity field with your flow solver
		call your_flow_solver(time) 
		
		!Update the LCS diagnostics using the new flow field
		call cfd2lcs_update(n,u,v,w,time)  
		
		timestep = timestep + 1
		time = time + DT
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

	subroutine your_partition_function()
		implicit none
		!----
		integer:: proc, guess_ni, guess_nj, guess_nk
		integer:: leftover_ni, leftover_nj, leftover_nk
		!----
		!Partition the domain into chunks for each processor
		!based on the number of processors specified in X,Y,Z direction.  Unequal sizes are permitted.
		!Actual grid coordinates are set in your_grid_function().
		!----

		if (myrank ==0)&
			write(*,'(a)') 'in your_partition_function...'

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
	end subroutine your_partition_function

	subroutine your_grid_function()
		implicit none
		!----
		integer:: i,j,k,ii,jj,kk
		real(LCSRP):: dx,dy,dz
		!----
		!Here set the grid coordinates x,y,z.
		!Just use a uniform Cartesian grid for now, arbitrary spacing is possible.
		!----
		if (myrank ==0)&
			write(*,'(a)') 'in your_grid_function...'

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
					z(ii,jj,kk) = 0.5*dz + real(k-1)*dz
				enddo
			enddo
		enddo

	end subroutine your_grid_function

	subroutine your_flow_solver(time)
		implicit none
		real(LCSRP),intent(in):: time
		!----
		integer:: i,j,k
		!----
		if (myrank ==0)&
			write(*,'(a)') 'in your_flow_solver...'
		
		!-----
		! Assumes periodicity of multiple 2pi in X,Y,Z
		!-----
		if (.NOT. allocated(u)) allocate(u(ni,nj,nk))
		if (.NOT. allocated(v)) allocate(v(ni,nj,nk))
		if (.NOT. allocated(w)) allocate(w(ni,nj,nk))
		do k =1,nk
		do j =1,nj
		do i =1,ni
			u(i,j,k) = (ABC_A+ABC_D*sin(time))*sin(z(i,j,k)) + ABC_C*cos(y(i,j,k))
			v(i,j,k) = ABC_B*sin(x(i,j,k)) + (ABC_A+ABC_D*sin(time))*cos(z(i,j,k))
			w(i,j,k) = ABC_C*sin(y(i,j,k)) + ABC_B*cos(x(i,j,k))
		enddo
		enddo
		enddo
	end subroutine your_flow_solver

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

end program abc_flow 
