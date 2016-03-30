!
!Copyright (C) 2015-2016, Justin R. Finn.  All rights reserved.
!libcfd2lcs is distributed is under the terms of the GNU General Public License
!
!
! A Simple working example to show how to call
! CFD2LCS.  Subroutines starting with "your_"
! are independent of the cfd2lcs functionality,
! and are used only to create a simple dataset
! for this example.
!
!
! FLOW FIELD:  ROMS data
!
! This example also contains a simple interface to read input data from netcdf files.
!

module netcdf_m
	use netcdf
	implicit none
	contains
	subroutine netcdf_read_chunk(fname,varname,n,offset,r1d,r2d,r3d,i1d,i2d,i3d)
		use netcdf
		implicit none
		include 'cfd2lcs_inc_sp.f90'
		!-----
		character(len=*):: fname,varname
		integer:: n(3),offset(3)
		real(LCSRP),optional:: r1d(1:n(1))
		real(LCSRP),optional:: r2d(1:n(1),1:n(2))
		real(LCSRP),optional:: r3d(1:n(1),1:n(2),1:n(3))
		integer,optional:: i1d(1:n(1))
		integer,optional:: i2d(1:n(1),1:n(2))
		integer,optional:: i3d(1:n(1),1:n(2),1:n(3))
		!-----
	  	integer :: ncid,varid,status
		real(LCSRP)::missing_r
		integer:: i,j,k
		!-----
		! Open the file.
		call check( nf90_open(fname, nf90_nowrite, ncid) )
		! Get the varids of the pressure and temperature netCDF variables.
		call check( nf90_inq_varid(ncid, trim(varname), varid) )
		! Read the data
		if(present(r1d)) then
			call check( nf90_get_var(ncid, varid, r1d, start = offset+1, count = n) )
			status = nf90_get_att(ncid, varid, "missing_value", missing_r)
			do i = 1,n(1)
				if(r1d(i) == missing_r) r1d(i) = 0.0_LCSRP
			enddo
		endif
		if(present(r2d)) then
			call check( nf90_get_var(ncid, varid, r2d, start = offset+1, count = n) )
			status = nf90_get_att(ncid, varid, "missing_value", missing_r)
			do j = 1,n(2)
			do i = 1,n(1)
				if(r2d(i,j) == missing_r) r2d(i,j) = 0.0_LCSRP
			enddo
			enddo
		endif
		if(present(r3d)) then
			call check( nf90_get_var(ncid, varid, r3d, start = offset+1, count = n) )
			status = nf90_get_att(ncid, varid, "missing_value", missing_r)
			do k = 1,n(3)
			do j = 1,n(2)
			do i = 1,n(1)
				if(r3d(i,j,k) == missing_r) r3d(i,j,k) = 0.0_LCSRP
			enddo
			enddo
			enddo
		endif
		if(present(i1d)) then
			call check( nf90_get_var(ncid, varid, i1d, start = offset+1, count = n) )
		endif
		if(present(i2d)) then
			call check( nf90_get_var(ncid, varid, i2d, start = offset+1, count = n) )
		endif
		if(present(i3d)) then
			call check( nf90_get_var(ncid, varid, i3d, start = offset+1, count = n) )
		endif
		! Close the file. This frees up any internal netCDF resources
		! associated with the file.
		call check( nf90_close(ncid) )
		contains
		subroutine check(status)
			integer, intent ( in) :: status
			if(status /= nf90_noerr) then
			print *, trim(nf90_strerror(status))
			stop "Stopped"
			end if
		end subroutine check
	end subroutine netcdf_read_chunk
end module netcdf_m

program roms_lcs
	use netcdf_m
	implicit none
	!-----
	include 'cfd2lcs_inc_sp.f90'  !This includes parameter definitions needed for cfd2lcs
	include 'mpif.h'
	!******BEGIN USER INPUT********************
	!-----
	!CONSTANTS
	!-----
	real(LCSRP),parameter:: PI = 4.0*atan(1.0)
	real(LCSRP),parameter:: R_EARTH = 6371000.0_LCSRP  !Radius of earth [m]
	real(LCSRP),parameter:: DEG2RAD = PI/180.0_LCSRP !Convert degrees to radians
	real(LCSRP),parameter:: D2S = real(60*60*24,LCSRP)  !Convert days to seconds

	!-----
	!Total number of grid points in each direction
	!-----
	integer, parameter:: NX = 330
	integer, parameter:: NY = 440
	integer, parameter:: NZ = 1
	!-----
	!"Simulation" parameters
	!-----
	integer,parameter:: START_DAY = 1 
	integer,parameter:: END_DAY = 58 
	real(LCSRP),parameter:: START_TIME = real(START_DAY)*D2S
	real(LCSRP),parameter:: END_TIME = real(END_DAY)*D2S
	real(LCSRP),parameter:: DT = 1.0_LCSRP*D2S
	real(LCSRP),parameter:: T = 7.0_LCSRP*D2S
	real(LCSRP),parameter:: H = 1.0_LCSRP*D2S
	integer,parameter:: RESOLUTION = 1 
	real(LCSRP),parameter:: CFL = 0.45
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
	real(LCSRP), allocatable:: deglat(:,:), deglon(:,:)
	integer(LCSIP),allocatable:: flag(:,:,:)
	integer:: n(3),offset(3)
	integer(LCSIP):: id_fwd,id_bkwd
	integer:: day = START_DAY !initialize to start
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
	write(*,'(a)') 'libcfd2lcs roms test'
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
	!Set the boundary conditions:
	!-----
	call your_bc_function()

	!-----
	!Initialize cfd2lcs for your data
	!-----
	n = (/ni,nj,nk/)  !number of grid points for THIS partition
	offset = (/offset_i,offset_j,offset_k/)  !Global offset of these grid points
	call cfd2lcs_init(mycomm,n,offset,x,y,z,flag)

	!-----
	!Initialize LCS diagnostics
	!-----
	call cfd2lcs_diagnostic_init(id_fwd,FTLE_FWD,RESOLUTION,T,H,'fwdFTLE')
	call cfd2lcs_diagnostic_init(id_bkwd,FTLE_BKWD,RESOLUTION,T,H,'bkwdFTLE')

	!-----
	!Set cfd2lcs options/parameters
	!-----
	call cfd2lcs_set_option('SYNCTIMER',LCS_FALSE)
	call cfd2lcs_set_option('DEBUG',LCS_FALSE)
	call cfd2lcs_set_option('WRITE_FLOWMAP',LCS_FALSE)
	call cfd2lcs_set_option('WRITE_BCFLAG',LCS_FALSE)
	call cfd2lcs_set_option('INCOMPRESSIBLE',LCS_FALSE)
	call cfd2lcs_set_option('AUX_GRID',LCS_FALSE)
	call cfd2lcs_set_option('INTEGRATOR',RK2)
	call cfd2lcs_set_option('INTERPOLATOR',TSE_LIMIT)
	call cfd2lcs_set_param('CFL', CFL)

	!-----
	!***Start of your flow solver timestepping loop***
	!-----
	timestep = 0
	time = START_TIME
	day = START_DAY
	do while (time <= END_TIME)

		if(myrank == 0) then
			write(*,'(a)') '------------------------------------------------------------------'
			write(*,'(a,i10.0,a,ES11.4,a,ES11.4)') 'STARTING TIMESTEP #',timestep,': time = ',time,', DT = ',DT
			write(*,'(a)') '------------------------------------------------------------------'
		endif

		!Produce the new velocity field with your flow solver
		call your_flow_solver()

		!Update the LCS diagnostics using the new flow field
		call cfd2lcs_update(n,u,v,w,time)

		timestep = timestep + 1
		time = time + DT
		day = day + 1
	enddo

	!-----
	!Cleanup
	!-----
	deallocate(x)
	deallocate(y)
	deallocate(z)
	deallocate(deglat)
	deallocate(deglon)
	deallocate(u)
	deallocate(v)
	deallocate(w)
	call cfd2lcs_finalize(ierr)
	call MPI_FINALIZE(ierr)

	contains

	subroutine your_partition_function()
		implicit none
		!----
		integer:: guess_ni, guess_nj, guess_nk
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
		integer:: i,j,k
		character(len=128):: ROMS_FILE = './inputData/romsGrid.nc'
		!----
		!----
		!Here set the grid coordinates x,y,z.
		!We load in lon/lat and the mask from the roms grid file:
		!----

		if (myrank ==0)&
			write(*,*) 'in your_grid_function...',trim(ROMS_FILE)

		!Read in lat, lon, and mask:
		allocate(x(1:ni,1:nj,1:nk))
		allocate(y(1:ni,1:nj,1:nk))
		allocate(z(1:ni,1:nj,1:nk))
		allocate(deglat(1:ni,1:nj))
		allocate(deglon(1:ni,1:nj))
		call netcdf_read_chunk(trim(ROMS_FILE),'lat',(/ni,nj,0/),(/offset_i,offset_j,0/),r2d=deglat)
		call netcdf_read_chunk(trim(ROMS_FILE),'lon',(/ni,nj,0/),(/offset_i,offset_j,0/),r2d=deglon)

		!Convert to Cartesian x,y,z and set the "flag" according to the roms mask
		do k =1,nk
		do j =1,nj
		do i =1,ni
			x(i,j,k) = R_EARTH*cos(DEG2RAD*deglat(i,j))*cos(DEG2RAD*deglon(i,j))
			y(i,j,k) = R_EARTH*cos(DEG2RAD*deglat(i,j))*sin(DEG2RAD*deglon(i,j))
			z(i,j,k) = R_EARTH*sin(DEG2RAD*deglat(i,j))
		enddo
		enddo
		enddo

	end subroutine your_grid_function

	subroutine your_bc_function()
		implicit none
		!----
		integer:: i,j
		integer,allocatable:: mask(:,:)
		character(len=128):: ROMS_FILE = './inputData/romsGrid.nc'
		!----
		if(myrank==0) &
			write(*,'(a)') 'in your_bc_function...'

		allocate(mask(1:ni,1:nj))
		allocate(flag(1:ni,1:nj,1:nk))

		!Read in the mask:
		call netcdf_read_chunk(trim(ROMS_FILE),'mask',(/ni,nj,0/),(/offset_i,offset_j,0/),i2d=mask)

		!Default is LCS_INTERNAL everywhere
		flag = LCS_INTERNAL

		!Set outflow everywhere on the (2D) domain exterior:
		if(offset_i==0) flag(1,:,:) = LCS_OUTFLOW
		if(offset_j==0) flag(:,1,:) = LCS_OUTFLOW
		if(offset_i+ni==NX) flag(ni,:,:) = LCS_OUTFLOW
		if(offset_j+nj==NY) flag(:,nj,:) = LCS_OUTFLOW

		!Set the flag based on the roms mask:
		do j = 1,nj
		do i = 1,ni
			if(mask(i,j) == 0 ) flag(i,j,:) = LCS_MASK
		enddo
		enddo

	end subroutine your_bc_function

	subroutine your_flow_solver()
		implicit none
		!----
		integer:: i,j,k
		!----
		real(LCSRP):: my_umax(3), my_umin(3),umax(3),umin(3)
		integer:: ierr
		real(LCSRP):: c(3)
		character(len=128):: ROMS_FILE
		real(LCSRP), allocatable:: utmp(:,:),vtmp(:,:)
		!----
		if (myrank ==0)&
			write(*,'(a)') 'in your_flow_solver...'

		if (.NOT. allocated(u)) allocate(u(ni,nj,nk))
		if (.NOT. allocated(v)) allocate(v(ni,nj,nk))
		if (.NOT. allocated(w)) allocate(w(ni,nj,nk))
		allocate(utmp(1:ni,1:nj))
		allocate(vtmp(1:ni,1:nj))

		!-----
		!Read in the u,v data (as 2d arrays)
		!set w=0 everywhere
		!-----
		if(day < 10) write(ROMS_FILE,'(a,i1.1,a)') './inputData/romsDay_',day,'.nc'
		if(day > 10) write(ROMS_FILE,'(a,i2.2,a)') './inputData/romsDay_',day,'.nc'
		if(myrank==0)write(*,*) 'Reading data from ',trim(ROMS_FILE)
		call netcdf_read_chunk(ROMS_FILE,'u',(/ni,nj,0/),(/offset_i,offset_j,0/),r2d=utmp)
		call netcdf_read_chunk(ROMS_FILE,'v',(/ni,nj,0/),(/offset_i,offset_j,0/),r2d=vtmp)
		u(1:ni,1:nj,1) = utmp(1:ni,1:nj)
		v(1:ni,1:nj,1) = vtmp(1:ni,1:nj)
		w = 0.0_LCSRP

		!-----
		!Convert from polar to cartesian
		!-----
		do k = 1,nk
		do j = 1,nj
		do i = 1,ni
			c(1:3) = sphere2cart(u(i,j,k),v(i,j,k),w(i,j,k),deglon(i,j),deglat(i,j))
			u(i,j,k) = c(1)
			v(i,j,k) = c(2)
			w(i,j,k) = c(3)
		enddo
		enddo
		enddo

		!-----
		!Check max/min
		!-----
		my_umax(1) = maxval(u); my_umin(1) = minval(u)
		my_umax(2) = maxval(v); my_umin(2) = minval(v)
		my_umax(3) = maxval(w); my_umin(3) = minval(w)
		call MPI_REDUCE(my_umax,umax,3,MPI_REAL,MPI_MAX,0,mycomm,ierr)
		call MPI_REDUCE(my_umin,umin,3,MPI_REAL,MPI_MIN,0,mycomm,ierr)
		if(myrank==0) then
			write(*,*)'myrank[',myrank,'] min/max [u]', umin(1),umax(1)
			write(*,*)'myrank[',myrank,'] min/max [v]', umin(2),umax(2)
			write(*,*)'myrank[',myrank,'] min/max [w]', umin(3),umax(3)
		endif

		deallocate(utmp)
		deallocate(vtmp)

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

	function sphere2cart(vaz,vel,vr,az,el)
		implicit none
		real(LCSRP):: sphere2cart(3)
		real(LCSRP):: vaz,vel,vr
		real(LCSRP):: az,el
		!-----
		real(LCSRP),parameter:: DEG2RAD = PI/180.0_LCSRP
		!-----
		!Convert spherical vector (vaz,vel,vr) to cartesian (x,y,z)
		!-----
		sphere2cart(1) = -vaz*sin(DEG2RAD*az) -vel*sin(DEG2RAD*el)*cos(DEG2RAD*az)&
			 + vr*cos(DEG2RAD*el)*cos(DEG2RAD*az)
		sphere2cart(2) = vaz*cos(DEG2RAD*az) -vel*sin(DEG2RAD*el)*sin(DEG2RAD*az)&
			 + vr*cos(DEG2RAD*el)*sin(DEG2RAD*az)
		sphere2cart(3) = vaz*0.0_LCSRP + vel*cos(DEG2RAD*el) + vr*sin(DEG2RAD*el)
	end function sphere2cart

end program roms_lcs
