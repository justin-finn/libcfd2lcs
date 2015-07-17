module io_m
	use mpi_m
	use data_m
	use hdf5
	implicit none
	!----
	!Read and write routine for r0,r1,r2 datatypes
	!All I/O operations completed using hdf5 library
	!----

	integer,parameter:: &
		IO_READ = 1, &
		IO_WRITE = 2, &
		IO_APPEND = 3

	integer,parameter:: &
		R0_DATA = 0, &
		R1_DATA = 1, &
		R2_DATA = 2, &
		GRID_DATA = 3

	character(len=25),parameter:: ACTION_STRING(3)= (/&
		' READING         : ', &
		' WRITING         : ', &
		' WRITING (APPEND): ' /)

	contains

	subroutine write_structured_data(fname,IO_ACTION,cart,r0,r1,r2)
     	IMPLICIT NONE
		!-----
		character(len=*):: fname
		integer:: IO_ACTION
		type(cart_t):: cart
		type(sr0_t),optional:: r0
		type(sr1_t),optional:: r1
		type(sr2_t),optional:: r2
		!-----
		integer:: NVAR, WORK_DATA
		character(len=32),allocatable:: dataname(:)
		character(len=32):: groupname
		real(WP) :: data (1:cart%ni,1:cart%nj,1:cart%nk)  ! Write buffer
		integer,parameter:: NDIM = 3  !all data considered 3 dimensional
		integer(HID_T) :: file_id       ! File identifier
		integer(HID_T) :: dset_id       ! Dataset identifier
		integer(HID_T) :: filespace     ! Dataspace identifier in file
		integer(HID_T) :: memspace      ! Dataspace identifier in memory
		integer(HID_T) :: plist_id      ! Property list identifier
		integer(HSIZE_T), DIMENSION(NDIM) :: global_size ! Dataset dimensions in the file
		integer(HSIZE_T), DIMENSION(NDIM) :: local_size ! Processor array dimension
		integer(HSSIZE_T), DIMENSION(NDIM) :: offset
		integer(HSIZE_T),  DIMENSION(NDIM) :: data_count = 1
		integer(HSIZE_T),  DIMENSION(NDIM) :: data_stride = 1
		integer :: error  ! Error flags
		integer :: info = MPI_INFO_NULL
		INTEGER(HID_T):: group_id
		integer:: i,j,k,ivar
		!-----

		!
		!Figure out what type of data we will dump, r0,r1,r2 or grid.
		!Only allow one type to be dumped in a single call to this routine
		!to dump multiple datasets to the same file, call with IO_ACTION = APPEND
		!
		if (present(r0)) then
			NVAR = 1
			WORK_DATA = R0_DATA
			write(groupname,'(a)') '/'
			allocate(dataname(1))
			write(dataname(1),'(a,a)')  trim(groupname),trim(r0%label)
			if(lcsrank==0) &
				write(*,'(a,a,a)') 'In write_structured_data... ',ACTION_STRING(IO_ACTION),trim(r0%label)
		elseif(present(r1)) then
			NVAR = 3
			WORK_DATA = R1_DATA
			write(groupname,'(a)')  trim(r1%label)
			allocate(dataname(3))
			write(dataname(1),'(a,a,a,a)')   trim(groupname),'/',trim(r1%label),'-X'
			write(dataname(2),'(a,a,a,a)')   trim(groupname),'/',trim(r1%label),'-Y'
			write(dataname(3),'(a,a,a,a)')   trim(groupname),'/',trim(r1%label),'-Z'
			if(lcsrank==0) &
				write(*,'(a,a,a)') 'In write_structured_data... ',ACTION_STRING(IO_ACTION),trim(r1%label)
		elseif(present(r2)) then
			NVAR = 9
			WORK_DATA = R2_DATA
			write(groupname,'(a)')  trim(r2%label)
			allocate(dataname(9))
			write(dataname(1),'(a,a,a,a)')   trim(groupname),'/',trim(r2%label),'-XX'
			write(dataname(2),'(a,a,a,a)')   trim(groupname),'/',trim(r2%label),'-XY'
			write(dataname(3),'(a,a,a,a)')   trim(groupname),'/',trim(r2%label),'-XZ'
			write(dataname(4),'(a,a,a,a)')   trim(groupname),'/',trim(r2%label),'-YX'
			write(dataname(5),'(a,a,a,a)')   trim(groupname),'/',trim(r2%label),'-YY'
			write(dataname(6),'(a,a,a,a)')   trim(groupname),'/',trim(r2%label),'-YZ'
			write(dataname(7),'(a,a,a,a)')   trim(groupname),'/',trim(r2%label),'-ZX'
			write(dataname(8),'(a,a,a,a)')   trim(groupname),'/',trim(r2%label),'-ZY'
			write(dataname(9),'(a,a,a,a)')   trim(groupname),'/',trim(r2%label),'-ZZ'
			if(lcsrank==0) &
				write(*,'(a,a,a)') 'In write_structured_data... ',ACTION_STRING(IO_ACTION),trim(r2%label)
		else
			NVAR = 3
			WORK_DATA = GRID_DATA
			write(groupname,'(a)')  'grid'
			allocate(dataname(3))
			write(dataname(1),'(a,a)')   trim(groupname),'/X'
			write(dataname(2),'(a,a)')   trim(groupname),'/Y'
			write(dataname(3),'(a,a)')   trim(groupname),'/Z'
			if(lcsrank==0) &
				write(*,'(a,a,a)') 'In write_structured_data... ',ACTION_STRING(IO_ACTION),'Grid'
		endif

		!
		! Set some sizes:
		!
		global_size(1) = cart%gni
		global_size(2) = cart%gnj
		global_size(3) = cart%gnk
		local_size(1) = cart%ni
		local_size(2) = cart%nj
		local_size(3) = cart%nk
		offset(1) = cart%offset_i
		offset(2) = cart%offset_j
		offset(3) = cart%offset_k

		!
		! Initialize HDF5 library and Fortran interfaces.
		!
		CALL h5open_f(error)

		!
		! Setup file access property list with parallel I/O access.
		!
		CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
		CALL h5pset_fapl_mpio_f(plist_id, lcscomm, info, error)

		!
		! Create/open the file collectively.
		!
		select case(IO_ACTION)
		case(IO_WRITE)
			CALL h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
		case(IO_READ)
			CALL h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
		case(IO_APPEND)
			CALL h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
		end select
		if(error/=0) then
			write(*,*) 'Error opening file', fname
			CFD2LCS_ERROR = 4
			return
		endif
		CALL h5pclose_f(plist_id, error)

		!
		! Create write group for the data (except for R0)
		!
		if(WORK_DATA /= R0_DATA) then
		if(IO_ACTION == IO_WRITE .OR. IO_ACTION == IO_APPEND) then
			CALL h5gcreate_f(file_id,trim(groupname),group_id,error)
			CALL h5gclose_f(group_id,error)
		endif
		endif

		!Loop through each variable and read/write the data
		do ivar = 1,NVAR

			!
			! Set the data:
			!
			if(IO_ACTION == IO_WRITE .OR. IO_ACTION == IO_APPEND) then
				select case(ivar)
				case(1)
					if (WORK_DATA == GRID_DATA) then
						do k =1,cart%nk
						do j =1,cart%nj
						do i =1,cart%ni
							data(i,j,k) = cart%x(i)
						enddo
						enddo
						enddo
					elseif(WORK_DATA == R0_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r0%r(1:cart%ni,1:cart%nj,1:cart%nk)
					elseif(WORK_DATA == R1_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r1%x(1:cart%ni,1:cart%nj,1:cart%nk)
					elseif(WORK_DATA == R2_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r2%xx(1:cart%ni,1:cart%nj,1:cart%nk)
					endif
				case(2)
					if (WORK_DATA == GRID_DATA) then
						do k =1,cart%nk
						do j =1,cart%nj
						do i =1,cart%ni
							data(i,j,k) = cart%y(j)
						enddo
						enddo
						enddo
					elseif(WORK_DATA == R0_DATA) then
						cycle
					elseif(WORK_DATA == R1_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r1%y(1:cart%ni,1:cart%nj,1:cart%nk)
					elseif(WORK_DATA == R2_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r2%xy(1:cart%ni,1:cart%nj,1:cart%nk)
					endif
				case(3)
					if (WORK_DATA == GRID_DATA) then
						do k =1,cart%nk
						do j =1,cart%nj
						do i =1,cart%ni
							data(i,j,k) = cart%z(k)
						enddo
						enddo
						enddo
					elseif(WORK_DATA == R0_DATA) then
						cycle
					elseif(WORK_DATA == R1_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r1%z(1:cart%ni,1:cart%nj,1:cart%nk)
					elseif(WORK_DATA == R2_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r2%xz(1:cart%ni,1:cart%nj,1:cart%nk)
					endif
				case(4)
					if(WORK_DATA == R2_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r2%yx(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case(5)
					if(WORK_DATA == R2_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r2%yy(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case(6)
					if(WORK_DATA == R2_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r2%yz(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case(7)
					if(WORK_DATA == R2_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r2%zx(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case(8)
					if(WORK_DATA == R2_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r2%zy(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case(9)
					if(WORK_DATA == R2_DATA) then
						data(1:cart%ni,1:cart%nj,1:cart%nk) = r2%zz(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case default
				cycle
				end select
			endif


			!
			! Create the data space for the  dataset.
			!
			CALL h5screate_simple_f(NDIM, global_size, filespace, error)
			CALL h5screate_simple_f(NDIM, local_size, memspace, error)

			!
			! Create chunked dataset.
			!
			CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
			CALL h5pset_chunk_f(plist_id, NDIM, local_size, error)
			select case(WP) !Handle single or double precision:
			case(4)
				if(IO_ACTION == IO_WRITE .OR. IO_ACTION == IO_APPEND) then
					CALL h5dcreate_f(file_id,trim(dataname(ivar)),H5T_NATIVE_REAL,filespace,dset_id,error,plist_id)
				endif
				if(IO_ACTION == IO_READ) then
					CALL h5dopen_f(file_id,trim(dataname(ivar)),dset_id,error)
				endif
			case(8)
				if(IO_ACTION == IO_WRITE .OR. IO_ACTION == IO_APPEND) then
					CALL h5dcreate_f(file_id,trim(dataname(ivar)),H5T_NATIVE_DOUBLE,filespace,dset_id,error,plist_id)
				endif
				if(IO_ACTION == IO_READ) then
					CALL h5dopen_f(file_id,trim(dataname(ivar)),dset_id,error)
				endif
			case default
				write(*,*) 'Error: bad WP'
				CFD2LCS_ERROR = 5
				return
			end select



			!
			! Select hyperslab in the file.
			!
			CALL h5dget_space_f(dset_id, filespace, error)
			CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, data_count, error, &
			                                data_stride, local_size)

			!
			! Create property list for collective dataset write
			!
			CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
			CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

			!
			! Read/Write the dataset collectively.
			!
			select case(WP) !Handle single or double precision:
			case(4)
				if(IO_ACTION == IO_WRITE .OR. IO_ACTION == IO_APPEND) then
					CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, data, global_size, error, &
			               file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
				endif
				if(IO_ACTION == IO_READ) then
					CALL h5dread_f(dset_id, H5T_NATIVE_REAL, data, global_size, error, &
			               file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
				endif
			case(8)
				if(IO_ACTION == IO_WRITE .OR. IO_ACTION == IO_APPEND) then
					CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, global_size, error, &
				               file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
				endif
				if(IO_ACTION == IO_READ) then
					CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, global_size, error, &
			               file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
				endif
			case default
				write(*,*) 'Error: bad WP'
				CFD2LCS_ERROR = 5
				return
			end select


			!
			! Copy from read buffer
			!
			if(IO_ACTION == IO_READ) then
				select case(ivar)
				case(1)
					if (WORK_DATA == GRID_DATA) then
						do k =1,cart%nk
						do j =1,cart%nj
						do i =1,cart%ni
							cart%x(k) = data(i,j,k)
						enddo
						enddo
						enddo
					elseif(WORK_DATA == R0_DATA) then
						r0%r(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					elseif(WORK_DATA == R1_DATA) then
						r1%x(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					elseif(WORK_DATA == R2_DATA) then
						r2%xx(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					endif
				case(2)
					if (WORK_DATA == GRID_DATA) then
						do k =1,cart%nk
						do j =1,cart%nj
						do i =1,cart%ni
							cart%y(k) = data(i,j,k)
						enddo
						enddo
						enddo
					elseif(WORK_DATA == R0_DATA) then
						cycle
					elseif(WORK_DATA == R1_DATA) then
						r1%y(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					elseif(WORK_DATA == R2_DATA) then
						r2%xy(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					endif
				case(3)
					if (WORK_DATA == GRID_DATA) then
						do k =1,cart%nk
						do j =1,cart%nj
						do i =1,cart%ni
							cart%z(k) = data(i,j,k)
						enddo
						enddo
						enddo
					elseif(WORK_DATA == R0_DATA) then
						cycle
					elseif(WORK_DATA == R1_DATA) then
						r1%z(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					elseif(WORK_DATA == R2_DATA) then
						r2%xz(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					endif
				case(4)
					if(WORK_DATA == R2_DATA) then
						r2%yx(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case(5)
					if(WORK_DATA == R2_DATA) then
						r2%yy(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case(6)
					if(WORK_DATA == R2_DATA) then
						r2%yz(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case(7)
					if(WORK_DATA == R2_DATA) then
						r2%zx(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case(8)
					if(WORK_DATA == R2_DATA) then
						r2%zy(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case(9)
					if(WORK_DATA == R2_DATA) then
						r2%zz(1:cart%ni,1:cart%nj,1:cart%nk) = data(1:cart%ni,1:cart%nj,1:cart%nk)
					else
						cycle
					endif
				case default
				cycle
				end select
			endif

			!
			! Close the dataspace/dataset.
			!
			CALL h5sclose_f(filespace, error) !close dataspace
			CALL h5sclose_f(memspace, error) !close dataspace
			CALL h5dclose_f(dset_id, error)  !close dataset

		enddo

		! Close the file
		CALL h5pclose_f(plist_id, error) !close property list
		CALL h5fclose_f(file_id, error) !close file
		CALL h5close_f(error) !close fortran interfaces and H5 library

	contains
	subroutine checkio(point)
		integer:: ierr
		integer:: point
		if(error/=0) then
			write(*,*) 'myrank[',lcsrank,'] error at point',point
			call mpi_barrier(lcscomm,ierr)
			stop
		else
			write(*,*) 'myrank[',lcsrank,'] ok at point',point
			call mpi_barrier(lcscomm,ierr)
		endif
	end subroutine checkio

	end subroutine write_structured_data
end module io_m

