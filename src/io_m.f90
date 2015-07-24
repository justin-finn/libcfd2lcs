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
		R2_DATA = 2

	character(len=25),parameter:: ACTION_STRING(3)= (/&
		' READING         : ', &
		' WRITING         : ', &
		' WRITING (APPEND): ' /)

	contains

	subroutine structured_io(fname,IO_ACTION,global_size,offset,r0,r1,r2)
     	IMPLICIT NONE
		!-----
		character(len=*):: fname
		integer:: IO_ACTION
		integer,dimension(3):: global_size
		integer,dimension(3):: offset
		type(sr0_t),optional:: r0
		type(sr1_t),optional:: r1
		type(sr2_t),optional:: r2
		!-----
		integer:: NVAR, WORK_DATA
		character(len=32),allocatable:: dataname(:)
		character(len=32):: groupname
		real(LCSRP),allocatable :: data (:,:,:)  ! Write buffer
		integer,parameter:: NDIM = 3  !all data considered 3 dimensional
		integer(HID_T) :: file_id       ! File identifier
		integer(HID_T) :: dset_id       ! Dataset identifier
		integer(HID_T) :: filespace     ! Dataspace identifier in file
		integer(HID_T) :: memspace      ! Dataspace identifier in memory
		integer(HID_T) :: plist_id      ! Property list identifier
		integer(HSIZE_T), DIMENSION(NDIM) :: local_size ! Processor array dimension
		integer(HSIZE_T),  DIMENSION(NDIM) :: data_count = 1
		integer(HSIZE_T),  DIMENSION(NDIM) :: data_stride = 1
		integer :: error  ! Error flags
		integer :: info = MPI_INFO_NULL
		INTEGER(HID_T):: group_id
		integer:: i,j,k,ivar
		integer:: ni,nj,nk
		!-----

		!
		!Figure out what type of data we will dump, r0,r1,r2 or grid.
		!Only allow one type to be dumped in a single call to this routine
		!to dump multiple datasets to the same file, call with IO_ACTION = APPEND
		!
		if (present(r0)) then
			NVAR = 1
			WORK_DATA = R0_DATA
			ni = r0%ni; nj = r0%nj; nk = r0%nk

			write(groupname,'(a)') '/'
			allocate(dataname(1))
			write(dataname(1),'(a,a)')  trim(groupname),trim(r0%label)
			if(lcsrank==0) &
				write(*,'(a,a,a)') 'In structured_io... ',ACTION_STRING(IO_ACTION),trim(r0%label)
		elseif(present(r1)) then
			NVAR = 3
			WORK_DATA = R1_DATA
			ni = r1%ni; nj = r1%nj; nk = r1%nk
			write(groupname,'(a)')  trim(r1%label)
			allocate(dataname(3))
			write(dataname(1),'(a,a,a,a)')   trim(groupname),'/',trim(r1%label),'-X'
			write(dataname(2),'(a,a,a,a)')   trim(groupname),'/',trim(r1%label),'-Y'
			write(dataname(3),'(a,a,a,a)')   trim(groupname),'/',trim(r1%label),'-Z'
			if(lcsrank==0) &
				write(*,'(a,a,a)') 'In structured_io... ',ACTION_STRING(IO_ACTION),trim(r1%label)
		elseif(present(r2)) then
			NVAR = 9
			WORK_DATA = R2_DATA
			ni = r2%ni; nj = r2%nj; nk = r2%nk
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
				write(*,'(a,a,a)') 'In structured_io... ',ACTION_STRING(IO_ACTION),trim(r2%label)
		else
			if(lcsrank==0) &
				write(*,'(a,a)') 'In structured_io...  No data present'
			return
		endif
		local_size(1) = ni
		local_size(2) = nj
		local_size(3) = nk

		!
		! Allocate read/write buffer
		!
		allocate(data(1:ni,1:nj,1:nk))

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
					if(WORK_DATA == R0_DATA) then
						data(1:ni,1:nj,1:nk) = r0%r(1:ni,1:nj,1:nk)
					elseif(WORK_DATA == R1_DATA) then
						data(1:ni,1:nj,1:nk) = r1%x(1:ni,1:nj,1:nk)
					elseif(WORK_DATA == R2_DATA) then
						data(1:ni,1:nj,1:nk) = r2%xx(1:ni,1:nj,1:nk)
					endif
				case(2)
					if(WORK_DATA == R0_DATA) then
						cycle
					elseif(WORK_DATA == R1_DATA) then
						data(1:ni,1:nj,1:nk) = r1%y(1:ni,1:nj,1:nk)
					elseif(WORK_DATA == R2_DATA) then
						data(1:ni,1:nj,1:nk) = r2%xy(1:ni,1:nj,1:nk)
					endif
				case(3)
					if(WORK_DATA == R0_DATA) then
						cycle
					elseif(WORK_DATA == R1_DATA) then
						data(1:ni,1:nj,1:nk) = r1%z(1:ni,1:nj,1:nk)
					elseif(WORK_DATA == R2_DATA) then
						data(1:ni,1:nj,1:nk) = r2%xz(1:ni,1:nj,1:nk)
					endif
				case(4)
					data(1:ni,1:nj,1:nk) = r2%yx(1:ni,1:nj,1:nk)
				case(5)
					data(1:ni,1:nj,1:nk) = r2%yy(1:ni,1:nj,1:nk)
				case(6)
					data(1:ni,1:nj,1:nk) = r2%yz(1:ni,1:nj,1:nk)
				case(7)
					data(1:ni,1:nj,1:nk) = r2%zx(1:ni,1:nj,1:nk)
				case(8)
					data(1:ni,1:nj,1:nk) = r2%zy(1:ni,1:nj,1:nk)
				case(9)
					data(1:ni,1:nj,1:nk) = r2%zz(1:ni,1:nj,1:nk)
				case default
					cycle
				end select
			endif


			!
			! Create the data space for the  dataset.
			!
			CALL h5screate_simple_f(NDIM, int(global_size,HSIZE_T), filespace, error)
			CALL h5screate_simple_f(NDIM, local_size, memspace, error)

			!
			! Create chunked dataset.
			!
			CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
			CALL h5pset_chunk_f(plist_id, NDIM, local_size, error)
			select case(LCSRP) !Handle single or double precision:
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
				write(*,*) 'Error: bad LCSRP'
				CFD2LCS_ERROR = 5
				return
			end select



			!
			! Select hyperslab in the file.
			!
			CALL h5dget_space_f(dset_id, filespace, error)
			CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, int(offset,HSSIZE_T), data_count, error, &
			                                data_stride, local_size)

			!
			! Create property list for collective dataset write
			!
			CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
			CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

			!
			! Read/Write the dataset collectively.
			!
			select case(LCSRP) !Handle single or double precision:
			case(4)
				if(IO_ACTION == IO_WRITE .OR. IO_ACTION == IO_APPEND) then
					CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, data, int(global_size,HSIZE_T), error, &
			               file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
				endif
				if(IO_ACTION == IO_READ) then
					CALL h5dread_f(dset_id, H5T_NATIVE_REAL, data, int(global_size,HSIZE_T), error, &
			               file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
				endif
			case(8)
				if(IO_ACTION == IO_WRITE .OR. IO_ACTION == IO_APPEND) then
					CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, int(global_size,HSIZE_T), error, &
				               file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
				endif
				if(IO_ACTION == IO_READ) then
					CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, int(global_size,HSIZE_T), error, &
			               file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
				endif
			case default
				write(*,*) 'Error: bad LCSRP'
				CFD2LCS_ERROR = 5
				return
			end select


			!
			! Copy from read buffer
			!
			if(IO_ACTION == IO_READ) then
				select case(ivar)
				case(1)
					if(WORK_DATA == R0_DATA) then
						r0%r(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
					elseif(WORK_DATA == R1_DATA) then
						r1%x(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
					elseif(WORK_DATA == R2_DATA) then
						r2%xx(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
					endif
				case(2)
					if(WORK_DATA == R0_DATA) then
						cycle
					elseif(WORK_DATA == R1_DATA) then
						r1%y(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
					elseif(WORK_DATA == R2_DATA) then
						r2%xy(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
					endif
				case(3)
					if(WORK_DATA == R0_DATA) then
						cycle
					elseif(WORK_DATA == R1_DATA) then
						r1%z(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
					elseif(WORK_DATA == R2_DATA) then
						r2%xz(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
					endif
				case(4)
					r2%yx(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
				case(5)
					r2%yy(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
				case(6)
					r2%yz(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
				case(7)
					r2%zx(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
				case(8)
					r2%zy(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
				case(9)
					r2%zz(1:ni,1:nj,1:nk) = data(1:ni,1:nj,1:nk)
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

		deallocate(data)

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

	end subroutine structured_io
end module io_m

