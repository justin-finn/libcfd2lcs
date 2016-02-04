module io_m
	use data_m
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

	character(len=4),parameter:: FILE_EXT = '.dat'

	contains

	subroutine write_lcs(lcs,time)
		implicit none
		!-----
		type(lcs_t):: lcs
		real(LCSRP):: time
		!-----
		integer:: gn(3),offset(3)
		character(len=128):: FMT1,fname
		integer(8):: findex
		!-----
		!Output a datafile containing the LCS diagnostic results
		!-----

		!-----
		!Generate the filename.
		!-----
		findex = nint(time/lcs%h,8)
		select case(findex)
			case(-999999999:-100000000)
				FMT1 = "(a,a,a,i10.9,a)"
			case(-99999999:-10000000)
				FMT1 = "(a,a,a,i9.8,a)"
			case(-9999999:-1000000)
				FMT1 = "(a,a,a,i8.7,a)"
			case(-999999:-100000)
				FMT1 = "(a,a,a,i7.6,a)"
			case(-99999:-10000)
				FMT1 = "(a,a,a,i6.5,a)"
			case(-9999:-1000)
				FMT1 = "(a,a,a,i5.4,a)"
			case(-999:-100)
				FMT1 = "(a,a,a,i4.3,a)"
			case(-99:-10)
				FMT1 = "(a,a,a,i3.2,a)"
			case(-9:-1)
				FMT1 = "(a,a,a,i2.1,a)"
			case(0:9)
				FMT1 = "(a,a,a,i1.1,a)"
			case(10:99)
				FMT1 = "(a,a,a,i2.2,a)"
			case(100:999)
				FMT1 = "(a,a,a,i3.3,a)"
			case(1000:9999)
				FMT1 = "(a,a,a,i4.4,a)"
			case(10000:99999)
				FMT1 = "(a,a,a,i5.5,a)"
			case(100000:999999)
				FMT1 = "(a,a,a,i6.6,a)"
			case(1000000:9999999)
				FMT1 = "(a,a,a,i7.7,a)"
			case(10000000:99999999)
				FMT1 = "(a,a,a,i8.8,a)"
			case(100000000:999999999)
				FMT1 = "(a,a,a,i9.9,a)"
			case default
				if(lcsrank==0)&
					write(*,*) 'ERROR: unsupported range for file index,',findex
				CFD2LCS_ERROR = 1
				return
		end select
		write(fname,'(a,i4.4,a,a,a,i10.9,a)')'./cfd2lcs_output/',lcsrank,'_',trim(lcs%label),'_',nint(time/lcs%h),FILE_EXT
		if(lcsrank==0)&
			write(*,*) 'In write_lcs...',trim(fname)

		!-----
		!Ouptut the LCS.  Data depends on diagnostic type.
		!-----
		select case (lcs%diagnostic)
			case(FTLE_FWD,FTLE_BKWD)
				gn = (/lcs%sgrid%gni,lcs%sgrid%gnj,lcs%sgrid%gnk/)
				offset = (/lcs%sgrid%offset_i,lcs%sgrid%offset_j,lcs%sgrid%offset_k/)
				call structured_io(trim(fname),IO_WRITE,gn,offset,r0=lcs%ftle)	!Append  the FTLE
			case(LP_TRACER)
				call unstructured_io(fname,IO_WRITE,r1=lcs%lp%xp)
				call unstructured_io(fname,IO_APPEND,r1=lcs%lp%up)
			case default
		end select

	end subroutine write_lcs


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
		character(len=LCS_NAMELEN),allocatable:: dataname(:)
		character(len=LCS_NAMELEN):: groupname
		real(LCSRP),allocatable :: data (:,:,:)  ! Write buffer
		integer,parameter:: NDIM = 3  !all data considered 3 dimensional
		integer:: i,j,k,ivar
		integer:: ni,nj,nk
		character(len=128):: junk
		integer:: funit
		integer,parameter:: funit_start = 1234
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


		!read/write the header
		funit = funit_start + lcsrank
		select case(IO_ACTION)
			case(IO_WRITE)
				open(funit,file=fname,status='replace',form='formatted')
				write(funit,'(a,a,a)') 'TITLE = "',trim(groupname),'"'
				write(funit,'(a)') 'VARIABLES ='
				do ivar = 1,NVAR
					write(funit,'(a,a,a)') '"',trim(dataname(ivar))	,'"'
				enddo
				write(funit,'(a)') 'ZONE'
				write(funit,'(a,a,a,i4.4,a)') 'T = "',trim(groupname),'_',lcsrank,'"'
				write(funit,'(a,i4.4,a,i4.4,a,i4.4)') 'I = ',ni,' J = ',nj, ' K = ',nk
				write(funit,'(a)') 'ZONETYPE = Ordered, DATAPACKING = BLOCK'
			case(IO_READ)
				open(funit,file=fname,action='read',form='formatted')
				read(funit,*) junk
				read(funit,*) junk
				do ivar = 1,NVAR
					read(funit,*) junk
				enddo
				read(funit,*) junk
				read(funit,*) junk
				read(funit,*) junk
				read(funit,*) junk
			case default
				if(lcsrank==0) &
				write(*,*) 'bad IO_ACTION.  Must be IO_WRITE or IO_READ'
				CFD2LCS_ERROR = 1
		end select
			

		!
		! Allocate read/write buffer
		!
		allocate(data(1:ni,1:nj,1:nk))

		!Loop through each variable and read/write the data
		do ivar = 1,NVAR

			!
			! Set the data:
			!
			!
			if(IO_ACTION == IO_WRITE) then
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

				!write data
				do k=1,nk
				do j=1,nj
				do i=1,ni
					write(funit,*) data(i,j,k)
				enddo
				enddo
				enddo
			endif

			!
			! Copy from read buffer
			!
			if(IO_ACTION == IO_READ) then
				!read data
				do k=1,nk
				do j=1,nj
				do i=1,ni
					read(funit,*) data(i,j,k)
				enddo
				enddo
				enddo
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

		enddo

		deallocate(data)
		close(funit)

	end subroutine structured_io

	subroutine unstructured_io(fname,IO_ACTION,r0,r1,r2)
     	IMPLICIT NONE
		!-----
		character(len=*):: fname
		integer:: IO_ACTION
		type(ur0_t),optional:: r0
		type(ur1_t),optional:: r1
		type(ur2_t),optional:: r2
		!-----
		
		if(lcsrank==0) &
				write(*,'(a,a,a)') 'In unstructured_io... NO ACTION (not supported):  Need to compile with HDF5 support'
	
	end subroutine unstructured_io

end module io_m

