!
! mpi module
!
module mpi_m
	use precision_m
	implicit none
	!----
	INCLUDE 'mpif.h'
	!----
	integer :: nprocs
	integer :: lcsrank
	integer :: lcscomm
	integer, parameter:: MPI_REAL_WP = WP
	integer,save:: CFD2LCS_ERROR
	!----

	contains

	subroutine init_lcs_mpi(cfdcomm)
		implicit none
		!----
		integer, intent(in):: cfdcomm
		!----
		integer:: proc,ierr,cfdrank,nprocs_cfd,maxrank
		!----

		!Error handling:
		CFD2LCS_ERROR = 0

		!Duplicate the MPI comm, so that we have our
		!own communications for LCS related stuff.
		call MPI_COMM_DUP(cfdcomm,lcscomm,ierr)
		call MPI_COMM_RANK(lcscomm,lcsrank,ierr)
		call MPI_COMM_SIZE(lcscomm,nprocs,ierr)
		if (lcsrank==0) &
			write(*,'(a,i6,a)') 'in init_lcs_mpi... Using ',nprocs, ' MPI processes.'

		!Check for any errors
		call cfd2lcs_error_check()

	end subroutine init_lcs_mpi







	subroutine cfd2lcs_error_check()
		implicit none
		!-----
		integer:: ierr,MAX_CFD2LCS_ERROR
		!-----
		!In the event of a cfd2lcs error, we dont necessarily want
		!to bring down the cfd solver.  So, flag an error instead:
		!-----

		if (CFD2LCS_ERROR == -1) then
			if (lcsrank==0) write(*,'(a)') &
				'FATAL CFD2LCS CFD2LCS_ERROR DETECTED, WILL NOT PERFORM LCS COMPUTATIONS'
			return
		endif

		!Check for any new errors:
		call MPI_ALLREDUCE(CFD2LCS_ERROR,MAX_CFD2LCS_ERROR,1,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
		if (MAX_CFD2LCS_ERROR > 0) then
			if(lcsrank == 0) then
				write(*,'(a)') '*******************************'
				write(*,'(a)') 'CFD2LCS_ERROR DETECTED IN CFD2LCS!'
			endif
			call MPI_BARRIER(lcscomm,ierr)
			if(CFD2LCS_ERROR > 0) write(*,'(a,i6,a)') 'rank[',lcsrank,']: ',trim(CFD2LCS_ERROR_MSG(CFD2LCS_ERROR))
			call MPI_BARRIER(lcscomm,ierr)
			if(lcsrank == 0) then
				write(*,'(a)') '*******************************'
			endif
			CFD2LCS_ERROR = -1
			call MPI_BARRIER(lcscomm,ierr)
		endif
	end subroutine cfd2lcs_error_check
	character(len=128) function CFD2LCS_ERROR_MSG(ERROR)
		integer:: ERROR
		select case (ERROR)
		case(1)
			CFD2LCS_ERROR_MSG = 'Trouble setting up communication pattern'
		case(2)
			CFD2LCS_ERROR_MSG = 'Bad CFD velocity array size'
		case(3)
			CFD2LCS_ERROR_MSG = 'Bad I/O'
		case(4)
			CFD2LCS_ERROR_MSG = 'Bad I/O ACTION'
		case default
			CFD2LCS_ERROR_MSG = 'Unknown Error'
		end select
	end function CFD2LCS_ERROR_MSG

end module mpi_m
