!
! mpi module
!
module mpi_m
	use data_m
	implicit none
	!----
	INCLUDE 'mpif.h'
	!----
	integer :: nprocs
	integer :: lcsrank
	integer :: lcscomm
	!----
	integer:: MPI_LCSRP  !Real precision for mpi

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

end module mpi_m
