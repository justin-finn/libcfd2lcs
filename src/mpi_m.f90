!
! mpi module
!
module mpi_m
	implicit none
	!----
	INCLUDE 'mpif.h'
	!----
	integer :: nprocs
	integer :: lcsrank
	integer :: lcscomm
	!----

	contains

	subroutine init_lcs_mpi(cfdcomm)
		implicit none
		!----
		integer, intent(in):: cfdcomm
		!----
		integer:: ierr
		!----
		!Duplicate the MPI comm, so that we have our
		!own communications for LCS related stuff.
		!----

		call MPI_COMM_DUP(cfdcomm,lcscomm,ierr)
		call MPI_COMM_RANK(lcscomm,lcsrank,ierr)
		call MPI_COMM_SIZE(lcscomm,nprocs,ierr)

		if (lcsrank==0) &
			write(*,'(a,i6,a)') 'in init_lcs_mpi... Using ',nprocs, ' MPI processes.'
	end subroutine init_lcs_mpi

end module mpi_m
