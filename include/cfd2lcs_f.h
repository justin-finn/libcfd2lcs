!********************************************************
!CFD2LCS Fortran Include file.
!Contains interface definitions needed for user-level API
!********************************************************

!
! Set the precision: (single =4, double = 8)
!	LCSRP = REAL precision
!	LCSIP = INTEGER precision
!
integer, parameter 	:: LCSRP = 4
integer, parameter 	:: LCSIP = 4
integer, parameter:: MPI_LCSRP = LCSRP
integer, parameter:: MPI_LCSIP = LCSIP

!
! Boundary condition flags:
!
integer(LCSIP),parameter:: &
	LCS_PERIODIC = 0, &
	LCS_INFLOW = 1, &
	LCS_OUTFLOW = 2, &
	LCS_WALL = 3
