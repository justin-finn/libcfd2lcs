!********************************************************
!CFD2LCS Fortran Include file.
!Contains interface definitions needed for user-level API
!********************************************************

!
! Set the precision: (single =4, double = 8)
!
integer, parameter :: LCS_PRECISION = 4

!
!	LCSRP = REAL precision
!	LCSIP = INTEGER precision
!
integer, parameter 	:: LCSRP = LCS_PRECISION
integer, parameter 	:: LCSIP = 4

!
! Boundary condition flags:
!
integer(LCSIP),parameter:: &
	LCS_PERIODIC = 0, &
	LCS_INFLOW = 1, &
	LCS_OUTFLOW = 2, &
	LCS_WALL = 3

!
! Set the Preallocate flag for communication buffers
!
logical,parameter:: PREALLOCATE_BUFFERS = .TRUE.
