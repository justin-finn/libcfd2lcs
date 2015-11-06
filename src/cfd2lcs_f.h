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
! Set the verbosity of output
!
logical, parameter:: LCS_VERBOSE = .FALSE.

!
! Boundary condition flags:
!
integer(LCSIP),parameter:: &
	LCS_PERIODIC = 0, &
	LCS_INFLOW = 1, &
	LCS_OUTFLOW = 2, &
	LCS_WALL = 3, &
	LCS_SLIP = 4

!
! Set the Preallocate flag for communication buffers
!
logical,parameter:: PREALLOCATE_BUFFERS = .TRUE.


!
! Define the different types of LCS diagnostics here
!
integer(LCSIP),parameter:: &
	FTLE_FWD = 0, &
	FTLE_BKWD = 1, &
	LP_TRACER = 2

!
! Define the grid options for computing ftle fields
!
integer(LCSIP),parameter:: &
	LCS_CFD_GRID = 0, &
	LCS_NATIVE_GRID = 1
