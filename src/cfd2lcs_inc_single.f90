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

!
! The default string length for lcs labels
!
integer,parameter:: LCS_NAMELEN = 32

!
! Boundary condition flags:
!
integer(LCSIP),parameter:: &
	LCS_INTERNAL = 0, &
	LCS_INFLOW = 1, &
	LCS_OUTFLOW = 2, &
	LCS_WALL = 3, &
	LCS_SLIP = 4, &
	LCS_2D = 5

!
! Define the different types of LCS diagnostics here
!
integer(LCSIP),parameter:: &
	FTLE_FWD = 0, &
	FTLE_BKWD = 1, &
	LP_TRACER = 2
	
!
!Integration methods:
!
integer,parameter:: &
	EULER = 0, &
	TRAPEZOIDAL = 1, &
	RK2 = 2, &
	RK3 = 3, &
	RK4 = 4

!
! Interpolation  methods
!
integer,parameter:: &
	NEAREST_NBR = 0, &
	LINEAR = 1, &
	QUADRATIC = 2, &
	CUBIC = 3, &
	IDW = 4, &
	TSE = 5
