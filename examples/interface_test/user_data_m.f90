module user_data_m
	implicit none
	include 'cfd2lcs_f.h'

	!-----
	!Domain dimensions
	!-----
	real(LCSRP), parameter:: PI = 4.0*atan(1.0)
	real(LCSRP), parameter:: LX = 2.0*PI
	real(LCSRP), parameter:: LY = 2.0*PI
	real(LCSRP), parameter:: LZ = 2.0*PI

	!-----
	!Total number of grid points in each direction
	!-----
	integer, parameter:: NX = 4
	integer, parameter:: NY = 32
	integer, parameter:: NZ = 64

	!-----
	!Boundary conditions for the domain exterior:
	!-----
	integer(LCSIP),parameter:: BC_IMIN = LCS_PERIODIC
	integer(LCSIP),parameter:: BC_JMIN = LCS_PERIODIC
	integer(LCSIP),parameter:: BC_KMIN = LCS_PERIODIC
	integer(LCSIP),parameter:: BC_IMAX = LCS_PERIODIC
	integer(LCSIP),parameter:: BC_JMAX = LCS_PERIODIC
	integer(LCSIP),parameter:: BC_KMAX = LCS_PERIODIC

	!-----
	!"Simulation" parameters
	!-----
	real(LCSRP),parameter:: DT = 1e-2
	real(LCSRP),parameter:: START_TIME = 0.0
	real(LCSRP),parameter:: END_TIME = 0.001

end module user_data_m


