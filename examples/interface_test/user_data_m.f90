module user_data_m
	implicit none
	include 'cfd2lcs_f.h'

	!-----
	!Domain dimensions
	!-----
	real(LCSRP), parameter:: PI = 4.0*atan(1.0)
	real(LCSRP), parameter:: LX = 2.0
	real(LCSRP), parameter:: LY = 1.0
	real(LCSRP), parameter:: LZ = 0.5

	!-----
	!Total number of grid points in each direction
	!-----
	integer, parameter:: NX = 512 
	integer, parameter:: NY = 256
	integer, parameter:: NZ = 1 

	!-----
	!Boundary conditions for the domain exterior:
	!-----
	integer(LCSIP),parameter:: BC_IMIN = LCS_WALL
	integer(LCSIP),parameter:: BC_JMIN = LCS_WALL
	integer(LCSIP),parameter:: BC_KMIN = LCS_WALL
	integer(LCSIP),parameter:: BC_IMAX = LCS_WALL
	integer(LCSIP),parameter:: BC_JMAX = LCS_WALL
	integer(LCSIP),parameter:: BC_KMAX = LCS_WALL

	!-----
	!"Simulation" parameters
	!-----
	real(LCSRP),parameter:: DT = 0.1
	real(LCSRP),parameter:: START_TIME = 0.0
	real(LCSRP),parameter:: END_TIME = 5 

end module user_data_m


