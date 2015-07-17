module precision_m
	implicit none

	integer, parameter, private :: SP = kind(1.0)
	integer, parameter, private :: DP = kind(1.0d0)

	! WP is the working precision...
	integer, parameter :: WP = SP

	real(WP), private :: sample_real
	real(WP), parameter :: MAX_REAL_WP = HUGE(sample_real)

	integer, private :: sample_int
	integer, parameter :: MAX_INTEGER = HUGE(sample_int)

end module precision_m
