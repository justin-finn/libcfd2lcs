/*
CFD2LCS C/C++ Include file.
Contains interface definitions needed for user-level API
To use:  #include "cfd2lcs_inc.h"
*/

#include <mpi.h>

/*
Set the precision here
*/
typedef float lcsdata_t;

/*
The string namelength
*/
#define LCS_NAMELEN 32

/*
Boundary Condition Flags
*/
#define LCS_INTERNAL 0
#define LCS_INFLOW 1
#define LCS_OUTFLOW 2
#define LCS_WALL 3
#define LCS_SLIP 4
#define LCS_2D 5
#define LCS_MASK 6

/*
Define the different types of LCS diagnostics here
*/
#define	FTLE_FWD 0
#define	FTLE_BKWD 1
#define	LP_TRACER 2

/*
Integration methods:
*/
#define	EULER 0
#define	TRAPEZOIDAL 1
#define	RK2 2
#define	RK3 3
#define	RK4 4

/*
Interpolation methods
*/
#define	NEAREST_NBR 0
#define	LINEAR 1
#define	QUADRATIC 2
#define	CUBIC 3
#define	IDW 4
#define	TSE 5

/*
Interface function prototypes:
These are available to the user.  Could put this
in the header, but I worry about multiple defines of mpi.h
*/
void cfd2lcs_init_c(
	MPI_Comm usercomm,
	int n[3],
	int offset[3],
	lcsdata_t *x,
	lcsdata_t *y,
	lcsdata_t *z,
	int *flag
);

int cfd2lcs_diagnostic_init_c(
	int lcs_type,
	int resolution,
	lcsdata_t t,
	lcsdata_t h,
	lcsdata_t rhop,
	lcsdata_t dp,
	char label[]
);

void cfd2lcs_update_c(
	int n[3],
	lcsdata_t *u,
	lcsdata_t *v,
	lcsdata_t *w,
	lcsdata_t time,
	lcsdata_t cfl
);

void cfd2lcs_diagnostic_destroy_c(
	int lcs_handle
);

void cfd2lcs_set_option_c(
	char option[],
	int val
);
