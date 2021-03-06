/*
Copyright (C) 2015-2016, Justin R. Finn.  All rights reserved.
libcfd2lcs is distributed is under the terms of the GNU General Public License
*/

/*
CFD2LCS C/C++ Interface file.
*/
#include <stdio.h>
#include <string.h>
#include <cfd2lcs_inc.h>

/*
Prototypes for CFD2LCS entry points (C calls to F90 library)
*/
void cfd2lcs_init_(
	int *usercomm,
	int *n,
	int *offset,
	lcsdata_t *x,
	lcsdata_t *y,
	lcsdata_t *z,
	int *flag
);

int cfd2lcs_diagnostic_init_(
	int *lcs_handle,
	int *lcs_type,
	int *resolution,
	lcsdata_t *t,
	lcsdata_t *h,
	char *label,
	int namelen
);

void cfd2lcs_update_(
	int *n,
	lcsdata_t *u,
	lcsdata_t *v,
	lcsdata_t *w,
	lcsdata_t *time
);

void cfd2lcs_diagnostic_destroy_(
	int *lcs_handle
);

void cfd2lcs_set_option_(
	char *option,
	int *val,
	int namelen
);

void cfd2lcs_set_param_(
	char *option,
	lcsdata_t *val,
	int namelen
);



/*
MAIN CFD2LCS C/C++ INTERFACE:
*/

//The default intialization.
void cfd2lcs_init_c(
	MPI_Comm usercomm,
	int n[3],
	int offset[3],
	lcsdata_t *x,
	lcsdata_t *y,
	lcsdata_t *z,
	int *flag
)
{
	//Convert mpicomm to Fortran type:
	MPI_Fint fortran_usercomm;
	fortran_usercomm=MPI_Comm_c2f(usercomm);

	//We just pass the 3 x,y,z vectors "as is" into the cfd2lcs arrays
	cfd2lcs_init_(&fortran_usercomm,&n[0],&offset[0],x,y,z,flag);
}

//Initialization of an lcs diagnostic:
int cfd2lcs_diagnostic_init_c(
	int lcs_type,
	int resolution,
	lcsdata_t t,
	lcsdata_t h,
	char label[]
)
{
	int lcs_handle;
	int namelen = strlen(label); // need to pass this explicitly to f90
	cfd2lcs_diagnostic_init_(&lcs_handle,&lcs_type,&resolution,&t,&h,&label[0],namelen);
	//Pass back the handle
	return lcs_handle;
}

//Update of lcs diagnostics
void cfd2lcs_update_c(
	int n[3],
	lcsdata_t *u,
	lcsdata_t *v,
	lcsdata_t *w,
	lcsdata_t time
)
{
	lcsdata_t *uu,*vv,*ww;
	int i;

	//We just pass the 3 x,y,z vectors "as is" into the cfd2lcs arrays
	cfd2lcs_update_(&n[0],u,v,w,&time);
}

//Destroy an existing lcs diagnostic
void cfd2lcs_diagnostic_destroy_c(
	int lcs_handle
){
	cfd2lcs_diagnostic_destroy_(&lcs_handle);
}

//Set a cfd2lcs user option
void cfd2lcs_set_option_c(
	char option[],
	int val
)
{
	int namelen = strlen(option); // need to pass this explicitly to f90
	cfd2lcs_set_option_(&option[0],&val,namelen);  //note, namelen should be last argument.
}

//Set a cfd2lcs user parameter
void cfd2lcs_set_param_c(
	char param[],
	lcsdata_t val
)
{
	int namelen = strlen(param); // need to pass this explicitly to f90
	cfd2lcs_set_param_(&param[0],&val,namelen);  //note, namelen should be last argument.
}
