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
	lcsdata_t *rhop,
	lcsdata_t *dp,
	char *label,
	int namelen
);

void cfd2lcs_update_(
	int *n,
	lcsdata_t *u,
	lcsdata_t *v,
	lcsdata_t *w,
	lcsdata_t *time,
	lcsdata_t *cfl
);

void cfd2lcs_diagnostic_destroy_(
	int *lcs_handle
);

void cfd2lcs_set_option_(
	char *option,
	int *val,
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
	void *x,
	void *y,
	void *z,
	void *flag,
	int datastride
)
{
	//Convert mpicomm to Fortran type:
	MPI_Fint fortran_usercomm;
	fortran_usercomm=MPI_Comm_c2f(usercomm);

	switch(datastride)
	{
		case LCS_3V:
			//We just pass the 3 x,y,z vectors "as is" into the cfd2lcs arrays
			cfd2lcs_init_(&fortran_usercomm,&n[0],&offset[0],x,y,z,flag);
			break;
		default:
			//JRF:  Is another treatment needed here???
			printf("ERROR:  Bad datastride: %d\n",datastride);
			break;
	}
}

//Initialization of an lcs diagnostic:
int cfd2lcs_diagnostic_init_c(
	int lcs_type,
	int resolution,
	lcsdata_t t,
	lcsdata_t h,
	lcsdata_t rhop,
	lcsdata_t dp,
	char label[]
)
{
	int lcs_handle;
	int namelen = strlen(label); // need to pass this explicitly to f90
	cfd2lcs_diagnostic_init_(&lcs_handle,&lcs_type,&resolution,&t,&h,&rhop,&dp,&label[0],namelen);
	//Pass back the handle
	return lcs_handle;
}

//Update of lcs diagnostics
void cfd2lcs_update_c(
	int n[3],
	void *u,
	void *v,
	void *w,
	lcsdata_t time,
	lcsdata_t cfl,
	int datastride
)
{
	lcsdata_t *uu,*vv,*ww;
	int i;
	int databytes = sizeof(lcsdata_t);

	switch(datastride)
	{
		case LCS_3V:
			//We just pass the 3 x,y,z vectors "as is" into the cfd2lcs arrays
			cfd2lcs_update_(&n[0],u,v,w,&time,&cfl);
			break;
		case LCS_1V_INTERLACED:
			//Assume user passes uvw into the "u" array
			//We want to separate this into three arrays, uu,vv,ww
			//to be passed to cfd2lcs_update.
			if(u==NULL) return;
			uu=malloc(n[0]*n[1]*n[2]*sizeof(lcsdata_t));
			vv=malloc(n[0]*n[1]*n[2]*sizeof(lcsdata_t));
			ww=malloc(n[0]*n[1]*n[2]*sizeof(lcsdata_t));

			for(i=0;i<n[0]*n[1]*n[2];i++){
			if(databytes==8)
				uu[i]=(lcsdata_t) (((double *) u)[i*datastride]);
			else if(databytes==4)
				uu[i]=(lcsdata_t) (((float *) u)[i*datastride]);
			else
				return;
			}
			for(i=0;i<n[0]*n[1]*n[2];i++){
				if(databytes==8)
					vv[i]=(lcsdata_t) (((double *) u)[(i+1)*datastride]);
				else if(databytes==4)
					vv[i]=(lcsdata_t) (((float *) u)[(i+1)*datastride]);
				else
					return;
			}
			for(i=0;i<n[0]*n[1]*n[2];i++){
				if(databytes==8)
					ww[i]=(lcsdata_t) (((double *) u)[(i+2)*datastride]);
				else if(databytes==4)
					ww[i]=(lcsdata_t) (((float *) u)[(i+2)*datastride]);
				else
					return;
			}
			cfd2lcs_update_(&n[0],uu,vv,ww,&time,&cfl);

		default:
			printf("ERROR:  Unrecongnized datastride\n");
			break;
	}
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
){
	int namelen = strlen(option); // need to pass this explicitly to f90
	cfd2lcs_set_option_(&option[0],&val,namelen);  //note, namelen should be last argument.
}
