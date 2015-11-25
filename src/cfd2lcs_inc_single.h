/*
CFD2LCS C/C++ Include file.
Contains interface definitions needed for user-level API
To use:  #include "cfd2lcs_inc.h"
*/
#include <string.h>

/*
Set the precision, and be sure that you match the precision
defined in cfd2lcs_inc.f90 (eventually we could have a compiler macro or something...)
*/
typedef float lcsdata_t;

/*
The string namelength
*/
#define LCS_NAMELEN 32

/*
Boundary Condition Flags
*/
#define LCS_PERIODIC 0
#define LCS_INFLOW 1
#define LCS_OUTFLOW 2
#define LCS_WALL 3
#define LCS_SLIP 4

/*
Define the different types of LCS diagnostics here
*/
#define	FTLE_FWD 0
#define	FTLE_BKWD 1
#define	LP_TRACER 2

/*
Define the different possible data layouts that we can accept.
These can be added as needed...
*/
#define LCS_3V 0
#define LCS_1V_INTERLACED 1



/*****************************************
Interface functions
******************************************/
//Prototypes for unctions to be called from C routine:
void cfd2lcs_init_c(int usercomm,int n[3], int offset[3], lcsdata_t *x,
	 lcsdata_t *y, lcsdata_t *z, int BC_LIST[6], lcsdata_t lperiodic[3], int datalayout);
void cfd2lcs_diagnostic_init_c(int lcs_type, int resolution, 
	lcsdata_t T, lcsdata_t h, lcsdata_t rhop, lcsdata_t dp, char label[]);
void cfd2lcs_update_c(int n[3],	lcsdata_t *u, lcsdata_t *v, lcsdata_t *w, 
	lcsdata_t time, int datalayout);
//Prototypes for F90 library functions
void cfd2lcs_init_(int *usercomm, int *n, int *offset, 
	float *x, float *y, float *z, int *BC_LIST, float *lperiodic);
void cfd2lcs_diagnostic_init_(int *lcs_handle, int *lcs_type, int *resolution,
	 float *T, float *h, float *rhop, float *dp, char *label, int namelen);
void cfd2lcs_update_(int *n,float *u, float *v, float *w, float *time);



/*
The default intialization.
We can eventually handle different data packing for grid coordinates
*/
void cfd2lcs_init_c(int usercomm,int n[3], int offset[3], 
	lcsdata_t *x, lcsdata_t *y, lcsdata_t *z, int BC_LIST[6], lcsdata_t lperiodic[3], int datalayout)
{
	switch(datalayout)
	{
		case LCS_3V:
			//We just pass the 3 x,y,z vectors "as is" into the cfd2lcs arrays
			cfd2lcs_init_(&usercomm,&n[0],&offset[0],x,y,z,&BC_LIST[0],&lperiodic[0]);
			break;		

		case LCS_1V_INTERLACED:
			//We will need to un-sort the data from 1 vector into separate x,y,z vectors
			//before passing.  But I think it is best to wait until Andrew and I can 
			//discuss how his data is stored before implementing this...
			printf("ERROR:  LCS_1V_INTERLACED not yet supported\n");
			break;		
	
		default:
			printf("ERROR:  Bad datalayout: %d\n",datalayout);
			break;		
	}		
}


/*
Initialize an LCS diagnositc:
Andrew, is it OK to return an integer here for the lcs_handle argument?
*/
void cfd2lcs_diagnostic_init_c(int lcs_type, int resolution, lcsdata_t T, lcsdata_t h, lcsdata_t rhop, lcsdata_t dp, char label[])
{
	int lcs_handle;
	int namelen = strlen(label); // need to pass this explicitly to f90 
	cfd2lcs_diagnostic_init_(&lcs_handle,&lcs_type,&resolution,&T,&h,&rhop,&dp,&label[0],namelen);
}


/*
Update the LCS diagnostics
Similarly to the init, we can allow for multiple data layouts here...
*/
void cfd2lcs_update_c(int n[3],	lcsdata_t *u, lcsdata_t *v, lcsdata_t *w, lcsdata_t time, int datalayout)
{
	switch(datalayout)
	{
		case LCS_3V:
			//We just pass the 3 x,y,z vectors "as is" into the cfd2lcs arrays
			cfd2lcs_update_(&n[0],u,v,w,&time);
			break;		

		case LCS_1V_INTERLACED:
			//We will need to un-sort the data from 1 vector into separate x,y,z vectors
			//before passing.  But I think it is best to wait until Andrew and I can 
			//discuss how his data is stored before implementing this...
			printf("ERROR:  LCS_1V_INTERLACED not yet supported\n");
			break;		
	
		default:
			printf("ERROR:  Bad datalayout: %d\n",datalayout);
			break;		
	}		
}
