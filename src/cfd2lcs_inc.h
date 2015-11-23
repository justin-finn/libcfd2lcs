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
//#define DOUBLE_PRECISION   // Uncomment this line for double precision
#ifdef DOUBLE_PRECISION
typedef double lcsdata_t;
#else
typedef float lcsdata_t;
#endif

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

/*
The default intialization.
We can eventually handle different data packing for grid coordinates
*/
void cfd2lcs_init_c(int usercomm,int n[3], int offset[3], 
	float *x, float *y, float *z, int BC_LIST[6], float lperiodic[3], int datalayout)
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
void cfd2lcs_diagnostic_init_c(int lcs_type, int resolution, float T, float h, float rhop, float dp, char label[])
{
	int lcs_handle;
	int namelen = strlen(label); // need to pass this explicitly to f90 
	cfd2lcs_diagnostic_init_(&lcs_handle,&lcs_type,&resolution,&T,&h,&rhop,&dp,&label[0],namelen);
}


/*
Update the LCS diagnostics
Similarly to the init, we can allow for multiple data layouts here...
*/
void cfd2lcs_update_c(int n[3],	float *u, float *v, float *w, float time, int datalayout)
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
