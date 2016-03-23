/*
A Simple working example to show how to call
CFD2LCS.  Subroutines starting with "your_"
are independent of the cfd2lcs functionality,
and are used only to create a simple dataset
for this example.

FLOW FIELD:  2D, time dependent double gyre flow
in the X-Y plane. User parameters included directly below.
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

//Uncomment only one of the following:
#include "cfd2lcs_inc_sp.h"  	// For Single Precision
//#include "cfd2lcs_inc_dp.h"  	// For Double Precision

////////////////////////////
///////BEGIN USER INPUT///////
////////////////////////////
//-----
//Domain dimensions
//-----
lcsdata_t LX = 2.0*M_PI;
lcsdata_t LY = 2.0*M_PI;
lcsdata_t LZ = 2.0*M_PI;
//-----
//Total number of grid points in each direction
//-----
int NX = 64;
int NY = 64;
int NZ = 64;
//-----
//"Simulation" parameters
//-----
lcsdata_t DT = 0.025;
lcsdata_t START_TIME = 0.0;
lcsdata_t END_TIME = 10.001;
lcsdata_t T = 10.0;
lcsdata_t H = 1.0;
int RESOLUTION = 0;
lcsdata_t CFL = 0.9;
//-----
//"ABC" parameters
//-----
#define ABC_A sqrt(3.0)
#define ABC_B sqrt(2.0)
#define ABC_C 1.0
#define ABC_D 0.0
//----
////////////////////////////
///////END USER INPUT///////
////////////////////////////

//
//Some other global variables
//
int nprocs,mycomm;
int nproc_x,nproc_y,nproc_z;
int myrank, myrank_i, myrank_j, myrank_k;
int ni, nj, nk, offset_i, offset_j, offset_k;
lcsdata_t *x, *y, *z, *u, *v, *w; //grid coordinates and velocities
int *flag; //boundary condition flag
int id_fwd,id_bkwd;

//-----
//These function set the rules for the relationship
//between 3D (i,j,k) indices of each processor and the processor rank.
//-----
int  rank2i(int rank,int npi)
{
	return((rank % npi)+1);
}
int  rank2j(int rank,int npi,int npj)
{
	return((rank/npi % npj)+1);
}
int  rank2k(int rank,int npi,int npj)
{
	return(rank/(npi*npj)+1);
}

//----
//Partition the domain into chunks for each processor
//based on the number of processors specified in X,Y,Z direction.  Unequal sizes are permitted.
//Actual grid coordinates are set in your_grid_function().
//----
void your_partition_function()
{
	//----
	int proc, guess_ni, guess_nj, guess_nk;
	int leftover_ni, leftover_nj, leftover_nk;

		if (myrank ==0){
			printf("In your_partition_function...\n");
		}

		//Check the number of procs:
		if (nprocs != nproc_x*nproc_y*nproc_z){
			if(myrank==0)printf("ERROR: Number of processors incompatible with partition\n");
			exit(1);
		}

		if (nproc_x > NX || nproc_y > NY || nproc_z > NZ){
			if(myrank==0)printf("ERROR: More processors than grid points...\n");
			exit(1);
		}

		//Set local proc coordinates:
		myrank_i = rank2i(myrank,nproc_x);
		myrank_j = rank2j(myrank,nproc_x,nproc_y);
		myrank_k = rank2k(myrank,nproc_x,nproc_y);

		//Set the array range on each processor, allowing for unequal sizes:
		guess_ni = NX/nproc_x;
		guess_nj = NY/nproc_y;
		guess_nk = NZ/nproc_z;
		leftover_ni = NX - (guess_ni*nproc_x);
		leftover_nj = NY - (guess_nj*nproc_y);
		leftover_nk = NZ - (guess_nk*nproc_z);
		if (myrank_i <=leftover_ni)
		{
			ni = guess_ni + 1;
			offset_i = (guess_ni+1)*(myrank_i-1);
		}
		else
		{
			ni = guess_ni;
			offset_i = (guess_ni+1)*(leftover_ni) + guess_ni*(myrank_i-leftover_ni-1);
		}

		if (myrank_j <= leftover_nj)
		{
			nj = guess_nj + 1;
			offset_j = (guess_nj+1)*(myrank_j-1);
		}
		else
		{
			nj = guess_nj;
			offset_j = (guess_nj+1)*(leftover_nj) + guess_nj*(myrank_j-leftover_nj-1);
		}

		if (myrank_k <= leftover_nk)
		{
			nk = guess_nk + 1;
			offset_k = (guess_nk+1)*(myrank_k-1);
		}
		else
		{
			nk = guess_nk;
			offset_k = (guess_nk+1)*(leftover_nk) + guess_nk*(myrank_k-leftover_nk-1);
		}
}

//----
//Set the grid coordinates x,y,z.
//Just use a uniform Cartesian grid for now, arbitrary spacing is possible.
//----
void your_grid_function()
{
	//----
	int i,j,k,ii,jj,kk;
	lcsdata_t dx,dy,dz;
	//----
	#define MAX(x, y) (((x) > (y)) ? (x) : (y))
	#define MIN(x, y) (((x) < (y)) ? (x) : (y))

	if (myrank ==0) printf("in your_grid_function...\n");

	dx = LX / (lcsdata_t)(MAX(NX,1));
	dy = LY / (lcsdata_t)(MAX(NY,1));
	dz = LZ / (lcsdata_t)(MAX(NZ,1));
	int	ind = 0;
	for( k = offset_k; k < offset_k+nk; k++)
	{
		for( j = offset_j; j < offset_j+nj; j++)
		{
			for( i = offset_i; i < offset_i+ni; i++)
			{
				x[ind] = 0.5*dx + (lcsdata_t)(i)*dx;
				y[ind] = 0.5*dy + (lcsdata_t)(j)*dy;
				z[ind] = 0.5*dz + (lcsdata_t)(k)*dz;
				flag[ind] = LCS_INTERNAL;  //Set  LCS_INTERNAL everywhere for tripply periodic
				ind++;
			}
		}
	}
}


//
//Set the velocity field
//
void set_velocity(lcsdata_t time)
{
	//----
	int i,j,k;
	lcsdata_t aoft,boft;
	lcsdata_t fofxt;
	lcsdata_t dfdx;
	//----

	if (myrank ==0) printf("in set_velocity...\n");

	//-----
	// Assumes periodicity of multiple 2pi in X,Y,Z
	//-----
	int	ind = 0;
	for( k = 0; k < nk; k++)
	{
		for( j = 0; j < nj; j++)
		{
			for( i = 0; i < ni; i++)
			{
				u[ind] = (ABC_A+ABC_D*sin(time))*sin(z[ind]) + ABC_C*cos(y[ind]);
				v[ind] = ABC_B*sin(x[ind]) + (ABC_A+ABC_D*sin(time))*cos(z[ind]);
				w[ind] = ABC_C*sin(y[ind]) + ABC_B*cos(x[ind]);
				ind ++;
			}
		}
	}
}


int main (argc, argv)
	int argc;
	char *argv[];
{
	int timestep;
	lcsdata_t time;
	int ierr;

	//-----
	//Initialize MPI
	//-----
	MPI_Init (&argc, &argv);	/* starts MPI */
	mycomm = MPI_COMM_WORLD;
	MPI_Comm_rank (mycomm, &myrank);	/* get current process id */
	MPI_Comm_size (mycomm, &nprocs);	/* get number of processes */

	//-----
	//Check the input on Rank 0 and broadcast input to everyone
	//Assume we want the last 3 arguments for NPX,NPY,NPZ
	//-----
	if(myrank==0)
	{
		if( argc >= 4 ) {
			nproc_x = atoi(argv[argc-3]);
			nproc_y = atoi(argv[argc-2]);
			nproc_z = atoi(argv[argc-1]);
		}
		else {
			printf("Error: must supply arguments for NPROCS_X, NPROCS_Y, NPROCS_Z\n");
			printf("Example usage:  mpirun -np 8 ./CFD2LCS_TEST 4 2 1\n");
		}
		printf("Will partition domain using: %d %d %d sub-domains\n",nproc_x,nproc_y,nproc_z);
	}
	ierr = MPI_Bcast(&nproc_x, 1,MPI_INTEGER,0,mycomm);
	ierr = MPI_Bcast(&nproc_y, 1,MPI_INTEGER,0,mycomm);
	ierr = MPI_Bcast(&nproc_z, 1,MPI_INTEGER,0,mycomm);


	//-----
	//Setup the parallel partition:
	//-----
	your_partition_function();

	//-----
	//Allocate space for data
	//-----
	x = malloc(ni*nj*nk*sizeof(lcsdata_t));
	y = malloc(ni*nj*nk*sizeof(lcsdata_t));
	z = malloc(ni*nj*nk*sizeof(lcsdata_t));
	u = malloc(ni*nj*nk*sizeof(lcsdata_t));
	v = malloc(ni*nj*nk*sizeof(lcsdata_t));
	w = malloc(ni*nj*nk*sizeof(lcsdata_t));
	flag = malloc(ni*nj*nk*sizeof(int));

	//-----
	//Set the grid points:
	//-----
	your_grid_function();

	//-----
	//Initialize cfd2lcs for your data
	//-----
	int n[3] = {ni,nj,nk};// number of grid points for THIS partition
	int offset[3]= {offset_i,offset_j,offset_k};// Global offset for these grid points
	cfd2lcs_init_c(mycomm,n,offset,x,y,z,flag);

	//-----
	//Initialize LCS diagnostics
	//-----
	char labelfwd[LCS_NAMELEN]="fwdFTLE";
	id_fwd = cfd2lcs_diagnostic_init_c(FTLE_FWD,RESOLUTION,T,H,labelfwd);
	char labelbkwd[LCS_NAMELEN]="bkwdFTLE";
	id_bkwd = cfd2lcs_diagnostic_init_c(FTLE_BKWD,RESOLUTION,T,H,labelbkwd);
	printf("fwd and bkwd id %d %d \n",id_fwd, id_bkwd);

	//-----
	//Set CFD2LCS options/parameters
	//-----
	cfd2lcs_set_option_c("SYNCTIMER",LCS_FALSE);
	cfd2lcs_set_option_c("DEBUG",LCS_FALSE);
	cfd2lcs_set_option_c("WRITE_FLOWMAP",LCS_FALSE);
	cfd2lcs_set_option_c("WRITE_BCFLAG",LCS_FALSE);
	cfd2lcs_set_option_c("INCOMPRESSIBLE",LCS_FALSE);
	cfd2lcs_set_option_c("AUX_GRID",LCS_FALSE);
	cfd2lcs_set_option_c("INTEGRATOR",RK3);
	cfd2lcs_set_option_c("INTERPOLATOR",LINEAR);
	cfd2lcs_set_param_c("CFL", CFL);

	//-----
	//***Start of your flow solver timestepping loop***
	//-----
	timestep = 0;
	time = START_TIME;
	set_velocity(time);
	do
	{
		if(myrank == 0)
		{
			printf("------------------------------------------------------------------\n");
			printf("STARTING TIMESTEP #%d: time = %f, DT= %f\n",timestep,time,DT);
			printf("------------------------------------------------------------------\n");
		}
		set_velocity(time);// !CFD Solve for the velocity field
		cfd2lcs_update_c(n,u,v,w,time);  //Update LCS fields
		timestep = timestep + 1;
		time = time + DT;
	}
	while (time <= END_TIME);

	//-----
	//Cleanup
	//-----
	free(x);
	free(y);
	free(z);
	free(u);
	free(v);
	free(w);
	free(flag);
	MPI_Finalize();
	return 0;
}
