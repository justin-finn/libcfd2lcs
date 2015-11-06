!
! Data structure definitions
!
module data_m
	implicit none
	!----
	INCLUDE 'cfd2lcs_f.h'
	INCLUDE 'mpif.h'
	!----

	!Some constants:
	integer,parameter:: MAX_NAMELEN = 32
	integer,parameter:: NGHOST_CFD = 1

	!Particle flags:
	integer,parameter:: &
			LP_IB = -1, &
			LP_UNKNOWN = -2, &
			LP_RECYCLE = -3, &
			LP_STICK = -4

	!Interpolation options
	integer,parameter::&
		ZERO_ORDER = 0, &
		GAUSSIAN_RBF = 1, &
		TRICUBIC = 2, &
		TRILINEAR = 3
	integer:: INTERPOLATION = GAUSSIAN_RBF
	!integer:: INTERPOLATION = ZERO_ORDER

	!MPI stuff:
	integer :: nprocs,lcsrank,lcscomm
	integer:: MPI_LCSRP  !Real precision for mpi

	!Real, structured rank 0 scalar (sr0_t)
	type sr0_t
		integer:: ni,nj,nk,ng
		character(len=MAX_NAMELEN):: label
		real(LCSRP), allocatable:: r(:,:,:)
	end type sr0_t

	!Real, structured Cartesian rank 1 vector (sr1_t)
	type sr1_t
		integer:: ni,nj,nk,ng
		character(len=MAX_NAMELEN):: label
		real(LCSRP), allocatable:: x(:,:,:)
		real(LCSRP), allocatable:: y(:,:,:)
		real(LCSRP), allocatable:: z(:,:,:)
		logical:: periodic_translate  !For translating x,y,z coordinates over periodic boundaries
	end type sr1_t

	!Real, structured Cartesian rank 2 tensor (sr2_t)
	type sr2_t
		integer:: ni,nj,nk,ng
		character(len=MAX_NAMELEN):: label
		real(LCSRP), allocatable:: xx(:,:,:)
		real(LCSRP), allocatable:: xy(:,:,:)
		real(LCSRP), allocatable:: xz(:,:,:)
		real(LCSRP), allocatable:: yx(:,:,:)
		real(LCSRP), allocatable:: yy(:,:,:)
		real(LCSRP), allocatable:: yz(:,:,:)
		real(LCSRP), allocatable:: zx(:,:,:)
		real(LCSRP), allocatable:: zy(:,:,:)
		real(LCSRP), allocatable:: zz(:,:,:)
	end type sr2_t

	!Real, unstructured rank 0 scalar (ur0_t)
	type ur0_t
		integer:: n
		character(len=MAX_NAMELEN):: label
		real(LCSRP), allocatable:: r(:)
	end type ur0_t

	!Real, unstructured rank 1 vector (ur1_t)
	type ur1_t
		integer:: n
		character(len=MAX_NAMELEN):: label
		real(LCSRP), allocatable:: x(:)
		real(LCSRP), allocatable:: y(:)
		real(LCSRP), allocatable:: z(:)
	end type ur1_t

	!Real, unstructured rank 2 tensor (ur2_t)
	type ur2_t
		integer:: n
		character(len=MAX_NAMELEN):: label
		real(LCSRP), allocatable:: xx(:)
		real(LCSRP), allocatable:: xy(:)
		real(LCSRP), allocatable:: xz(:)
		real(LCSRP), allocatable:: yx(:)
		real(LCSRP), allocatable:: yy(:)
		real(LCSRP), allocatable:: yz(:)
		real(LCSRP), allocatable:: zx(:)
		real(LCSRP), allocatable:: zy(:)
		real(LCSRP), allocatable:: zz(:)
	end type ur2_t

	!Integer, unstructured rank 0 scalar (ui0_t)
	type ui0_t
		integer:: n
		character(len=MAX_NAMELEN):: label
		integer(LCSIP), allocatable:: i(:)
	end type ui0_t

	!Integer, unstructured rank 1 vector (ui1_t)
	type ui1_t
		integer:: n
		character(len=MAX_NAMELEN):: label
		integer(LCSIP), allocatable:: x(:)
		integer(LCSIP), allocatable:: y(:)
		integer(LCSIP), allocatable:: z(:)
	end type ui1_t


	!Structured Comms (scomm_t):
	type scomm_t
		character(len=MAX_NAMELEN):: label
		integer:: connectivity !Either FACE_CONN or MAX_CONN
		integer:: datatype !Either R0_COMM, R1_COMM, R2_COMM
		integer :: nbr_rank(-1:1,-1:1,-1:1) !Rank of processor in each direction
		integer :: flag(-1:1,-1:1,-1:1) !tells us what to do in each direction
		integer :: checker(-1:1,-1:1,-1:1) !Checkerboard pattern for each comm direction
		integer :: pack_start(-1:1,-1:1,-1:1), unpack_start(-1:1,-1:1,-1:1)  !First index into pack/unpack buffers
		integer :: n_pack(-1:1,-1:1,-1:1), n_unpack(-1:1,-1:1,-1:1)  !Number of pack/unpack elements in each dir
		integer :: pack_list_min(-1:1,-1:1,-1:1,1:3),pack_list_max(-1:1,-1:1,-1:1,1:3) !3D Range for packing buffers
		integer :: unpack_list_min(-1:1,-1:1,-1:1,1:3),unpack_list_max(-1:1,-1:1,-1:1,1:3) !3D Range for unpacking buffers
		real(LCSRP), allocatable :: pack_buffer(:), unpack_buffer(:)  !Exchange buffers
		integer:: pack_bufsize,unpack_bufsize
		real(LCSRP):: periodic_shift(-1:1,-1:1,-1:1,1:3)  !For shifting coordinates across periodic boundaries
	end type scomm_t

	!Structured grid (sgrid_t):
	type sgrid_t
		character(len=MAX_NAMELEN):: label

		!Dimensions
		integer:: ni, nj, nk, ng
		integer:: gni, gnj, gnk
		integer:: offset_i, offset_j, offset_k

		!Boundary condition flags:
		integer:: global_bc_list(6) !for external domain boundary
		integer:: bc_list(6)  !for processor boundary
		real(LCSRP):: lperiodic(3)

		!Communicators
		type(scomm_t):: scomm_face_r0
		type(scomm_t):: scomm_max_r0
		type(scomm_t):: scomm_face_r1
		type(scomm_t):: scomm_max_r1
		type(scomm_t):: scomm_face_r2
		type(scomm_t):: scomm_max_r2

		!Data
		type(sr1_t):: grid  !Cartesian grid coordinates
	end type sgrid_t
	integer:: NSGRID
	type(sgrid_t),allocatable,target:: sgrid_c(:) !Collection of NSGRID sgrid structures

	!Lagrangian Particles
	type lp_t
		character(len=MAX_NAMELEN):: label
		integer:: np !Number of particles (on each proc)
		real(LCSRP):: dp  !Diameter
		real(LCSRP):: rhop !Density
		type(ur1_t):: xp !Position
		type(ur1_t):: up !Velocity
		type(ui1_t):: no !Nearest node index (i,j,k)
		type(ui0_t):: proc0 !Origin proc
		type(ui0_t):: no0 !Origin node
		type(ui0_t):: flag !Multipurpose flag
	end type lp_t
	integer:: NLP
	type(lp_t),allocatable,target:: lp_c(:)  !Collection of NLP lp structures

	!lcs:
	type lcs_t
		integer(LCSIP):: id  !A unique integer identifier to allow the user to modify aspects of each lcs
		character(len=MAX_NAMELEN):: label
		integer(LCSIP):: diagnostic
		real(LCSRP):: T ! Integration time
		real(LCSRP):: h ! Visualization timestep
		type(sgrid_t),pointer :: sgrid  !pointer to the structured grid
		type(lp_t),pointer:: lp !pointer to Lagrangian particles if used
		type(sr1_t):: fm  !Flow Map
	end type lcs_t
	integer:: NLCS
	type(lcs_t),allocatable,target:: lcs_c(:)  !Collection of NLCS lcs structures

	!The CFD side data:
	type scfd_t
		character(len=MAX_NAMELEN):: label
		real(LCSRP):: t_n, t_np1
		type(sgrid_t),pointer:: sgrid
		type(sr1_t):: u_n	!old velocity field
		type(sr1_t):: u_np1  !latest velocity field
	end type scfd_t
	type(scfd_t):: scfd

	!Error handling:
	integer:: CFD2LCS_ERROR


	contains
	!These function set the rules for the relationship
	!between 3D (i,j,k) indices and 1D (l) index for structured grids.
	integer function ijk2l(i,j,k,ni,nj)
		implicit none
		integer::i,j,k !current index
		integer::ni,nj !num pts in i,j directions
		ijk2l= (k-1)*ni*nj + (j-1)*ni +i
	end function ijk2l
	!JRF: Return  i,j,k, from l:
	integer function l2i(l,ni)
		implicit none
		integer::l !current index
		integer::ni !num pts in i direction
	    l2i= mod(l,ni)
	end function l2i
	integer function l2j(l,ni,nj)
		implicit none
		integer::l !current index
		integer::ni,nj !num pts in i,j directions
		l2j = mod((l)/ni,nj)
	end function l2j
	integer function l2k(l,ni,nj)
		implicit none
		integer::l !current index
		integer::ni,nj !num pts in i,j directions
		l2k = (l)/(ni*nj)
	end function l2k




end module data_m
