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

	!Structured CFD data (scfd_t):
	type scfd_t
		character(len=MAX_NAMELEN):: label

		!Dimensions
		integer:: ni, nj, nk, ng
		integer:: gni, gnj, gnk
		integer:: offset_i, offset_j, offset_k

		!Boundary condition flags:
		integer:: bc_list(6)

		!Communicators
		type(scomm_t):: scomm_face_r0
		type(scomm_t):: scomm_max_r0
		type(scomm_t):: scomm_face_r1
		type(scomm_t):: scomm_max_r1
		type(scomm_t):: scomm_face_r2
		type(scomm_t):: scomm_max_r2

		!Data
		type(sr1_t):: grid  !Cartesian grid coordinates
		type(sr1_t):: u		!Current velocity field

	end type scfd_t

	!Lagrangian Particles
	type lp_t
		character(len=MAX_NAMELEN):: label
		integer:: np
		real(LCSRP):: dp  !Diameter
		real(LCSRP):: rhop !Density
		type(ur1_t):: xp !Position
		type(ur1_t):: up !Velocity
	end type lp_t


	!lcs:
	type lcs_t
		character(len=MAX_NAMELEN):: label
		integer:: diagnostic

		real(LCSRP):: T ! Integration time
		real(LCSRP):: h ! Visualization timestep

		type(scfd_t),pointer :: flow
		type(lp_t),pointer:: lp
		type(sr1_t),pointer:: fm_fwd
		type(sr1_t),pointer:: fm_bkwd
	end type lcs_t


	!The CFD side data:
	type(scfd_t):: scfd

	!The LCS diagnostics:
	integer:: NLCS
	type(lcs_t),allocatable:: lcs(:)

	!Error handling:
	integer:: CFD2LCS_ERROR

	contains



end module data_m
