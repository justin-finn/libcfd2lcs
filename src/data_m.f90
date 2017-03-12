!
!Copyright (C) 2015-2016, Justin R. Finn.  All rights reserved.
!libcfd2lcs is distributed is under the terms of the GNU General Public License
!
!
! Data structure definitions
!
module data_m
      implicit none
      !----
      INCLUDE 'cfd2lcs_inc.f90'
      INCLUDE 'mpif.h'
      !----

      !!!!!!! USER ACCESSIBLE OPTIONS !!!!!!!!
      !
      ! Set the verbosity of output
      !
      logical:: LCS_VERBOSE = .FALSE.

      !
      ! Syncronize the cpu timer (with possible slowdown for mpi_barrier calls)
      !
      logical:: SYNC_TIMER = .FALSE.

      !
      ! Write all the sgrids during initialization:
      !
      logical:: DEBUG_SGRID = .TRUE.

      !
      ! Zero any negative FTLE values (assuming incompressible flow)
      !
      logical:: INCOMPRESSIBLE = .FALSE.

      !
      ! Auxillary grid method (See Farazmand & Haller, 2012)
      !
      logical:: AUX_GRID = .FALSE.

      !
      ! Write the flowmap along with the lcs diagnostic:
      !
      logical:: FLOWMAP_IO = .TRUE.
      
      !
      ! Write the velocity at T_0 along with the lcs diagnostic:
      !
      logical:: VELOCITY_IO = .TRUE.
      
      !
      ! Write the bcflag along with the lcs diagnostic:
      !
      logical:: BCFLAG_IO = .FALSE.

      !
      ! Particle integration scheme
      !
      integer:: INTEGRATOR = RK2

      !
      ! Particle interpolation scheme
      !
      integer:: INTERPOLATOR = LINEAR
      
      !
      ! Number of particles to integrate before updating the user
      !
      integer:: N_UPDATE = 10000000
      

      !!!!!!! USER ACCESSIBLE PARAMS: !!!!!!!!
      
      !The cfl number used for particle integration (Fwd & Bkwd)
      !This gets passed by the user in cfd2lcs_update
      real(LCSRP):: CFL_MAX = 0.5_LCSRP

      !
      !For injecting a blob of tracers:
      !
      real(LCSRP):: TRACER_INJECT_X = huge(1.0_LCSRP)
      real(LCSRP):: TRACER_INJECT_Y = huge(1.0_LCSRP)
      real(LCSRP):: TRACER_INJECT_Z = huge(1.0_LCSRP)
      real(LCSRP):: TRACER_INJECT_RADIUS = 0.0_LCSRP
      
      !Spacing of the auxillary grid (relative to the lcs diagnostic grid):
      real(LCSRP):: AUX_GRID_SCALING = 0.2_LCSRP

      !!!!!!! END USER ACCESSIBLE OPTIONS/PARAMS !!!!!!!!


      !
      ! Name of the output and temp directories
      !
      character(len=32),parameter:: OUTPUT_DIR = 'cfd2lcs_output'
      character(len=32),parameter:: TEMP_DIR = 'cfd2lcs_tmp'

      !----
      !Some constants:
      !----
      integer,parameter:: NGHOST_CFD = 1
      integer,parameter:: NMAX_STRUCT = 1000

      !Particle flags:
      integer,parameter:: &
                  LP_IB = -1, &
                  LP_UNKNOWN = -2, &
                  LP_RECYCLE = -3, &
                  LP_STICK = -4

      !Integration directions
      integer,parameter:: FWD = 1,&
                                    BKWD = -1,&
                                    IGNORE = 0

      !i,j,k Cartesian offsets
      integer,parameter:: N_NBR = 27
      integer(LCSIP),parameter:: NBR_OFFSET(3,N_NBR) = reshape((/&
               0, 0,     0, &  !self
              -1,   0,   0, &  !face
               0,  -1,   0, &  !face
               0,   0,  -1, &  !face
               1,   0,   0, &  !face
               0,   1,   0, &  !face
               0,   0,   1, &  !face
              -1,  -1,  -1, &
               0,  -1,  -1, &
               1,  -1,  -1, &
              -1,   0,  -1, &
               1,   0,  -1, &
              -1,   1,  -1, &
               0,   1,  -1, &
               1,   1,  -1, &
              -1,  -1,   0, &
               1,  -1,   0, &
              -1,   1,   0, &
               1,   1,   0, &
              -1,  -1,   1, &
               0,  -1,   1, &
               1,  -1,   1, &
              -1,   0,   1, &
               1,   0,   1, &
              -1,   1,   1, &
               0,   1,   1, &
               1,   1,   1 &
            /),(/3,N_NBR/))


      !----
      !MPI stuff:
      !----
      integer :: nprocs,lcsrank,lcscomm
      integer:: MPI_LCSRP  !Real precision for mpi
      !Structured Comms (scomm_t):
      type scomm_t
            character(len=LCS_NAMELEN):: label
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

      !----
      !Data Structures:
      !----
      !Real, structured rank 0 scalar (sr0_t)
      type sr0_t
            integer:: ni,nj,nk,ng
            character(len=LCS_NAMELEN):: label
            real(LCSRP), allocatable:: r(:,:,:)
      end type sr0_t

      !Real, structured Cartesian rank 1 vector (sr1_t)
      type sr1_t
            integer:: ni,nj,nk,ng
            character(len=LCS_NAMELEN):: label
            real(LCSRP), allocatable:: x(:,:,:)
            real(LCSRP), allocatable:: y(:,:,:)
            real(LCSRP), allocatable:: z(:,:,:)
            logical:: periodic_translate  !For translating x,y,z coordinates over periodic boundaries
      end type sr1_t

      !Real, structured Cartesian rank 2 tensor (sr2_t)
      type sr2_t
            integer:: ni,nj,nk,ng
            character(len=LCS_NAMELEN):: label
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

      !Integer, structured rank 0 scalar (si0_t)
      type si0_t
            integer:: ni,nj,nk,ng
            character(len=LCS_NAMELEN):: label
            integer(LCSIP), allocatable:: i(:,:,:)
      end type si0_t


      !Integer, structured Cartesian rank 1 vector (sr1_t)
      type si1_t
            integer:: ni,nj,nk,ng
            character(len=LCS_NAMELEN):: label
            integer(LCSIP), allocatable:: x(:,:,:)
            integer(LCSIP), allocatable:: y(:,:,:)
            integer(LCSIP), allocatable:: z(:,:,:)
            logical:: periodic_translate  !For translating x,y,z coordinates over periodic boundaries
      end type si1_t


      !Real, unstructured rank 0 scalar (ur0_t)
      type ur0_t
            integer:: n
            character(len=LCS_NAMELEN):: label
            real(LCSRP), allocatable:: r(:)
      end type ur0_t

      !Real, unstructured rank 1 vector (ur1_t)
      type ur1_t
            integer:: n
            character(len=LCS_NAMELEN):: label
            real(LCSRP), allocatable:: x(:)
            real(LCSRP), allocatable:: y(:)
            real(LCSRP), allocatable:: z(:)
      end type ur1_t

      !Real, unstructured rank 2 tensor (ur2_t)
      type ur2_t
            integer:: n
            character(len=LCS_NAMELEN):: label
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
            character(len=LCS_NAMELEN):: label
            integer(LCSIP), allocatable:: i(:)
      end type ui0_t

      !Integer, unstructured rank 1 vector (ui1_t)
      type ui1_t
            integer:: n
            character(len=LCS_NAMELEN):: label
            integer(LCSIP), allocatable:: x(:)
            integer(LCSIP), allocatable:: y(:)
            integer(LCSIP), allocatable:: z(:)
      end type ui1_t

      !Structured grid (sgrid_t):
      type sgrid_t
            character(len=LCS_NAMELEN):: label

            !Dimensions
            integer:: ni, nj, nk, ng
            integer:: gni, gnj, gnk
            integer:: offset_i, offset_j, offset_k

            !Communicators
            type(scomm_t):: scomm_face_r0
            type(scomm_t):: scomm_max_r0
            type(scomm_t):: scomm_face_r1
            type(scomm_t):: scomm_max_r1
            type(scomm_t):: scomm_face_r2
            type(scomm_t):: scomm_max_r2

            !Data
            type(sr1_t):: grid  !Cartesian grid coordinates

            !Least squares gradient wts (for non-rectilinear grids)
            logical:: rectilinear
            integer:: nbr_f,nbr_l
            type(sr1_t),allocatable:: lsgw(:) !Least squares gradient weights

            !Generic boundary conditions:
            type(si0_t):: bcflag !A user indicator for each node
            logical:: periodic_i, periodic_j, periodic_k
            type(sr1_t):: bcnorm !normal for each node

            !periodic_shifts and bounding boxes:
            real(LCSRP),allocatable::  bb(:,:)  !0:nproc-1,1:6
            real(LCSRP),allocatable::  ps(:,:,:,:,:)  !0:nproc-1,1:27,1:3

      end type sgrid_t
      integer:: NSGRID
      type(sgrid_t),target:: sgrid_c(NMAX_STRUCT) !Collection of NSGRID sgrid structures

      !Lagrangian Particles
      type lp_t
            integer(LCSIP):: id  !A unique integer identifier to allow the user to modify aspects of each lcs
            character(len=LCS_NAMELEN):: label
            integer:: np !Number of particles (on each proc)
            integer:: npall !Number of particles (across all procs)
            integer:: direction !Integration direction (FWD or BKWD)
            logical:: recursive_tracking
            real(LCSRP):: lifetime  !time in the domain
            real(LCSRP):: dt_factor !used to compute cfl-like condition when integrating bkwd time flowmap
            type(ur1_t):: xp !Position
            type(ur1_t):: up !Velocity
            type(ur1_t):: dx !Net displacement (needed to handle periodic domains)
            type(ui1_t):: no !Nearest node index (i,j,k) wrt the lcs%sgrid
            type(ui1_t):: no_scfd !Nearest node index (i,j,k) wrt the scfd%sgrid
            type(ui0_t):: proc0 !Origin proc
            type(ui0_t):: no0 !Origin node
            type(ui0_t):: flag !Multipurpose flag
            type(sgrid_t),pointer :: sgrid  !pointer to the structured grid that this lp is associated with
            type(sr1_t):: fm  !Flow Map (On sgrid)
            type(sr1_t):: ugrid !Velocity field (On sgrid)
      end type lp_t
      integer:: NLP
      type(lp_t),target:: lp_c(NMAX_STRUCT)  !Collection of NLP lp structures

      !lcs:
      type lcs_t
            integer(LCSIP):: id  !A unique integer identifier to allow the user to modify aspects of each lcs
            character(len=LCS_NAMELEN):: label
            integer(LCSIP):: diagnostic
            real(LCSRP):: T ! Integration time
            real(LCSRP):: h ! Visualization timestep
            type(sgrid_t),pointer :: sgrid  !pointer to the structured grid
            type(lp_t),pointer:: lp !pointer to Lagrangian particles if used
            type(sr0_t):: ftle  !FTLE field
            !Auxillary grids:
            type(lp_t),pointer:: lpX0, lpX1, lpY0, lpY1, lpZ0, lpZ1
      end type lcs_t
      integer:: NLCS
      type(lcs_t),target:: lcs_c(NMAX_STRUCT)  !Collection of NLCS lcs structures

      !The CFD side data:  There can only be one of these, and we can make it globally available
      type scfd_t
            character(len=LCS_NAMELEN):: label
            real(LCSRP):: t_n, t_np1
            type(sgrid_t),pointer:: sgrid
            type(sr1_t):: u_n !old velocity field
            type(sr1_t):: u_np1  !latest velocity field
            type(sr1_t):: delta  !characteristic length (used for cfl)
      end type scfd_t
      type(scfd_t):: scfd

      !Error handling:
      integer:: CFD2LCS_ERROR = 0  !Initialize to zero here:

      !CPU timing:
      real(LCSRP),save:: t_start_global, t_start_update, t_finish_update
      real(LCSRP),save:: cpu_total_sim,cpu_total_lcs,cpu_fwd,cpu_bkwd,cpu_reconstruct,cpu_io,cpu_lpmap,cpu_ftle,cpu_error
      integer(8),save:: integrations_fwd, integrations_bkwd,integrations_fwd_c,integrations_bkwd_c
      real(LCSRP),save:: this_cpu_fwd, this_cpu_bkwd, cpu_fwd_c,cpu_bkwd_c

      contains
      !These function set the rules for the relationship
      !between 3D (i,j,k) indices and 1D (l) index for structured grids.
      integer function lcs_ijk2l(i,j,k,ni,nj)
            integer::i,j,k !current index
            integer::ni,nj !num pts in i,j directions
            lcs_ijk2l= (k-1)*ni*nj + (j-1)*ni +i
      end function lcs_ijk2l
      integer function l2i(l,ni)
            implicit none
            integer::l !current index
            integer::ni !num pts in i direction
          l2i= mod(l-1,ni)+1
      end function l2i
      integer function l2j(l,ni,nj)
            implicit none
            integer::l !current index
            integer::ni,nj !num pts in i,j directions
            l2j = mod((l-1)/ni,nj)+1
      end function l2j
      integer function l2k(l,ni,nj)
            implicit none
            integer::l !current index
            integer::ni,nj !num pts in i,j directions
            l2k = (l-1)/(ni*nj)+1
      end function l2k

      !A Simple timer function:
      real(LCSRP) function cputimer(mpicomm,sync)
            implicit none
            !-----
            integer:: mpicomm
            logical:: sync
            !-----
            integer:: ierr
            integer:: time_array(8)
            if(sync)then
                  call MPI_BARRIER(mpicomm,ierr)
            endif
            call date_and_time(values=time_array)
            cputimer = time_array (5) * 3600.0_LCSRP + time_array (6) * 60.0_LCSRP &
           + time_array (7) + 0.001_LCSRP * time_array (8)
      end

end module data_m
