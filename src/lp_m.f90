!
!Copyright (C) 2015-2016, Justin R. Finn.  All rights reserved.
!libcfd2lcs is distributed is under the terms of the GNU General Public License
!
module lp_m
      use data_m
      use unstructured_m
      use structured_m
      implicit none
      !Basic routines for the lp structure.
      contains

      subroutine init_lp(lp,sgrid,direction,dt_factor,label)
            !----
            type(lp_t),pointer:: lp
            character(len=*):: label
            type(sgrid_t),pointer:: sgrid
            integer:: direction
            real(LCSRP):: dt_factor
            integer::ierr
            !----
            integer:: i,j,k,ip
            !----
            !Initialize a new set of lagrangian particles with certain properties
            !Initial positions will be determined by the sgrid passed in.
            !----

            if (lcsrank==0)&
                  write(*,'(a,a)') 'in init_lp...',trim(label)

            !----
            !Point to the next available LP
            !----
            NLP = NLP + 1
            if(NLP > NMAX_STRUCT) then
                  if(lcsrank==0)then
                        write(*,*) 'ERROR: NLP exceeds NMAX_STRUCT.',NLP,NMAX_STRUCT
                        write(*,*) 'Increase NMAX_STRUCT in data_m.f90 and recompile libcfd2lcs'
                  endif
                  CFD2LCS_ERROR = 1
                  return
            endif
            lp => lp_c(NLP)
            lp%label = trim(label)
            lp%np = sgrid%ni*sgrid%nj*sgrid%nk
            lp%sgrid => sgrid
            lp%dt_factor = dt_factor
            call MPI_ALLREDUCE(lp%np,lp%npall,1,MPI_INTEGER,MPI_SUM,lcscomm,ierr)
            lp%id= NLP
            lp%lifetime = 0.0_LCSRP

            !initialize
            call init_ur1(lp%xp,lp%np,'XP')
            call init_ur1(lp%up,lp%np,'UP')
            call init_ur1(lp%dx,lp%np,'DX')
            call init_ui1(lp%no,lp%np,'NODE')
            call init_ui1(lp%no_scfd,lp%np,'SCFD_NODE')
            call init_ui0(lp%no0,lp%np,'NODE0')
            call init_ui0(lp%proc0,lp%np,'PROC0')
            call init_ui0(lp%flag,lp%np,'FLAG')

            select case(direction)
                  case(FWD)
                        lp%direction = FWD
                  case(BKWD)
                        lp%direction = BKWD
                  case default
                        write(*,*) 'ERROR:  BAD DIRECTION:', direction
                        CFD2LCS_ERROR = 1
            end select

            ip = 0
            do k = 1,sgrid%nk
            do j = 1,sgrid%nj
            do i = 1,sgrid%ni
                  ip = ip + 1
                  lp%xp%x(ip) = sgrid%grid%x(i,j,k)
                  lp%xp%y(ip) = sgrid%grid%y(i,j,k)
                  lp%xp%z(ip) = sgrid%grid%z(i,j,k)
                  lp%proc0%i(ip) = lcsrank
                  lp%no0%i(ip) = lcs_ijk2l(i,j,k,sgrid%grid%ni,sgrid%grid%nj) !Cartesian ordering
                  lp%no%x(ip) =  i
                  lp%no%y(ip) =  j
                  lp%no%z(ip) =  k
                  !particles start at unknown scfd node, but you can make an educated guess:
                  lp%no_scfd%x(ip) = min(max(nint(real(i*scfd%sgrid%ni)/real(sgrid%grid%ni)),1),scfd%sgrid%ni) !Guess
                  lp%no_scfd%y(ip) = min(max(nint(real(j*scfd%sgrid%nj)/real(sgrid%grid%nj)),1),scfd%sgrid%nj) !Guess
                  lp%no_scfd%z(ip) = min(max(nint(real(k*scfd%sgrid%nk)/real(sgrid%grid%nk)),1),scfd%sgrid%nk) !Guess
                  lp%flag%i(ip) = LP_UNKNOWN
            enddo
            enddo
            enddo

            !Make sure we do a recursive search
            lp%recursive_tracking = .true.

            !Initalize the flowmap on the structured grid:
            call init_sr1(lp%fm,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,trim(lp%label)//'-FM',translate=.false.)
            lp%fm%x = 0.0_LCSRP
            lp%fm%y = 0.0_LCSRP
            lp%fm%z = 0.0_LCSRP
            
            if (VELOCITY_IO) then
                  call init_sr1(lp%ugrid,sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,trim(lp%label)//'-UGRID',translate=.false.)
                  lp%ugrid%x = 0.0_LCSRP
                  lp%ugrid%y = 0.0_LCSRP
                  lp%ugrid%z = 0.0_LCSRP
            endif

            if(lcsrank==0)write(*,'(a,a,a,ES18.4,a)') &
                  'LP: ',trim(lp%label),' has ',real(lp%npall,LCSRP), ' particles across all procs'

      end subroutine init_lp

      subroutine init_lp_tracer(lp,sgrid,direction,dt_factor,label)
            !----
            type(lp_t),pointer:: lp
            character(len=*):: label
            type(sgrid_t),pointer:: sgrid
            integer:: direction
            real(LCSRP):: dt_factor
            !----
            integer:: i,j,k,ip
            real(LCSRP):: r,dx,dy,dz
            integer::ierr
            !----
            !Initialize a new set of lagrangian particles with certain properties
            !Initial positions will be determined by the sgrid passed in and the 
            !Tracer injection properties
            !----

            if (lcsrank==0)&
                  write(*,'(a,a)') 'in init_lp_tracer...',trim(label)

            !----
            !Point to the next available LP
            !----
            NLP = NLP + 1
            if(NLP > NMAX_STRUCT) then
                  if(lcsrank==0)then
                        write(*,*) 'ERROR: NLP exceeds NMAX_STRUCT.',NLP,NMAX_STRUCT
                        write(*,*) 'Increase NMAX_STRUCT in data_m.f90 and recompile libcfd2lcs'
                  endif
                  CFD2LCS_ERROR = 1
                  return
            endif
            lp => lp_c(NLP)
            lp%label = trim(label)
            lp%sgrid => sgrid
            lp%dt_factor = dt_factor
            lp%id= NLP
            
            !Figure out np
            lp%np = 0
            do k = 1,sgrid%nk
            do j = 1,sgrid%nj
            do i = 1,sgrid%ni
                  dx = sgrid%grid%x(i,j,k) - TRACER_INJECT_X
                  dy = sgrid%grid%y(i,j,k) - TRACER_INJECT_Y
                  dz = sgrid%grid%z(i,j,k) - TRACER_INJECT_Z
                  r = sqrt(dx*dx+dy*dy+dz*dz)
                  if (r <= TRACER_INJECT_RADIUS) lp%np = lp%np+1
            enddo
            enddo
            enddo
            call MPI_ALLREDUCE(lp%np,lp%npall,1,MPI_INTEGER,MPI_SUM,lcscomm,ierr)

            !initialize
            call init_ur1(lp%xp,lp%np,'XP')
            call init_ur1(lp%up,lp%np,'UP')
            call init_ur1(lp%dx,lp%np,'DX')
            call init_ui1(lp%no,lp%np,'NODE')
            call init_ui1(lp%no_scfd,lp%np,'SCFD_NODE')
            call init_ui0(lp%no0,lp%np,'NODE0')
            call init_ui0(lp%proc0,lp%np,'PROC0')
            call init_ui0(lp%flag,lp%np,'FLAG')

            select case(direction)
                  case(FWD)
                        lp%direction = FWD
                  case(BKWD)
                        lp%direction = BKWD
                  case default
                        write(*,*) 'ERROR:  BAD DIRECTION:', direction
                        CFD2LCS_ERROR = 1
            end select

            ip = 0
            do k = 1,sgrid%nk
            do j = 1,sgrid%nj
            do i = 1,sgrid%ni
                  dx = sgrid%grid%x(i,j,k) - TRACER_INJECT_X
                  dy = sgrid%grid%y(i,j,k) - TRACER_INJECT_Y
                  dz = sgrid%grid%z(i,j,k) - TRACER_INJECT_Z
                  r = sqrt(dx*dx+dy*dy+dz*dz)
                  if (r > TRACER_INJECT_RADIUS .or. ip == lp%np) cycle

                  ip = ip + 1
                  lp%xp%x(ip) = sgrid%grid%x(i,j,k)
                  lp%xp%y(ip) = sgrid%grid%y(i,j,k)
                  lp%xp%z(ip) = sgrid%grid%z(i,j,k)
                  lp%proc0%i(ip) = lcsrank
                  lp%no0%i(ip) = lcs_ijk2l(i,j,k,sgrid%grid%ni,sgrid%grid%nj) !Cartesian ordering
                  lp%no%x(ip) =  i
                  lp%no%y(ip) =  j
                  lp%no%z(ip) =  k
                  !particles start at unknown scfd node, but you can make an educated guess:
                  lp%no_scfd%x(ip) = min(max(nint(real(i*scfd%sgrid%ni)/real(sgrid%grid%ni)),1),scfd%sgrid%ni) !Guess
                  lp%no_scfd%y(ip) = min(max(nint(real(j*scfd%sgrid%nj)/real(sgrid%grid%nj)),1),scfd%sgrid%nj) !Guess
                  lp%no_scfd%z(ip) = min(max(nint(real(k*scfd%sgrid%nk)/real(sgrid%grid%nk)),1),scfd%sgrid%nk) !Guess
                  lp%flag%i(ip) = LP_UNKNOWN
            enddo
            enddo
            enddo

            !Make sure we do a recursive search
            lp%recursive_tracking = .true.

            !Note, we dont initialize the lp%fm or lp%ugrid here.
            
            
            if(lcsrank==0)write(*,'(a,a,a,ES18.4,a)') &
                  'LP: ',trim(lp%label),' has ',real(lp%npall,LCSRP), ' particles across all procs'

      end subroutine init_lp_tracer

      subroutine reset_lp(lp)
            !----
            type(lp_t),pointer:: lp
            !----
            type(sr1_t),pointer:: grid
            integer:: i,j,k,ip
            !----
            !Reset Lagrangian particles to a grid
            !----

            if(.NOT. associated(lp)) then
                  write(*,'(a)') 'ERROR:  trying to reset lp that is not associated...'
                  CFD2LCS_ERROR = 1
                  return
            else
                  if (lcsrank==0 .AND. LCS_VERBOSE)&
                        write(*,*) 'in reset_lp...',trim(lp%label)
            endif

            !Reset to the grid
            grid => lp%sgrid%grid

            !Point to the new lp and initialize
            lp%np = grid%ni*grid%nj*grid%nk
            call init_ur1(lp%xp,lp%np,'XP')
            call init_ur1(lp%up,lp%np,'UP')
            call init_ur1(lp%dx,lp%np,'DX')
            call init_ui1(lp%no,lp%np,'NODE')
            call init_ui1(lp%no_scfd,lp%np,'NODE_SCFD')
            call init_ui0(lp%no0,lp%np,'NODE0')
            call init_ui0(lp%proc0,lp%np,'PROC0')
            call init_ui0(lp%flag,lp%np,'FLAG')

            ip = 0
            do k = 1,grid%nk
            do j = 1,grid%nj
            do i = 1,grid%ni
                  ip = ip + 1
                  lp%xp%x(ip) = grid%x(i,j,k)
                  lp%xp%y(ip) = grid%y(i,j,k)
                  lp%xp%z(ip) = grid%z(i,j,k)
                  lp%proc0%i(ip) = lcsrank
                  lp%no0%i(ip) = lcs_ijk2l(i,j,k,grid%ni,grid%nj) !Cartesian ordering
                  lp%no%x(ip) =i
                  lp%no%y(ip) =j
                  lp%no%z(ip) =k
                  !particles start at unknown scfd node, but you can make an educated guess:
                  lp%no_scfd%x(ip) = min(max(nint(real(i*scfd%sgrid%ni)/real(grid%ni)),1),scfd%sgrid%ni) !Guess
                  lp%no_scfd%y(ip) = min(max(nint(real(j*scfd%sgrid%nj)/real(grid%nj)),1),scfd%sgrid%nj) !Guess
                  lp%no_scfd%z(ip) = min(max(nint(real(k*scfd%sgrid%nk)/real(grid%nk)),1),scfd%sgrid%nk) !Guess
                  lp%flag%i(ip) = LP_UNKNOWN
            enddo
            enddo
            enddo

            !Make sure we do a recursive search
            lp%recursive_tracking = .true.

            !Reset the flowmap and velocity on the structured grid:
            if(allocated(lp%fm%x)) lp%fm%x = 0.0_LCSRP
            if(allocated(lp%fm%y)) lp%fm%y = 0.0_LCSRP
            if(allocated(lp%fm%z)) lp%fm%z = 0.0_LCSRP
            if(allocated(lp%ugrid%x)) lp%ugrid%x = 0.0_LCSRP
            if(allocated(lp%ugrid%y)) lp%ugrid%y = 0.0_LCSRP
            if(allocated(lp%ugrid%z)) lp%ugrid%z = 0.0_LCSRP

      end subroutine reset_lp


      subroutine resize_lp(lp,np)
            implicit none
            !-----
            type(lp_t):: lp
            integer:: np
            !-----
            !Dynamically resize lp
            !-----

            if(lcsrank==0 .AND. LCS_VERBOSE)&
                  write(*,*) 'in resize_lp...',trim(lp%label)

            !Check the size of the arrays and resize if needed:
            call resize_ur1(lp%xp,np)
            call resize_ur1(lp%up,np)
            call resize_ur1(lp%dx,np)
            call resize_ui1(lp%no,np)
            call resize_ui1(lp%no_scfd,np)
            call resize_ui0(lp%proc0,np)
            call resize_ui0(lp%no0,np)
            call resize_ui0(lp%flag,np)

            !Set the new np
            lp%np = np

      end subroutine resize_lp

      subroutine reorder_lp(lp)
            implicit none
            !-------
            type(lp_t):: lp
            !-------
            integer:: ip, new_np
            integer:: recycle_count
            !-------
            !Clean up any holes in the lp structure
            !This needs to be called after processors
            !exchange lp's, or if LP's leave the domain,
            !and before unpacking any new ones.
            !It is up to the user to resize afterwards.
            !-------

            if (lcsrank==0 .AND. LCS_VERBOSE) &
                  write(*,*) 'in reorder_lp...',trim(lp%label)

            new_np = 0
            recycle_count = 0
            do ip = 1,lp%np
                  if (lp%flag%i(ip) == LP_RECYCLE) then
                        recycle_count = recycle_count+1
                  else
                        new_np = new_np+1
                        if (new_np/=ip) then
                              lp%xp%x(new_np)         = lp%xp%x(ip)
                              lp%xp%y(new_np)         = lp%xp%y(ip)
                              lp%xp%z(new_np)         = lp%xp%z(ip)
                              lp%up%x(new_np)         = lp%up%x(ip)
                              lp%up%y(new_np)         = lp%up%y(ip)
                              lp%up%z(new_np)         = lp%up%z(ip)
                              lp%dx%x(new_np)         = lp%dx%x(ip)
                              lp%dx%y(new_np)         = lp%dx%y(ip)
                              lp%dx%z(new_np)         = lp%dx%z(ip)
                              lp%no%x(new_np)         = lp%no%x(ip)
                              lp%no%y(new_np)         = lp%no%y(ip)
                              lp%no%z(new_np)         = lp%no%z(ip)
                              lp%no_scfd%x(new_np)    = lp%no_scfd%x(ip)
                              lp%no_scfd%y(new_np)    = lp%no_scfd%y(ip)
                              lp%no_scfd%z(new_np)    = lp%no_scfd%z(ip)
                              lp%proc0%i(new_np)            = lp%proc0%i(ip)
                              lp%no0%i(new_np)        = lp%no0%i(ip)
                              lp%flag%i(new_np)             = lp%flag%i(ip)
                        endif
                  endif
            enddo

            !Note, no resize done here, but need to set n for all members of lp.
            lp%np = new_np
            lp%xp%n = new_np
            lp%up%n = new_np
            lp%dx%n = new_np
            lp%no%n = new_np
            lp%no_scfd%n = new_np
            lp%proc0%n = new_np
            lp%no0%n = new_np
            lp%flag%n = new_np

      end subroutine reorder_lp

      subroutine destroy_lp(lp)
            implicit none
            !-------
            type(lp_t),pointer:: lp
            !-------
            !Remove the lp from the lp collection (lp_c)
            !-------

            if (lcsrank==0)&
                  write(*,*) 'in destroy_lp...',trim(lp%label)

            lp%direction = IGNORE
            lp%label = 'LP_NOT_USED'
            lp%np = 0
            lp%recursive_tracking = .true.
            lp%dt_factor = 1.0_LCSRP

            call destroy_ur1(lp%xp)
            call destroy_ur1(lp%up)
            call destroy_ur1(lp%dx)
            call destroy_ui1(lp%no)
            call destroy_ui1(lp%no_scfd)
            call destroy_ui0(lp%no0)
            call destroy_ui0(lp%proc0)
            call destroy_ui0(lp%flag)
            call destroy_sr1(lp%fm)
            call destroy_sr1(lp%ugrid)

            if(associated(lp%sgrid)) nullify(lp%sgrid)
            if(associated(lp)) nullify(lp)

      end subroutine destroy_lp

      subroutine project_lp2plane(lp,sgrid,no)
            implicit none
            type(lp_t):: lp
            type(sgrid_t):: sgrid
            type(ui1_t):: no
            !----
            type(sr1_t):: v1,v2,norm
            type(ur1_t):: d2p,normp
            type(ur0_t):: nmag
            integer::ni,nj,nk,ng,np,ip
            integer::i,j,k
            !----
            !This routine projects the particles into the (local) plane of data.
            !If gni,gnj,gnk>1, we dont need to do this:
            !----
            if(sgrid%gni >1 .and. sgrid%gnj > 1 .and. sgrid%gnk > 1) return

            if(lcsrank==0 .and. LCS_VERBOSE) &
                  write(*,*) 'in project_lp2plane...'

            ni = sgrid%ni
            nj = sgrid%nj
            nk = sgrid%nk
            ng = sgrid%ng
            np = lp%np

            call init_sr1(v1,ni,nj,nk,ng,'v1',translate=.false.)
            call init_sr1(v2,ni,nj,nk,ng,'v2',translate=.false.)
            call init_sr1(norm,ni,nj,nk,ng,'norm',translate=.false.)
            call init_ur1(d2p,lp%np,'D2P')
            call init_ur1(normp,lp%np,'NORMP')
            call init_ur0(nmag,lp%np,'nmag')

            !Compute the normal to the plane
            !by cross product of two local vectors
            !Also, ensure that the lp node index is correct:
            if(sgrid%gni == 1) then
                  v1%x(1:ni,1:nj,1:nk) = sgrid%grid%x(1:ni,2:nj+1,1:nk) - sgrid%grid%x(1:ni,0:nj-1,1:nk)
                  v1%y(1:ni,1:nj,1:nk) = sgrid%grid%y(1:ni,2:nj+1,1:nk) - sgrid%grid%y(1:ni,0:nj-1,1:nk)
                  v1%z(1:ni,1:nj,1:nk) = sgrid%grid%z(1:ni,2:nj+1,1:nk) - sgrid%grid%z(1:ni,0:nj-1,1:nk)
                  v2%x(1:ni,1:nj,1:nk) = sgrid%grid%x(1:ni,1:nj,2:nk+1) - sgrid%grid%x(1:ni,1:nj,0:nk-1)
                  v2%y(1:ni,1:nj,1:nk) = sgrid%grid%y(1:ni,1:nj,2:nk+1) - sgrid%grid%y(1:ni,1:nj,0:nk-1)
                  v2%z(1:ni,1:nj,1:nk) = sgrid%grid%z(1:ni,1:nj,2:nk+1) - sgrid%grid%z(1:ni,1:nj,0:nk-1)
                  lp%no%x(1:lp%np) = 1
                  lp%no_scfd%x(1:lp%np) = 1
            endif
            if(sgrid%gnj == 1) then
                  v1%x(1:ni,1:nj,1:nk) = sgrid%grid%x(2:ni+1,1:nj,1:nk) - sgrid%grid%x(0:ni-1,1:nj,1:nk)
                  v1%y(1:ni,1:nj,1:nk) = sgrid%grid%y(2:ni+1,1:nj,1:nk) - sgrid%grid%y(0:ni-1,1:nj,1:nk)
                  v1%z(1:ni,1:nj,1:nk) = sgrid%grid%z(2:ni+1,1:nj,1:nk) - sgrid%grid%z(0:ni-1,1:nj,1:nk)
                  v2%x(1:ni,1:nj,1:nk) = sgrid%grid%x(1:ni,1:nj,2:nk+1) - sgrid%grid%x(1:ni,1:nj,0:nk-1)
                  v2%y(1:ni,1:nj,1:nk) = sgrid%grid%y(1:ni,1:nj,2:nk+1) - sgrid%grid%y(1:ni,1:nj,0:nk-1)
                  v2%z(1:ni,1:nj,1:nk) = sgrid%grid%z(1:ni,1:nj,2:nk+1) - sgrid%grid%z(1:ni,1:nj,0:nk-1)
                  lp%no%y(1:lp%np) = 1
                  lp%no_scfd%y(1:lp%np) = 1
            endif
            if(sgrid%gnk == 1) then
                  v1%x(1:ni,1:nj,1:nk) = sgrid%grid%x(2:ni+1,1:nj,1:nk) - sgrid%grid%x(0:ni-1,1:nj,1:nk)
                  v1%y(1:ni,1:nj,1:nk) = sgrid%grid%y(2:ni+1,1:nj,1:nk) - sgrid%grid%y(0:ni-1,1:nj,1:nk)
                  v1%z(1:ni,1:nj,1:nk) = sgrid%grid%z(2:ni+1,1:nj,1:nk) - sgrid%grid%z(0:ni-1,1:nj,1:nk)
                  v2%x(1:ni,1:nj,1:nk) = sgrid%grid%x(1:ni,2:nj+1,1:nk) - sgrid%grid%x(1:ni,0:nj-1,1:nk)
                  v2%y(1:ni,1:nj,1:nk) = sgrid%grid%y(1:ni,2:nj+1,1:nk) - sgrid%grid%y(1:ni,0:nj-1,1:nk)
                  v2%z(1:ni,1:nj,1:nk) = sgrid%grid%z(1:ni,2:nj+1,1:nk) - sgrid%grid%z(1:ni,0:nj-1,1:nk)
                  lp%no%z(1:lp%np) = 1
                  lp%no_scfd%z(1:lp%np) = 1
            endif
            norm%x = v1%y*v2%z - v1%z*v2%y
            norm%y = v1%z*v2%x - v1%x*v2%z
            norm%z = v1%x*v2%y - v1%y*v2%x

            !get the unit normal and distance to each plane for every particle
            do ip = 1,np
                  !find the nearest in-bounds node
                  i = max(min(no%x(ip),ni),1)
                  j = max(min(no%y(ip),nj),1)
                  k = max(min(no%z(ip),nk),1)
                  normp%x(ip) = norm%x(i,j,k)
                  normp%y(ip) = norm%y(i,j,k)
                  normp%z(ip) = norm%z(i,j,k)
                  d2p%x(ip) = sgrid%grid%x(i,j,k)
                  d2p%y(ip) = sgrid%grid%y(i,j,k)
                  d2p%z(ip) = sgrid%grid%z(i,j,k)
            enddo

            !Make normp a unit vector
            nmag%r = sqrt(normp%x*normp%x+normp%y*normp%y+normp%z*normp%z)
            normp%x = normp%x/nmag%r
            normp%y = normp%y/nmag%r
            normp%z = normp%z/nmag%r

            !distance to particle
            d2p%x(1:np) = lp%xp%x(1:np) - d2p%x(1:np)
            d2p%y(1:np) = lp%xp%y(1:np) - d2p%y(1:np)
            d2p%z(1:np) = lp%xp%z(1:np) - d2p%z(1:np)

            !store d2p dot normp
            nmag%r = d2p%x*normp%x + d2p%y*normp%y + d2p%z*normp%z

            !Ensure sign of normp and d2p are alligned
            do ip =1,np
                  if(nmag%r(ip) < 0.0_LCSRP) then
                        nmag%r(ip) = -nmag%r(ip)
                        normp%x(ip) = -normp%x(ip)
                        normp%y(ip) = -normp%y(ip)
                        normp%z(ip) = -normp%z(ip)
                  endif
            enddo

            !correct position
            lp%xp%x(1:np) = lp%xp%x(1:np) - nmag%r(1:np) * normp%x(1:np)
            lp%xp%y(1:np) = lp%xp%y(1:np) - nmag%r(1:np) * normp%y(1:np)
            lp%xp%z(1:np) = lp%xp%z(1:np) - nmag%r(1:np) * normp%z(1:np)

            call destroy_sr1(v1)
            call destroy_sr1(v2)
            call destroy_sr1(norm)
            call destroy_ur1(d2p)
            call destroy_ur1(normp)
            call destroy_ur0(nmag)

      end subroutine project_lp2plane

end module lp_m
