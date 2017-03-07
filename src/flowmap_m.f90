!
!Copyright (C) 2015-2016, Justin R. Finn.  All rights reserved.
!libcfd2lcs is distributed is under the terms of the GNU General Public License
!
module flowmap_m
      use data_m
      use io_m
      use structured_m
      use lp_m
      use comms_m
      use lp_tracking_m
      use lp_motion_m
      implicit none

      contains

      subroutine update_flowmap_sl(lp,sgrid,fm,flow)
            use comms_m
            implicit none
            !-----
            type(lp_t):: lp
            type(sgrid_t):: sgrid
            type(sr1_t):: fm
            type(scfd_t):: flow
            !-----
            real(LCSRP):: dt,lp_time,dt_factor
            type(ur1_t):: lpgrid
            type(ur1_t):: fmp
            type(ur0_t):: mask
            integer::n_subcycle, subcycle, ip
            integer:: npall,ierr
            !-----
            !Perform a semi-lagrangian update of the flow map
            !Subcycling is used, if necessary, to ensure that
            !the particle CFL number (wrt the flow grid and LCS grid)
            !is not greater than CFL_MAX
            !-----
            call MPI_REDUCE(lp%np,npall,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)

            if(lcsrank==0 .AND. LCS_VERBOSE) &
                  write(*,*) 'in update_flowmap_sl... ', trim(fm%label), ' NP =',npall

            !-----
            !Initialize some temporary structures for integration:
            !-----
            call init_ur0(mask,lp%np,'MASK')
            call init_ur1(fmp,lp%np,'FM_P')
            call init_ur1(lpgrid,lp%np,'ParticleGrid')
            lpgrid%x(1:lp%np)=lp%xp%x(1:lp%np)
            lpgrid%y(1:lp%np)=lp%xp%y(1:lp%np)
            lpgrid%z(1:lp%np)=lp%xp%z(1:lp%np)

            !-----
            !Set dt, and remember to account for different lcs spacing
            !-----
            call set_dt(dt,n_subcycle,flow,lp)

            !-----
            !Set a mask variable:
            !-----
            do ip = 1,lp%np
                  if(flow%sgrid%bcflag%i(lp%no_scfd%x(ip),lp%no_scfd%y(ip),lp%no_scfd%z(ip)) /= LCS_MASK) then
                        mask%r(ip) = 1.0_LCSRP
                  endif
            enddo

            !-----
            !Advance the particles through the flow timestep
            !-----
            lp_time = flow%t_np1
            do subcycle = 1,n_subcycle
                  if(LCS_VERBOSE) N_UPDATE = 1 
                  if(lcsrank==0 .AND. mod(subcycle, max(N_UPDATE/lp%np,1))==0)&
                        write(*,*) ' Starting SL subcycle',subcycle

                  !integrate
                  call integrate_lp(flow,lp,lp_time,dt)

                  !Interpolate the flow map to the particles at t+dt,
                  call interp_s2u_r1(lp,sgrid,lp%no,fmp,fm) !Interp flow map to lp

                  !Update the flow map (stored in displacement form, relative to the fixed grid).
                  do ip = 1,lp%np
                        fm%x(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip)) = (fmp%x(ip) + (lp%xp%x(ip)-lpgrid%x(ip)))*mask%r(ip)
                        fm%y(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip)) = (fmp%y(ip) + (lp%xp%y(ip)-lpgrid%y(ip)))*mask%r(ip)
                        fm%z(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip)) = (fmp%z(ip) + (lp%xp%z(ip)-lpgrid%z(ip)))*mask%r(ip)
                  enddo

                  !Set Flow map BC:
                  call exchange_sdata(sgrid%scomm_max_r1,r1=fm)
                  call set_flowmap_bc(sgrid,fm)

                  !relocate particles back to original grid
                  lp%xp%x(1:lp%np)=lpgrid%x(1:lp%np)
                  lp%xp%y(1:lp%np)=lpgrid%y(1:lp%np)
                  lp%xp%z(1:lp%np)=lpgrid%z(1:lp%np)
            enddo

            !cleanup
            call destroy_ur1(lpgrid)
            call destroy_ur1(fmp)
            call destroy_ur0(mask)

      end subroutine update_flowmap_sl

      subroutine write_flowmap_substep(lp,T,h)
            implicit none
            !-----
            type(lp_t),pointer:: lp
            real(LCSRP):: T,h
            integer:: gn(3), offset(3)
            character (len=128):: fname,fnameX0,fnameX1,fnameY0,fnameY1,fnameZ0,fnameZ1
            integer:: step
            !-----
            !Write a time h flowmap to disk
            !-----

            if(lcsrank==0 .AND. LCS_VERBOSE) &
                  write(*,*) 'in write_flowmap_substep... ',trim(lp%fm%label)

            step = nint( mod(scfd%t_np1,T)/h)
            gn = (/lp%sgrid%gni,lp%sgrid%gnj,lp%sgrid%gnk/)
            offset = (/lp%sgrid%offset_i,lp%sgrid%offset_j,lp%sgrid%offset_k/)

            if(FILE_EXT == '.h5') then
                  write(fname,'(a,a,a,a,a,i4.4,a)')&
                  './',trim(TEMP_DIR),'/',trim(lp%fm%label),'_',step,FILE_EXT
            endif
            if(FILE_EXT == '.dat') then
                  write(fname,'(a,a,a,i4.4,a,a,a,i4.4,a)')&
                  './',trim(TEMP_DIR),'/',lcsrank,'_',trim(lp%fm%label),'_',step,FILE_EXT
            endif
            call structured_io(trim(fname),IO_WRITE,gn,offset,r1=lp%fm)
            if(lcsrank==0 .AND. LCS_VERBOSE) &
                  write(*,*) ' Time h flowmap saved to: ', trim(fname)

      end subroutine write_flowmap_substep

      subroutine reconstruct_flowmap(lp,T,h,success)
            use sgrid_m
            implicit none
            type(lp_t),pointer:: lp
            real(LCSRP):: T,h
            logical:: success
            !-----
            integer:: i, step, nstep, inc
            integer:: gn(3), offset(3)
            character(len=128):: fname
            type(sr1_t):: fm_tmp
            logical:: file_exists
            integer:: ni,nj,nk
            integer:: ip,np
            real:: t1,t0
            !-----
            !Reconstruct the time T flow map from N time h substeps
            !-----

            if(lcsrank==0) &
                  write(*,'(a,a)') 'Reconstructing flowmap: ',trim(lp%fm%label)

            !-----
            !Initialize temporary data:
            !-----
            gn = (/lp%sgrid%gni,lp%sgrid%gnj,lp%sgrid%gnk/)
            offset = (/lp%sgrid%offset_i,lp%sgrid%offset_j,lp%sgrid%offset_k/)
            call init_sr1(fm_tmp,lp%sgrid%ni,lp%sgrid%nj,lp%sgrid%nk,lp%sgrid%ng,'fm_tmp',.FALSE.)

            !-----
            !Reset the lp grid, which will be used to gather the time T flow map
            !-----
            call reset_lp(lp)
            ni = lp%sgrid%ni; nj = lp%sgrid%nj; nk = lp%sgrid%nk;
            do ip = 1,lp%np
                  lp%no%x(ip) =  l2i(lp%no0%i(ip),ni)
                  lp%no%y(ip) =  l2j(lp%no0%i(ip),ni,nj)
                  lp%no%z(ip) =  l2k(lp%no0%i(ip),ni,nj)
            enddo

            !-----
            !Loop through each of the substeps and construct the time T flowmap:
            !-----
            step = nint( mod(scfd%t_np1,T)/h)
            nstep = nint(T/h)
            if (lp%direction == FWD) then
                  inc = 1
                  step = step + inc
            else
                  inc = -1
            endif
            !Set the lp direction to FWD, and correct after reconstruct
            lp%direction = FWD

            success = .true.
            do i = 1,nstep
                  !Load this substep (or skip, if there is no data yet:
                  if (step >= nstep) step = 0
                  if (step < 0) step = nstep -1
                  if(FILE_EXT == '.h5') then
                        write(fname,'(a,a,a,a,a,i4.4,a)')'./',trim(TEMP_DIR),'/',trim(lp%fm%label),'_',step,FILE_EXT
                  endif
                  if(FILE_EXT == '.dat') then
                        write(fname,'(a,a,a,i4.4,a,a,a,i4.4,a)')'./'&
                        ,trim(TEMP_DIR),'/',lcsrank,'_',trim(lp%fm%label),'_',step,FILE_EXT
                  endif
                  INQUIRE(FILE=trim(fname), EXIST=file_exists)
                  if (.NOT. file_exists) then
                        !next step
                        step = step + inc
                        success = .false.
                        cycle
                  endif
                  !if(lcsrank==0 .AND. LCS_VERBOSE) &
                  if(lcsrank==0 .AND. LCS_VERBOSE) &
                        write(*,*) i,'about to open file: ', trim(fname)
                  call structured_io(trim(fname),IO_READ,gn,offset,r1=lp%fm)  !read the flowmap

                  t0 = cputimer(lcscomm,SYNC_TIMER)

                  !Set Flow map BC:
                  call exchange_sdata(lp%sgrid%scomm_max_r1,r1=lp%fm)
                  call set_flowmap_bc(lp%sgrid,lp%fm)

                  !Project into the plane (if 2d) and handle boundary conditions:
                  !Note, we always use the grid associated with the lp, not the flow here.
                  !This is to avoid issues with Losing particles during reconstruct
                  call project_lp2plane(lp,lp%sgrid,lp%no)
                  call set_lp_bc(lp,lp%sgrid,lp%no)
                  
                  !Interp flow map to particles (store in lp%up)
                  call interp_s2u_r1(lp,lp%sgrid,lp%no,lp%up,lp%fm)

                  !Update particle positions and flow map
                  np = lp%np
                  lp%xp%x(1:np) = lp%xp%x(1:np) + lp%up%x(1:np)
                  lp%xp%y(1:np) = lp%xp%y(1:np) + lp%up%y(1:np)
                  lp%xp%z(1:np) = lp%xp%z(1:np) + lp%up%z(1:np)
                  lp%dx%x(1:np) = lp%dx%x(1:np) + lp%up%x(1:np)
                  lp%dx%y(1:np) = lp%dx%y(1:np) + lp%up%y(1:np)
                  lp%dx%z(1:np) = lp%dx%z(1:np) + lp%up%z(1:np)

                  !Exchange particles
                  call exchange_lp_alltoall(lp,lp%sgrid)

                  !Track to the grid:  Note, we want the LP grid, not the scfd grid
                  !in order to evaluate the flow map.
                  call track_lp2node(lp,lp%sgrid,lp%no)

                  !next step
                  step = step + inc

                  t1 = cputimer(lcscomm,SYNC_TIMER)
                  cpu_reconstruct = cpu_reconstruct + max(t1-t0,0.0)
            enddo

            !map back to the origin:
            call exchange_lpmap(lp)

            !Set Flow map BC:
            call exchange_sdata(lp%sgrid%scomm_max_r1,r1=lp%fm)
            call set_flowmap_bc(lp%sgrid,lp%fm)

            !cleanup
            call destroy_sr1(fm_tmp)
            if (inc == -1) lp%direction = BKWD  !reset direction

      end subroutine reconstruct_flowmap

      subroutine set_flowmap_bc(sgrid,fm)
            implicit none
            !-----
            type(sr1_t):: fm
            type(sgrid_t):: sgrid
            !-----
            integer:: i,j,k,ni,nj,nk,ng
            integer:: i_b,j_b,k_b
            integer:: i_ib,j_ib,k_ib
            !-----
            !Set the flowmap boundary conditions.
            !-----
            if(lcsrank==0 .and. LCS_VERBOSE) &
                  write(*,*) 'In set_flowmap_bc... ', trim(fm%label)

            ni = fm%ni
            nj = fm%nj
            nk = fm%nk
            ng = fm%ng

            !First pass: handle any mask conditions  or inflow outflow, where we
            !want the flowmap to be zero in the IB Nodes:
            do k = 1,nk
            do j = 1,nj
            do i = 1,ni
                  select case(sgrid%bcflag%i(i,j,k))
                  case(LCS_MASK,LCS_INFLOW,LCS_OUTFLOW)
                        fm%x(i,j,k) = 0.0_LCSRP
                        fm%y(i,j,k) = 0.0_LCSRP
                        fm%z(i,j,k) = 0.0_LCSRP
                  case default
                        !Do nothing...
                  end select
            enddo
            enddo
            enddo

            !Second pass:  set the ghost/fake values according to the desired condition
            do k = 1-ng,nk+ng
            do j = 1-ng,nj+ng
            do i = 1-ng,ni+ng
                  if(i>=1 .and. j>=1 .and. k>=1 .and. i<=ni .and. j<=nj .and. k<=nk) cycle
                  select case(sgrid%bcflag%i(i,j,k))
                  case(LCS_INTERNAL)
                        cycle
                  case(LCS_WALL,LCS_MASK,LCS_INFLOW,LCS_OUTFLOW)
                        !Zero the flow map in any ghost/fake that is a LCS_WALL
                        !This corresponds with the Lagrangian Stick condition.
                        fm%x(i,j,k) = 0.0_LCSRP
                        fm%y(i,j,k) = 0.0_LCSRP
                        fm%z(i,j,k) = 0.0_LCSRP
                  case(LCS_SLIP,LCS_2D)
                        !Set 0 gradient WRT to the In-bounds direction at ghost/fake nodes:
                        i_b = max(min(i,ni),1)
                        j_b = max(min(j,nj),1)
                        k_b = max(min(k,nk),1)
                        fm%x(i,j,k) = fm%x(i_b,j_b,k_b)
                        fm%y(i,j,k) = fm%y(i_b,j_b,k_b)
                        fm%z(i,j,k) = fm%z(i_b,j_b,k_b)
                  case default
                        write(*,*) 'lcsrank[',lcsrank,'] ERROR: Unknown bcflag:',sgrid%bcflag%i(i,j,k)
                        CFD2LCS_ERROR=1
                  end select
            enddo
            enddo
            enddo

      end subroutine set_flowmap_bc

end module flowmap_m
