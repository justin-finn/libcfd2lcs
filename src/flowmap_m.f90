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

	subroutine update_flowmap_sl(lcs,flow)
		use comms_m
		implicit none
		!-----
		type(lcs_t),pointer:: lcs
		type(scfd_t):: flow
		!-----
		type(lp_t),pointer:: lp
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
		call MPI_REDUCE(lcs%lp%np,npall,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)

		if(lcsrank==0) &
			write(*,*) 'in update_flowmap_sl... ', trim(lcs%fm%label), ' NP =',npall

		!-----
		!Initialize some temporary structures for integration:
		!-----
		lp => lcs%lp
		call init_ur0(mask,lp%np,'MASK')
		call init_ur1(fmp,lp%np,'FM_P')
		call init_ur1(lpgrid,lp%np,'ParticleGrid')
		lpgrid%x(1:lp%np)=lp%xp%x(1:lp%np)
		lpgrid%y(1:lp%np)=lp%xp%y(1:lp%np)
		lpgrid%z(1:lp%np)=lp%xp%z(1:lp%np)

		!-----
		!Set dt, and remember to account for different lcs spacing
		!-----
		dt_factor = 1.0_LCSRP/real(max(1+lcs%resolution, 1))
		call set_dt(dt,n_subcycle,flow,lp,dt_factor)

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
			if(lcsrank==0 .AND. mod(subcycle, max(N_UPDATE/lp%np,1))==0)&
				write(*,*) ' Starting SL subcycle',subcycle

			!integrate
			call integrate_lp(flow,lp,lp_time,dt)

			!Interpolate the flow map to the particles at t+dt,
			call interp_s2u_r1(lp,lcs%sgrid,lp%no,fmp,lcs%fm) !Interp flow map to lp

			!Update the flow map (stored in displacement form, relative to the fixed grid).
			do ip = 1,lp%np
				lcs%fm%x(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip)) = (fmp%x(ip) + (lp%xp%x(ip)-lpgrid%x(ip)))*mask%r(ip)
				lcs%fm%y(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip)) = (fmp%y(ip) + (lp%xp%y(ip)-lpgrid%y(ip)))*mask%r(ip)
				lcs%fm%z(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip)) = (fmp%z(ip) + (lp%xp%z(ip)-lpgrid%z(ip)))*mask%r(ip)
			enddo

			!Set Flow map BC:
			call exchange_sdata(lcs%sgrid%scomm_max_r1,r1=lcs%fm)
			call set_flowmap_bc(lcs%sgrid,lcs%fm)

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

	subroutine write_flowmap_substep(lcs)
		implicit none
		!-----
		type(lcs_t):: lcs
		integer:: gn(3), offset(3)
		character (len=128):: fname
		integer:: step
		!-----
		!Write a time h flowmap to disk
		!-----

		if(lcsrank==0 .AND. LCS_VERBOSE) &
			write(*,*) 'in write_flowmap_substep... ',trim(lcs%fm%label)

		step = nint( mod(scfd%t_np1,lcs%T)/lcs%h)
		gn = (/lcs%sgrid%gni,lcs%sgrid%gnj,lcs%sgrid%gnk/)
		offset = (/lcs%sgrid%offset_i,lcs%sgrid%offset_j,lcs%sgrid%offset_k/)
		if(FILE_EXT == '.h5') then
			write(fname,'(a,a,a,a,a,i4.4,a)'),'./',trim(TEMP_DIR),'/',trim(lcs%label),'_',step,FILE_EXT
		endif
		if(FILE_EXT == '.dat') then
			write(fname,'(a,a,a,i4.4,a,a,a,i4.4,a)'),'./',trim(TEMP_DIR),'/',lcsrank,'_',trim(lcs%label),'_',step,FILE_EXT
		endif
		call structured_io(trim(fname),IO_WRITE,gn,offset,r1=lcs%fm)

		if(lcsrank==0 .AND. LCS_VERBOSE) &
			write(*,*) ' Time h flowmap saved to: ', trim(fname)

	end subroutine write_flowmap_substep

	subroutine reconstruct_flowmap(lcs)
		implicit none
		type(lcs_t):: lcs
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
			write(*,*) 'in reconstruct_flowmap... ',trim(lcs%fm%label)

		!-----
		!Initialize temporary data:
		!-----
		gn = (/lcs%sgrid%gni,lcs%sgrid%gnj,lcs%sgrid%gnk/)
		offset = (/lcs%sgrid%offset_i,lcs%sgrid%offset_j,lcs%sgrid%offset_k/)
		call init_sr1(fm_tmp,lcs%sgrid%ni,lcs%sgrid%nj,lcs%sgrid%nk,lcs%sgrid%ng,'fm_tmp',.FALSE.)

		!-----
		!Reset the lp grid, which will be used to gather the time T flow map
		!-----
		call reset_lp(lcs%lp,lcs%sgrid%grid)
		ni = lcs%sgrid%ni; nj = lcs%sgrid%nj; nk = lcs%sgrid%nk;
		do ip = 1,lcs%lp%np
			lcs%lp%no%x(ip) =  l2i(lcs%lp%no0%i(ip),ni)
			lcs%lp%no%y(ip) =  l2j(lcs%lp%no0%i(ip),ni,nj)
			lcs%lp%no%z(ip) =  l2k(lcs%lp%no0%i(ip),ni,nj)
		enddo

		!-----
		!Loop through each of the substeps and construct the time T flowmap:
		!-----
		step = nint( mod(scfd%t_np1,lcs%T)/lcs%h)
		nstep = nint(lcs%T/lcs%h)
		if (lcs%lp%direction == FWD) then
			inc = 1
			step = step + inc
		else
			inc = -1
		endif
		!Set the lp direction to FWD, and correct after reconstruct
		lcs%lp%direction = FWD

		do i = 1,nstep
			!Load this substep (or skip, if there is no data yet:
			if (step >= nstep) step = 0
			if (step < 0) step = nstep -1
			if(FILE_EXT == '.h5') then
				write(fname,'(a,a,a,a,a,i4.4,a)'),'./',trim(TEMP_DIR),'/',trim(lcs%label),'_',step,FILE_EXT
			endif
			if(FILE_EXT == '.dat') then
				write(fname,'(a,a,a,i4.4,a,a,a,i4.4,a)'),'./',trim(TEMP_DIR),'/',lcsrank,'_',trim(lcs%label),'_',step,FILE_EXT
			endif
			INQUIRE(FILE=trim(fname), EXIST=file_exists)
			if (.NOT. file_exists) then
				!next step
				step = step + inc
				cycle
			endif
			if(lcsrank==0 .AND. LCS_VERBOSE) &
				write(*,*) i,'about to open file: ', trim(fname)
			call structured_io(trim(fname),IO_READ,gn,offset,r1=lcs%fm)	!read the flowmap

			t0 = cputimer(lcscomm,SYNC_TIMER)

			!Set Flow map BC:
			call exchange_sdata(lcs%sgrid%scomm_max_r1,r1=lcs%fm)
			call set_flowmap_bc(lcs%sgrid,lcs%fm)

			!Set the particle bc: JRF do you want to use the lcs node or the scfd node?
			!call set_lp_bc(lcs%lp,scfd%sgrid,lcs%lp%no_scfd)
			call set_lp_bc(lcs%lp,lcs%sgrid,lcs%lp%no)

			!project onto 2d face (if applicable)
			call project_lp2plane(lcs%lp,lcs%sgrid,lcs%lp%no)

			!Interp flow map to particles (store in lp%up)
			call interp_s2u_r1(lcs%lp,lcs%sgrid,lcs%lp%no,lcs%lp%up,lcs%fm)

			!Update particle positions and flow map
			np = lcs%lp%np
			lcs%lp%xp%x(1:np) = lcs%lp%xp%x(1:np) + lcs%lp%up%x(1:np)
			lcs%lp%xp%y(1:np) = lcs%lp%xp%y(1:np) + lcs%lp%up%y(1:np)
			lcs%lp%xp%z(1:np) = lcs%lp%xp%z(1:np) + lcs%lp%up%z(1:np)
			lcs%lp%dx%x(1:np) = lcs%lp%dx%x(1:np) + lcs%lp%up%x(1:np)
			lcs%lp%dx%y(1:np) = lcs%lp%dx%y(1:np) + lcs%lp%up%y(1:np)
			lcs%lp%dx%z(1:np) = lcs%lp%dx%z(1:np) + lcs%lp%up%z(1:np)

			!Exchange particles
			call exchange_lp_alltoall(lcs%lp,lcs%sgrid)

			!Track to the grid
			call track_lp2node(lcs%lp,lcs%sgrid,lcs%lp%no)

			!next step
			step = step + inc

			t1 = cputimer(lcscomm,SYNC_TIMER)
			cpu_reconstruct = cpu_reconstruct + max(t1-t0,0.0)
		enddo

		!map back to the origin:
		call exchange_lpmap(lcs%lp,lcs%fm)

		!Set Flow map BC:
		call exchange_sdata(lcs%sgrid%scomm_max_r1,r1=lcs%fm)
		call set_flowmap_bc(lcs%sgrid,lcs%fm)

		!cleanup
		call destroy_sr1(fm_tmp)
		if (inc == -1) lcs%lp%direction = BKWD  !reset direction

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
