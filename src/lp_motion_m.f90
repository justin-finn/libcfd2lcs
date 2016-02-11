module lp_motion_m
	use data_m
	use lp_tracking_m
	use structured_m
	use unstructured_m
	implicit none



	contains

	subroutine update_flowmap_sl(lcs,flow)
		use flowmap_m
		use comms_m
		implicit none
		!-----
		type(lcs_t),pointer:: lcs
		type(scfd_t):: flow
		!-----
		type(lp_t),pointer:: lp
		logical:: DONE
		real(LCSRP):: dt,lp_time,dt_factor
		type(ur1_t):: lpgrid
		type(ur1_t):: fmp
		integer:: subcycle, ip
		integer::i,j,k
		integer:: npall,ierr
		type(ui1_t):: tmp_no
		!-----
		!Perform a semi-lagrangian update of the flow map
		!Subcycling is used, if necessary, to ensure that
		!the particle CFL number (wrt the flow grid and LCS grid)
		!is not greater than CFL_MAX
		!-----
		call MPI_REDUCE(lcs%lp%np,npall,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)

		if(lcsrank==0) &
			write(*,*) 'in update_flowmap_sl...', trim(lcs%fm%label), '(',npall,' particles)'

		!-----
		!Initialize some temporary structures for integration:
		!-----
		lp => lcs%lp
		call init_ur1(lpgrid,lp%np,'ParticleGrid')
		call init_ur1(fmp,lp%np,'FM_P')
		lpgrid%x(1:lp%np)=lp%xp%x(1:lp%np)
		lpgrid%y(1:lp%np)=lp%xp%y(1:lp%np)
		lpgrid%z(1:lp%np)=lp%xp%z(1:lp%np)
		if (lcs%resolution /=0) then
			call init_ui1(tmp_no,lp%np,'TMP_NODE')
			tmp_no%x(1:lp%np) = lp%no%x(1:lp%np)
			tmp_no%y(1:lp%np) = lp%no%y(1:lp%np)
			tmp_no%z(1:lp%np) = lp%no%z(1:lp%np)
		endif

		!-----
		!Set dt, and remember to account for different lcs spacing
		!-----
		dt_factor = 1.0_LCSRP/real(max(1+lcs%resolution, 1))
		call set_dt(dt,flow,lp,dt_factor)

		!-----
		!Advance the particles through the flow timestep
		!-----
		subcycle= 0
		DONE = .FALSE.
		lp_time = flow%t_np1
		do while (.NOT. DONE)
			subcycle = subcycle + 1
			if(lcsrank==0 .AND. LCS_VERBOSE )&
				write(*,*) ' Starting SL subcycle',subcycle

			!Check if this should be the last subcycle
			if(lp_time + dt <= flow%t_n) then
				dt = flow%t_n-lp_time
				DONE = .TRUE.
			endif
			if(dt >= 0.0_LCSRP) exit

			!integrate
			!if the resolution of the scfd and lcs grids are not the same, be careful:
			if (lcs%resolution ==0) then
				call integrate_lp(flow,lp,lp_time,dt)
			else
				lp%no%x(1:lp%np) = lcs%scfd_node%x(1:lp%np)
				lp%no%y(1:lp%np) = lcs%scfd_node%y(1:lp%np)
				lp%no%z(1:lp%np) = lcs%scfd_node%z(1:lp%np)
				call integrate_lp(flow,lp,lp_time,dt)
				lp%no%x(1:lp%np) = tmp_no%x(1:lp%np)
				lp%no%y(1:lp%np) = tmp_no%y(1:lp%np)
				lp%no%z(1:lp%np) = tmp_no%z(1:lp%np)
			endif


			!Interpolate the flow map to the particles at t+dt,
			call interp_s2u_r1(lp,lcs%sgrid,fmp,lcs%fm) !Interp flow map to lp

			!Update the flow map (stored in displacement form, relative to the fixed grid).
			do ip = 1,lp%np
				lcs%fm%x(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip)) = fmp%x(ip) + (lp%xp%x(ip)-lpgrid%x(ip))
				lcs%fm%y(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip)) = fmp%y(ip) + (lp%xp%y(ip)-lpgrid%y(ip))
				lcs%fm%z(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip)) = fmp%z(ip) + (lp%xp%z(ip)-lpgrid%z(ip))
			enddo
			call exchange_sdata(lcs%sgrid%scomm_max_r1,r1=lcs%fm)

			!relocate particles back to original grid
			lp%xp%x(1:lp%np)=lpgrid%x(1:lp%np)
			lp%xp%y(1:lp%np)=lpgrid%y(1:lp%np)
			lp%xp%z(1:lp%np)=lpgrid%z(1:lp%np)
		enddo

		!-----
		!Handle backward flow map bc
		!-----
		call set_flowmap_bc(lcs%fm, lcs%sgrid)
		
		!cleanup
		call destroy_ur1(lpgrid)
		call destroy_ur1(fmp)
		call destroy_ui1(tmp_no)

	end subroutine update_flowmap_sl

	subroutine update_lp(lp,flow)
		use comms_m
		implicit none
		!-----
		type(lp_t),pointer:: lp
		type(scfd_t):: flow
		!-----
		logical:: DONE
		real(LCSRP):: dt,lp_time
		integer:: subcycle,npall,ierr
		!-----
		!Update the LP positions and velocities.
		!Subcycling is used, if necessary, to ensure that
		!the particle CFL number (wrt the flow grid) is not
		!greater than CFL_MAX
		!-----
		call MPI_REDUCE(lp%np,npall,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)
		if(lcsrank==0) &
			write(*,*) 'in update_lp...', trim(lp%label), '(',npall,' particles)'

		!-----
		!Decide on the subcycling timestep
		!-----
		call set_dt(dt,flow,lp,1.0_LCSRP)

		!-----
		!Advance the particles through the flow timestep
		!-----
		subcycle= 0
		DONE = .FALSE.
		lp_time = flow%t_n
		do while (.NOT. DONE)
			subcycle = subcycle + 1
			if(lcsrank==0 .AND. LCS_VERBOSE )&
				write(*,*) ' Starting LP subcycle',subcycle

			!Check if this should be the last subcycle
			if(lp_time + dt >= flow%t_np1) then
				dt = flow%t_np1-lp_time
				DONE = .TRUE.
			endif
			if(dt <= 0.0_LCSRP) exit

			!Integrate for new positions/velocities
			call integrate_lp(flow,lp,lp_time,dt)

			!Track lp to nearest node for forward integration
			call track_lp2node(lp,flow%sgrid)
		enddo

	end subroutine update_lp

	subroutine set_dt(dt,flow,lp,dt_factor)
		implicit none
		!-----
		real(LCSRP):: dt
		type(scfd_t):: flow
		type(lp_t):: lp
		real(LCSRP):: dt_factor
		!-----
		integer:: i,j,k,ni,nj,nk,ng
		type(sr1_t),pointer:: grid
		integer:: ierr
		type(sr0_t):: cfl
		integer:: my_subcycle,subcycle
		real(LCSRP):: dt_f
		!-----
		!Set the subcycling timestep based on several criteria:
		!1. dt_p <= dt_f
		!2. CFL < CFL_MAX
		!3. CFL_P < CFL_MAX (TODO, eventually, when we have inertial particles)
		!-----
real(LCSRP):: dtf(1:flow%sgrid%ni,1:flow%sgrid%nj,1:flow%sgrid%nk)
real(LCSRP):: mydt
logical:: debug=.false.

		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in set_dt...'

		!brevity...
		ni = flow%sgrid%ni
		nj = flow%sgrid%nj
		nk = flow%sgrid%nk
		ng = flow%sgrid%ng
		grid => flow%sgrid%grid
if(debug) then
		!Rectilinear CFL restriction:
		dtf(1:ni,1:nj,1:nk) = 0.5*CFL_MAX/(&
		 abs(flow%u_np1%x(1:ni,1:nj,1:nk)/(grid%x(2:ni+1,1:nj,1:nk) - grid%x(0:ni-1,1:nj,1:nk))) &
		+abs(flow%u_np1%y(1:ni,1:nj,1:nk)/(grid%y(1:ni,2:nj+1,1:nk) - grid%y(1:ni,0:nj-1,1:nk))) &
		+abs(flow%u_np1%z(1:ni,1:nj,1:nk)/(grid%z(1:ni,1:nj,2:nk+1) - grid%z(1:ni,1:nj,0:nk-1))) &
		)
		
		!Set dt based on most restrictive condition, across all procs.
		mydt = min(flow%t_np1-flow%t_n,minval(dtf))
		call MPI_ALLREDUCE(mydt,dt,1,MPI_LCSRP,MPI_MIN,lcscomm,ierr)

		!Allow user to pass a factor to account for
		!a grid spacing that is different from the flow grid (for the sl update of bkwd flow map)
		dt = dt*dt_factor

		!Negative dt for bkwd integration
		if(lp%direction==BKWD) dt = -1.0_LCSRP*dt

		if(lcsrank==0 .AND. abs(dt) > 0.0_LCSRP) &
			write(*,'(a,ES11.4,a,i5)')  '  particle DT = ', dt,&
				&', N-subcycle = ',ceiling(abs((flow%t_np1-flow%t_n)/dt))
		
else

		!General CFL calculation:
		!First time through, we need to compute characteristic dimensions for each node
		call init_sr0(cfl,ni,nj,nk,ng,'CFL')
		if(.NOT. allocated(flow%delta%x)) then
			call init_sr1(flow%delta,ni,nj,nk,ng,'TMP',translate=.false.)
			flow%delta%x = huge(1.0)
			flow%delta%y = huge(1.0)
			flow%delta%z = huge(1.0)
			do k = 1,nk
			do j = 1,nj
			do i = 1,ni
				flow%delta%x(i,j,k) = 0.5_LCSRP*abs(&
				maxval(grid%x(i-1:i+1,j-1:j+1,k-1:k+1))-minval(grid%x(i-1:i+1,j-1:j+1,k-1:k+1)))
				flow%delta%y(i,j,k) = 0.5_LCSRP*abs(&
				maxval(grid%y(i-1:i+1,j-1:j+1,k-1:k+1))-minval(grid%y(i-1:i+1,j-1:j+1,k-1:k+1)))
				flow%delta%z(i,j,k) = 0.5_LCSRP*abs(&
				maxval(grid%z(i-1:i+1,j-1:j+1,k-1:k+1))-minval(grid%z(i-1:i+1,j-1:j+1,k-1:k+1)))
			enddo
			enddo
			enddo
		endif
		dt_f = flow%t_np1-flow%t_n
		cfl%r = dt_f*(abs(flow%u_np1%x/flow%delta%x)+abs(flow%u_np1%y/flow%delta%y)+abs(flow%u_np1%z/flow%delta%z))

		!Compute required number of subcycles, and take max across all procs	
		my_subcycle = max(ceiling(maxval(cfl%r)/(CFL_MAX*dt_factor)),1)
		call MPI_ALLREDUCE(my_subcycle,subcycle,1,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
		dt = dt_f / real(subcycle,LCSRP)
		
		!Negative dt for bkwd integration
		if(lp%direction==BKWD) dt = -1.0_LCSRP*dt

		if(lcsrank==0 .AND. abs(dt) > 0.0_LCSRP) &
			write(*,'(a,ES11.4,a,i5)')  '  particle DT = ', dt,', N-subcycle = ',subcycle
		
		call destroy_sr0(cfl)
endif

	end subroutine set_dt

	subroutine integrate_lp(flow,lp,t,dt)
		implicit none
		!-----
		type(scfd_t):: flow
		type(lp_t):: lp
		real(LCSRP),intent(inout):: t
		real(LCSRP),intent(in):: dt
		!-----
		type(sr1_t):: u
		type(ur1_t):: xp0
		real(LCSRP):: t0,t1,ct
		!-----
		!Integrate the particle equation of motion from  t to t + dt
		!Use the known velocity at time level n and np1
		!-----

		!Some temporary structures:
		call init_sr1(u,flow%u_n%ni,flow%u_n%nj,flow%u_n%nk,flow%u_n%ng,'utmp',.FALSE.)
		call init_ur1(xp0,lp%np,'xp0')
		xp0%x(1:lp%np) = lp%xp%x(1:lp%np)
		xp0%y(1:lp%np) = lp%xp%y(1:lp%np)
		xp0%z(1:lp%np) = lp%xp%z(1:lp%np)

		t0 = flow%t_n
		t1 = flow%t_np1
		select case(INTEGRATOR)
		case (EULER) !1st order
			!Advance based on beginning of subcycle
			ct = (t-t0)/(t1-t0)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%up,u) !Interp u
			!advance
			lp%xp%x(1:lp%np) = lp%xp%x(1:lp%np) + dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = lp%xp%y(1:lp%np) + dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = lp%xp%z(1:lp%np) + dt*lp%up%z(1:lp%np)

		case(TRAPEZOIDAL) !2nd order
			!Advance based on midpoint of the subcycle
			ct = (t+0.5_LCSRP*dt-t0)/(t1-t0)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%up,u) !Interp u
			lp%xp%x(1:lp%np) = lp%xp%x(1:lp%np) + dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = lp%xp%y(1:lp%np) + dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = lp%xp%z(1:lp%np) + dt*lp%up%z(1:lp%np)

		case(RK2) !2nd order
			!advance to t_np1h
			ct = (t-t0)/(t1-t0)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%up,u) !Interp u
			lp%xp%x(1:lp%np) = lp%xp%x(1:lp%np) + 0.5_LCSRP*dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = lp%xp%y(1:lp%np) + 0.5_LCSRP*dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = lp%xp%z(1:lp%np) + 0.5_LCSRP*dt*lp%up%z(1:lp%np)
			!Compute velocity at xnp1h,tnp1h
			ct = (t+0.5_LCSRP*dt-t0)/(t1-t0)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%up,u) !Interp u
			!advance to tnp1
			lp%xp%x(1:lp%np) = xp0%x(1:lp%np) + dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = xp0%y(1:lp%np) + dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = xp0%z(1:lp%np) + dt*lp%up%z(1:lp%np)
		case (RK3)
			!TODO
			write(*,*) 'ERROR:  RK3 not yet implemented'
			CFD2LCS_ERROR = 1
			return
		case (RK4)
			!TODO
			write(*,*) 'ERROR:  RK4 not yet implemented'
			CFD2LCS_ERROR = 1
			return
		end select

		!Increment time:
		t = t + dt

		!Increment dx:
		lp%dx%x(1:lp%np) = lp%dx%x(1:lp%np) + (lp%xp%x(1:lp%np)-xp0%x(1:lp%np))
		lp%dx%y(1:lp%np) = lp%dx%y(1:lp%np) + (lp%xp%y(1:lp%np)-xp0%y(1:lp%np))
		lp%dx%z(1:lp%np) = lp%dx%z(1:lp%np) + (lp%xp%z(1:lp%np)-xp0%z(1:lp%np))

		!cleanup
		call destroy_ur1(xp0)
		call destroy_sr1(u)
	end subroutine integrate_lp

end module lp_motion_m
