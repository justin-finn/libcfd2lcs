module lp_motion_m
	use data_m
	use lp_tracking_m
	implicit none

	real(LCSRP),parameter:: CFL_MAX = 0.5_LCSRP

	contains

	subroutine update_lp(lp,flow)
		use comms_m
		implicit none
		!-----
		type(lp_t),pointer:: lp
		type(scfd_t):: flow
		!-----
		logical:: DONE
		real(LCSRP):: dt,lp_time
		integer:: subcycle
		!-----
		!Update the LP positions and velocities
		!based on the current scfd data. 
		!Subcycling is used, if necessary, to ensure that
		!the particle CFL number (wrt the flow grid) is not
		!greater than CFL_MAX
		!-----
		if(lcsrank==0) &
			write(*,*) 'in update_lp...', trim(lp%label)

		!-----
		!Decide on the subcycling timestep
		!-----
		call set_dt()
	
		!-----
		!Advance the particles through the flow timestep	
		!-----
		subcycle= 0
		DONE = .FALSE.
		lp_time = flow%t_n
		do while (.NOT. DONE)
			
			!Check if this should be the last subcycle
			if(lp_time + dt >= flow%t_np1) then
				dt = flow%t_np1-lp_time
				DONE = .TRUE.
			endif
			
			!Integrate for new positions/velocities
			if(lcsrank==0)&
				write(*,*) ' Starting LP subcycle',subcycle+1
			call integrate_lp()
			
			!Track lp to nearest node 
			call track_lp2node(lp,flow%sgrid)

			!Update time
			lp_time = lp_time + dt
			subcycle = subcycle + 1

		enddo

		contains
		subroutine set_dt()
			implicit none
			!-----
			real(LCSRP):: dtf(1:flow%sgrid%ni,1:flow%sgrid%nj,1:flow%sgrid%nk)
			real(LCSRP):: dtp(1:lp%np)
			integer:: ni,nj,nk
			type(sr1_t),pointer:: grid
			real(LCSRP):: mydt
			integer:: ierr
			!-----
			if(lcsrank==0 .AND. LCS_VERBOSE)&
				write(*,*) 'in set_dt...'
 
			!brevity...
			ni = flow%sgrid%ni
			nj = flow%sgrid%nj
			nk = flow%sgrid%nk
			grid => flow%sgrid%grid

			!CFL restriction based on fluid velocity:
			dtf(1:ni,1:nj,1:nk) = 0.5*CFL_MAX/(&
			 abs(flow%u_np1%x(1:ni,1:nj,1:nk)/(grid%x(2:ni+1,1:nj,1:nk) - grid%x(0:ni-1,1:nj,1:nk))) &
			+abs(flow%u_np1%y(1:ni,1:nj,1:nk)/(grid%y(1:ni,2:nj+1,1:nk) - grid%y(1:ni,0:nj-1,1:nk))) &
			+abs(flow%u_np1%z(1:ni,1:nj,1:nk)/(grid%z(1:ni,1:nj,2:nk+1) - grid%z(1:ni,1:nj,0:nk-1))) &
			)
			
			!TODO: CFL restriction based on particle velocity (for inertial particles)
			
			!Set dt based on most restrictive condition, across all procs.
			mydt = min(flow%t_np1-flow%t_n,minval(dtf))
			call MPI_ALLREDUCE(mydt,dt,1,MPI_LCSRP,MPI_MIN,lcscomm,ierr)
			
			if(lcsrank==0) &
				write(*,'(a,ES11.4,a,i5)')  'Particle DT = ', dt,', N-subcycle = ',ceiling((flow%t_np1-flow%t_n)/dt)
			
		end subroutine set_dt

		subroutine integrate_lp()
			implicit none
			!-----
			!Here, we integrate the particle equation of motion
			!using the known velocity at time level n and np1
			!-----
	
			!Euler (1st order)
			call interp_s2u_r1(lp,flow%sgrid%grid,lp%up,flow%u_n) !Interp u 
			lp%xp%x(1:lp%np) = lp%xp%x(1:lp%np) + dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = lp%xp%y(1:lp%np) + dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = lp%xp%z(1:lp%np) + dt*lp%up%z(1:lp%np)

			!Trapezoidal (1st order)

			!RK2

			!RK4

		end subroutine integrate_lp


















	end subroutine update_lp
end module lp_motion_m
