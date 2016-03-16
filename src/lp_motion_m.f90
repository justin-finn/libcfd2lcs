module lp_motion_m
	use data_m
	use lp_m
	use lp_tracking_m
	use structured_m
	use unstructured_m
	implicit none

	contains

	subroutine update_lp(lp,flow)
		use comms_m
		implicit none
		!-----
		type(lp_t),pointer:: lp
		type(scfd_t):: flow
		!-----
		real(LCSRP):: dt,lp_time
		integer:: n_subcycle,subcycle,ierr
		!-----
		!Update the LP positions and velocities.
		!Subcycling is used, if necessary, to ensure that
		!the particle CFL number (wrt the flow grid) is not
		!greater than CFL_MAX
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE) &
			write(*,*) 'in update_lp... '

		!-----
		!Decide on the subcycling timestep
		!-----
		call set_dt(dt,n_subcycle,flow,lp)

		!-----
		!Advance the particles through the flow timestep
		!-----
		lp_time = flow%t_n
		do subcycle = 1,n_subcycle
			if(LCS_VERBOSE) N_UPDATE = 1 
			if(lcsrank==0 .AND. mod(subcycle, max(N_UPDATE/max(lp%np,1),1))==0)&
				write(*,*) ' Starting LP subcycle',subcycle

			!Integrate for new positions/velocities
			call integrate_lp(flow,lp,lp_time,dt)

			!Track lp to nearest node for forward integration
			call track_lp2node(lp,flow%sgrid,lp%no_scfd)
		enddo

	end subroutine update_lp

	subroutine set_dt(dt,n_subcycle,flow,lp)
		implicit none
		!-----
		real(LCSRP):: dt
		integer:: n_subcycle
		type(scfd_t):: flow
		type(lp_t):: lp
		!-----
		integer:: ierr
		type(sr0_t):: cfl
		integer:: my_subcycle
		real(LCSRP):: dt_f,this_cfl
		integer:: ni,nj,nk,ng,npall
		character(len=64)::myfmt
		integer:: int1,int2
		!-----
		!Set the subcycling timestep based on several criteria:
		!1. dt_p <= dt_f
		!2. CFL < CFL_MAX
		!3. CFL_P < CFL_MAX (TODO, eventually, when we have inertial particles)
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in set_dt...'

		ni = flow%sgrid%ni
		nj = flow%sgrid%nj
		nk = flow%sgrid%nk
		ng = flow%sgrid%ng

		!General CFL calculation:
		!Note, we actually store 1/delta in flow%delta
		call init_sr0(cfl,flow%sgrid%ni,flow%sgrid%nj,flow%sgrid%nk,flow%sgrid%ng,'CFL')

		dt_f = flow%t_np1-flow%t_n
		cfl%r = 0.0_LCSRP
		cfl%r = cfl%r + abs(flow%u_np1%x * flow%delta%x)
		cfl%r = cfl%r + abs(flow%u_np1%y * flow%delta%y)
		cfl%r = cfl%r + abs(flow%u_np1%z * flow%delta%z)
		cfl%r = cfl%r*dt_f

		!Compute required number of subcycles, and take max across all procs
		this_cfl = maxval(cfl%r(1:ni,1:nj,1:nk))
		my_subcycle = max(ceiling(this_cfl/(CFL_MAX*lp%dt_factor)),1)
		call MPI_ALLREDUCE(my_subcycle,n_subcycle,1,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
		dt = dt_f / real(n_subcycle,LCSRP)

		!Negative dt for bkwd integration
		if(lp%direction==BKWD) dt = -1.0_LCSRP*dt

		!Some info:
		call MPI_REDUCE(lp%np,npall,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)
		if(lcsrank==0 .AND. abs(dt) > 0.0_LCSRP) then
			int1 = ceiling(log10(real(npall)))+1
			int2 = ceiling(log10(real(n_subcycle)))+1
			if(int1>10) then
				write(myfmt,'(a,i2,a,i1,a)') '(a,a,a,i',int1,',a,ES11.4,a,i',int2,')'
			else
				write(myfmt,'(a,i1,a,i1,a)') '(a,a,a,i',int1,',a,ES11.4,a,i',int2,')'
			endif
			write(*,trim(myfmt))  'LP: ',trim(lp%label),' (NP=',npall,'): DT = ',dt,', N-subcycle = ',n_subcycle

			if(lp%direction==FWD) integrations_fwd = integrations_fwd+npall*n_subcycle
			if(lp%direction==BKWD) integrations_bkwd = integrations_bkwd+npall*n_subcycle
		endif

		!cleanup
		call destroy_sr0(cfl)

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
		integer:: ip
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
			ct = max(min((t-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u) !Interp u
			!advance
			lp%xp%x(1:lp%np) = lp%xp%x(1:lp%np) + dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = lp%xp%y(1:lp%np) + dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = lp%xp%z(1:lp%np) + dt*lp%up%z(1:lp%np)

		case(TRAPEZOIDAL) !2nd order
			!Advance based on midpoint of the subcycle
			ct = max(min((t+0.5_LCSRP*dt-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u) !Interp u
			lp%xp%x(1:lp%np) = lp%xp%x(1:lp%np) + dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = lp%xp%y(1:lp%np) + dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = lp%xp%z(1:lp%np) + dt*lp%up%z(1:lp%np)

		case(RK2) !2nd order
			!advance to t_np1h
			!ct = (t-t0)/(t1-t0)
			ct = max(min((t-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u) !Interp u
			lp%xp%x(1:lp%np) = lp%xp%x(1:lp%np) + 0.5_LCSRP*dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = lp%xp%y(1:lp%np) + 0.5_LCSRP*dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = lp%xp%z(1:lp%np) + 0.5_LCSRP*dt*lp%up%z(1:lp%np)
			!Compute velocity at xnp1h,tnp1h
			!ct = (t+0.5_LCSRP*dt-t0)/(t1-t0)
			ct = max(min((t+0.5_LCSRP*dt-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u) !Interp u
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
		lp%lifetime = lp%lifetime + abs(dt)

		!Project into the plane (if 2d) and handle boundary conditions:
		!NB: Here, we use the grid associated with the flow for bc enforcement.
		call project_lp2plane(lp,scfd%sgrid,lp%no_scfd)
		call set_lp_bc(lp,scfd%sgrid,lp%no_scfd)
		
		!Make sure we stick if needed:
		do ip = 1,lp%np
			if(lp%flag%i(ip) == LP_STICK) then
				lp%xp%x(ip) = xp0%x(ip)
				lp%xp%y(ip) = xp0%y(ip)
				lp%xp%z(ip) = xp0%z(ip)
			endif
		enddo

		!Increment dx:
		lp%dx%x(1:lp%np) = lp%dx%x(1:lp%np) + (lp%xp%x(1:lp%np)-xp0%x(1:lp%np))
		lp%dx%y(1:lp%np) = lp%dx%y(1:lp%np) + (lp%xp%y(1:lp%np)-xp0%y(1:lp%np))
		lp%dx%z(1:lp%np) = lp%dx%z(1:lp%np) + (lp%xp%z(1:lp%np)-xp0%z(1:lp%np))

		!cleanup
		call destroy_ur1(xp0)
		call destroy_sr1(u)
	end subroutine integrate_lp

end module lp_motion_m
