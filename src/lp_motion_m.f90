!
!Copyright (C) 2015-2016, Justin R. Finn.  All rights reserved.
!libcfd2lcs is distributed is under the terms of the GNU General Public License
!
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
		integer:: ni,nj,nk,ng
		integer(8):: npall
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
		call MPI_REDUCE(int(lp%np,8),npall,1,MPI_INTEGER8,MPI_SUM,0,lcscomm,ierr)
		if(lcsrank==0 .AND. abs(dt) > 0.0_LCSRP) then
			int1 = ceiling(log10(real(npall)))+1
			int2 = ceiling(log10(real(n_subcycle)))+1
			if(int1>=10) then
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
		type(ur1_t):: xn, k, xnp1
		real(LCSRP):: t0,t1,ct
		integer:: ip
		real(LCSRP), parameter:: TWO       = 2.0_LCSRP
		real(LCSRP), parameter:: ONEHALF   = 1.0_LCSRP/2.0_LCSRP
		real(LCSRP), parameter:: ONETHIRD  = 1.0_LCSRP/3.0_LCSRP
		real(LCSRP), parameter:: TWOTHIRD  = 2.0_LCSRP/3.0_LCSRP
		real(LCSRP), parameter:: ONESIXTH  = 1.0_LCSRP/6.0_LCSRP
		!-----
		!Integrate the particle equation of motion from  t to t + dt
		!Use the known velocity at time level n and np1
		!Available methods:  EULER, TRAPEZOIDAL, RK2, RK3, RK4
		!-----

		!Some temporary structures:
		call init_sr1(u,flow%u_n%ni,flow%u_n%nj,flow%u_n%nk,flow%u_n%ng,'utmp',.FALSE.)
		call init_ur1(xn,lp%np,'xn')
		xn%x(1:lp%np) = lp%xp%x(1:lp%np)
		xn%y(1:lp%np) = lp%xp%y(1:lp%np)
		xn%z(1:lp%np) = lp%xp%z(1:lp%np)
	
		!For higher order, we need extra storage:
		if(INTEGRATOR == RK3 .OR. INTEGRATOR == RK4) then	
			call init_ur1(k,lp%np,'k')
			call init_ur1(xnp1,lp%np,'xnp1')
		endif

		t0 = flow%t_n
		t1 = flow%t_np1
		select case(INTEGRATOR)
		case (EULER) !1st order
			!Advance based on beginning of subcycle
			ct = max(min((t-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u)
			!advance
			lp%xp%x(1:lp%np) = lp%xp%x(1:lp%np) + dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = lp%xp%y(1:lp%np) + dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = lp%xp%z(1:lp%np) + dt*lp%up%z(1:lp%np)

		case(TRAPEZOIDAL) !2nd order
			!Advance based on midpoint of the subcycle
			ct = max(min((t+ONEHALF*dt-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u)
			lp%xp%x(1:lp%np) = lp%xp%x(1:lp%np) + dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = lp%xp%y(1:lp%np) + dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = lp%xp%z(1:lp%np) + dt*lp%up%z(1:lp%np)

		case(RK2) !2nd order
			!advance to t_np1h
			ct = max(min((t-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u)
			lp%xp%x(1:lp%np) = lp%xp%x(1:lp%np) + ONEHALF*dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = lp%xp%y(1:lp%np) + ONEHALF*dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = lp%xp%z(1:lp%np) + ONEHALF*dt*lp%up%z(1:lp%np)
			!Compute velocity at xnp1h,tnp1h
			ct = max(min((t+ONEHALF*dt-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u)
			!advance to tnp1
			lp%xp%x(1:lp%np) = xn%x(1:lp%np) + dt*lp%up%x(1:lp%np)
			lp%xp%y(1:lp%np) = xn%y(1:lp%np) + dt*lp%up%y(1:lp%np)
			lp%xp%z(1:lp%np) = xn%z(1:lp%np) + dt*lp%up%z(1:lp%np)

		case (RK3)
			!IC
			xnp1%x = xn%x
			xnp1%y = xn%y
			xnp1%z = xn%z
			!Evaluate K1
			ct = max(min((t-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u)
			k%x(1:lp%np) = dt*lp%up%x(1:lp%np)  !k1
			k%y(1:lp%np) = dt*lp%up%y(1:lp%np)  !k1
			k%z(1:lp%np) = dt*lp%up%z(1:lp%np)  !k1
			xnp1%x(1:lp%np) = xnp1%x(1:lp%np)+ONESIXTH*k%x(1:lp%np)	
			xnp1%y(1:lp%np) = xnp1%y(1:lp%np)+ONESIXTH*k%y(1:lp%np)	
			xnp1%z(1:lp%np) = xnp1%z(1:lp%np)+ONESIXTH*k%z(1:lp%np)
			!Evaluate K2
			!NB: Perform half update of using k1, compute k2, then update again.
			ct = max(min((t+ONEHALF*dt-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			lp%xp%x(1:lp%np) = xn%x(1:lp%np) + ONEHALF*k%x(1:lp%np)
			lp%xp%y(1:lp%np) = xn%y(1:lp%np) + ONEHALF*k%y(1:lp%np)
			lp%xp%z(1:lp%np) = xn%z(1:lp%np) + ONEHALF*k%z(1:lp%np)
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u)
			lp%xp%x(1:lp%np) = xn%x(1:lp%np) - k%x(1:lp%np)
			lp%xp%y(1:lp%np) = xn%y(1:lp%np) - k%y(1:lp%np)
			lp%xp%z(1:lp%np) = xn%z(1:lp%np) - k%z(1:lp%np)
			k%x(1:lp%np) = dt*lp%up%x(1:lp%np)  !k2
			k%y(1:lp%np) = dt*lp%up%y(1:lp%np)  !k2
			k%z(1:lp%np) = dt*lp%up%z(1:lp%np)  !k2
			lp%xp%x(1:lp%np) = xn%x(1:lp%np) + TWO*k%x(1:lp%np)
			lp%xp%y(1:lp%np) = xn%y(1:lp%np) + TWO*k%y(1:lp%np)
			lp%xp%z(1:lp%np) = xn%z(1:lp%np) + TWO*k%z(1:lp%np)
			xnp1%x(1:lp%np) = xnp1%x(1:lp%np)+TWOTHIRD*k%x(1:lp%np)	
			xnp1%y(1:lp%np) = xnp1%y(1:lp%np)+TWOTHIRD*k%y(1:lp%np)	
			xnp1%z(1:lp%np) = xnp1%z(1:lp%np)+TWOTHIRD*k%z(1:lp%np)
			!Evaluate K3.  Note particle already at correct position
			ct = max(min((t+dt-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u)
			k%x(1:lp%np) = dt*lp%up%x(1:lp%np)  !k3
			k%y(1:lp%np) = dt*lp%up%y(1:lp%np)  !k3
			k%z(1:lp%np) = dt*lp%up%z(1:lp%np)  !k3
			!Set the final position at tnp1:
			lp%xp%x(1:lp%np) = xnp1%x(1:lp%np)+ONESIXTH*k%x(1:lp%np)	
			lp%xp%y(1:lp%np) = xnp1%y(1:lp%np)+ONESIXTH*k%y(1:lp%np)	
			lp%xp%z(1:lp%np) = xnp1%z(1:lp%np)+ONESIXTH*k%z(1:lp%np)	

		case (RK4) !4th order
			!IC
			xnp1%x = xn%x
			xnp1%y = xn%y
			xnp1%z = xn%z
			!Evaluate K1
			ct = max(min((t-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u)
			k%x(1:lp%np) = dt*lp%up%x(1:lp%np)  !k1
			k%y(1:lp%np) = dt*lp%up%y(1:lp%np)  !k1
			k%z(1:lp%np) = dt*lp%up%z(1:lp%np)  !k1
			xnp1%x(1:lp%np) = xnp1%x(1:lp%np)+ONESIXTH*k%x(1:lp%np)	
			xnp1%y(1:lp%np) = xnp1%y(1:lp%np)+ONESIXTH*k%y(1:lp%np)	
			xnp1%z(1:lp%np) = xnp1%z(1:lp%np)+ONESIXTH*k%z(1:lp%np)
			!Evaluate K2
			ct = max(min((t+ONEHALF*dt-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			lp%xp%x(1:lp%np) = xn%x(1:lp%np) + ONEHALF*k%x(1:lp%np)
			lp%xp%y(1:lp%np) = xn%y(1:lp%np) + ONEHALF*k%y(1:lp%np)
			lp%xp%z(1:lp%np) = xn%z(1:lp%np) + ONEHALF*k%z(1:lp%np)
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u)
			k%x(1:lp%np) = dt*lp%up%x(1:lp%np)  !k2
			k%y(1:lp%np) = dt*lp%up%y(1:lp%np)  !k2
			k%z(1:lp%np) = dt*lp%up%z(1:lp%np)  !k2
			xnp1%x(1:lp%np) = xnp1%x(1:lp%np)+ONETHIRD*k%x(1:lp%np)	
			xnp1%y(1:lp%np) = xnp1%y(1:lp%np)+ONETHIRD*k%y(1:lp%np)	
			xnp1%z(1:lp%np) = xnp1%z(1:lp%np)+ONETHIRD*k%z(1:lp%np)
			!Evaluate K3
			ct = max(min((t+ONEHALF*dt-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			lp%xp%x(1:lp%np) = xn%x(1:lp%np) + ONEHALF*k%x(1:lp%np)
			lp%xp%y(1:lp%np) = xn%y(1:lp%np) + ONEHALF*k%y(1:lp%np)
			lp%xp%z(1:lp%np) = xn%z(1:lp%np) + ONEHALF*k%z(1:lp%np)
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u)
			k%x(1:lp%np) = dt*lp%up%x(1:lp%np)  !k3
			k%y(1:lp%np) = dt*lp%up%y(1:lp%np)  !k3
			k%z(1:lp%np) = dt*lp%up%z(1:lp%np)  !k3
			xnp1%x(1:lp%np) = xnp1%x(1:lp%np)+ONETHIRD*k%x(1:lp%np)	
			xnp1%y(1:lp%np) = xnp1%y(1:lp%np)+ONETHIRD*k%y(1:lp%np)	
			xnp1%z(1:lp%np) = xnp1%z(1:lp%np)+ONETHIRD*k%z(1:lp%np)
			!Evaluate K4
			ct = max(min((t+dt-t0)/(t1-t0),1.0_LCSRP),0.0_LCSRP)
			u%x = (1.0_LCSRP-ct)*flow%u_n%x + ct*flow%u_np1%x
			u%y = (1.0_LCSRP-ct)*flow%u_n%y + ct*flow%u_np1%y
			u%z = (1.0_LCSRP-ct)*flow%u_n%z + ct*flow%u_np1%z
			lp%xp%x(1:lp%np) = xn%x(1:lp%np) + k%x(1:lp%np)
			lp%xp%y(1:lp%np) = xn%y(1:lp%np) + k%y(1:lp%np)
			lp%xp%z(1:lp%np) = xn%z(1:lp%np) + k%z(1:lp%np)
			call interp_s2u_r1(lp,flow%sgrid,lp%no_scfd,lp%up,u)
			k%x(1:lp%np) = dt*lp%up%x(1:lp%np)  !k4
			k%y(1:lp%np) = dt*lp%up%y(1:lp%np)  !k4
			k%z(1:lp%np) = dt*lp%up%z(1:lp%np)  !k4
			!Set the final position at tnp1:
			lp%xp%x(1:lp%np) = xnp1%x(1:lp%np)+ONESIXTH*k%x(1:lp%np)	
			lp%xp%y(1:lp%np) = xnp1%y(1:lp%np)+ONESIXTH*k%y(1:lp%np)	
			lp%xp%z(1:lp%np) = xnp1%z(1:lp%np)+ONESIXTH*k%z(1:lp%np)	
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
				lp%xp%x(ip) = xn%x(ip)
				lp%xp%y(ip) = xn%y(ip)
				lp%xp%z(ip) = xn%z(ip)
			endif
		enddo

		!Increment dx:
		lp%dx%x(1:lp%np) = lp%dx%x(1:lp%np) + (lp%xp%x(1:lp%np)-xn%x(1:lp%np))
		lp%dx%y(1:lp%np) = lp%dx%y(1:lp%np) + (lp%xp%y(1:lp%np)-xn%y(1:lp%np))
		lp%dx%z(1:lp%np) = lp%dx%z(1:lp%np) + (lp%xp%z(1:lp%np)-xn%z(1:lp%np))

		!cleanup
		call destroy_ur1(xn)
		call destroy_sr1(u)
		if(INTEGRATOR == RK3 .OR. INTEGRATOR == RK4) then	
			call destroy_ur1(xnp1)
			call destroy_ur1(k)
		endif
	end subroutine integrate_lp

end module lp_motion_m
