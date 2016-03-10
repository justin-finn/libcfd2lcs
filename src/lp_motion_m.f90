module lp_motion_m
	use data_m
	use lp_tracking_m
	use structured_m
	use unstructured_m
	implicit none

	integer,parameter:: N_UPDATE = 1000000

	contains

	subroutine update_lp(lp,flow)
		use comms_m
		implicit none
		!-----
		type(lp_t),pointer:: lp
		type(scfd_t):: flow
		!-----
		real(LCSRP):: dt,lp_time
		integer:: n_subcycle,subcycle,npall,ierr
		!-----
		!Update the LP positions and velocities.
		!Subcycling is used, if necessary, to ensure that
		!the particle CFL number (wrt the flow grid) is not
		!greater than CFL_MAX
		!-----

		call MPI_REDUCE(lp%np,npall,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)
		if(lcsrank==0) &
			write(*,*) 'in update_lp... ',trim(lp%label), ' NP =',npall

		!-----
		!Decide on the subcycling timestep
		!-----
		call set_dt(dt,n_subcycle,flow,lp,1.0_LCSRP)

		!-----
		!Advance the particles through the flow timestep
		!-----
		lp_time = flow%t_n
		do subcycle = 1,n_subcycle
			if(lcsrank==0 .AND. mod(subcycle, max(N_UPDATE/lp%np,1))==0)&
				write(*,*) ' Starting LP subcycle',subcycle

			!Integrate for new positions/velocities
			call integrate_lp(flow,lp,lp_time,dt)

			!Track lp to nearest node for forward integration
			call track_lp2node(lp,flow%sgrid,lp%no_scfd)
		enddo

	end subroutine update_lp

	subroutine set_dt(dt,n_subcycle,flow,lp,dt_factor)
		implicit none
		!-----
		real(LCSRP):: dt
		integer:: n_subcycle
		type(scfd_t):: flow
		type(lp_t):: lp
		real(LCSRP):: dt_factor
		!-----
		integer:: i,j,k,ni,nj,nk,ng
		type(sr1_t),pointer:: grid
		integer:: ierr
		type(sr0_t):: cfl
		integer:: my_subcycle
		real(LCSRP):: dt_f,this_cfl
		integer:: im1,jm1,km1,ip1,jp1,kp1
		!-----
		!Set the subcycling timestep based on several criteria:
		!1. dt_p <= dt_f
		!2. CFL < CFL_MAX
		!3. CFL_P < CFL_MAX (TODO, eventually, when we have inertial particles)
		!-----

		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in set_dt...'

		!brevity...
		ni = flow%sgrid%ni
		nj = flow%sgrid%nj
		nk = flow%sgrid%nk
		ng = flow%sgrid%ng
		grid => flow%sgrid%grid

		!General CFL calculation:
		!First time through, we need to compute characteristic dimensions for each node
		!Note, we actually store 1/delta in flow%delta
		call init_sr0(cfl,ni,nj,nk,ng,'CFL')
		if(.NOT. allocated(flow%delta%x)) then
			call init_sr1(flow%delta,ni,nj,nk,ng,'TMP',translate=.false.)
			do k = 1,nk
			do j = 1,nj
			do i = 1,ni
				!set range
				im1 = max(i-1,1)
				jm1 = max(j-1,1)
				km1 = max(k-1,1)
				ip1 = min(i+1,ni)
				jp1 = min(j+1,nj)
				kp1 = min(k+1,nk)

				flow%delta%x(i,j,k) = 0.5_LCSRP*abs(&
				maxval(grid%x(im1:ip1,jm1:jp1,km1:kp1))-minval(grid%x(im1:ip1,jm1:jp1,km1:kp1)))
				flow%delta%y(i,j,k) = 0.5_LCSRP*abs(&
				maxval(grid%y(im1:ip1,jm1:jp1,km1:kp1))-minval(grid%y(im1:ip1,jm1:jp1,km1:kp1)))
				flow%delta%z(i,j,k) = 0.5_LCSRP*abs(&
				maxval(grid%z(im1:ip1,jm1:jp1,km1:kp1))-minval(grid%z(im1:ip1,jm1:jp1,km1:kp1)))
				!protect against 0 dist, store 1/dx,1/dy,1/dz
				if(flow%delta%x(i,j,k) > 0.0_LCSRP) then
					flow%delta%x(i,j,k) = 1.0_LCSRP / flow%delta%x(i,j,k)
				else
					flow%delta%x(i,j,k)= 0.0_LCSRP
				endif
				if(flow%delta%y(i,j,k) > 0.0_LCSRP) then
					flow%delta%y(i,j,k) = 1.0_LCSRP / flow%delta%y(i,j,k)
				else
					flow%delta%y(i,j,k)= 0.0_LCSRP
				endif
				if(flow%delta%z(i,j,k) > 0.0_LCSRP) then
					flow%delta%z(i,j,k) = 1.0_LCSRP / flow%delta%z(i,j,k)
				else
					flow%delta%z(i,j,k)= 0.0_LCSRP
				endif
			enddo
			enddo
			enddo
		endif
		dt_f = flow%t_np1-flow%t_n
		cfl%r = 0.0_LCSRP
		cfl%r = cfl%r + abs(flow%u_np1%x * flow%delta%x)
		cfl%r = cfl%r + abs(flow%u_np1%y * flow%delta%y)
		cfl%r = cfl%r + abs(flow%u_np1%z * flow%delta%z)
		cfl%r = cfl%r*dt_f

		!Compute required number of subcycles, and take max across all procs
		this_cfl = maxval(cfl%r)
		my_subcycle = max(ceiling(this_cfl/(CFL_MAX*dt_factor)),1)
		call MPI_ALLREDUCE(my_subcycle,n_subcycle,1,MPI_INTEGER,MPI_MAX,lcscomm,ierr)
		dt = dt_f / real(n_subcycle,LCSRP)

		!Negative dt for bkwd integration
		if(lp%direction==BKWD) dt = -1.0_LCSRP*dt

		if(lcsrank==0 .AND. abs(dt) > 0.0_LCSRP) &
			write(*,'(a,ES11.4,a,i5,a,ES11.4)')  '  particle DT = ', dt,&
			', N-subcycle = ',n_subcycle,', CFL = ',this_cfl

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
	
		!Make sure we stick if needed:	
		do ip = 1,lp%np
			if(lp%flag%i(ip) == LP_STICK) then
				lp%xp%x(ip) = xp0%x(ip)
				lp%xp%y(ip) = xp0%y(ip)
				lp%xp%z(ip) = xp0%z(ip)
			endif
		enddo

		!Project into the plane (if 2d):
		call project_lp2plane(lp,flow%sgrid,lp%no_scfd)

		!Handle particle boundary conditions
		!ONLY if direction = FWD??  JRF
		call set_lp_bc(lp,flow%sgrid,lp%no_scfd)

		!Increment dx:
		lp%dx%x(1:lp%np) = lp%dx%x(1:lp%np) + (lp%xp%x(1:lp%np)-xp0%x(1:lp%np))
		lp%dx%y(1:lp%np) = lp%dx%y(1:lp%np) + (lp%xp%y(1:lp%np)-xp0%y(1:lp%np))
		lp%dx%z(1:lp%np) = lp%dx%z(1:lp%np) + (lp%xp%z(1:lp%np)-xp0%z(1:lp%np))

		!cleanup
		call destroy_ur1(xp0)
		call destroy_sr1(u)
	end subroutine integrate_lp


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
		!ip = lp%np  !JRF:  Potentially sinister bug here previously where ip was set instead of np
		np = lp%np  !JRF: this is what it should be


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

end module lp_motion_m
