module lp_tracking_m
	use data_m
	use structured_m
	use unstructured_m
	use comms_m
	use bspline_oo_module
	use sgrid_m
	implicit none

	contains

	recursive subroutine track_lp2node(lp,sgrid,no)
		implicit none
		!-----
		type(lp_t):: lp
		type(sgrid_t):: sgrid
		type(ui1_t)::no
		!-----
		integer:: ip,ni,nj,nk,ng
		integer:: inbr,i,j,k,ierr
		real(LCSRP),allocatable:: rsq_min(:),rsq(:)
		real(LCSRP),allocatable:: xg(:),yg(:),zg(:)
		integer,allocatable:: ioff(:),joff(:),koff(:)
		logical:: my_do_recursion, do_recursion
		integer,save::  ntrack=0
		integer:: max_track
		integer::my_nmissing=0
		!-----
		!Find the nearest i,j,k node for each LP
		!-----

		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in track_lp2no... ', trim(lp%label),' => ',trim(sgrid%label),ntrack

		!-----
		!Check if a recursive search is needed:
		!-----
		if(lp%recursive_tracking) then
			do_recursion = .true.
		else
			do_recursion = .false.
		endif
		ntrack = ntrack+1

		!Brevity...
		ni = sgrid%grid%ni
		nj = sgrid%grid%nj
		nk = sgrid%grid%nk
		ng = sgrid%grid%ng

		!allocate
		allocate(rsq(lp%np))
		allocate(rsq_min(lp%np))
		allocate(xg(lp%np))
		allocate(yg(lp%np))
		allocate(zg(lp%np))
		allocate(ioff(lp%np))
		allocate(joff(lp%np))
		allocate(koff(lp%np))

		!-----
		!Using the last known particle node, do a local neighborhood search to
		!update the node based on current particle position
		!-----
		rsq_min(1:lp%np) = huge(0.0_LCSRP)
		do inbr = 1,N_NBR

			!Dont search for neighbors in 3d if grid is only 2d:
			if(sgrid%gni==1 .and. NBR_OFFSET(1,inbr) /=0) cycle
			if(sgrid%gnj==1 .and. NBR_OFFSET(2,inbr) /=0) cycle
			if(sgrid%gnk==1 .and. NBR_OFFSET(3,inbr) /=0) cycle

			!gather all the local grid points to a buffer:
			do ip = 1,lp%np
				!JRF: MIN/MAX are to protect for backward fm when nearest cfd node
				!could be a ghost.
				i = max(min(no%x(ip)+NBR_OFFSET(1,inbr),ni+ng),1-ng)
				j = max(min(no%y(ip)+NBR_OFFSET(2,inbr),nj+ng),1-ng)
				k = max(min(no%z(ip)+NBR_OFFSET(3,inbr),nk+ng),1-ng)
				xg(ip) = sgrid%grid%x(i,j,k)
				yg(ip) = sgrid%grid%y(i,j,k)
				zg(ip) = sgrid%grid%z(i,j,k)
			enddo
			!Compute the offset to the closest node.  These loops should vectorize:
			do ip = 1,lp%np
				rsq(ip) = &
					+(xg(ip)-lp%xp%x(ip))*(xg(ip)-lp%xp%x(ip)) &
					+(yg(ip)-lp%xp%y(ip))*(yg(ip)-lp%xp%y(ip)) &
					+(zg(ip)-lp%xp%z(ip))*(zg(ip)-lp%xp%z(ip))
			enddo
			where(rsq < rsq_min)
				rsq_min = rsq
				ioff = NBR_OFFSET(1,inbr)
				joff = NBR_OFFSET(2,inbr)
				koff = NBR_OFFSET(3,inbr)
			end where
		enddo
		!update the node
		no%x(1:lp%np) = no%x(1:lp%np)+ioff(1:lp%np)
		no%y(1:lp%np) = no%y(1:lp%np)+joff(1:lp%np)
		no%z(1:lp%np) = no%z(1:lp%np)+koff(1:lp%np)


		!-----
		!Check if we will need to call again.
		!All processors must call together.
		!-----
		if(lp%recursive_tracking) then
			!Count how many not found
			my_nmissing=0
			do ip = 1,lp%np
				if(ioff(ip) /=0 .or. joff(ip) /=0 .or. koff(ip)/=0)then
				my_nmissing= my_nmissing+1
				endif
			enddo
			if(LCS_VERBOSE)&
				write(*,*) 'lcsrank[',lcsrank,'] before iteration',ntrack, 'missing ',my_nmissing,'lp'

			!ensure we dont end up in an infinite tracking loop
			max_track = 4*max(sgrid%ni,sgrid%nj,sgrid%nk)
			if(ntrack>max_track) then
				do ip = 1,lp%np
					if( ioff(ip)/=0 .or. joff(ip)/=0 .or. koff(ip) /=0) then
						write(*,*) 'lcsrank[',lcsrank,'] cant find ip=',ip,'proc,node',lp%proc0%i(ip),lp%no0%i(ip)
						write(*,*) 'lcsrank[',lcsrank,'] location =',lp%xp%x(ip),lp%xp%y(ip),lp%xp%z(ip)
						write(*,*) 'lcsrank[',lcsrank,'] last node =',no%x(ip),no%y(ip),no%z(ip)
						write(*,*) 'lcsrank[',lcsrank,'] ioff,joff,koff =',ioff(ip),joff(ip),koff(ip)
					endif
				enddo
				write(*,*) 'ERROR: lcsrank[',lcsrank,'] cannot locate all particles after',ntrack,' tracking iterations.'
				CFD2LCS_ERROR = 1
				return
			endif
			!If particle has not moved, we can stop.
			if(any(abs(ioff)>0) .OR. any(abs(joff)>0) .OR. any(abs(koff)>0))then
				my_do_recursion = .true.
			else
				my_do_recursion = .false.
			endif
			call MPI_ALLREDUCE(my_do_recursion,do_recursion,1,MPI_LOGICAL,MPI_LOR,lcscomm,ierr)
		else
			ntrack = 0
		endif

		!-----
		!If this is a forward trajectory exchange particles
		!-----
		if(lp%direction == FWD) then
			call exchange_lpdata(lp,sgrid,no)
		endif

		!-----
		!If needed, cleanup memory and call again:
		!-----
		if(do_recursion) then
			deallocate(rsq)
			deallocate(rsq_min)
			deallocate(xg)
			deallocate(yg)
			deallocate(zg)
			deallocate(ioff)
			deallocate(joff)
			deallocate(koff)
			call track_lp2node(lp,sgrid,no)
		else
			if(lcsrank==0 .and. lp%recursive_tracking) &
				write(*,*) ' Tracking all particles required ',ntrack,' iterations'
			ntrack = 0

			!Can now turn off recursive tracking (for everyone)
			lp%recursive_tracking = .false.
		endif

	end subroutine track_lp2node

	subroutine set_lp_bc(lp,sgrid,no)
		implicit none
		!-----
		type(lp_t):: lp
		type(sgrid_t):: sgrid
		type(ui1_t):: no  !the node index associated with sgrid
		!-----
		integer:: ip,ierr
		integer:: nbounce, nstick, my_nbounce, my_nstick
		integer:: i,j,k
		integer:: ni,nj,nk,ng
		real(LCSRP):: norm(3), p2no(3),dist
		!-----
		!Handle particle boundary conditions
		!-----

		ni = sgrid%ni
		nj = sgrid%nj
		nk = sgrid%nk
		ng = sgrid%ng

		my_nbounce = 0
		my_nstick = 0
		do ip = 1,lp%np

			!Enusre an IB node:
			i = max(min(no%x(ip),ni),1)
			j = max(min(no%y(ip),nj),1)
			k = max(min(no%z(ip),nk),1)
			no%x(ip) = i
			no%y(ip) = j
			no%z(ip) = k

			select case(sgrid%bcflag%i(i,j,k))
				case(LCS_INTERNAL,LCS_2D)
					cycle !do nothing (for LCS_2D, we handle in project_lp2face)
				case(LCS_MASK)
					!Project the particle back to the boundary:
					p2no(1) = lp%xp%x(ip) -sgrid%grid%x(i,j,k)
					p2no(2) = lp%xp%y(ip) -sgrid%grid%y(i,j,k)
					p2no(3) = lp%xp%z(ip) -sgrid%grid%z(i,j,k)
					norm(1) = sgrid%bcnorm%x(i,j,k)
					norm(2) = sgrid%bcnorm%y(i,j,k)
					norm(3) = sgrid%bcnorm%z(i,j,k)
					dist = dot_product(p2no,norm)
					if(dist > 0.0_LCSRP) then
						lp%xp%x(ip) = lp%xp%x(ip) - norm(1)*dist
						lp%xp%y(ip) = lp%xp%y(ip) - norm(2)*dist
						lp%xp%z(ip) = lp%xp%z(ip) - norm(3)*dist
						my_nbounce = my_nbounce + 1
					endif

				case(LCS_SLIP,LCS_WALL)
					!Project the particle back to the boundary:
					if(i>1 .and. j>1 .and. k>1 .and. i<ni .and. j<nj .and. k<nk) cycle
					p2no(1) = lp%xp%x(ip) -sgrid%grid%x(i,j,k)
					p2no(2) = lp%xp%y(ip) -sgrid%grid%y(i,j,k)
					p2no(3) = lp%xp%z(ip) -sgrid%grid%z(i,j,k)
					norm(1) = sgrid%bcnorm%x(i,j,k)
					norm(2) = sgrid%bcnorm%y(i,j,k)
					norm(3) = sgrid%bcnorm%z(i,j,k)
					dist = dot_product(p2no,norm)
					if(dist > 0.0_LCSRP) then
						lp%xp%x(ip) = lp%xp%x(ip) - norm(1)*dist
						lp%xp%y(ip) = lp%xp%y(ip) - norm(2)*dist
						lp%xp%z(ip) = lp%xp%z(ip) - norm(3)*dist
						my_nbounce = my_nbounce + 1
					endif

				case(LCS_INFLOW,LCS_OUTFLOW)
					!Stick the particle to the nearest IB node
					if(i>1 .and. j>1 .and. k>1 .and. i<ni .and. j<nj .and. k<nk) cycle
					lp%flag%i(ip) = LP_STICK
					lp%xp%x(ip) = sgrid%grid%x(i,j,k)
					lp%xp%y(ip) = sgrid%grid%y(i,j,k)
					lp%xp%z(ip) = sgrid%grid%z(i,j,k)
					my_nstick = my_nstick+1
				case default
					write(*,*) 'lcsrank[',lcsrank,'] ERROR: Unknown bcflag:',sgrid%bcflag%i(i,j,k)
					CFD2LCS_ERROR=1
			end select
		enddo

		!count bc encounters:
		if(LCS_VERBOSE) then
			call MPI_REDUCE(my_nstick,nstick,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)
			if (lcsrank==0) write(*,*) 'handled',nstick,' stick bc'
			call MPI_REDUCE(my_nbounce,nbounce,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)
			if (lcsrank==0) write(*,*) 'handled',nbounce,' bounce bc'
		endif

	end subroutine set_lp_bc

	subroutine interp_s2u_r1(lp,sgrid,no,ur1,sr1)
		use gradient_m
		implicit none
		!-----
		type(lp_t):: lp
		type(sgrid_t):: sgrid
		type(ui1_t):: no  !lp node index associated with sgrid
		type(ur1_t):: ur1
		type(sr1_t):: sr1
		!-----
		integer:: ip,np,ni,nj,nk,ng
		integer:: this_interpolator
		!Interpolation using B-splines forrectilinear grids:
		type(bspline_3d):: bsp_x,bsp_y,bsp_z
		real(LCSRP),allocatable:: xg(:),yg(:),zg(:)
		integer:: iflag,idx,idy,idz
		integer:: order_x,order_y,order_z
		integer:: INTERPOLATION_ORDER
		!Trilinear for rectilinear grids:
		integer:: i0,j0,k0,i1,j1,k1
		type(ur1_t):: t,mt,x0,x1
		type(ur1_t):: f0,f1,f2,f3,f4,f5,f6,f7
		!Taylor series for non-rectilinear grids:(with and without limiter)
		type(sr2_t):: gradsr1
		type(ur1_t):: r1node, delta
		type(ur2_t):: gradsr1_p
		type(ur1_t):: r1min,r1max
		!-----
		!Interpolate structured data to unstructured pts
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in interp_s2u_r1... ',trim(sr1%label),' => ',trim(ur1%label)

		!brevity:
		ni = sgrid%ni
		nj = sgrid%nj
		nk = sgrid%nk
		ng = sgrid%ng
		np = lp%np

		!-----
		!For non-rectilinear grids, we have to use the TSE or TSE_LIMIT schemes
		!(until you implement something better...)
		!-----
		this_interpolator = INTERPOLATOR
		if(.NOT. sgrid%rectilinear) then
			if(this_interpolator /= TSE_LIMIT) then
				this_interpolator = TSE
			endif
		endif

		select case(this_interpolator)

		case(NEAREST_NBR)
			!-----
			!Zeroth order, nearest node interp:
			!-----
			do ip = 1,np
				ur1%x(ip) = sr1%x(no%x(ip),no%y(ip),no%z(ip))
				ur1%y(ip) = sr1%y(no%x(ip),no%y(ip),no%z(ip))
				ur1%z(ip) = sr1%z(no%x(ip),no%y(ip),no%z(ip))
			enddo

		case(LINEAR)
			!-----
			!First order, trilinear interpolation: (2nd order accurate)
			!-----
			allocate(xg(1-ng:ni+ng))
			allocate(yg(1-ng:nj+ng))
			allocate(zg(1-ng:nk+ng))
			call init_ur1(x0,np,'X0')
			call init_ur1(x1,np,'X0')
			call init_ur1(t,np,'T')
			call init_ur1(mt,np,'1MINUST')
			call init_ur1(f0,np,'F0')
			call init_ur1(f1,np,'F1')
			call init_ur1(f2,np,'F2')
			call init_ur1(f3,np,'F3')
			call init_ur1(f4,np,'F4')
			call init_ur1(f5,np,'F5')
			call init_ur1(f6,np,'F6')
			call init_ur1(f7,np,'F7')

			xg(1-ng:ni+ng) =sgrid%grid%x(1-ng:ni+ng,1,1)
			yg(1-ng:nj+ng) =sgrid%grid%y(1,1-ng:nj+ng,1)
			zg(1-ng:nk+ng) =sgrid%grid%z(1,1,1-ng:nk+ng)

			!Special consideration for 2D:  JRF:  DO you really want/need to do this?
			if(sgrid%gni == 1) then
				xg(0) = xg(1) - 1.0;
				xg(2) = xg(1) + 1.0;
			endif
			if(sgrid%gnj == 1) then
				yg(0) = yg(1) - 1.0;
				yg(2) = yg(1) + 1.0;
			endif
			if(sgrid%gnk == 1) then
				zg(0) = zg(1) - 1.0;
				zg(2) = zg(1) + 1.0;
			endif


			!Gather scattered data into vectors
			do ip = 1,np
				if(lp%xp%x(ip) > xg(no%x(ip))) then
					i0 = no%x(ip); i1 = no%x(ip)+1
				else
					i0 = no%x(ip)-1; i1 = no%x(ip)
				endif
				if(lp%xp%y(ip) > yg(no%y(ip))) then
					j0 = no%y(ip); j1 = no%y(ip)+1
				else
					j0 = no%y(ip)-1; j1 = no%y(ip)
				endif
				if(lp%xp%z(ip) > zg(no%z(ip))) then
					k0 = no%z(ip); k1 = no%z(ip)+1
				else
					k0 = no%z(ip)-1; k1 = no%z(ip)
				endif
				x0%x(ip)= xg(i0)
				x1%x(ip)= xg(i1)
				x0%y(ip)= yg(j0)
				x1%y(ip)= yg(j1)
				x0%z(ip)= zg(k0)
				x1%z(ip)= zg(k1)
				!node 0
				f0%x(ip) = sr1%x(i0,j0,k0)
				f0%y(ip) = sr1%y(i0,j0,k0)
				f0%z(ip) = sr1%z(i0,j0,k0)
				!node 1
				f1%x(ip) = sr1%x(i1,j0,k0)
				f1%y(ip) = sr1%y(i1,j0,k0)
				f1%z(ip) = sr1%z(i1,j0,k0)
				!node 2
				f2%x(ip) = sr1%x(i0,j1,k0)
				f2%y(ip) = sr1%y(i0,j1,k0)
				f2%z(ip) = sr1%z(i0,j1,k0)
				!node 3
				f3%x(ip) = sr1%x(i1,j1,k0)
				f3%y(ip) = sr1%y(i1,j1,k0)
				f3%z(ip) = sr1%z(i1,j1,k0)
				!node 4
				f4%x(ip) = sr1%x(i0,j0,k1)
				f4%y(ip) = sr1%y(i0,j0,k1)
				f4%z(ip) = sr1%z(i0,j0,k1)
				!node 5
				f5%x(ip) = sr1%x(i1,j0,k1)
				f5%y(ip) = sr1%y(i1,j0,k1)
				f5%z(ip) = sr1%z(i1,j0,k1)
				!node 6
				f6%x(ip) = sr1%x(i0,j1,k1)
				f6%y(ip) = sr1%y(i0,j1,k1)
				f6%z(ip) = sr1%z(i0,j1,k1)
				!node 7
				f7%x(ip) = sr1%x(i1,j1,k1)
				f7%y(ip) = sr1%y(i1,j1,k1)
				f7%z(ip) = sr1%z(i1,j1,k1)
			enddo

			!Now the computations (vectorized)
			t%x(1:np) = (lp%xp%x(1:np) - x0%x(1:np)) / (x1%x(1:np)-x0%x(1:np))
			t%y(1:np) = (lp%xp%y(1:np) - x0%y(1:np)) / (x1%y(1:np)-x0%y(1:np))
			t%z(1:np) = (lp%xp%z(1:np) - x0%z(1:np)) / (x1%z(1:np)-x0%z(1:np))
			mt%x = 1.0_LCSRP-t%x
			mt%y = 1.0_LCSRP-t%y
			mt%z = 1.0_LCSRP-t%z
			ur1%x =	mt%x*mt%y*mt%z*f0%x + t%x*mt%y*mt%z*f1%x + mt%x*t%y*mt%z*f2%x + &
				t%x*t%y*mt%z*f3%x +	mt%x*mt%y*t%z*f4%x + t%x*mt%y*t%z*f5%x + &
				mt%x*t%y*t%z*f6%x + t%x*t%y*t%z*f7%x
			ur1%y =	mt%x*mt%y*mt%z*f0%y + t%x*mt%y*mt%z*f1%y + mt%x*t%y*mt%z*f2%y + &
				t%x*t%y*mt%z*f3%y +	mt%x*mt%y*t%z*f4%y + t%x*mt%y*t%z*f5%y + &
				mt%x*t%y*t%z*f6%y + t%x*t%y*t%z*f7%y
			ur1%z =	mt%x*mt%y*mt%z*f0%z + t%x*mt%y*mt%z*f1%z + mt%x*t%y*mt%z*f2%z + &
				t%x*t%y*mt%z*f3%z +	mt%x*mt%y*t%z*f4%z + t%x*mt%y*t%z*f5%z + &
				mt%x*t%y*t%z*f6%z + t%x*t%y*t%z*f7%z

			!cleanup
			deallocate(xg)
			deallocate(yg)
			deallocate(zg)
			call destroy_ur1(x0)
			call destroy_ur1(x1)
			call destroy_ur1(t)
			call destroy_ur1(mt)
			call destroy_ur1(f0)
			call destroy_ur1(f1)
			call destroy_ur1(f2)
			call destroy_ur1(f3)
			call destroy_ur1(f4)
			call destroy_ur1(f5)
			call destroy_ur1(f6)
			call destroy_ur1(f7)

		case(QUADRATIC,CUBIC)

			if(this_interpolator == QUADRATIC) INTERPOLATION_ORDER = 2
			if(this_interpolator == CUBIC) INTERPOLATION_ORDER = 3

			!-----
			!Polynomial splines of arbitrary order.
			!uses the bspline-fortran library.
			!Note, we pass INTERPOLATION_ORDER +1 to the library, corresponding to polynomial degre +1
			!-----
			allocate(xg(1-ng:ni+ng))
			allocate(yg(1-ng:nj+ng))
			allocate(zg(1-ng:nk+ng))
			xg(1-ng:ni+ng) =sgrid%grid%x(1-ng:ni+ng,1,1)
			yg(1-ng:nj+ng) =sgrid%grid%y(1,1-ng:nj+ng,1)
			zg(1-ng:nk+ng) =sgrid%grid%z(1,1,1-ng:nk+ng)
			!Special consideration for 2D:
			if(sgrid%gni == 1) then
				xg(0) = xg(1) - 1.0;
				xg(2) = xg(1) + 1.0;
			endif
			if(sgrid%gnj == 1) then
				yg(0) = yg(1) - 1.0;
				yg(2) = yg(1) + 1.0;
			endif
			if(sgrid%gnk == 1) then
				zg(0) = zg(1) - 1.0;
				zg(2) = zg(1) + 1.0;
			endif
			idx = 0; idy=0; idz=0; iflag = 0
			order_x = MAX(MIN(INTERPOLATION_ORDER+1,ni+2*ng-1),2)
			order_y = MAX(MIN(INTERPOLATION_ORDER+1,nj+2*ng-1),2)
			order_z = MAX(MIN(INTERPOLATION_ORDER+1,nk+2*ng-1),2)
			call bsp_x%initialize(xg,yg,zg,sr1%x,order_x,order_y,order_z,iflag)
			call bsp_y%initialize(xg,yg,zg,sr1%y,order_x,order_y,order_z,iflag)
			call bsp_z%initialize(xg,yg,zg,sr1%z,order_x,order_y,order_z,iflag)
			!interpolate
			do ip = 1,np
				call bsp_x%evaluate(lp%xp%x(ip),lp%xp%y(ip),lp%xp%z(ip),idx,idy,idz,ur1%x(ip),iflag)
				call bsp_y%evaluate(lp%xp%x(ip),lp%xp%y(ip),lp%xp%z(ip),idx,idy,idz,ur1%y(ip),iflag)
				call bsp_z%evaluate(lp%xp%x(ip),lp%xp%y(ip),lp%xp%z(ip),idx,idy,idz,ur1%z(ip),iflag)
			enddo
			!cleanup
			call bsp_x%destroy
			call bsp_y%destroy
			call bsp_z%destroy
			deallocate(xg)
			deallocate(yg)
			deallocate(zg)

		case(IDW)
			!-----
			!Inverse Distance Weighting
			!-----
			call interp_s2u_r1_idw(lp,sgrid%grid,ur1,sr1)

		case(TSE)
			!-----
			!Taylor series expansion about nearest node (second order)
			!Assumes that the ghost vale of sr1 are already updated
			!-----
			!allocate
			call init_sr2(gradsr1,ni,nj,nk,ng,'GradPhi')   !Need to change to accomodate non-rectilinear grids
			call init_ur2(gradsr1_p,np,'GradPhiAtP')
			call init_ur1(r1node,np,'UNODE')
			call init_ur1(delta,np,'DELTA')

			!compute grad(phi)
			if(sgrid%rectilinear) then
				call grad_sr1(sgrid,sr1,gradsr1)
			else
				call grad_sr1_ls(sgrid,sr1,gradsr1)
			endif

			!Gather
			do ip = 1,np
				r1node%x(ip) = sr1%x(no%x(ip),no%y(ip),no%z(ip))
				r1node%y(ip) = sr1%y(no%x(ip),no%y(ip),no%z(ip))
				r1node%z(ip) = sr1%z(no%x(ip),no%y(ip),no%z(ip))
				delta%x(ip) = sgrid%grid%x(no%x(ip),no%y(ip),no%z(ip))
				delta%y(ip) = sgrid%grid%y(no%x(ip),no%y(ip),no%z(ip))
				delta%z(ip) = sgrid%grid%z(no%x(ip),no%y(ip),no%z(ip))

				gradsr1_p%xx(ip) = gradsr1%xx(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%xy(ip) = gradsr1%xy(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%xz(ip) = gradsr1%xz(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%yx(ip) = gradsr1%yx(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%yy(ip) = gradsr1%yy(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%yz(ip) = gradsr1%yz(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%zx(ip) = gradsr1%zx(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%zy(ip) = gradsr1%zy(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%zz(ip) = gradsr1%zz(no%x(ip),no%y(ip),no%z(ip))
			enddo

			!Compute phi_p = phi_no + grad(phi)|_no . (x_p-x_no)
			delta%x(1:np) = lp%xp%x(1:np) - delta%x(1:np)
			delta%y(1:np) = lp%xp%y(1:np) - delta%y(1:np)
			delta%z(1:np) = lp%xp%z(1:np) - delta%z(1:np)
			ur1%x(1:np) = r1node%x(1:np) + delta%x(1:np)*gradsr1_p%xx(1:np) &
			+ delta%y(1:np)*gradsr1_p%xy(1:np)+ delta%z(1:np)*gradsr1_p%xz(1:np)
			ur1%y(1:np) = r1node%y(1:np) + delta%x(1:np)*gradsr1_p%yx(1:np) &
			+ delta%y(1:np)*gradsr1_p%yy(1:np)+ delta%z(1:np)*gradsr1_p%yz(1:np)
			ur1%z(1:np) = r1node%z(1:np) + delta%x(1:np)*gradsr1_p%zx(1:np) &
			+ delta%y(1:np)*gradsr1_p%zy(1:np)+ delta%z(1:np)*gradsr1_p%zz(1:np)

			call destroy_ur1(r1node)
			call destroy_ur1(delta)
			call destroy_ur2(gradsr1_p)
			call destroy_sr2(gradsr1)

		case(TSE_LIMIT)
			!-----
			!Taylor series expansion about nearest node (second order)
			!Assumes that the ghost vale of sr1 are already updated
			!Impose a limiter so that the interpolated value does not go out of bounds
			!-----
			!allocate
			call init_sr2(gradsr1,ni,nj,nk,ng,'GradPhi')   !Need to change to accomodate non-rectilinear grids
			call init_ur2(gradsr1_p,np,'GradPhiAtP')
			call init_ur1(r1node,np,'UNODE')
			call init_ur1(delta,np,'DELTA')
			call init_ur1(r1min,np,'MIN_LIMIT')
			call init_ur1(r1max,np,'MAX_LIMIT')

			!compute grad(phi)
			if(sgrid%rectilinear) then
				call grad_sr1(sgrid,sr1,gradsr1)
			else
				call grad_sr1_ls(sgrid,sr1,gradsr1)
			endif

			!Gather
			do ip = 1,np
				r1node%x(ip) = sr1%x(no%x(ip),no%y(ip),no%z(ip))
				r1node%y(ip) = sr1%y(no%x(ip),no%y(ip),no%z(ip))
				r1node%z(ip) = sr1%z(no%x(ip),no%y(ip),no%z(ip))
				delta%x(ip) = sgrid%grid%x(no%x(ip),no%y(ip),no%z(ip))
				delta%y(ip) = sgrid%grid%y(no%x(ip),no%y(ip),no%z(ip))
				delta%z(ip) = sgrid%grid%z(no%x(ip),no%y(ip),no%z(ip))

				gradsr1_p%xx(ip) = gradsr1%xx(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%xy(ip) = gradsr1%xy(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%xz(ip) = gradsr1%xz(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%yx(ip) = gradsr1%yx(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%yy(ip) = gradsr1%yy(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%yz(ip) = gradsr1%yz(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%zx(ip) = gradsr1%zx(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%zy(ip) = gradsr1%zy(no%x(ip),no%y(ip),no%z(ip))
				gradsr1_p%zz(ip) = gradsr1%zz(no%x(ip),no%y(ip),no%z(ip))

				!limiters:
				r1min%x(ip) = minval(sr1%x(no%x(ip)-1:no%x(ip)+1, no%y(ip)-1:no%y(ip)+1, no%z(ip)-1:no%z(ip)+1))
				r1min%y(ip) = minval(sr1%y(no%x(ip)-1:no%x(ip)+1, no%y(ip)-1:no%y(ip)+1, no%z(ip)-1:no%z(ip)+1))
				r1min%z(ip) = minval(sr1%z(no%x(ip)-1:no%x(ip)+1, no%y(ip)-1:no%y(ip)+1, no%z(ip)-1:no%z(ip)+1))
				r1max%x(ip) = maxval(sr1%x(no%x(ip)-1:no%x(ip)+1, no%y(ip)-1:no%y(ip)+1, no%z(ip)-1:no%z(ip)+1))
				r1max%y(ip) = maxval(sr1%y(no%x(ip)-1:no%x(ip)+1, no%y(ip)-1:no%y(ip)+1, no%z(ip)-1:no%z(ip)+1))
				r1max%z(ip) = maxval(sr1%z(no%x(ip)-1:no%x(ip)+1, no%y(ip)-1:no%y(ip)+1, no%z(ip)-1:no%z(ip)+1))
	
			enddo

			!Compute phi_p = phi_no + grad(phi)|_no . (x_p-x_no)
			delta%x(1:np) = lp%xp%x(1:np) - delta%x(1:np)
			delta%y(1:np) = lp%xp%y(1:np) - delta%y(1:np)
			delta%z(1:np) = lp%xp%z(1:np) - delta%z(1:np)
			ur1%x(1:np) = r1node%x(1:np) + delta%x(1:np)*gradsr1_p%xx(1:np) &
			+ delta%y(1:np)*gradsr1_p%xy(1:np)+ delta%z(1:np)*gradsr1_p%xz(1:np)
			ur1%y(1:np) = r1node%y(1:np) + delta%x(1:np)*gradsr1_p%yx(1:np) &
			+ delta%y(1:np)*gradsr1_p%yy(1:np)+ delta%z(1:np)*gradsr1_p%yz(1:np)
			ur1%z(1:np) = r1node%z(1:np) + delta%x(1:np)*gradsr1_p%zx(1:np) &
			+ delta%y(1:np)*gradsr1_p%zy(1:np)+ delta%z(1:np)*gradsr1_p%zz(1:np)

			!limit
			ur1%x(1:np) = max(ur1%x(1:np),r1min%x(1:np))
			ur1%y(1:np) = max(ur1%y(1:np),r1min%y(1:np))
			ur1%z(1:np) = max(ur1%z(1:np),r1min%z(1:np))
			ur1%x(1:np) = min(ur1%x(1:np),r1max%x(1:np))
			ur1%y(1:np) = min(ur1%y(1:np),r1max%y(1:np))
			ur1%z(1:np) = min(ur1%z(1:np),r1max%z(1:np))

			call destroy_ur1(r1node)
			call destroy_ur1(delta)
			call destroy_ur2(gradsr1_p)
			call destroy_sr2(gradsr1)
			call destroy_ur1(r1min)
			call destroy_ur1(r1max)

		case default
			write(*,*) 'ERROR:  Unrecognized interpolator:',this_interpolator
			CFD2LCS_ERROR = 1
			return
		end select

	end subroutine interp_s2u_r1

	subroutine interp_s2u_r1_idw(lp,grid,ur1,sr1)
		implicit none
		!-----
		type(lp_t):: lp
		type(sr1_t):: grid
		type(ur1_t):: ur1
		type(sr1_t):: sr1
		!-----
		!IDW: based on (possibly) 27 nearest nodes:
		!-----
		real(LCSRP):: tmp1(1:lp%np,1:N_NBR)
		real(LCSRP):: tmp2(1:lp%np,1:N_NBR)
		real(LCSRP):: tmp3(1:lp%np,1:N_NBR)
		real(LCSRP):: sumwts(1:lp%np)
		real(LCSRP):: fx(1:lp%np,1:N_NBR)
		real(LCSRP):: fy(1:lp%np,1:N_NBR)
		real(LCSRP):: fz(1:lp%np,1:N_NBR)
		integer:: i,j,k,ip,inbr,np
		!-----
		real(LCSRP),parameter:: SMALL = 100.0_LCSRP / huge(LCSRP)
		real(LCSRP),parameter:: P = 2.0_LCSRP

		if(lcsrank==0 )&
			write(*,*) 'in_interp_s2u_r1_idw... ERROR:  This is crap, why are you using this interp?'
		CFD2LCS_ERROR = 0
		return

		!brevity...
		np =lp%np

		!Gather all the data:
		do inbr = 1,N_NBR
			do ip = 1,np
				i = lp%no%x(ip)+NBR_OFFSET(1,inbr)
				j = lp%no%y(ip)+NBR_OFFSET(2,inbr)
				k = lp%no%z(ip)+NBR_OFFSET(3,inbr)
				tmp1(ip,inbr) = grid%x(i,j,k)
				tmp2(ip,inbr) = grid%y(i,j,k)
				tmp3(ip,inbr) = grid%z(i,j,k)
				fx(ip,inbr) = sr1%x(i,j,k)
				fy(ip,inbr) = sr1%y(i,j,k)
				fz(ip,inbr) = sr1%z(i,j,k)
			enddo
		enddo

		!compute distance (tmp1) and check if we are on a node
		do inbr = 1,N_NBR
				tmp1(:,inbr) = lp%xp%x(1:np) - tmp1(:,inbr) !dx
				tmp2(:,inbr) = lp%xp%y(1:np) - tmp2(:,inbr) !dy
				tmp3(:,inbr) = lp%xp%z(1:np) - tmp3(:,inbr) !dz
				!tmp1(:,inbr) = sqrt(tmp1(:,inbr)*tmp1(:,inbr) + tmp2(:,inbr)*tmp2(:,inbr) + tmp3(:,inbr)*tmp3(:,inbr)) !dist
				tmp1(:,inbr) = (tmp1(:,inbr)*tmp1(:,inbr) + tmp2(:,inbr)*tmp2(:,inbr) + tmp3(:,inbr)*tmp3(:,inbr)) !dist^2
				tmp1(:,inbr) = max(tmp1(:,inbr),SMALL)
		enddo

		!compute the weights (tmp2)
		do inbr = 1,N_NBR
			!tmp2(:,inbr) = 1.0_LCSRP/(tmp1(:,inbr))**P
			tmp2(:,inbr) = 1.0_LCSRP/(tmp1(:,inbr))  !For P=2
		enddo

		!Compute sumwts and then normalize
		sumwts = 0.0_LCSRP
		do inbr = 1,N_NBR
			sumwts(:) = sumwts(:) + tmp2(:,inbr)
		enddo
		do inbr = 1,N_NBR
			tmp2(:,inbr) = tmp2(:,inbr) / sumwts(:)
		enddo

		!Compute the values at the particles:
		ur1%x = 0.0_LCSRP
		ur1%y = 0.0_LCSRP
		ur1%z = 0.0_LCSRP
		do inbr = 1,N_NBR
			ur1%x(1:np) = ur1%x(1:np) + tmp2(:,inbr)*fx(:,inbr)
			ur1%y(1:np) = ur1%y(1:np) + tmp2(:,inbr)*fy(:,inbr)
			ur1%z(1:np) = ur1%z(1:np) + tmp2(:,inbr)*fz(:,inbr)
		enddo

	end subroutine interp_s2u_r1_idw

end module lp_tracking_m
!!DEBUG:
!if(i0 < 0 .or. i1 > ni+1 .or. j0<0 .or. j1>nj+1 .or. k0<0 .or. k1>nk+1) then
!write(*,*) '[',lcsrank,'] WTF: ip, nox,noy,noz',ip,no%x(ip),no%y(ip),no%z(ip)
!write(*,*) '[',lcsrank,'] WTF: ip, i0,j0,k0',ip,i0,j0,k0
!write(*,*) '[',lcsrank,'] WTF: ip, i1,j1,k1',ip,i1,j1,k1
!write(*,*) '[',lcsrank,'] WTF: ip, xp,yp,zp',ip,lp%xp%x(ip),lp%xp%y(ip),lp%xp%z(ip)
!write(*,*) '[',lcsrank,'] WTF: ip, xg,yg,zg',ip,xg(no%x(ip)),yg(no%y(ip)),zg(no%z(ip))
!write(*,*) '[',lcsrank,'] WTF: ip, node bcflag',ip,sgrid%bcflag%i(no%x(ip),no%y(ip),no%z(ip))
!endif

