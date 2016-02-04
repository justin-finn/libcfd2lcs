module lp_tracking_m
	use data_m
	implicit none

	!i,j,k Cartesian offsets
	integer,parameter:: N_NBR = 27
	integer(LCSIP),parameter:: NBR_OFFSET(3,N_NBR) = reshape((/&
		   0,	0,	 0, &  !self
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

	contains

	recursive subroutine track_lp2node(lp,sgrid)
		use comms_m
		implicit none
		!-----
		type(lp_t):: lp
		type(sgrid_t):: sgrid 
		!-----
		integer:: ip,ni,nj,nk,ng,ijk(3)
		real(LCSRP):: xp,yp,zp
		integer:: inbr,i,j,k,ierr
		real(LCSRP),allocatable:: rsq_min(:),rsq(:)
		real(LCSRP),allocatable:: xg(:),yg(:),zg(:)
		integer,allocatable:: ioff(:),joff(:),koff(:)
		integer:: offx, offy,offz
		logical:: recursive_check,my_recursive_check, my_do_recursion, do_recursion
		integer,save::  ntrack=0
		integer:: max_track
		!-----
		!Find the nearest i,j,k node for each LP
		!-----

		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in track_lp2no... ', trim(lp%label),' => ',trim(sgrid%label)

		!-----
		!Check if a recursive search is needed:
		!-----
		ntrack = ntrack+1
		if(any(lp%flag%i(1:lp%np)==LP_UNKNOWN))then
			my_recursive_check = .TRUE.
		else
			my_recursive_check = .FALSE.
		endif
		call MPI_ALLREDUCE(my_recursive_check,recursive_check,1,MPI_LOGICAL,MPI_LOR,lcscomm,ierr)

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
			!gather all the local grid points to a buffer:
			do ip = 1,lp%np
				i = lp%no%x(ip)+NBR_OFFSET(1,inbr)
				j = lp%no%y(ip)+NBR_OFFSET(2,inbr)
				k = lp%no%z(ip)+NBR_OFFSET(3,inbr)
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
		lp%no%x(1:lp%np) = lp%no%x(1:lp%np)+ioff(1:lp%np)
		lp%no%y(1:lp%np) = lp%no%y(1:lp%np)+joff(1:lp%np)
		lp%no%z(1:lp%np) = lp%no%z(1:lp%np)+koff(1:lp%np)

		!Some things we only do for FWD lp
		if(lp%direction == FWD) then
			!-----
			!Set lagrangian particle boundary condition
			!-----
			call set_lp_bc(lp,sgrid)

			!-----
			!Exchange particles that have crossed proc. boundaries
			!-----
			call exchange_lpdata(lp,sgrid)
		endif

		!-----
		!Check if we need to call again.
		!All processors must call together.
		!-----
		if(recursive_check) then
			!ensure we dont end up in an infinite tracking loop
			max_track = 2*max(sgrid%ni,sgrid%nj,sgrid%nk)
			if(ntrack>max_track) then
				write(*,*) 'ERROR: lcsrank[',lcsrank,'] cannot locate all particles after',ntrack,' tracking iterations.'
				CFD2LCS_ERROR = 1
				lp%flag%i(1:lp%np) = LP_IB
				return
			endif

			!If particle has not moved, we can stop.
			if(any(abs(ioff)>0) .OR. any(abs(joff)>0) .OR. any(abs(koff)>0))then
				my_do_recursion = .true.
			else
				my_do_recursion = .false.
			endif
			call MPI_ALLREDUCE(my_do_recursion,do_recursion,1,MPI_LOGICAL,MPI_LOR,lcscomm,ierr)

			if(do_recursion) then
				!cleanup memory and call again
				deallocate(rsq)
				deallocate(rsq_min)
				deallocate(xg)
				deallocate(yg)
				deallocate(zg)
				deallocate(ioff)
				deallocate(joff)
				deallocate(koff)
				call track_lp2node(lp,sgrid)
			else
				if(lcsrank==0) &
					write(*,*) ' Tracking all particles required ',ntrack,' iterations'
				lp%flag%i(1:lp%np) = LP_IB
				ntrack = 0
			endif
		else
			ntrack = 0
		endif

	end subroutine track_lp2node

	subroutine set_lp_bc(lp,sgrid)
		implicit none
		!-----
		type(lp_t):: lp
		type(sgrid_t):: sgrid
		!-----
		integer:: ip,ierr
		integer:: nbounce, nstick, my_nbounce, my_nstick
		!-----
		!Handle particles which encounter a non-periodic boundary
		!LCS_WALL, LCS_SLIP, LCS_INFLOW:  For now, just project the particle back to the IB node.
		! 	TODO: eventually, you will want a nice specular collision.
		!LCS_OUTFLOW:  Make the particle stick to the boundary TODO
		!-----

		my_nbounce = 0
		my_nstick = 0
		do ip = 1,lp%np
			!x0
			if(lp%no%x(ip) < 1 .AND. sgrid%bc_list(1)/= LCS_PERIODIC) then
				select case(sgrid%bc_list(1))
					case(LCS_WALL,LCS_SLIP,LCS_INFLOW,LCS_OUTFLOW)
						lp%no%x(ip) = 1
						lp%xp%x(ip) =  sgrid%grid%x(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
						my_nbounce = my_nbounce+1
						if(sgrid%bc_list(6)==LCS_OUTFLOW) then
							my_nstick = my_nstick+1
							lp%flag%i(ip) = LP_STICK
						endif
					case default
				end select
			endif
			!y0
			if(lp%no%y(ip) < 1 .AND. sgrid%bc_list(2)/= LCS_PERIODIC) then
				select case(sgrid%bc_list(2))
					case(LCS_WALL,LCS_SLIP,LCS_INFLOW,LCS_OUTFLOW)
						lp%no%y(ip) = 1
						lp%xp%y(ip) = sgrid%grid%y(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
						my_nbounce = my_nbounce+1
						if(sgrid%bc_list(6)==LCS_OUTFLOW) then
							my_nstick = my_nstick+1
							lp%flag%i(ip) = LP_STICK
						endif
					case default
				end select
			endif
			!z0
			if(lp%no%z(ip) < 1 .AND. sgrid%bc_list(3)/= LCS_PERIODIC) then
				select case(sgrid%bc_list(3))
					case(LCS_WALL,LCS_SLIP,LCS_INFLOW,LCS_OUTFLOW)
						lp%no%z(ip) = 1
						lp%xp%z(ip) = sgrid%grid%z(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
						my_nbounce = my_nbounce+1
						if(sgrid%bc_list(6)==LCS_OUTFLOW) then
							my_nstick = my_nstick+1
							lp%flag%i(ip) = LP_STICK
						endif
					case default
				end select
			endif

			!x1
			if(lp%no%x(ip) > sgrid%grid%ni .AND. sgrid%bc_list(4)/= LCS_PERIODIC) then
				select case(sgrid%bc_list(4))
					case(LCS_WALL,LCS_SLIP,LCS_INFLOW,LCS_OUTFLOW)
						lp%no%x(ip) = sgrid%grid%ni
						lp%xp%x(ip) = sgrid%grid%x(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
						my_nbounce = my_nbounce+1
						if(sgrid%bc_list(6)==LCS_OUTFLOW) then
							my_nstick = my_nstick+1
							lp%flag%i(ip) = LP_STICK
						endif
					case default
				end select
			endif
			!y1
			if(lp%no%y(ip) > sgrid%grid%nj .AND. sgrid%bc_list(5)/= LCS_PERIODIC) then
				select case(sgrid%bc_list(5))
					case(LCS_WALL,LCS_SLIP,LCS_INFLOW,LCS_OUTFLOW)
						lp%no%y(ip) = sgrid%grid%nj
						lp%xp%y(ip) = sgrid%grid%y(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
						my_nbounce = my_nbounce+1
						if(sgrid%bc_list(6)==LCS_OUTFLOW) then
							my_nstick = my_nstick+1
							lp%flag%i(ip) = LP_STICK
						endif
					case default
				end select
			endif
			!z1
			if(lp%no%z(ip) > sgrid%grid%nk .AND. sgrid%bc_list(6)/= LCS_PERIODIC) then
				select case(sgrid%bc_list(6))
					case(LCS_WALL,LCS_SLIP,LCS_INFLOW,LCS_OUTFLOW)
						lp%no%z(ip) = sgrid%grid%nk
						lp%xp%z(ip) = sgrid%grid%z(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
						my_nbounce = my_nbounce+1
						if(sgrid%bc_list(6)==LCS_OUTFLOW) then
							my_nstick = my_nstick+1
							lp%flag%i(ip) = LP_STICK
						endif
					case default
				end select
			endif
		enddo

		!count bc encounters:
		if(LCS_VERBOSE) then
			call MPI_REDUCE(my_nstick,nstick,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)
			if (lcsrank==0) write(*,*) 'handled',nstick,' stick bc'
			call MPI_REDUCE(my_nbounce,nbounce,1,MPI_INTEGER,MPI_SUM,0,lcscomm,ierr)
			if (lcsrank==0) write(*,*) 'handled',nbounce,' bounce bc'
		endif

	end subroutine set_lp_bc

	subroutine interp_s2u_r1(lp,grid,ur1,sr1)
		use bspline_oo_module
		use unstructured_m
		implicit none
		!-----
		type(lp_t):: lp
		type(sr1_t):: grid
		type(ur1_t):: ur1
		type(sr1_t):: sr1
		!-----
		!Interpolation using B-splines for 3D orthogonal grids
		type(bspline_3d):: bsp_x,bsp_y,bsp_z
		real(LCSRP),allocatable:: xg(:),yg(:),zg(:)
		integer:: iflag,idx,idy,idz
		integer:: order_x,order_y,order_z
		!Trilinear:
		integer:: ip
		integer:: i0,j0,k0,i1,j1,k1 
		type(ur1_t):: t,mt,x0,x1
		type(ur1_t):: f0,f1,f2,f3,f4,f5,f6,f7
		!-----
		!Interpolate structured data to unstructured pts
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in interp_s2u_r1... ',trim(sr1%label),' => ',trim(ur1%label)

		select case(INTERPOLATION_ORDER)
		case( : 0)
			!-----
			!Zeroth order, nearest node interp:
			!-----
			do ip = 1,lp%np
				ur1%x(ip) = sr1%x(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
				ur1%y(ip) = sr1%y(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
				ur1%z(ip) = sr1%z(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
			enddo

		case(1)
			!-----
			!First order, trilinear interpolation: (2nd order accurate)
			!-----
			allocate(xg(1-grid%ng:grid%ni+grid%ng))
			allocate(yg(1-grid%ng:grid%nj+grid%ng))
			allocate(zg(1-grid%ng:grid%nk+grid%ng))
			call init_ur1(x0,lp%np,'X0')
			call init_ur1(x1,lp%np,'X0')
			call init_ur1(t,lp%np,'T')
			call init_ur1(mt,lp%np,'1MINUST')
			call init_ur1(f0,lp%np,'F0')
			call init_ur1(f1,lp%np,'F1')
			call init_ur1(f2,lp%np,'F2')
			call init_ur1(f3,lp%np,'F3')
			call init_ur1(f4,lp%np,'F4')
			call init_ur1(f5,lp%np,'F5')
			call init_ur1(f6,lp%np,'F6')
			call init_ur1(f7,lp%np,'F7')
			
			xg(1-grid%ng:grid%ni+grid%ng) =grid%x(1-grid%ng:grid%ni+grid%ng,1,1)
			yg(1-grid%ng:grid%nj+grid%ng) =grid%y(1,1-grid%ng:grid%nj+grid%ng,1)
			zg(1-grid%ng:grid%nk+grid%ng) =grid%z(1,1,1-grid%ng:grid%nk+grid%ng)
			
			!Gather scattered data into vectors
			do ip = 1,lp%np
				if(lp%xp%x(ip) > xg(lp%no%x(ip))) then
					i0 = lp%no%x(ip); i1 = lp%no%x(ip)+1
				else
					i0 = lp%no%x(ip)-1; i1 = lp%no%x(ip)
				endif
				if(lp%xp%y(ip) > yg(lp%no%y(ip))) then
					j0 = lp%no%y(ip); j1 = lp%no%y(ip)+1
				else
					j0 = lp%no%y(ip)-1; j1 = lp%no%y(ip)
				endif
				if(lp%xp%z(ip) > zg(lp%no%z(ip))) then
					k0 = lp%no%z(ip); k1 = lp%no%z(ip)+1
				else
					k0 = lp%no%z(ip)-1; k1 = lp%no%z(ip)
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
			t%x(1:lp%np) = (lp%xp%x(1:lp%np) - x0%x(1:lp%np)) / (x1%x(1:lp%np)-x0%x(1:lp%np))
			t%y(1:lp%np) = (lp%xp%y(1:lp%np) - x0%y(1:lp%np)) / (x1%y(1:lp%np)-x0%y(1:lp%np))
			t%z(1:lp%np) = (lp%xp%z(1:lp%np) - x0%z(1:lp%np)) / (x1%z(1:lp%np)-x0%z(1:lp%np))
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

		case(2 : )
			!-----
			!Polynomial splines of arbitrary order.
			!uses the bspline-fortran library.
			!Note, we pass INTERPOLATION_ORDER +1 to the library, corresponding to polynomial degre +1
			!-----
			allocate(xg(1-grid%ng:grid%ni+grid%ng))
			allocate(yg(1-grid%ng:grid%nj+grid%ng))
			allocate(zg(1-grid%ng:grid%nk+grid%ng))
			xg(1-grid%ng:grid%ni+grid%ng) =grid%x(1-grid%ng:grid%ni+grid%ng,1,1)
			yg(1-grid%ng:grid%nj+grid%ng) =grid%y(1,1-grid%ng:grid%nj+grid%ng,1)
			zg(1-grid%ng:grid%nk+grid%ng) =grid%z(1,1,1-grid%ng:grid%nk+grid%ng)
			idx = 0; idy=0; idz=0; iflag = 0
			order_x = MAX(MIN(INTERPOLATION_ORDER+1,grid%ni+2*grid%ng-1),2)
			order_y = MAX(MIN(INTERPOLATION_ORDER+1,grid%nj+2*grid%ng-1),2)
			order_z = MAX(MIN(INTERPOLATION_ORDER+1,grid%nk+2*grid%ng-1),2)
			call bsp_x%initialize(xg,yg,zg,sr1%x,order_x,order_y,order_z,iflag)
			call bsp_y%initialize(xg,yg,zg,sr1%y,order_x,order_y,order_z,iflag)
			call bsp_z%initialize(xg,yg,zg,sr1%z,order_x,order_y,order_z,iflag)
			!interpolate
			do ip = 1,lp%np
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

		end select

	end subroutine interp_s2u_r1

end module lp_tracking_m
