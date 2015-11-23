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

		!Bspline order:
		integer,parameter:: BSP_ORDER = 2
		real(LCSRP),parameter:: IDW_EXPONENT = 2.0_LCSRP

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
		implicit none
		!-----
		type(lp_t):: lp
		type(sr1_t):: grid
		type(ur1_t):: ur1
		type(sr1_t):: sr1
		!-----
		!Interpolation using RBF
		integer:: ip,inbr,ierr
		integer,allocatable:: ig(:),jg(:),kg(:)
		real(LCSRP),allocatable::wts(:,:)
		!Interpolation using B-splines for 3D orthogonal grids
		type(bspline_3d):: bsp_x,bsp_y,bsp_z
		real(LCSRP),allocatable:: xb(:),yb(:),zb(:)
		integer:: iflag,idx,idy,idz
		integer:: order_x,order_y,order_z
		!-----
		!Interpolate structured data to unstructured pts
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in interp_s2u_r1... ',trim(sr1%label),' => ',trim(ur1%label)

		select case(INTERPOLATION)
		case(ZERO_ORDER)
			!-----
			!Simplest possible nearest node interp:
			!-----
			do ip = 1,lp%np
				ur1%x(ip) = sr1%x(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
				ur1%y(ip) = sr1%y(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
				ur1%z(ip) = sr1%z(lp%no%x(ip),lp%no%y(ip),lp%no%z(ip))
			enddo

		case(GAUSSIAN_RBF,IDW)
			!-----
			!Gaussian Radial Basis Function:
			!-----
			allocate(ig(1:lp%np))
			allocate(jg(1:lp%np))
			allocate(kg(1:lp%np))
			allocate(wts(1:lp%np,1:N_NBR))

			if(INTERPOLATION==GAUSSIAN_RBF)	call rbf_wts()
			if(INTERPOLATION==IDW)	call idw_wts()
			ur1%x = 0.0_LCSRP
			ur1%y = 0.0_LCSRP
			ur1%z = 0.0_LCSRP
			do inbr = 1,N_NBR
				ig(1:lp%np) = lp%no%x(1:lp%np) + NBR_OFFSET(1,inbr)
				jg(1:lp%np) = lp%no%y(1:lp%np) + NBR_OFFSET(2,inbr)
				kg(1:lp%np) = lp%no%z(1:lp%np) + NBR_OFFSET(3,inbr)
				do ip = 1,lp%np
					ur1%x(ip) = ur1%x(ip) + sr1%x(ig(ip),jg(ip),kg(ip))*wts(ip,inbr)
					ur1%y(ip) = ur1%y(ip) + sr1%y(ig(ip),jg(ip),kg(ip))*wts(ip,inbr)
					ur1%z(ip) = ur1%z(ip) + sr1%z(ig(ip),jg(ip),kg(ip))*wts(ip,inbr)
				enddo
			enddo

			deallocate(ig)
			deallocate(jg)
			deallocate(kg)
			deallocate(wts)

		case(BSPLINE)

			!Initialize
			allocate(xb(1-grid%ng:grid%ni+grid%ng))
			allocate(yb(1-grid%ng:grid%nj+grid%ng))
			allocate(zb(1-grid%ng:grid%nk+grid%ng))
			xb(1-grid%ng:grid%ni+grid%ng) =grid%x(1-grid%ng:grid%ni+grid%ng,1,1)
			yb(1-grid%ng:grid%nj+grid%ng) =grid%y(1,1-grid%ng:grid%nj+grid%ng,1)
			zb(1-grid%ng:grid%nk+grid%ng) =grid%z(1,1,1-grid%ng:grid%nk+grid%ng)
			idx = 0; idy=0; idz=0; iflag = 0
			order_x = MIN(BSP_ORDER,grid%ni+2*grid%ng-1)
			order_y = MIN(BSP_ORDER,grid%nj+2*grid%ng-1)
			order_z = MIN(BSP_ORDER,grid%nk+2*grid%ng-1)
			call bsp_x%initialize(xb,yb,zb,sr1%x,order_x,order_y,order_z,iflag)
			call bsp_y%initialize(xb,yb,zb,sr1%y,order_x,order_y,order_z,iflag)
			call bsp_z%initialize(xb,yb,zb,sr1%z,order_x,order_y,order_z,iflag)
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
			deallocate(xb)
			deallocate(yb)
			deallocate(zb)

		case default
			write(*,*) 'ERROR:  UNKNOWN INTERPOLATION',INTERPOLATION
			CFD2LCS_ERROR =1
			return
		end select

		contains
		subroutine rbf_wts()
			implicit none
			!-----
			real(LCSRP):: deltasq(1:lp%np),sumwts(1:lp%np)
			real(LCSRP):: xg(1:lp%np),yg(1:lp%np),zg(1:lp%np)
			integer:: i,j,k
			!-----
			!Compute weights for radial basis function interpolation
			!-----
			do inbr = 1,N_NBR
				!gather all the local grid points to a buffer:
				do ip = 1,lp%np
					i = lp%no%x(ip)+NBR_OFFSET(1,inbr)
					j = lp%no%y(ip)+NBR_OFFSET(2,inbr)
					k = lp%no%z(ip)+NBR_OFFSET(3,inbr)
					xg(ip) = grid%x(i,j,k)
					yg(ip) = grid%y(i,j,k)
					zg(ip) = grid%z(i,j,k)

					!Set characteristic length scale of the RBF= Delta^2
					if(inbr==1) then!self
						deltasq(ip) = (0.5_LCSRP*min(&
							grid%x(i+1,j,k)-grid%x(i-1,j,k),&
							grid%y(i,j+1,k)-grid%y(i,j-1,k),&
							grid%z(i,j,k+1)-grid%z(i,j,k-1)))**2
					endif
				enddo
				do ip = 1,lp%np
					wts(ip,inbr) = &
						+(xg(ip)-lp%xp%x(ip))*(xg(ip)-lp%xp%x(ip)) &
						+(yg(ip)-lp%xp%y(ip))*(yg(ip)-lp%xp%y(ip)) &
						+(zg(ip)-lp%xp%z(ip))*(zg(ip)-lp%xp%z(ip))  !save r^2 in wts
				enddo


			enddo

			!-----
			!Compute, normalize and save the interpolation weights
			!-----
			sumwts = 0.0_LCSRP
			do inbr = 1,N_NBR
				wts(1:lp%np,inbr) = exp(-wts(1:lp%np,inbr)/deltasq(1:lp%np))
				sumwts(1:lp%np) = sumwts(1:lp%np) + wts(1:lp%np,inbr)
			enddo
			do inbr = 1,N_NBR
				wts(1:lp%np,inbr) = wts(1:lp%np,inbr)/sumwts(1:lp%np)
			enddo
		end subroutine rbf_wts

		subroutine idw_wts()
			implicit none
			!-----
			real(LCSRP):: sumwts(1:lp%np)
			real(LCSRP):: xg(1:lp%np),yg(1:lp%np),zg(1:lp%np)
			integer:: i,j,k
			real(LCSRP),parameter:: P = IDW_EXPONENT*0.5_LCSRP
			real(LCSRP),parameter:: SMALL = 1e-10
			!-----
			!Compute weights for radial basis function interpolation
			!-----
			do inbr = 1,N_NBR
				!gather all the local grid points to a buffer:
				do ip = 1,lp%np
					i = lp%no%x(ip)+NBR_OFFSET(1,inbr)
					j = lp%no%y(ip)+NBR_OFFSET(2,inbr)
					k = lp%no%z(ip)+NBR_OFFSET(3,inbr)
					xg(ip) = grid%x(i,j,k)
					yg(ip) = grid%y(i,j,k)
					zg(ip) = grid%z(i,j,k)
				enddo
				do ip = 1,lp%np
					wts(ip,inbr) = 1.0_LCSRP/(&
						+(xg(ip)-lp%xp%x(ip))*(xg(ip)-lp%xp%x(ip)) &
						+(yg(ip)-lp%xp%y(ip))*(yg(ip)-lp%xp%y(ip)) &
						+(zg(ip)-lp%xp%z(ip))*(zg(ip)-lp%xp%z(ip))+SMALL)**P
				enddo
			enddo

			!-----
			!Normalize
			!-----
			sumwts = 0.0_LCSRP
			do inbr = 1,N_NBR
				sumwts(1:lp%np) = sumwts(1:lp%np) + wts(1:lp%np,inbr)
			enddo
			do inbr = 1,N_NBR
				wts(1:lp%np,inbr) = wts(1:lp%np,inbr)/sumwts(1:lp%np)
			enddo

		end subroutine idw_wts

	end subroutine interp_s2u_r1

end module lp_tracking_m
