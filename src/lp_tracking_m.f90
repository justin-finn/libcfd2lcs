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

	subroutine track_lp2node(lp,sgrid)
		use comms_m
		implicit none
		!-----
		type(lp_t):: lp
		type(sgrid_t):: sgrid
		!-----
		integer:: ip,ni,nj,nk,ng,np,ijk(3)
		real(LCSRP):: xp,yp,zp
		integer:: inbr,i,j,k
		real(LCSRP):: rsq3d(1:sgrid%grid%ni,1:sgrid%grid%nj,1:sgrid%grid%nk)
		real(LCSRP):: rsq_min(1:lp%np),rsq(1:lp%np)
		real(LCSRP):: xg(1:lp%np),yg(1:lp%np),zg(1:lp%np)
		integer:: ioff(1:lp%np),joff(1:lp%np),koff(1:lp%np)
		integer:: offx, offy,offz
		!-----
		!Find the nearest i,j,k node for each LP
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in track_lp2no... ', trim(lp%label),' => ',trim(sgrid%label)

		!Brevity...
		ni = sgrid%grid%ni
		nj = sgrid%grid%nj
		nk = sgrid%grid%nk
		ng = sgrid%grid%ng

		!-----
		!Use a brute force loop if needed to
		!locate all particles to the nearest in-bound node
		!try to vectorize over grid points
		!-----
		do ip = 1,lp%np
			if(lp%flag%i(ip)==LP_UNKNOWN) then
				xp = lp%xp%x(ip)
				yp = lp%xp%y(ip)
				zp = lp%xp%z(ip)
				rsq3d(1:ni,1:nj,1:nk) = &
					(sgrid%grid%x(1:ni,1:nj,1:nk)-xp)*(sgrid%grid%x(1:ni,1:nj,1:nk)-xp) &
				+	(sgrid%grid%y(1:ni,1:nj,1:nk)-yp)*(sgrid%grid%y(1:ni,1:nj,1:nk)-yp) &
				+	(sgrid%grid%z(1:ni,1:nj,1:nk)-zp)*(sgrid%grid%z(1:ni,1:nj,1:nk)-zp)
				ijk = minloc(rsq3d)
				lp%no%x(ip) = ijk(1)
				lp%no%y(ip) = ijk(2)
				lp%no%z(ip) = ijk(3)
				lp%flag%i(ip) = LP_IB
			endif
		enddo

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

		!Set lagrangian particle boundary condition
		call set_lp_bc(lp,sgrid)
		
		!Exchange particles that have crossed proc. boundaries
		call exchange_lpdata(lp,sgrid)

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
		implicit none
		!-----
		type(lp_t):: lp
		type(sr1_t):: grid
		type(ur1_t):: ur1
		type(sr1_t):: sr1
		!-----
		integer:: ip,inbr,ierr
		integer:: ig(1:lp%np),jg(1:lp%np),kg(1:lp%np)
		real(LCSRP)::wts(1:lp%np,1:N_NBR)
		!-----
		!Interpolate structured data to unstructured pts
		!-----
		if(lcsrank==0 .AND. LCS_VERBOSE)&
			write(*,*) 'in interp_sgrid2lp... ',trim(sr1%label),' => ',trim(ur1%label)

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
		case(GAUSSIAN_RBF)
			!-----
			!Gaussian Radial Basis Function:
			!-----
			call rbf_wts()
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
		case(TRICUBIC)
			!TODO:
			write(*,*) 'ERROR:  TRICUBIC NOT YET IMPLEMENTED'
			CFD2LCS_ERROR =1
			return
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
			deltasq = 0.0_LCSRP
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
					wts(ip,inbr) = &
						+(xg(ip)-lp%xp%x(ip))*(xg(ip)-lp%xp%x(ip)) &
						+(yg(ip)-lp%xp%y(ip))*(yg(ip)-lp%xp%y(ip)) &
						+(zg(ip)-lp%xp%z(ip))*(zg(ip)-lp%xp%z(ip))  !save r^2 in wts
				enddo
				deltasq(1:lp%np) = deltasq(1:lp%np) + wts(1:lp%np,inbr)
			enddo
			
			!-----
			!Set the length scale based on the local grid
			!-----
			deltasq = (1.0_LCSRP/3.0_LCSRP) * deltasq / real(N_NBR,LCSRP)  !this is the interp length scale
			
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
	end subroutine interp_s2u_r1

end module lp_tracking_m
