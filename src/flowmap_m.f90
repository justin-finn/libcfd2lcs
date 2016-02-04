module flowmap_m
	use data_m
	use io_m
	use structured_m
	implicit none

	contains

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

		if(lcsrank==0) &
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

		if(lcsrank==0)&
			write(*,*) ' Time h flowmap saved to: ', trim(fname)

	end subroutine write_flowmap_substep

	subroutine reconstruct_flowmap(lcs)
		use lp_m
		use comms_m
		use lp_tracking_m
		implicit none
		type(lcs_t):: lcs
		!-----
		integer:: i, step, nstep, inc
		integer:: gn(3), offset(3)
		character(len=128):: fname
		type(sr1_t):: fm_tmp
		logical:: file_exists
		integer:: ni,nj,nk,ng
		integer:: ierr
		integer:: ip,np
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
			if(lcsrank==0)&
				write(*,*) i,'about to open file: ', trim(fname)
			call structured_io(trim(fname),IO_READ,gn,offset,r1=lcs%fm)	!read the flowmap

			!Set BC:
			call exchange_sdata(lcs%sgrid%scomm_face_r1,r1=lcs%fm)
			call set_flowmap_bc(lcs%fm,lcs%sgrid)
			call set_lp_bc(lcs%lp,lcs%sgrid)

			!Interp flow map to particles (store in lp%up)
			call interp_s2u_r1(lcs%lp,lcs%sgrid%grid,lcs%lp%up,lcs%fm)

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
			call track_lp2node(lcs%lp,lcs%sgrid)

			!next step
			step = step + inc
		enddo

		!map back to the origin:
		call exchange_lpmap(lcs%lp,lcs%fm)

		!cleanup
		call destroy_sr1(fm_tmp)
		if (inc == -1) lcs%lp%direction = BKWD  !reset direction 

	end subroutine reconstruct_flowmap

	subroutine set_flowmap_bc(fm, sgrid)
		implicit none
		!-----
		type(sr1_t):: fm
		type(sgrid_t):: sgrid
		!-----
		integer:: i,j,k,ni,nj,nk,ng
		!-----

		if(lcsrank==0) &
			write(*,*) 'In set_flowmap_bc... '

		ni = fm%ni
		nj = fm%nj
		nk = fm%nk
		ng = fm%ng

		!
		! X0
		!
		select case(sgrid%bc_list(1))
			case(LCS_OUTFLOW) !No slip...
				do i = 0,1-ng
					fm%x(i,:,:) = 0.0
					fm%y(i,:,:) = 0.0
					fm%z(i,:,:) = 0.0
				enddo
			case(LCS_INFLOW,LCS_WALL,LCS_SLIP)
				do i = 0,1-ng
					fm%x(i,:,:) = fm%x(1,:,:)
					fm%y(i,:,:) = fm%y(1,:,:)
					fm%z(i,:,:) = fm%z(1,:,:)
				enddo
			case default
				!do nothing...
		end select
		!
		! X1
		!
		select case(sgrid%bc_list(4))
			case(LCS_OUTFLOW) !No slip...
				do i = ni+1,ni+ng
					fm%x(i,:,:) = 0.0
					fm%y(i,:,:) = 0.0
					fm%z(i,:,:) = 0.0
				enddo
			case(LCS_INFLOW,LCS_WALL,LCS_SLIP) !zero gradient
				do i = ni+1,ni+ng
					fm%x(i,:,:) = fm%x(ni,:,:)
					fm%y(i,:,:) = fm%y(ni,:,:)
					fm%z(i,:,:) = fm%z(ni,:,:)
				enddo
			case default
				!do nothing...
		end select
		!
		! Y0
		!
		select case(sgrid%bc_list(2))
			case(LCS_OUTFLOW) !No slip...
				do j = 0,1-ng
					fm%x(:,j,:) = 0.0
					fm%y(:,j,:) = 0.0
					fm%z(:,j,:) = 0.0
				enddo
			case(LCS_INFLOW,LCS_WALL,LCS_SLIP)
				do j = 0,1-ng
					fm%x(:,j,:) = fm%x(:,1,:)
					fm%y(:,j,:) = fm%y(:,1,:)
					fm%z(:,j,:) = fm%z(:,1,:)
				enddo
			case default
				!do nothing...
		end select
		!
		! Y1
		!
		select case(sgrid%bc_list(5))
			case(LCS_OUTFLOW) !No slip...
				do j = nj+1,nj+ng
					fm%x(:,j,:) = 0.0
					fm%y(:,j,:) = 0.0
					fm%z(:,j,:) = 0.0
				enddo
			case(LCS_INFLOW,LCS_WALL,LCS_SLIP) !zero gradient
				do j = nj+1,nj+ng
					fm%x(:,j,:) = fm%x(:,nj,:)
					fm%y(:,j,:) = fm%y(:,nj,:)
					fm%z(:,j,:) = fm%z(:,nj,:)
				enddo
			case default
				!do nothing...
		end select
		!
		! Z0
		!
		select case(sgrid%bc_list(3))
			case(LCS_OUTFLOW) !No slip...
				do k = 0,1-ng
					fm%x(:,:,k) = 0.0
					fm%y(:,:,k) = 0.0
					fm%z(:,:,k) = 0.0
				enddo
			case(LCS_INFLOW,LCS_WALL,LCS_SLIP)
				do k = 0,1-ng
					fm%x(:,:,k) = fm%x(:,:,1)
					fm%y(:,:,k) = fm%y(:,:,1)
					fm%z(:,:,k) = fm%z(:,:,1)
				enddo
			case default
				!do nothing...
		end select
		!
		! Z1
		!
		select case(sgrid%bc_list(6))
			case(LCS_OUTFLOW) !No slip...
				do k = nk+1,nk+ng
					fm%x(:,:,k) = 0.0
					fm%y(:,:,k) = 0.0
					fm%z(:,:,k) = 0.0
				enddo
			case(LCS_INFLOW,LCS_WALL,LCS_SLIP) !zero gradient
				do k = nk+1,nk+ng
					fm%x(:,:,k) = fm%x(:,:,nk)
					fm%y(:,:,k) = fm%y(:,:,nk)
					fm%z(:,:,k) = fm%z(:,:,nk)
				enddo
			case default
				!do nothing...
		end select

	end subroutine set_flowmap_bc

end module flowmap_m
