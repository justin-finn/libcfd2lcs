!Top level, user interface module.
subroutine cfd2lcs_init(cfdcomm,n,offset,x,y,z,flag)
	use io_m
	use sgrid_m
	implicit none
	!----
	integer(LCSIP):: cfdcomm
	integer(LCSIP):: n(3),offset(3)
	real(LCSRP):: x(1:n(1),1:n(2),1:n(3))
	real(LCSRP):: y(1:n(1),1:n(2),1:n(3))
	real(LCSRP):: z(1:n(1),1:n(2),1:n(3))
	integer(LCSIP):: flag(1:n(1),1:n(2),1:n(3))
	!----
	integer:: error,success1,success2,ierr
	!----
	if(CFD2LCS_ERROR /= 0) return

	!Initialize data counters
	NLCS = 0
	NLP = 0
	NSGRID = 0
	integrations_fwd = 0
	integrations_bkwd = 0
	integrations_fwd_c = 0
	integrations_bkwd_c = 0

	!init the mpi
	call init_lcs_mpi(cfdcomm)

	if(lcsrank ==0)&
		write(*,*) 'in cfd2lcs_init...'

	!Make sure the required output and tmp directories exist
	if(lcsrank ==0) then
		call system("mkdir -p ./"//trim(OUTPUT_DIR), success1)
		call system("mkdir -p ./"//trim(TEMP_DIR), success2)
		if (success1 /= 0 .OR. success2 /=0) then
			write(*,*) 'ERROR:  cfd2lcs cannot create required output directories'
			CFD2LCS_ERROR = 1
		endif
	endif

	!Init the default structured cfd storage (scfd) :
	scfd%label = 'CFD_DATA'
	call init_sgrid(scfd%sgrid,'CFD_GRID',n,offset,x,y,z,flag)
	call init_sr1(scfd%u_n,scfd%sgrid%ni,scfd%sgrid%nj,scfd%sgrid%nk,scfd%sgrid%ng,'U_N',translate=.false.)
	call init_sr1(scfd%u_np1,scfd%sgrid%ni,scfd%sgrid%nj,scfd%sgrid%nk,scfd%sgrid%ng,'U_NP1',translate=.false.)
	call compute_delta_xyz(scfd%sgrid,scfd%delta)

	!Check:
	call cfd2lcs_error_check(error)

end subroutine cfd2lcs_init

subroutine cfd2lcs_update(n,ux,uy,uz,time)
	use data_m
	use io_m
	use comms_m
	use sgrid_m
	use lp_motion_m
	use lcs_m
	use lp_m
	use flowmap_m
	implicit none
	!----
	integer:: n(3)
	real(LCSRP):: ux(1:n(1),1:n(2),1:n(3))
	real(LCSRP):: uy(1:n(1),1:n(2),1:n(3))
	real(LCSRP):: uz(1:n(1),1:n(2),1:n(3))
	real(LCSRP), intent(in):: time
	!----
	integer:: error,ierr
	integer:: ilp,ilcs
	type(lp_t),pointer:: lp
	type(lcs_t),pointer:: lcs
	logical,save:: FIRST_CALL = .true.
	real:: t2,t3,t0,t1
	logical:: fm_complete
	!----
	if(CFD2LCS_ERROR /= 0) return

	t_start_update = cputimer(lcscomm,SYNC_TIMER)

	if(lcsrank ==0)then
		write(*,'(a)') 	'------------------------------------------------------------------'
		write(*,'(a,ES11.4,a)')	'-----libcfd2lcs update: T=',time,'-----------------------------'
	endif

	!-----
	!Check we got an arrays of the correct size:
	!-----
	if	( scfd%sgrid%ni/=n(1) .OR. scfd%sgrid%nj /= n(2) .OR. scfd%sgrid%nk /=n(3)) then
		write(*,'(a,i6,a)') 'rank[',lcsrank,'] received velocity array of incorrect dimension'
		write(*,'(a,i6,a,i4,i4,i4,a)') 'rank[',lcsrank,'] [ni,nj,nk]= [',n(1),n(2),n(3),']'
		write(*,'(a,i6,a,i4,i4,i4,a)') 'rank[',lcsrank,'] sgrid[ni,nj,nk]= [',scfd%sgrid%ni,scfd%sgrid%nj,scfd%sgrid%nk,']'
		CFD2LCS_ERROR = 1
		return
	endif

	!-----
	!Set the new velocity, update ghosts and fakes:
	!-----
	if(FIRST_CALL) then
		scfd%t_n = time  !allows us to not start at t=0
		scfd%t_np1 = time  !allows us to not start at t=0
		FIRST_CALL = .FALSE.
	else
		scfd%t_n = scfd%t_np1
		scfd%t_np1 = time
	endif
	scfd%u_n = scfd%u_np1  !Shift down the velocity field from np1 => n
	scfd%u_np1%x(1:n(1),1:n(2),1:n(3)) = ux(1:n(1),1:n(2),1:n(3))
	scfd%u_np1%y(1:n(1),1:n(2),1:n(3)) = uy(1:n(1),1:n(2),1:n(3))
	scfd%u_np1%z(1:n(1),1:n(2),1:n(3)) = uz(1:n(1),1:n(2),1:n(3))
	call exchange_sdata(scfd%sgrid%scomm_max_r1,r1=scfd%u_np1)
	call set_velocity_bc(scfd%sgrid,scfd%u_np1)

	!-----
	! Update each LP set:
	!-----
	do ilp = 1, NLP
		lp => lp_c(ilp)
		if(lp%direction == FWD) then
			t2 = cputimer(lcscomm,SYNC_TIMER)
			call update_lp(lp,scfd)
			t3 = cputimer(lcscomm,SYNC_TIMER)
			cpu_fwd = cpu_fwd + max(t3-t2,0.0)
			this_cpu_fwd = this_cpu_fwd + max(t3-t2,0.0)
		endif
		if(lp%direction == BKWD) then
			t2 = cputimer(lcscomm,SYNC_TIMER)
			call update_flowmap_sl(lp,lp%sgrid,lp%fm,scfd)
			t3 = cputimer(lcscomm,SYNC_TIMER)
			cpu_bkwd = cpu_bkwd + max(t3-t2,0.0)
			this_cpu_bkwd = this_cpu_bkwd + max(t3-t2,0.0)
		endif
	enddo

	!-----
	! Update each lcs diagnostic:
	! If this timestep corresponds to a flowmap substep inerval, then compute/output the diagnostic
	!-----
	do ilcs = 1, NLCS
		lcs => lcs_c(ilcs)
		if(lcs%id<0) cycle !only for active lcs
		if( int(scfd%t_np1/lcs%h) /= int(scfd%t_n/lcs%h) ) then

			select case(lcs%diagnostic)
				case(FTLE_FWD,FTLE_BKWD)
					!-----
					!Map the forward time particles back to their original grid:
					!-----
					if(lcs%diagnostic==FTLE_FWD) then
						call exchange_lpmap(lcs%lp)
						if(AUX_GRID) then
							call exchange_lpmap(lcs%lpX0)
							call exchange_lpmap(lcs%lpY0)
							call exchange_lpmap(lcs%lpZ0)
							call exchange_lpmap(lcs%lpX1)
							call exchange_lpmap(lcs%lpY1)
							call exchange_lpmap(lcs%lpZ1)
						endif
					end if
					!-----
					!Write temp files for the scfd map substep
					!-----
					call write_flowmap_substep(lcs%lp,lcs%T,lcs%h)
					if(AUX_GRID) then
						call write_flowmap_substep(lcs%lpX0,lcs%T,lcs%h)
						call write_flowmap_substep(lcs%lpY0,lcs%T,lcs%h)
						call write_flowmap_substep(lcs%lpZ0,lcs%T,lcs%h)
						call write_flowmap_substep(lcs%lpX1,lcs%T,lcs%h)
						call write_flowmap_substep(lcs%lpY1,lcs%T,lcs%h)
						call write_flowmap_substep(lcs%lpZ1,lcs%T,lcs%h)
					endif
					!-----
					!Reconstruct the time T flowmap from time h substeps
					!-----
					call reconstruct_flowmap(lcs%lp,lcs%T,lcs%h,fm_complete)
					if(AUX_GRID) then
						call reconstruct_flowmap(lcs%lpX0,lcs%T,lcs%h,fm_complete)
						call reconstruct_flowmap(lcs%lpY0,lcs%T,lcs%h,fm_complete)
						call reconstruct_flowmap(lcs%lpZ0,lcs%T,lcs%h,fm_complete)
						call reconstruct_flowmap(lcs%lpX1,lcs%T,lcs%h,fm_complete)
						call reconstruct_flowmap(lcs%lpY1,lcs%T,lcs%h,fm_complete)
						call reconstruct_flowmap(lcs%lpZ1,lcs%T,lcs%h,fm_complete)
					endif
					!-----
					!Compute the FTLE
					!-----
					call compute_ftle(lcs)
					!-----
					!Write the time T LCS
					!-----
					!if(fm_complete) then
					call write_lcs(lcs,scfd%t_np1)
					!endif
					!-----
					!Reset the scfd maps for FTLE type diagnostics
					!Need to track particles to scfd grid when you do:
					!-----
					call reset_lp(lcs%lp)
					call track_lp2node(lcs%lp,scfd%sgrid,lcs%lp%no_scfd)
					if(AUX_GRID) then
						call reset_lp(lcs%lpX0)
						call reset_lp(lcs%lpY0)
						call reset_lp(lcs%lpZ0)
						call reset_lp(lcs%lpX1)
						call reset_lp(lcs%lpY1)
						call reset_lp(lcs%lpZ1)
						call track_lp2node(lcs%lpX0,scfd%sgrid,lcs%lpX0%no_scfd)
						call track_lp2node(lcs%lpY0,scfd%sgrid,lcs%lpY0%no_scfd)
						call track_lp2node(lcs%lpZ0,scfd%sgrid,lcs%lpZ0%no_scfd)
						call track_lp2node(lcs%lpX1,scfd%sgrid,lcs%lpX1%no_scfd)
						call track_lp2node(lcs%lpY1,scfd%sgrid,lcs%lpY1%no_scfd)
						call track_lp2node(lcs%lpZ1,scfd%sgrid,lcs%lpZ1%no_scfd)
					endif

				case(LP_TRACER)
					call write_lcs(lcs,scfd%t_np1)
					!Kill the tracers at the end of their lifetime (T):
					if( lcs%lp%lifetime > lcs%T)  then
						call destroy_lp(lcs%lp)
						call destroy_lcs(lcs)
					endif

				case default

			end select
		endif
	enddo

	!-----
	!Check
	!-----
	call cfd2lcs_error_check(error)

	!-----
	!Spit out some info:
	!-----
	t_finish_update = cputimer(lcscomm,SYNC_TIMER)
	call cfd2lcs_info()

end subroutine cfd2lcs_update

subroutine cfd2lcs_diagnostic_init(lcs_handle,lcs_type,resolution,T,h,label)
	use data_m
	use sgrid_m
	use lp_m
	use lp_tracking_m
	use lcs_m
	implicit none
	!----
	integer(LCSIP),intent(out):: lcs_handle
	integer(LCSIP),intent(in):: lcs_type
	integer(LCSIP),intent(in):: resolution
	real(LCSRP),intent(in):: T
	real(LCSRP),intent(in):: h
	character(len=*),intent(in):: label
	!----
	type(lcs_t),pointer:: lcs
	integer:: error
	integer:: ilp,ilcs
	!----
	!Initialize an lcs diagnostic.
	!----
	if(CFD2LCS_ERROR /= 0) return

	if(lcsrank ==0)&
		write(*,*) 'in cfd2lcs_diagnostic_init... ',trim(label)

	!----
	!Point to the next available lcs and set this up
	!----
	NLCS = NLCS + 1
	if(NLCS > NMAX_STRUCT) then
		if(lcsrank==0)then
			write(*,*) 'ERROR: NLCS exceeds NMAX_STRUCT.',NLCS,NMAX_STRUCT
			write(*,*) 'Increase NMAX_STRUCT in data_m.f90 and recompile libcfd2lcs'
		endif
		CFD2LCS_ERROR = 1
		return
	endif
	lcs => lcs_c(NLCS)
	lcs%id = NLCS
	lcs%diagnostic = lcs_type
	lcs%label = trim(label)
	lcs%T = T
	lcs%h = h

	!-----
	!Define the grid for this LCS diagnostic
	!resolution = 0 :	We are using the CFD grid for the LCS calculations
	!resolution = 1 :	Add/Remove gride points from existing CFD grid
	!-----
	if(resolution == 0) then
		if(lcsrank ==0)&
			write(*,*) 'Using Native CFD grid'
		lcs%sgrid => scfd%sgrid
	else
		if(lcsrank ==0)&
			write(*,*) 'New Grid with resolution factor',resolution
		call new_sgrid_from_sgrid(lcs%sgrid,scfd%sgrid,trim(lcs%label)//'-grid',resolution)
	endif

	!-----
	!Figure out what we are dealing with and initialize appropriately:
	!-----
	select case(lcs%diagnostic)

		case(FTLE_FWD)
			if(lcsrank ==0)&
				write(*,*) 'FWD Time FTLE:  Name: ',(lcs%label)
			call init_ftle(lcs,label,FWD,resolution)

		case(FTLE_BKWD)
			if(lcsrank ==0)&
				write(*,*) 'BKWD Time FTLE:  Name: ',(lcs%label)
			call init_ftle(lcs,label,BKWD,resolution)

		case(LP_TRACER)
			if(lcsrank ==0)&
				write(*,*) 'Lagrangian Particle Tracers:  Name: ',(lcs%label)
			call init_lp_tracer(lcs%lp,lcs%sgrid,FWD,1.0_LCSRP,trim(label)//'-LP')
			!Check that some tracers actually got injected:
			if (lcs%lp%npall ==0) then
				call destroy_lp(lcs%lp)
				call destroy_lcs(lcs)
				lcs_handle = -1
				return
			endif
			call track_lp2node(lcs%lp,scfd%sgrid,lcs%lp%no_scfd) !Track lp to the cfd grid
		
		case default
			if(lcsrank ==0)&
				write(*,'(a)') 'ERROR, bad specification for lcs_type.&
					& Options are: FTLE_FWD, FTLE_BKWD, LP_TRACER'
			CFD2LCS_ERROR = 1

	end select

	!-----
	!Pass back the id in lcs_handle
	!-----
	lcs_handle = lcs%id

	!-----
	!Check
	!-----
	call cfd2lcs_error_check(error)

end subroutine cfd2lcs_diagnostic_init

subroutine cfd2lcs_diagnostic_destroy(lcs_handle)
	use data_m
	implicit none
	!-----
	integer:: lcs_handle
	!-----
	type(lcs_t),pointer:: lcs
	integer:: ilcs
	integer:: error
	!-----
	if(CFD2LCS_ERROR /= 0) return

	do ilcs = 1, NLCS
		lcs => lcs_c(ilcs)
		if(lcs%id==lcs_handle) then
			if(lcsrank==0) &
				write(*,*) 'Destroying LCS: ', trim(lcs%label)
				!TODO:
		endif
	enddo
	
	!-----
	!Check
	!-----
	call cfd2lcs_error_check(error)

end subroutine cfd2lcs_diagnostic_destroy

subroutine cfd2lcs_finalize()
	use data_m
	use sgrid_m
	use lcs_m
	use lp_m
	implicit none
	!-----
	integer:: idata
	type(sgrid_t),pointer:: sgrid
	type(lp_t),pointer:: lp
	type(lcs_t),pointer:: lcs
	integer:: error
	!-----
	if(CFD2LCS_ERROR /= 0) return

	if(lcsrank ==0)&
		write(*,'(a)') 'in cfd2lcs_finalize...'

	!Deallocate all LCS
	do idata = 1,NLCS
		lcs => lcs_c(idata)
		call destroy_lcs(lcs)
	enddo

	!Deallocate all lp
	do idata = 1,NLP
		lp => lp_c(idata)
		call destroy_lp(lp)
	enddo

	!Deallocate all sgrid
	do idata = 1,NSGRID
		sgrid => sgrid_c(idata)
		call destroy_sgrid(sgrid)
	enddo

	!Cleanup scfd
	call destroy_sr1(scfd%u_n)
	call destroy_sr1(scfd%u_np1)
	scfd%label = 'Unused CFD data'
	
	!-----
	!Check
	!-----
	call cfd2lcs_error_check(error)

end subroutine cfd2lcs_finalize

subroutine cfd2lcs_set_option(option,val)
	use data_m
	implicit none
	!-----
	character(len=*):: option
	integer:: val
	!-----
	character(len=32):: str
	logical:: warn = .false.
	integer:: error
	!-----
	!Parse any options passed by the user
	!-----
	if(CFD2LCS_ERROR /= 0) return

	if(lcsrank ==0 .AND. LCS_VERBOSE)&
		write(*,'(a)') 'in cfd2lcs_set_option... '

	select case(trim(option))
		!Sync timer:
		case("SYNCTIMER")
			select case(val)
				case(LCS_TRUE)
					SYNC_TIMER = .TRUE.
					str= 'TRUE'
				case(LCS_FALSE)
					SYNC_TIMER = .FALSE.
					str= 'FALSE'
				case default
					warn = .true.
			end select
		
		!Debug configuration:
		case("DEBUG")
			select case(val)
				case(LCS_TRUE)
					LCS_VERBOSE = .TRUE.
					DEBUG_SGRID = .TRUE.
					str= 'TRUE'
				case(LCS_FALSE)
					LCS_VERBOSE = .FALSE.
					DEBUG_SGRID = .FALSE.
					str= 'FALSE'
				case default
					warn = .true.
			end select
		
		!Write Flowmap:
		case("WRITE_FLOWMAP")
			select case(val)
				case(LCS_TRUE)
					FLOWMAP_IO = .TRUE.
					str= 'TRUE'
				case(LCS_FALSE)
					FLOWMAP_IO = .FALSE.
					str= 'FALSE'
				case default
					warn = .true.
			end select
		
		!Write Bcflag:
		case("WRITE_BCFLAG")
			select case(val)
				case(LCS_TRUE)
					BCFLAG_IO = .TRUE.
					str= 'TRUE'
				case(LCS_FALSE)
					BCFLAG_IO = .FALSE.
					str= 'FALSE'
				case default
					warn = .true.
			end select
		
		!Compressibility:
		case("INCOMPRESSIBLE")
			select case(val)
				case(LCS_TRUE)
					INCOMPRESSIBLE = .TRUE.
					str= 'TRUE'
				case(LCS_FALSE)
					INCOMPRESSIBLE = .FALSE.
					str= 'FALSE'
				case default
					warn = .true.
			end select
		
		!Auxillary Grids (for CG Def tensor)
		case("AUX_GRID")
			select case(val)
				case(LCS_TRUE)
					AUX_GRID = .TRUE.
					str= 'TRUE'
				case(LCS_FALSE)
					AUX_GRID = .FALSE.
					str= 'FALSE'
				case default
					warn = .true.
			end select
		
		!Update frequency
		case("UPDATE_FREQ")
			N_UPDATE = val
			write(str,'(i32)') val

		!Interpolation
		case("INTERPOLATOR")
			select case(val)
				case(NEAREST_NBR)
					 str= 'NEAREST_NBR'
				case(LINEAR)
					 str= 'LINEAR'
				case(QUADRATIC)
					 str= 'QUADRATIC'
				case(CUBIC)
					 str= 'CUBIC'
				case(TSE)
					 str= 'TSE'
				case(TSE_LIMIT)
					 str= 'TSE_LIMIT'
				case default
					warn = .true.	
			end select
			INTERPOLATOR = val

		!Integration
		case("INTEGRATOR")
			select case(val)
				case(EULER)
					 str= 'EULER'
				case(TRAPEZOIDAL)
					 str= 'TRAPEZOIDAL'
				case(RK2)
					 str= 'RK2'
				case(RK3)
					 str= 'RK3'
				case(RK4)
					 str= 'RK4'
				case default
					warn = .true.	
			end select
			INTEGRATOR = val
		case default
			warn = .true.	
	end select
		
	!Tell the user what was set	
	if(lcsrank==0) then
		if(warn) then
			write(*,'(a,a,a,i10)') 'ERROR, bad libcfd2lcs option/val combination: ',trim(option),',',val
			CFD2LCS_ERROR = 1
		else
			write(*,'(a,a,a,a)') 'Setting libcfd2lcs ',trim(option),' = ',trim(str)
		endif
	endif
	
	!-----
	!Check
	!-----
	call cfd2lcs_error_check(error)

end subroutine cfd2lcs_set_option

subroutine cfd2lcs_set_param(param,val)
	use data_m
	implicit none
	!-----
	character(len=*)::param 
	real(LCSRP):: val
	!-----
	character(len=32):: str
	logical:: warn = .false.
	integer:: error
	!-----
	!Parse any parameters (real valued)  passed by the user
	!-----
	if (CFD2LCS_ERROR/=0) return

	if(lcsrank ==0 .AND. LCS_VERBOSE)&
		write(*,'(a)') 'in cfd2lcs_set_param... '

	select case(trim(param))
		!MAX CFL:
		case("CFL")
			CFL_MAX = val

		!Tracer injection X coordinate
		case("TRACER_INJECT_X")
			TRACER_INJECT_X = val

		!Tracer injection Y coordinate
		case("TRACER_INJECT_Y")
			TRACER_INJECT_Y = val

		!Tracer injection Z coordinate
		case("TRACER_INJECT_Z")
			TRACER_INJECT_Z = val

		!Tracer injection radius
		case("TRACER_INJECT_RADIUS")
			TRACER_INJECT_RADIUS = val

		case default
			warn = .true.	
	end select	
	
	!Tell the user what was set	
	if(lcsrank==0) then
		if(warn) then
			write(*,'(a,a)') 'ERROR, unknown libcfd2lcs parameter: ',trim(param)
			CFD2LCS_ERROR = 1
		else
			write(*,'(a,a,a,ES18.4)') 'Setting libcfd2lcs ',trim(param),' = ',val
		endif
	endif
	
	!-----
	!Check
	!-----
	call cfd2lcs_error_check(error)

end subroutine cfd2lcs_set_param

subroutine cfd2lcs_error_check(error)
	use data_m
	implicit none
	!-----
	integer:: error
	integer:: ierr,MAX_CFD2LCS_ERROR
	real(LCSRP):: t0,t1
	!-----
	!In the event of a cfd2lcs error, we dont necessarily want
	!to bring down the cfd solver.  So, flag an error instead:
	!-----
	t0 = cputimer(lcscomm,SYNC_TIMER)
	error= 0
	call MPI_ALLREDUCE(CFD2LCS_ERROR,MAX_CFD2LCS_ERROR,1,MPI_INTEGER,MPI_SUM,lcscomm,ierr)
	if (MAX_CFD2LCS_ERROR /= 0) then
		if (lcsrank==0) write(*,'(a)') &
			'FATAL CFD2LCS_ERROR DETECTED, WILL NOT PERFORM LCS COMPUTATIONS'
		CFD2LCS_ERROR = 1
		error = 1
	endif
	t0 = cputimer(lcscomm,SYNC_TIMER)
	cpu_error = cpu_error + max(t1-t0,0.0)
end subroutine cfd2lcs_error_check

subroutine cfd2lcs_info()
	use data_m
	implicit none
	integer:: ierr
	logical,save:: FIRST_CALL = .TRUE.
	!-----
	!Output cpu times and other efficiency related information:
	!-----
	if(CFD2LCS_ERROR /= 0) return
	
	!Initialize CPU timing
	if(FIRST_CALL) then
		t_start_global = cputimer(lcscomm,SYNC_TIMER)
		if(SYNC_TIMER) then
			call MPI_BCAST(t_start_global,1,MPI_REAL,0,lcscomm,ierr) !sync
		endif
		cpu_total_sim = 0.0
		cpu_total_lcs= 0.0
		cpu_fwd = 0.0
		cpu_bkwd = 0.0
		cpu_reconstruct = 0.0
		cpu_io = 0.0
		cpu_lpmap = 0.0
		cpu_ftle = 0.0
		cpu_error = 0.0
		
		this_cpu_fwd = 0.0
		this_cpu_bkwd =0.0
		cpu_fwd_c = 0.0
		cpu_bkwd_c =0.0
	
		FIRST_CALL = .FALSE.	
		return
	endif

	!Update...	
	cpu_total_lcs = cpu_total_lcs + max(t_finish_update-t_start_update,0.0)
	cpu_total_sim = cpu_total_sim + max(t_finish_update-t_start_global,0.0)
	integrations_fwd_c = integrations_fwd_c + integrations_fwd
	integrations_bkwd_c = integrations_bkwd_c + integrations_bkwd
	cpu_fwd_c = cpu_fwd_c+this_cpu_fwd
	cpu_bkwd_c = cpu_bkwd_c+this_cpu_bkwd
	t_start_global = t_finish_update
	if(SYNC_TIMER) then
		call MPI_BCAST(t_start_global,1,MPI_REAL,0,lcscomm,ierr) !sync
	endif
	if(lcsrank==0) then
		if(SYNC_TIMER)then
		write(*,'(a)') 	'-----libcfd2lcs CPU overhead (MPI synchronized)-------------------'
		else
		write(*,'(a)') 	'-----libcfd2lcs CPU overhead (NOT synchronized, rank 0 only)------'
		endif
		write(*,'(a)') 	'|      TASK             |  WALL CLOCK [SEC] |     % OF TOTAL     |'
		write(*,'(a,ES18.4,a,a,a)') &
						'|Application, Total:    |',&
		real(cpu_total_sim), ' | ', '         -         |'
		write(*,'(a,ES18.4,a,F17.4,a)') &
						'|cfd2lcs, Total:        |',&
		real(cpu_total_lcs), ' | ', cpu_total_lcs/cpu_total_sim*100.0_LCSRP,'% |'
		write(*,'(a,ES18.4,a,F17.4,a)') &
						'|cfd2lcs, FWD Advect:   |',&
		real(cpu_fwd), ' | ', cpu_fwd/cpu_total_sim*100.0_LCSRP,'% |'
		write(*,'(a,ES18.4,a,F17.4,a)') &
						'|cfd2lcs, BKWD Advect:  |',&
		real(cpu_bkwd), ' | ', cpu_bkwd/cpu_total_sim*100.0_LCSRP,'% |'
		write(*,'(a,ES18.4,a,F17.4,a)') &
						'|cfd2lcs, Reconstruct:  |',&
		real(cpu_reconstruct), ' | ', cpu_reconstruct/cpu_total_sim*100.0_LCSRP,'% |'
		write(*,'(a,ES18.4,a,F17.4,a)') &
						'|cfd2lcs, LP Re-map:    |',&
		real(cpu_lpmap), ' | ', cpu_lpmap/cpu_total_sim*100.0_LCSRP,'% |'
		write(*,'(a,ES18.4,a,F17.4,a)') &
						'|cfd2lcs, I/O:          |',&
		real(cpu_io), ' | ', cpu_io/cpu_total_sim*100.0_LCSRP,'% |'
		write(*,'(a,ES18.4,a,F17.4,a)') &
						'|cfd2lcs, LCS:          |',&
		real(cpu_ftle), ' | ', cpu_ftle/cpu_total_sim*100.0_LCSRP,'% |'
		write(*,'(a,ES18.4,a,F17.4,a)') &
						'|cfd2lcs, Error Check:  |',&
		real(cpu_error), ' | ', cpu_error/cpu_total_sim*100.0_LCSRP,'% |'
		write(*,'(a,ES18.4,a,F17.4,a)') &
						'|cfd2lcs, Other:        |',&
		real(cpu_total_lcs-(cpu_io+cpu_reconstruct+cpu_bkwd+cpu_fwd+cpu_lpmap+cpu_ftle+cpu_error)), ' | ',&
		real(cpu_total_lcs-(cpu_io+cpu_reconstruct+cpu_bkwd+cpu_fwd+cpu_lpmap+cpu_ftle+cpu_error))/cpu_total_sim*100.0_LCSRP, '% |'
		write(*,'(a)') 	'| PARTICLE INTEGRATIONS |  NUMBER [N/STEP]  |     RATE [N/SEC]   |'
		write(*,'(a,ES18.4,a,ES18.4,a)') &
						'|Forward (this step):   |',real(integrations_fwd),' | ', real(integrations_fwd,LCSRP)/this_cpu_fwd,' |'
		write(*,'(a,ES18.4,a,ES18.4,a)') &
						'|Backward (this step):  |',real(integrations_bkwd),' | ', real(integrations_bkwd,LCSRP)/this_cpu_bkwd,' |'
		write(*,'(a,ES18.4,a,ES18.4,a)') &
						'|Forward (cumulative):  |',real(integrations_fwd_c),' | ', real(integrations_fwd_c,LCSRP)/cpu_fwd_c,' |'
		write(*,'(a,ES18.4,a,ES18.4,a)') &
						'|Backward (cumulative): |',real(integrations_bkwd_c),' | ', real(integrations_bkwd_c,LCSRP)/cpu_bkwd_c,' |'
		write(*,'(a)') 	'------------------------------------------------------------------'
		integrations_fwd = 0
		integrations_bkwd = 0
		this_cpu_fwd = 0.0_LCSRP
		this_cpu_bkwd = 0.0_LCSRP
	endif

end subroutine cfd2lcs_info
