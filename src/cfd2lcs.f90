!Top level, user interface module.
subroutine cfd2lcs_init(cfdcomm,n,offset,x,y,z,BC_LIST,lperiodic)
	use scfd_m
	use mpi_m
	implicit none
	!----
	integer:: cfdcomm
	integer:: n(3),offset(3)
	real(LCSRP):: x(1:n(1),1:n(2),1:n(3))
	real(LCSRP):: y(1:n(1),1:n(2),1:n(3))
	real(LCSRP):: z(1:n(1),1:n(2),1:n(3))
	integer(LCSIP),dimension(6):: BC_LIST
	real(LCSRP):: lperiodic(3)
	!----
	integer:: error
	!----

	!Error handling:
	CFD2LCS_ERROR = 0

	!init the mpi
	call init_lcs_mpi(cfdcomm)

	if(lcsrank ==0)&
		write(*,'(a)') 'in cfd2lcs_init...'

	!Init the default cfd storage structure:
	call init_scfd(scfd,'CFD Data',n,offset,x,y,z,BC_LIST,lperiodic)

	!Check:
	call cfd2lcs_error_check(error)

end subroutine cfd2lcs_init

subroutine cfd2lcs_update(n,ux,uy,uz,time)
	use data_m
	use io_m
	use comms_m
	use scfd_m
	implicit none
	!----
	integer:: n(3)
	real(LCSRP):: ux(1:n(1),1:n(2),1:n(3))
	real(LCSRP):: uy(1:n(1),1:n(2),1:n(3))
	real(LCSRP):: uz(1:n(1),1:n(2),1:n(3))
	real(LCSRP), intent(in):: time
	!----
	integer:: gn(3)
	integer:: offset(3)
	integer:: error

	type(sr2_t):: gradu
	!----

	if(CFD2LCS_ERROR /= 0) return

	if(lcsrank ==0)&
		write(*,'(a)') 'in cfd2lcs_update...'

	!Check we got an arrays of the correct size:
	if	( scfd%ni/=n(1) .OR. scfd%nj /= n(2) .OR. scfd%nk /=n(3)) then
		write(*,'(a,i6,a)') 'rank[',lcsrank,'] received velocity array of incorrect dimension'
		write(*,'(a,i6,a,i4,i4,i4,a)') 'rank[',lcsrank,'] [ni,nj,nk]= [',n(1),n(2),n(3),']'
		write(*,'(a,i6,a,i4,i4,i4,a)') 'rank[',lcsrank,'] scfd[ni,nj,nk]= [',scfd%ni,scfd%nj,scfd%nk,']'
		CFD2LCS_ERROR = 1
	endif

	!Set the velocity:
	scfd%u%x(1:n(1),1:n(2),1:n(3)) = ux(1:n(1),1:n(2),1:n(3))
	scfd%u%y(1:n(1),1:n(2),1:n(3)) = uy(1:n(1),1:n(2),1:n(3))
	scfd%u%z(1:n(1),1:n(2),1:n(3)) = uz(1:n(1),1:n(2),1:n(3))

	!Compute grad(U):
	call init_sr2(gradu,scfd%ni,scfd%nj,scfd%nk,scfd%ng,'Grad U')
	call exchange_sdata(scfd%scomm_max_r1,r1=scfd%u)
	call grad_sr1(scfd%ni,scfd%nj,scfd%nk,scfd%ng,scfd%grid,scfd%u,gradu)

	!Test the I/O
	gn = (/scfd%gni,scfd%gnj,scfd%gnk/)
	offset = (/scfd%offset_i,scfd%offset_j,scfd%offset_k/)
	call structured_io('./dump/iotest.h5',IO_WRITE,gn,offset,r1=scfd%u)		!Write U
	call structured_io('./dump/iotest.h5',IO_APPEND,gn,offset,r1=scfd%grid)	!Append the grid
	call structured_io('./dump/iotest.h5',IO_APPEND,gn,offset,r2=gradu)	!Append grad(U)

	!cleanup
	call destroy_sr2(gradu)

	call cfd2lcs_error_check(error)

end subroutine cfd2lcs_update

subroutine cfd2lcs_finalize()
	use data_m
	use mpi_m
	use scfd_m
	implicit none

	if(lcsrank ==0)&
		write(*,'(a)') 'in cfd2lcs_finalize...'

	call destroy_scfd(scfd)
end subroutine cfd2lcs_finalize

subroutine cfd2lcs_error_check(error)
	use data_m
	use mpi_m
	implicit none
	!-----
	integer:: error
	integer:: ierr,MAX_CFD2LCS_ERROR
	!-----
	!In the event of a cfd2lcs error, we dont necessarily want
	!to bring down the cfd solver.  So, flag an error instead:
	!-----
	error= 0
	call MPI_ALLREDUCE(CFD2LCS_ERROR,MAX_CFD2LCS_ERROR,1,MPI_INTEGER,MPI_SUM,lcscomm,ierr)
	if (MAX_CFD2LCS_ERROR /= 0) then
		if (lcsrank==0) write(*,'(a)') &
			'FATAL CFD2LCS_ERROR DETECTED, WILL NOT PERFORM LCS COMPUTATIONS'
		CFD2LCS_ERROR = 1
		error = 1
	endif
end subroutine cfd2lcs_error_check

