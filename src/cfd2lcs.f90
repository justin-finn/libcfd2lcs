!Top level, user interface module.
subroutine cfd2lcs_init(cfdcomm,ni,nj,nk,offset_i,offset_j,offset_k,x,y,z)
	use data_m
	use mpi_m
	implicit none
	!----
	integer:: cfdcomm
	integer:: ni,nj,nk,offset_i,offset_j,offset_k
	real(WP):: x(1:ni,1:nj,1:nk)
	real(WP):: y(1:ni,1:nj,1:nk)
	real(WP):: z(1:ni,1:nj,1:nk)
	!----
	integer:: nmax(3),my_nmax(3)
	integer:: ierr
	!----

	!init the mpi
	call init_lcs_mpi(cfdcomm)

	if(lcsrank ==0)&
		write(*,'(a)') 'in cfd2lcs_init...'

	!Init the default cfd storage structure:
	call init_scfd(scfd,'CFD Data',ni,nj,nk,offset_i,offset_j,offset_k,x,y,z)


	call cfd2lcs_error_check()

end subroutine cfd2lcs_init

subroutine cfd2lcs_update(ni,nj,nk,ux,uy,uz,time)
	use data_m
	use mpi_m
	use io_m
	implicit none
	!----
	integer:: ni,nj,nk
	real(WP), intent(in):: ux(1:ni,1:nj,1:nk)
	real(WP), intent(in):: uy(1:ni,1:nj,1:nk)
	real(WP), intent(in):: uz(1:ni,1:nj,1:nk)
	real(WP), intent(in):: time
	!----
	integer:: gn(3)
	integer:: offset(3)
	!----

	if(CFD2LCS_ERROR /= 0) return

	if(lcsrank ==0)&
		write(*,'(a)') 'in cfd2lcs_update...'

	!Check we got an arrays of the correct size:
	if	( scfd%ni/=ni .OR. scfd%nj /= nj .OR. scfd%nk /=nk) then
		write(*,'(a,i6,a)') 'rank[',lcsrank,'] received velocity array of incorrect dimension'
		write(*,'(a,i6,a,i4,i4,i4,a)') 'rank[',lcsrank,'] [ni,nj,nk]= [',ni,nj,nk,']'
		write(*,'(a,i6,a,i4,i4,i4,a)') 'rank[',lcsrank,'] scfd[ni,nj,nk]= [',scfd%ni,scfd%nj,scfd%nk,']'
		CFD2LCS_ERROR = 2
	endif

	!Set the velocity:
	scfd%u%x(1:ni,1:nj,1:nk) = ux(1:ni,1:nj,1:nk)
	scfd%u%y(1:ni,1:nj,1:nk) = uy(1:ni,1:nj,1:nk)
	scfd%u%z(1:ni,1:nj,1:nk) = uz(1:ni,1:nj,1:nk)


	!Test the I/O
	gn = (/scfd%gni,scfd%gnj,scfd%gnk/)
	offset = (/scfd%offset_i,scfd%offset_j,scfd%offset_k/)
	call structured_io('./dump/iotest.h5',IO_WRITE,gn,offset,r1=scfd%u)		!Write U
	call structured_io('./dump/iotest.h5',IO_APPEND,gn,offset,r1=scfd%grid)	!Append the grid



	call cfd2lcs_error_check()
end subroutine cfd2lcs_update

subroutine cfd2lcs_finalize()
	use data_m
	use mpi_m
	implicit none

	if(lcsrank ==0)&
		write(*,'(a)') 'in cfd2lcs_finalize...'

	call destroy_scfd(scfd)
end subroutine cfd2lcs_finalize
