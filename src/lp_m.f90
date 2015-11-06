module lp_m
	use data_m
	use unstructured_m
	implicit none
	!Basic routines for the lp structure.	
	contains

	subroutine init_lp(lp,label,rhop,dp,grid)
		!----
		type(lp_t),pointer:: lp
		character(len=*):: label
		real(LCSRP) rhop,dp
		type(sr1_t):: grid
		!----
		integer:: i,j,k,ip
		integer:: ilp, ilcs
		integer, allocatable::lpptr(:)
		type(lp_t),allocatable::lp_c_tmp(:)
		!----
		!Initialize a new set of lagrangian particles with certain properties
		!Initial positions will be determined by the grid passed in.
		!----

		if (lcsrank==0)&
			write(*,*) 'in init_lp...',trim(label)

		!----
		!Add a new item to the lp array
		!Careful to preserve ptrs from LCS
		!----
		if(NLP == 0 ) then
			NLP = NLP + 1
			allocate(lp_c(NLP))
		else
			!Save ptrs from lcs
			allocate(lpptr(1:NLCS))
			lpptr = -1
			do ilcs = 1,NLCS
				do ilp = 1,NLP
					if (associated(lcs_c(ilcs)%lp,lp_c(ilp))) then
						lpptr(ilcs) = ilp
					endif
				enddo
			enddo

			!expand array of structures
			allocate(lp_c_tmp(NLP))
			lp_c_tmp = lp_c
			deallocate(lp_c)
			allocate(lp_c(NLP+1))
			lp_c(1:NLP) = lp_c_tmp(1:NLP)

			!fix old ptrs
			do ilcs = 1,NLCS
				if(lpptr(ilcs) > 0) then
					lcs_c(ilcs)%lp => lp_c(lpptr(ilcs))
				endif
			enddo

			NLP = NLP + 1
		endif

		!Point to the new lp and initialize
		lp => lp_c(NLP)
		lp%label = trim(label)
		lp%np = grid%ni*grid%nj*grid%nk
		lp%rhop = rhop
		lp%dp = dp
		call init_ur1(lp%xp,lp%np,'XP')
		call init_ur1(lp%up,lp%np,'UP')
		call init_ui1(lp%no,lp%np,'NODE')
		call init_ui0(lp%no0,lp%np,'NODE0')
		call init_ui0(lp%proc0,lp%np,'PROC0')
		call init_ui0(lp%flag,lp%np,'FLAG')
		ip = 0
		do k = 1,grid%nk
		do j = 1,grid%nj
		do i = 1,grid%ni
			ip = ip + 1
			lp%xp%x(ip) = grid%x(i,j,k)
			lp%xp%y(ip) = grid%y(i,j,k)
			lp%xp%z(ip) = grid%z(i,j,k)
			lp%proc0%i(ip) = lcsrank
			lp%no0%i(ip) = ijk2l(i,j,k,grid%ni,grid%nj)  !a guess that works if the grid is the same
			lp%flag%i(ip) = LP_UNKNOWN
		enddo
		enddo
		enddo

		!cleanup
		if(allocated(lpptr))deallocate(lpptr)
	end subroutine init_lp

	subroutine resize_lp(lp,np)
		implicit none
		!-----
		type(lp_t):: lp
		integer:: np
		!-----
		!Dynamically resize lp
		!-----

		if(LCS_VERBOSE)&
			write(*,*) 'in resize_lp...',trim(lp%label)

		!Check the size of the arrays and resize if needed:
		call resize_ur1(lp%xp,np)
		call resize_ur1(lp%up,np)
		call resize_ui1(lp%no,np)
		call resize_ui0(lp%proc0,np)
		call resize_ui0(lp%no0,np)
		call resize_ui0(lp%flag,np)

		!Set the new np
		lp%np = np
	
	end subroutine resize_lp

	subroutine reorder_lp(lp)
		implicit none
		!-------
		type(lp_t):: lp
		!-------
		integer:: ip, new_np
		integer:: recycle_count
		!-------
		!Clean up any holes in the lp structure
		!This needs to be called after processors
		!exchange lp's, or if LP's leave the domain,
		!and before unpacking any new ones.
		!It is up to the user to resize afterwards.
		!-------

		if (lcsrank==0) &
			write(*,*) 'in reorder_lp...',trim(lp%label)

		new_np = 0
		recycle_count = 0
		do ip = 1,lp%np
			if (lp%flag%i(ip) == LP_RECYCLE) then
				recycle_count = recycle_count+1
			else
				new_np = new_np+1
				if (new_np/=ip) then
					lp%xp%x(new_np) 		= lp%xp%x(ip)
					lp%xp%y(new_np) 		= lp%xp%y(ip)
					lp%xp%z(new_np) 		= lp%xp%z(ip)
					lp%up%x(new_np) 		= lp%up%x(ip)
					lp%up%y(new_np) 		= lp%up%y(ip)
					lp%up%z(new_np) 		= lp%up%z(ip)
					lp%no%x(new_np) 		= lp%no%x(ip)
					lp%no%y(new_np) 		= lp%no%y(ip)
					lp%no%z(new_np) 		= lp%no%z(ip)
					lp%proc0%i(new_np) 		= lp%proc0%i(ip)
					lp%no0%i(new_np) 		= lp%no0%i(ip)
					lp%flag%i(new_np) 		= lp%flag%i(ip)
				endif
			endif
		enddo

		!Note, no resize done here, but need to set n for all members of lp.
		lp%np = new_np
		lp%xp%n = new_np
		lp%up%n = new_np
		lp%no%n = new_np
		lp%proc0%n = new_np
		lp%no0%n = new_np
		lp%flag%n = new_np

	end subroutine reorder_lp

end module lp_m
