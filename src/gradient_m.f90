!
!Copyright (C) 2015-2016, Justin R. Finn.  All rights reserved.
!libcfd2lcs is distributed is under the terms of the GNU General Public License
!
module gradient_m
	use data_m
	use structured_m
	implicit none
	contains

	subroutine grad_sr1(sgrid,sr1,grad)
		implicit none
		!----
		type(sgrid_t),intent(in):: sgrid
		type(sr1_t),intent(in):: sr1
		type(sr2_t),intent(inout):: grad
		!----
		integer:: i,j,k
		!----

		if (lcsrank==0 .AND. LCS_VERBOSE) &
			write(*,*) 'in grad_sr1... ',trim(sr1%label),' => ',trim(grad%label)

		!Check that the grid is rectilinear:
		if(.NOT. sgrid%rectilinear) then
			if(lcsrank==0 .and. LCS_VERBOSE) &
				write(*,*) 'WARNING, in grad_sr1: sgrid is not rectilinear, calling grad_sr1_ls instead'
			call grad_sr1_ls(sgrid,sr1,grad)
			return
		endif

		!2nd order central scheme:
		do k = 1,sgrid%nk
		do j = 1,sgrid%nj
		do i = 1,sgrid%ni
			grad%xx(i,j,k) = (sr1%x(i+1,j,k)-sr1%x(i-1,j,k)) / (sgrid%grid%x(i+1,j,k) - sgrid%grid%x(i-1,j,k))
			grad%xy(i,j,k) = (sr1%x(i,j+1,k)-sr1%x(i,j-1,k)) / (sgrid%grid%y(i,j+1,k) - sgrid%grid%y(i,j-1,k))
			grad%xz(i,j,k) = (sr1%x(i,j,k+1)-sr1%x(i,j,k-1)) / (sgrid%grid%z(i,j,k+1) - sgrid%grid%z(i,j,k-1))
			grad%yx(i,j,k) = (sr1%y(i+1,j,k)-sr1%y(i-1,j,k)) / (sgrid%grid%x(i+1,j,k) - sgrid%grid%x(i-1,j,k))
			grad%yy(i,j,k) = (sr1%y(i,j+1,k)-sr1%y(i,j-1,k)) / (sgrid%grid%y(i,j+1,k) - sgrid%grid%y(i,j-1,k))
			grad%yz(i,j,k) = (sr1%y(i,j,k+1)-sr1%y(i,j,k-1)) / (sgrid%grid%z(i,j,k+1) - sgrid%grid%z(i,j,k-1))
			grad%zx(i,j,k) = (sr1%z(i+1,j,k)-sr1%z(i-1,j,k)) / (sgrid%grid%x(i+1,j,k) - sgrid%grid%x(i-1,j,k))
			grad%zy(i,j,k) = (sr1%z(i,j+1,k)-sr1%z(i,j-1,k)) / (sgrid%grid%y(i,j+1,k) - sgrid%grid%y(i,j-1,k))
			grad%zz(i,j,k) = (sr1%z(i,j,k+1)-sr1%z(i,j,k-1)) / (sgrid%grid%z(i,j,k+1) - sgrid%grid%z(i,j,k-1))
		enddo
		enddo
		enddo

		!1D Checks (for rectilinear):
		if(sgrid%gni==1) then
			grad%xx = 0.0_LCSRP
			grad%yx = 0.0_LCSRP
			grad%zx = 0.0_LCSRP
		endif
		if(sgrid%gnj==1) then
			grad%xy = 0.0_LCSRP
			grad%yy = 0.0_LCSRP
			grad%zy = 0.0_LCSRP
		endif
		if(sgrid%gnk==1) then
			grad%xz = 0.0_LCSRP
			grad%yz = 0.0_LCSRP
			grad%zz = 0.0_LCSRP
		endif

	end subroutine grad_sr1

	subroutine grad_sr0(sgrid,sr0,grad)
		implicit none
		!----
		type(sgrid_t),intent(in):: sgrid
		type(sr0_t),intent(in):: sr0
		type(sr1_t),intent(inout):: grad
		!----
		integer:: i,j,k
		!----

		if (lcsrank==0 .AND. LCS_VERBOSE) &
			write(*,*) 'in grad_sr0... ',trim(sr0%label),' => ',trim(grad%label)

		!Check that the grid is rectilinear:
		if(.NOT. sgrid%rectilinear) then
			if(lcsrank==0 .and. LCS_VERBOSE) &
				write(*,*) 'WARNING, in grad_sr0: sgrid is not rectilinear, calling grad_sr0_ls instead'
			call grad_sr0_ls(sgrid,sr0,grad)
			return
		endif

		!2nd order central scheme:
		do k = 1,sgrid%nk
		do j = 1,sgrid%nj
		do i = 1,sgrid%ni
			grad%x(i,j,k) = (sr0%r(i+1,j,k)-sr0%r(i-1,j,k)) / (sgrid%grid%x(i+1,j,k) - sgrid%grid%x(i-1,j,k))
			grad%y(i,j,k) = (sr0%r(i,j+1,k)-sr0%r(i,j-1,k)) / (sgrid%grid%y(i,j+1,k) - sgrid%grid%y(i,j-1,k))
			grad%z(i,j,k) = (sr0%r(i,j,k+1)-sr0%r(i,j,k-1)) / (sgrid%grid%z(i,j,k+1) - sgrid%grid%z(i,j,k-1))
		enddo
		enddo
		enddo

		!1D Checks (for rectilinear):
		if(sgrid%gni==1) then
			grad%x = 0.0_LCSRP
		endif
		if(sgrid%gnj==1) then
			grad%y = 0.0_LCSRP
		endif
		if(sgrid%gnk==1) then
			grad%z = 0.0_LCSRP
		endif

	end subroutine grad_sr0


	subroutine compute_lsgw(sgrid,full_conn)
		implicit none
		!-----
		type(sgrid_t):: sgrid
		logical:: full_conn
		!-----
		integer:: i,j,k,ii,jj,kk,nbr
		real(lcsrp):: sum_wdx2,sum_wdy2,sum_wdz2,sum_wdxdy,sum_wdxdz,sum_wdydz,weight,dx(3),denom
		character(len=32):: label
		!-----
		! Compute and save the weights used for least-squares gradient
		! estimation for non-rectilinear grids.
		!-----
		if (lcsrank == 0) &
			write(*,*) 'in calc_lsg_weights...', trim(sgrid%label),full_conn

		!cleanup old wts, if any:
		call destroy_lsgw(sgrid)

		!Determine the nbr range depending on the desired connectivity
		if(full_conn) then
			sgrid%nbr_f= 2
			sgrid%nbr_l= 27
		else
			sgrid%nbr_f=2
			sgrid%nbr_l=7
		endif

		!Allocate space for the weights:
		allocate(sgrid%lsgw(sgrid%nbr_f:sgrid%nbr_l))
		do nbr = sgrid%nbr_f,sgrid%nbr_l
			write(label,'(a,i2.2)') 'LSG_WTS_',nbr
			call init_sr1(sgrid%lsgw(nbr),sgrid%ni,sgrid%nj,sgrid%nk,sgrid%ng,trim(label),translate=.false.)
		enddo

		do k= 1,sgrid%nk
		do j= 1,sgrid%nj
		do i= 1,sgrid%ni
			sum_wdx2 = 0.0_LCSRP
			sum_wdy2 = 0.0_LCSRP
			sum_wdz2 = 0.0_LCSRP
			sum_wdxdy = 0.0_LCSRP
			sum_wdxdz = 0.0_LCSRP
			sum_wdydz = 0.0_LCSRP

			do nbr = sgrid%nbr_f,sgrid%nbr_l
				ii = i+NBR_OFFSET(1,nbr)
				jj = j+NBR_OFFSET(2,nbr)
				kk = k+NBR_OFFSET(3,nbr)

				dx(1) = sgrid%grid%x(ii,jj,kk) - sgrid%grid%x(i,j,k)
				dx(2) = sgrid%grid%y(ii,jj,kk) - sgrid%grid%y(i,j,k)
				dx(3) = sgrid%grid%z(ii,jj,kk) - sgrid%grid%z(i,j,k)
				!Inverse distance weight:
				weight = 1.0_LCSRP/sum(dx(1:3)**2)

				sum_wdx2 = sum_wdx2 + weight*dx(1)**2
				sum_wdy2 = sum_wdy2 + weight*dx(2)**2
				sum_wdz2 = sum_wdz2 + weight*dx(3)**2
				sum_wdxdy = sum_wdxdy + weight*dx(1)*dx(2)
				sum_wdxdz = sum_wdxdz + weight*dx(1)*dx(3)
				sum_wdydz = sum_wdydz + weight*dx(2)*dx(3)
			enddo

			denom = 2.0_LCSRP*sum_wdxdy*sum_wdxdz*sum_wdydz + &
				sum_wdx2*sum_wdy2*sum_wdz2 - &
				sum_wdx2*sum_wdydz**2 - &
				sum_wdy2*sum_wdxdz**2 - &
				sum_wdz2*sum_wdxdy**2

			do nbr = sgrid%nbr_f,sgrid%nbr_l
				ii = i+NBR_OFFSET(1,nbr)
				jj = j+NBR_OFFSET(2,nbr)
				kk = k+NBR_OFFSET(3,nbr)
				dx(1) = sgrid%grid%x(ii,jj,kk) - sgrid%grid%x(i,j,k)
				dx(2) = sgrid%grid%y(ii,jj,kk) - sgrid%grid%y(i,j,k)
				dx(3) = sgrid%grid%z(ii,jj,kk) - sgrid%grid%z(i,j,k)

				!Inverse distance weight:
				weight = 1.0_LCSRP/sum(dx(1:3)**2)
				! x
				sgrid%lsgw(nbr)%x(i,j,k) = weight*( &
					(sum_wdy2*sum_wdz2-sum_wdydz**2)*dx(1) + &
					(sum_wdxdz*sum_wdydz-sum_wdxdy*sum_wdz2)*dx(2) + &
					(sum_wdxdy*sum_wdydz-sum_wdxdz*sum_wdy2)*dx(3) )/denom
				! y
				sgrid%lsgw(nbr)%y(i,j,k) = weight*( &
					(sum_wdxdz*sum_wdydz-sum_wdxdy*sum_wdz2)*dx(1) + &
					(sum_wdx2*sum_wdz2-sum_wdxdz**2)*dx(2) + &
					(sum_wdxdy*sum_wdxdz-sum_wdydz*sum_wdx2)*dx(3) )/denom
				! z
				sgrid%lsgw(nbr)%z(i,j,k) = weight*( &
					(sum_wdxdy*sum_wdydz-sum_wdxdz*sum_wdy2)*dx(1) + &
					(sum_wdxdy*sum_wdxdz-sum_wdydz*sum_wdx2)*dx(2) + &
					(sum_wdx2*sum_wdy2-sum_wdxdy**2)*dx(3) )/denom
			end do
		enddo
		enddo
		enddo

	end subroutine compute_lsgw

	subroutine destroy_lsgw(sgrid)
		implicit none
		!-----
		type(sgrid_t):: sgrid
		!-----
		integer:: nbr
		!-----

		if(.NOT. allocated(sgrid%lsgw)) return

		do nbr = sgrid%nbr_f,sgrid%nbr_l
			call destroy_sr1(sgrid%lsgw(nbr))
		enddo
		if(allocated(sgrid%lsgw)) then
			deallocate(sgrid%lsgw)
		endif

		sgrid%nbr_f = 0
		sgrid%nbr_l = 0

	end subroutine destroy_lsgw

	subroutine grad_sr1_ls(sgrid,sr1,grad)
		implicit none
		!----
		type(sgrid_t),intent(in):: sgrid
		type(sr1_t),intent(in):: sr1
		type(sr2_t),intent(inout):: grad
		!----
		integer:: i,j,k
		integer:: ni,nj,nk,ng
		integer:: nbr
		type(sr1_t):: tmp
		!----

		if (lcsrank==0 .AND. LCS_VERBOSE) &
			write(*,*) 'in grad_sr1_ls... ',trim(sr1%label),' => ',trim(grad%label)

		grad%xx = 0.0_LCSRP
		grad%xy = 0.0_LCSRP
		grad%xz = 0.0_LCSRP
		grad%yx = 0.0_LCSRP
		grad%yy = 0.0_LCSRP
		grad%yz = 0.0_LCSRP
		grad%zx = 0.0_LCSRP
		grad%zy = 0.0_LCSRP
		grad%zz = 0.0_LCSRP

		ni = sgrid%ni
		nj = sgrid%nj
		nk = sgrid%nk
		ng = sgrid%ng
		call init_sr1(tmp,ni,nj,nk,ng,'TMP',translate=.false.)
		!These should all vectorize (confirmed with gfortran)
		do nbr = sgrid%nbr_f,sgrid%nbr_l
			i = NBR_OFFSET(1,nbr)
			j = NBR_OFFSET(2,nbr)
			k = NBR_OFFSET(3,nbr)
			tmp%x(1:ni,1:nj,1:nk) = sr1%x(1+i:ni+i, 1+j:nj+j, 1+k:nk+k)
			tmp%y(1:ni,1:nj,1:nk) = sr1%y(1+i:ni+i, 1+j:nj+j, 1+k:nk+k)
			tmp%z(1:ni,1:nj,1:nk) = sr1%z(1+i:ni+i, 1+j:nj+j, 1+k:nk+k)
			tmp%x = tmp%x-sr1%x
			grad%xx = grad%xx + tmp%x*sgrid%lsgw(nbr)%x
			grad%xy = grad%xy + tmp%x*sgrid%lsgw(nbr)%y
			grad%xz = grad%xz + tmp%x*sgrid%lsgw(nbr)%z
			tmp%y = tmp%y-sr1%y
			grad%yx = grad%yx + tmp%y*sgrid%lsgw(nbr)%x
			grad%yy = grad%yy + tmp%y*sgrid%lsgw(nbr)%y
			grad%yz = grad%yz + tmp%y*sgrid%lsgw(nbr)%z
			tmp%z = tmp%z-sr1%z
			grad%zx = grad%zx + tmp%z*sgrid%lsgw(nbr)%x
			grad%zy = grad%zy + tmp%z*sgrid%lsgw(nbr)%y
			grad%zz = grad%zz + tmp%z*sgrid%lsgw(nbr)%z
		enddo
		call destroy_sr1(tmp)

	end subroutine grad_sr1_ls

	subroutine grad_sr0_ls(sgrid,sr0,grad)
		implicit none
		!----
		type(sgrid_t),intent(in):: sgrid
		type(sr0_t),intent(in):: sr0
		type(sr1_t),intent(inout):: grad
		!----
		integer:: i,j,k
		integer:: ni,nj,nk,ng
		integer:: nbr
		type(sr0_t):: tmp
		!----

		if (lcsrank==0 .AND. LCS_VERBOSE) &
			write(*,*) 'in grad_sr0_ls... ',trim(sr0%label),' => ',trim(grad%label)

		grad%x = 0.0_LCSRP
		grad%y = 0.0_LCSRP
		grad%z = 0.0_LCSRP

		ni = sgrid%ni
		nj = sgrid%nj
		nk = sgrid%nk
		ng = sgrid%ng
		call init_sr0(tmp,ni,nj,nk,ng,'TMP')
		!These should all vectorize (confirmed with gfortran)
		do nbr = sgrid%nbr_f,sgrid%nbr_l
			i = NBR_OFFSET(1,nbr)
			j = NBR_OFFSET(2,nbr)
			k = NBR_OFFSET(3,nbr)
			tmp%r(1:ni,1:nj,1:nk) = sr0%r(1+i:ni+i, 1+j:nj+j, 1+k:nk+k)
			tmp%r = tmp%r-sr0%r
			grad%x = grad%x + tmp%r*sgrid%lsgw(nbr)%x
			grad%y = grad%y + tmp%r*sgrid%lsgw(nbr)%y
			grad%z = grad%z + tmp%r*sgrid%lsgw(nbr)%z
		enddo
		call destroy_sr0(tmp)

	end subroutine grad_sr0_ls

end module gradient_m
