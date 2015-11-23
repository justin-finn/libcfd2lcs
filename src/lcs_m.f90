module lcs_m
	use data_m
	use structured_m
	use comms_m
	implicit none

	contains

	subroutine compute_ftle(lcs)
		implicit none
		!-----
		type(lcs_t):: lcs
		!-----
		type(sr2_t):: cg,gradfm
		integer:: ni,nj,nk,ng
		integer:: i,j,k
		integer,parameter::LWMAX=1000
		real(LCSRP):: A(3,3), W(3), WORK(LWMAX)
		integer:: INFO,LWORK
		!-----

		if(lcsrank==0)&
			write(*,*) 'in compute_ftle... ', trim(lcs%fm%label),' => ',trim(lcs%ftle%label)

		!brevity
		ni = lcs%fm%ni
		nj = lcs%fm%nj
		nk = lcs%fm%nk
		ng = lcs%fm%ng

		!Compute grad(fm)
		call init_sr2(gradfm,ni,nj,nk,ng,'GradFM')
		call init_sr2(cg,ni,nj,nk,ng,'CGtensor')
		call exchange_sdata(lcs%sgrid%scomm_face_r1,r1=lcs%fm)
		call grad_sr1(ni,nj,nk,ng,lcs%sgrid%grid,lcs%fm,gradfm)

		!Compute the right Cauchy-Green deformation tensor
		cg%xx = gradfm%xx*gradfm%xx + gradfm%yx*gradfm%yx + gradfm%zx*gradfm%zx
		cg%xy = gradfm%xx*gradfm%xy + gradfm%yx*gradfm%yy + gradfm%zx*gradfm%zy
		cg%xz = gradfm%xx*gradfm%xz + gradfm%yx*gradfm%yz + gradfm%zx*gradfm%zz
		cg%yx = gradfm%xy*gradfm%xx + gradfm%yy*gradfm%yx + gradfm%zy*gradfm%zx
		cg%yy = gradfm%xy*gradfm%xy + gradfm%yy*gradfm%yy + gradfm%zy*gradfm%zy
		cg%yz = gradfm%xy*gradfm%xz + gradfm%yy*gradfm%yz + gradfm%zy*gradfm%zz
		cg%zx = gradfm%xz*gradfm%xx + gradfm%yz*gradfm%yx + gradfm%zz*gradfm%zx
		cg%zy = gradfm%xz*gradfm%xy + gradfm%yz*gradfm%yy + gradfm%zz*gradfm%zy
		cg%zz = gradfm%xz*gradfm%xz + gradfm%yz*gradfm%yz + gradfm%zz*gradfm%zz

		!For each node, solve the eigensystem using LAPACK
		do k = 1,nk
		do j = 1,nj
		do i = 1,ni
			A(1,1) = cg%xx(i,j,k)
			A(1,2) = cg%xy(i,j,k)
			A(1,3) = cg%xz(i,j,k)
			A(2,1) = cg%yx(i,j,k)
			A(2,2) = cg%yy(i,j,k)
			A(2,3) = cg%yz(i,j,k)
			A(3,1) = cg%zx(i,j,k)
			A(3,2) = cg%zy(i,j,k)
			A(3,3) = cg%zz(i,j,k)
			select case(LCSRP)
			case(4)
				if(i == 1 .AND. j==1 .AND. k==1) then
					LWORK= -1
					call SSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
					LWORK = MIN(LWMAX,INT(WORK(1)))
				endif
				call SSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
				lcs%ftle%r(i,j,k) = log(maxval(W))
			case(8)
				if(i == 1 .AND. j==1 .AND. k==1) then
					LWORK= -1
					call DSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
					LWORK = MIN(LWMAX,INT(WORK(1)))
				endif
				call DSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
				lcs%ftle%r(i,j,k) = log(maxval(W))
			end select
		enddo
		enddo
		enddo

		call destroy_sr2(gradfm)
		call destroy_sr2(cg)

	end subroutine compute_ftle

end module lcs_m
