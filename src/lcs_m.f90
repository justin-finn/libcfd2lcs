!
!Copyright (C) 2015-2016, Justin R. Finn.  All rights reserved.
!libcfd2lcs is distributed is under the terms of the GNU General Public License
!
module lcs_m
      use data_m
      use sgrid_m
      use comms_m
      use gradient_m
      use lp_m
      use lp_tracking_m
      implicit none

      contains

      subroutine init_ftle(lcs,label,direction,resolution)
            implicit none
            !-----
            type(lcs_t),pointer:: lcs
            character(len=*),intent(in):: label
            integer:: direction, resolution
            !-----
            integer:: ni,nj,nk,ng
            real(LCSRP):: dt_factor
            !-----
            !Initialize the FTLE type lcs
            !-----

            !brevity:
            ni = lcs%sgrid%ni
            nj = lcs%sgrid%nj
            nk = lcs%sgrid%nk
            ng = lcs%sgrid%ng

            !The ftle:
            call init_sr0(lcs%ftle,ni,nj,nk,ng,trim(lcs%label)//'-FTLE')

            !Set the timestepping factor:
            !This helps decide the dt for the semi-lagrangian update of the flowmap
            if(direction == BKWD) then
                  dt_factor = 1.0_LCSRP/real(max(1+resolution, 1))
            else
                  dt_factor = 1.0_LCSRP
            endif
                  
            !Initialize the particles on the lcs grid:
            call init_lp(lcs%lp,lcs%sgrid,direction,dt_factor,trim(label)//'-LP')
            call track_lp2node(lcs%lp,scfd%sgrid,lcs%lp%no_scfd)

            !If you are using an Auxilllary grid, you will create 6 new grids of particles:
            if(AUX_GRID) then
                  call create_aux_grid(lcs)
            endif

      end subroutine init_ftle

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
            real(LCSRP),parameter:: SMALL = 100.0_LCSRP*tiny(1.0_LCSRP)
            real(LCSRP):: t0,t1
            !-----

            if(lcsrank==0)&
                  write(*,'(a,a)') 'Computing FTLE: ' ,trim(lcs%ftle%label)
                  
            t0 = cputimer(lcscomm,SYNC_TIMER)

            !brevity
            ni = lcs%sgrid%ni
            nj = lcs%sgrid%nj
            nk = lcs%sgrid%nk
            ng = lcs%sgrid%ng

            !Compute grad(fm)
            call init_sr2(gradfm,ni,nj,nk,ng,'GradFM')
            call init_sr2(cg,ni,nj,nk,ng,'CGtensor')
            if(AUX_GRID) then
                  !Careful to avoid div by zero:
                  do k = 1,nk
                  do j = 1,nj
                  do i = 1,ni
                        if( abs(lcs%lpX1%sgrid%grid%x(i,j,k)-lcs%lpX0%sgrid%grid%x(i,j,k))>SMALL) then
                              gradfm%xx(i,j,k) = (lcs%lpX1%fm%x(i,j,k)-lcs%lpX0%fm%x(i,j,k))&
                               / (lcs%lpX1%sgrid%grid%x(i,j,k)-lcs%lpX0%sgrid%grid%x(i,j,k))
                              gradfm%yx(i,j,k) = (lcs%lpX1%fm%y(i,j,k)-lcs%lpX0%fm%y(i,j,k))&
                               / (lcs%lpX1%sgrid%grid%x(i,j,k)-lcs%lpX0%sgrid%grid%x(i,j,k))
                              gradfm%zx(i,j,k) = (lcs%lpX1%fm%z(i,j,k)-lcs%lpX0%fm%z(i,j,k))&
                               / (lcs%lpX1%sgrid%grid%x(i,j,k)-lcs%lpX0%sgrid%grid%x(i,j,k))
                        endif             
                        if( abs(lcs%lpY1%sgrid%grid%y(i,j,k)-lcs%lpY0%sgrid%grid%y(i,j,k))>SMALL) then
                              gradfm%xy(i,j,k) = (lcs%lpY1%fm%x(i,j,k)-lcs%lpY0%fm%x(i,j,k))&
                               / (lcs%lpY1%sgrid%grid%y(i,j,k)-lcs%lpY0%sgrid%grid%y(i,j,k))
                              gradfm%yy(i,j,k) = (lcs%lpY1%fm%y(i,j,k)-lcs%lpY0%fm%y(i,j,k))&
                               / (lcs%lpY1%sgrid%grid%y(i,j,k)-lcs%lpY0%sgrid%grid%y(i,j,k))
                              gradfm%zy(i,j,k) = (lcs%lpY1%fm%z(i,j,k)-lcs%lpY0%fm%z(i,j,k))&
                               / (lcs%lpY1%sgrid%grid%y(i,j,k)-lcs%lpY0%sgrid%grid%y(i,j,k))
                        endif             
                        if( abs(lcs%lpZ1%sgrid%grid%z(i,j,k)-lcs%lpZ0%sgrid%grid%z(i,j,k))>SMALL) then
                              gradfm%xz(i,j,k) = (lcs%lpZ1%fm%x(i,j,k)-lcs%lpZ0%fm%x(i,j,k))&
                               / (lcs%lpZ1%sgrid%grid%z(i,j,k)-lcs%lpZ0%sgrid%grid%z(i,j,k))
                              gradfm%yz(i,j,k) = (lcs%lpZ1%fm%y(i,j,k)-lcs%lpZ0%fm%y(i,j,k))&
                               / (lcs%lpZ1%sgrid%grid%z(i,j,k)-lcs%lpZ0%sgrid%grid%z(i,j,k))
                              gradfm%zz(i,j,k) = (lcs%lpZ1%fm%z(i,j,k)-lcs%lpZ0%fm%z(i,j,k))&
                               / (lcs%lpZ1%sgrid%grid%z(i,j,k)-lcs%lpZ0%sgrid%grid%z(i,j,k))
                        endif             
                  enddo
                  enddo
                  enddo
            else
                  call exchange_sdata(lcs%sgrid%scomm_face_r1,r1=lcs%lp%fm)
                  if(lcs%sgrid%rectilinear) then
                        call grad_sr1(lcs%sgrid,lcs%lp%fm,gradfm)
                  else
                        call grad_sr1_ls(lcs%sgrid,lcs%lp%fm,gradfm)
                  endif
            endif

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
                        lcs%ftle%r(i,j,k) = log(max(sqrt(W(3)),tiny(1.0_LCSRP)))/lcs%t
                  case(8)
                        if(i == 1 .AND. j==1 .AND. k==1) then
                              LWORK= -1
                              call DSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
                              LWORK = MIN(LWMAX,INT(WORK(1)))
                        endif
                        call DSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
                        lcs%ftle%r(i,j,k) = log(max(sqrt(W(3)),tiny(1.0_LCSRP)))/lcs%t
                  end select
            enddo
            enddo
            enddo

            !Zero any negative values of FTLE:
            if(INCOMPRESSIBLE) then
                  lcs%ftle%r = max(lcs%ftle%r,0.0_LCSRP)
            endif

            call destroy_sr2(gradfm)
            call destroy_sr2(cg)
                  
            t1 = cputimer(lcscomm,SYNC_TIMER)
            cpu_ftle = cpu_ftle + max(t1-t0,0.0)

      end subroutine compute_ftle

      subroutine destroy_lcs(lcs)
            implicit none
            type(lcs_t),pointer :: lcs
            !-----
            lcs%id = -1
            lcs%label = 'LCS_NOT_USED'
            lcs%diagnostic = -1
            lcs%T = 0.0_LCSRP
            lcs%h = 0.0_LCSRP

            if(associated(lcs%sgrid))     nullify(lcs%sgrid)
            if(associated(lcs%lp)) nullify(lcs%lp)
            if(associated(lcs%lpX0)) nullify(lcs%lpX0)
            if(associated(lcs%lpY0)) nullify(lcs%lpY0)
            if(associated(lcs%lpZ0)) nullify(lcs%lpZ0)
            if(associated(lcs%lpX1)) nullify(lcs%lpX1)
            if(associated(lcs%lpY1)) nullify(lcs%lpY1)
            if(associated(lcs%lpZ1)) nullify(lcs%lpZ1)

            call destroy_sr0(lcs%ftle)

            if(associated(lcs)) nullify(lcs)
      end subroutine destroy_lcs
      
      subroutine create_aux_grid(lcs)
            implicit none
            !-----
            type(lcs_t):: lcs
            !-----
            integer:: ni,nj,nk
            integer:: nbr,i,j,k,ii,jj,kk,ind
            type(lp_t),pointer:: dummyX0,dummyY0,dummyZ0,dummyX1,dummyY1,dummyZ1,thisLP
            integer:: n(3),offset(3)
            real(LCSRP):: x(1:lcs%sgrid%ni,1:lcs%sgrid%nj,1:lcs%sgrid%nk)
            real(LCSRP):: y(1:lcs%sgrid%ni,1:lcs%sgrid%nj,1:lcs%sgrid%nk)
            real(LCSRP):: z(1:lcs%sgrid%ni,1:lcs%sgrid%nj,1:lcs%sgrid%nk)
            real(LCSRP):: aux_delta(1:lcs%sgrid%ni,1:lcs%sgrid%nj,1:lcs%sgrid%nk)
            integer(LCSIP):: flag(1:lcs%sgrid%ni,1:lcs%sgrid%nj,1:lcs%sgrid%nk)
            type(sgrid_t),pointer:: sgridX0,sgridX1,sgridY0,sgridY1,sgridZ0,sgridZ1
            real(LCSRP):: nbrsum,distsum
            real(LCSRP):: dx,dy,dz
            !-----
            !Create an auxillary grid of particles around the lcs%lp, to allow for increased
            !accuracy when computing the CG def tensor.
            !See "Computing Lagrangian coherent structures from their variational theory"
            !by Farazmand & Haller, Chaos 2012, for more detail.
            !-----
            
            !brevity:
            ni = lcs%sgrid%ni
            nj = lcs%sgrid%nj
            nk = lcs%sgrid%nk
            n = (/lcs%sgrid%ni,lcs%sgrid%nj,lcs%sgrid%nk/)
            offset = (/lcs%sgrid%offset_i,lcs%sgrid%offset_j,lcs%sgrid%offset_k/)

            !-----
            !Compute the auxillary grid offsets as a fraction of the mean distance to
            !each of the in-bounds neighbors:
            !-----
            aux_delta = 0.0_LCSRP
            do k = 1,nk
            do j = 1,nj
            do i = 1,ni
                  nbrsum = 0.0_LCSRP
                  distsum = 0.0_LCSRP
                  do nbr = 2,7  !use only faces (NBR 2-7)
                        ii = i+NBR_OFFSET(1,nbr)
                        jj = j+NBR_OFFSET(2,nbr)
                        kk = k+NBR_OFFSET(3,nbr)
                        if(ii>0 .and. jj>0 .and. kk>0 .and. ii<=ni .and. jj<=nj.and. kk<=nk) then
                              nbrsum = nbrsum+1.0_LCSRP
                              dx = lcs%sgrid%grid%x(i,j,k) - lcs%sgrid%grid%x(ii,jj,kk)
                              dy = lcs%sgrid%grid%y(i,j,k) - lcs%sgrid%grid%y(ii,jj,kk)
                              dz = lcs%sgrid%grid%z(i,j,k) - lcs%sgrid%grid%z(ii,jj,kk)
                              distsum = distsum + sqrt(dx*dx+dy*dy+dz*dz)
                        endif
                  enddo
                  aux_delta(i,j,k) = distsum/nbrsum*AUX_GRID_SCALING
            enddo
            enddo
            enddo

            !-----
            !Initialize some dummy lp structures based on the original lp:
            !-----
            call init_lp(dummyX0,lcs%lp%sgrid,lcs%lp%direction,lcs%lp%dt_factor,'DUMMYX0')
            call init_lp(dummyY0,lcs%lp%sgrid,lcs%lp%direction,lcs%lp%dt_factor,'DUMMYY0')
            call init_lp(dummyZ0,lcs%lp%sgrid,lcs%lp%direction,lcs%lp%dt_factor,'DUMMYZ0')
            call init_lp(dummyX1,lcs%lp%sgrid,lcs%lp%direction,lcs%lp%dt_factor,'DUMMYX1')
            call init_lp(dummyY1,lcs%lp%sgrid,lcs%lp%direction,lcs%lp%dt_factor,'DUMMYY1')
            call init_lp(dummyZ1,lcs%lp%sgrid,lcs%lp%direction,lcs%lp%dt_factor,'DUMMYZ1')
            
            !-----
            !Shift the particles in the dummy structures:
            !-----
            do k = 1,nk
            do j = 1,nj
            do i = 1,ni
                  ind = lcs_ijk2l(i,j,k,ni,nj)
                  dummyX0%xp%x(ind) = lcs%lp%xp%x(ind) - aux_delta(i,j,k)
                  dummyY0%xp%y(ind) = lcs%lp%xp%y(ind) - aux_delta(i,j,k)
                  dummyZ0%xp%z(ind) = lcs%lp%xp%z(ind) - aux_delta(i,j,k)
                  dummyX1%xp%x(ind) = lcs%lp%xp%x(ind) + aux_delta(i,j,k)
                  dummyY1%xp%y(ind) = lcs%lp%xp%y(ind) + aux_delta(i,j,k)
                  dummyZ1%xp%z(ind) = lcs%lp%xp%z(ind) + aux_delta(i,j,k)
            enddo
            enddo
            enddo

            !-----
            !Project these particles into the data plane (if 2d)
            !-----
            call project_lp2plane(dummyX0,lcs%lp%sgrid,dummyX0%no)
            call project_lp2plane(dummyY0,lcs%lp%sgrid,dummyY0%no)
            call project_lp2plane(dummyZ0,lcs%lp%sgrid,dummyZ0%no)
            call project_lp2plane(dummyX1,lcs%lp%sgrid,dummyX1%no)
            call project_lp2plane(dummyY1,lcs%lp%sgrid,dummyY1%no)
            call project_lp2plane(dummyZ1,lcs%lp%sgrid,dummyZ1%no)

            !-----
            !Now, create new sgrid based on these x,y,z points:
            !-----
            !X0
            thisLP => dummyX0
            do k = 1,nk
            do j = 1,nj
            do i = 1,ni
                  ind = lcs_ijk2l(i,j,k,ni,nj)
                  x(i,j,k) = thisLP%xp%x(ind) 
                  y(i,j,k) = thisLP%xp%y(ind) 
                  z(i,j,k) = thisLP%xp%z(ind) 
                  flag(i,j,k) = lcs%sgrid%bcflag%i(i,j,k)
            enddo
            enddo
            enddo
            call init_sgrid(sgridX0,trim(lcs%label)//'-AUXGRID-X0',n,offset,x,y,z,flag)
            !Y0
            thisLP => dummyY0
            do k = 1,nk
            do j = 1,nj
            do i = 1,ni
                  ind = lcs_ijk2l(i,j,k,ni,nj)
                  x(i,j,k) = thisLP%xp%x(ind) 
                  y(i,j,k) = thisLP%xp%y(ind) 
                  z(i,j,k) = thisLP%xp%z(ind) 
                  flag(i,j,k) = lcs%sgrid%bcflag%i(i,j,k)
            enddo
            enddo
            enddo
            call init_sgrid(sgridY0,trim(lcs%label)//'-AUXGRID-Y0',n,offset,x,y,z,flag)
            !Z0
            thisLP => dummyZ0
            do k = 1,nk
            do j = 1,nj
            do i = 1,ni
                  ind = lcs_ijk2l(i,j,k,ni,nj)
                  x(i,j,k) = thisLP%xp%x(ind) 
                  y(i,j,k) = thisLP%xp%y(ind) 
                  z(i,j,k) = thisLP%xp%z(ind) 
                  flag(i,j,k) = lcs%sgrid%bcflag%i(i,j,k)
            enddo
            enddo
            enddo
            call init_sgrid(sgridZ0,trim(lcs%label)//'-AUXGRID-Z0',n,offset,x,y,z,flag)
            !X1
            thisLP => dummyX1
            do k = 1,nk
            do j = 1,nj
            do i = 1,ni
                  ind = lcs_ijk2l(i,j,k,ni,nj)
                  x(i,j,k) = thisLP%xp%x(ind) 
                  y(i,j,k) = thisLP%xp%y(ind) 
                  z(i,j,k) = thisLP%xp%z(ind) 
                  flag(i,j,k) = lcs%sgrid%bcflag%i(i,j,k)
            enddo
            enddo
            enddo
            call init_sgrid(sgridX1,trim(lcs%label)//'-AUXGRID-X1',n,offset,x,y,z,flag)
            !Y1
            thisLP => dummyY1
            do k = 1,nk
            do j = 1,nj
            do i = 1,ni
                  ind = lcs_ijk2l(i,j,k,ni,nj)
                  x(i,j,k) = thisLP%xp%x(ind) 
                  y(i,j,k) = thisLP%xp%y(ind) 
                  z(i,j,k) = thisLP%xp%z(ind) 
                  flag(i,j,k) = lcs%sgrid%bcflag%i(i,j,k)
            enddo
            enddo
            enddo
            call init_sgrid(sgridY1,trim(lcs%label)//'-AUXGRID-Y1',n,offset,x,y,z,flag)
            !Z1
            thisLP => dummyZ1
            do k = 1,nk
            do j = 1,nj
            do i = 1,ni
                  ind = lcs_ijk2l(i,j,k,ni,nj)
                  x(i,j,k) = thisLP%xp%x(ind) 
                  y(i,j,k) = thisLP%xp%y(ind) 
                  z(i,j,k) = thisLP%xp%z(ind) 
                  flag(i,j,k) = lcs%sgrid%bcflag%i(i,j,k)
            enddo
            enddo
            enddo
            call init_sgrid(sgridZ1,trim(lcs%label)//'-AUXGRID-Z1',n,offset,x,y,z,flag)

            !-----
            !Initialize the particles on the 6 auxillary grids:
            !And, track the lp to the scfd grid
            !-----
            call init_lp(lcs%lpX0,sgridX0,lcs%lp%direction,lcs%lp%dt_factor,trim(lcs%label)//'-LP0')
            call init_lp(lcs%lpY0,sgridY0,lcs%lp%direction,lcs%lp%dt_factor,trim(lcs%label)//'-LP1')
            call init_lp(lcs%lpZ0,sgridZ0,lcs%lp%direction,lcs%lp%dt_factor,trim(lcs%label)//'-LP2')
            call init_lp(lcs%lpX1,sgridX1,lcs%lp%direction,lcs%lp%dt_factor,trim(lcs%label)//'-LP3')
            call init_lp(lcs%lpY1,sgridY1,lcs%lp%direction,lcs%lp%dt_factor,trim(lcs%label)//'-LP4')
            call init_lp(lcs%lpZ1,sgridZ1,lcs%lp%direction,lcs%lp%dt_factor,trim(lcs%label)//'-LP5')
            call track_lp2node(lcs%lpX0,scfd%sgrid,lcs%lpX0%no_scfd)
            call track_lp2node(lcs%lpY0,scfd%sgrid,lcs%lpY0%no_scfd)
            call track_lp2node(lcs%lpZ0,scfd%sgrid,lcs%lpZ0%no_scfd)
            call track_lp2node(lcs%lpX1,scfd%sgrid,lcs%lpX1%no_scfd)
            call track_lp2node(lcs%lpY1,scfd%sgrid,lcs%lpY1%no_scfd)
            call track_lp2node(lcs%lpZ1,scfd%sgrid,lcs%lpZ1%no_scfd)

            !-----
            !Destroy the dummy structures
            !-----
            call destroy_lp(dummyX0)
            call destroy_lp(dummyY0)
            call destroy_lp(dummyZ0)
            call destroy_lp(dummyX1)
            call destroy_lp(dummyY1)
            call destroy_lp(dummyZ1)

      end subroutine create_aux_grid

end module lcs_m
