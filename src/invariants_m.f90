      module invariants_m
      use data_m
      use comms_m
      use structured_m
      use gradient_m
      implicit none

      real(LCSRP), parameter:: eigth = 1./8.
      real(LCSRP), parameter:: twenyforth = 1./24.
      real(LCSRP), parameter:: coef = 27./4.
      integer, parameter:: LWMAX = 1000

      contains

      subroutine compute_invariants(lcs)
         type(lcs_t):: lcs
         !-----
         type(sr2_t):: gradu
         integer:: ni,nj,nk,ng
         integer:: i,j,k
         real(LCSRP) :: d1u1,d1u2,d1u3,d2u1,d2u2,d2u3,d3u1,d3u2,d3u3
         real(LCSRP) :: S11,S12,S13,S21,S22,S23,S31,S32,S33
         real(LCSRP) :: W11,W12,W13,W21,W22,W23,W31,W32,W33
         real(LCSRP) :: INV(5), D
         real(LCSRP):: A(3,3), W(3), WORK(LWMAX)
         integer:: INFO,LWORK
         !-----

         !Allocate space
         ni = lcs%sgrid%ni
         nj = lcs%sgrid%nj
         nk = lcs%sgrid%nk
         ng = lcs%sgrid%ng
         call init_sr2(gradu,ni,nj,nk,ng,'GradU')
         
         !Compute Grad(u)
         call exchange_sdata(lcs%sgrid%scomm_face_r1,r1=lcs%lp%ugrid)
         if(lcs%sgrid%rectilinear) then
            call grad_sr1(lcs%sgrid,lcs%lp%ugrid,gradu)
         else
            call grad_sr1_ls(lcs%sgrid,lcs%lp%ugrid,gradu)
         endif
      
         
         !Loop over grid points and compute invariants
         do k = 1,nk
         do j = 1,nj
         do i = 1,ni
            d1u1 = gradu%xx(i,j,k)
            d2u1 = gradu%xy(i,j,k) ! d2u1 is du1/dx2 
            d3u1 = gradu%xz(i,j,k)
            d1u2 = gradu%yx(i,j,k)
            d2u2 = gradu%yy(i,j,k)
            d3u2 = gradu%yz(i,j,k)
            d1u3 = gradu%zx(i,j,k)
            d2u3 = gradu%zy(i,j,k)
            d3u3 = gradu%zz(i,j,k)

            S11 = d1u1 + d1u1
            S12 = d1u2 + d2u1
            S13 = d1u3 + d3u1
            S21 = d2u1 + d1u2
            S22 = d2u2 + d2u2
            S23 = d2u3 + d3u2
            S31 = d3u1 + d1u3
            S32 = d3u2 + d2u3
            S33 = d3u3 + d3u3
            
            W11 = d1u1 - d1u1
            W12 = d1u2 - d2u1
            W13 = d1u3 - d3u1
            W21 = d2u1 - d1u2
            W22 = d2u2 - d2u2
            W23 = d2u3 - d3u2
            W31 = d3u1 - d1u3
            W32 = d3u2 - d2u3
            W33 = d3u3 - d3u3

            !!!!!!! QS
            INV(1)=-eigth*(                     &
                             S11*S11 + S12*S21 + S13*S31 &
                           + S21*S12 + S22*S22 + S23*S32 &
                           + S31*S13 + S32*S23 + S33*S33 &
                                  )
            !!!!!!! RS   
            INV(2)=-twenyforth*(                 &
                  S11*S11*S11 + S11*S12*S21 + S11*S13*S31 &
                + S21*S11*S12 + S21*S12*S22 + S21*S13*S32 &
                + S31*S11*S13 + S31*S12*S23 + S31*S13*S33 &
                + S12*S21*S11 + S12*S22*S21 + S12*S23*S31 &
                + S13*S31*S11 + S13*S32*S21 + S13*S33*S31 &
                + S22*S22*S22 + S22*S21*S12 + S22*S23*S32 &
                + S23*S31*S13 + S23*S32*S22 + S23*S33*S32 &
                + S32*S22*S23 + S32*S21*S13 + S32*S23*S33 &
                + S33*S33*S33 + S33*S23*S32 + S33*S13*S31 &
                                        )
            !!!!!!! Q
            INV(3)= INV(1) - eigth*( &
                        W11*W11 + W12*W21 + W13*W31      &
                      + W21*W12 + W22*W22 + W23*W32      &
                      + W31*W13 + W32*W23 + W33*W33      &
                                                       )
            !!!!!! R
            INV(4)= INV(2) - eigth*(    &
                  W11*W11*S11 + W11*W12*S21 + W11*W13*S31 &
                + W21*W11*S12 + W21*W12*S22 + W21*W13*S32 &
                + W31*W11*S13 + W31*W12*S23 + W31*W13*S33 &
                + W12*W21*S11 + W12*W22*S21 + W12*W23*S31 &
                + W13*W31*S11 + W13*W32*S21 + W13*W33*S31 &
                + W22*W22*S22 + W22*W21*S12 + W22*W23*S32 &
                + W23*W31*S13 + W23*W32*S22 + W23*W33*S32 &
                + W32*W22*S23 + W32*W21*S13 + W32*W23*S33 &
                + W33*W33*S33 + W33*W23*S32 + W33*W13*S31 &
                                                     )
            !!!!!! Qw=Q-Qs
            INV(5)=  INV(3) - INV(1)

            ! the Q-criterium: The invariant Q quantifies the balance 
            ! between dissipation and enstrophy, hence Q>0 is a vortex &
            ! and Q<0 is a straining node
            IF (INV(3)>0) THEN
              lcs%inv%Q%r(i,j,k) = 1.0_LCSRP
            ELSE
              lcs%inv%Q%r(i,j,k) = 0.0_LCSRP 
            ENDIF

            !!!!! In the joint pdf (Q,R), the curve D=Q^3+27/4*R^2 
            !separates regions of real eigenvalues from region of complexe
            !eigenvalues i.e. straining nodes (D<0) from vortices (D>0). 
            !Besides invariant R gives information on whether we face a 
            !stretching structure or a compression structure. Thus, we cut 
            !(Q,R) in 4 zones that define whether the surrounding flow 
            !structure is a Vortex compression (lcs%inv%D%r=1) or Vortex stretching
            !(=2) or Sheet structure (=3) or Tube structure (=4). The D-criterai
            !from Jeong and Hussain 95 is recovered by only considering
            !the criteria D<0 or >0
            D=INV(3)**3 + coef*INV(4)**2
            IF (D.gt.0.0.AND.INV(4).gt.0.0) THEN
              lcs%inv%D%r(i,j,k) = 1.0_LCSRP
            ELSEIF (D.gt.0.0.AND.INV(4).le.0.0) THEN
              lcs%inv%D%r(i,j,k) = 2.0_LCSRP
            ELSEIF (D.le.0.0.AND.INV(4).gt.0.0) THEN
              lcs%inv%D%r(i,j,k) = 3.0_LCSRP
            ELSEIF (D.le.0.0.AND.INV(4).le.0.0) THEN
              lcs%inv%D%r(i,j,k) = 4.0_LCSRP
            ENDIF

            ! The invariant RS informs on whether we face a contraction 
            !or expansion
            IF (INV(2).gt.0.0) THEN
               lcs%inv%C%r(i,j,k) = 1.0_LCSRP
            ELSE
               lcs%inv%C%r(i,j,k) = 0.0_LCSRP
            ENDIF

            !!! vorticity number from Truesdell 54 (see Ooi J. Fluid Mech. (1999)
            ! vol. 381, pp. 141–174) telles whether we face Irrotational
            !dissipation (large H), Vortex tube (small H) and Vortex sheets (H=1)
            lcs%inv%H%r(i,j,k) = sqrt(-INV(5)/INV(1))

            !!!!! Eigenvalues, we only keep the information whether lambda2 of S^2+OMEGA^2 is negative (vortex) or not
            !         symmetric tensor S^2+OMEGA^2 (see Jeong & Hussain 95 JFM), we only compute the upper part
            !                     S11 S12 S13
            !                     S21 S22 S23
            !                     S31 S32 S33
            !        S11 S12 S13  ...
            !        S21 S22 S23
            !        S31 S32 S33
            A(1,1) = S11**2   + S21*S12 + S31*S13 + W11**2   + W21*W12 + W31*W13
            A(1,2) = S12*S11 + S22*S12 + S32*S13 + W12*W11 + W22*W12 + W32*W13
            A(1,3) = S13*S11 + S23*S12 + S33*S13 + W13*W11 + W23*W12 + W33*W13
            A(2,1) = 0.
            A(2,2) = S12*S21 + S22**2   + S32*S23 + W12*W21 + W22**2   + W32*W23
            A(2,3) = S13*S21 + S23*S22 + S33*S23 + W13*W21 + W23*W22 + W33*W23
            A(3,1) = 0.
            A(3,2) = 0.
            A(3,3) = 0.
            select case(LCSRP)
               case(4)
                  if(i == 1 .AND. j==1 .AND. k==1) then
                     LWORK= -1
                     call SSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
                     LWORK = MIN(LWMAX,INT(WORK(1)))
                  endif
                  call SSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
               case(8)
                     if(i == 1 .AND. j==1 .AND. k==1) then
                     LWORK= -1
                     call DSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
                     LWORK = MIN(LWMAX,INT(WORK(1)))
                  endif
                  call DSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
            end select
            call sort3(W)
            IF(W(2)<0) then
               lcs%inv%L2%r(i,j,k) = 1.0_LCSRP
            ELSE
               lcs%inv%L2%r(i,j,k) = 0.0_LCSRP
            ENDIF

         enddo
         enddo
         enddo

         !Cleanup
         call destroy_sr2(gradu)
      
      end subroutine compute_invariants

      subroutine sort3(x)
            implicit none
            !------
            real(LCSRP), intent(inout):: x(3)
            real(LCSRP):: a,b,c,tmp
            !------
            
            !Sort the three real 'X' values in ascending order
            
            a = x(1); b=x(2); c=x(3);
            if (b < a) then
                  tmp = b
                  b = a
                  a = tmp
            endif
            if (c < b) then
                  tmp = c 
                  c = b
                  b = tmp
            endif
            if (b < a) then
                  tmp = b
                  b = a
                  a = tmp
            endif
            x(1)=a; x(2)=b; x(3)=c;

      end subroutine sort3
      
      subroutine sort3pair(x,y)
            implicit none
            !------
            real(LCSRP), intent(inout):: x(3),y(3)
            real(LCSRP):: a,b,c,tmp
            real(LCSRP):: ay,by,cy,tmpy
            !------
            
            !Sort the three real 'X' values in ascending order
            !Also, reorder the 'Y' values accordingly
            
            a=x(1); b=x(2); c=x(3);
            ay=y(1); by=y(2); cy=y(3);
            if (b < a) then
                  tmp = b
                  tmpy = by
                  b = a
                  by = ay
                  a = tmp
                  ay = tmpy
            endif
            if (c < b) then
                  tmp = c 
                  tmpy = cy 
                  c = b
                  cy = by
                  b = tmp
                  by = tmpy
            endif
            if (b < a) then
                  tmp = b
                  tmpy = by
                  b = a
                  by = ay
                  a = tmp
                  ay = tmpy
            endif
            x(1)=a; x(2)=b; x(3)=c;
            y(1)=ay; y(2)=by; y(3)=cy;
      end subroutine sort3pair

      end module invariants_m

