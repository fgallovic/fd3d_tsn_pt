!-------------------------------------------------------
!   Main routine for the FD dynamic simulation
!-------------------------------------------------------
      subroutine fd3d()
      USE medium_com
      USE displt_com
      USE strfld_com
      USE fd3dparam_com
      USE friction_com
      USE source_com
      IMPLICIT NONE

      real    :: time,tabs,friction,tmax,xmax,ymax,numer,denom,veltest,dd
      integer :: incrack(nxt,nzt),broken(nxt,nzt),nsurf
      real    :: U1OUT,SLIPRATEOUT(nxt,nzt),gliss(nxt,nzt)
      REAL    :: CPUT1,CPUT2
      REAL :: maxvel,maxvelsave
      INTEGER i,j,it,k

!---------------------------
! Write down the input
!---------------------------
       if (ioutput.eq.1) then
         open(95,file='result/vmodel.inp')
         open(96,file='result/friction.inp')
         do k = 1,nzt
           do i = 1,nxt
             write(95,*) mu1(i,nysc,k)
             write(96,'(4E13.5)') strinix(i,k),peak_xz(i,k),Dc(i,k),peak_xz(i,k)/normstress(k)
            enddo
         enddo
         close(95)
         close(96)
       endif

!-------------------------------------------------------
!   initialize arrays
!-------------------------------------------------------
      u1=0.;v1=0.;w1=0.
      xx=0.;yy=0.;zz=0.;xy=0.;yz=0.;xz=0.
      ruptime=0.;slip=0.;rise=0.;schange=0.
      broken=0;incrack=0;gliss=0.
      sliprate=0.

!-------------------------------------------------------
!     Loop over time
!-------------------------------------------------------

      CALL CPU_TIME(CPUT1)
      maxvelsave=0.
      
      !$ACC DATA COPYIN (LAM1,MU1,D1) &
      !$ACC      COPYIN (U1,V1,W1) COPYIN (XX,YY,ZZ,XY,YZ,XZ) &
      !$ACC      COPYIN (BROKEN,gliss,STRINIX,PEAK_XZ,DC) &
      !$ACC      COPY (INCRACK,RUPTIME,SLIP,SCHANGE)

      do it = 1,ntfd
        time = it*dt

         SLIPRATEOUT(:,:)=0.
!$ACC DATA COPY (SLIPRATEOUT)
        
!-------------------------------------------------------------
!   Apply 4th-order differencing to interior particle velocities
!             from  2:nxt;2:nyt;2:nzt
!-------------------------------------------------------------
         call dvel(nxt,nysc,nzt,dt,dh)

!-------------------------------------------------------------
!   Compute absorbing boundary conditions in velocities
!-------------------------------------------------------------
         call abc(nxt,nysc,nzt,dt,dh)

!-------------------------------------------------------------
!   Compute velocities of 2nd order accuracy at boundary
!-------------------------------------------------------------
         call bnd2d(nxt,nysc,nzt,dt,dh)

!Stress free
         call fuvw(nxt,nyt,nzt-2)

!   Symetrie rychlosti podle plochy zlomu (y = nysc)
!----------------------------------------------------------------

        !$ACC PARALLEL DEFAULT (PRESENT)
        !$ACC LOOP GANG
         do i = 1,nxt
           !$ACC LOOP VECTOR
           do k = 1,nzt
             u1(i,nysc+1, k) = - u1(i,nysc-1+1, k)
             v1(i,nysc+1, k) = v1(i,nysc-1, k)
             w1(i,nysc+1, k) = - w1(i,nysc-1+1, k)
             u1(i,nysc+2, k) = - u1(i,nysc-2+1, k)
             v1(i,nysc+2, k) = v1(i,nysc-2, k)
             w1(i,nysc+2, k) = - w1(i,nysc-2+1, k)
           end do
         end do
        !$ACC END PARALLEL

!-------------------------------------------------------------
!   4th-order differencing of pressure
!-------------------------------------------------------------
         call dstres(nxt,nysc,nzt,dh,dt)
!-------------------------------------------------------------
!   Compute stress tensor of 2nd order accuracy at boundary
!-------------------------------------------------------------
         call strbnd(nxt,nysc,nzt,dh,dt)
         
!----------------------------------------------------------------
!   Symetrie napeti podle plochy zlomu (y = nysc)
!----------------------------------------------------------------

        !$ACC PARALLEL DEFAULT (PRESENT)
        !$ACC LOOP GANG
         do i = 1,nxt
           !$ACC LOOP VECTOR
           do k = 1,nzt
             xx(i,nysc+1, k) = - xx(i,nysc-1+1, k)
             yy(i,nysc+1, k) = - yy(i,nysc-1+1, k)
             zz(i,nysc+1, k) = - zz(i,nysc-1+1, k)
             xy(i,nysc+1, k) = xy(i,nysc-1, k)
             xz(i,nysc+1, k) = - xz( i,nysc-1+1, k)
             yz(i,nysc+1, k) = yz(i,nysc-1, k)
             xx(i,nysc+2, k) = - xx(i,nysc-2+1, k)
             yy(i,nysc+2, k) = - yy(i,nysc-2+1, k)
             zz(i,nysc+2, k) = - zz(i,nysc-2+1, k)
             xy(i,nysc+2, k) = xy(i,nysc-2, k)
             xz(i,nysc+2, k) = - xz( i,nysc-2+1, k)
             yz(i,nysc+2, k) = yz(i,nysc-2, k)
           end do
         end do
        !$ACC END PARALLEL


!----------------------------------------------------------------
! CRACK BOUNDARY CONDITIONS (applied right after stress update)
!----------------------------------------------------------------
!----------------------------------------------------------------
! The boundary conditions are applied now in ABSOLUTE STRESS,
! not in relative stress, otherwise we get all sort of problems
!----------------------------------------------------------------
        !$ACC PARALLEL DEFAULT (PRESENT)
        !$ACC LOOP GANG
        do k = 3,nzt-2
          !$ACC LOOP VECTOR
          do i = 3,nxt-2
            !-------------------------------
            ! Only for non-strike-slip fault
            !-------------------------------
            !tabsx = xy(i,nysc,k) + strinix(i,k)
            !tabsz = yz(i,nysc,k) + striniy(i,k)
            !tabs  = sqrt(tabsx**2+tabsy**2)

#if defined DIPSLIP
              tabs = yz(i,nysc,k) + strinix(i,k)
#else
              tabs = xy(i,nysc,k) + strinix(i,k)
#endif

            if (tabs.gt.peak_xz(i,k)) then
              broken(i,k)  = 1
              incrack(i,k) = 1
              if (ruptime(i,k).eq.0.) ruptime(i,k) = time
            endif
          enddo
        enddo
        !$ACC END PARALLEL

!-------------------------------------------------------------
!  Apply Boundary Conditions
!-------------------------------------------------------------
        !$ACC PARALLEL DEFAULT (PRESENT)
        !$ACC LOOP GANG
        do k = 3,nzt-2
          !$ACC LOOP VECTOR
          do i = 3,nxt-2
!-------------------------------------------------------------
!   Do boundary conditions only for points that are considered
!   broken, i.e. those that are slipping
!-------------------------------------------------------------
            if (incrack(i,k).eq.1) then
              if (broken(i,k).eq.1) then
#if defined DIPSLIP
                  veltest=w1(i,nysc+1,k)
#else
                  veltest=u1(i,nysc+1,k)
#endif

                if (veltest.lt.0.) then
                  broken(i,k)=0;gliss(i,k)=0.          !INSTANTANOUS HEALING
                  
#if defined DIPSLIP
                  w1(i,nysc+1,k) = 0.
                  w1(i,nysc,k) = 0.
!                  w1(i,nysc-1,k) = 0.
#else
                  u1(i,nysc+1,k) = 0.
                  u1(i,nysc,k) = 0.
!                  u1(i,nysc-1,k) = 0.
#endif
                  v1(i,nysc,k) = 0.
                else
!----------------------------------------------------------------------------
!   Apply friction law to absolute slip, absolute sliprate, absolute traction
!   Classical slip weakening law
!----------------------------------------------------------------------------
                  if (Dc(i,k)<=0.) then
                    dd = 0.
                  else
                    if (gliss(i,k).le.Dc(i,k)) then
                      dd = 1.0 - gliss(i,k)/Dc(i,k)
                    else
                      dd = 0.0
                    endif
                  endif
                  friction     = peak_xz(i,k)*dd
#if defined DIPSLIP
                  yz(i,nysc,k) = friction - strinix(i,k)  ! THIN B.C.
                  U1OUT=W1(I,NYSC+1,K)
#else
!                  if(xy(i,nysc,k)+ strinix(i,k) > friction)then !VERZE S PODMINKOU NA PREKROCENI NAPETI
!                    xy(i,nysc,k) = friction - strinix(i,k)  ! THIN B.C.
!                  else
!                    u1(i,nysc+1,k) = 0.
!                    u1(i,nysc,k) = 0.
!                    broken(i,k)=0;gliss(i,k)=0.
!                  endif
!                  U1OUT=U1(I,NYSC+1,K)

                  xy(i,nysc,k) = friction - strinix(i,k)  ! THIN B.C.
                  U1OUT=U1(I,NYSC+1,K)
#endif
                  slip(i,k) = slip(i,k) + 2.0*u1out*dt
                  gliss(i,k) = gliss(i,k) + 2.0*u1out*dt
                  sliprateout(i,k) = 2.*u1out
                endif
              endif
            endif
#if defined DIPSLIP
          SCHANGE(I,K)=YZ(I,NYSC,K)
#else
          SCHANGE(I,K)=XY(I,NYSC,K)
#endif
          enddo
        enddo
        !$ACC END PARALLEL

!Stress free
         call fres(nxt,nyt,nzt-2)
         
        !$ACC END DATA
        sliprate(:,:,it)=sliprateout(:,:)
        if(mod(it,int(1./dt))==0)then
          maxvel=maxval(sliprateout(:,:))
          write(*,*)'Checkuji',time,maxvel
          if(maxvel>maxvelsave)maxvelsave=maxvel
          if (maxvel<=0.01*maxvelsave)exit
        endif
        
      enddo ! --- End of the time loop

      !$ACC END DATA

      CALL CPU_TIME(CPUT2)
      PRINT *,'CPU TIME OF TIME LOOP: ',CPUT2-CPUT1

      
!-------------------
! Open output files:
!-------------------
      if(ioutput.eq.1) then
!       OPEN(23, file='result/ssx3d.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
!       OPEN(24, file='result/disp.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
        OPEN(25, file='result/sliprate.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
!        WRITE(25) sliprate(:,:,:)
        do it = 1,ntfd
!          do k = 1,nzt
!            do i = 1,nxt
!               write(23,'(1g12.4)') schange(i,k)
!               write(24,'(1g12.4)') slip(i,k)
!               write(25,'(1g12.4)') sliprate(i,k,it)
!               WRITE(23) schange(I,K)
!               WRITE(24) slip(I,K)
!               WRITE(25) sliprate(I,K,it)
               WRITE(25) sliprate(1:nxt,1:nzt,it)
!            enddo
!          enddo
        enddo
!        do k = 1,nzt
!          do i = 1,nxt
!             write(23,'(1g12.4)') schange(i,k)
!             write(24,'(1g12.4)') slip(i,k)
!             WRITE(23) schange(I,K)
!             WRITE(24) slip(I,K)
!          enddo
!        enddo
         close(25)
       endif
      
      do it=1,ntfd           !Saving slip rate at a point on the fault
        time = (it-1)*dt
        write(388,*)time,sliprate(nxt/3,nzt/2,it)
      enddo

      tmax            = -1.
      output_param(1) =  0.
      nsurf           =  0
      output_param(2) =  0.
      numer           =  0.
      denom           =  0.
      do k = 1,nzt
        do i = 1,nxt
          ! --- Average rupture speed:
          if (ruptime(i,k).gt.tmax) then
            xmax            = dh/2. + (i-1)*dh
            ymax            = dh/2. + (nzt-(k-1))*dh
!            dist            = sqrt((xmax-estk)**2 + (ymax-edip)**2) !nefunguje, nezna hypocentrum
            tmax            = ruptime(i,k)
!            output_param(1) = dist/tmax !nefunguje, nezna hypocentrum
          endif

          ! --- Surface of rupture:
          if (incrack(i,k).eq.1) nsurf = nsurf + 1
          output_param(3) = nsurf * (dh*dh)

          ! --- Seismic moment:
          output_param(2) = output_param(2) + slip(i,k)*mu1(i,nysc+1,k)*(dh*dh)

           ! --- Rise time:
           if (maxval(sliprate(i,k,:)).eq.0) then
             rise(i,k) = 0.
           elseif (ruptime(i,k).ne.0.) then
             rise(i,k) = slip(i,k)/maxval(sliprate(i,k,:))
           else
             rise(i,k) = 0.
           endif

           ! --- Stress drop:
           if (ruptime(i,k).ne.0.) then
             numer = numer + schange(i,k)*slip(i,k)
             denom = denom + slip(i,k)
           endif
           output_param(4) = -(numer/denom)
         enddo
       enddo
       output_param(5) = (1./2.)**sum(peak_xz*Dc)/dble(nxt*nzt)
       output_param(6) = (1./2.)*output_param(4)*(output_param(2)/(mu_mean*output_param(3)))

!---------------------------
! Write down the output
!---------------------------
       if (ioutput.eq.1) then
         open(96,file='result/risetime.res')
         open(97,file='result/ruptime.res')
         open(98,file='result/slip.res')
         open(99,file='result/stressdrop.res')
         do k = 1,nzt
           do i = 1,nxt
             write(96,*) rise(i,k)
             write(97,*) ruptime(i,k)
             write(98,*) slip(i,k)
             write(99,*) schange(i,k)
           enddo
         enddo
         close(96)
         close(97)
         close(98)
         close(99)
       endif

      return
      end
