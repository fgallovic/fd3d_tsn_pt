      MODULE displt_com
        real,allocatable,dimension(:,:,:):: u1,v1,w1
      END MODULE

      MODULE strfld_com
        real,allocatable,dimension(:,:,:):: xx,yy,zz,xy,yz,xz
      END MODULE

!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine dvel(nxt,nyt,nzt,dt,dh)
!----------------------------------------------------------
!     4th order finite-difference of velocity components
!     nxt   nodal points in x dir  (integer)(sent)
!     nyt   nodal points in y dir  (integer)(sent)
!     nzt   nodal points in z dir  (integer)(sent)
!     dh    spatial discretization (real)   (sent)
!     dt    temporal discretization(real)   (sent)
!----------------------------------------------------------
      integer :: nxt,nyt,nzt
      real    :: dh,dt
!----------------------------------------------------------
!     Find displacement fields at time t+1/2
!----------------------------------------------------------
      !call uxx1(3,nxt-2,3,nyt-2,3,nzt-2,dh,dt)
      !call vyy1(3,nxt-3,3,nyt-3,3,nzt-2,dh,dt)
      !call wzz1(3,nxt-3,3,nyt-2,3,nzt-3,dh,dt)
      !ve smeru y se pocita az ke zlomu  (y=nysc)
      call uxx1(3,nxt-2,3,nyt,3,nzt-2,dh,dt)
      call vyy1(3,nxt-3,3,nyt,3,nzt-2,dh,dt)
      call wzz1(3,nxt-3,3,nyt,3,nzt-2,dh,dt)

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine uxx1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     4nd order finite-difference of u1
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------

      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dh,dt,dth,c1,c2,d

      dth = dt/dh
      c1  = 9./8.
      c2  = -1./24.
!----------------------------------------------------------
!     Find u-displacement fields at time t+1/2
!----------------------------------------------------------
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe

        d         = d1(i,j,k)
        u1(i,j,k) = u1(i,j,k) + (dth/d)*(    &
           c1*(xx(i,j,k)   - xx(i-1,j,k)) +  &
           c2*(xx(i+1,j,k) - xx(i-2,j,k)) +  &
           c1*(xy(i,j,k)   - xy(i,j-1,k)) +  &
           c2*(xy(i,j+1,k) - xy(i,j-2,k)) +  &
           c1*(xz(i,j,k)   - xz(i,j,k-1)) +  &
           c2*(xz(i,j,k+1) - xz(i,j,k-2)))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine vyy1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     4nd order finite-difference of v1
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dh,dt,dth,c1,c2,d

      dth = dt/dh
      c1  = 9./8.
      c2  = -1./24.
!----------------------------------------------------------
!     Find v-displacement fields at time t+1/2
!----------------------------------------------------------
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe

        d         = d1(i,j,k)
        v1(i,j,k) = v1(i,j,k) + (dth/d)*(    &
           c1*(xy(i+1,j,k) - xy(i,j,k))   +  &
           c2*(xy(i+2,j,k) - xy(i-1,j,k)) +  &
           c1*(yy(i,j+1,k) - yy(i,j,k))   +  &
           c2*(yy(i,j+2,k) - yy(i,j-1,k)) +  &
           c1*(yz(i,j,k)   - yz(i,j,k-1)) +  &
           c2*(yz(i,j,k+1) - yz(i,j,k-2)))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine wzz1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     4nd order finite-difference of w1
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dh,dt,dth,c1,c2,d

      dth = dt/dh
      c1  = 9./8.
      c2  = -1./24.
!----------------------------------------------------------
!     Find w-displacement fields at time t+1/2
!----------------------------------------------------------
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe

        d         = d1(i,j,k)
        w1(i,j,k) = w1(i,j,k) + (dth/d)*(    &
           c1*(xz(i+1,j,k) - xz(i,j,k))   +  &
           c2*(xz(i+2,j,k) - xz(i-1,j,k)) +  &
           c1*(yz(i,j,k)   - yz(i,j-1,k)) +  &
           c2*(yz(i,j+1,k) - yz(i,j-2,k)) +  &
           c1*(zz(i,j,k+1) - zz(i,j,k))   +  &
           c2*(zz(i,j,k+2) - zz(i,j,k-1)))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine abc(nxt,nyt,nzt,dt,dh)
!----------------------------------------------------------
!     absorbing boundary condition
!     nxt   nodal points in x dir  (integer)(sent)
!     nyt   nodal points in y dir  (integer)(sent)
!     nzt   nodal points in z dir  (integer)(sent)
!     dh    spatial discretization (real)   (sent)
!     dt    temporal discretization(real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxt,nyt,nzt
      real    :: dh,dt,dth,d,xl,xm,dp1,dp,ds1,ds
      integer :: nxm,nxm2,nym,nym2,nzm,nzm2

      nxm  = nxt - 1
      nxm2 = nxt - 2
      nym  = nyt - 1
      nym2 = nyt - 2

      nzm  = nzt - 1
      nzm2 = nzt - 2
      dth  = dt/dh

!----------------------------------------------------------
!     Apply abc at i=1 and i=nxt
!----------------------------------------------------------
!     Composente u1

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do k = 1,nzt
      !$ACC LOOP VECTOR
      do j = 1,nyt

        d   = d1(1,j,k)
        xl  = lam1(1,j,k)
        xm  = mu1(1,j,k)
        dp1 = xl + 2.*xm
        dp  = dth*sqrt(dp1/d)
        u1(1,j,k) = u1(1,j,k) + dp*(u1(2,j,k)-u1(1,j,k))
        d   = d1(nxt,j,k)
        xl  = lam1(nxt,j,k)
        xm  = mu1(nxt,j,k)
        dp1 = xl + 2.*xm
        dp  = dth*sqrt(dp1/d)
        u1(nxt,j,k) = u1(nxt,j,k) - dp*(u1(nxt,j,k)-u1(nxm,j,k))

      enddo
      enddo
      !$ACC END PARALLEL

!     Composente v1

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do k = 1,nzt
      !$ACC LOOP VECTOR
      !do j = 1,nym
      do j = 1,nyt
        d   = d1(1,j,k)
        xm  = mu1(1,j,k)
        ds1 = xm/d
        ds  = dth*sqrt(ds1)
        v1(1,j,k) = v1(1,j,k) + ds*(v1(2,j,k)-v1(1,j,k))
        d   = d1(nxm,j,k)
        xm  = mu1(nxm,j,k)
        ds1 = xm/d
        ds  = dth*sqrt(ds1)
        v1(nxm,j,k) = v1(nxm,j,k) - ds*(v1(nxm,j,k)-v1(nxm2,j,k))

      enddo
      enddo
      !$ACC END PARALLEL

!     Composente w1

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do k = 1,nzm
      !$ACC LOOP VECTOR
      do j = 1,nyt

        d   = d1(1,j,k)
        xm  = mu1(1,j,k)
        ds1 = xm/d
        ds  = dth*sqrt(ds1)
        w1(1,j,k) = w1(1,j,k) + ds*(w1(2,j,k)-w1(1,j,k))
        d   = d1(nxm,j,k)
        xm  = mu1(nxm,j,k)
        ds1 = xm/d
        ds  = dth*sqrt(ds1)
        w1(nxm,j,k) = w1(nxm,j,k) - ds*(w1(nxm,j,k)-w1(nxm2,j,k))

      enddo
      enddo
      !$ACC END PARALLEL

!----------------------------------------------------------
!     Apply abc at j=1 and j=nyt
!----------------------------------------------------------
!     Composente u1

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do k = 1,nzt
      !$ACC LOOP VECTOR
      do i = 2,nxm

        d   = d1(i,1,k)
        xm  = mu1(i,1,k)
        ds1 = xm/d
        ds  = dth*sqrt(ds1)
        u1(i,1,k) = u1(i,1,k) + ds*(u1(i,2,k)-u1(i,1,k))
        !d   = d1(i,nyt,k)
        !xm  = mu1(i,nyt,k)
        !ds1 = xm/d
        !ds  = dth*sqrt(ds1)
        !u1(i,nyt,k) = u1(i,nyt,k) - ds*(u1(i,nyt,k)-u1(i,nym,k))

      enddo
      enddo
      !$ACC END PARALLEL

!     Composente v1

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do k = 1,nzt
      !$ACC LOOP VECTOR
      do i = 2,nxm2

        d   = d1(i,1,k)
        xl  = lam1(i,1,k)
        xm  = mu1(i,1,k)
        dp1 = xl+2.*xm
        dp  = dth*sqrt(dp1/d)
        v1(i,1,k) = v1(i,1,k) + dp*(v1(i,2,k)-v1(i,1,k))
        !d   = d1(i,nym,k)      !Vypocet je symetricky v y
        !xl  = lam1(i,nym,k)
        !xm  = mu1(i,nym,k)
        !dp1 = xl+2.*xm
        !dp  = dth*sqrt(dp1/d)
        !v1(i,nym,k) = v1(i,nym,k) - dp*(v1(i,nym,k)-v1(i,nym2,k))

      enddo
      enddo
      !$ACC END PARALLEL

!     Composente w1

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do k = 1,nzm
      !$ACC LOOP VECTOR
      do i = 2,nxm2

        d   = d1(i,1,k)
        xm  = mu1(i,1,k)
        ds1 = xm/d
        ds  = dth*sqrt(ds1)
        w1(i,1,k) = w1(i,1,k) + ds*(w1(i,2,k)-w1(i,1,k))
        !d   = d1(i,nyt,k)  !Vypocet je symetricky v y
        !xm  = mu1(i,nyt,k)
        !ds1 = xm/d
        !ds  = dth*sqrt(ds1)
        !w1(i,nyt,k) = w1(i,nyt,k) - ds*(w1(i,nyt,k)-w1(i,nym,k))

      enddo
      enddo
      !$ACC END PARALLEL
!----------------------------------------------------------
!     Apply abc at z=1 and z=nzt
!----------------------------------------------------------
!     Composente u1

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      !do j = 2,nym
      do j = 2,nyt
      !$ACC LOOP VECTOR
      do i = 2,nxm

        d   = d1(i,j,1)
        xm  = mu1(i,j,1)
        ds1 = xm/d
        ds  = dth*sqrt(ds1)
        u1(i,j,1) = u1(i,j,1) + ds*(u1(i,j,2)-u1(i,j,1))
!        d   = d1(i,j,nzt)
!        xm  = mu1(i,j,nzt)
!        ds1 = xm/d
!        ds  = dth*sqrt(ds1)
!        u1(i,j,nzt) = u1(i,j,nzt) - ds*(u1(i,j,nzt)-u1(i,j,nzm)) !proc to tady je? nad volnym povrchem (nzt-2) nula

      enddo
      enddo
      !$ACC END PARALLEL

!     Composente v1

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      !do j = 2,nym2
      do j = 2,nyt
      !$ACC LOOP VECTOR
      do i = 2,nxm2

        d   = d1(i,j,1)
        xm  = mu1(i,j,1)
        ds1 = xm/d
        ds  = dth*sqrt(ds1)
        v1(i,j,1) = v1(i,j,1) + ds*(v1(i,j,2)-v1(i,j,1))
!        d   = d1(i,j,nzt)
!        xm  = mu1(i,j,nzt)
!        ds1 = xm/d
!        ds  = dth*sqrt(ds1)
!        v1(i,j,nzt) = v1(i,j,nzt) - ds*(v1(i,j,nzt)-v1(i,j,nzm))!proc to tady je? nad volnym povrchem (nzt-2) nula

      enddo
      enddo
      !$ACC END PARALLEL

!     Composente w1

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      !do j = 2,nym
      do j = 2,nyt
      !$ACC LOOP VECTOR
      do i = 2,nxm2

        d   = d1(i,j,1)
        xl  = lam1(i,j,1)
        xm  = mu1(i,j,1)
        dp1 = xl+2.*xm
        dp  = dth*sqrt(dp1/d)
        w1(i,j,1) = w1(i,j,1) + dp*(w1(i,j,2)-w1(i,j,1))
!        d   = d1(i,j,nzt)
!        xl  = lam1(i,j,nzt)
!        xm  = mu1(i,j,nzt)
!        dp1 = xl+2.*xm
!        dp  = dth*sqrt(dp1/d)
!        w1(i,j,nzm) = w1(i,j,nzm) - dp*(w1(i,j,nzm)-w1(i,j,nzm2))!proc to tady je? nad volnym povrchem (nzt-2) nula

      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine bnd2d(nxt,nyt,nzt,dt,dh)
!----------------------------------------------------------
!     bnd2d finds 2nd-order differencing of wave eq at bnd
!     to obtain the velocity values
!     nxt   nodal points in x dir          (integer)(sent)
!     nyt   nodal points in y dir          (integer)(sent)
!     nzt   nodal points in z dir          (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      integer :: nxt,nyt,nzt
      real    :: dh,dt
!----------------------------------------------------------
!     Find displacement fields at time t+1/2 at x=2 and x=nx-1 by 2nd
!     order differences
!----------------------------------------------------------
      !call uxx0(2,2,2,nyt-1,2,nzt-1,dh,dt)
      call uxx0(2,2,2,nyt,2,nzt-2,dh,dt)
      !call uxx0(nxt-1,nxt-1,2,nyt-1,2,nzt-1,dh,dt)
      call uxx0(nxt-1,nxt-1,2,nyt,2,nzt-2,dh,dt)

      !call vyy0(2,2,2,nyt-2,2,nzt-1 ,dh,dt)
      call vyy0(2,2,2,nyt,2,nzt-2 ,dh,dt)

      !call vyy0(nxt-2,nxt-2,2,nyt-2,2,nzt-1,dh,dt)
      call vyy0(nxt-2,nxt-2,2,nyt,2,nzt-2,dh,dt)

      !call wzz0(2,2,2,nyt-1,2,nzt-2,dh,dt)
      call wzz0(2,2,2,nyt,2,nzt-2,dh,dt)
      !call wzz0(nxt-2,nxt-2,2,nyt-1,2,nzt-2,dh,dt)
      call wzz0(nxt-2,nxt-2,2,nyt,2,nzt-2,dh,dt)
!----------------------------------------------------------
!     Find displacement fields at time t+1/2 at y=2 and y=ny-1
!----------------------------------------------------------
      call uxx0(3,nxt-2,2,2,2,nzt-2,dh,dt)
      !call uxx0(3,nxt-2,nyt-1,nyt-1,2,nzt-1,dh,dt)

      call vyy0(3,nxt-3,2,2,2,nzt-2,dh,dt)
      !call vyy0(3,nxt-3,nyt-2,nyt-2,2,nzt-1,dh,dt)

      call wzz0(3,nxt-3,2,2,2,nzt-2,dh,dt)
      !call wzz0(3,nxt-3,nyt-1,nyt-1,2,nzt-2,dh,dt)
!----------------------------------------------------------
!     Find displacement fields at time t+1/2 at z=2
!----------------------------------------------------------
      !call uxx0(3,nxt-2,3,nyt-2,2,2,dh,dt)
      call uxx0(3,nxt-2,3,nyt,2,2,dh,dt)

      !call vyy0(3,nxt-3,3,nyt-3,2,2,dh,dt)
      call vyy0(3,nxt-3,3,nyt,2,2,dh,dt)

      !call wzz0(3,nxt-3,3,nyt-2,2,2,dh,dt)
      call wzz0(3,nxt-3,3,nyt,2,2,dh,dt)

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine uxx0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     2nd order finite-difference of u1
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dt,dh,dth,d

      dth = dt/dh
!----------------------------------------------------------
!     Find u-displacement fields at time t+1/2
!----------------------------------------------------------
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe

        d         = d1(i,j,k)
        u1(i,j,k) = u1(i,j,k) + (dth/d)*(  &
           (xx(i,j,k) - xx(i-1,j,k)) +     &
           (xy(i,j,k) - xy(i,j-1,k)) +     &
           (xz(i,j,k) - xz(i,j,k-1)))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine vyy0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     2nd order finite-difference of v1
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dt,dh,dth,d

      dth = dt/dh
!----------------------------------------------------------
!     Find v-displacement fields at time t+1/2
!----------------------------------------------------------
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe

        d         = d1(i,j,k)
        v1(i,j,k) = v1(i,j,k) + (dth/d)*(  &
           (xy(i+1,j,k) - xy(i,j,k)) +     &
           (yy(i,j+1,k) - yy(i,j,k)) +     &
           (yz(i,j,k)   - yz(i,j,k-1)))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine wzz0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     2nd order finite-difference of w1
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dt,dh,dth,d

      dth = dt/dh
!----------------------------------------------------------
!     Find w-displacement fields at time t+1/2
!----------------------------------------------------------
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe

        d         = d1(i,j,k)
        w1(i,j,k) = w1(i,j,k) + (dth/d)*(  &
          (xz(i+1,j,k) - xz(i,j,k))   +    &
          (yz(i,j,k)   - yz(i,j-1,k)) +    &
          (zz(i,j,k+1) - zz(i,j,k)))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine dstres(nxt,nyt,nzt,dh,dt)
!----------------------------------------------------------
!     4th order finite-difference of stress components
!     nxt   nodal points in x dir          (integer)(sent)
!     nyt   nodal points in y dir          (integer)(sent)
!     nzt   nodal points in z dir          (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxt,nyt,nzt
      real    :: dh,dt,dth,c1,c2,xl,xm,a,b,xm1,xm2,xmu

      dth = dt/dh
      c1  = 9./8.
      c2  = -1./24.
!----------------------------------------------------------
!     Find displacement fields at time t+1/2
!----------------------------------------------------------
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = 3,(nzt-2)
      !do j = 3,(nyt-2)
      !do j = 3,(nyt-2)
      do j = 3,nyt
      !$ACC LOOP VECTOR
      do i = 2,(nxt-2)

        xl = lam1(i,j,k)
        xm = mu1(i,j,k)
        a  = xl + 2.*xm
        b  = xl

!       Find xx stress

        xx(i,j,k) = xx(i,j,k)                      +  &
            dth*a*(c1*(u1(i+1,j,k) - u1(i,j,k))    +  &
                   c2*(u1(i+2,j,k) - u1(i-1,j,k))) +  &
            dth*b*(c1*(v1(i,j,k)   - v1(i,j-1,k))  +  &
                   c2*(v1(i,j+1,k) - v1(i,j-2,k))  +  &
                   c1*(w1(i,j,k)   - w1(i,j,k-1))  +  &
                   c2*(w1(i,j,k+1) - w1(i,j,k-2)))

!     Find yy stress

        yy(i,j,k) = yy(i,j,k)                      +  &
            dth*a*(c1*(v1(i,j,k)   - v1(i,j-1,k))  +  &
                   c2*(v1(i,j+1,k) - v1(i,j-2,k))) +  &
            dth*b*(c1*(u1(i+1,j,k) - u1(i,j,k))    +  &
                   c2*(u1(i+2,j,k) - u1(i-1,j,k))  +  &
                   c1*(w1(i,j,k)   - w1(i,j,k-1))  +  &
                   c2*(w1(i,j,k+1) - w1(i,j,k-2)))

!     Find zz stress

         zz(i,j,k) = zz(i,j,k)                     +  &
            dth*a*(c1*(w1(i,j,k)   - w1(i,j,k-1))  +  &
                   c2*(w1(i,j,k+1) - w1(i,j,k-2))) +  &
            dth*b*(c1*(u1(i+1,j,k) - u1(i,j,k))    +  &
                   c2*(u1(i+2,j,k) - u1(i-1,j,k))  +  &
                   c1*(v1(i,j,k)   - v1(i,j-1,k))  +  &
                   c2*(v1(i,j+1,k) - v1(i,j-2,k)))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

!----------------------------------------------------------
!     shear stresses
!----------------------------------------------------------

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = 2,(nzt-2)
      !do j = 2,(nyt-2)
      do j = 2,nyt
      !$ACC LOOP VECTOR
      do i = 3,(nxt-2)

        dth = .5*dt/dh
        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j+1,k)
        xmu = xm1+xm2

!       Find xy stress

        xy(i,j,k) = xy(i,j,k)                      +  &
          dth*xmu*(c1*(u1(i,j+1,k) - u1(i,j,k))    +  &
                   c2*(u1(i,j+2,k) - u1(i,j-1,k))  +  &
                   c1*(v1(i,j,k)   - v1(i-1,j,k))  +  &
                   c2*(v1(i+1,j,k) - v1(i-2,j,k)))

        dth = .5*dt/dh
        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j,k+1)
        xmu = xm1+xm2

!       Find xz stress

        xz(i,j,k) = xz(i,j,k)                      +  &
          dth*xmu*(c1*(u1(i,j,k+1) - u1(i,j,k))    +  &
                   c2*(u1(i,j,k+2) - u1(i,j,k-1))  +  &
                   c1*(w1(i,j,k)   - w1(i-1,j,k))  +  &
                   c2*(w1(i+1,j,k) - w1(i-2,j,k)))

        dth = .5*dt/dh
        xm1 = mu1(i,j,k)
        xm2 = mu1(i+1,j+1,k+1)
        xmu = xm1+xm2

!       Find yz stress

        yz(i,j,k) = yz(i,j,k)                      +  &
          dth*xmu*(c1*(v1(i,j,k+1) - v1(i,j,k))    +  &
                   c2*(v1(i,j,k+2) - v1(i,j,k-1))  +  &
                   c1*(w1(i,j+1,k) - w1(i,j,k))    +  &
                   c2*(w1(i,j+2,k) - w1(i,j-1,k)))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      !call sxy1(3,nxt-2,2,nyt-2,1,1,dh,dt)
      !call sxy1(3,nxt-2,2,nyt-2,nzt-1,nzt,dh,dt)

      call sxy1(3,nxt-2,2,nyt,1,1,dh,dt)
      !call sxy1(3,nxt-2,2,nyt,nzt-1,nzt,dh,dt) !proc to tady je? nad volnym povrchem (nzt-2) nula

      !call sxz1(3,nxt-2,1,1,2,nzt-2,dh,dt)
      !call sxz1(3,nxt-2,nyt-1,nyt,2,nzt-2,dh,dt)

      call sxz1(3,nxt-2,1,1,2,nzt-2,dh,dt)
      !call sxz1(3,nxt-2,nyt-1,nyt,2,nzt-2,dh,dt)

      !call syz1(1,2,2,nyt-2,2,nzt-2,dh,dt)
      !call syz1(nxt-1,nxt-1,2,nyt-2,2,nzt-2,dh,dt)
      call syz1(1,2,2,nyt,2,nzt-2,dh,dt)
      call syz1(nxt-1,nxt-1,2,nyt,2,nzt-2,dh,dt)

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine sxy1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     4th order finite-difference of xy
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dh,dt,dth,c1,c2,xm1,xm2,xmu

      dth = .5*dt/dh
      c1  = 9./8.
      c2  = -1./24.

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe
!----------------------------------------------------------
!     Find xy stress
!----------------------------------------------------------
         xm1 = mu1(i,j,k)
         xm2 = mu1(i,j+1,k)
         xmu = xm1+xm2
         xy(i,j,k) = xy(i,j,k)                      +  &
           dth*xmu*(c1*(u1(i,j+1,k) - u1(i,j,k))    +  &
                    c2*(u1(i,j+2,k) - u1(i,j-1,k))  +  &
                    c1*(v1(i,j,k)   - v1(i-1,j,k))  +  &
                    c2*(v1(i+1,j,k) - v1(i-2,j,k)))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine sxz1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     4th order finite-difference of xz
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dh,dt,dth,c1,c2,xm1,xm2,xmu

      dth = .5*dt/dh
      c1  = 9./8.
      c2  = -1./24.

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe
!----------------------------------------------------------
!     Find xz stress
!----------------------------------------------------------
        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j,k+1)
        xmu = xm1+xm2
        xz(i,j,k) = xz(i,j,k)                      +  &
          dth*xmu*(c1*(u1(i,j,k+1) - u1(i,j,k))    +  &
                   c2*(u1(i,j,k+2) - u1(i,j,k-1))  +  &
                   c1*(w1(i,j,k)   - w1(i-1,j,k))  +  &
                   c2*(w1(i+1,j,k) - w1(i-2,j,k)))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine syz1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     4th order finite-difference of yz
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dh,dt,dth,c1,c2,xm1,xm2,xmu

      dth = .5*dt/dh
      c1  = 9./8.
      c2  = -1./24.

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe
!----------------------------------------------------------
!     Find yz stress
!----------------------------------------------------------
        xm1 = mu1(i,j,k)
        xm2 = mu1(i+1,j+1,k+1)
        xmu = xm1+xm2
        yz(i,j,k) = yz(i,j,k)                        +  &
            dth*xmu*(c1*(v1(i,j,k+1) - v1(i,j,k))    +  &
                     c2*(v1(i,j,k+2) - v1(i,j,k-1))  +  &
                     c1*(w1(i,j+1,k) - w1(i,j,k))    +  &
                     c2*(w1(i,j+2,k) - w1(i,j-1,k)))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine strbnd(nxt,nyt,nzt,dh,dt)
!----------------------------------------------------------
!     2th order finite-difference of stresses at boundaries
!     nxt   nodal points in x dir  (integer)(sent)
!     nyt   nodal points in y dir  (integer)(sent)
!     nzt   nodal points in z dir  (integer)(sent)
!     dh    spatial discretization (real)   (sent)
!     dt    temporal discretization(real)   (sent)
!----------------------------------------------------------
      integer :: nxt,nyt,nzt
      real    :: dh,dt
!----------------------------------------------------------
!     compute stresses at bnd with 2nd order difference operator
!----------------------------------------------------------
      call sxx(1,1,2,nyt,2,nzt-2,dh,dt)
      call sxx(nxt-1,nxt-1,2,nyt,2,nzt-2,dh,dt)

      call sxx(2,nxt-2,2,2,2,nzt-2,dh,dt)
      !call sxx(2,nxt-2,nyt-1,nyt,2,nzt,dh,dt)

      !call sxx(2,nxt-2,3,nyt-2,2,2,dh,dt)
      !call sxx(2,nxt-2,3,nyt-2,nzt-1,nzt,dh,dt)
      call sxx(2,nxt-2,3,nyt,2,2,dh,dt)

      !call sxx(2,nxt-2,3,nyt,nzt-1,nzt,dh,dt)
!----------------------------------------------------------
!     xy stress on x=2 plane and x=nxt-1
!----------------------------------------------------------
      !call sxy0(2,2,1,nyt-1,1,nzt,dh,dt)
      !call sxy0(nxt-1,nxt,1,nyt-1,1,nzt,dh,dt)
      call sxy0(2,2,1,nyt,1,nzt-2,dh,dt)
      call sxy0(nxt-1,nxt,1,nyt,1,nzt-2,dh,dt)
!----------------------------------------------------------
!     xy stress on y=2 plane and y=nyt-1
!----------------------------------------------------------
      call sxy0(3,nxt-2,1,1,1,nzt-2,dh,dt)
      !call sxy0(3,nxt-2,nyt-1,nyt-1,1,nzt,dh,dt)
!----------------------------------------------------------
!     xz stress on z=1, x=2, x=nxt-1 planes
!----------------------------------------------------------
      call sxz0(2,nxt,1,nyt,1,1,dh,dt)
      !call sxz0(2,nxt,1,nyt,nzt-1,nzt-1,dh,dt)
      call sxz0(2,2,1,nyt,2,nzt-2,dh,dt)
      call sxz0(nxt-1,nxt,1,nyt,2,nzt-2,dh,dt)
!----------------------------------------------------------
!     yz stress on z=1, y=1, y=nyt-1 planes
!----------------------------------------------------------
      !call syz0(1,nxt-1,1,nyt-1,1,1,dh,dt)
      !call syz0(1,nxt-1,1,nyt-1,nzt-1,nzt-1,dh,dt)
      call syz0(1,nxt-1,1,nyt,1,1,dh,dt)
      !call syz0(1,nxt-1,1,nyt,nzt-1,nzt-1,dh,dt)
      call syz0(1,nxt-1,1,1,2,nzt-2,dh,dt)

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine sxx(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     STRESS FREE BOUNDARY CONDITION AT k=nzt
!----------------------------------------------------------
!----------------------------------------------------------
!     2th order finite-difference of normal stresses
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dh,dt,dth,xl,xm,a,b

      dth = dt/dh
!----------------------------------------------------------
!     Find displacement fields at time t+1/2
!----------------------------------------------------------
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe

        xl = lam1(i,j,k)
        xm = mu1(i,j,k)
        a  = xl + 2.*xm
        b  = xl

!       Find xx stress

        xx(i,j,k) = xx(i,j,k)               +  &
          dth*a*(u1(i+1,j,k) - u1(i,j,k))   +  &
          dth*b*(v1(i,j,k)   - v1(i,j-1,k)  +  &
                 w1(i,j,k)   - w1(i,j,k-1))

!       Find yy stress

        yy(i,j,k) = yy(i,j,k)               +  &
          dth*a*(v1(i,j,k)   - v1(i,j-1,k)) +  &
          dth*b*(u1(i+1,j,k) - u1(i,j,k)    +  &
                 w1(i,j,k)   - w1(i,j,k-1))

!       Find zz stress

        zz(i,j,k) = zz(i,j,k)               +  &
          dth*a*(w1(i,j,k)   - w1(i,j,k-1)) +  &
          dth*b*(u1(i+1,j,k) - u1(i,j,k)    +  &
                 v1(i,j,k)   - v1(i,j-1,k))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine sxy0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     2th order finite-difference of xy
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dh,dt,dth,xm1,xm2,xmu

      dth = .5*dt/dh

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe

!       Find xy stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j+1,k)
        xmu = xm1+xm2

        xy(i,j,k) = xy(i,j,k) + dth*xmu*(u1(i,j+1,k) - u1(i,j,k) + v1(i,j,k) - v1(i-1,j,k))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine sxz0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     2th order finite-difference of xz
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dh,dt,dth,xm1,xm2,xmu

      dth = .5*dt/dh

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k=nzb,nze
      do j=nyb,nye
      !$ACC LOOP VECTOR
      do i=nxb,nxe

!       Find xz stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j,k+1)
        xmu = xm1+xm2
        xz(i,j,k) = xz(i,j,k) + dth*xmu*(u1(i,j,k+1) - u1(i,j,k) + w1(i,j,k) - w1(i-1,j,k))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine syz0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
!----------------------------------------------------------
!     2th order finite-difference of yz
!     nxb   starting point for FD in x dir (integer)(sent)
!     nxe   ending point for FD in x dir   (integer)(sent)
!     nyb   starting point for FD in y dir (integer)(sent)
!     nye   ending point for FD in y dir   (integer)(sent)
!     nzb   starting point for FD in z dir (integer)(sent)
!     nze   ending point for FD in z dir   (integer)(sent)
!     dh    spatial discretization         (real)   (sent)
!     dt    temporal discretization        (real)   (sent)
!----------------------------------------------------------
      USE medium_com
      USE displt_com
      USE strfld_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dh,dt,dth,xm1,xm2,xmu

      dth = .5*dt/dh

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k=nzb,nze
      do j=nyb,nye
      !$ACC LOOP VECTOR
      do i=nxb,nxe

!        Find yz stress

         xm1 = mu1(i,j,k)
         xm2 = mu1(i+1,j+1,k+1)
         xmu = xm1+xm2

         yz(i,j,k) = yz(i,j,k) + dth*xmu*(v1(i,j,k+1) - v1(i,j,k) + w1(i,j+1,k) - w1(i,j,k))

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
    end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------

     subroutine fuvw(nxt,nyt,nzt)
      USE medium_com
      USE displt_com

!     free-surface B.C. for velocities

!     nxt   nodal points in x dir (integer)(sent)
!     nyt   nodal points in y dir (integer)(sent)
!     nzt   nodal points in z dir (integer)(sent)

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do j=1,nyt
      !$ACC LOOP VECTOR
         do i=2,nxt-1
            u1(i,j,nzt+1)=u1(i,j,nzt)-(w1(i,j,nzt)-w1(i-1,j,nzt))
         enddo
      enddo
      !$ACC END PARALLEL
      
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do j=2,nyt-2
      !$ACC LOOP VECTOR
         do i=2,nxt-1
            v1(i,j,nzt+1)=v1(i,j,nzt)-(w1(i,j+1,nzt)-w1(i,j,nzt))
         enddo
      enddo
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do j=2,nyt-2
      !$ACC LOOP VECTOR
         do i=2,nxt-2
            xl=lam1(i,j,nzt+1)
            xm=mu1(i,j,nzt+1)
            a=2.*xl
            b=xl+2.*xm
            w1(i,j,nzt+1)=w1(i,j,nzt)-(a/b)*(u1(i+1,j,nzt+1)-u1(i,j,nzt+1)+v1(i,j,nzt+1)-v1(i,j-1,nzt+1))
         enddo
      enddo
      !$ACC END PARALLEL

    end
    
    
    subroutine fres(nxt,nyt,nzt)
      USE medium_com
      USE displt_com
      USE strfld_com
!     free-surface B.C. for stresses

!     nxt   nodal points in x dir  (integer)(sent)
!     nyt   nodal points in y dir  (integer)(sent)
!     nzt   nodal points in z dir  (integer)(sent)
!     dh    spatial discretization (real)   (sent)
!     dt    temporal discretization(real)   (sent)

      nyf = nyt - 2
      nxf = nxt - 2
      nzm = nzt-1

      nzt2=nzt+2
      nzt1=nzt+1
      nztm=nzt-1
      nztm2=nzt-2

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do j=1,nyt
      !$ACC LOOP VECTOR
        do i=1,nxt

          zz(i,j,nzt1) = -zz(i,j,nzt)
          zz(i,j,nzt2) = -zz(i,j,nztm)
          xz(i,j,nzt1) = -xz(i,j,nztm)
          xz(i,j,nzt2) = -xz(i,j,nztm2)
          yz(i,j,nzt1) = -yz(i,j,nztm)
          yz(i,j,nzt2) = -yz(i,j,nztm2)
!     zero yz and xz at free-surface.
          xz(i,j,nzt) = 0.
          yz(i,j,nzt) = 0.

        enddo
      enddo
      !$ACC END PARALLEL
      
    end 
    