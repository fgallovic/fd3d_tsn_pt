      subroutine dvel(nxt,nyt,nzt,dt,dh)
         USE pml_com
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

      call uxx1(nabc+1,nxt-nabc,nabc+1,nyt,nabc+1,nzt-nfs,dh,dt)
      call vyy1(nabc+1,nxt-nabc,nabc+1,nyt,nabc+1,nzt-nfs,dh,dt)
      call wzz1(nabc+1,nxt-nabc,nabc+1,nyt,nabc+1,nzt-nfs,dh,dt)

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
      subroutine bnd2d(nxt,nyt,nzt,dt,dh)
         USE pml_com
         USE traction_com
         USE medium_com
         USE displt_com
         USE strfld_com
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
!     Find displacement fields at time t+1/2 at y=nysc-1 2nd
!     order differences
!----------------------------------------------------------

      call uxx0a (nabc+1, nxt-nabc, nyt+1, nyt+1, nabc+1, nzt-nfs,dh, dt)

      call vyy0a (nabc+1, nxt-nabc, nyt+1, nyt+1, nabc+1, nzt-nfs,dh, dt)

      call wzz0a (nabc+1, nxt-nabc, nyt+1, nyt+1, nabc+1, nzt-nfs,dh, dt)
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
      subroutine uxx0a(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
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
      USE traction_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dt,dh,dth,d, pdx

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
        pdx = (xx(i,j,k) - xx(i-1,j,k)) +     &
              (xy(i,j,k) - xy(i,j-1,k)) +     &
              (xz(i,j,k) - xz(i,j,k-1))
        au1(i,k) = pdx - au1(i,k)
        d         = d1(i,j,k)
        u1(i,j,k) = u1(i,j,k) + (dth/d)*(pdx + damp_s*au1(i,k) )
        au1(i,k) = pdx
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
      subroutine vyy0a(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
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
      USE traction_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dt,dh,dth,d, pdy

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
        pdy = (xy(i+1,j,k) - xy(i,j,k)) +     &
              (yy(i,j+1,k) - yy(i,j,k)) +     &
              (yz(i,j,k)   - yz(i,j,k-1))
        d         = d1(i,j,k)

        av1(i,k) = pdy - av1(i,k)
        v1(i,j,k) = v1(i,j,k) + (dth/d)*( pdy + damp_s*av1(i,k))
        av1(i,k) = pdy
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
      subroutine wzz0a(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
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
      USE traction_com

      integer :: nxb,nxe,nyb,nye,nzb,nze
      real    :: dt,dh,dth,d, pdz

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
        pdz = (xz(i+1,j,k) - xz(i,j,k))   +    &
              (yz(i,j,k)   - yz(i,j-1,k)) +    &
              (zz(i,j,k+1) - zz(i,j,k))
        d         = d1(i,j,k)
        aw1(i,k)= pdz - aw1(i,k)
        w1(i,j,k) = w1(i,j,k) + (dth/d)*( pdz + damp_s*aw1(i,k))
        aw1(i,k)= pdz
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
      USE pml_com
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
      do k = nabc+1,(nzt-nfs)
      do j = nabc+1,nyt
      !$ACC LOOP VECTOR
      do i = nabc+1,(nxt-nabc)

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
!----------------------------------------------------------
!     shear stresses
!----------------------------------------------------------

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j+1,k)
        xmu = 0.5*( xm1+xm2)



!       Find xy stress

        xy(i,j,k) = xy(i,j,k)                      +  &
          dth*xmu*(c1*(u1(i,j+1,k) - u1(i,j,k))    +  &
                   c2*(u1(i,j+2,k) - u1(i,j-1,k))  +  &
                   c1*(v1(i,j,k)   - v1(i-1,j,k))  +  &
                   c2*(v1(i+1,j,k) - v1(i-2,j,k)))

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j,k+1)
        xmu = 0.5*( xm1+xm2)

!       Find xz stress

        xz(i,j,k) = xz(i,j,k)                      +  &
          dth*xmu*(c1*(u1(i,j,k+1) - u1(i,j,k))    +  &
                   c2*(u1(i,j,k+2) - u1(i,j,k-1))  +  &
                   c1*(w1(i,j,k)   - w1(i-1,j,k))  +  &
                   c2*(w1(i+1,j,k) - w1(i-2,j,k)))

        xm1 = mu1(i,j,k)
        xm2 = mu1(i+1,j+1,k+1)
        xmu = 0.5*( xm1+xm2)

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
         USE pml_com
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
!     compute all stresses at y = nysc-1 with 2nd order difference operator
!----------------------------------------------------------

      call sxx (nabc+1,nxt-nabc, nyt+1, nyt+1, nabc+1,nzt-nfs,dh,dt)
      call sxy0(nabc+1,nxt-nabc, nyt+1, nyt+1, nabc+1,nzt-nfs,dh,dt)
      call sxz0(nabc+1,nxt-nabc, nyt+1, nyt+1, nabc+1,nzt-nfs,dh,dt)
      call syz0(nabc+1,nxt-nabc, nyt+1, nyt+1, nabc+1,nzt-nfs,dh,dt)

      return
      end
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
      subroutine sxx(nxb,nxe,nyb,nye,nzb,nze,dh,dt)
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
      do j=2,nyt-1		!mozna nyt 
      !$ACC LOOP VECTOR
         do i=2,nxt-1
            v1(i,j,nzt+1)=v1(i,j,nzt)-(w1(i,j+1,nzt)-w1(i,j,nzt))
         enddo
      enddo
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do j=2,nyt !mozna nyt 
      !$ACC LOOP VECTOR
         do i=2,nxt-2
            xl=lam1(i,j,nzt+1)
            xm=mu1(i,j,nzt+1)
            a=2.*xl
            b=xl+2.*xm
            w1(i,j,nzt+1)=w1(i,j,nzt-1)-(a/b)*(u1(i+1,j,nzt)-u1(i,j,nzt)+u1(i+1,j,nzt+1)-u1(i,j,nzt+1) &
                          +v1(i,j,nzt)-v1(i,j-1,nzt)+v1(i,j,nzt+1)-v1(i,j-1,nzt+1))
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

    subroutine tasu1(nxt,nyt,nzt,dt,dh)

      USE medium_com
      USE displt_com
      USE strfld_com
      USE traction_com
      USE pml_com
!u component at the fault

!     nxt   nodal points in x dir  (integer)(sent)
!     nyt   position of fault  (integer)(sent)
!     nzt   nodal points in z dir  (integer)(sent)
!     dh    spatial discretization (real)   (sent)
!     dt    temporal discretization(real)   (sent)

      real    ::  dh, dt, d, dth
      integer :: nxt, nyt, nzt

      dth = dt/dh

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do k = nabc+1,nzt-nfs
      !$ACC LOOP VECTOR
          do i = nabc+1,nxt-nabc
            d         = d1(i,nyt,k)
            u1(i,nyt,k) = u1(i,nyt,k) + (dth/d)*(2*(tx(i,k) + RFx(i,k)))
          enddo
      enddo
      !$ACC END PARALLEL

    end

    subroutine tasw1(nxt,nyt,nzt,dt,dh)

      USE medium_com
      USE displt_com
      USE strfld_com
      USE traction_com
      USE pml_com
!w component at the fault

!     nxt   nodal points in x dir  (integer)(sent)
!     nyt   position of fault  (integer)(sent)
!     nzt   nodal points in z dir  (integer)(sent)
!     dh    spatial discretization (real)   (sent)
!     dt    temporal discretization(real)   (sent)
      real    ::  dh, dt, d, dth
      integer :: nxt, nyt, nzt

      dth = dt/dh
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do k = nabc+1,nzt-nfs
      !$ACC LOOP VECTOR
          do i = nabc+1,nxt-nabc
            d         = d1(i,nyt,k)
            w1(i,nyt,k) = w1(i,nyt,k) + (dth/d)*(2*(tz(i,k) + RFz(i,k)))
          enddo
      enddo
      !$ACC END PARALLEL

    end

    subroutine tasxz(nxt,nyt,nzt,dt,dh)

      USE medium_com
      USE displt_com
      USE strfld_com
      USE traction_com
      USE pml_com
!xz component at the fault

!     nxt   nodal points in x dir  (integer)(sent)
!     nyt   position of fault  (integer)(sent)
!     nzt   nodal points in z dir  (integer)(sent)
!     dh    spatial discretization (real)   (sent)
!     dt    temporal discretization(real)   (sent)
      real dh, dt
      integer :: nxt, nyt, nzt

      call sxz1(nabc+1,nxt-nabc, nyt, nyt, nabc+1,nzt-nfs,dh,dt)

    end

    subroutine tasii(nxt,nyt,nzt,dt,dh)

      USE medium_com
      USE displt_com
      USE strfld_com
      USE traction_com
      USE pml_com
!w component at the fault

!     nxt   nodal points in x dir  (integer)(sent)
!     nyt   position of fault  (integer)(sent)
!     nzt   nodal points in z dir  (integer)(sent)
!     dh    spatial discretization (real)   (sent)
!     dt    temporal discretization(real)   (sent)
      real    ::  dh, dt, d, dth, diff1, diff2,c1,c2
      integer :: nxt, nyt, nzt
      c1  = 9./8.
      c2  = -1./24.
      dth = dt/dh
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
      do k = nabc+1,nzt-nfs
      !$ACC LOOP VECTOR
          do i = nabc+1,nxt-nabc
            xl = lam1(i,nyt,k)
            xm = mu1(i,nyt,k)
            a  = xl + 2.*xm
            b  = xl

	    diff1=c1*(u1(i+1,nyt,k) - u1(i,nyt,k)) + c2*(u1(i+2,nyt,k) - u1(i-1,nyt,k)) 
            diff3=c1*(w1(i,nyt,k)   - w1(i,nyt,k-1)) + c2*(w1(i,nyt,k+1) - w1(i,nyt,k-2))
		!diff1=(u1(i+1,nyt,k) - u1(i,nyt,k))
		!diff3=(w1(i,nyt,k)   - w1(i,nyt,k-1))
		
		
            v1t(i,k)=v1(i,nyt-1,k) - b*(diff1 + diff3)/(2*a)

            xx(i,nyt,k) = xx(i,nyt,k)             +  &
            dth*a*(diff1)   +  &
            dth*b*(2*(v1t(i,k)   - v1(i,nyt-1,k)) +  &
            diff3)

            yy(i,nyt,k) = yy(i,nyt,k)             +  &
            dth*a*2*(v1t(i,k)   - v1(i,nyt-1,k))  +  &
            dth*b*(diff1    +  &
            diff3)

            zz(i,nyt,k) = zz(i,nyt,k)             +  &
            dth*a*(diff3) +  &
            dth*b*(diff1    +  &
            2*(v1t(i,k)   - v1(i,nyt-1,k)))



          enddo
      enddo
      !$ACC END PARALLEL

    end

     subroutine uxxa(nxb,nxe,nyb,nye,nzb,nze,dh,dt, omegax, omegay, omegaz, p1,p2,p3)
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
      USE pml_com

      integer :: nxb,nxe,nyb,nye,nzb,nze,i2,j2,k2
      real    :: dt,dh,dth,d, pt
      real    :: omegax(nxe-nxb+1),  omegay(nye-nyb+1), omegaz(nze-nzb+1)
      real    :: p1(nxe-nxb+1,nye-nyb+1,nze-nzb+1),p2(nxe-nxb+1,nye-nyb+1,nze-nzb+1),p3(nxe-nxb+1,nye-nyb+1,nze-nzb+1)
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
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k

        pt = p1(i2,j2,k2)*(1.-omegax(i2))  + (dt/d)*(xx(i,j,k) - xx(i-1,j,k))/dh
        p1(i2,j2,k2) = pt/(1.+omegax(i2))
		
        pt = p2(i2,j2,k2)*(1.-omegay(j2))  + (dt/d)*(xy(i,j,k) - xy(i,j-1,k))/dh
        p2(i2,j2,k2) = pt/(1.+omegay(j2))
		
        pt = p3(i2,j2,k2)*(1.-omegaz(k2))  + (dt/d)*(xz(i,j,k) - xz(i,j,k-1))/dh
        p3(i2,j2,k2) = pt/(1.+omegaz(k2))

        u1(i,j,k) = p1(i2,j2,k2) + p2(i2,j2,k2) + p3(i2,j2,k2)
		
		! pt = p1(1-nxb+i,1-nyb+j,1-nzb+k)*(1.-omegax(1-nxb+i))  + (dt/d)*(xx(i,j,k) - xx(i-1,j,k))/dh
        ! p1(1-nxb+i,1-nyb+j,1-nzb+k) = pt/(1.+omegax(1-nxb+i))
		
        ! pt = p2(1-nxb+i,1-nyb+j,1-nzb+k)*(1.-omegay(1-nyb+j))  + (dt/d)*(xy(i,j,k) - xy(i,j-1,k))/dh
        ! p2(1-nxb+i,1-nyb+j,1-nzb+k) = pt/(1.+omegay(1-nyb+j))
		
        ! pt = p3(1-nxb+i,1-nyb+j,1-nzb+k)*(1.-omegaz(1-nzb+k))  + (dt/d)*(xz(i,j,k) - xz(i,j,k-1))/dh
        ! p3(1-nxb+i,1-nyb+j,1-nzb+k) = pt/(1.+omegaz(1-nzb+k))

        ! u1(i,j,k) = p1(1-nxb+i,1-nyb+j,1-nzb+k) + p2(1-nxb+i,1-nyb+j,1-nzb+k) + p3(1-nxb+i,1-nyb+j,1-nzb+k)

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
      subroutine vyya(nxb,nxe,nyb,nye,nzb,nze,dh,dt, omegax, omegay, omegaz, p1, p2, p3)
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
      USE pml_com
      integer :: nxb,nxe,nyb,nye,nzb,nze,i2,j2,k2
      real    :: dt,dh,dth,d,pt, cl
      real    :: omegax(nxe-nxb+1),  omegay(nye-nyb+1), omegaz(nze-nzb+1)
      real    :: p1(nxe-nxb+1,nye-nyb+1,nze-nzb+1),p2(nxe-nxb+1,nye-nyb+1,nze-nzb+1),p3(nxe-nxb+1,nye-nyb+1,nze-nzb+1)

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
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
		
        pt = p1(i2,j2,k2)*(1.-omegax(i2))  + (dt/d)*(xy(i+1,j,k) - xy(i,j,k))/dh
        p1(i2,j2,k2) = pt/(1.+omegax(i2))

        pt = p2(i2,j2,k2)*(1.-omegay(j2))  + (dt/d)*(yy(i,j+1,k) - yy(i,j,k))/dh
        p2(i2,j2,k2) = pt/(1.+omegay(j2))

        pt = p3(i2,j2,k2)*(1.-omegaz(k2))  + (dt/d)*(yz(i,j,k) - yz(i,j,k-1))/dh
        p3(i2,j2,k2) = pt/(1.+omegaz(k2))

        v1(i,j,k) = p1(i2,j2,k2) + p2(i2,j2,k2) + p3(i2,j2,k2)
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
      subroutine wzza(nxb,nxe,nyb,nye,nzb,nze,dh,dt, omegax, omegay, omegaz,p1,p2,p3)
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
      USE pml_com
      integer :: nxb,nxe,nyb,nye,nzb,nze,i2,j2,k2
      real    :: dt,dh,dth,d,pt
      real    :: omegax(nxe-nxb+1),  omegay(nye-nyb+1), omegaz(nze-nzb+1)
      real    :: p1(nxe-nxb+1,nye-nyb+1,nze-nzb+1),p2(nxe-nxb+1,nye-nyb+1,nze-nzb+1),p3(nxe-nxb+1,nye-nyb+1,nze-nzb+1)
     
!----------------------------------------------------------
!     Find w-displacement fields at time t+1/2
!----------------------------------------------------------
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
        d         = d1(i,j,k)
        pt = p1(i2,j2,k2)*(1.-omegax(i2))  + (dt/d)*(xz(i+1,j,k) - xz(i,j,k))/dh
        p1(i2,j2,k2) = pt/(1.+omegax(i2))

        pt = p2(i2,j2,k2)*(1.-omegay(j2))  + (dt/d)*(yz(i,j,k) - yz(i,j-1,k))/dh
        p2(i2,j2,k2) = pt/(1.+omegay(j2))

        pt = p3(i2,j2,k2)*(1.-omegaz(k2))  + (dt/d)*(zz(i,j,k+1) - zz(i,j,k))/dh
        p3(i2,j2,k2) = pt/(1.+omegaz(k2))

        w1(i,j,k) = p1(i2,j2,k2) + p2(i2,j2,k2) + p3(i2,j2,k2)
      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      return
      end
     subroutine sxxa(nxb,nxe,nyb,nye,nzb,nze,dh,dt, omegax, omegay, omegaz, px1, px2, px3, py1,py2,py3,pz1,pz2,pz3)
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
      USE pml_com
      integer :: nxb,nxe,nyb,nye,nzb,nze,i2,j2,k2
      real    :: dh,dt,dth,xl,xm,a,b, pt, omt
      real    :: omegax(nxe-nxb+1),  omegay(nye-nyb+1), omegaz(nze-nzb+1), diff1, diff2, diff3, cl1, cl2, cl3
      real    :: px1(nxe-nxb+1,nye-nyb+1,nze-nzb+1),px2(nxe-nxb+1,nye-nyb+1,nze-nzb+1),px3(nxe-nxb+1,nye-nyb+1,nze-nzb+1)
      real    :: py1(nxe-nxb+1,nye-nyb+1,nze-nzb+1),py2(nxe-nxb+1,nye-nyb+1,nze-nzb+1),py3(nxe-nxb+1,nye-nyb+1,nze-nzb+1)
      real    :: pz1(nxe-nxb+1,nye-nyb+1,nze-nzb+1),pz2(nxe-nxb+1,nye-nyb+1,nze-nzb+1),pz3(nxe-nxb+1,nye-nyb+1,nze-nzb+1)

!----------------------------------------------------------
!     Find displacement fields at time t+1/2
!----------------------------------------------------------
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
        xl = lam1(i,j,k)
        xm = mu1(i,j,k)
        a  = xl + 2.*xm
        b  = xl
		
		diff1=u1(i+1,j,k) - u1(i,j,k)
		diff2=v1(i,j,k) - v1(i,j-1,k)
		diff3=w1(i,j,k) - w1(i,j,k-1)
		
!       Find xx stress
        
        pt = px1(i2,j2,k2)*(1.-omegax(i2)) + dt*a*diff1/dh
        px1(i2,j2,k2) = pt/(1.+omegax(i2))

        pt = px2(i2,j2,k2)*(1.-omegay(j2)) + dt*b*diff2/dh
        px2(i2,j2,k2) = pt/(1.+omegay(j2))

        pt = px3(i2,j2,k2)*(1.-omegaz(k2)) + dt*b*diff3/dh
        px3(i2,j2,k2) = pt/(1.+omegaz(k2))

        xx(i,j,k)= px1(i2,j2,k2) + px2(i2,j2,k2) + px3(i2,j2,k2)

!       Find yy stress

        pt = py1(i2,j2,k2)*(1.-omegax(i2)) + dt*b*diff1/dh
        py1(i2,j2,k2) = pt/(1.+omegax(i2))

        pt = py2(i2,j2,k2)*(1.-omegay(j2)) + dt*a*diff2/dh
        py2(i2,j2,k2) = pt/(1.+omegay(j2))

        pt = py3(i2,j2,k2)*(1.-omegaz(k2)) + dt*b*diff3/dh
        py3(i2,j2,k2) = pt/(1.+omegaz(k2))

        yy(i,j,k)= py1(i2,j2,k2) + py2(i2,j2,k2) + py3(i2,j2,k2)
!       Find zz stress
        pt = pz1(i2,j2,k2)*(1.-omegax(i2)) + dt*b*diff1/dh
        pz1(i2,j2,k2) = pt/(1.+omegax(i2))

        pt = pz2(i2,j2,k2)*(1.-omegay(j2)) + dt*b*diff2/dh
        pz2(i2,j2,k2) = pt/(1.+omegay(j2))

        pt = pz3(i2,j2,k2)*(1.-omegaz(k2)) + dt*a*diff3/dh
        pz3(i2,j2,k2) = pt/(1.+omegaz(k2))

        zz(i,j,k)= pz1(i2,j2,k2) + pz2(i2,j2,k2) + pz3(i2,j2,k2)

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
      subroutine sxya(nxb,nxe,nyb,nye,nzb,nze,dh,dt, omegax, omegay, omegaz,p1,p2)
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
      USE pml_com
      integer :: nxb,nxe,nyb,nye,nzb,nze,i2,j2,k2
      real    :: dh,dt,dth,xm1,xm2,xmu, pt
      real    :: omegax(nxe-nxb+1),  omegay(nye-nyb+1), omegaz(nze-nzb+1)
      real    :: p1(nxe-nxb+1,nye-nyb+1,nze-nzb+1),p2(nxe-nxb+1,nye-nyb+1,nze-nzb+1)

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!       Find xy stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j+1,k)
        xmu = (xm1+xm2)/2.

        pt = p1(i2,j2,k2)*(1.-omegay(j2)) + dt*xmu*(u1(i,j+1,k) - u1(i,j,k))/dh
        p1(i2,j2,k2) = pt/(1.+omegay(j2))

        pt = p2(i2,j2,k2)*(1.-omegax(i2)) + dt*xmu*(v1(i,j,k) - v1(i-1,j,k))/dh
        p2(i2,j2,k2) = pt/(1.+omegax(i2))

        xy(i,j,k)= p1(i2,j2,k2) + p2(i2,j2,k2)

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
      subroutine sxza(nxb,nxe,nyb,nye,nzb,nze,dh,dt, omegax, omegay, omegaz,p1,p2)
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
      USE pml_com
      integer :: nxb,nxe,nyb,nye,nzb,nze,i2,j2,k2
      real    :: dh,dt,dth,xm1,xm2,xmu, pt
      real    :: omegax(nxe-nxb+1),  omegay(nye-nyb+1), omegaz(nze-nzb+1)
      real    :: p1(nxe-nxb+1,nye-nyb+1,nze-nzb+1),p2(nxe-nxb+1,nye-nyb+1,nze-nzb+1)


      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k=nzb,nze
      do j=nyb,nye
      !$ACC LOOP VECTOR
      do i=nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!       Find xz stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j,k+1)
        xmu = (xm1+xm2)/2.

        pt = p1(i2,j2,k2)*(1.-omegaz(k2)) + dt*xmu*(u1(i,j,k+1) - u1(i,j,k))/dh
        p1(i2,j2,k2) = pt/(1.+omegaz(k2))

        pt = p2(i2,j2,k2)*(1.-omegax(i2)) + dt*xmu*(w1(i,j,k) - w1(i-1,j,k))/dh
        p2(i2,j2,k2) = pt/(1.+omegax(i2))

        xz(i,j,k)= p1(i2,j2,k2) + p2(i2,j2,k2)

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
      subroutine syza(nxb,nxe,nyb,nye,nzb,nze,dh,dt, omegax, omegay, omegaz,p1,p2)
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
      USE pml_com
      integer :: nxb,nxe,nyb,nye,nzb,nze,i2,j2,k2
      real    :: dh,dt,xm1,xm2,xmu,pt
      real    ::  omegax(nxe-nxb+1),  omegay(nye-nyb+1), omegaz(nze-nzb+1)
      real    :: p1(nxe-nxb+1,nye-nyb+1,nze-nzb+1),p2(nxe-nxb+1,nye-nyb+1,nze-nzb+1)


      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k=nzb,nze
      do j=nyb,nye
      !$ACC LOOP VECTOR
      do i=nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!        Find yz stress

         xm1 = mu1(i,j,k)
         xm2 = mu1(i+1,j+1,k+1)
         xmu = (xm1+xm2)/2.

        pt = p1(i2,j2,k2)*(1.-omegaz(k2)) + dt*xmu*(v1(i,j,k+1) - v1(i,j,k))/dh
        p1(i2,j2,k2) = pt/(1.+omegaz(k2))

        pt = p2(i2,j2,k2)*(1.-omegay(j2)) + dt*xmu*(w1(i,j+1,k) - w1(i,j,k))/dh
        p2(i2,j2,k2) = pt/(1.+omegay(j2))

        yz(i,j,k)= p1(i2,j2,k2) + p2(i2,j2,k2)

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

    !Perfectly matched layers

     subroutine pml_uvw (nxt,nyt,nzt,dt,dh)

      USE medium_com
      USE displt_com
      USE strfld_com
      USE traction_com
      USE pml_com

      real    ::  dh, dt, d, dth,pt
      integer :: nxt, nyt, nzt, nxb, nxe, nyb, nye, nzb, nze
	  integer :: i,j,k,i2,j2,k2

      nxb=2
      nxe=nabc
      nyb=nabc+1
      nye=nyt-1
      nzb=nabc+1
      nze=nzt-nfs

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
		do j = nyb,nye+1
      !$ACC LOOP VECTOR
          do i = nxb,nxe
            k2=1-nzb+k
			j2=1-nyb+j
			i2=1-nxb+i
			d = d1(i,j,k)
					
			pt = (u11(i2,j2,k2)*(1.-omegax1(i2))  + (dt/d)*(xx(i,j,k) - xx(i-1,j,k))/dh)
			u11(i2,j2,k2) = pt/(1.+omegax1(i2))
		
			pt = (u12(i2,j2,k2)*(1.-omegay1(j2))  + (dt/d)*(xy(i,j,k) - xy(i,j-1,k))/dh)
			u12(i2,j2,k2) = pt/(1.+omegay1(j2))
		
			pt = (u13(i2,j2,k2)*(1.-omegaz1(k2))  + (dt/d)*(xz(i,j,k) - xz(i,j,k-1))/dh)
			u13(i2,j2,k2) = pt/(1.+omegaz1(k2))
		
			u1(i,j,k) = u11(i2,j2,k2) + u12(i2,j2,k2) + u13(i2,j2,k2)
			
			
			pt = w11(i2,j2,k2)*(1.-omegaxS1(i2))  + (dt/d)*(xz(i+1,j,k) - xz(i,j,k))/dh
			w11(i2,j2,k2) = pt/(1.+omegaxS1(i2))

			pt = w12(i2,j2,k2)*(1.-omegay1(j2))  + (dt/d)*(yz(i,j,k) - yz(i,j-1,k))/dh
			w12(i2,j2,k2) = pt/(1.+omegay1(j2))

			pt = w13(i2,j2,k2)*(1.-omegaz1(k2))  + (dt/d)*(zz(i,j,k+1) - zz(i,j,k))/dh
			w13(i2,j2,k2) = pt/(1.+omegaz1(k2))

			w1(i,j,k) = w11(i2,j2,k2) + w12(i2,j2,k2) + w13(i2,j2,k2)
			
          enddo
        enddo
      enddo
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
		do j = nyb,nye
      !$ACC LOOP VECTOR
			do i = nxb,nxe
            k2=1-nzb+k
			j2=1-nyb+j
			d         = d1(i,j,k)
			i2=1-nxb+i
		
			pt = v11(i2,j2,k2)*(1.-omegaxS1(i2))  + (dt/d)*(xy(i+1,j,k) - xy(i,j,k))/dh
			v11(i2,j2,k2) = pt/(1.+omegaxS1(i2))

			pt = v12(i2,j2,k2)*(1.-omegay1(j2))  + (dt/d)*(yy(i,j+1,k) - yy(i,j,k))/dh
			v12(i2,j2,k2) = pt/(1.+omegay1(j2))

			pt = v13(i2,j2,k2)*(1.-omegaz1(k2))  + (dt/d)*(yz(i,j,k) - yz(i,j,k-1))/dh
			v13(i2,j2,k2) = pt/(1.+omegaz1(k2))

			v1(i,j,k) = v11(i2,j2,k2) + v12(i2,j2,k2) + v13(i2,j2,k2)
			enddo
		enddo
      enddo
      !$ACC END PARALLEL

    !  call uxxa(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt, omegax1, omegay1, omegaz1, u11, u12, u13)
   !   call vyya(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegaxS1b,omegay1b,omegaz1b,v11,v12,v13)
   !   call wzza(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt, omegaxS1, omegay1, omegaz1, w11, w12, w13)


!     x near nxt

      nxb=nxt-nabc+1
      nxe=nxt-1
      nyb=nabc+1
      nye=nyt-1
      nzb=nabc+1
      nze=nzt-nfs

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
		do j = nyb,nye+1
      !$ACC LOOP VECTOR
            do i = nxb,nxe
            k2=1-nzb+k
			j2=1-nyb+j

			d = d1(i,j,k)
			i2=1-nxb+i
		
			pt = (u21(i2,j2,k2)*(1.-omegax2(i2))  + (dt/d)*(xx(i,j,k) - xx(i-1,j,k))/dh)
			u21(i2,j2,k2) = pt/(1.+omegax2(i2))
		
			pt = (u22(i2,j2,k2)*(1.-omegay2(j2))  + (dt/d)*(xy(i,j,k) - xy(i,j-1,k))/dh)
			u22(i2,j2,k2) = pt/(1.+omegay2(j2))
		
			pt = (u23(i2,j2,k2)*(1.-omegaz2(k2))  + (dt/d)*(xz(i,j,k) - xz(i,j,k-1))/dh)
			u23(i2,j2,k2) = pt/(1.+omegaz2(k2))
		
			u1(i,j,k) = u21(i2,j2,k2) + u22(i2,j2,k2) + u23(i2,j2,k2)
			
			
			pt = w21(i2,j2,k2)*(1.-omegaxS2(i2))  + (dt/d)*(xz(i+1,j,k) - xz(i,j,k))/dh
			w21(i2,j2,k2) = pt/(1.+omegaxS2(i2))

			pt = w22(i2,j2,k2)*(1.-omegay2(j2))  + (dt/d)*(yz(i,j,k) - yz(i,j-1,k))/dh
			w22(i2,j2,k2) = pt/(1.+omegay2(j2))

			pt = w23(i2,j2,k2)*(1.-omegaz2(k2))  + (dt/d)*(zz(i,j,k+1) - zz(i,j,k))/dh
			w23(i2,j2,k2) = pt/(1.+omegaz2(k2))

			w1(i,j,k) = w21(i2,j2,k2) + w22(i2,j2,k2) + w23(i2,j2,k2)
			
			enddo
        enddo
      enddo
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
		do j = nyb,nye
      !$ACC LOOP VECTOR
			do i = nxb,nxe
            k2=1-nzb+k
			j2=1-nyb+j
			d         = d1(i,j,k)
			i2=1-nxb+i
		
			pt = v21(i2,j2,k2)*(1.-omegaxS2(i2))  + (dt/d)*(xy(i+1,j,k) - xy(i,j,k))/dh
			v21(i2,j2,k2) = pt/(1.+omegaxS2(i2))

			pt = v22(i2,j2,k2)*(1.-omegay2(j2))  + (dt/d)*(yy(i,j+1,k) - yy(i,j,k))/dh
			v22(i2,j2,k2) = pt/(1.+omegay2(j2))

			pt = v23(i2,j2,k2)*(1.-omegaz2(k2))  + (dt/d)*(yz(i,j,k) - yz(i,j,k-1))/dh
			v23(i2,j2,k2) = pt/(1.+omegaz2(k2))

			v1(i,j,k) = v21(i2,j2,k2) + v22(i2,j2,k2) + v23(i2,j2,k2)
			enddo
		enddo
      enddo
      !$ACC END PARALLEL
   !   call uxxa(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt,omegax2, omegay2, omegaz2, u21, u22, u23)
   !   call vyya(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegaxS2b,omegay2b,omegaz2b,v21, v22, v23)
   !   call wzza(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt,omegaxS2, omegay2, omegaz2, w21, w22, w23)

!     y near 1

      nxb=2
      nxe=nxt-1
      nyb=2
      nye=nabc
      nzb=nabc+1
      nze=nzt-nfs
	  
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
		do j = nyb,nye
      !$ACC LOOP VECTOR
          do i = nxb,nxe
            k2=1-nzb+k
			j2=1-nyb+j

			d = d1(i,j,k)
			i2=1-nxb+i
		
			pt = (u31(i2,j2,k2)*(1.-omegax3(i2))  + (dt/d)*(xx(i,j,k) - xx(i-1,j,k))/dh)
			u31(i2,j2,k2) = pt/(1.+omegax3(i2))
		
			pt = (u32(i2,j2,k2)*(1.-omegay3(j2))  + (dt/d)*(xy(i,j,k) - xy(i,j-1,k))/dh)
			u32(i2,j2,k2) = pt/(1.+omegay3(j2))
		
			pt = (u33(i2,j2,k2)*(1.-omegaz3(k2))  + (dt/d)*(xz(i,j,k) - xz(i,j,k-1))/dh)
			u33(i2,j2,k2) = pt/(1.+omegaz3(k2))
		
			u1(i,j,k) = u31(i2,j2,k2) + u32(i2,j2,k2) + u33(i2,j2,k2)
			
					
			pt = v31(i2,j2,k2)*(1.-omegaxS3(i2))  + (dt/d)*(xy(i+1,j,k) - xy(i,j,k))/dh
			v31(i2,j2,k2) = pt/(1.+omegaxS3(i2))

			pt = v32(i2,j2,k2)*(1.-omegay3(j2))  + (dt/d)*(yy(i,j+1,k) - yy(i,j,k))/dh
			v32(i2,j2,k2) = pt/(1.+omegay3(j2))

			pt = v33(i2,j2,k2)*(1.-omegaz3(k2))  + (dt/d)*(yz(i,j,k) - yz(i,j,k-1))/dh
			v33(i2,j2,k2) = pt/(1.+omegaz3(k2))

			v1(i,j,k) = v31(i2,j2,k2) + v32(i2,j2,k2) + v33(i2,j2,k2)

			
			pt = w31(i2,j2,k2)*(1.-omegaxS3(i2))  + (dt/d)*(xz(i+1,j,k) - xz(i,j,k))/dh
			w31(i2,j2,k2) = pt/(1.+omegaxS3(i2))

			pt = w32(i2,j2,k2)*(1.-omegay3(j2))  + (dt/d)*(yz(i,j,k) - yz(i,j-1,k))/dh
			w32(i2,j2,k2) = pt/(1.+omegay3(j2))

			pt = w33(i2,j2,k2)*(1.-omegaz3(k2))  + (dt/d)*(zz(i,j,k+1) - zz(i,j,k))/dh
			w33(i2,j2,k2) = pt/(1.+omegaz3(k2))

			w1(i,j,k) = w31(i2,j2,k2) + w32(i2,j2,k2) + w33(i2,j2,k2)
			
			enddo
		enddo	  
      enddo
      !$ACC END PARALLEL

     ! call uxxa(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegax3, omegay3, omegaz3, u31, u32, u33)
     ! call vyya(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegaxS3, omegayS3, omegaz3, v31, v32, v33)
      !call wzza(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegaxS3, omegay3, omegaz3, w31, w32, w33)

!     z near 1

      nxb=2
      nxe=nxt-1
      nyb=2
      nye=nyt-1
      nzb=2
      nze=nabc

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
		do j = nyb,nye+1
      !$ACC LOOP VECTOR
          do i = nxb,nxe
            k2=1-nzb+k
			j2=1-nyb+j

			d = d1(i,j,k)
			i2=1-nxb+i
		
			pt = (u41(i2,j2,k2)*(1.-omegax4(i2))  + (dt/d)*(xx(i,j,k) - xx(i-1,j,k))/dh)
			u41(i2,j2,k2) = pt/(1.+omegax4(i2))
		
			pt = (u42(i2,j2,k2)*(1.-omegay4(j2))  + (dt/d)*(xy(i,j,k) - xy(i,j-1,k))/dh)
			u42(i2,j2,k2) = pt/(1.+omegay4(j2))
		
			pt = (u43(i2,j2,k2)*(1.-omegaz4(k2))  + (dt/d)*(xz(i,j,k) - xz(i,j,k-1))/dh)
			u43(i2,j2,k2) = pt/(1.+omegaz4(k2))
		
			u1(i,j,k) = u41(i2,j2,k2) + u42(i2,j2,k2) + u43(i2,j2,k2)
			
			
			pt = w41(i2,j2,k2)*(1.-omegaxS4(i2))  + (dt/d)*(xz(i+1,j,k) - xz(i,j,k))/dh
			w41(i2,j2,k2) = pt/(1.+omegaxS4(i2))

			pt = w42(i2,j2,k2)*(1.-omegay4(j2))  + (dt/d)*(yz(i,j,k) - yz(i,j-1,k))/dh
			w42(i2,j2,k2) = pt/(1.+omegay4(j2))

			pt = w43(i2,j2,k2)*(1.-omegaz4(k2))  + (dt/d)*(zz(i,j,k+1) - zz(i,j,k))/dh
			w43(i2,j2,k2) = pt/(1.+omegaz4(k2))

			w1(i,j,k) = w41(i2,j2,k2) + w42(i2,j2,k2) + w43(i2,j2,k2)
			
			enddo
        enddo
      enddo
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
		do j = nyb,nye
      !$ACC LOOP VECTOR
			do i = nxb,nxe
            k2=1-nzb+k
			j2=1-nyb+j
			d         = d1(i,j,k)
			i2=1-nxb+i
		
			pt = v41(i2,j2,k2)*(1.-omegaxS4(i2))  + (dt/d)*(xy(i+1,j,k) - xy(i,j,k))/dh
			v41(i2,j2,k2) = pt/(1.+omegaxS4(i2))

			pt = v42(i2,j2,k2)*(1.-omegay4(j2))  + (dt/d)*(yy(i,j+1,k) - yy(i,j,k))/dh
			v42(i2,j2,k2) = pt/(1.+omegay4(j2))

			pt = v43(i2,j2,k2)*(1.-omegaz4(k2))  + (dt/d)*(yz(i,j,k) - yz(i,j,k-1))/dh
			v43(i2,j2,k2) = pt/(1.+omegaz4(k2))

			v1(i,j,k) = v41(i2,j2,k2) + v42(i2,j2,k2) + v43(i2,j2,k2)
			enddo
		enddo
      enddo
      !$ACC END PARALLEL
      !call uxxa(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt,omegax4, omegay4, omegaz4, u41, u42, u43)
      !call vyya(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegaxS4b, omegayS4b, omegaz4b,v41,v42,v43)
      !call wzza(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt, omegaxS4, omegay4, omegazS4, w41, w42, w43)

     end

     subroutine pml_xyz (nxt,nyt,nzt,dt,dh)

      USE medium_com
      USE displt_com
      USE strfld_com
      USE traction_com
      USE pml_com

      real    :: dh, dt,a,b,diff1,diff2,diff3,xl,xm,pt
      integer :: nxt, nyt, nzt, nxb, nxe, nyb, nye, nzb, nze


!     near x=1
      nxb=2
      nxe=nabc
      nyb=nabc+1
      nye=nyt-1
      nzb=nabc+1
      nze=nzt-nfs

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k = nzb,nze
      do j = nyb,nye+1
      !$ACC LOOP VECTOR
      do i = nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
        xl = lam1(i,j,k)
        xm = mu1(i,j,k)
        a  = xl + 2.*xm
        b  = xl
		
		diff1=u1(i+1,j,k) - u1(i,j,k)
		diff2=v1(i,j,k) - v1(i,j-1,k)
		diff3=w1(i,j,k) - w1(i,j,k-1)
		
!       Find xx stress
        
        pt = xx11(i2,j2,k2)*(1.-omegaxS1(i2)) + dt*a*diff1/dh
        xx11(i2,j2,k2) = pt/(1.+omegaxS1(i2))

        pt = xx12(i2,j2,k2)*(1.-omegay1(j2)) + dt*b*diff2/dh
        xx12(i2,j2,k2) = pt/(1.+omegay1(j2))

        pt = xx13(i2,j2,k2)*(1.-omegaz1(k2)) + dt*b*diff3/dh
        xx13(i2,j2,k2) = pt/(1.+omegaz1(k2))

        xx(i,j,k)= xx11(i2,j2,k2) + xx12(i2,j2,k2) + xx13(i2,j2,k2)

!       Find yy stress

        pt = yy11(i2,j2,k2)*(1.-omegaxS1(i2)) + dt*b*diff1/dh
        yy11(i2,j2,k2) = pt/(1.+omegaxS1(i2))

        pt = yy12(i2,j2,k2)*(1.-omegay1(j2)) + dt*a*diff2/dh
        yy12(i2,j2,k2) = pt/(1.+omegay1(j2))

        pt = yy13(i2,j2,k2)*(1.-omegaz1(k2)) + dt*b*diff3/dh
        yy13(i2,j2,k2) = pt/(1.+omegaz1(k2))

        yy(i,j,k)= yy11(i2,j2,k2) + yy12(i2,j2,k2) + yy13(i2,j2,k2)
!       Find zz stress
        pt = zz11(i2,j2,k2)*(1.-omegaxS1(i2)) + dt*b*diff1/dh
        zz11(i2,j2,k2) = pt/(1.+omegaxS1(i2))

        pt = zz12(i2,j2,k2)*(1.-omegay1(j2)) + dt*b*diff2/dh
        zz12(i2,j2,k2) = pt/(1.+omegay1(j2))

        pt = zz13(i2,j2,k2)*(1.-omegaz1(k2)) + dt*a*diff3/dh
        zz13(i2,j2,k2) = pt/(1.+omegaz1(k2))

        zz(i,j,k)= zz11(i2,j2,k2) + zz12(i2,j2,k2) + zz13(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL
	  
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!       Find xy stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j+1,k)
        xmu = (xm1+xm2)/2.

        pt = xy11(i2,j2,k2)*(1.-omegay1(j2)) + dt*xmu*(u1(i,j+1,k) - u1(i,j,k))/dh
        xy11(i2,j2,k2) = pt/(1.+omegay1(j2))

        pt = xy12(i2,j2,k2)*(1.-omegax1(i2)) + dt*xmu*(v1(i,j,k) - v1(i-1,j,k))/dh
        xy12(i2,j2,k2) = pt/(1.+omegax1(i2))

        xy(i,j,k)= xy11(i2,j2,k2) + xy12(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL
	  
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k=nzb,nze
      do j=nyb,nye+1
      !$ACC LOOP VECTOR
      do i=nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!       Find xz stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j,k+1)
        xmu = (xm1+xm2)/2.

        pt = xz11(i2,j2,k2)*(1.-omegaz1(k2)) + dt*xmu*(u1(i,j,k+1) - u1(i,j,k))/dh
        xz11(i2,j2,k2) = pt/(1.+omegaz1(k2))

        pt = xz12(i2,j2,k2)*(1.-omegax1(i2)) + dt*xmu*(w1(i,j,k) - w1(i-1,j,k))/dh
        xz12(i2,j2,k2) = pt/(1.+omegax1(i2))

        xz(i,j,k)= xz11(i2,j2,k2) + xz12(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL
	  
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k=nzb,nze
      do j=nyb,nye
      !$ACC LOOP VECTOR
      do i=nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!        Find yz stress

         xm1 = mu1(i,j,k)
         xm2 = mu1(i+1,j+1,k+1)
         xmu = (xm1+xm2)/2.

        pt = yz11(i2,j2,k2)*(1.-omegaz1(k2)) + dt*xmu*(v1(i,j,k+1) - v1(i,j,k))/dh
        yz11(i2,j2,k2) = pt/(1.+omegaz1(k2))

        pt = yz12(i2,j2,k2)*(1.-omegay1(j2)) + dt*xmu*(w1(i,j+1,k) - w1(i,j,k))/dh
        yz12(i2,j2,k2) = pt/(1.+omegay1(j2))

        yz(i,j,k)= yz11(i2,j2,k2) + yz12(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL
	  
     ! call sxxa(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt,omegaxS1, omegay1, omegaz1, xx11,xx12,xx13,yy11,yy12,yy13,zz11,zz12,zz13)
     ! call sxya(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegax1b,omegay1b,omegaz1b, xy11,xy12)
     ! call sxza(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt,omegax1, omegay1, omegaz1, xz11,xz12)
     ! call syza(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegaxS1b,omegay1b,omegaz1b, yz11,yz12)


!      xz11(:,:,nze-nzb+1) = 0.
!      xz12(:,:,nze-nzb+1) = 0.
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
	  do i=1,nxe-nxb+1
      !$ACC LOOP VECTOR
        do j=1,nye-nyb+1+1
          xz11(i,j,nze-nzb+1) = 0.
          xz12(i,j,nze-nzb+1) = 0.
        enddo
      enddo
      !$ACC END PARALLEL

      
!      yz11(:,:,nze-nzb+1) = 0.
!      yz12(:,:,nze-nzb+1) = 0.
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
	  do i=1,nxe-nxb+1
      !$ACC LOOP VECTOR
        do j=1,nye-nyb+1
          yz11(i,j,nze-nzb+1) = 0.
          yz12(i,j,nze-nzb+1) = 0.
        enddo
      enddo
      !$ACC END PARALLEL
      
!     near x=nxt

      nxb=nxt-nabc+1
      nxe=nxt-1
      nyb=nabc+1
      nye=nyt-1
      nzb=nabc+1
      nze=nzt-nfs
	  
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k = nzb,nze
      do j = nyb,nye+1
      !$ACC LOOP VECTOR
      do i = nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
        xl = lam1(i,j,k)
        xm = mu1(i,j,k)
        a  = xl + 2.*xm
        b  = xl
		
		diff1=u1(i+1,j,k) - u1(i,j,k)
		diff2=v1(i,j,k) - v1(i,j-1,k)
		diff3=w1(i,j,k) - w1(i,j,k-1)
		
!       Find xx stress
        
        pt = xx21(i2,j2,k2)*(1.-omegaxS2(i2)) + dt*a*diff1/dh
        xx21(i2,j2,k2) = pt/(1.+omegaxS2(i2))

        pt = xx22(i2,j2,k2)*(1.-omegay2(j2)) + dt*b*diff2/dh
        xx22(i2,j2,k2) = pt/(1.+omegay2(j2))

        pt = xx23(i2,j2,k2)*(1.-omegaz2(k2)) + dt*b*diff3/dh
        xx23(i2,j2,k2) = pt/(1.+omegaz2(k2))

        xx(i,j,k)= xx21(i2,j2,k2) + xx22(i2,j2,k2) + xx23(i2,j2,k2)

!       Find yy stress

        pt = yy21(i2,j2,k2)*(1.-omegaxS2(i2)) + dt*b*diff1/dh
        yy21(i2,j2,k2) = pt/(1.+omegaxS2(i2))

        pt = yy22(i2,j2,k2)*(1.-omegay2(j2)) + dt*a*diff2/dh
        yy22(i2,j2,k2) = pt/(1.+omegay2(j2))

        pt = yy23(i2,j2,k2)*(1.-omegaz2(k2)) + dt*b*diff3/dh
        yy23(i2,j2,k2) = pt/(1.+omegaz2(k2))

        yy(i,j,k)= yy21(i2,j2,k2) + yy22(i2,j2,k2) + yy23(i2,j2,k2)
!       Find zz stress
        pt = zz21(i2,j2,k2)*(1.-omegaxS2(i2)) + dt*b*diff1/dh
        zz21(i2,j2,k2) = pt/(1.+omegaxS2(i2))

        pt = zz22(i2,j2,k2)*(1.-omegay2(j2)) + dt*b*diff2/dh
        zz22(i2,j2,k2) = pt/(1.+omegay2(j2))

        pt = zz23(i2,j2,k2)*(1.-omegaz2(k2)) + dt*a*diff3/dh
        zz23(i2,j2,k2) = pt/(1.+omegaz2(k2))

        zz(i,j,k)= zz21(i2,j2,k2) + zz22(i2,j2,k2) + zz23(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL
	  
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!       Find xy stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j+1,k)
        xmu = (xm1+xm2)/2.

        pt = xy21(i2,j2,k2)*(1.-omegay2(j2)) + dt*xmu*(u1(i,j+1,k) - u1(i,j,k))/dh
        xy21(i2,j2,k2) = pt/(1.+omegay2(j2))

        pt = xy22(i2,j2,k2)*(1.-omegax2(i2)) + dt*xmu*(v1(i,j,k) - v1(i-1,j,k))/dh
        xy22(i2,j2,k2) = pt/(1.+omegax2(i2))

        xy(i,j,k)= xy21(i2,j2,k2) + xy22(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL
	  
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k=nzb,nze
      do j=nyb,nye+1
      !$ACC LOOP VECTOR
      do i=nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!       Find xz stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j,k+1)
        xmu = (xm1+xm2)/2.

        pt = xz21(i2,j2,k2)*(1.-omegaz2(k2)) + dt*xmu*(u1(i,j,k+1) - u1(i,j,k))/dh
        xz21(i2,j2,k2) = pt/(1.+omegaz2(k2))

        pt = xz22(i2,j2,k2)*(1.-omegax2(i2)) + dt*xmu*(w1(i,j,k) - w1(i-1,j,k))/dh
        xz22(i2,j2,k2) = pt/(1.+omegax2(i2))

        xz(i,j,k)= xz21(i2,j2,k2) + xz22(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k=nzb,nze
      do j=nyb,nye
      !$ACC LOOP VECTOR
      do i=nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!        Find yz stress

         xm1 = mu1(i,j,k)
         xm2 = mu1(i+1,j+1,k+1)
         xmu = (xm1+xm2)/2.

        pt = yz21(i2,j2,k2)*(1.-omegaz2(k2)) + dt*xmu*(v1(i,j,k+1) - v1(i,j,k))/dh
        yz21(i2,j2,k2) = pt/(1.+omegaz2(k2))

        pt = yz22(i2,j2,k2)*(1.-omegay2(j2)) + dt*xmu*(w1(i,j+1,k) - w1(i,j,k))/dh
        yz22(i2,j2,k2) = pt/(1.+omegay2(j2))

        yz(i,j,k)= yz21(i2,j2,k2) + yz22(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      !call sxxa(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt,omegaxS2, omegay2, omegaz2,xx21,xx22,xx23,yy21,yy22,yy23,zz21,zz22,zz23)
      !call sxya(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegax2b,omegay2b,omegaz2b, xy21,xy22)
      !call sxza(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt,omegax2, omegay2, omegaz2,xz21,xz22)
      !call syza(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegaxS2b,omegay2b,omegaz2b, yz21,yz22)

!      xz21(:,:,nze-nzb+1) = 0.
!      xz22(:,:,nze-nzb+1) = 0.
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
	  do i=1,nxe-nxb+1
      !$ACC LOOP VECTOR
        do j=1,nye-nyb+1+1
          xz21(i,j,nze-nzb+1) = 0.
          xz22(i,j,nze-nzb+1) = 0.
        enddo
      enddo
      !$ACC END PARALLEL

      
!      yz21(:,:,nze-nzb+1) = 0.
!      yz22(:,:,nze-nzb+1) = 0.
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
	  do i=1,nxe-nxb+1
      !$ACC LOOP VECTOR
        do j=1,nye-nyb+1
          yz21(i,j,nze-nzb+1) = 0.
          yz22(i,j,nze-nzb+1) = 0.
        enddo
      enddo
      !$ACC END PARALLEL

!     near y=1

      nxb=2
      nxe=nxt-1
      nyb=2
      nye=nabc
      nzb=nabc+1
      nze=nzt-nfs
	  
      !call sxxa(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegaxS3, omegay3, omegaz3,xx31,xx32,xx33,yy31,yy32,yy33,zz31,zz32,zz33)
      !call sxya(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegax3, omegayS3, omegaz3, xy31,xy32)
      !call sxza(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegax3, omegay3, omegaz3,xz31,xz32)
      !call syza(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegaxS3, omegayS3, omegaz3, yz31,yz32)

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
        xl = lam1(i,j,k)
        xm = mu1(i,j,k)
        a  = xl + 2.*xm
        b  = xl
		
		diff1=u1(i+1,j,k) - u1(i,j,k)
		diff2=v1(i,j,k) - v1(i,j-1,k)
		diff3=w1(i,j,k) - w1(i,j,k-1)
		
!       Find xx stress
        
        pt = xx31(i2,j2,k2)*(1.-omegaxS3(i2)) + dt*a*diff1/dh
        xx31(i2,j2,k2) = pt/(1.+omegaxS3(i2))

        pt = xx32(i2,j2,k2)*(1.-omegay3(j2)) + dt*b*diff2/dh
        xx32(i2,j2,k2) = pt/(1.+omegay3(j2))

        pt = xx33(i2,j2,k2)*(1.-omegaz3(k2)) + dt*b*diff3/dh
        xx33(i2,j2,k2) = pt/(1.+omegaz3(k2))

        xx(i,j,k)= xx31(i2,j2,k2) + xx32(i2,j2,k2) + xx33(i2,j2,k2)

!       Find yy stress

        pt = yy31(i2,j2,k2)*(1.-omegaxS3(i2)) + dt*b*diff1/dh
        yy31(i2,j2,k2) = pt/(1.+omegaxS3(i2))

        pt = yy32(i2,j2,k2)*(1.-omegay3(j2)) + dt*a*diff2/dh
        yy32(i2,j2,k2) = pt/(1.+omegay3(j2))

        pt = yy33(i2,j2,k2)*(1.-omegaz3(k2)) + dt*b*diff3/dh
        yy33(i2,j2,k2) = pt/(1.+omegaz3(k2))

        yy(i,j,k)= yy31(i2,j2,k2) + yy32(i2,j2,k2) + yy33(i2,j2,k2)
!       Find zz stress
        pt = zz31(i2,j2,k2)*(1.-omegaxS3(i2)) + dt*b*diff1/dh
        zz31(i2,j2,k2) = pt/(1.+omegaxS3(i2))

        pt = zz32(i2,j2,k2)*(1.-omegay3(j2)) + dt*b*diff2/dh
        zz32(i2,j2,k2) = pt/(1.+omegay3(j2))

        pt = zz33(i2,j2,k2)*(1.-omegaz3(k2)) + dt*a*diff3/dh
        zz33(i2,j2,k2) = pt/(1.+omegaz3(k2))

        zz(i,j,k)= zz31(i2,j2,k2) + zz32(i2,j2,k2) + zz33(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!       Find xy stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j+1,k)
        xmu = (xm1+xm2)/2.

        pt = xy31(i2,j2,k2)*(1.-omegayS3(j2)) + dt*xmu*(u1(i,j+1,k) - u1(i,j,k))/dh
        xy31(i2,j2,k2) = pt/(1.+omegayS3(j2))

        pt = xy32(i2,j2,k2)*(1.-omegax3(i2)) + dt*xmu*(v1(i,j,k) - v1(i-1,j,k))/dh
        xy32(i2,j2,k2) = pt/(1.+omegax3(i2))

        xy(i,j,k)= xy31(i2,j2,k2) + xy32(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL
	  
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k=nzb,nze
      do j=nyb,nye
      !$ACC LOOP VECTOR
      do i=nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!       Find xz stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j,k+1)
        xmu = (xm1+xm2)/2.

        pt = xz31(i2,j2,k2)*(1.-omegaz3(k2)) + dt*xmu*(u1(i,j,k+1) - u1(i,j,k))/dh
        xz31(i2,j2,k2) = pt/(1.+omegaz3(k2))

        pt = xz32(i2,j2,k2)*(1.-omegax3(i2)) + dt*xmu*(w1(i,j,k) - w1(i-1,j,k))/dh
        xz32(i2,j2,k2) = pt/(1.+omegax3(i2))

        xz(i,j,k)= xz31(i2,j2,k2) + xz32(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
      do k=nzb,nze
      do j=nyb,nye
      !$ACC LOOP VECTOR
      do i=nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!        Find yz stress

         xm1 = mu1(i,j,k)
         xm2 = mu1(i+1,j+1,k+1)
         xmu = (xm1+xm2)/2.

        pt = yz31(i2,j2,k2)*(1.-omegaz3(k2)) + dt*xmu*(v1(i,j,k+1) - v1(i,j,k))/dh
        yz31(i2,j2,k2) = pt/(1.+omegaz3(k2))

        pt = yz32(i2,j2,k2)*(1.-omegayS3(j2)) + dt*xmu*(w1(i,j+1,k) - w1(i,j,k))/dh
        yz32(i2,j2,k2) = pt/(1.+omegayS3(j2))

        yz(i,j,k)= yz31(i2,j2,k2) + yz32(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL

!      xz31(:,:,nze-nzb+1) = 0.
!      xz32(:,:,nze-nzb+1) = 0.
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
	  do i=1,nxe-nxb+1
      !$ACC LOOP VECTOR
        do j=1,nye-nyb+1
          xz31(i,j,nze-nzb+1) = 0.
          xz32(i,j,nze-nzb+1) = 0.
        enddo
      enddo
      !$ACC END PARALLEL
      
!      yz31(:,:,nze-nzb+1) = 0.
!      yz32(:,:,nze-nzb+1) = 0.
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
	  do i=1,nxe-nxb+1
      !$ACC LOOP VECTOR
        do j=1,nye-nyb+1
          yz31(i,j,nze-nzb+1) = 0.
          yz32(i,j,nze-nzb+1) = 0.
        enddo
      enddo
      !$ACC END PARALLEL

!     near z=1
      nxb=2
      nxe=nxt-1
      nyb=2
      nye=nyt-1
      nzb=2
      nze=nabc
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k = nzb,nze
      do j = nyb,nye+1
      !$ACC LOOP VECTOR
      do i = nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
        xl = lam1(i,j,k)
        xm = mu1(i,j,k)
        a  = xl + 2.*xm
        b  = xl
		
		diff1=u1(i+1,j,k) - u1(i,j,k)
		diff2=v1(i,j,k) - v1(i,j-1,k)
		diff3=w1(i,j,k) - w1(i,j,k-1)
		
!       Find xx stress
        
        pt = xx41(i2,j2,k2)*(1.-omegaxS4(i2)) + dt*a*diff1/dh
        xx41(i2,j2,k2) = pt/(1.+omegaxS4(i2))

        pt = xx42(i2,j2,k2)*(1.-omegay4(j2)) + dt*b*diff2/dh
        xx42(i2,j2,k2) = pt/(1.+omegay4(j2))

        pt = xx43(i2,j2,k2)*(1.-omegaz4(k2)) + dt*b*diff3/dh
        xx43(i2,j2,k2) = pt/(1.+omegaz4(k2))

        xx(i,j,k)= xx41(i2,j2,k2) + xx42(i2,j2,k2) + xx43(i2,j2,k2)

!       Find yy stress

        pt = yy41(i2,j2,k2)*(1.-omegaxS4(i2)) + dt*b*diff1/dh
        yy41(i2,j2,k2) = pt/(1.+omegaxS4(i2))

        pt = yy42(i2,j2,k2)*(1.-omegay4(j2)) + dt*a*diff2/dh
        yy42(i2,j2,k2) = pt/(1.+omegay4(j2))

        pt = yy43(i2,j2,k2)*(1.-omegaz4(k2)) + dt*b*diff3/dh
        yy43(i2,j2,k2) = pt/(1.+omegaz4(k2))

        yy(i,j,k)= yy41(i2,j2,k2) + yy42(i2,j2,k2) + yy43(i2,j2,k2)
!       Find zz stress
        pt = zz41(i2,j2,k2)*(1.-omegaxS4(i2)) + dt*b*diff1/dh
        zz41(i2,j2,k2) = pt/(1.+omegaxS4(i2))

        pt = zz42(i2,j2,k2)*(1.-omegay4(j2)) + dt*b*diff2/dh
        zz42(i2,j2,k2) = pt/(1.+omegay4(j2))

        pt = zz43(i2,j2,k2)*(1.-omegaz4(k2)) + dt*a*diff3/dh
        zz43(i2,j2,k2) = pt/(1.+omegaz4(k2))

        zz(i,j,k)= zz41(i2,j2,k2) + zz42(i2,j2,k2) + zz43(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL
	  
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k = nzb,nze
      do j = nyb,nye
      !$ACC LOOP VECTOR
      do i = nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!       Find xy stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j+1,k)
        xmu = (xm1+xm2)/2.

        pt = xy41(i2,j2,k2)*(1.-omegayS4(j2)) + dt*xmu*(u1(i,j+1,k) - u1(i,j,k))/dh
        xy41(i2,j2,k2) = pt/(1.+omegayS4(j2))

        pt = xy42(i2,j2,k2)*(1.-omegax4(i2)) + dt*xmu*(v1(i,j,k) - v1(i-1,j,k))/dh
        xy42(i2,j2,k2) = pt/(1.+omegax4(i2))

        xy(i,j,k)= xy41(i2,j2,k2) + xy42(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL
	  
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k=nzb,nze
      do j=nyb,nye+1
      !$ACC LOOP VECTOR
      do i=nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!       Find xz stress

        xm1 = mu1(i,j,k)
        xm2 = mu1(i,j,k+1)
        xmu = (xm1+xm2)/2.

        pt = xz41(i2,j2,k2)*(1.-omegazS4(k2)) + dt*xmu*(u1(i,j,k+1) - u1(i,j,k))/dh
        xz41(i2,j2,k2) = pt/(1.+omegazS4(k2))

        pt = xz42(i2,j2,k2)*(1.-omegax4(i2)) + dt*xmu*(w1(i,j,k) - w1(i-1,j,k))/dh
        xz42(i2,j2,k2) = pt/(1.+omegax4(i2))

        xz(i,j,k)= xz41(i2,j2,k2) + xz42(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL
	  
      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG COLLAPSE (2)
	  do k=nzb,nze
      do j=nyb,nye
      !$ACC LOOP VECTOR
      do i=nxb,nxe
		i2=1-nxb+i
		j2=1-nyb+j
		k2=1-nzb+k
!        Find yz stress

         xm1 = mu1(i,j,k)
         xm2 = mu1(i+1,j+1,k+1)
         xmu = (xm1+xm2)/2.

        pt = yz41(i2,j2,k2)*(1.-omegazS4(k2)) + dt*xmu*(v1(i,j,k+1) - v1(i,j,k))/dh
        yz41(i2,j2,k2) = pt/(1.+omegazS4(k2))

        pt = yz42(i2,j2,k2)*(1.-omegayS4(j2)) + dt*xmu*(w1(i,j+1,k) - w1(i,j,k))/dh
        yz42(i2,j2,k2) = pt/(1.+omegayS4(j2))

        yz(i,j,k)= yz41(i2,j2,k2) + yz42(i2,j2,k2)

      enddo
      enddo
      enddo
      !$ACC END PARALLEL
     ! call sxxa(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt,omegaxS4, omegay4, omegaz4, xx41,xx42,xx43,yy41,yy42,yy43,zz41,zz42,zz43)
     ! call sxya(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegax4b, omegayS4b, omegaz4b, xy41,xy42)
     ! call sxza(nxb,nxe,nyb,nye+1,nzb,nze,dh,dt,omegax4, omegay4, omegazS4, xz41,xz42)
     ! call syza(nxb,nxe,nyb,nye,nzb,nze,dh,dt,omegaxS4b, omegayS4b, omegazS4b, yz41,yz42)


     end

