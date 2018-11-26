    MODULE fd3dparam_com
      integer :: nxt,nyt,nzt,ntfd,nysc
      real    :: dh,dt
    END MODULE

    MODULE medium_com
      real,allocatable,dimension(:,:,:):: lam1,mu1,d1
      real:: mu_mean
    END MODULE

    MODULE friction_com
      USE fd3dparam_com
      real,allocatable,dimension(:,:):: strinix,peak_xz,Dc
      real:: dip,topdepth
      real,parameter:: pi=3.1415926535
      CONTAINS
      FUNCTION normstress(j)
      IMPLICIT NONE
      real:: normstress
      integer:: j
#if defined DIPSLIP
      normstress=max(1.e5,8520.*(dh*real(nzt-j-2)*sin(dip/180.*pi)+topdepth*1000.))
#else
      normstress=max(1.e5,16200.*(dh*real(nzt-j-2)*sin(dip/180.*pi)+topdepth*1000.))
#endif
      END FUNCTION
    END MODULE

    MODULE source_com
      REAL,ALLOCATABLE,DIMENSION(:,:):: ruptime,rise,slip,schange
      real,allocatable,dimension(:,:,:):: sliprate
      real    :: output_param(6)
      integer :: ioutput
    END MODULE


    SUBROUTINE fd3d_init()
      USE medium_com
      USE displt_com
      USE strfld_com
      USE fd3dparam_com
      USE friction_com
      USE source_com
      IMPLICIT NONE

!--------------------
! Read the input file
!--------------------
      write(*,*)'Reading FD3D parameters...'
      open(11, file='inputfd3d.dat', status='old')
      read(11,*) nxt,nyt,nzt
      nysc=nyt-2
      read(11,*) dh
      read(11,*) ntfd
      read(11,*) dt
      read(11,*) dip,topdepth
      close(11)
      topdepth=-2.*dh*sin(dip/180.*pi)/1000.
write(*,*)'Changing topdepth to ',topdepth
!----------------------------
! Allocate FD module arrays
!----------------------------
      allocate(lam1(nxt,nyt,nzt),mu1(nxt,nyt,nzt),d1(nxt,nyt,nzt))
      allocate(u1(nxt,nyt,nzt),v1(nxt,nyt,nzt),w1(nxt,nyt,nzt))
      allocate(xx(nxt,nyt,nzt),yy(nxt,nyt,nzt),zz(nxt,nyt,nzt),xy(nxt,nyt,nzt),yz(nxt,nyt,nzt),xz(nxt,nyt,nzt))
      allocate(strinix(nxt,nzt),peak_xz(nxt,nzt),Dc(nxt,nzt))
      allocate(ruptime(nxt,nzt),slip(nxt,nzt),rise(nxt,nzt),schange(nxt,nzt),sliprate(nxt,nzt,ntfd))

      strinix=0.;peak_xz=0.;Dc=0.

!------------------------------------------------------------
! Read the velocity model
! Be careful: the velocity model for the FD is upside down
!------------------------------------------------------------
      CALL readcrustalmodel(dip,topdepth)

    END


   SUBROUTINE readcrustalmodel(dip,topdepth)
    USE medium_com
    USE fd3dparam_com
    IMPLICIT NONE
    real*8,parameter:: PI=3.1415926535
    real dip,topdepth
    real    :: vpe(2),vse(2),den(2),CFL,dum,dd,vpp,vss
    real,allocatable,dimension(:):: vp,vs,depth,rho
    INTEGER ndepth,j,k

    vpe(2)  = 0.
    vpe(1)  = 1.0E+10
    vse(2)  = 0.
    vse(1)  = 1.0E+10
    den(2)  = 0.
    den(1)  = 1.0E+10
    mu_mean = 0.
    open(10, file='crustal.dat', status='old')
    read(10,*)
    read(10,*)
    read(10,*)ndepth
    allocate(depth(ndepth),vp(ndepth),vs(ndepth),rho(ndepth))
    read(10,*)
    read(10,*)
    do k=1,ndepth
      read(10,*)depth(k),vp(k),vs(k),rho(k)
    enddo
    depth=depth*1.e3;vp=vp*1.e3;vs=vs*1.e3;rho=rho*1.e3
    close(10)
!    open(10,FILE='model.new.dat')
    do k=nzt,1,-1
      dum=topdepth*1.e3+(dh*real(nzt-k+1)-dh/2.)*sin(dip/180.d0*PI)
      if(dum>depth(ndepth))then
        vpp=vp(ndepth)
        vss=vs(ndepth)
        dd=rho(ndepth)
      else
        do j=2,ndepth
          if(dum<depth(j))exit
        enddo
        vpp=vp(j-1)
        vss=vs(j-1)
        dd=rho(j-1)
      endif
!      write(10,*)vpp,vss,dd
      if (vpp.gt.vpe(2)) vpe(2) = vpp
      if (vpp.lt.vpe(1)) vpe(1) = vpp
      if (vss.gt.vse(2)) vse(2) = vss
      if (vss.lt.vse(1)) vse(1) = vss
      if (dd.gt.den(2)) den(2) = dd
      if (dd.lt.den(1)) den(1) = dd
      mu_mean = mu_mean + vss*vss*dd
      lam1(1:nxt,1:nyt,k) = dd*(vpp**2-2.*vss**2)
      mu1(1:nxt,1:nyt,k)  = dd*vss**2
      d1(1:nxt,1:nyt,k)   = dd
    enddo
!    close(10)
    mu_mean = (mu_mean/nzt)
!    write(*,*)mu_mean
    deallocate(depth,vp,vs,rho)

!-------------------------------------------------------
!     Make sure the simulation is numerically stable
!     CFL < 0.25
!-------------------------------------------------------
    CFL = vpe(2)*dt/dh
    if (CFL.gt.0.25) then
      print *,'Your simulation is numerically unstable'
      stop
    endif
    END


