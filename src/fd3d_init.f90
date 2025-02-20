!-------------------------------------------------------
! Initialization and communication modules for the finite difference code
!-------------------------------------------------------
! Authors: Jan Premus and Frantisek Gallovic (8/2019)
! Charles University in Prague, Faculty of Mathematics and Physics

! This code is published under the GNU General Public License. To any
! licensee is given permission to modify the work, as well as to copy
! and redistribute the work or any derivative version. Still we would
! like to kindly ask you to acknowledge the authors and don't remove
! their names from the code. This code is distributed in the hope
! that it will be useful, but WITHOUT ANY WARRANTY.
! ------------------------------------------------------

    MODULE fd3dparam_com
      integer :: nxt,nyt,nzt,ntfd,nysc
      real    :: dh,dt
    END MODULE

    MODULE medium_com
      real,allocatable,dimension(:,:,:):: lam1,mu1,d1
      real,allocatable,dimension(:,:):: vpls
      real:: mu_mean,constnormstress
    END MODULE
    
    MODULE pml_com
      integer :: nabc,nfs  ! number of fringe points
      real :: omegaM_pml
      real,allocatable, dimension (:) :: omega_pml,omegaR_pml,omega_pmlM,omegaR_pmlM
      real,allocatable, dimension (:,:,:):: u11,u12,u13,u21,u22,u23,u31,u32,u33,u41,u42,u43
      real,allocatable, dimension (:,:,:):: v11,v12,v13,v21,v22,v23,v31,v32,v33,v41,v42,v43
      real,allocatable, dimension (:,:,:):: w11,w12,w13,w21,w22,w23,w31,w32,w33,w41,w42,w43
      real,allocatable, dimension (:,:,:):: xx11,xx12,xx13,xx21,xx22,xx23,xx31,xx32,xx33,xx41,xx42,xx43
      real,allocatable, dimension (:,:,:):: yy11,yy12,yy13,yy21,yy22,yy23,yy31,yy32,yy33,yy41,yy42,yy43
      real,allocatable, dimension (:,:,:):: zz11,zz12,zz13,zz21,zz22,zz23,zz31,zz32,zz33,zz41,zz42,zz43
      real,allocatable, dimension (:,:,:):: xy11,xy12,xy21,xy22,xy31,xy32,xy41,xy42
      real,allocatable, dimension (:,:,:):: xz11,xz12,xz21,xz22,xz31,xz32,xz41,xz42
      real,allocatable, dimension (:,:,:):: yz11,yz12,yz21,yz22,yz31,yz32,yz41,yz42
      real,allocatable, dimension (:):: omegax1,omegax2,omegax3,omegax4
      real,allocatable, dimension (:):: omegay1,omegay2,omegay3,omegay4
      real,allocatable, dimension (:):: omegaz1,omegaz2,omegaz3,omegaz4
      real,allocatable, dimension (:):: omegaxS1,omegaxS2,omegaxS3,omegaxS4
      real,allocatable, dimension (:):: omegayS3,omegayS4
      real,allocatable, dimension (:):: omegazS4 
#if defined FSPACE
      real,allocatable, dimension (:,:,:):: u51,u52,u53,v51,v52,v53,w51,w52,w53
	  real,allocatable, dimension (:,:,:):: xx51,xx52,xx53,yy51,yy52,yy53,zz51,zz52,zz53
	  real,allocatable, dimension (:,:,:):: xy51,xy52,xz51,xz52,yz51,yz52
	  real,allocatable, dimension (:):: omegax5,omegay5,omegaz5,omegaxS5,omegayS5,omegazS5 
#endif	  
    END MODULE
    
    MODULE friction_com
      USE fd3dparam_com
      USE pml_com
      
      real,allocatable,dimension(:,:):: uZ,wX
      real,allocatable,dimension(:,:):: striniZ,striniX, T0X, T0Z
      real,allocatable,dimension(:,:):: Dc
      real,allocatable,dimension(:,:):: tabsX,tabsZ
#if defined FVW
      real v0, hx0, hz0, perturb, RR2, TT2, strinixI
      real,allocatable,dimension(:,:):: a, ba, psi, vw, f0, fw, Sn, uini,wini
      real,allocatable,dimension(:,:):: aX, baX, psiX,vwX, f0X, fwX, SnX
      real,allocatable,dimension(:,:):: aZ, baZ, psiZ, vwZ, f0Z, fwZ, SnZ
#endif
      real,allocatable,dimension(:,:):: peak_xz,dyn_xz,coh
      real,allocatable,dimension(:,:):: peakX,DcX,dynX
      real,allocatable,dimension(:,:):: peakZ,DcZ,dynZ
      
      real:: dip
      real,parameter:: pi=3.1415926535	 
      
      CONTAINS
      
      FUNCTION normstress(j)
      USE medium_com, only: constnormstress
      IMPLICIT NONE
      real:: normstress
      integer:: j
#if defined FSPACE
      normstress=constnormstress
#else
#if defined DIPSLIP
      normstress=max(1.e5,8520.*dh*real(nzt-nfs-j)*sin(dip/180.*pi))
      !normstress=min(18.*dh*real(nzt-nfs-j)*sin(dip/180.*pi)/1.e3,100.);normstress=1.e6*max(1.,normstress)
#else
      normstress=max(1.e5,16200.*dh*real(nzt-nfs-j)*sin(dip/180.*pi))
      !normstress=min(18.*dh*real(nzt-nfs-j)*sin(dip/180.*pi)/1.e3,100.);normstress=1.e6*max(1.,normstress)
       !normstress=min(18.*dh*real(nzt-nfs-j)*sin(dip/180.*pi)/1.e3,60.);normstress=1.e6*max(1.,normstress)
#endif
#endif
      END FUNCTION
      
    END MODULE

    MODULE source_com
      REAL,ALLOCATABLE,DIMENSION(:,:):: ruptime,rise,slipZ,schangeX,schangeZ,sliptime,slipX
      real    :: output_param(6)
      integer :: ioutput
      integer:: Nstations
      integer,allocatable:: staX(:), staY(:), staZ(:)
      real,allocatable :: seisU(:), seisV(:), seisW(:)
      real,allocatable :: seissurfU(:,:), seissurfV(:,:), seissurfW(:,:)
	  real :: waveT
	  real :: Eg, Er
	 ! REAL,allocatable :: waveU(:,:,:),waveV(:,:,:), waveW(:,:,:)
	  
    END MODULE
    
    MODULE displt_com
      real,allocatable,dimension(:,:,:):: u1,v1,w1
    END MODULE

    MODULE strfld_com
      real,allocatable,dimension(:,:,:):: xx,yy,zz,xy,yz,xz
    END MODULE

    MODULE traction_com
      real,allocatable,dimension(:,:):: tx,tz,v1t,avdx,avdz,RFx,RFz, T0XI, T0ZI  ! x and z component of traction, v component of velocity at fault
      real,allocatable,dimension(:,:):: au1,av1,aw1  !attenuation of velocities near fault
      real :: damp_s
    END MODULE

    MODULE SlipRates_com
      INTEGER nSR,NL,NW
      REAL dL,dW,dtseis,SRdur
      REAL M0,Mw
      REAL,allocatable,dimension(:):: MSRX,MSRZ,MomentRate
      REAL,allocatable,dimension(:,:):: muSource
    END MODULE
	
	MODULE PostSeismic_com
	
      REAL,allocatable,dimension(:):: MSX, MSZ, TS
	  REAL,allocatable,dimension(:,:):: gpssyntN, gpssyntE, gpssyntZ
	  REAL,allocatable,dimension(:,:):: gpsrealN, gpsrealE, gpsrealZ, gpsrealT
	  REAL,allocatable,dimension(:,:):: gpsgfN, gpsgfE, gpsgfZ, gpssigma
	  integer,allocatable,dimension(:)::gpsrealTN
	  REAL :: T1S,T2S,SigmaGPS,SigmaGPS2, gpsweight, VRGPS
	  integer :: NTS, NGPS, igps, NTSrv
	  
	END MODULE

    SUBROUTINE fd3d_init()
      USE medium_com
      USE fd3dparam_com
      USE friction_com
      USE source_com
      USE pml_com
      USE traction_com
      USE SlipRates_com
#if defined FVW
      USE RATESTATE,only:gkoef,MDIS,NDIS,slipOUT,f0PS,baPS,fwPS,vwPS,aPS,dcPS,snPS
#endif
      IMPLICIT NONE
      integer nxtT, nytT, nztT
      integer i
      real pml_vp,pml_fact  

!--------------------
! Read the input file
!--------------------

      write(*,*)'Reading FD3D parameters...'
      open(11, file='inputfd3d.dat', status='old')
      read(11,*) nxtT,nytT,nztT
      read(11,*) dh
      read(11,*) ntfd
      read(11,*) dt
      read(11,*) dip
      read(11,*) nabc, pml_vp,pml_fact   !(pml_fact=-(N+1)*log(0.001), see Komatitsch and Martin, 2007, Geophysics 72)
      read(11,*) damp_s

#if defined FSPACE
	  nfs=nabc
      read(11,*) constnormstress
#else
      nfs=2 ! Number of layers above free surface 
#endif	
      nxt=nxtT+2*nabc
      nyt=nytT+nabc
      nzt=nztT+nabc+nfs
      nysc=nyt
      omegaM_pml=pml_fact*pml_vp/(2.*dh*(nabc-1))

      read(11,*) Nstations
      if(Nstations>0) then
        allocate(staX(Nstations),staY(Nstations),staZ(Nstations),seisU(Nstations),seisV(Nstations),seisW(Nstations),seissurfU(nxt,nyt),seissurfV(nxt,nyt),seissurfW(nxt,nyt))
        do i=1,Nstations
          read(11,*) staX(i),staY(i),staZ(i)
          staX(i)=nabc+staX(i)
          staY(i)=nabc+staY(i)
          staZ(i)=nabc+staZ(i)
        enddo
      endif 
!	  read(11,*) waveT
      close(11)
 
!----------------------------
! Allocate FD module arrays
!----------------------------
      allocate(lam1(nxt,nyt,nzt),mu1(nxt,nyt,nzt),d1(nxt,nyt,nzt), vpls(nxt,nzt))
      allocate(uZ(nxt,nzt),wX(nxt,nzt),tabsX(nxt,nzt),tabsZ(nxt,nzt))
      
#if defined FVW
      allocate(dc(nxt,nzt),dcX(nxt,nzt),dcZ(nxt,nzt))
      allocate(striniZ(nxt,nzt),striniX(nxt,nzt), a(nxt,nzt), ba(nxt,nzt))
	  
      allocate(psi(nxt,nzt),vw(nxt,nzt),t0X(nxt,nzt),t0Z(nxt,nzt))
      allocate(aX(nxt,nzt),baX(nxt,nzt),psiX(nxt,nzt),vwX(nxt,nzt))
      allocate(aZ(nxt,nzt),baZ(nxt,nzt),psiZ(nxt,nzt),vwZ(nxt,nzt))
	  allocate(f0(nxt,nzt),fw(nxt,nzt),Sn(nxt,nzt))
	  allocate(f0X(nxt,nzt),fwX(nxt,nzt),SnX(nxt,nzt))
	  allocate(f0Z(nxt,nzt),fwZ(nxt,nzt),SnZ(nxt,nzt))
	  allocate(uini(nxt,nzt),wini(nxt,nzt),t0Xi(nxt,nzt),t0Zi(nxt,nzt))
	  perturb=0.
	  gkoef=4
	  MDIS=((nxt-2*nabc-1)/gkoef)+1    
	  NDIS=((nzt-nfs-nabc-1)/gkoef)+1 
	  allocate(slipOUT(2,MDIS*NDIS))
	  allocate(f0PS(MDIS,NDIS),baPS(MDIS,NDIS),fwPS(MDIS,NDIS),vwPS(MDIS,NDIS),aPS(MDIS,NDIS),dcPS(MDIS,NDIS),snPS(MDIS,NDIS))
#else

      allocate(striniZ(nxt,nzt),striniX(nxt,nzt),peak_xz(nxt,nzt),Dc(nxt,nzt),dyn_xz(nxt,nzt))
      allocate(peakX(nxt,nzt),T0X(nxt,nzt),DcX(nxt,nzt),dynX(nxt,nzt))
      allocate(peakZ(nxt,nzt),T0Z(nxt,nzt),DcZ(nxt,nzt),dynZ(nxt,nzt))

#endif
      
      allocate(ruptime(nxt,nzt),slipZ(nxt,nzt),slipX(nxt,nzt),rise(nxt,nzt),schangeZ(nxt,nzt),schangeX(nxt,nzt),sliptime(nxt,nzt))
      allocate(muSource(nxt,nzt),coh(nxt,nzt))

      striniX=0.; striniZ=0.; peak_xz=0.; Dc=0.

!------------------------------------------------------------
! Read the velocity model
! Be careful: the velocity model for the FD is upside down
!------------------------------------------------------------
      CALL readcrustalmodel(dip)

    END SUBROUTINE


    SUBROUTINE readcrustalmodel(dip)
    USE medium_com
    USE fd3dparam_com
    USE pml_com
    use mod_pgamisf, only : v30 !obtain v30 for gmpe calculation
    use SlipRates_com
    IMPLICIT NONE
    real*8,parameter:: PI=3.1415926535
    real dip, d_zone, vlow_zone,vlowMax_zone, vlowMin_zone 
    real    :: vpe(2),vse(2),den(2),CFL,dum,dd,vpp,vss
    real,allocatable,dimension(:):: vp,vs,depth,rho
    INTEGER ndepth,j,k,i

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
    do k=nzt,1,-1
      dum=(dh*real(nzt-nfs-k)+dh/2.)*sin(dip/180.d0*PI)
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
      if (dum<=30.) then
        v30=vss !k from surface, dum from "0", rewrites v30 if there is layer<=30
      else
        if(k==nzt) v30=vss ! take first value
      endif
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
	  vpls(1:nxt,k)=2/(dd*vss)
    enddo
    mu_mean = (mu_mean/nzt)
!   write(*,*)mu_mean
    deallocate(depth,vp,vs,rho)

!damaged fault zone
    muSource(:,:)=mu1(:,nysc,:)
 !   d_zone=500 ! fault zone width
 !   vlowMax_zone=0.3 !velocity reduction at the fault
 !   vlowMin_zone=0.1 !velocity reduction at the fault
 !   vlowMax_zone=0.0 !no velocity reduction
 !   vlowMin_zone=0.0 !no velocity reduction
!	  j=0
 !   do while (j*dh<=d_zone)
 !     do k=1,nzt-2
 !       do i=1,nxt
!         lam1(i,nyt-j,k)=(1-vlow_zone)**2*lam1(i,nyt-j,k)
  !        vlow_zone= vlowMax_zone - (vlowMax_zone-vlowMin_zone)*(j*dh/d_zone)**2
  !        mu1(i,nyt-j,k)=(1-vlow_zone)**2*mu1(i,nyt-j,k)
  !      enddo
  !    enddo
!	    j=j+1
!	  enddo
	  
!-------------------------------------------------------
!     Make sure the simulation is numerically stable
!     CFL < 0.25
!-------------------------------------------------------
    CFL = vpe(2)*dt/dh
    if (CFL.gt.0.25) then
      print *,'Your simulation is numerically unstable', CFL
      stop
    endif
    END SUBROUTINE
    
#if defined FVW
    
    SUBROUTINE forwardspecialTPV103()
!   Setup of dynamic parameters for TPV103 benchmark
    USE inversion_com
    USE fd3dparam_com
    USE friction_com
    USE pml_com
    IMPLICIT NONE

    integer i,j,k
    real dum
    real wfx,wfz, trans, fringe, f0ini, v0ini, a0ini, bini, d0h, fwini, vw0ini, vwdeltaini, vini,RRini
    real Tini, Snini, Psiini, BX1, BZ1, hx, hz, hx0t, hz0t

    open(244,FILE='scecmodel.dat')
    read (244,*) wfx,wfz, hx0t, hz0t, trans, fringe  !polovina sirky zlomu v kmetrech,poloha hypocentra, okraj v metrech
    read (244,*) f0ini ! 
    read (244,*) v0ini ! 
    read (244,*) a0ini !
    read (244,*) bini !
    read (244,*) d0h ! 
    read (244,*) fwini ! 
    read (244,*) vw0ini !
    read (244,*) vwdeltaini !
    read (244,*) vini !
    read (244,*) Tini ! 
    read (244,*) Snini ! 
    read (244,*) perturb
    read (244,*) RRini
    read (244,*) TT2
    striniZ=0.
    close(244)
    
    RR2 = RRini*RRini
    Sn=Snini
    uini=vini
    wini=0.
    f0=f0ini
    v0=v0ini
    strinixI=Tini
    striniX=Tini
    fw=fwini
    
    hx0 = hx0t + fringe + real(nabc)*dh
    hz0 = hz0t + fringe + real(nabc)*dh
    print*,hx0,hz0
    striniZ=0.
    do k=1,nzt
      do i = 1,nxt

        Dc(i,k)=d0h
        hx = real(i)*dh
        hz = real(k)*dh
        if ( abs(hx - hx0) >= wfx + trans) then
          BX1 = 0.
        elseif (( abs(hx - hx0) > wfx ) .AND. (abs(hx - hx0) < wfx + trans)) then
          BX1 = 0.5*(1+tanh(trans/(abs(hx - hx0) - wfx - trans) + trans/(abs(hx - hx0) - wfx)))
        else
          BX1 = 1.
        endif
        
        if ( abs(hz - hz0) >= wfz + trans) then
          BZ1 = 0.
        elseif (( abs(hz - hz0) > wfz ) .AND. (abs(hz - hz0) < wfz + trans)) then
          BZ1 = 0.5*(1+tanh(trans/(abs(hz - hz0) - wfz - trans) + trans/(abs(hz - hz0) - wfz)))
        else
          BZ1 = 1.
        endif
        
        !print *, BX*BZ
        a(i,k) = a0ini + a0ini*(1.-BX1*BZ1)
	ba(i,k)=bini-a(i,k)
        vw(i,k) = vw0ini + vwdeltaini*(1.-BX1*BZ1) 
        !strinix(i,k) = Sn*(f0-(b-a(i,k))*log(2*vini/v0ini))
        psi(i,k) = a(i,k)*(log((2*v0ini/(2*vini))) + log(sinh(striniX(i,k)/(a(i,k)*Sn(i,k)))))

      enddo
    enddo
   
    do k=1,nzt
      do i = 1,nxt

        hx = real(i)*dh +dh/2.
        hz = real(k)*dh +dh/2.
        if ( abs(hx - hx0) >= wfx + trans) then
          BX1 = 0.
        elseif (( abs(hx - hx0) > wfx ) .AND. (abs(hx - hx0) < wfx + trans)) then
          BX1 = 0.5*(1+tanh(trans/(abs(hx - hx0) - wfx - trans) + trans/(abs(hx - hx0) - wfx)))
        else
          BX1 = 1.
        endif
        
        if ( abs(hz - hz0) >= wfz + trans) then
          BZ1 = 0.
        elseif (( abs(hz - hz0) > wfz ) .AND. (abs(hz - hz0) < wfz + trans)) then
          BZ1 = 0.5*(1+tanh(trans/(abs(hz - hz0) - wfz - trans) + trans/(abs(hz - hz0) - wfz)))	
        else
          BZ1 = 1.
        endif
        
        !print *, BX*BZ
        aX(i,k) = a0ini + a0ini*(1.-BX1*BZ1)
        vwX(i,k) = vw0ini + vwdeltaini*(1.-BX1*BZ1) 
        !strinix(i,k) = Sn*(f0-(b-a(i,k))*log(2*vini/v0ini))
        psiX(i,k)=aX(i,k)*(log((2*v0/(2*uini(i,k)))) + log(sinh(striniX(i,k)/(aX(i,k)*Sn(i,k)))))
        baX(i,k)=bini-aX(i,k)
      enddo
    enddo

    END SUBROUTINE
	
	SUBROUTINE forwardspecialTPV104()
!   Setup of dynamic parameters for TPV104 benchmark
    USE inversion_com
    USE fd3dparam_com
    USE friction_com
    USE pml_com
    IMPLICIT NONE

    integer i,j,k
    real dum
    real wfx,wfz, trans, fringe, f0ini, v0ini, a0ini, bini, d0h, fwini, vw0ini, vwdeltaini, vini,RRini
    real Tini, Snini, Psiini, BX1, BZ1, hx, hz, hx0t, hz0t

    open(244,FILE='scecmodel.dat')
    read (244,*) wfx,wfz, hx0t, hz0t, trans, fringe  !polovina sirky zlomu v kmetrech,poloha hypocentra, okraj v metrech
    read (244,*) f0ini ! 
    read (244,*) v0ini ! 
    read (244,*) a0ini !
    read (244,*) bini !
    read (244,*) d0h ! 
    read (244,*) fwini ! 
    read (244,*) vw0ini !
    read (244,*) vwdeltaini !
    read (244,*) vini !
    read (244,*) Tini ! 
    read (244,*) Snini ! 
    read (244,*) perturb
    read (244,*) RRini
    read (244,*) TT2
    striniZ=0.
    close(244)
    
    RR2 = RRini*RRini
    Sn=Snini
    uini=vini
    wini=0.
    f0=f0ini
    v0=v0ini
    strinixI=Tini
    fw=fwini
    
    hx0 = hx0t + fringe + real(nabc)*dh
    hz0 = hz0t + fringe + real(nabc)*dh
    print*,hx0,hz0
    striniZ=0.
    do k=1,nzt
      do i = 1,nxt
        striniX(i,k)=Tini
        Dc(i,k)=d0h
        hx = real(i)*dh
        hz = real(k)*dh
        if ( abs(hx - hx0) >= wfx + trans) then
          BX1 = 0.
        elseif (( abs(hx - hx0) > wfx ) .AND. (abs(hx - hx0) < wfx + trans)) then
          BX1 = 0.5*(1+tanh(trans/(abs(hx - hx0) - wfx - trans) + trans/(abs(hx - hx0) - wfx)))
        else
          BX1 = 1.
        endif
        
        if ( -(hz - hz0) >= wfz + trans) then
          BZ1 = 0.
        elseif (( -(hz - hz0) > wfz ) .AND. (-(hz - hz0) < wfz + trans)) then
          BZ1 = 0.5*(1+tanh(trans/(abs(hz - hz0) - wfz - trans) + trans/(abs(hz - hz0) - wfz)))
        else
          BZ1 = 1.
        endif
        
        !print *, BX*BZ
        a(i,k) = a0ini + a0ini*(1.-BX1*BZ1)
        ba(i,k) = bini -a(i,k)
        vw(i,k) = vw0ini + vwdeltaini*(1.-BX1*BZ1) 
        !strinix(i,k) = Sn*(f0-(b-a(i,k))*log(2*vini/v0ini))
        psi(i,k) = a(i,k)*(log((2*v0ini/(2*vini))) + log(sinh(striniX(i,k)/(a(i,k)*Sn(i,k)))))

      enddo
    enddo
   
    do k=1,nzt
      do i = 1,nxt

        hx = real(i)*dh +dh/2.
        hz = real(k)*dh +dh/2.
        if ( abs(hx - hx0) >= wfx + trans) then
          BX1 = 0.
        elseif (( abs(hx - hx0) > wfx ) .AND. (abs(hx - hx0) < wfx + trans)) then
          BX1 = 0.5*(1+tanh(trans/(abs(hx - hx0) - wfx - trans) + trans/(abs(hx - hx0) - wfx)))
        else
          BX1 = 1.
        endif
        
        if ( -(hz - hz0) >= wfz + trans) then
          BZ1 = 0.
        elseif (( -(hz - hz0) > wfz ) .AND. (-(hz - hz0) < wfz + trans)) then
          BZ1 = 0.5*(1+tanh(trans/(abs(hz - hz0) - wfz - trans) + trans/(abs(hz - hz0) - wfz)))	
        else
          BZ1 = 1.
        endif
        
        !print *, BX*BZ
        aX(i,k) = a0ini + a0ini*(1.-BX1*BZ1)
        baX(i,k) = bini - aX(i,k)
        vwX(i,k) = vw0ini + vwdeltaini*(1.-BX1*BZ1) 
        !strinix(i,k) = Sn*(f0-(b-a(i,k))*log(2*vini/v0ini))
        psiX(i,k)=aX(i,k)*(log((2*v0/(2*uini(i,k)))) + log(sinh(striniX(i,k)/(aX(i,k)*Sn(i,k)))))

      enddo
    enddo

    END SUBROUTINE

#else

    SUBROUTINE forwardspecialTPV5()
!   Setup of dynamic parameters for TPV5 benchmark
    USE inversion_com
    USE fd3dparam_com
    USE friction_com
    USE medium_com
    IMPLICIT NONE
    REAL,PARAMETER:: x0=15.e3,z0=6.e3,a=10.e3,b=4.e3,phi=0.,xn=17.e3,zn=6.e3,rn=1.5e3
    real hx0, hz0, h1x0, h1z0, h2x0, h2z0, hdelta, T0, T0n, sn, mus, mud, d0h
    real x,z,rr,DL,DW, muso, T0h1, T0h2
    integer i,j,k,no0,j2
    real dum
    real d_zone, vlow_zone

    open(244,FILE='scecmodel.dat')
    read (244,*) no0, hx0, hz0, h1x0, h1z0, h2x0, h2z0  !okraj,stred nukleace, stred leve heterogenity, stred prave heterogenity
    read (244,*) hdelta ! sirka nukleacni zony
    read (244,*) T0, T0n, T0h1, T0h2 ! predpeti, predpeti v nukleacni zone
    read (244,*) sn ! normalove napeti
    read (244,*) mus, muso !staticke treni, staticke treni na okraji
    read (244,*) mud !dynamicke treni
    read (244,*) d0h ! kriticky slip
    read (244,*) d_zone, vlow_zone ! sirka zlomove zony, pokles elastickeho modulu ve zlomove zone
    close(244)
    striniZ=0.
    coh=0.
    do k=1,nzt-2
        do i = 1,nxt
        striniX(i,k)=T0

          if ((((real(i)-1.)*dh-hx0>= -hdelta/2.0) .and. ((real(i)-1.)*dh-hx0 <= hdelta/2.0)) &
            .and. (((real(k)-1.)*dh-hz0>= -hdelta/2.0) .and. ((real(k)-1.)*dh-hz0 <= hdelta/2.0))) then
            striniX(i,k)=T0n
          endif
          if ((((real(i)-1.)*dh-h1x0>= -hdelta/2.0) .and. ((real(i)-1.)*dh-h1x0 <= hdelta/2.0)) &
            .and. (((real(k)-1.)*dh-h1z0>= -hdelta/2.0) .and. ((real(k)-1.)*dh-h1z0 <= hdelta/2.0))) then
            striniX(i,k)=T0h1
          endif

          if ((((real(i)-1.)*dh-h2x0>= -hdelta/2.0) .and. ((real(i)-1.)*dh-h2x0 <= hdelta/2.0)) &
            .and. (((real(k)-1.)*dh-h2z0>= -hdelta/2.0) .and. ((real(k)-1.)*dh-h2z0 <= hdelta/2.0))) then
            striniX(i,k)=T0h2
          endif

          peak_xz(i,k)=sn*muso
          dyn_xz(i,k) = sn*mud
          Dc(i,k)=d0h
        enddo
    enddo
    
    do k=no0+1,nzt-2
      do i = no0+1,nxt-no0
        peak_xz(i,k)=sn*mus
      enddo
    enddo
    
    do i=1,nxt
      peak_xz(i,nzt-1)=peak_xz(i,nzt-2)
      dyn_xz(i,nzt-1)=dyn_xz(i,nzt-2)
      Dc(i,nzt-1)=Dc(i,nzt-2)
      striniX(i,nzt-1)=striniX(i,nzt-2)
      peak_xz(i,nzt)=peak_xz(i,nzt-3)
      dyn_xz(i,nzt)=dyn_xz(i,nzt-3)
      Dc(i,nzt)=Dc(i,nzt-3)
      striniX(i,nzt)=striniX(i,nzt-3)
    enddo

    END SUBROUTINE

    SUBROUTINE forwardspecialTPV9()
!   Setup of dynamic parameters for TPV9 benchmark
    USE inversion_com
    USE fd3dparam_com
    USE friction_com
    USE medium_com
    USE strfld_com
    IMPLICIT NONE
    REAL,PARAMETER:: x0=15.e3,z0=6.e3,a=10.e3,b=4.e3,phi=0.,xn=17.e3,zn=6.e3,rn=1.5e3
    real hx0, hz0, h1x0, h1z0, h2x0, h2z0, hdelta, T0, T0n, sn, mus, mud, d0h
    real x,z,rr,DL,DW, muso, T0h1, T0h2,sigma_n
    integer i,j,k,no0,j2
    real dum
    real d_zone, vlow_zone

    open(244,FILE='scecmodel.dat')
    read (244,*) no0, hx0, hz0, h1x0, h1z0, h2x0, h2z0  !okraj,stred nukleace, stred leve heterogenity, stred prave heterogenity
    read (244,*) hdelta ! sirka nukleacni zony
    read (244,*) T0, T0n, T0h1, T0h2 ! predpeti, predpeti v nukleacni zone
    read (244,*) sn ! normalove napeti
    read (244,*) mus, muso !staticke treni, staticke treni na okraji
    read (244,*) mud !dynamicke treni
    read (244,*) d0h ! kriticky slip
    read (244,*) d_zone, vlow_zone ! sirka zlomove zony, pokles elastickeho modulu ve zlomove zone
    close(244)

    coh=1.e6

    striniX=0.
    do k=2,nzt-2
      do i = 1,nxt
        sigma_n=sn*dh*(real(nzt-2-k))
        peak_xz(i,k)=sigma_n*muso
        dyn_xz(i,k) = sigma_n*mud
        striniZ(i,k)=T0*sigma_n
        
       if ((((real(i)-1.)*dh-hx0>= -hdelta/2.0) .and. ((real(i)-1.)*dh-hx0 <= hdelta/2.0)) &
         .and. (((real(k)-1.)*dh-hz0>= -hdelta/2.0) .and. ((real(k)-1.)*dh-hz0 <= hdelta/2.0))) then
!       if ((((real(i))*dh-hx0>= -hdelta/2.0) .and. ((real(i))*dh-hx0 <= hdelta/2.0)) &
!         .and. (((real(k))*dh-hz0>= -hdelta/2.0) .and. ((real(k))*dh-hz0 <= hdelta/2.0))) then
           striniZ(i,k)=T0n*mus*sigma_n+coh(i,k)
       endif

       Dc(i,k)=d0h
      enddo
    enddo

    do k=no0+1,nzt-2
      do i = no0+1,nxt-no0
        sigma_n=sn*dh*(real(nzt-2-k))
        peak_xz(i,k)=sigma_n*mus
      enddo
    enddo

    do i = 1,nxt
      peak_xz(i,nzt-2)=0.25*peak_xz(i,nzt-3)
      dyn_xz(i,nzt-2) = 0.25*dyn_xz(i,nzt-3)
      striniZ(i,nzt-2)=0.25*striniZ(i,nzt-3)
      Dc(i,nzt-2)=d0h
    enddo

    do i=1,nxt
      peak_xz(i,nzt-1)=peak_xz(i,nzt-2)
      dyn_xz(i,nzt-1)=dyn_xz(i,nzt-2)
      Dc(i,nzt-1)=Dc(i,nzt-2)
      striniZ(i,nzt-1)=-striniZ(i,nzt-2)
      peak_xz(i,nzt)=peak_xz(i,nzt-3)
      dyn_xz(i,nzt)=dyn_xz(i,nzt-3)
      Dc(i,nzt)=Dc(i,nzt-3)
      striniZ(i,nzt)=-striniZ(i,nzt-3)
    enddo

    END SUBROUTINE
    
    SUBROUTINE forwardspecialTPV8()
!   Setup of dynamic parameters for TPV8 benchmark
    
    USE inversion_com
    USE fd3dparam_com
    USE friction_com
    USE medium_com
    USE strfld_com
    IMPLICIT NONE
    REAL,PARAMETER:: x0=15.e3,z0=6.e3,a=10.e3,b=4.e3,phi=0.,xn=17.e3,zn=6.e3,rn=1.5e3
    real hx0, hz0, h1x0, h1z0, h2x0, h2z0, hdelta, T0, T0n, sn, mus, mud, d0h
    real x,z,rr,DL,DW, muso, T0h1, T0h2,sigma_n
    integer i,j,k,no0,j2
    real dum
    real d_zone, vlow_zone

    open(244,FILE='scecmodel.dat')
    read (244,*) no0, hx0, hz0, h1x0, h1z0, h2x0, h2z0  !okraj,stred nukleace, stred leve heterogenity, stred prave heterogenity
    read (244,*) hdelta ! sirka nukleacni zony
    read (244,*) T0, T0n, T0h1, T0h2 ! predpeti, predpeti v nukleacni zone
    read (244,*) sn ! normalove napeti
    read (244,*) mus, muso !staticke treni, staticke treni na okraji
    read (244,*) mud !dynamicke treni
    read (244,*) d0h ! kriticky slip
    read (244,*) d_zone, vlow_zone ! sirka zlomove zony, pokles elastickeho modulu ve zlomove zone
    close(244)

    coh=1.e6

    striniZ=0.
    do k=2,nzt-2
      do i = 1,nxt
        sigma_n=sn*dh*(real(nzt-2-k)+0.5)
        peak_xz(i,k)=sigma_n*muso
        dyn_xz(i,k) = sigma_n*mud
        striniX(i,k)=T0*sigma_n
            
        if ((((real(i)-1.)*dh-hx0>= -hdelta/2.0) .and. ((real(i)-1.)*dh-hx0 <= hdelta/2.0)) &
          .and. (((real(k)-1.)*dh-hz0>= -hdelta/2.0) .and. ((real(k)-1.)*dh-hz0 <= hdelta/2.0))) then
            striniX(i,k)=T0n*mus*sigma_n+coh(i,k)
        endif

        Dc(i,k)=d0h
      enddo
    enddo

    do k=no0+1,nzt-2
      do i = no0+1,nxt-no0
        sigma_n=sn*dh*(real(nzt-2-k)+0.5)
        peak_xz(i,k)=sigma_n*mus
      enddo
    enddo

    do i=1,nxt
      peak_xz(i,nzt-1)=peak_xz(i,nzt-2)
      dyn_xz(i,nzt-1)=dyn_xz(i,nzt-2)
      Dc(i,nzt-1)=Dc(i,nzt-2)
      strinix(i,nzt-1)=-strinix(i,nzt-2)
      peak_xz(i,nzt)=peak_xz(i,nzt-3)
      dyn_xz(i,nzt)=dyn_xz(i,nzt-3)
      Dc(i,nzt)=Dc(i,nzt-3)
      strinix(i,nzt)=-strinix(i,nzt-3)
    enddo
    END SUBROUTINE
#endif

#if defined FVW
    subroutine inversion_modeltofd3d() ! inversion from model control points to fd3d grid
    use inversion_com
    use fd3dparam_com
    use friction_com
    use pml_com
	use medium_com
    implicit none
    real,dimension(:),allocatable:: xintrpl,yintrpl,xnew,ynew

    integer i,k,ii,kk
    real DW,DL,dum,xs,zs,t,u, hx, hz, rr, FXZ
	hx0=nucl(1)
	hz0=nucl(2)
	RR2=nucl(3)
	TT2=nucl(4)
	perturb=nucl(5)

    v0 =1.e-6
    coh=0.	
	wini=0.
! Bilinear interpolation
    DL=dh*(nxt-2*nabc-1)/real(NLI-1)
    DW=dh*(nzt-nfs-nabc-1)/real(NWI-1)
    do k=nabc+1,nzt-nfs
      ZS=dh*(k-1-nabc)
      kk=min(NWI-1,int(ZS/DW)+1)
      u=min(1.,(ZS-DW*(kk-1))/DW)
      do i=nabc+1,nxt-nabc
        XS=dh*(i-1-nabc)
        ii=min(NLI-1,int(XS/DL)+1)
        t=min(1.,(XS-DL*(ii-1))/DL) 
        Dc(i,k)     = (1.-t)*(1.-u)*DcI(ii,kk)+t*(1.-u)*DcI(ii+1,kk)+t*u*DcI(ii+1,kk+1)+(1.-t)*u*DcI(ii,kk+1)
        Sn(i,k)     = normstress(k)
        a(i,k)      = (1.-t)*(1.-u)*aI(ii,kk)+t*(1.-u)*aI(ii+1,kk)+t*u*aI(ii+1,kk+1)+(1.-t)*u*aI(ii,kk+1)
        ba(i,k)     = (1.-t)*(1.-u)*baI(ii,kk)+t*(1.-u)*baI(ii+1,kk)+t*u*baI(ii+1,kk+1)+(1.-t)*u*baI(ii,kk+1)
        !psi(i,k)    = (1.-t)*(1.-u)*psiI(ii,kk)+t*(1.-u)*psiI(ii+1,kk)+t*u*psiI(ii+1,kk+1)+(1.-t)*u*psiI(ii,kk+1)
        f0(i,k)     = (1.-t)*(1.-u)*f0I(ii,kk)+t*(1.-u)*f0I(ii+1,kk)+t*u*f0I(ii+1,kk+1)+(1.-t)*u*f0I(ii,kk+1)
        fw(i,k)     = (1.-t)*(1.-u)*fwI(ii,kk)+t*(1.-u)*fwI(ii+1,kk)+t*u*fwI(ii+1,kk+1)+(1.-t)*u*fwI(ii,kk+1)
        vw(i,k)     = (1.-t)*(1.-u)*vwI(ii,kk)+t*(1.-u)*vwI(ii+1,kk)+t*u*vwI(ii+1,kk+1)+(1.-t)*u*vwI(ii,kk+1)

#if defined DIPSLIP
		wini(i,k)   = (1.-t)*(1.-u)*viniI(ii,kk)+t*(1.-u)*viniI(ii+1,kk)+t*u*viniI(ii+1,kk+1)+(1.-t)*u*viniI(ii,kk+1)
        striniZ(i,k)= (1.-t)*(1.-u)*T0I(ii,kk)+t*(1.-u)*T0I(ii+1,kk)+t*u*T0I(ii+1,kk+1)+(1.-t)*u*T0I(ii,kk+1)
        psi(i,k) = a(i,k)*(log((2*v0/(2*wini(i,k)))) + log(sinh(striniZ(i,k)/(a(i,k)*Sn(i,k)))))
        striniX(i,k)= 0.
#else
		uini(i,k)   = (1.-t)*(1.-u)*viniI(ii,kk)+t*(1.-u)*viniI(ii+1,kk)+t*u*viniI(ii+1,kk+1)+(1.-t)*u*viniI(ii,kk+1)
        striniX(i,k)= (1.-t)*(1.-u)*T0I(ii,kk)+t*(1.-u)*T0I(ii+1,kk)+t*u*T0I(ii+1,kk+1)+(1.-t)*u*T0I(ii,kk+1)
		psi(i,k) = a(i,k)*(log((2*v0/(2*uini(i,k)))) + log(sinh(striniX(i,k)/(a(i,k)*Sn(i,k)))))
        striniZ(i,k)= 0.
#endif
      enddo
    enddo
	k=nabc

      do i=nabc+1,nxt-nabc
	  
        Dc(i,k)     = Dc(i,k+1)
        Sn(i,k)     = Sn(i,k+1)
        a(i,k)      = a(i,k+1)
        ba(i,k)     = ba(i,k+1)
        psi(i,k)    = psi(i,k+1)
        f0(i,k)     = f0(i,k+1)
        fw(i,k)     = fw(i,k+1)
        vw(i,k)     = vw(i,k+1)
		wini(i,k)   = wini(i,k+1)
		uini(i,k)   = uini(i,k+1)
#if defined DIPSLIP
        striniZ(i,k)= striniZ(i,k+1)
        !psi(i,k) = a(i,k)*(log((2*v0/(2*wini))) + log(sinh(striniZ(i,k)/(a(i,k)*Sn(i,k)))))
        striniX(i,k)= 0.
#else
        striniX(i,k)= striniX(i,k+1)
        !psi(i,k) = a(i,k)*(log((2*v0/(2*uini(i,k)))) + log(sinh(striniX(i,k)/(a(i,k)*Sn(i,k)))))
        striniZ(i,k)= 0.
#endif
      enddo
	
	do k=nabc,nzt-nfs
          i=nabc
          Dc(i,k)     = Dc(i+1,k)
          Sn(i,k)     = Sn(i+1,k)
          a(i,k)      = a(i+1,k)
          ba(i,k)     = ba(i+1,k)
          psi(i,k)    = psi(i+1,k)
          f0(i,k)     = f0(i+1,k)
          fw(i,k)     = fw(i+1,k)
          vw(i,k)     = vw(i+1,k)
		  wini(i,k)   = wini(i+1,k)
		  uini(i,k)   = uini(i+1,k)
#if defined DIPSLIP
          striniZ(i,k)= striniZ(i+1,k)
          !psi(i,k) = a(i,k)*(log((2*v0/(2*wini))) + log(sinh(striniZ(i,k)/(a(i,k)*Sn(i,k)))))
          striniX(i,k)= 0.
#else
          striniX(i,k)= striniX(i+1,k)
          !psi(i,k) = a(i,k)*(log((2*v0/(2*uini(i,k)))) + log(sinh(striniX(i,k)/(a(i,k)*Sn(i,k)))))
          striniZ(i,k)= 0.
#endif
    enddo	
    	
	hx0=hx0+real(nabc)*dh
	hz0=hz0+real(nabc)*dh
	
    END SUBROUTINE
	
    SUBROUTINE readinversionresult()
    USE inversion_com
    USE fd3dparam_com
    USE friction_com
    IMPLICIT NONE
    real x,z,rr,DL,DW
    integer i,j
    real dum

    open(244,FILE='forwardmodel.dat')
    read(244,*)dum,dum, nucl(1:5), T0I(:,:),aI(:,:), baI(:,:), psiI(:,:), f0I(:,:), fwI(:,:), DcI(:,:), vwI(:,:), viniI(:,:)
    close(244)
    
    CALL inversion_modeltofd3d()
    
    END SUBROUTINE

#else

    SUBROUTINE readinversionresult()
    USE inversion_com
    USE fd3dparam_com
    USE friction_com
    IMPLICIT NONE
    real x,z,rr,DL,DW
    integer i,j
    real dum
    
    open(244,FILE='forwardmodel.dat')
    read(244,*)dum,dum,T0I(:,:),TsI(:,:),DcI(:,:)
    close(244)
    
    CALL inversion_modeltofd3d()
    
    END SUBROUTINE
    SUBROUTINE inversion_modeltofd3d() ! inversion from model control points to fd3d grid
    USE inversion_com
    USE fd3dparam_com
    USE friction_com
    USE pml_com
    IMPLICIT NONE
    real,dimension(:),allocatable:: xintrpl,yintrpl,xnew,ynew
    real,parameter:: dyn=0.4  !dynamic friction coefficient
    integer i,k,ii,kk
    real DW,DL,dum,xs,zs,t,u

! Bilinear interpolation
    DL=dh*(nxt-2*nabc-1)/real(NLI-1)
    DW=dh*(nzt-nfs-nabc-1)/real(NWI-1)
    do k=nabc+1,nzt-nfs
      ZS=dh*(k-1-nabc)
      kk=min(NWI-1,int(ZS/DW)+1)
      u=min(1.,(ZS-DW*(kk-1))/DW)
      do i=nabc+1,nxt-nabc
        XS=dh*(i-1-nabc)
        ii=min(NLI-1,int(XS/DL)+1)
        t=min(1.,(XS-DL*(ii-1))/DL) 
        Dc(i,k)     =(1.-t)*(1.-u)*DcI(ii,kk)+t*(1.-u)*DcI(ii+1,kk)+t*u*DcI(ii+1,kk+1)+(1.-t)*u*DcI(ii,kk+1)
        dyn_xz(i,k)=dyn*normstress(k)
        peak_xz(i,k)=((1.-t)*(1.-u)*TsI(ii,kk)+t*(1.-u)*TsI(ii+1,kk)+t*u*TsI(ii+1,kk+1)+(1.-t)*u*TsI(ii,kk+1))*normstress(k)+dyn_xz(i,k)
#if defined DIPSLIP
        striniZ(i,k)=(1.-t)*(1.-u)*T0I(ii,kk)+t*(1.-u)*T0I(ii+1,kk)+t*u*T0I(ii+1,kk+1)+(1.-t)*u*T0I(ii,kk+1)+dyn_xz(i,k)
        striniX(i,k)=0.
#else
        striniX(i,k)=(1.-t)*(1.-u)*T0I(ii,kk)+t*(1.-u)*T0I(ii+1,kk)+t*u*T0I(ii+1,kk+1)+(1.-t)*u*T0I(ii,kk+1)+dyn_xz(i,k)
        striniZ(i,k)=0.
#endif
      enddo
    enddo
    coh=0.5e6

    END SUBROUTINE
    
#endif
