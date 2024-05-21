!-------------------------------------------------------
! Main routine for the FD dynamic simulation and enforcement of
! the fault boundary condition.
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
! Preprocessor macros: OpenACC directives
#define _ACC_PARALLEL        !$acc parallel default (present)
#define _ACC_LOOP_COLLAPSE_2 !$acc loop collapse (2)
#define _ACC_LOOP_COLLAPSE_3 !$acc loop collapse (3)
#define _ACC_LOOP            !$acc loop
#define _ACC_END_PARALLEL    !$acc end parallel
! ------------------------------------------------------
      MODULE outputs_com


      real,allocatable,dimension (:,:):: sliprateoutX,sliprateoutZ
#if defined FVW
      real,allocatable,dimension (:,:):: psiout
#endif
      END MODULE



      SUBROUTINE fd3d()
      USE medium_com
      USE displt_com
      USE strfld_com
      USE fd3dparam_com
      USE friction_com
      USE source_com
      USE traction_com
      USE pml_com
	  USE outputs_com
      USE SlipRates_com
      use PostSeismic_com
	  use ieee_arithmetic
      IMPLICIT NONE

      real    :: time,friction,tmax,xmax,ymax,numer,denom,veltest,dd,slipmax
      real    :: pdx, pdz,tabs
      real    :: u1out
      real    :: CPUT1,CPUT2
      real    :: maxvelX,maxvelZ,maxvelsave,stoptimecheck
      real    :: tint, tint2, eraddt !temporary variables, time integrated part of radiated energy
      real    :: dht, ek, es, ef, c1, c2
      integer :: i,j,it,k, nxe, nxb, nyb, nye, nzb, nze
      integer :: ifrom,ito,jfrom,jto,kk,ii,jj
      integer :: stopnexttime
      integer, allocatable :: nkk(:)
      real    :: rup_tresh, rv, cz, efracds, alphakoef, schangef
      real,allocatable,dimension (:,:):: distX,distZ
      real,allocatable,dimension (:):: erad, efrac, timek
      real,allocatable,dimension (:,:,:)::slipt
!     real,allocatable,dimension (:,:,:):: waveU,waveV,waveW
#if defined FVW
      real    :: fss, flv, psiss, dpsi,  sr
      real    :: FXZ, GT, hx, hz, rr,AA,BB
#endif
!-------------------------------------------------------
!   initialize arrays
!-------------------------------------------------------
      
      allocate(v1(nxt,nyt,nzt))
      allocate(u1(nxt,nyt,nzt))
      allocate(w1(nxt,nyt,nzt))
      allocate(xx(nxt,nyt,nzt),yy(nxt,nyt,nzt),zz(nxt,nyt,nzt),xy(nxt,nyt,nzt),yz(nxt,nyt,nzt),xz(nxt,nyt,nzt))
      allocate(tx(nxt,nzt),tz(nxt,nzt),v1t(nxt,nzt),avdx(nxt,nzt),avdz(nxt,nzt),RFx(nxt,nzt),RFz(nxt,nzt))
      allocate(sliprateoutX(nxt,nzt),sliprateoutZ(nxt,nzt),distX(nxt,nzt),distZ(nxt,nzt))
      allocate(omega_pml(nabc-1), omegaR_pml(nabc-1),omega_pmlM(nabc-1), omegaR_pmlM(nabc-1))
      allocate(au1(nxt,nzt),av1(nxt,nzt),aw1(nxt,nzt))
      allocate(erad(nSR), efrac(nSR), slipt(nxt,nzt,nSR),timek(nSR))
!   allocate(waveU(nxt,nyt,nzt),waveV(nxt,nyt,nzt),waveW(nxt,nyt,nzt))
#if defined FVW
	  allocate(psiout(nxt,nzt))
#endif
      u1=0.; v1=0.; w1=0.
      xx=0.; yy=0.; zz=0.; xy=0.; yz=0.; xz=0.
      ruptime=1.e4; rise=0.; sliptime=1.e4
      tx=0.; tz=0.; v1t=0.
      uZ=0.; wX=0.
      avdx = 0.; avdz = 0.
      RFx = 0.; RFz = 0.
      au1=0.; av1=0.; aw1=0 
      MSRX=0.; MSRZ=0.; MomentRate=0.
      allocate(nkk(nSR))
      nkk=0
      stopnexttime=0
      distX=0.; distZ=0.
      rup_tresh=1e-3 !	Rate treshold for rupture time calculation
      c1  = 9./8. !	4th order central FD formula parameters	
      c2  = -1./24.
      tabsX=0.; tabsZ=0.
      sliprateoutZ=0.; sliprateoutX=0.
      SCHANGEZ=0.; SCHANGEX=0.
	  timek=0.
      tabs=0.; 
      slipX=0.; slipZ=0.
      efrac=0.;erad=0.;eraddt=0.;efracds=0.;schangef=0.;slipt=0.
!     waveU=0.; waveV=0.; waveW=0.
      dht = dh/dt
      if(Nstations>0) then
        seisU=0.;seisV=0.;seisW=0.
        seissurfU=0.;seissurfV=0.;seissurfW=0.
      endif
      call init_pml()
      call interp_fric()

#if defined FVW
      tabsX=1.
      tabsZ=1.
	  T0XI=T0X
	  T0ZI=T0Z
	  RFx=1.
	  RFz=1.
      psiout=0.
	  
!   Nucleation for fast velocity  weakening friction
!#if defined DIPSLIP
!    do k = nabc-1,nzt-2
!      do i = nabc-1,nxt-nabc
!        hx = real(i)*dh
!        hz = real(k)*dh
!        rr = (hx-hx0)**2 + (hz-hz0)**2
!        if (rr<RR2) then
!          FXZ = 1. !Alternatively smooth nucl. FXZ=exp(rr/(rr-RR2))
!        else
!          FXZ = 0.
!        endif
!          T0Z(i,k) = T0Z(i,k)+perturb*FXZ!SnZ(i,k)*(f0Z(i,k) + perturb*FXZ) !
!      enddo
!    enddo
!#else

!Fromer version of nucleation used in the Napa inversion!!!!
  !  do k = nabc-1,nzt-2	
  !    do i = nabc-1,nxt-nabc
  !      hx = real(i)*dh
  !      hz = real(k)*dh
  !      rr = (hx-hx0)**2 + (hz-hz0)**2
  !      if (rr<RR2) then
  !        FXZ=1.! !Alternatively smooth nucl. FXZ=exp(rr/(rr-RR2))
  !      else
  !        FXZ = 0.
  !      endif
  !        T0X(i,k) = T0X(i,k)+ SnX(i,k)*perturb*FXZ  !SnX(i,k)*(f0X(i,k) + perturb*FXZ)
  !    enddo
  !  enddo
!#endif
	
#endif



!---------------------------
! Write down the input
!---------------------------
      if (ioutput.eq.1) then
#if defined FVW
        open(95,file='result/vmodel.inp')
        open(96,file='result/friction.inp')
        do k = nabc+1,nzt-nfs
          do i = nabc+1,nxt-nabc
            write(95,*) mu1(i,nysc,k)
            write(96,'(100E13.5)') real(T0X(i,k),4),&
			real(T0Z(i,k),4), &
			real(aZ(i,k),4),&
			real(baZ(i,k),4),&
			real(vwZ(i,k),4),&
			real(SnZ(i,k),4),&
			real(DCZ(i,k),4),&
			real(log10(uini(i,k)),4),&
			real(psiZ(i,k),4),&
			real(f0Z(i,k),4),&
			real(DCZ(i,k)*log(v0/uini(i,k)),4),&
			real(DCZ(i,k)*exp((psiZ(i,k)-f0Z(i,k))/(baZ(i,k)+aZ(i,k)))/v0,4),&
			real(uini(i,k)*exp((psiZ(i,k)-f0Z(i,k))/(baZ(i,k)+aZ(i,k)))/v0,4)
          enddo
        enddo
        close(95)
        close(96)
#else
        open(95,file='result/vmodel.inp')
        open(96,file='result/friction.inp')
        do k = nabc+1,nzt-nfs
          do i = nabc+1,nxt-nabc
            write(95,*) mu1(i,nysc,k)
            write(96,'(5E13.5)') striniX(i,k),striniZ(i,k),peak_xz(i,k),Dc(i,k),peak_xz(i,k)/normstress(k)
          enddo
        enddo
        close(95)
        close(96)
#endif
        if (nstations>0) then
          open(31,file='result/seisout.dat')
          open(32,file='result/seisoutU.surface.gnuplot.dat')
          open(33,file='result/seisoutV.surface.gnuplot.dat')
          open(34,file='result/seisoutW.surface.gnuplot.dat')
        endif 
      endif
!-------------------------------------------------------
!     Loop over time
!-------------------------------------------------------
      if(ioutput.eq.1) then
#if defined FVW
        OPEN(24, file='result/psi.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
#endif
        OPEN(25, file='result/sliprateZ.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
        OPEN(27, file='result/sliprateX.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
        OPEN(26, file='result/shearstressZ.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
        OPEN(28, file='result/shearstressX.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
        OPEN(41, file='result/waveformX.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
        OPEN(42, file='result/waveformY.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
        OPEN(43, file='result/waveformZ.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')	
      endif
      
      CALL CPU_TIME(CPUT1)
      maxvelsave=0.
      stoptimecheck=min(1.,dt*real(ntfd)/10.)

      !$ACC DATA COPYIN (LAM1,MU1,D1) &
      !$ACC      COPYIN (U1,V1,W1) COPYIN (XX,YY,ZZ,XY,YZ,XZ) &
      !$ACC      COPYIN (tx,tz,v1t,avdx,avdz,RFx,RFz,uZ,wX) &
      !$ACC      COPYIN (au1,av1,aw1) &
      !$ACC      COPYIN (omega_pml,omegaR_pml,omega_pmlM,omegaR_pmlM) &
      !$ACC      COPYIN (u11,u12,u13,u21,u22,u23,u31,u32,u33,u41,u42,u43) &
      !$ACC      COPYIN (v11,v12,v13,v21,v22,v23,v31,v32,v33,v41,v42,v43) &
      !$ACC      COPYIN (w11,w12,w13,w21,w22,w23,w31,w32,w33,w41,w42,w43) &
      !$ACC      COPYIN (xx11,xx12,xx13,xx21,xx22,xx23,xx31,xx32,xx33,xx41,xx42,xx43) &
      !$ACC      COPYIN (yy11,yy12,yy13,yy21,yy22,yy23,yy31,yy32,yy33,yy41,yy42,yy43) &
      !$ACC      COPYIN (zz11,zz12,zz13,zz21,zz22,zz23,zz31,zz32,zz33,zz41,zz42,zz43) &
      !$ACC      COPYIN (xy11,xy12,xy21,xy22,xy31,xy32,xy41,xy42) &
      !$ACC      COPYIN (xz11,xz12,xz21,xz22,xz31,xz32,xz41,xz42) &
      !$ACC      COPYIN (yz11,yz12,yz21,yz22,yz31,yz32,yz41,yz42) &
      !$ACC      COPYIN (omegax1,omegax2,omegax3,omegax4) &
      !$ACC      COPYIN (omegay1,omegay2,omegay3,omegay4) &
      !$ACC      COPYIN (omegaz1,omegaz2,omegaz3,omegaz4) &
      !$ACC      COPYIN (omegaxS1,omegaxS2,omegaxS3,omegaxS4) &
      !$ACC      COPYIN (omegayS3,omegayS4,omegazS4) &
      !$ACC      COPYIN (dyn_xz,striniZ,striniX,peak_xz,Dc,coh,tabsX,tabsZ) &
      !$ACC      COPYIN (peakX, dynX, DcX, peakZ, dynZ, DcZ,staX,staY,staZ, distX, distZ) &
#if defined FVW
      !$ACC      COPYIN (aX,baX,psiX,vwX,SnX,f0X,fwX)&
      !$ACC      COPYIN (aZ,baZ,psiZ,vwZ,SnZ,f0Z,fwZ)&
      !$ACC      COPYIN (uini,wini,t0xi,t0zi)&
#endif
#if defined FSPACE
      !$ACC      COPYIN (u51,u52,u53,v51,v52,v53,w51,w52,w53)&
      !$ACC      COPYIN (xx51,xx52,xx53,yy51,yy52,yy53,zz51,zz52,zz53)&
      !$ACC      COPYIN (xy51,xy52,xz51,xz52,yz51,yz52)&
      !$ACC      COPYIN (omegaxS5,omegayS5,omegazS5) &
      !$ACC      COPYIN (omegax5,omegay5,omegaz5) &
#endif
      !$ACC      COPY (ruptime,sliptime,rise)
	  
      do it = 1,ntfd
	  
        time = dt*it
        sliprateoutX(1:nxt,1:nzt)=0.
        sliprateoutZ(1:nxt,1:nzt)=0.
        schangeZ(1:nxt,1:nzt)=0.
        schangeX(1:nxt,1:nzt)=0.

#if defined FVW
        psiout(1:nxt,1:nzt) = 0.
        !$ACC DATA COPY (T0X,T0Z,sliprateoutX,sliprateoutZ,schangeZ,schangeX,seisU,seisV,seisW,seissurfU,seissurfV,seissurfW,psiout,slipX,slipZ)
#else
        !$ACC DATA COPY (T0X,T0Z,sliprateoutX,sliprateoutZ,schangeZ,schangeX,seisU,seisV,seisW,seissurfU,seissurfV,seissurfW,slipX,slipZ)
#endif
!-------------------------------------------------------------
!   Velocity tick
!-------------------------------------------------------------
        call dvel(nxt,nyt-2,nzt,dt,dh)   !Apply 4th-order differencing to interior particle velocities
        call bnd2d(nxt,nyt,nzt,dt,dh)    !Compute velocities of 2nd order accuracy near fault and apply symmetry in normal component
        if (it>1) then
		call tasu1(nxt,nysc,nzt,dt,dh)   !Compute velocities at fault boundary
        call tasw1(nxt,nysc,nzt,dt,dh)
		endif
        call pml_uvw (nxt,nyt,nzt,dt,dh) !Absorbing condition
#if !defined FSPACE
        call fuvw(nxt,nyt,nzt-nfs)       !Stress free condition
#endif		
         
!-------------------------------------------------------------
!   Stress tick
!-------------------------------------------------------------
        call dstres(nxt,nyt-2,nzt,dh,dt) !4th-order differencing of stress
        call strbnd(nxt,nyt,nzt,dh,dt)   !Compute stress tensor of 2nd order accuracy near fault and apply symmetry in normal components
        call tasxz(nxt,nysc,nzt,dt,dh)   !Compute stress components at fault boundary
        call tasii(nxt,nysc,nzt,dt,dh)
        call pml_xyz (nxt,nyt,nzt,dt,dh) !Absorbing condition 
#if !defined FSPACE
        call fres(nxt,nyt,nzt-nfs)       !Stress free condition
#endif	
!----------------------------------------------------------------
! Fault BOUNDARY CONDITION applied in ABSOLUTE STRESS,
! not in relative stress, otherwise we get all sort of problems
!----------------------------------------------------------------
        !traction calculation
        _ACC_PARALLEL
        _ACC_LOOP_COLLAPSE_2
        do k = nabc+1,nzt-nfs-1
          do i = nabc+1,nxt-nabc
            !pdz = ((xz(i+1,nysc,k) - xz(i,nysc,k))/(2.*2.) + (zz(i,nysc,k+1) - zz(i,nysc,k))/(2.*2.) - yz(i,nysc-1,k))
            pdz = (c1*(xz(i+1,nysc,k) - xz(i,nysc,k))   +  c2*(xz(i+2,nysc,k) - xz(i-1,nysc,k)))/(2.) + &
              (c1*(zz(i,nysc,k+1) - zz(i,nysc,k))   +  c2*(zz(i,nysc,k+2) - zz(i,nysc,k-1)))/(2.) - &
              yz(i,nysc-1,k)
            avdz(i,k)= damp_s*(pdz - avdz(i,k))
            RFz(i,k) = pdz + avdz(i,k)
            tz(i,k) = -RFz(i,k) - 0.5*d1(i,nysc,k)*dht*w1(i,nysc,k)
            avdz(i,k) = pdz

            !pdx = ((xx(i,nysc,k) - xx(i-1,nysc,k))/(2.*2.) + (xz(i,nysc,k) - xz(i,nysc,k-1))/(2.*2.) - xy(i,nysc-1,k))
            pdx = (c1*(xx(i,nysc,k)   - xx(i-1,nysc,k)) +  c2*(xx(i+1,nysc,k) - xx(i-2,nysc,k)))/(2.) + &
              (c1*(xz(i,nysc,k)   - xz(i,nysc,k-1)) +  c2*(xz(i,nysc,k+1) - xz(i,nysc,k-2)))/(2.) - &
              xy(i,nysc-1,k)
            avdx(i,k)= damp_s*(pdx - avdx(i,k))
            RFx(i,k) = pdx + avdx(i,k)
            tx(i,k) = -RFx(i,k) - 0.5*d1(i,nysc,k)*dht*u1(i,nysc,k)
            avdx(i,k) = pdx
          enddo
        enddo
        _ACC_END_PARALLEL
        
        !Traction calculation near free surface
        k=nzt-nfs
#if defined FSPACE
        _ACC_PARALLEL
        _ACC_LOOP
      do i = nabc+1,nxt-nabc
		  pdz = (c1*(xz(i+1,nysc,k) - xz(i,nysc,k))   +  c2*(xz(i+2,nysc,k) - xz(i-1,nysc,k)))/(2.) + &
		    (c1*(zz(i,nysc,k+1) - zz(i,nysc,k))   +  c2*(zz(i,nysc,k+2) - zz(i,nysc,k-1)))/(2.) - &
		    yz(i,nysc-1,k)
		  avdz(i,k)= damp_s*(pdz - avdz(i,k))
		  RFz(i,k) = pdz + avdz(i,k)
		  tz(i,k) = -RFz(i,k) - 0.5*d1(i,nysc,k)*dht*w1(i,nysc,k)
		  avdz(i,k) = pdz
		  
          pdx = (c1*(xx(i,nysc,k)   - xx(i-1,nysc,k)) +  c2*(xx(i+1,nysc,k) - xx(i-2,nysc,k)))/(2.) + &
		    (c1*(xz(i,nysc,k)   - xz(i,nysc,k-1)) +  c2*(xz(i,nysc,k+1) - xz(i,nysc,k-2)))/(2.) - &
		    xy(i,nysc-1,k)
		  avdx(i,k)= damp_s*(pdx - avdx(i,k))
		  RFx(i,k) = pdx + avdx(i,k)
		  tx(i,k) = -RFx(i,k) - 0.5*d1(i,nysc,k)*dht*u1(i,nysc,k)
		  avdx(i,k) = pdx
	    enddo
		_ACC_END_PARALLEL
#else
        _ACC_PARALLEL
        _ACC_LOOP
        do i = nabc+1,nxt-nabc
          pdz = ((xz(i+1,nysc,k-1) - xz(i,nysc,k-1))/(2.*4) + (zz(i,nysc,k+1) - zz(i,nysc,k))/(2.*2.) - yz(i,nysc-1,k-1)/(4.))
          avdz(i,k)= damp_s*(pdz - avdz(i,k))
          RFz(i,k) = pdz + avdz(i,k)
          tz(i,k) = -RFz(i,k) - 0.5*d1(i,nysc,k)*dht*w1(i,nysc,k)/2.
          avdz(i,k) = pdz
          
          !pdx = ((xx(i,nysc,k) - xx(i-1,nysc,k))/(2.*2.) + (xz(i,nysc,k) - xz(i,nysc,k-1))/(2.*2.) - xy(i,nysc-1,k))
          pdx = (c1*(xx(i,nysc,k)   - xx(i-1,nysc,k)) +  c2*(xx(i+1,nysc,k) - xx(i-2,nysc,k)))/(2.) + &
            (c1*(xz(i,nysc,k)   - xz(i,nysc,k-1)) +  c2*(xz(i,nysc,k+1) - xz(i,nysc,k-2)))/(2.) - &
            xy(i,nysc-1,k)
          avdx(i,k)= damp_s*(pdx - avdx(i,k))
          RFx(i,k) = pdx + avdx(i,k)
          tx(i,k) = -RFx(i,k) - 0.5*d1(i,nysc,k)*dht*u1(i,nysc,k)
          avdx(i,k) = pdx
        enddo 
        _ACC_END_PARALLEL
        
        _ACC_PARALLEL
        _ACC_LOOP
        do i = nabc+1,nxt-nabc
          tx(i,nzt-1)=+tx(i,nzt-2)
          tz(i,nzt-1)=-tz(i,nzt-3)
		  RFx(i,nzt-1)=RFx(i,nzt-2)
		!  RFz(i,nzt-1)=-RFz(i,nzt-3)
		
        enddo
        _ACC_END_PARALLEL
#endif
!-------------------------------------------------------------
!  Apply Boundary Conditions
!-------------------------------------------------------------
        !Traction interpolation in staggered positions
        _ACC_PARALLEL
        _ACC_LOOP_COLLAPSE_2
        do k = nabc+1,nzt-nfs-1
          do i = nabc+1,nxt-nabc-1
            tint=(tx(i,k)+T0X(i,k)+tx(i+1,k)+T0X(i+1,k)+tx(i,k+1)+T0X(i,k+1)+tx(i+1,k+1)+T0X(i+1,k+1))/4.
            tabsZ(i,k) = sqrt(tint**2 + (tz(i,k)+T0Z(i,k))**2)
            tint=(tz(i,k)+T0Z(i,k)+tz(i-1,k)+T0Z(i-1,k)+tz(i,k-1)+T0Z(i,k-1)+tz(i-1,k-1)+T0Z(i-1,k-1))/4.
            tabsX(i,k) = sqrt((tx(i,k)+T0X(i,k))**2 + (tint)**2)
	        uZ(i,k)=(U1(I,NYSC,K)+U1(I+1,NYSC,K)+U1(I,NYSC,K+1)+U1(I+1,NYSC,K+1))/4.
	        wX(i,k)=(w1(I,NYSC,K)+w1(I-1,NYSC,K)+w1(I,NYSC,K-1)+w1(I-1,NYSC,K-1))/4.
          enddo
        enddo
        _ACC_END_PARALLEL
		
		i=nxt-nabc 		
		_ACC_PARALLEL
        _ACC_LOOP
        do k = nabc+1,nzt-nfs-1
            tint=(tx(i,k)+T0X(i,k)+tx(i,k+1)+T0X(i,k+1))/2.
            tabsZ(i,k) = sqrt(tint**2 + (tz(i,k)+T0Z(i,k))**2)
            tint=(tz(i,k)+T0Z(i,k)+tz(i-1,k)+T0Z(i-1,k)+tz(i,k-1)+T0Z(i,k-1)+tz(i-1,k-1)+T0Z(i-1,k-1))/4.
            tabsX(i,k) = sqrt((tx(i,k)+T0X(i,k))**2 + (tint)**2)
	        uZ(i,k)=(U1(I,NYSC,K)+U1(I,NYSC,K+1))/2.
	        wX(i,k)=(w1(I,NYSC,K)+w1(I-1,NYSC,K)+w1(I,NYSC,K-1)+w1(I-1,NYSC,K-1))/4.
        enddo
        _ACC_END_PARALLEL
		
        k=nzt-nfs

	    _ACC_PARALLEL
        _ACC_LOOP
          do i = nabc+1,nxt-nabc+1
            tint=(tx(i,k)+T0X(i,k)+tx(i+1,k)+T0X(i+1,k))/2.
            tabsZ(i,k) = sqrt(tint**2 + (tz(i,k)+T0Z(i,k))**2)
            tint=(tz(i,k)+T0Z(i,k)+tz(i-1,k)+T0Z(i-1,k)+tz(i,k-1)+T0Z(i,k-1)+tz(i-1,k-1)+T0Z(i-1,k-1))/4.
            tabsX(i,k) = sqrt((tx(i,k)+T0X(i,k))**2 + (tint)**2)
	        uZ(i,k)=(U1(I,NYSC,K)+U1(I+1,NYSC,K))/2.
	        wX(i,k)=(w1(I,NYSC,K)+w1(I-1,NYSC,K)+w1(I,NYSC,K-1)+w1(I-1,NYSC,K-1))/4.
          enddo
        _ACC_END_PARALLEL
		
        i=nxt-nabc
	    k=nzt-nfs
        tint=(tx(i,k)+T0X(i,k))
        tabsZ(i,k) = sqrt(tint**2 + (tz(i,k)+T0Z(i,k))**2)
        tint=(tz(i,k)+T0Z(i,k)+tz(i-1,k)+T0Z(i-1,k)+tz(i,k-1)+T0Z(i,k-1)+tz(i-1,k-1)+T0Z(i-1,k-1))/4.
        tabsX(i,k) = sqrt((tx(i,k)+T0X(i,k))**2 + (tint)**2)
	    uZ(i,k)=U1(I,NYSC,K)
	    wX(i,k)=(w1(I,NYSC,K)+w1(I-1,NYSC,K)+w1(I,NYSC,K-1)+w1(I-1,NYSC,K-1))/4.
		  
        _ACC_PARALLEL
        _ACC_LOOP_COLLAPSE_2
        do k = nabc+1,nzt-nfs
          do i = nabc+1,nxt-nabc
            tabs=tabsZ(i,k)
            u1out=-sqrt(W1(I,NYSC,K)**2+uZ(i,k)**2)

#if defined FVW
            
            sr=sqrt((2*w1(i,nyt,k)-wini(i,k))**2+(2*uZ(i,k)-uini(i,k))**2)                          
            flv = f0Z(i,k) - baZ(i,k)*log(sr/v0)
            fss = fwZ(i,k) + (flv - fwZ(i,k))/((1. + (sr/vwZ(i,k))**8)**(1./8.))
            psiss = aZ(i,k)*(log(sinh(fss/aZ(i,k))) + log(2*v0/(sr))) 

            psiZ(i,k)=(psiZ(i,k)-psiss)*exp(-sr*dt/DcZ(i,k)) + psiss
            friction  = SnZ(i,k) * aZ(i,k)*asinh(sr*exp(psiZ(i,k)/aZ(i,k))/(2*v0)) +coh(i,k)
            distZ(i,k) = distZ(i,k)  - 2*u1out*dt
            schangeZ(I,K) = (tz(i,k) + t0Z(i,k))*friction/tabs! - t0Z(i,k)
            psiout(i,k)=psiZ(i,k)
			if (abs(2*u1out)>rup_tresh) then
                if (ruptime(i,k).eq.1.e4) ruptime(i,k) = time				
            endif


#else
            if (distZ(i,k).le.DcZ(i,k)) then
              friction = peakZ(i,k) * (1.0 - distZ(i,k)/DcZ(i,k)) + dynZ(i,k)*distZ(i,k)/DcZ(i,k) + coh(i,k)
            else
              friction = dynZ(i,k) + coh(i,k)
            endif
            
            if (tabs>=friction.and.time<=SRdur) then
              distZ(i,k) = distZ(i,k)  - 2*u1out*dt
              tz(i,k) =  (tz(i,k) + T0Z(i,k))*friction/tabs - T0Z(i,k)

              if (-2*u1out>rup_tresh) then
                if (ruptime(i,k).ne.1.e4) rise(i,k) = time
                if (ruptime(i,k).eq.1.e4) ruptime(i,k) = time
              endif
            endif
#endif

            if ((sliptime(i,k)==1.e4).AND.(distZ(i,k)>Dc(i,k))) sliptime(i,k)=time
          enddo
        enddo
        _ACC_END_PARALLEL
        
        _ACC_PARALLEL
        _ACC_LOOP_COLLAPSE_2
        do k = nabc+1,nzt-nfs
          do i = nabc+1,nxt-nabc
            tabs=tabsX(i,k)
            u1out=-sqrt(wX(i,k)**2+U1(I,NYSC,K)**2)
            
#if defined FVW
            sr=sqrt((2*wX(i,k)-wini(i,k))**2+(2*u1(i,nyt,k)-uini(i,k))**2)
            flv = f0X(i,k) - baX(i,k)*log(sr/v0)
            fss = fwX(i,k) + (flv - fwX(i,k))/((1. + (sr/vwX(i,k))**8)**(1./8.))
            psiss = aX(i,k)*(log(sinh(fss/aX(i,k))) + log(2*v0/(sr))) 
            psiX(i,k)=(psiX(i,k)-psiss)*exp(-sr*dt/DcX(i,k)) + psiss
            friction  = SnX(i,k) * aX(i,k)*asinh(sr*exp(psiX(i,k)/aX(i,k))/(2*v0)) 	+ coh(i,k)		
            distX(i,k) = distX(i,k)  - 2*u1out*dt
            schangeX(I,K) = (tx(i,k) + t0X(i,k))*friction/tabs! - t0X(i,k)
            
#else
            if (distX(i,k).le.Dc(i,k)) then
              friction = peakX(i,k) * (1.0 - distX(i,k)/DcX(i,k)) + dynX(i,k)*distX(i,k)/DcX(i,k) + coh(i,k)
            else
              friction = dynX(i,k) + coh(i,k)
            endif
            
            if (tabs>=friction.and.time<=SRdur) then
              distX(i,k) = distX(i,k)  - 2*u1out*dt
              tx(i,k) = (tx(i,k) + T0X(i,k))*friction/tabs - T0X(i,k)
            endif
#endif 

          enddo
        enddo
        _ACC_END_PARALLEL

        _ACC_PARALLEL
        _ACC_LOOP_COLLAPSE_2
        do k = nabc+1,nzt-nfs-1!+1
          do i = nabc+1,nxt-nabc-1!+1
#if defined FVW			
            sliprateoutX(i,k) = (-2.*U1(I,NYSC,K)-2.*U1(I+1,NYSC,K)-2.*U1(I,NYSC,K+1)-2.*U1(I+1,NYSC,K+1))/4.
	        sliprateoutZ(i,k) = - 2.*W1(I,NYSC,K)
            slipZ(i,k)=slipZ(i,k)+sliprateoutZ(i,k)*dt
            slipX(i,k)=slipX(i,k)+sliprateoutX(i,k)*dt
#else
            SCHANGEZ(I,K) = tz(i,k) + T0Z(i,k)
            sliprateoutZ(i,k) = - 2.*W1(I,NYSC,K)
            SCHANGEX(I,K) = (tx(i,k)+T0X(i,k)+tx(i+1,k)+T0X(i+1,k)+tx(i,k+1)+T0X(i,k+1)+tx(i+1,k+1)+T0X(i+1,k+1))/4.
            sliprateoutX(i,k) = (-2.*U1(I,NYSC,K)-2.*U1(I+1,NYSC,K)-2.*U1(I,NYSC,K+1)-2.*U1(I+1,NYSC,K+1))/4.
            slipZ(i,k)=slipZ(i,k)+sliprateoutZ(i,k)*dt
            slipX(i,k)=slipX(i,k)+sliprateoutX(i,k)*dt
#endif 
            efracds=efracds+schangeX(i,k)*sliprateoutX(i,k)*dt+schangeZ(i,k)*sliprateoutZ(i,k)*dt
            eraddt=eraddt+(schangeX(i,k)-t0X(i,k))*sliprateoutX(i,k)*dt+(schangeZ(i,k)-t0Z(i,k))*sliprateoutZ(i,k)*dt
		  enddo
        enddo
        _ACC_END_PARALLEL
		
        i=nxt-nabc
        _ACC_PARALLEL
        _ACC_LOOP
        do k = nabc+1,nzt-nfs-1!+1
#if defined FVW			
            sliprateoutX(i,k) = (-2.*U1(I,NYSC,K)-2.*U1(I,NYSC,K+1))/2.
	        sliprateoutZ(i,k) = - 2.*W1(I,NYSC,K)
            slipZ(i,k)=slipZ(i,k)+sliprateoutZ(i,k)*dt
            slipX(i,k)=slipX(i,k)+sliprateoutX(i,k)*dt
#else
            SCHANGEZ(I,K) = tz(i,k) + T0Z(i,k)
            sliprateoutZ(i,k) = - 2.*W1(I,NYSC,K)
            SCHANGEX(I,K) = (tx(i,k)+T0X(i,k)+tx(i,k+1)+T0X(i,k+1))/2.
            sliprateoutX(i,k) = (-2.*U1(I,NYSC,K)-2.*U1(I,NYSC,K+1))/2.
            slipZ(i,k)=slipZ(i,k)+sliprateoutZ(i,k)*dt
            slipX(i,k)=slipX(i,k)+sliprateoutX(i,k)*dt
#endif 
            efracds=efracds+schangeX(i,k)*sliprateoutX(i,k)*dt+schangeZ(i,k)*sliprateoutZ(i,k)*dt
            eraddt=eraddt+(schangeX(i,k)-t0X(i,k))*sliprateoutX(i,k)*dt+(schangeZ(i,k)-t0Z(i,k))*sliprateoutZ(i,k)*dt
        enddo
        _ACC_END_PARALLEL

        k=nzt-nfs!+1

	_ACC_PARALLEL
        _ACC_LOOP
          do i = nabc+1,nxt-nabc+1
            
#if defined FVW
            sliprateoutX(i,k) = (-2.*U1(I,NYSC,K)-2.*U1(I+1,NYSC,K))/2.
	        sliprateoutZ(i,k) = - 2.*W1(I,NYSC,K)
#else
            SCHANGEZ(I,K) = tz(i,k) + T0Z(i,k)
            sliprateoutZ(i,k) = - 2.*W1(I,NYSC,K)
            SCHANGEX(I,K) = (tx(i,k)+T0X(i,k)+tx(i+1,k)+T0X(i+1,k))/2.
            sliprateoutX(i,k) = (-2.*U1(I,NYSC,K)-2.*U1(I+1,NYSC,K))/2.
#endif 

          enddo
        _ACC_END_PARALLEL
		
        if(ioutput.eq.1) then
        if (Nstations>0) then
          _ACC_PARALLEL
          _ACC_LOOP
          do i=1,Nstations
            seisU(i)=seisU(i)+u1(staX(i),staY(i),staZ(i))
            seisV(i)=seisV(i)+v1(staX(i),staY(i),staZ(i))
            seisW(i)=seisW(i)+w1(staX(i),staY(i),staZ(i))
          enddo
          _ACC_END_PARALLEL
          _ACC_PARALLEL
          _ACC_LOOP_COLLAPSE_2
		  do j=1,nyt
            do i=1,nxt
              seissurfU(i,j)=seissurfU(i,j)+u1(i,j,nzt-nfs)
              seissurfV(i,j)=seissurfV(i,j)+v1(i,j,nzt-nfs)
              seissurfW(i,j)=seissurfW(i,j)+w1(i,j,nzt-nfs)
            enddo
          enddo
          _ACC_END_PARALLEL
        endif
        endif
		
        ! if(ioutput.eq.1) then
		! if ((time .GE. waveT) .AND. (waveT .NE. 0.)) then
!		! _ACC_PARALLEL
!        ! _ACC_LOOP_COLLAPSE_3
!		   do k=1,nzt
!		   do j=1,nyt
!          do i=1,nxt
!             waveU(i,j,k)=u1(i,j,k)
!             waveV(i,j,k)=v1(i,j,k)
!             waveW(i,j,k)=w1(i,j,k)
!           enddo
!		   enddo
!		   enddo
!        ! _ACC_END_PARALLEL
!         endif
!         endif  
    
#if defined TPV103
!   Smooth nucleation for the tpv104 benchmark
        if ((time-dt/2.) <= TT2) then
        
          if ((time-dt/2.) < TT2) then
            GT=exp(((time-dt/2.)-TT2)**2/((time-dt/2.)*((time-dt/2.)-2*TT2)))
          else
            GT=1.
          endif
          _ACC_PARALLEL
          _ACC_LOOP_COLLAPSE_2
          do k = nabc-1,nzt-nfs
            do i = nabc-1,nxt-nabc
              hx = real(i)*dh
              hz = real(k)*dh
              rr = (hx-hx0)**2 + (hz-hz0)**2
              if (rr<RR2-dh) then
                FXZ=exp(rr/(rr-RR2-dh))
              else
                FXZ = 0.
              endif
              T0X(i,k) = strinixI + perturb*FXZ!*GT

            enddo
          enddo
          _ACC_END_PARALLEL
        endif
#endif
#if defined TPV104
!   Smooth nucleation for the tpv104 benchmark
        if ((time-dt/2.) <= TT2) then
        
          if ((time-dt/2.) < TT2) then
            GT=exp(((time-dt/2.)-TT2)**2/((time-dt/2.)*((time-dt/2.)-2*TT2)))
          else
            GT=1.
          endif
          _ACC_PARALLEL
          _ACC_LOOP_COLLAPSE_2
          do k = nabc-1,nzt-nfs
            do i = nabc-1,nxt-nabc
              hx = real(i)*dh
              hz = real(k)*dh
              rr = (hx-hx0)**2 + (hz-hz0)**2
              if (rr<RR2-dh) then
                FXZ=exp(rr/(rr-RR2-dh))
              else
                FXZ = 0.
              endif

              T0X(i,k) = T0XI(i,k) + SnX(i,k)*perturb*FXZ*GT

            enddo
          enddo
          _ACC_END_PARALLEL
        endif
#endif

#if defined FVW
!Smooth nucleation for inversion

		if ((time) <= TT2) then
        
          if ((time) < TT2) then
            GT=exp(((time)-TT2)**2/((time)*((time)-2*TT2)))
          else
            GT=1.
          endif
_ACC_PARALLEL
_ACC_LOOP_COLLAPSE_2
          do k = nabc-1,nzt-2
            do i = nabc-1,nxt-nabc
              hx = real(i)*dh
              hz = real(k)*dh
              rr = (hx-hx0)**2 + (hz-hz0)**2
              if (rr<RR2) then
                FXZ=exp(rr/(rr-RR2))!FXZ=1.
              else
                FXZ = 0.
              endif
#if defined DIPSLIP
				T0Z(i,k) = T0ZI(i,k)+ SnZ(i,k)*perturb*FXZ*GT 
#else
                T0X(i,k) = T0XI(i,k)+ SnX(i,k)*perturb*FXZ*GT  
#endif
            enddo
          enddo
_ACC_END_PARALLEL
        endif

#endif

!$ACC END DATA
        if(mod(it,ceiling(stoptimecheck/dt))==0)then
          maxvelZ=maxval(sliprateoutZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs))
          maxvelX=maxval(sliprateoutX(nabc+1:nxt-nabc,nabc+1:nzt-nfs))
          write(*,*)'Time: ',time,'Slip rate max: ',maxvelX,maxvelZ
        endif
        if(ioutput.eq.1.and.Nstations>0)then
          if(mod(it,ceiling(dtseis/dt))==0)then
            write (31,'(1000000E13.5)')time,(seisU(i)/float(ceiling(dtseis/dt)),seisV(i)/float(ceiling(dtseis/dt)),seisW(i)/float(ceiling(dtseis/dt)),i=1,Nstations)
            seisU=0.;seisV=0.;seisW=0.
            do i=nabc+1,nxt-nabc,2
              write(32,'(1000000E13.5)')(seissurfU(i,j)/float(ceiling(dtseis/dt)),j=nabc+1,nyt,2)
              write(33,'(1000000E13.5)')(seissurfV(i,j)/float(ceiling(dtseis/dt)),j=nabc+1,nyt,2)
              write(34,'(1000000E13.5)')(seissurfW(i,j)/float(ceiling(dtseis/dt)),j=nabc+1,nyt,2)
            enddo
            write(32,*);write(32,*)
            write(33,*);write(33,*)
            write(34,*);write(34,*)
            seissurfU=0.;seissurfV=0.;seissurfW=0.;
          endif
        endif

        k=int(real(it-1)*dt/dtseis)+1
        nkk(k)=nkk(k)+1

        if(k<=nSR)then
		timek(k)=time
        if(ioutput.eq.1) then
#if defined FVW
          WRITE(24) psiout(nabc+1:nxt-nabc,nabc+1:nzt-nfs)
#endif
#if defined DIPSLIP
          WRITE(25) sliprateoutZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs)
          WRITE(26) SCHANGEZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs)
#else
          WRITE(27) sliprateoutX(nabc+1:nxt-nabc,nabc+1:nzt-nfs)
          WRITE(28) SCHANGEX(nabc+1:nxt-nabc,nabc+1:nzt-nfs)
#endif
        endif
 
        do j=1,NW
          jto=max(1,int(dW/dh*j))+1+nabc
          jfrom=min(jto,int(dW/dh*(j-1))+1)+1+nabc
          do i=1,NL
            ifrom=int(dL/dh*(i-1))+1+nabc
            ito=int(dL/dh*i)+nabc
            kk=((j-1)*NL+i-1)*nSR+k
            MSRX(kk)=MSRX(kk)+sum(sliprateoutX(ifrom:ito,jfrom:jto))/(dble((ito-ifrom+1)*(jto-jfrom+1)))!*(dtseis/dt))
            MSRZ(kk)=MSRZ(kk)+sum(sliprateoutZ(ifrom:ito,jfrom:jto))/(dble((ito-ifrom+1)*(jto-jfrom+1)))!*(dtseis/dt))
          enddo
        enddo
        
        do j = nabc+1,nzt-nfs
          do i = nabc+1,nxt-nabc
            MomentRate(k)=MomentRate(k)+sqrt(sliprateoutX(i,j)**2+sliprateoutZ(i,j)**2)*muSource(i,j)*dh*dh!/(dtseis/dt)
          enddo
        enddo
        !tint=0.
        tint2=0.
        do j = nabc+1,nzt-nfs
          do i = nabc+1,nxt-nabc
             !tint=tint+sqrt(schangeX(i,k)**2+schangeZ(i,k)**2)*sqrt(slipX(i,k)**2+slipZ(i,k)**2)
             tint2=tint2+0.5*((schangeX(i,j)-t0X(i,j))*slipX(i,j)+(schangeZ(i,j)-t0Z(i,j))*slipZ(i,j))
             slipt(i,j,k)=sqrt(slipX(i,j)**2+slipZ(i,j)**2)
          enddo
        enddo
        efrac(k)=efracds
        erad(k)=tint2-eraddt
		
        if(mod(it,ceiling(stoptimecheck/dt))==0)then
          if(stopnexttime==1)exit
#if defined DIPSLIP
          if(maxvelZ<1.e-7)stopnexttime=1
          if(maxvelZ>maxvelsave)maxvelsave=maxvelZ
          if(maxvelZ<=0.01*maxvelsave)stopnexttime=1
!          stopnexttime=0  !Don't stop
#else
          if(maxvelX<1.e-7)stopnexttime=1
          if(maxvelX>maxvelsave)maxvelsave=maxvelX
          if(maxvelX<=0.01*maxvelsave)stopnexttime=1
!          stopnexttime=0  !Don't stop
#endif
        endif
        
        endif
! output waveforms
!		if ((time .GE. waveT) .AND. (waveT .NE. 0.)) then
!			waveT=0.;
!			WRITE(41) waveU(1:nxt,1:nyt,1:nzt)
!			WRITE(42) waveV(1:nxt,1:nyt,1:nzt)
!			WRITE(43) waveW(1:nxt,1:nyt,1:nzt)
!		endif

      enddo ! --- End of the time loop
       
      !$ACC END DATA
      
      
      !divide by number of dt timesteps in each dtseis time interval
      do k=1,nSR
        if (nkk(k)>0) then
          MomentRate(k)=MomentRate(k)/nkk(k)
          do j=1,NW
            do i=1,NL
              kk=((j-1)*NL+i-1)*nSR+k
              MSRX(kk)=MSRX(kk)/nkk(k)
              MSRZ(kk)=MSRZ(kk)/nkk(k)
            enddo
          enddo
        endif           
      enddo
      
!===============================================================================================================================================
!Postprocessing
!===============================================================================================================================================

	  if (igps==1) then
#if defined FVW	  
	    call QDYN3D()
        !MSX=0.
		!MSZ=0.    
        !do k=1,nSr 
        !    kk=0
        !    do j=1,NW
        !        do i=1,NL
		!          kk=((j-1)*NL+i-1)*nSR!kk=kk+1
        !          MSX(kk)=MSX(kk)+MSRX(kk+(k-1)*NL*NW)*dtseis
	!			  MSZ(kk)=MSZ(kk)+MSRZ(kk+(k-1)*NL*NW)*dtseis
    !            enddo
    !        enddo
    !    enddo
		
#else
        MSX=0.
		MSZ=0.
        do k=1,NSr
            kk=0
            do j=1,NW
                do i=1,NL
		          kk=kk+1
                  MSX(kk)=MSX(kk)+MSRX(kk+(k-1)*NL*NW)*dtseis
				  MSZ(kk)=MSZ(kk)+MSRZ(kk+(k-1)*NL*NW)*dtseis
                enddo
            enddo
        enddo
#endif	 
      endif	 
	  
#if defined FVW	  
      if(ioutput.eq.1) then
	  	do k=1,NSr
        do j = nabc+1,nzt-nfs
          do i = nabc+1,nxt-nabc
			 if ((slipt(i,j,k) >= 0.95*maxval(slipt(i,j,:))).and.(rise(i,j)==0.)) then
				rise(i,j)=timek(k)-ruptime(i,j)
			 endif
          enddo
        enddo		
		enddo
	  endif
#else
	rise=rise-ruptime
#endif
      SCHANGEZ(:,:)=SCHANGEZ(:,:)-T0Z(:,:)   !stress drop
      SCHANGEX(:,:)=SCHANGEX(:,:)-T0X(:,:)
     ! schangef=sum(distZ*sqrt((schangeX+T0X)**2+(schangeZ+T0Z)**2))

      deallocate(u1,v1,w1)
      deallocate(xx,yy,zz,xy,yz,xz)
      deallocate(tx,tz,v1t,avdx,avdz, RFx,RFz)
      deallocate(sliprateoutX,sliprateoutZ,distX,distZ)
      deallocate(omega_pml,omegaR_pml,omega_pmlM,omegaR_pmlM,au1,av1,aw1)
      deallocate(omegax1,omegay1,omegaz1,omegaxS1)    
      deallocate(u11,u12,u13,v11,v12,v13,w11,w12,w13)
      deallocate(xx11,xx12,xx13,yy11,yy12,yy13,zz11,zz12,zz13,xz11,xz12,xy11,xy12,yz11,yz12)
      deallocate(omegax2,omegay2,omegaz2,omegaxS2)      
      deallocate(u21,u22,u23,v21,v22,v23,w21,w22,w23)
      deallocate(xx21,xx22,xx23,yy21,yy22,yy23,zz21,zz22,zz23,xz21,xz22,xy21,xy22,yz21,yz22)  
      deallocate(omegax3,omegay3,omegaz3,omegaxS3,omegayS3) 
      deallocate(u31,u32,u33,v31,v32,v33,w31,w32,w33)
      deallocate(xx31,xx32,xx33,yy31,yy32,yy33,zz31,zz32,zz33,xz31,xz32,xy31,xy32,yz31,yz32)  
      deallocate(omegax4,omegay4,omegaz4,omegaxS4,omegayS4,omegazS4)
      deallocate(u41,u42,u43,v41,v42,v43,w41,w42,w43)
      deallocate(xx41,xx42,xx43,yy41,yy42,yy43,zz41,zz42,zz43,xz41,xz42,xy41,xy42,yz41,yz42)    
#if defined FSPACE	  
      deallocate(omegax5,omegay5,omegaz5,omegaxS5,omegayS5,omegazS5)  
      deallocate(u51,u52,u53,v51,v52,v53,w51,w52,w53)
      deallocate(xx51,xx52,xx53,yy51,yy52,yy53,zz51,zz52,zz53,xz51,xz52,xy51,xy52,yz51,yz52)    
#endif	  
#if defined FVW
      deallocate (psiout,t0xi,t0zi)
#endif

      CALL CPU_TIME(CPUT2)
      PRINT *,'CPU TIME OF TIME LOOP: ',CPUT2-CPUT1

!-------------------
! Open output files:
!-------------------
      if(ioutput.eq.1) then
#if defined FVW
        close(24)
#endif
        close(25)
        close(26)
        close(27)
        close(28)
        close(31)
        close(32);close(33);close(34)
      endif
      
      tmax            = -1.
      output_param(1) =  0.
      output_param(2) =  0.
      numer           =  0.
      denom           =  0.
      schangef        =  0.
      do k = nabc+1,nzt-nfs	
        do i = nabc+1,nxt-nabc
          ! --- Seismic moment:
          output_param(2) = output_param(2) + sqrt(slipX(i,k)**2 + slipZ(i,k)**2)*muSource(i,k)*(dh*dh)
          ! --- Stress drop:
          
          numer = numer + schangeZ(i,k)*slipZ(i,k)+schangeX(i,k)*slipX(i,k)	!correct only for pure dipslip/strikeslip
          denom = denom + sqrt(slipX(i,k)**2 + slipZ(i,k)**2)
          if (denom .ne. 0.0) then
            output_param(4) = -(numer/denom)
          else
            output_param(4) = 0.0
          endif
          schangef=schangef+slipX(i,k)*(schangeX(i,k)+T0X(i,k))+slipZ(i,k)*(schangeZ(i,k)+T0Z(i,k))
        enddo
      enddo
    ! output_param(5) = (1./2.)**sum(peak_xz*Dc)/dble((nxt-2*nabc)*(nzt-nfs-nabc))
    ! output_param(6) = (1./2.)*output_param(4)*(output_param(2)/(mu_mean*output_param(3)))
      M0=output_param(2)
      Mw=(log10(M0)-9.1)/1.5

      do k=1,nSR
        !write(398,*)dt*(k-1),efrac(k),(efrac(k)-schangef)*dh**2
        if ((efrac(k)-schangef)*dh**2>0.) Eg=(efrac(k)-schangef)*dh**2
        if (erad(k)*dh**2>0.) Er=erad(k)*dh**2		   
      enddo
      output_param(5)=Eg
      output_param(6)=Er

!---------------------------
! Write down the output
!---------------------------

      if (ioutput.eq.1) then
        open(96,file='result/risetime.res')
        open(97,file='result/ruptime.res')
#if defined DIPSLIP
        open(94,file='result/slipZ.res')
        open(95,file='result/stressdropZ.res')
#else
        open(98,file='result/slipX.res')
        open(99,file='result/stressdropX.res')
#endif
        open(297,FILE='mtildeX.dat')
        open(298,FILE='mtildeZ.dat')
        open(299,FILE='mtildemomentrate.dat')

        do k=1,NL*NW*nSR
           write(297,'(1E13.5)')MSRX(k) 
           write(298,'(1E13.5)')MSRZ(k)	
        enddo
        do k=1,nSR
           write(299,*)dtseis*(k-1),Momentrate(k)
        enddo

		if(igps==1)then
          open(295,FILE='mtildeT.dat')
          open(296,FILE='mtildeS.dat')
          do k=1,NL*NW*NTS
#if defined DIPSLIP
            write (296,*) MSZ(k)
#else
            write (296,*) MSX(k)
#endif
          enddo
          do k=1,NTS
            write(295,*) TS(k)
          enddo
          close(295)
          close(296)
        endif

        close(297)
        close(298)
        close(299)

        do k = nabc+1,nzt-nfs
          do i = nabc+1,nxt-nabc
            write(96,*) rise(i,k)
            write(97,*) ruptime(i,k)
#if defined DIPSLIP
            write(94,*) slipZ(i,k)
            write(95,*) schangeZ(i,k)
#else
            write(98,*) slipX(i,k)
            write(99,*) schangeX(i,k)
#endif
          enddo
        enddo
        close(96)
        close(97)
        close(98)
        close(99)
        close(94)
        close(95)

        open(501,file='result/contour.res')
        write(501,*) 'j k t'
        do k = nabc+1,nzt-nfs
          do i = nabc+1,nxt-nabc
            write (501,*) (real(i-nabc-1) - real(nxt-2*nabc)/2.)*dh,real(nzt-k-nabc-1)*dh,ruptime(i,k) 
          enddo
        enddo
        close(501)
        open(502,file='result/czone.res')
        open(503,file='result/rvel.res')
        tint=0.
        tint2=0.
        do k = nabc+1,nzt-nfs
          do i = nabc+1,nxt-nabc
            cz=0.
            rv=0.
            if ((ruptime(i,k).ne.1.e4).and.(ruptime(i+1,k).ne.1.e4).and.(ruptime(i,k+1).ne.1.e4) &
            .and.(ruptime(i-1,k).ne.1.e4).and.(ruptime(i,k-1).ne.1.e4)) then
              rv = (sqrt((ruptime(i+1,k)-ruptime(i,k))**2+(ruptime(i,k+1)-ruptime(i,k))**2) &
                +sqrt((ruptime(i,k)-ruptime(i-1,k))**2+(ruptime(i,k)-ruptime(i,k-1))**2))
              if (rv.ne.0.) then
                rv=2*dh/(sqrt((ruptime(i+1,k)-ruptime(i,k))**2+(ruptime(i,k+1)-ruptime(i,k))**2) &
                  +sqrt((ruptime(i,k)-ruptime(i-1,k))**2+(ruptime(i,k)-ruptime(i,k-1))**2))
                tint=tint+rv*sqrt(slipX(i,k)**2+slipZ(i,k)**2)
                tint2=tint2+sqrt(slipX(i,k)**2+slipZ(i,k)**2)
              else
                rv = 0.
              endif
              if (sliptime(i,k).ne.1.e4) then
                cz=rv*(sliptime(i,k)-ruptime(i,k))
              endif
            endif
            write (502,*) cz
            write (503,*) rv
          enddo
        enddo
        output_param(1)=tint/tint2
        close(502)
        close(503)
      endif

#if defined DIPSLIP
      slipmax=maxval(slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs))
      output_param(3)=count(slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs)>0.05*slipmax)*dh*dh
#else
      slipmax=maxval(slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs))
      output_param(3)=count(slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs)>0.05*slipmax)*dh*dh
#endif

      deallocate(efrac,erad,slipt,timek)
      deallocate(nkk)

      END SUBROUTINE
