      MODULE displt_com
        real,allocatable,dimension(:,:,:):: u1,v1,w1
      END MODULE

      MODULE strfld_com
        real,allocatable,dimension(:,:,:):: xx,yy,zz,xy,yz,xz
      END MODULE

    MODULE traction_com
      real,allocatable,dimension(:,:):: tx,tz,v1t,avdx,avdz,RFx,RFz  ! x and z component of traction, v component of velocity at fault
      real,allocatable,dimension(:,:):: au1,av1,aw1  !attenuation of velocities near fault
      real :: damp_s
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
    END MODULE
    
    MODULE SlipRates_com
      INTEGER nSR,NL,NW
      REAL dL,dW,dtseis
      REAL M0,Mw
      REAL,allocatable,dimension(:):: MSR,MomentRate

    END MODULE
    
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
      USE traction_com
      USE pml_com
      USE SlipRates_com
      IMPLICIT NONE

      real    :: time,tabs,friction,tmax,xmax,ymax,numer,denom,veltest,dd
      integer :: incrack(nxt,nzt),broken(nxt,nzt),nsurf
      real    :: pdx, pdz
      real    :: U1OUT,SLIPRATEOUT(nxt,nzt)
      real    :: CPUT1,CPUT2
	  REAL :: maxvel,maxvelsave
      real    :: dht, ek, es, ef, c1, c2
      integer :: i,j,it,k, nxe, nxb, nyb, nye, nzb, nze
      real,allocatable,dimension(:)::etot_out, epot_out, ekin_out, efault_out
      integer :: ifrom,ito,jfrom,jto,kk
	  real :: rup_tresh, rv, cz
	  
      rup_tresh=1.e-3
      
      c1  = 9./8.
      c2  = -1./24.

!---------------------------
! Write down the input
!---------------------------
       if (ioutput.eq.1) then
         open(95,file='result/vmodel.inp')
         open(96,file='result/friction.inp')
         do k = nabc+1,nzt-nfs
           do i = nabc+1,nxt-nabc
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
      allocate(u1(nxt,nyt,nzt),v1(nxt,nyt,nzt),w1(nxt,nyt,nzt))
      allocate(xx(nxt,nyt,nzt),yy(nxt,nyt,nzt),zz(nxt,nyt,nzt),xy(nxt,nyt,nzt),yz(nxt,nyt,nzt),xz(nxt,nyt,nzt))
      !pml initialization
      allocate(tx(nxt,nzt),tz(nxt,nzt),v1t(nxt,nzt),avdx(nxt,nzt),avdz(nxt,nzt), RFx(nxt,nzt),RFz(nxt,nzt))
      allocate(etot_out(ntfd), epot_out(ntfd), ekin_out(ntfd), efault_out(ntfd))
      allocate(omega_pml(nabc-1), omegaR_pml(nabc-1),omega_pmlM(nabc-1), omegaR_pmlM(nabc-1), au1(nxt,nzt),av1(nxt,nzt),aw1(nxt,nzt))
      u1=0.;v1=0.;w1=0.
      xx=0.;yy=0.;zz=0.;xy=0.;yz=0.;xz=0.
      ruptime=1.e4;slip=0.;rise=0.;sliptime=1.e4
      broken=0;incrack=0
      tx=0.;tz=0.;v1t=0.;gliss=0.0
      avdx = 0.; avdz = 0.; RFx = 0.; RFz = 0.
      au1=0.; av1=0.;aw1=0 
      MSR=0.;MomentRate=0.
      
      nxb=2
      nxe=nabc
      nyb=nabc+1
      nye=nyt-1
      nzb=nabc+1
      nze=nzt-nfs
      allocate (omegax1(nxe-nxb+1),omegay1(nye-nyb+1+1),omegaz1(nze-nzb+1))      
      allocate (omegaxS1(nxe-nxb+1))    
      allocate (u11(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),u12(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),u13(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (v11(nxe-nxb+1,nye-nyb+1,nze-nzb+1),v12(nxe-nxb+1,nye-nyb+1,nze-nzb+1),v13(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (w11(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),w12(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),w13(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (xx11(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),xx12(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),xx13(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (yy11(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),yy12(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),yy13(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (zz11(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),zz12(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),zz13(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (xz11(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),xz12(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (xy11(nxe-nxb+1,nye-nyb+1,nze-nzb+1),xy12(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (yz11(nxe-nxb+1,nye-nyb+1,nze-nzb+1),yz12(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      
      u11=0.;u12=0.;u13=0.;v11=0.;v12=0.;v13=0.;w11=0.;w12=0.;w13=0.
      xx11=0.;xx12=0.;xx13=0.;yy11=0.;yy12=0.;yy13=0.;zz11=0.;zz12=0.;zz13=0.
      xy11=0.;xy12=0.;xz11=0.;xz12=0.;yz11=0.;yz12=0.
      
      nxb=nxt-nabc+1
      nxe=nxt-1
      nyb=nabc+1
      nye=nyt-1
      nzb=nabc+1
      nze=nzt-nfs
      allocate (omegax2(nxe-nxb+1),omegay2(nye-nyb+1+1),omegaz2(nze-nzb+1))  
      allocate (omegaxS2(nxe-nxb+1))      
      allocate (u21(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),u22(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),u23(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (v21(nxe-nxb+1,nye-nyb+1,nze-nzb+1),v22(nxe-nxb+1,nye-nyb+1,nze-nzb+1),v23(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (w21(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),w22(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),w23(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (xx21(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),xx22(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),xx23(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (yy21(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),yy22(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),yy23(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (zz21(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),zz22(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),zz23(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (xz21(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),xz22(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (xy21(nxe-nxb+1,nye-nyb+1,nze-nzb+1),xy22(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (yz21(nxe-nxb+1,nye-nyb+1,nze-nzb+1),yz22(nxe-nxb+1,nye-nyb+1,nze-nzb+1))  
      
      
      u21=0.;u22=0.;u23=0.;v21=0.;v22=0.;v23=0.;w21=0.;w22=0.;w23=0.
      xx21=0.;xx22=0.;xx23=0.;yy21=0.;yy22=0.;yy23=0.;zz21=0.;zz22=0.;zz23=0.
      xy21=0.;xy22=0.;xz21=0.;xz22=0.;yz21=0.;yz22=0.
      
      nxb=2
      nxe=nxt-1
      nyb=2
      nye=nabc
      nzb=nabc+1
      nze=nzt-nfs
      allocate (omegax3(nxe-nxb+1),omegay3(nye-nyb+1),omegaz3(nze-nzb+1)) 
      allocate (omegaxS3(nxe-nxb+1),omegayS3(nye-nyb+1)) 
      allocate (u31(nxe-nxb+1,nye-nyb+1,nze-nzb+1),u32(nxe-nxb+1,nye-nyb+1,nze-nzb+1),u33(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (v31(nxe-nxb+1,nye-nyb+1,nze-nzb+1),v32(nxe-nxb+1,nye-nyb+1,nze-nzb+1),v33(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (w31(nxe-nxb+1,nye-nyb+1,nze-nzb+1),w32(nxe-nxb+1,nye-nyb+1,nze-nzb+1),w33(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (xx31(nxe-nxb+1,nye-nyb+1,nze-nzb+1),xx32(nxe-nxb+1,nye-nyb+1,nze-nzb+1),xx33(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (yy31(nxe-nxb+1,nye-nyb+1,nze-nzb+1),yy32(nxe-nxb+1,nye-nyb+1,nze-nzb+1),yy33(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (zz31(nxe-nxb+1,nye-nyb+1,nze-nzb+1),zz32(nxe-nxb+1,nye-nyb+1,nze-nzb+1),zz33(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (xz31(nxe-nxb+1,nye-nyb+1,nze-nzb+1),xz32(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (xy31(nxe-nxb+1,nye-nyb+1,nze-nzb+1),xy32(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (yz31(nxe-nxb+1,nye-nyb+1,nze-nzb+1),yz32(nxe-nxb+1,nye-nyb+1,nze-nzb+1))  
      
      u31=0.;u32=0.;u33=0.;v31=0.;v32=0.;v33=0.;w31=0.;w32=0.;w33=0.
      xx31=0.;xx32=0.;xx33=0.;yy31=0.;yy32=0.;yy33=0.;zz31=0.;zz32=0.;zz33=0.
      xy31=0.;xy32=0.;xz31=0.;xz32=0.;yz31=0.;yz32=0.
      
      nxb=2
      nxe=nxt-1
      nyb=2
      nye=nyt-1
      nzb=2
      nze=nabc
      allocate (omegax4(nxe-nxb+1),omegay4(nye-nyb+1+1),omegaz4(nze-nzb+1))  
      allocate (omegaxS4(nxe-nxb+1),omegayS4(nye-nyb+1+1),omegazS4(nze-nzb+1))  
      allocate (u41(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),u42(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),u43(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (v41(nxe-nxb+1,nye-nyb+1,nze-nzb+1),v42(nxe-nxb+1,nye-nyb+1,nze-nzb+1),v43(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (w41(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),w42(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),w43(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (xx41(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),xx42(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),xx43(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (yy41(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),yy42(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),yy43(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (zz41(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),zz42(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),zz43(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (xz41(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1),xz42(nxe-nxb+1,nye-nyb+1+1,nze-nzb+1))
      allocate (xy41(nxe-nxb+1,nye-nyb+1,nze-nzb+1),xy42(nxe-nxb+1,nye-nyb+1,nze-nzb+1))
      allocate (yz41(nxe-nxb+1,nye-nyb+1,nze-nzb+1),yz42(nxe-nxb+1,nye-nyb+1,nze-nzb+1))    

      u41=0.;u42=0.;u43=0.;v41=0.;v42=0.;v43=0.;w41=0.;w42=0.;w43=0.
      xx41=0.;xx42=0.;xx43=0.;yy41=0.;yy42=0.;yy43=0.;zz41=0.;zz42=0.;zz43=0.
      xy41=0.;xy42=0.;xz41=0.;xz42=0.;yz41=0.;yz42=0.

     ! pi = 4.0*atan(1.0)

      do i = 1,nabc-1
        omega_pml(i) = 0.5*dt*omegaM_pml * (real(i)/real((nabc-1)))**4
        omegaR_pml(nabc-i) = omega_pml(i)
        omega_pmlM(i) = 0.5*dt*omegaM_pml * ((real(i)-0.5)/real((nabc-1)))**4
        omegaR_pmlM(nabc-i) = omega_pmlM(i)
       ! print*, omega_pml(i), omega_pmlM(i)
      end do
      
      omegax1=0.;omegax2=0.;omegax3=0.;omegax4=0.
      omegay1=0.;omegay2=0.;omegay3=0.;omegay4=0.
      omegaz1=0.;omegaz2=0.;omegaz3=0.;omegaz4=0.
      omegaxS1=0.;omegaxS2=0.;omegaxS3=0.;omegaxS4=0.
      omegayS3=0.;omegayS4=0.
      omegazS4=0.
      do i=1,nabc-1
        omegax1(i)= omegaR_pml(i)
        omegaxS1(i)= omegaR_pmlM(i)
        
        omegax2(i)= omega_pml(i)
        omegaxS2(i)= omega_pmlM(i)
        
        omegay3(i)= omegaR_pml(i)
        omegax3(i)= omegaR_pml(i)
        omegax3(nxt-i-1)= omega_pmlM(i)   
        omegayS3(i)= omegaR_pmlM(i)
        omegaxS3(i)= omegaR_pmlM(i)
        omegaxS3(nxt-i-1)= omega_pml(i)  
        
        omegaz4(i)= omegaR_pml(i)
        omegax4(i)= omegaR_pml(i)
        omegax4(nxt-i-1)= omega_pmlM(i)
        omegay4(i)= omegaR_pml(i)
        
        omegazS4(i)= omegaR_pmlM(i)
        omegaxS4(i)= omegaR_pmlM(i)
        omegaxS4(nxt-i-1)= omega_pml(i)
        omegayS4(i)= omegaR_pmlM(i)
      enddo
      
!-------------------------------------------------------
!     Loop over time
!-------------------------------------------------------

      if(ioutput.eq.1) then
        OPEN(25, file='result/sliprate.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
        OPEN(26, file='result/shearstress.res',FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
      endif
      CALL CPU_TIME(CPUT1)
      maxvelsave=0.

      !$ACC DATA COPYIN (LAM1,MU1,D1) &
      !$ACC      COPYIN (U1,V1,W1) COPYIN (XX,YY,ZZ,XY,YZ,XZ) &
      !$ACC      COPYIN (tx,tz,v1t,avdx,avdz,RFx,RFz) &
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
      !$ACC      COPYIN (BROKEN,dyn_xz,gliss,STRINIX,PEAK_XZ,DC) &
      !$ACC      COPY (INCRACK,RUPTIME,SLIPTIME,SLIP)
      dht = dh/dt
      do it = 1,ntfd
        time = it*dt
         SLIPRATEOUT(:,:)=0.
         SCHANGE(:,:)=0.
!$ACC DATA COPY (SLIPRATEOUT,SCHANGE)

!-------------------------------------------------------------
!   Apply 4th-order differencing to interior particle velocities
!             from  2:nxt;2:nyt-2;2:nzt
!-------------------------------------------------------------
         call dvel(nxt,nyt-2,nzt,dt,dh)

!-------------------------------------------------------------
!   Compute velocities of 2nd order accuracy near fault
!-------------------------------------------------------------
         call bnd2d(nxt,nyt-2,nzt,dt,dh)
!-------------------------------------------------------------
!   Compute velocities at fault boundary
!-------------------------------------------------------------
         call tasu1(nxt,nysc,nzt,dt,dh)
         call tasw1(nxt,nysc,nzt,dt,dh)

!Stress free and absorbing conditions
         call pml_uvw (nxt,nyt,nzt,dt,dh)
         
		 call fuvw(nxt,nyt,nzt-nfs)

      !$ACC PARALLEL DEFAULT (PRESENT)
      !$ACC LOOP GANG
         do i = 1,nxt
      !$ACC LOOP VECTOR
            do k = 1,nzt
               v1(i,nyt,k)=v1(i,nyt-1,k)
            enddo
         enddo
      !$ACC END PARALLEL
         
         
!-------------------------------------------------------------
!   4th-order differencing of pressure
!-------------------------------------------------------------
         call dstres(nxt,nyt-2,nzt,dh,dt)
!-------------------------------------------------------------
!   Compute stress tensor of 2nd order accuracy near fault
!-------------------------------------------------------------
         call strbnd(nxt,nyt-2,nzt,dh,dt)
!-------------------------------------------------------------
!   Compute stress components at fault boundary
!-------------------------------------------------------------
        call tasxz(nxt,nysc,nzt,dt,dh)
        call tasii(nxt,nysc,nzt,dt,dh)

!Stress free and absorbing conditions
        call pml_xyz (nxt,nyt,nzt,dt,dh)
       
		call fres(nxt,nyt,nzt-nfs)   !TADY PADA GPU!!!!

        !$ACC PARALLEL DEFAULT (PRESENT)
        !$ACC LOOP GANG
        do i = 1,nxt
          !$ACC LOOP VECTOR
          do k = 1,nzt
            xy(i,nyt,k)=xy(i,nyt-1,k)
            yz(i,nyt,k)=yz(i,nyt-1,k)
          enddo
        enddo
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
        do k = nabc+1,nzt-nfs
          !$ACC LOOP VECTOR
          do i = nabc+1,nxt-nabc
		!pdz = (+(xz(i+1,nysc,k) - xz(i,nysc,k))/2 + (zz(i,nysc,k+1) - zz(i,nysc,k))/2 - yz(i,nysc-1,k))
             pdz = (c1*(xz(i+1,nysc,k) - xz(i,nysc,k))   +  c2*(xz(i+2,nysc,k) - xz(i-1,nysc,k)))/2 + &
	         (c1*(zz(i,nysc,k+1) - zz(i,nysc,k))   +  c2*(zz(i,nysc,k+2) - zz(i,nysc,k-1)))/2 - &
			 yz(i,nysc-1,k)
			
            avdz(i,k)= damp_s*(pdz - avdz(i,k))
            RFz(i,k) = pdz + avdz(i,k)
            tz(i,k) = -RFz(i,k) - 0.5*d1(i,nysc,k)*dht*w1(i,nysc,k)
            avdz(i,k) = pdz

            !pdx = (+(xx(i,nysc,k) - xx(i-1,nysc,k))/2 + (xz(i,nysc,k) - xz(i,nysc,k-1))/2 - xy(i,nysc-1,k))
			 pdx = (c1*(xx(i,nysc,k)   - xx(i-1,nysc,k)) +  c2*(xx(i+1,nysc,k) - xx(i-2,nysc,k)))/2 + &
			 (c1*(xz(i,nysc,k)   - xz(i,nysc,k-1)) +  c2*(xz(i,nysc,k+1) - xz(i,nysc,k-2)))/2 - &
			 xy(i,nysc-1,k)
            avdx(i,k)= damp_s*(pdx - avdx(i,k))
            RFx(i,k) = pdx + avdx(i,k)
            tx(i,k) = -RFx(i,k) - 0.5*d1(i,nysc,k)*dht*u1(i,nysc,k)
            avdx(i,k) = pdx
          enddo
        enddo
        !$ACC END PARALLEL
!-------------------------------------------------------------
!  Apply Boundary Conditions
!-------------------------------------------------------------
        !$ACC PARALLEL DEFAULT (PRESENT)
        !$ACC LOOP GANG
        do k = nabc+1,nzt-nfs
          !$ACC LOOP VECTOR
          do i = nabc+1,nxt-nabc


#if defined DIPSLIP
            tabs = tz(i,k) + strinix(i,k)
            u1out=W1(I,NYSC,K)
#else
            tabs = tx(i,k) + strinix(i,k)
            u1out=U1(I,NYSC,K)
#endif

!----------------------------------------------------------------------------
!   Apply friction law to absolute slip, absolute sliprate, absolute traction
!   Classical slip weakening law
!----------------------------------------------------------------------------

            if (tabs.gt.peak_xz(i,k)) then
              broken(i,k)  = 1
              incrack(i,k) = 1
            endif

            if (broken(i,k) .eq. 1) then
              if (Dc(i,k)<=0.) then
                dd = 0.
              else
                if (gliss(i,k).le.Dc(i,k)) then
                  friction = peak_xz(i,k) * (1.0 - gliss(i,k)/Dc(i,k)) + dyn_xz(i,k)*gliss(i,k)/Dc(i,k)
                else
                  friction = dyn_xz(i,k)
                endif
              endif

              if ((ruptime(i,k).eq.1.e4).AND.(- 2*u1out>rup_tresh)) ruptime(i,k) = time
              
#if defined DIPSLIP
              !if (u1out .ge. 0.0) then
              if (tabs .ge. friction) then
                slip(i,k) = slip(i,k)  - 2*u1out*dt
                gliss(i,k) = gliss(i,k)  - 2*u1out*dt
                sliprateout(i,k) = - 2.*u1out
                tz(i,k) = + friction - strinix(i,k)

              !else
              !  broken(i,k) = 0
              !  gliss(i,k) = 0.0
              endif
#else
              !if (u1out .ge. 0.0) then
              if (tabs .ge. friction) then
                slip(i,k) = slip(i,k)  - 2*u1out*dt
                gliss(i,k) = gliss(i,k)  - 2*u1out*dt
                sliprateout(i,k) = - 2.*u1out
                tx(i,k) = + friction - strinix(i,k)

              !else
              !  broken(i,k) = 0
              !  gliss(i,k) = 0.0
              endif
#endif

              if ((sliptime(i,k)==1.e4).AND.(slip(i,k)>Dc(i,k))) sliptime(i,k)=time
            endif
#if defined DIPSLIP
          SCHANGE(I,K) = tz(i,k) + strinix(i,k)
#else
          SCHANGE(I,K) = tx(i,k) + strinix(i,k)

#endif
          enddo
        enddo
        !$ACC END PARALLEL

!Stress free
        call fres(nxt,nyt,nzt-nfs)
     !   energy
         ! ek=0.
         ! es=0.
		 ! ef= 0.
         ! do k = nabc+1,nzt-nfs
             ! do j = nabc+1,nyt-1
	    	 ! do i = nabc+1,nxt-nabc                
                     ! ek = ek + d1(i,j,k)*(u1(i,j,k)**2+v1(i,j,k)**2+w1(i,j,k)**2)
                     ! es = es + xx(i,j,k)*(u1(i+1,j,k)-u1(i,j,k)) + yy(i,j,k)*(v1(i,j,k)-v1(i,j-1,k)) + zz(i,j,k)*(w1(i,j,k)-w1(i,j,k-1)) + &
                         ! xy(i,j,k)*(u1(i,j+1,k)-u1(i,j,k) + v1(i,j,k) - v1(i-1,j,k)) +   &
                         ! xz(i,j,k)*(u1(i,j,k+1) - u1(i,j,k) + w1(i,j,k) - w1(i-1,j,k)) + &
                         ! yz(i,j,k)*(v1(i,j,k+1) - v1(i,j,k) + w1(i,j+1,k) - w1(i,j,k))
                 ! enddo
             ! enddo
         ! enddo
         ! do k = nabc+1,nzt-nfs
	    	 ! do i = nabc+1,nxt-nabc 
		     ! ef = ef + (schange(i,k)-dyn_xz(i,k))*slip(i,k)
                 ! enddo
         ! enddo
		
         ! etot_out(it)=ek + es/dh + ef*dh**2
         ! epot_out(it)=ek
		 ! ekin_out(it)=es/dh
		 ! efault_out(it)=ef*dh**2

        !$ACC END DATA

      if(ioutput.eq.1) then
        WRITE(25) sliprateout(nabc+1:nxt-nabc,nabc+1:nzt-nfs)
        WRITE(26) SCHANGE(nabc+1:nxt-nabc,nabc+1:nzt-nfs)
      endif
        
      k=int(real(it-1)*dt/dtseis)+1
      if(k<1.or.k>nSR)write(*,*)'CHYBA!',k
      do j=1,NW
        jto=max(1,int(dW/dh*j))+1+nabc
        jfrom=min(jto,int(dW/dh*(j-1))+1)+1+nabc
        do i=1,NL
          ifrom=int(dL/dh*(i-1))+1+nabc
          ito=int(dL/dh*i)+nabc
          kk=((j-1)*NL+i-1)*nSR+k
          MSR(kk)=MSR(kk)+sum(sliprateout(ifrom:ito,jfrom:jto))/dble((ito-ifrom+1)*(jto-jfrom+1)*(dtseis/dt))
        enddo
      enddo
      MomentRate(k)=MomentRate(k)+sum(sliprateout(:,:)*mu1(:,nysc,:))*dh*dh/(dtseis/dt)
    
        if(mod(it,int(1./dt))==0)then
          maxvel=maxval(sliprateout(nabc+1:nxt-nabc,nabc+1:nzt-nfs))
          write(*,*)'Time: ',time,'Slip rate max: ',maxvel
          if(maxvel>maxvelsave)maxvelsave=maxvel
          if (maxvel<=0.01*maxvelsave)exit
        endif

      enddo ! --- End of the time loop
    
      !$ACC END DATA

      SCHANGE(:,:)=SCHANGE(:,:)-strinix(:,:)   !stress drop
      
      deallocate(u1,v1,w1)
      deallocate(xx,yy,zz,xy,yz,xz)
      deallocate(tx,tz,v1t,avdx,avdz, RFx,RFz)

      deallocate(omega_pml,omegaR_pml,omega_pmlM,omegaR_pmlM,au1,av1,aw1)
      deallocate (omegax1,omegay1,omegaz1,omegaxS1)    
      deallocate (u11,u12,u13,v11,v12,v13,w11,w12,w13)
      deallocate (xx11,xx12,xx13,yy11,yy12,yy13,zz11,zz12,zz13,xz11,xz12,xy11,xy12,yz11,yz12)
      deallocate (omegax2,omegay2,omegaz2,omegaxS2)      
      deallocate (u21,u22,u23,v21,v22,v23,w21,w22,w23)
      deallocate (xx21,xx22,xx23,yy21,yy22,yy23,zz21,zz22,zz23,xz21,xz22,xy21,xy22,yz21,yz22)  
      deallocate (omegax3,omegay3,omegaz3,omegaxS3,omegayS3) 
      deallocate (u31,u32,u33,v31,v32,v33,w31,w32,w33)
      deallocate (xx31,xx32,xx33,yy31,yy32,yy33,zz31,zz32,zz33,xz31,xz32,xy31,xy32,yz31,yz32)  
      deallocate (omegax4,omegay4,omegaz4,omegaxS4,omegayS4,omegazS4)
      deallocate (u41,u42,u43,v41,v42,v43,w41,w42,w43)
      deallocate (xx41,xx42,xx43,yy41,yy42,yy43,zz41,zz42,zz43,xz41,xz42,xy41,xy42,yz41,yz42)    

      CALL CPU_TIME(CPUT2)
      PRINT *,'CPU TIME OF TIME LOOP: ',CPUT2-CPUT1


!-------------------
! Open output files:
!-------------------
      if(ioutput.eq.1) then
        close(25)
        close(26)
        do it = 1,ntfd
           !write(410) etot_out(it),ekin_out(it),epot_out(it),efault_out(it)
        enddo
         !close(410)
      endif
      deallocate(etot_out,epot_out,ekin_out,efault_out)

      tmax            = -1.
      output_param(1) =  0.
      nsurf           =  0
      output_param(2) =  0.
      numer           =  0.
      denom           =  0.
      do k = nabc+1,nzt-nfs	
        do i = nabc+1,nxt-nabc
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
          output_param(2) = output_param(2) + slip(i,k)*mu1(i,nysc,k)*(dh*dh)

           ! --- Rise time:  MUSI SE PREDEFINOVAT
!           if (maxval(sliprate(i,k,:)).eq.0) then
!             rise(i,k) = 0.
!           elseif (ruptime(i,k).ne.0.) then
!             rise(i,k) = slip(i,k)/maxval(sliprate(i,k,:))
!           else
!             rise(i,k) = 0.
!           endif

           ! --- Stress drop:
           if (ruptime(i,k).ne.0.) then
             numer = numer + schange(i,k)*slip(i,k)
             denom = denom + slip(i,k)
           endif
           if (denom .ne. 0.0) then
             output_param(4) = -(numer/denom)
           else
             output_param(4) = 0.0
           endif

         enddo
       enddo
       output_param(5) = (1./2.)**sum(peak_xz*Dc)/dble((nxt-2*nabc)*(nzt-nfs-nabc))
       output_param(6) = (1./2.)*output_param(4)*(output_param(2)/(mu_mean*output_param(3)))
       M0=output_param(2)

!---------------------------
! Write down the output
!---------------------------
       if (ioutput.eq.1) then
!         open(96,file='result/risetime.res')
         open(97,file='result/ruptime.res')
         open(98,file='result/slip.res')
         open(99,file='result/stressdrop.res')
         do k = nabc+1,nzt-nfs
           do i = nabc+1,nxt-nabc
!             write(96,*) rise(i,k)
             write(97,*) ruptime(i,k)
             write(98,*) slip(i,k)
             write(99,*) schange(i,k)
           enddo
         enddo
!         close(96)
         close(97)
         close(98)
         close(99)

         open(501,file='result/contour.res')
         write(501,*) 'j k t'
         do k = nabc+1,nzt-nfs
           do i = nabc+1,nxt-nabc
             write (501,*) (real(i-nabc-1) - real(nxt-2*nabc)/2.)*dh,real(k-nabc-1)*dh,ruptime(i,k) 
           enddo
         enddo		
         close(501)
         open(502,file='result/czone.res')
         open(503,file='result/rvel.res')
         do k = nabc+1,nzt-nfs
           do i = nabc+1,nxt-nabc
             cz=0.
             rv=0.
             if ((ruptime(i,k).ne.1.e4).and.(ruptime(i+1,k).ne.1.e4).and.(ruptime(i,k+1).ne.1.e4).and.(ruptime(i-1,k).ne.1.e4).and.(ruptime(i,k-1).ne.1.e4)) then
               if (sliptime(i,k).ne.1.e4) then
                 rv = (sqrt((ruptime(i+1,k)-ruptime(i,k))**2+(ruptime(i,k+1)-ruptime(i,k))**2)+sqrt((ruptime(i,k)-ruptime(i-1,k))**2+(ruptime(i,k)-ruptime(i,k-1))**2))
                 if (rv.ne.0.) then
                   rv=2*dh/(sqrt((ruptime(i+1,k)-ruptime(i,k))**2+(ruptime(i,k+1)-ruptime(i,k))**2)+sqrt((ruptime(i,k)-ruptime(i-1,k))**2+(ruptime(i,k)-ruptime(i,k-1))**2))
                 else
                   rv = 0.
                 endif
                 cz=rv*(sliptime(i,k)-ruptime(i,k))
               endif
             endif
             write (502,*) cz
             write (503,*) rv
           enddo
         enddo
         close(502)
         close(503)
       endif

       END
