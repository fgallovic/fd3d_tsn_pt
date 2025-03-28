!-------------------------------------------------------
! Processing of Bayesian posterior samples from for fd3d_tsn_pt, considering multiple sampls*.dat a removing burn-in phase
!-------------------------------------------------------
! Authors: Frantisek Gallovic, Jan Premus, �. Valentov� (10/2024)
! Charles University in Prague, Faculty of Mathematics and Physics

! This code is published under the GNU General Public License. To any
! licensee is given permission to modify the work, as well as to copy
! and redistribute the work or any derivative version. Still we would
! like to kindly ask you to acknowledge the authors and don't remove
! their names from the code. This code is distributed in the hope
! that it will be useful, but WITHOUT ANY WARRANTY.
! ------------------------------------------------------

    MODULE interp_pb
      INTEGER NLI,NWI
    END MODULE

    MODULE fd3dparam_pb
      integer :: nxt,nyt,nzt,ntfd,nysc
      real    :: dh,dt
    END MODULE

    MODULE mask_pb
      REAL,ALLOCATABLE,DIMENSION(:,:):: mask
    END MODULE
    
    MODULE friction_pb
      USE fd3dparam_pb
      real,allocatable,dimension(:,:,:):: strinix,peak_xz,Dc
      real,parameter:: pi=3.1415926535
      REAL dip
      CONTAINS

      FUNCTION normstress(j)
      IMPLICIT NONE
      real:: normstress
      integer:: j
#if defined FSPACE
!      normstress=100.e6
      normstress=2.8e9
!      normstress=21.e9
#else
#if defined DIPSLIP
      normstress=max(1.e5,8520.*dh*real(nzt-j)*sin(dip/180.*pi))
      !normstress=min(18.*dh*real(nzt-j)*sin(dip/180.*pi)/1.e3,100.);normstress=1.e6*max(1.,normstress)
#else
      normstress=max(1.e5,16200.*dh*real(nzt-j)*sin(dip/180.*pi))
      !normstress=min(18.*dh*real(nzt-j)*sin(dip/180.*pi)/1.e3,100.);normstress=1.e6*max(1.,normstress)
      !normstress=min(18.*dh*real(nzt-j)*sin(dip/180.*pi)/1.e3,100.);normstress=1.e6*max(1.,normstress)
      !normstress=min(18.*dh*real(nzt-j)*sin(dip/180.*pi)/1.e3,60.);normstress=1.e6*max(1.,normstress)
#endif
#endif
      END FUNCTION
    END MODULE

    MODULE medium_pb
      real,allocatable,dimension(:,:,:):: lam1,mu1,d1
      real:: mu_mean
    END MODULE

    PROGRAM ProcessBayes
    USE fd3dparam_pb
    USE friction_pb
    USE interp_pb
    USE medium_pb
    USE mask_pb
    IMPLICIT NONE
    INTEGER,PARAMETER:: NMAX=1e6
    INTEGER,PARAMETER:: nper=56    !Three lines related to the smoothing of the moment rate spectra
    REAL,PARAMETER:: fc1=0.5286,fc2=24.5279
    REAL per(NMAX)
    INTEGER,ALLOCATABLE,DIMENSION(:):: accepted,iknow,kchain
    REAL,ALLOCATABLE,DIMENSION(:):: normalstress,VRs,misfits,meansd,meansl,duration,nuclsize,EG,ER,RE,meanoverstress,M0,meanDc,meanStrengthExcess,meanslip,rupturearea,meanruptvel,meanstrength,EGrate,VRgps,x0,z0
    REAL,ALLOCATABLE,DIMENSION(:,:,:):: DcA,TsA,T0A,SEA,ruptime1,slip1,rise1,schange1,es1,strengthexcess1,rupvel1,EGfault
    REAL,ALLOCATABLE,DIMENSION(:,:):: dum11,dum12,dum13,dum21,dum22,dum23,dum24,ms1
    real, allocatable, dimension(:) :: MRate,fc
    real, allocatable, dimension(:,:) :: MomentRate,MomentSpec
    INTEGER, allocatable, dimension(:) :: indx
    REAL*8,ALLOCATABLE,DIMENSION(:,:) :: dumFFT,strinixSpatialmean,peak_xzSpatialmean,Dcspatialmean
    COMPLEX*16,ALLOCATABLE :: dumFFTq(:),DcFFTc(:,:,:),strinixFFTc(:,:,:),peak_xzFFTc(:,:,:)!,DcFFTq(:,:)
    INTEGER nxtfft,nztfft
    REAL kx,ky,krad,strinixpwr,peak_xzpwr,Dcpwr,strinixpeak_xzxcorr,strinixDcxcorr,peak_xzDcxcorr
    REAL M0dum,EGdum,ERdum,VRgpsdum
    REAL bestmisfit,misfitaccept,dum,coh,dyn,dumarr(8)
    REAL vr,mf,x,z,slipmax
    REAL M0toMw,Mw,maxmr,dur1,dur2
    complex, allocatable :: cmrate(:)
    INTEGER NM,NTOT,ndoub,nchains,kall
    INTEGER i,j,k,kk,ml(1),ncent,ichain,ierr
    integer :: nsr,np
    real :: T,SRdur,dtseis,dfseis
    character(50) :: fname
!    real :: burn=0.1
!    integer :: idown=5
    real :: burn=0.2,slipfact=0.05
    integer :: idown=1

!--------------------
! Read the input file
!--------------------
    write(*,*)'Reading FD3D parameters...'
    open(11, file='inputfd3d.dat', status='old')
    read(11,*) nxt,nyt,nzt
    read(11,*) dh
    read(11,*) ntfd
    read(11,*) dt
    read(11,*) dip
    close(11)

    open(10,file='input.dat',action='read')
    read(10,*)
    read(10,*) !nfmax
    read(10,*)
    read(10,*) T,SRdur !,T1,T2
    read(10,*)
    read(10,*) !artifDT
    read(10,*)
    read(10,*) !NRseis
    read(10,*)
    read(10,*) !NL,NW
    read(10,*)
    read(10,*) !M0aprior
    read(10,*)
    read(10,*) !strike,dip
    read(10,*)
    read(10,*) !hypodepth
    read(10,*)
    read(10,*) !leng,widt
    read(10,*)
    read(10,*) !epicL,epicW
    read(10,*)
    read(10,*) np
    close(10)

    dtseis=T/real(np)
    dfseis=1./T
!    nSR=ceiling(real(ntfd)/(dtseis/dt))
    nSR=ceiling(SRdur/dtseis)
    allocate(MRate(nSR))

    allocate(lam1(nxt,nyt,nzt),mu1(nxt,nyt,nzt),d1(nxt,nyt,nzt))

    CALL readcrustalmodel(dip)

    open(10,FILE='inputinv.dat')
    read(10,*)
    read(10,*)NLI,NWI
    
    ALLOCATE(misfits(NMAX),VRs(NMAX),accepted(NMAX),iknow(NMAX))
    ALLOCATE(dum11(NLI,NWI),dum12(NLI,NWI),dum13(NLI,NWI),dum21(nxt,nzt),dum22(nxt,nzt),dum23(nxt,nzt),dum24(nxt,nzt))
    ALLOCATE(normalstress(nzt))
    do i=1,nzt
       normalstress(i)=normstress(i)
    enddo

!------ Learn about misfits
    accepted=0;iknow=0
    write(*,*)'Set number of chains to read sampls files. If set to 0, will read the joint sampls.dat'
    read(*,*) nchains

    if (nchains==0) then    !READ FROM A SINGLE FILE
      k=0
      open(101,FILE='sampls.dat',FORM='UNFORMATTED',ACCESS='STREAM')
10    read(101,END=11,ERR=11)mf,vr,dum11(:,:),dum12(:,:),dum13(:,:),dum21(:,:),dum22(:,:),dum23(:,:),dum24(:,:),MRate(:),dumarr(1:5)
      k=k+1
      misfits(k)=mf
      if(k==NMAX)stop 'Increase dimension!'
      goto 10
11    close(101)
      NTOT=k
      bestmisfit=minval(misfits(1:NTOT))
      write(*,*)'Total number of models: ',NTOT
      write(*,*)'Best misfit: ',bestmisfit
!      misfitaccept=maxval(misfits(1:ntot))
!      misfitaccept=bestmisfit-log(0.02) !Probability threashold (2% for real data)
!      misfitaccept=bestmisfit-log(0.05) !Probability threashold (5% for real data)
!      misfitaccept=bestmisfit-log(0.01) !Probability threashold (1% for real data)
      misfitaccept=bestmisfit-log(0.001) !Probability threashold (1%% for inv1)
!     misfitaccept=bestmisfit-log(0.00001) !Probability threashold (.1%% for pga)
!     misfitaccept=bestmisfit+20. !Honzuv napad
      allocate(T0A(NLI,NWI,NTOT),TsA(NLI,NWI,NTOT),DcA(NLI,NWI,NTOT))
      open(101,FILE='sampls.dat',FORM='UNFORMATTED',ACCESS='STREAM')
      do i=1,NTOT
        read(101)mf,vr,dum11(:,:),dum12(:,:),dum13(:,:),dum21(:,:),dum22(:,:),dum23(:,:),dum24(:,:),MRate(:),dumarr(1:5)
        do k=1,i-1
          if (all(abs((dum11(:,:)-T0A(:,:,k)))<1.e-5) .and. all(abs((dum12(:,:)-TsA(:,:,k)))<1.e-5) .and. all(abs((dum13(:,:)-DcA(:,:,k)))<1.e-5))then
            iknow(i)=1
            exit
          endif
        enddo
        if(mf<=misfitaccept.and.iknow(i)==0)accepted(i)=1
        T0A(:,:,i)=dum11(:,:)
        TsA(:,:,i)=dum12(:,:)
        DcA(:,:,i)=dum13(:,:)
      enddo
      NM=sum(accepted)
      ndoub=sum(iknow)
      write(*,*)'Removed doubles:',ndoub
      write(*,*)'Number of accepted models: ',NM
      DEALLOCATE(T0A,TsA,DcA)
      DEALLOCATE(misfits,VRs)
      close(101)
      
    else    !READ FROM MULTIPLE FILES
      kall=0
      allocate(kchain(nchains))
      do ichain=1,nchains
        write(fname,'(a,i3.3)') 'sampls',ichain
        write(*,*)fname
        open(101,FILE=fname,FORM='UNFORMATTED',ACCESS='STREAM')
        ierr=0
        k=0
        do while(ierr==0)
          read(101,iostat=ierr)mf,vr,dum11(:,:),dum12(:,:),dum13(:,:),dum21(:,:),dum22(:,:),dum23(:,:),dum24(:,:),MRate(:),dumarr(1:5)
          k=k+1
          kall=kall+1
          misfits(kall)=mf
        enddo
        close(101)
        kchain(ichain)=k-1
        kall=kall-1
      enddo
      NTOT=kall !=sum(kchain(:))
      bestmisfit=minval(misfits(1:NTOT))
      write(*,*)'Total number of models: ',NTOT
      write(*,*)'Best misfit: ',bestmisfit
      !reduce samples by the burn-in phase and by idown on each chain
      kall=0;kk=0
      do ichain=1,nchains
        do j=1,kchain(ichain) !NTOT
           kk=kk+1
           if (j>nint(burn*kchain(ichain)) .and. mod(kk,idown)==0)kall=kall+1
        enddo
      enddo
      NTOT=kall
!      NTOT=ceiling(sum(nint((1.-burn)*kchain(:)))/real(idown))
      write(*,*)'Total number of models after reduction: ',NTOT
!      misfitaccept=maxval(misfits(1:ntot))
!      misfitaccept=bestmisfit-log(0.02) !Probability threashold (2% for real data)
!      misfitaccept=bestmisfit-log(0.05) !Probability threashold (5% for real data)
!      misfitaccept=bestmisfit-log(0.01) !Probability threashold (1% for real data)
      misfitaccept=bestmisfit-log(0.001) !Probability threashold (1%% for inv1)
!     misfitaccept=bestmisfit-log(0.00001) !Probability threashold (.1%% for pga)
!     misfitaccept=bestmisfit+20. !Honzuv napad
      misfitaccept=1.e10
      allocate(T0A(NLI,NWI,NTOT),TsA(NLI,NWI,NTOT),DcA(NLI,NWI,NTOT))
      i=0
      kk=0
      iknow=0
      ndoub=0
      do ichain=1,nchains
        write(fname,'(a,i3.3)') 'sampls',ichain
        open(101,FILE=fname,FORM='UNFORMATTED',ACCESS='STREAM')
        do j=1,kchain(ichain) !NTOT
          kk=kk+1
          read(101)mf,vr,dum11(:,:),dum12(:,:),dum13(:,:),dum21(:,:),dum22(:,:),dum23(:,:),dum24(:,:),MRate(:),dumarr(1:5)
          if (j>nint(burn*kchain(ichain)) .and. mod(kk,idown)==0) then
            i=i+1
            do k=1,i-1
              if (all(abs((dum11(:,:)-T0A(:,:,k)))<1.e-5) .and. all(abs((dum12(:,:)-TsA(:,:,k)))<1.e-5) .and. all(abs((dum13(:,:)-DcA(:,:,k)))<1.e-5))then
                iknow(kk)=1
                exit
              endif
            enddo
            if(mf<=misfitaccept.and.iknow(kk)==0)accepted(kk)=1
            T0A(:,:,i)=dum11(:,:)
            TsA(:,:,i)=dum12(:,:)
            DcA(:,:,i)=dum13(:,:)
          endif
        enddo
        close(101)
      enddo
      NM=sum(accepted)
      ndoub=sum(iknow)
      write(*,*)'Removed doubles:',ndoub
      write(*,*)'Number of accepted models: ',NM
      DEALLOCATE(T0A,TsA,DcA)
      DEALLOCATE(misfits,VRs)
    endif
      
!------ Read accepted models
    allocate(strinix(nxt,nzt,NM),peak_xz(nxt,nzt,NM),Dc(nxt,nzt,NM))
    ALLOCATE(misfits(NM),VRs(NM),meansd(NM),meansl(NM),duration(NM),nuclsize(NM),EG(NM),ER(NM),RE(NM),meanoverstress(NM),M0(NM),EGrate(NM),VRgps(NM))
    ALLOCATE(meanDc(NM),meanStrengthExcess(NM),meanslip(NM),rupturearea(NM),meanruptvel(NM),meanstrength(NM),x0(NM),z0(NM))
    allocate(DcA(NLI,NWI,NM),TsA(NLI,NWI,NM),T0A(NLI,NWI,NM),SEA(NLI,NWI,NM))
    allocate(ruptime1(nxt,nzt,NM),slip1(nxt,nzt,NM),rise1(nxt,nzt,NM),schange1(nxt,nzt,NM),es1(nxt,nzt,NM),ms1(nxt,nzt),strengthexcess1(nxt,nzt,NM),rupvel1(nxt,nzt,NM),EGfault(nxt,nzt,NM))
    allocate(MomentRate(nSr,NM),MomentSpec(nper,NM),indx(NM),fc(NM))

    if (nchains==0) then    !READ FROM A SINGLE FILE
      open(101,FILE='sampls.dat',FORM='UNFORMATTED',ACCESS='STREAM')
      k=0
      do i=1,NTOT
        read(101)mf,vr,dum11(:,:),dum12(:,:),dum13(:,:),dum21(:,:),dum22(:,:),dum23(:,:),dum24(:,:),Mrate(:),M0dum,EGdum,ERdum,dum,VRgpsdum
        if(accepted(i)==1)then
          k=k+1
          misfits(k)=mf
          VRs(k)=vr
          T0A(:,:,k)=dum11(:,:)
          TsA(:,:,k)=dum12(:,:)
          DcA(:,:,k)=dum13(:,:)
          do kk=1,NWI
            SEA(:,kk,k)=TsA(:,kk,k)*normalstress((kk-1)*((nzt)/(NWI-1))+1)-T0A(:,kk,k)
          enddo
          ruptime1(:,:,k)=dum21(:,:)
          slip1(:,:,k)=dum22(:,:)
          rise1(:,:,k)=dum23(:,:)
          schange1(:,:,k)=dum24(:,:)
          MomentRate(:,k)=Mrate(:)
          M0(k)=M0dum
          Eg(k)=Egdum
          Er(k)=Erdum
          VRgps(k)=VRgpsdum
        endif
      enddo
      close(101)
    else
      k=0
      i=0
      do ichain=1,nchains
        write(fname,'(a,i3.3)') 'sampls',ichain
        open(101,FILE=fname,FORM='UNFORMATTED',ACCESS='STREAM')
        do j=1,kchain(ichain)
          read(101)mf,vr,dum11(:,:),dum12(:,:),dum13(:,:),dum21(:,:),dum22(:,:),dum23(:,:),dum24(:,:),Mrate(:),M0dum,EGdum,ERdum,dum,VRgpsdum
          i=i+1
          if(accepted(i)==1)then
            k=k+1
            misfits(k)=mf
            VRs(k)=vr
            T0A(:,:,k)=dum11(:,:)
            TsA(:,:,k)=dum12(:,:)
            DcA(:,:,k)=dum13(:,:)
            do kk=1,NWI
              SEA(:,kk,k)=TsA(:,kk,k)*normalstress((kk-1)*((nzt)/(NWI-1))+1)-T0A(:,kk,k)
            enddo
            ruptime1(:,:,k)=dum21(:,:)
            slip1(:,:,k)=dum22(:,:)
            rise1(:,:,k)=dum23(:,:)
            schange1(:,:,k)=dum24(:,:)
            MomentRate(:,k)=Mrate(:)
            M0(k)=M0dum
            Eg(k)=Egdum
            Er(k)=Erdum
            VRgps(k)=VRgpsdum
          endif
        enddo
        close(101)
      enddo
     
    endif
! Save best and worst accepted models
    open(201,FILE='forwardmodel.best.dat')
    ml(:)=minloc(misfits(:))
    write(201,'(10000E13.5)')misfits(ml(1)),VRs(ml(1)),T0A(:,:,ml(1)),TsA(:,:,ml(1)),DcA(:,:,ml(1))
    close(201)
    open(201,FILE='forwardmodel.worst.dat')
    ml(:)=maxloc(misfits(:))
    write(201,'(10000E13.5)')misfits(ml(1)),VRs(ml(1)),T0A(:,:,ml(1)),TsA(:,:,ml(1)),DcA(:,:,ml(1))
    close(201)
    
    open(231,FILE='forwardmodel.lowoverstress.dat')
    open(232,FILE='forwardmodel.smallestnucl.dat')

coh=0.5e6

    allocate(cmrate(np))
    do k=1,NM
      CALL interpolate(T0A(:,:,k),strinix(:,:,k))
      CALL interpolate(TSA(:,:,k),ms1(:,:))
      CALL interpolate(DcA(:,:,k),Dc(:,:,k))
      do i=1,nzt
        peak_xz(:,i,k)=ms1(:,i)*normalstress(i)+coh
      enddo
      es1(:,:,k)=(peak_xz(:,:,k)-strinix(:,:,k))/max(1.,strinix(:,:,k))
      slipmax=maxval(slip1(:,:,k))

!Fracture energy along the fault
      EGfault(:,:,k)=0.5*peak_xz(:,:,k)*(Dc(:,:,k)-max(0.,Dc(:,:,k)-slip1(:,:,k))*(Dc(:,:,k)-slip1(:,:,k))/Dc(:,:,k))
      write(*,*)Eg(k),sum(EGfault(:,:,k))*dh*dh
      
!calculate nucleation center
      x0(k)=0.
      z0(k)=0.
      ncent=0
      !find  center
      do j=1,nzt
        z=dh*(real(j)-0.5)
        do i=1,nxt
          x=dh*(real(i)-0.5)
          if (peak_xz(i,j,k)<=strinix(i,j,k)) then
            x0(k)=x0(k)+x
            z0(k)=z0(k)+z
            ncent=ncent+1
          endif
        enddo
      enddo
      x0(k)=x0(k)/ncent
      z0(k)=z0(k)/ncent

      strengthexcess1(:,:,k)=strinix(:,:,k)-peak_xz(:,:,k)
!      nuclsize(k)=dh*dh*COUNT(strengthexcess1(:,:,k)>=1.e5)/1.e6
!      meanoverstress(k)=sum(strengthexcess1(:,:,k),strengthexcess1(:,:,k)>=1.e5)*dh*dh/1.e6/nuclsize(k)/1.e6
      nuclsize(k)=dh*dh*COUNT(strengthexcess1(:,:,k)>=0.)/1.e6
      meanoverstress(k)=sum(strengthexcess1(:,:,k),strengthexcess1(:,:,k)>=0.)*dh*dh/1.e6/nuclsize(k)/1.e6
      meanstrength(k)=sum(peak_xz(:,:,k)*slip1(:,:,k))/sum(slip1(:,:,k))

    if(meanoverstress(k)<2.)write(231,'(10000E13.5)')misfits(k),VRs(k),T0A(:,:,k),TsA(:,:,k),DcA(:,:,k)
    if(nuclsize(k)<2.)write(232,'(10000E13.5)')misfits(k),VRs(k),T0A(:,:,k),TsA(:,:,k),DcA(:,:,k)

!dependence of slip on Dc
      if(mod(k,10)==0)then
        do j=1,nzt,10
          do i=1,nxt,10
            if(slip1(i,j,k)>slipfact*slipmax)write(428,*)slip1(i,j,k),Dc(i,j,k),strengthexcess1(i,j,k)
          enddo
        enddo
        write(428,*);write(428,*)
      endif
      
!Rupture velocity
      do j=2,nzt-1
        do i=2,nxt-1
          if ((ruptime1(i,j,k).ne.1.e4).and.(ruptime1(i+1,j,k).ne.1.e4).and.(ruptime1(i,j+1,k).ne.1.e4) &
          .and.(ruptime1(i-1,j,k).ne.1.e4).and.(ruptime1(i,j-1,k).ne.1.e4)) then
            dum = (sqrt((ruptime1(i+1,j,k)-ruptime1(i,j,k))**2+(ruptime1(i,j+1,k)-ruptime1(i,j,k))**2) &
              +sqrt((ruptime1(i,j,k)-ruptime1(i-1,j,k))**2+(ruptime1(i,j,k)-ruptime1(i,j-1,k))**2))
            if (dum.ne.0.) then
              rupvel1(i,j,k)=2*dh/dum
            else
              rupvel1(i,j,k) = 0.
            endif
          endif
        enddo
      enddo
      !meanruptvel(k)=1./(sum(dum21(2:nxt-1,2:nzt-1)*slip1(2:nxt-1,2:nzt-1,k),ruptime1(2:nxt-1,2:nzt-1,k)>1.)/sum(slip1(2:nxt-1,2:nzt-1,k),ruptime1(2:nxt-1,2:nzt-1,k)>1.))/1.e3
      meanruptvel(k)=sum(rupvel1(2:nxt-1,2:nzt-1,k)*slip1(2:nxt-1,2:nzt-1,k))/sum(slip1(2:nxt-1,2:nzt-1,k))/1.e3
      
!      EG(k)=sum(peak_xz(:,:)*(Dc(:,:)-(Dc(:,:)-slip1(:,:,k))/Dc(:,:)*max(0.,Dc(:,:)-slip1(:,:,k))))/2.*dh*dh   ! Dissipated breakdown work (protoze je tam ten min)
!      ER(k)=sum((strinix(:,:)+peak_xz(:,:)*max(0.,1.-slip1(:,:,k)/Dc(:,:)))*slip1(:,:,k)-peak_xz(:,:)*(Dc(:,:)-(max(0.,Dc(:,:)-slip1(:,:,k)))**2/Dc(:,:)))/2.*dh*dh
      RE(k)=ER(k)/(ER(k)+EG(k)) !Radiation efficiency

      meansd(k)=-sum(schange1(:,:,k)*slip1(:,:,k))/sum(slip1(:,:,k))/1.e6

!      EGrate(k)=sum((peak_xz(1:nxt,1:nzt-2)*(Dc(1:nxt,1:nzt-2)-(Dc(1:nxt,1:nzt-2)-slip1(1:nxt,1:nzt-2,k))/Dc(1:nxt,1:nzt-2)*max(0.,Dc(1:nxt,1:nzt-2)-slip1(1:nxt,1:nzt-2,k))))*slip1(1:nxt,1:nzt-2,k))/2./sum(slip1(1:nxt,1:nzt-2,k))   ! Dissipated breakdown work rate
!      EGrate(k)=EG(k)/sum(slip1(:,:,k))*real(nxt*nzt)

      meansl(k)=sum(strinix(:,:,k)/peak_xz(:,:,k)*slip1(:,:,k))/sum(slip1(:,:,k)) !Stress level (Eq. 4 in Ripperger et al., 2007)
      
      duration(k)=maxval(ruptime1(:,:,k),slip1(:,:,k)>slipfact*slipmax)-minval(ruptime1(:,:,k),slip1(:,:,k)>slipfact*slipmax)
      
!      M0(k)=sum(slip1(:,:,k)*mu1(:,nyt,:))*dh*dh

      meanDc(k)=sum(Dc(:,:,k)*slip1(:,:,k))/sum(slip1(:,:,k))
      
!      meanslip(k)=sum(slip1(:,:,k)*slip1(:,:,k))/sum(slip1(:,:,k))
      meanslip(k)=sum(slip1(:,:,k))/count(slip1(:,:,k)>0.001)

      rupturearea(k)=count(slip1(:,:,k)>slipfact*slipmax)*dh*dh
      
      EGrate(k)=EG(k)/rupturearea(k)
      
      meanStrengthExcess(k)=-sum(strengthexcess1(:,:,k)*slip1(:,:,k),strengthexcess1(:,:,k)<1.e5)/sum(slip1(:,:,k),strengthexcess1(:,:,k)<1.e5)/1.e6

      cmrate=0.
      cmrate(1:nSr)=MomentRate(1:nSr,k)
      call four1(cmrate,np,-1)
      cmrate=cmrate*dtseis
      call findfc(cmrate(:),m0(k),np,dtseis,fc(k))
      do j=1,np/2+1
        cmrate(j)=-cmrate(j)*(2.*pi*dfseis*(j-1))**2
      enddo
      CALL smoothspectrum(np,nper,dfseis,fc1,fc2,abs(cmrate),per(:),MomentSpec(:,k))
      maxMR=maxval(MomentRate(1:nSr,k))
      i=0
      ierr=0
      do while (ierr==0)
        i=i+1
        if (Momentrate(i,k)>0.1*maxMR) ierr=1
      enddo
      ierr=0
      j=nSr+1
      do while (ierr==0)
        j=j-1
        if (momentrate(j,k)>0.1*maxMR) ierr=1
      enddo
      dur1=(j-i)*dtseis
      dur2=2*m0(k)/maxMR
    enddo

    CALL indexx(NM,exp(-(misfits(:)-bestmisfit)),indx)
    open(201,FILE='processBayes.dat')
    do j=NM,1,-1
      k=indx(j)
!                               1         2       3         4            5        6    7      8       9             10          11      12             13               14        15     16      17             18            19            20       21
      write(201,'(100E13.5)')misfits(k),VRs(k),meansd(k),duration(k),nuclsize(k),EG(k),ER(k),RE(k),meansl(k),meanoverstress(k),M0(k),meanDc(k),meanStrengthExcess(k),meanslip(k),x0(k),z0(k),rupturearea(k),meanruptvel(k),meanstrength(k),Egrate(k),fc(k)
    enddo
    close(201)

    write(*,*)'Scalar moment (x 1e17): ',sum(M0(:)/1.e17)/NM,'+-',sqrt(sum((M0(:)/1.e17)**2)/NM-(sum(M0(:)/1.e17)/NM)**2)
    write(*,*)'Mean stress drop: ',sum(meansd(:))/NM,'+-',sqrt(sum(meansd(:)**2)/NM-(sum(meansd(:))/NM)**2)
    write(*,*)'Fracture energy: ',sum(EG(:))/NM,'+-',sqrt(sum(EG(:)**2)/NM-(sum(EG(:))/NM)**2)
    write(*,*)'Fracture energy rate: ',sum(EGrate(:))/NM,'+-',sqrt(sum(EGrate(:)**2)/NM-(sum(EGrate(:))/NM)**2)
    write(*,*)'Radiated energy: ',sum(ER(:))/NM,'+-',sqrt(sum(ER(:)**2)/NM-(sum(ER(:))/NM)**2)
    write(*,*)'Mean Dc: ',sum(meanDc(:))/NM,'+-',sqrt(sum(meanDc(:)**2)/NM-(sum(meanDc(:))/NM)**2)
    write(*,*)'Radiation efficiency: ',sum(ER(:)/(ER(:)+EG(:)))/NM,'+-',sqrt(sum((ER(:)/(ER(:)+EG(:)))**2)/NM-(sum(ER(:)/(ER(:)+EG(:)))/NM)**2)
    write(*,*)'Scaled energy: ',sum(ER(:)/M0(:))/NM,'+-',sqrt(sum((ER(:)/M0(:))**2)/NM-(sum(ER(:)/M0(:))/NM)**2)
    write(*,*)'Mean slip: ',sum(meanslip(:))/NM,'+-',sqrt(sum(meanslip(:)**2)/NM-(sum(meanslip(:))/NM)**2)
    
    open(201,FILE='processBayes.slipmodels.dat')
    open(202,FILE='processBayes.strengthexcess.dat')
    do j=1,nzt
      do i=1,nxt
        write(201,'(10000E13.5)')(slip1(i,j,indx(k)),k=NM,1,-1)
        write(202,'(10000E13.5)')(strengthexcess1(i,j,indx(k))/1.e6,k=NM,1,-1)
      enddo
    enddo
    close(201)
    close(202)

    open(201,FILE='processBayes.MomentRates.dat')
    open(202,FILE='processBayes.MomentSpectra.dat')
    do k=1,NM
      do j=1,nSR
        write(201,'(10000E13.5)')dtseis*(j-1),MomentRate(j,indx(k)),exp(-(misfits(indx(k))-bestmisfit))
      enddo
      do j=1,nper
        write(202,'(10000E13.5)')per(j),MomentSpec(j,indx(k)),exp(-(misfits(indx(k))-bestmisfit))
      enddo
      write(201,*);write(201,*)
      write(202,*);write(202,*)
    enddo
    close(201)
    close(202)

    allocate(mask(nxt,nzt))
    mask=1.
    !Optional: Apply mask from file mask.dat, which can be obtained by processBayes.slipmodels.sh:
    !open(201,FILE='mask.dat')
    !do j=1,nzt
    !  read(201,*)(mask(i,j),i=1,nxt)
    !enddo
    !close(201)
    
    open(201,FILE='processBayes.meansigma.dat')
    CALL meansigma(T0A(:,:,:),NM)   !Removed conversion to MPa, so the plotting scripts must be modified accordingly
    CALL meansigma(TSA(:,:,:),NM)
    CALL meansigma(SEA(:,:,:),NM)   !Removed conversion to MPa, so the plotting scripts must be modified accordingly
    CALL meansigma(DcA(:,:,:),NM)
    CALL meansigma2(es1(:,:,:),NM)
!    CALL meansigma(SEA(:,:,:)/max(1.,T0A(:,:,:)),NM)  !Dava jen uzky pruh, hodne ta hodnota asi lita.
    close(201)

    open(201,FILE='processBayes.meansigma2.dat')
    CALL meansigma2(slip1(:,:,:),NM)
    CALL meansigma2(schange1(:,:,:),NM)
    CALL meansigma2(rise1(:,:,:),NM)
    CALL meansigma2(ruptime1(:,:,:),NM)
    CALL meansigma2(rupvel1(:,:,:),NM)
    CALL meansigma2(Egfault(:,:,:),NM)
    close(201)

!FFT of dynamic parameters
Dc=log(Dc)   !POZOR!!
    nxtfft=2**int(log(dble(nxt))/log(2.d0)+1)
    nztfft=2**int(log(dble(nzt))/log(2.d0)+1)
    allocate(strinixFFTc(nxtfft/2,nztfft,NM),peak_xzFFTc(nxtfft/2,nztfft,NM),DcFFTc(nxtfft/2,nztfft,NM))
    allocate(dumFFT(nxtfft,nztfft),dumFFTq(nztfft))
    allocate(strinixSpatialmean(nxt,nzt),peak_xzSpatialmean(nxt,nzt),DcSpatialmean(nxt,nzt))
    strinixSpatialmean=0.d0;peak_xzSpatialmean=0.d0;DcSpatialmean=0.d0
    strinixpwr=0.;peak_xzpwr=0.;Dcpwr=0.
    strinixpeak_xzxcorr=0.;strinixDcxcorr=0.;peak_xzDcxcorr=0.
    do k=1,NM
      strinixSpatialmean(:,:)=strinixSpatialmean(:,:)+strinix(:,:,k)/1.e6/real(NM)
      peak_xzSpatialmean(:,:)=peak_xzSpatialmean(:,:)+peak_xz(:,:,k)/1.e6/real(NM)
      DcSpatialmean(:,:)=DcSpatialmean(:,:)+Dc(:,:,k)/real(NM)
    enddo
    do k=1,NM
      dumFFT=0.d0;dumFFT(1:nxt,1:nzt)=(strinix(1:nxt,1:nzt,k)/1.e6-strinixSpatialmean(1:nxt,1:nzt))*slip1(1:nxt,1:nzt,k)/sqrt(sum(slip1(:,:,k)**2)/dble(nxt*nzt))
      strinixpwr=strinixpwr+sum((strinix(1:nxt,1:nzt,k)/1.e6-strinixSpatialmean(1:nxt,1:nzt))**2*slip1(1:nxt,1:nzt,k)**2)/sum(slip1(1:nxt,1:nzt,k)**2)/real(NM)
!      if(k==100)then
!        do j=1,nzt
!          write(399,'(10000E13.5)')(strinix(i,j,k)/1.e6-strinixSpatialmean(i,j),i=1,nxt)
!        enddo
!      endif      
      CALL rlft3(dumFFT,dumFFTq,nxtfft,nztfft,1,1)
      do i=1,nxtfft/2
        strinixFFTc(i,:,k)=cmplx(dumFFT(2*i-1,:),dumFFT(2*i,:))*(dh*dh)/sqrt(dble(nxt*nzt)*dh*dh)
      enddo
!      dumFFT=0.d0;dumFFT(1:nxt,1:nzt)=(peak_xz(1:nxt,1:nzt,k)/1.e6-peak_xzSpatialmean(1:nxt,1:nzt))*slip1(1:nxt,1:nzt,k)/sqrt(sum(slip1(:,:,k)**2)/dble(nxt*nzt))  !peak stress
!      peak_xzpwr=peak_xzpwr+sum((peak_xz(1:nxt,1:nzt,k)/1.e6-peak_xzSpatialmean(1:nxt,1:nzt))**2*slip1(1:nxt,1:nzt,k)**2)/sum(slip1(1:nxt,1:nzt,k)**2)/real(NM)   !peak stress
      dumFFT=0.d0;dumFFT(1:nxt,1:nzt)=(peak_xz(1:nxt,1:nzt,k)/1.e6-peak_xzSpatialmean(1:nxt,1:nzt)-(strinix(1:nxt,1:nzt,k)/1.e6-strinixSpatialmean(1:nxt,1:nzt)))*slip1(1:nxt,1:nzt,k)/sqrt(sum(slip1(:,:,k)**2)/dble(nxt*nzt))   !strength excess
      peak_xzpwr=peak_xzpwr+sum((peak_xz(1:nxt,1:nzt,k)/1.e6-peak_xzSpatialmean(1:nxt,1:nzt)-(strinix(1:nxt,1:nzt,k)/1.e6-strinixSpatialmean(1:nxt,1:nzt)))**2*slip1(1:nxt,1:nzt,k)**2)/sum(slip1(1:nxt,1:nzt,k)**2)/real(NM)   !strength excess
      CALL rlft3(dumFFT,dumFFTq,nxtfft,nztfft,1,1)
      do i=1,nxtfft/2
        peak_xzFFTc(i,:,k)=cmplx(dumFFT(2*i-1,:),dumFFT(2*i,:))*(dh*dh)/sqrt(dble(nxt*nzt)*dh*dh)
      enddo
      dumFFT=0.d0;dumFFT(1:nxt,1:nzt)=(Dc(1:nxt,1:nzt,k)-DcSpatialmean(1:nxt,1:nzt))*slip1(1:nxt,1:nzt,k)/sqrt(sum(slip1(:,:,k)**2)/dble(nxt*nzt))
      Dcpwr=Dcpwr+sum((Dc(1:nxt,1:nzt,k)-DcSpatialmean(1:nxt,1:nzt))**2*slip1(1:nxt,1:nzt,k)**2)/sum(slip1(1:nxt,1:nzt,k)**2)/real(NM)
      CALL rlft3(dumFFT,dumFFTq,nxtfft,nztfft,1,1)
      do i=1,nxtfft/2
        DcFFTc(i,:,k)=cmplx(dumFFT(2*i-1,:),dumFFT(2*i,:))*(dh*dh)/sqrt(dble(nxt*nzt)*dh*dh)
      enddo
      !DcFFTq(:,k)=dumFFTq(:,k)*(dh*dh)/sqrt(dble(nxt*nzt)*dh*dh)   !neglecting Nyqist frequency in DcFFTq
      !strinixpeak_xzxcorr=strinixpeak_xzxcorr+sum((strinix(1:nxt,1:nzt,k)/1.e6-strinixSpatialmean(1:nxt,1:nzt))*(peak_xz(1:nxt,1:nzt,k)/1.e6-peak_xzSpatialmean(1:nxt,1:nzt))*slip1(1:nxt,1:nzt,k)**2)/sum(slip1(1:nxt,1:nzt,k)**2)/real(NM)   !peak stress
      strinixpeak_xzxcorr=strinixpeak_xzxcorr+sum((strinix(1:nxt,1:nzt,k)/1.e6-strinixSpatialmean(1:nxt,1:nzt))*(peak_xz(1:nxt,1:nzt,k)/1.e6-peak_xzSpatialmean(1:nxt,1:nzt)-(strinix(1:nxt,1:nzt,k)/1.e6-strinixSpatialmean(1:nxt,1:nzt)))*slip1(1:nxt,1:nzt,k)**2)/sum(slip1(1:nxt,1:nzt,k)**2)/real(NM)   !strength excess
      strinixDcxcorr=strinixDcxcorr+sum((strinix(1:nxt,1:nzt,k)/1.e6-strinixSpatialmean(1:nxt,1:nzt))*(Dc(1:nxt,1:nzt,k)-DcSpatialmean(1:nxt,1:nzt))*slip1(1:nxt,1:nzt,k)**2)/sum(slip1(1:nxt,1:nzt,k)**2)/real(NM)
      peak_xzDcxcorr=peak_xzDcxcorr+sum((peak_xz(1:nxt,1:nzt,k)/1.e6-peak_xzSpatialmean(1:nxt,1:nzt))*(Dc(1:nxt,1:nzt,k)-DcSpatialmean(1:nxt,1:nzt))*slip1(1:nxt,1:nzt,k)**2)/sum(slip1(1:nxt,1:nzt,k)**2)/real(NM)
    enddo
    deallocate(dumFFT,dumFFTq)
    write(*,*)'strinixpeak_xzxcorr,strinixDcxcorr,peak_xzDcxcorr: ',strinixpeak_xzxcorr/sqrt(strinixpwr)/sqrt(peak_xzpwr),strinixDcxcorr/sqrt(strinixpwr)/sqrt(Dcpwr),peak_xzDcxcorr/sqrt(peak_xzpwr)/sqrt(Dcpwr)
    write(*,*)'strinixstd,peak_xzstd,Dcstd: ',sqrt(strinixpwr),sqrt(peak_xzpwr),sqrt(Dcpwr)
    open(438,FILE='paramspectra.dat')
    strinixpwr=0.;peak_xzpwr=0.;Dcpwr=0.
    do j=1,nztfft
      if(j<=nztfft/2+1)then
        ky=1./dble(nztfft)/dh*real(j-1)
      else
        ky=-1./dble(nztfft)/dh*real(nztfft-j+1)
      endif
      do i=1,nxtfft/2    !neglecting Nyqist frequency in DcFFTq, otherwise +1
        kx=1./dble(nxtfft)/dh*real(i-1)
        krad=sqrt(kx**2+ky**2)
        !write(438,*)krad,abs(DcFFTc(i,j,1))
        write(438,'(4E13.5)')krad,sum(abs(strinixFFTc(i,j,:))**2)/real(NM),sum(abs(peak_xzFFTc(i,j,:))**2)/real(NM),sum(abs(DcFFTc(i,j,:))**2)/real(NM)
        strinixpwr=strinixpwr+2.*sum(abs(strinixFFTc(i,j,:))**2)/real(NM)/dble(nxtfft*nztfft)/dh**2
        peak_xzpwr=peak_xzpwr+2.*sum(abs(peak_xzFFTc(i,j,:))**2)/real(NM)/dble(nxtfft*nztfft)/dh**2
        Dcpwr=Dcpwr+2.*sum(abs(DcFFTc(i,j,:))**2)/real(NM)/dble(nxtfft*nztfft)/dh**2
      enddo
    enddo
    write(*,*)sqrt(strinixpwr),sqrt(peak_xzpwr),sqrt(Dcpwr)
    close(438)
    
    END
    
    
    SUBROUTINE meansigma(arr,NM)
    USE fd3dparam_pb
    USE interp_pb
    USE mask_pb
    IMPLICIT NONE
    INTEGER NM
    real arr(NLI,NWI,NM)
    REAL,ALLOCATABLE,DIMENSION(:,:):: meanA,sigmaA,mean,sigma
    INTEGER i,j
    
    allocate(meanA(NLI,NWI),sigmaA(NLI,NWI),mean(nxt,nzt),sigma(nxt,nzt))
    do j=1,NWI
      do i=1,NLI
        meanA(i,j)=sum(arr(i,j,1:NM))/real(NM)
        sigmaA(i,j)=sqrt(sum(arr(i,j,1:NM)**2)/real(NM)-meanA(i,j)**2)
!        sigmaA(i,j)=(maxval(arr(i,j,:))-minval(arr(i,j,:)))    !Histogram maximum width
      enddo
    enddo
    CALL interpolate(meanA(:,:),mean(:,:))
    CALL interpolate(sigmaA(:,:),sigma(:,:)) 
    do j=1,nzt
      write(201,'(10000E13.5)')(mean(i,j)*mask(i,j),i=1,nxt)
    enddo
    write(201,*);write(201,*)
    do j=1,nzt
      write(201,'(10000E13.5)')(sigma(i,j)*2.*mask(i,j),i=1,nxt)
!      write(201,'(10000E13.5)')(sigma(i,j)/mean(i,j),i=1,nxt)  !Relative sigma
    enddo
    write(201,*);write(201,*)
    deallocate(meanA,sigmaA,mean,sigma)
    
    END
    
    
    SUBROUTINE meansigma2(arr,NM)
    USE fd3dparam_pb
    IMPLICIT NONE
    INTEGER NM
    real arr(nxt,nzt,NM)
    REAL,ALLOCATABLE,DIMENSION(:,:):: mean,sigma
    INTEGER i,j
    
    allocate(mean(nxt,nzt),sigma(nxt,nzt))
    do j=1,nzt
      do i=1,nxt
        mean(i,j)=sum(arr(i,j,1:NM))/real(NM)
        sigma(i,j)=sqrt(sum(arr(i,j,1:NM)**2)/real(NM)-mean(i,j)**2)
!        sigma(i,j)=(maxval(arr(i,j,:))-minval(arr(i,j,:)))    !Histogram maximum width
      enddo
    enddo
    do j=1,nzt
      write(201,'(10000E13.5)')(mean(i,j),i=1,nxt)
    enddo
    write(201,*);write(201,*)
    do j=1,nzt
      write(201,'(10000E13.5)')(sigma(i,j)*2.,i=1,nxt)
!      write(201,'(10000E13.5)')(sigma(i,j)/mean(i,j),i=1,nxt)  !Relative sigma
    enddo
    write(201,*);write(201,*)
    deallocate(mean,sigma)
    
    END
    
subroutine findfc(mrspectr,m0,np,dtseis,fc)
implicit none
complex :: mrspectr(1:np)
real :: m0,dtseis
real :: fc
integer :: n1,np,ipoint,npoints,i
integer :: fmin,fmax
real :: fprop,df
real*8 :: misf,bestmisf
real :: fi(1:np),teor(1:np)
df=1./(dtseis*np)
fmin=max(1,floor(0.02/df))
fmax=min(np/2,ceiling(5./df))
n1=np/2+1
do i=1,n1
  fi(i)=(i-1)*df
enddo
bestmisf=1.e30
do ipoint=fmin,fmax,2
  fprop=df*(ipoint-1)
  teor(:)=m0/(1.+(fi(:)/fprop)**2.)
  misf=0.5*sum((log(abs(mrspectr(fmin:fmax))/teor(fmin:fmax)))**2)  
  if (misf<bestmisf) then
    bestmisf=misf
    fc=fprop
  endif
enddo
end subroutine    


      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
    END

    
    FUNCTION M0toMw(m0)
    implicit none
    real :: m0,m0tomw
    m0tomw=(log10(m0)-9.1)/1.5
    END FUNCTION 
    

    SUBROUTINE smoothspectrum(Nf,Nfsmooth,df,flo,fro,spec,freqaxis,smoothspec)
    !Smoothing spectrum by Konno & Omachi, 1998 BSSA, method
    IMPLICIT NONE
    REAL,PARAMETER:: b=40.
    INTEGER Nf,Nfsmooth
    REAL spec(Nf),WB(Nf/2+1,Nfsmooth),freqaxis(Nfsmooth),smoothspec(Nfsmooth),flo,fro
    REAL freq,df
    INTEGER i,j
    do j=1,Nfsmooth
      freqaxis(j)=10.**((log10(fro)-log10(flo))/real(Nfsmooth-1)*real(j-1)+log10(flo))
      WB(:,j)=0.
      do i=2,Nf/2+1
        freq=df*(i-1)
        if(freq.ne.freqaxis(j))then
          WB(i,j)=(sin(b*log10(freq/freqaxis(j)))/b/log10(freq/freqaxis(j)))**4
        else
          WB(i,j)=1.
        endif
      enddo
    enddo
    do j=1,Nfsmooth
      smoothspec(j)=sum(abs(spec(1:Nf/2+1))*WB(1:Nf/2+1,j))/sum(WB(1:Nf/2+1,j))
    enddo        
    END
    
    
    SUBROUTINE interpolate(arrin,arrout)     ! Bilinear interpolation
    USE interp_pb
    USE fd3dparam_pb
    IMPLICIT NONE
    REAL arrin(NLI,NWI),arrout(nxt,nzt)
    REAL DL,DW,ZS,XS,t,u
    INTEGER i,k,ii,kk

    DL=dh*(nxt-1)/real(NLI-1)
    DW=dh*(nzt-1)/real(NWI-1)
    do k=1,nzt
      ZS=dh*(k-1)
      kk=min(NWI-1,int(ZS/DW)+1)
      u=min(1.,(ZS-DW*(kk-1))/DW)
      do i=1,nxt
        XS=dh*(i-1)
        ii=min(NLI-1,int(XS/DL)+1)
        t=min(1.,(XS-DL*(ii-1))/DL)
        arrout(i,k)=(1.-t)*(1.-u)*arrin(ii,kk)+t*(1.-u)*arrin(ii+1,kk)+t*u*arrin(ii+1,kk+1)+(1.-t)*u*arrin(ii,kk+1)
      enddo
    enddo

    END
    

   SUBROUTINE readcrustalmodel(dip)
    USE medium_pb
    USE fd3dparam_pb
    IMPLICIT NONE
    real*8,parameter:: PI=3.1415926535
    real dip
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
    do k=nzt,1,-1
      dum=(dh*real(nzt-k)+dh/2.)*sin(dip/180.d0*PI)
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
      lam1(:,:,k) = dd*(vpp**2-2.*vss**2)
      mu1(:,:,k)  = dd*vss**2
      d1(:,:,k)   = dd
    enddo
!    close(10)
    mu_mean = (mu_mean/nzt)
!    write(*,*)mu_mean
    deallocate(depth,vp,vs,rho)

    END


      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
    END


      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
      INTEGER isign,nn1,nn2,nn3
      COMPLEX*16 data(nn1/2,nn2,nn3),speq(nn2,nn3)
      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      COMPLEX*16 c1,c2,h1,h2,w
      c1=dcmplx(0.5d0,0.0d0)
      c2=dcmplx(0.0d0,-0.5d0*isign)
      theta=6.28318530717959d0/dble(isign*nn1)
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      nn(1)=nn1/2
      nn(2)=nn2
      nn(3)=nn3
      if(isign.eq.1)then
        call fourn(data,nn,3,isign)
        do 12 i3=1,nn3
          do 11 i2=1,nn2
            speq(i2,i3)=data(1,i2,i3)
11        continue
12      continue
      endif
      do 15 i3=1,nn3
        j3=1
        if (i3.ne.1) j3=nn3-i3+2
        wr=1.0d0
        wi=0.0d0
        do 14 i1=1,nn1/4+1
          j1=nn1/2-i1+2
          do 13 i2=1,nn2
            j2=1
            if (i2.ne.1) j2=nn2-i2+2
            if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
            else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
            endif
13        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          w=dcmplx(dble(wr),dble(wi))
14      continue
15    continue
      if(isign.eq.-1)then
        call fourn(data,nn,3,isign)
      endif
      return
      END

    
      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      DOUBLE PRECISION data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=dble(wr)*data(k2)-dble(wi)*data(k2+1)
                tempi=dble(wr)*data(k2+1)+dble(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END
