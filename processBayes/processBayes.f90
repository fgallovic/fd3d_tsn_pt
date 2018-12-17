    MODULE interp_pb
      INTEGER NLI,NWI
    END MODULE

    MODULE fd3dparam_pb
      integer :: nxt,nyt,nzt,ntfd,nysc
      real    :: dh,dt
    END MODULE

    MODULE friction_pb
      USE fd3dparam_pb
      real,allocatable,dimension(:,:):: strinix,peak_xz,Dc
      real,parameter:: pi=3.1415926535
      REAL dip
      CONTAINS
      FUNCTION normstress(j)
      IMPLICIT NONE
      real:: normstress
      integer:: j
#if defined DIPSLIP
      normstress=max(1.e5,8520.*dh*real(nzt-j)*sin(dip/180.*pi))
#else
      normstress=max(1.e5,16200.*dh*real(nzt-j)*sin(dip/180.*pi))
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
    IMPLICIT NONE
    INTEGER,PARAMETER:: NMAX=1e6
    REAL,ALLOCATABLE,DIMENSION(:):: normalstress,VRs,misfits,meansd,meansl,duration,nuclsize,EG,ER,RE,meanoverstress,M0,meanDc,meanStrengthExcess,meanslip,rupturearea,meanruptvel,meanstrength
    REAL,ALLOCATABLE,DIMENSION(:,:,:):: DcA,TsA,T0A,SEA,ruptime1,slip1,rise1,schange1,es1,strengthexcess1
    REAL,ALLOCATABLE,DIMENSION(:,:):: dum11,dum12,dum13,dum21,dum22,dum23,dum24,ms1

    REAL bestmisfit,misfitaccept,dum
    REAL vr,mf,x0,z0,x,z,slipmax
    INTEGER NM,NTOT
    INTEGER i,j,k,ml(1),ncent
    
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

    allocate(lam1(nxt,nyt,nzt),mu1(nxt,nyt,nzt),d1(nxt,nyt,nzt))
    allocate(strinix(nxt,nzt),peak_xz(nxt,nzt),Dc(nxt,nzt))

    CALL readcrustalmodel(dip)

    open(10,FILE='inputinv.dat')
    read(10,*)
    read(10,*)NLI,NWI
    
    ALLOCATE(misfits(NMAX),VRs(NMAX))
    ALLOCATE(dum11(NLI,NWI),dum12(NLI,NWI),dum13(NLI,NWI),dum21(nxt,nzt),dum22(nxt,nzt),dum23(nxt,nzt),dum24(nxt,nzt))
    open(719,FILE='normalstressprofile.dat')
    ALLOCATE(normalstress(nzt))
    do i=1,nzt
      read(719,*)dum,normalstress(i)
    enddo
    normalstress(:)=normalstress(:)*1.e6
    close(719)

!------ Learn about misfits
    k=0
    open(101,FILE='sampls.dat',FORM='UNFORMATTED',ACCESS='STREAM')
10  read(101,END=11,ERR=11)mf,vr,dum11(:,:),dum12(:,:),dum13(:,:),dum21(:,:),dum22(:,:),dum23(:,:),dum24(:,:)
    k=k+1
    misfits(k)=mf
    if(k==NMAX)stop 'Increase dimension!'
    goto 10
11  close(101)
    NTOT=k
    write(*,*)'Total number of models: ',NTOT
    bestmisfit=minval(misfits(1:NTOT))
    write(*,*)'Best misfit: ',bestmisfit
    
!    misfitaccept=maxval(misfits(1:ntot))
!    misfitaccept=bestmisfit-log(0.02) !Probability threashold (2% for real data)
    misfitaccept=bestmisfit-log(0.001) !Probability threashold (1%% for inv1)
!   misfitaccept=bestmisfit-log(0.0001) !Probability threashold (.1%% for pga)
    
    NM=0
    do i=1,k
      if(misfits(i)<=misfitaccept)NM=NM+1
    enddo
    write(*,*)'Number of accepted models: ',NM
    DEALLOCATE(misfits,VRs)

!------ Read accepted models
    ALLOCATE(misfits(NM),VRs(NM),meansd(NM),meansl(NM),duration(NM),nuclsize(NM),EG(NM),ER(NM),RE(NM),meanoverstress(NM),M0(NM))
    ALLOCATE(meanDc(NM),meanStrengthExcess(NM),meanslip(NM),rupturearea(NM),meanruptvel(NM),meanstrength(NM))
    allocate(DcA(NLI,NWI,NM),TsA(NLI,NWI,NM),T0A(NLI,NWI,NM),SEA(NLI,NWI,NM))
    allocate(ruptime1(nxt,nzt,NM),slip1(nxt,nzt,NM),rise1(nxt,nzt,NM),schange1(nxt,nzt,NM),es1(nxt,nzt,NM),ms1(nxt,nzt),strengthexcess1(nxt,nzt,NM))
    open(101,FILE='sampls.dat',FORM='UNFORMATTED',ACCESS='STREAM')
    k=0
    do i=1,NTOT
      read(101)mf,vr,dum11(:,:),dum12(:,:),dum13(:,:),dum21(:,:),dum22(:,:),dum23(:,:),dum24(:,:)
      if(mf<=misfitaccept)then
        k=k+1
        misfits(k)=mf
        VRs(k)=vr
        T0A(:,:,k)=dum11(:,:)
        TsA(:,:,k)=dum12(:,:)
        DcA(:,:,k)=dum13(:,:)
        
        do j=1,NWI
          SEA(:,j,k)=TsA(:,j,k)*normalstress((j-1)*((nzt)/(NWI-1))+1)-T0A(:,j,k)
        enddo
        
        ruptime1(:,:,k)=dum21(:,:)
        slip1(:,:,k)=dum22(:,:)
        rise1(:,:,k)=dum23(:,:)
        schange1(:,:,k)=dum24(:,:)
      endif
    enddo
    close(101)

! Save best and worst accepted models
    open(201,FILE='forwardmodel.best.dat')
    ml(:)=minloc(misfits(:))
    write(201,'(10000E13.5)')misfits(ml(1)),VRs(ml(1)),T0A(:,:,ml(1)),TsA(:,:,ml(1)),DcA(:,:,ml(1))
    close(201)
    open(201,FILE='forwardmodel.worst.dat')
    ml(:)=maxloc(misfits(:))
    write(201,'(10000E13.5)')misfits(ml(1)),VRs(ml(1)),T0A(:,:,ml(1)),TsA(:,:,ml(1)),DcA(:,:,ml(1))
    close(201)
    
    open(201,FILE='processBayes.dat')
    do k=1,NM
      strinix=0.
      ms1=0.
      peak_xz=0.
      CALL interpolate(T0A(:,:,k),strinix(:,:))
      CALL interpolate(TSA(:,:,k),ms1(:,:))
      CALL interpolate(DcA(:,:,k),Dc(:,:))
      do i=1,nzt
        peak_xz(:,i)=ms1(:,i)*normalstress(i)
      enddo
      es1(:,:,k)=(peak_xz(:,:)-strinix(:,:))/max(1.,strinix(:,:))
      slipmax=maxval(slip1(:,:,k))

!calculate nucleation center
      x0=0.
      z0=0.
      ncent=0
      !find  center
      do j=1,nzt
        z=dh*(real(j)-0.5)
        do i=1,nxt
          x=dh*(real(i)-0.5)
          if (peak_xz(i,j)<=strinix(i,j)) then
            x0=x0+x
            z0=z0+z
            ncent=ncent+1
          endif
        enddo
      enddo
      x0=x0/ncent
      z0=z0/ncent

      strengthexcess1(:,:,k)=strinix(:,:)-peak_xz(:,:)
      nuclsize(k)=dh*dh*COUNT(strengthexcess1(:,:,k)>=1.e5)/1.e6
      meanoverstress(k)=sum(strengthexcess1(:,:,k),strengthexcess1(:,:,k)>=1.e5)*dh*dh/1.e6/nuclsize(k)/1.e6
      meanstrength(k)=sum(peak_xz(:,:)*slip1(:,:,k))/sum(slip1(:,:,k))

!dependence of slip on Dc
      if(mod(k,10)==0)then
        do j=1,nzt,10
          do i=1,nxt,10
            if(slip1(i,j,k)>0.05*slipmax)write(428,*)slip1(i,j,k),Dc(i,j),strengthexcess1(i,j,k)
          enddo
        enddo
        write(428,*);write(428,*)
      endif
      
!Rupture velocity
      do j=2,nzt-3
        do i=2,nxt-1
          dum21(i,j)=sqrt((ruptime1(i+1,j,k)-ruptime1(i-1,j,k))**2+(ruptime1(i,j+1,k)-ruptime1(i,j-1,k))**2)/2./dh  !slowness
        enddo
      enddo
      meanruptvel(k)=1./(sum(dum21(2:nxt-1,2:nzt-3)*slip1(2:nxt-1,2:nzt-3,k),ruptime1(2:nxt-1,2:nzt-3,k)>1.)/sum(slip1(2:nxt-1,2:nzt-3,k),ruptime1(2:nxt-1,2:nzt-3,k)>1.))/1.e3
      
      EG(k)=sum(peak_xz(:,:)*(Dc(:,:)-(Dc(:,:)-slip1(:,:,k))/Dc(:,:)*max(0.,Dc(:,:)-slip1(:,:,k))))/2.*dh*dh   ! Dissipated breakdown work (protoze je tam ten min)
      ER(k)=sum(strinix(:,:)*slip1(:,:,k)-peak_xz(:,:)*(Dc(:,:)-(max(0.,Dc(:,:)-slip1(:,:,k)))**2/Dc(:,:)))/2.*dh*dh
      RE(k)=ER(k)/(ER(k)+EG(k)) !Radiation efficiency
            
      meansd(k)=-sum(schange1(:,:,k)*slip1(:,:,k))/sum(slip1(:,:,k))/1.e6
     
      meansl(k)=sum(strinix(:,:)/peak_xz(:,:)*slip1(:,:,k))/sum(slip1(:,:,k)) !Stress level (Eq. 4 in Ripperger et al., 2007)
      
      duration(k)=maxval(ruptime1(:,:,k),slip1(:,:,k)>0.05*slipmax)-minval(ruptime1(:,:,k),slip1(:,:,k)>0.05*slipmax)
      
      M0(k)=sum(slip1(:,:,k)*mu1(:,nyt,:))*dh*dh

      meanDc(k)=sum(Dc(:,:)*slip1(:,:,k))/sum(slip1(:,:,k))
      
      meanslip(k)=sum(slip1(:,:,k)*slip1(:,:,k))/sum(slip1(:,:,k))
      
      rupturearea(k)=count(slip1(:,:,k)>0.05*slipmax)*dh*dh
      
      meanStrengthExcess(k)=-sum(strengthexcess1(:,:,k)*slip1(:,:,k),strengthexcess1(:,:,k)<1.e5)/sum(slip1(:,:,k),strengthexcess1(:,:,k)<1.e5)/1.e6

!                               1         2       3         4            5        6    7      8       9             10          11      12             13               14       15 16      17             18            19
      write(201,'(100E13.5)')misfits(k),VRs(k),meansd(k),duration(k),nuclsize(k),EG(k),ER(k),RE(k),meansl(k),meanoverstress(k),M0(k),meanDc(k),meanStrengthExcess(k),meanslip(k),x0,z0,rupturearea(k),meanruptvel(k),meanstrength(k)
    enddo
    close(201)

    write(*,*)'Mean stress drop: ',sum(meansd(:))/NM,'+-',sqrt(sum(meansd(:)**2)/NM-(sum(meansd(:))/NM)**2)
    write(*,*)'Fracture energy: ',sum(EG(:))/NM,'+-',sqrt(sum(EG(:)**2)/NM-(sum(EG(:))/NM)**2)
    write(*,*)'Radiated energy: ',sum(ER(:))/NM,'+-',sqrt(sum(ER(:)**2)/NM-(sum(ER(:))/NM)**2)
    write(*,*)'Mean Dc: ',sum(meanDc(:))/NM,'+-',sqrt(sum(meanDc(:)**2)/NM-(sum(meanDc(:))/NM)**2)
    write(*,*)'Mean slip: ',sum(meanslip(:))/NM,'+-',sqrt(sum(meanslip(:)**2)/NM-(sum(meanslip(:))/NM)**2)
    
    
    open(201,FILE='processBayes.slipmodels.dat')
    open(202,FILE='processBayes.strengthexcess.dat')
    do j=1,nzt
      do i=1,nxt
        write(201,'(10000E13.5)')slip1(i,j,1:NM)
        write(202,'(10000E13.5)')strengthexcess1(i,j,1:NM)/1.e6
      enddo
    enddo
    close(201)
    close(202)

    open(201,FILE='smaz.txt')
    do j=1,nzt
      write(201,'(10000E13.5)')(strengthexcess1(i,j,10),i=1,nxt)
    enddo
    close(201)
    
    open(201,FILE='processBayes.meansigma.dat')
    CALL meansigma(T0A(:,:,:)/1.e6,NM)
    CALL meansigma(TSA(:,:,:),NM)
    CALL meansigma(SEA(:,:,:)/1.e6,NM)
    CALL meansigma(DcA(:,:,:),NM)
    CALL meansigma2(es1(:,:,:),NM)
!    CALL meansigma(SEA(:,:,:)/max(1.,T0A(:,:,:)),NM)  !Dava jen uzky pruh, hodne ta hodnota asi lita.
    close(201)
    
    open(201,FILE='processBayes.meansigma2.dat')
    CALL meansigma2(slip1(:,:,:),NM)
    CALL meansigma2(schange1(:,:,:),NM)
    CALL meansigma2(rise1(:,:,:),NM)
    CALL meansigma2(ruptime1(:,:,:),NM)
    close(201)
    
    END
    
    
    SUBROUTINE meansigma(arr,NM)
    USE fd3dparam_pb
    USE interp_pb
    IMPLICIT NONE
    INTEGER NM
    real arr(NLI,NWI,NM)
    REAL,ALLOCATABLE,DIMENSION(:,:):: meanA,sigmaA,mean,sigma
    INTEGER i,j
    
    allocate(meanA(NLI,NWI),sigmaA(NLI,NWI),mean(nxt,nzt),sigma(nxt,nzt))
    do j=1,NWI
      do i=1,NLI
        meanA(i,j)=sum(arr(i,j,1:NM))/real(NM)
        sigmaA(i,j)=sqrt(sum(arr(i,j,:)**2)/real(NM)-meanA(i,j)**2)
!        sigmaA(i,j)=(maxval(arr(i,j,:))-minval(arr(i,j,:)))    !Histogram maximum width
      enddo
    enddo
    CALL interpolate(meanA(:,:),mean(:,:))
    CALL interpolate(sigmaA(:,:),sigma(:,:)) 
    do j=1,nzt
      write(201,'(10000E13.5)')(mean(i,j),i=1,nxt)
    enddo
    write(201,*);write(201,*)
    do j=1,nzt
      write(201,'(10000E13.5)')(sigma(i,j)*2.,i=1,nxt)
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
        sigma(i,j)=sqrt(sum(arr(i,j,:)**2)/real(NM)-mean(i,j)**2)
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
    
    
    SUBROUTINE interpolate(arrin,arrout)     ! Bilinear interpolation
    USE interp_pb
    USE fd3dparam_pb
    IMPLICIT NONE
    REAL arrin(NLI,NWI),arrout(nxt,nzt)
    REAL DL,DW,ZS,XS,t,u
    INTEGER i,k,ii,kk

    DL=dh*nxt/real(NLI-1)
    DW=dh*nzt/real(NWI-1)
    do k=1,nzt
      ZS=dh*(k-1)+dh/2.
      kk=int(ZS/DW)+1
      u=(ZS-DW*(kk-1))/DW
      do i=1,nxt
        XS=dh*(i-1)+dh/2.
        ii=int(XS/DL)+1
        t=(XS-DL*(ii-1))/DL
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
      dum=(dh*real(nzt-k)+dh/2.)*sin(dip/180.d0*PI)    ! TADY SE TO MUSI OPRAVIT!
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


