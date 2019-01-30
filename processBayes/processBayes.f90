    MODULE interp_com
      INTEGER NLI,NWI
    END MODULE
    
    PROGRAM ProcessBayes
    USE fd3dparam_com
    USE friction_com
    USE interp_com
    USE medium_com
    IMPLICIT NONE
    INTEGER,PARAMETER:: NMAX=1e6
    REAL,ALLOCATABLE,DIMENSION(:):: normalstress,VRs,misfits,meansd,meansl,duration,nuclsize,EG,ER,RE,meanoverstress,M0,meanDc,meanStrengthExcess,meanslip,rupturearea,meanruptvel,meanstrength,EGrate
    REAL,ALLOCATABLE,DIMENSION(:,:,:):: DcA,TsA,T0A,SEA,ruptime1,slip1,rise1,schange1,es1,strengthexcess1,rupvel1
    REAL,ALLOCATABLE,DIMENSION(:,:):: dum11,dum12,dum13,dum21,dum22,dum23,dum24,ms1

    REAL bestmisfit,misfitaccept,dum
    REAL vr,mf,x0,z0,x,z,slipmax
    INTEGER NM,NTOT
    INTEGER i,j,k,ml(1),ncent
    
    CALL fd3d_init()
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
    misfitaccept=bestmisfit-log(0.02) !Probability threashold (2% for real data)
!    misfitaccept=bestmisfit-log(0.001) !Probability threashold (1%% for inv1)
!   misfitaccept=bestmisfit-log(0.0001) !Probability threashold (.1%% for pga)
    
    NM=0
    do i=1,k
      if(misfits(i)<=misfitaccept)NM=NM+1
    enddo
    write(*,*)'Number of accepted models: ',NM
    DEALLOCATE(misfits,VRs)

!------ Read accepted models
    ALLOCATE(misfits(NM),VRs(NM),meansd(NM),meansl(NM),duration(NM),nuclsize(NM),EG(NM),ER(NM),RE(NM),meanoverstress(NM),M0(NM),EGrate(NM))
    ALLOCATE(meanDc(NM),meanStrengthExcess(NM),meanslip(NM),rupturearea(NM),meanruptvel(NM),meanstrength(NM))
    allocate(DcA(NLI,NWI,NM),TsA(NLI,NWI,NM),T0A(NLI,NWI,NM),SEA(NLI,NWI,NM))
    allocate(ruptime1(nxt,nzt,NM),slip1(nxt,nzt,NM),rise1(nxt,nzt,NM),schange1(nxt,nzt,NM),es1(nxt,nzt,NM),ms1(nxt,nzt),strengthexcess1(nxt,nzt,NM),rupvel1(nxt,nzt,NM))
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
          SEA(:,j,k)=TsA(:,j,k)*normalstress((j-1)*((nzt-2)/(NWI-1))+1)-T0A(:,j,k)
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
      do j=1,nzt-2
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
      nuclsize(k)=dh*dh*COUNT(strengthexcess1(1:nxt,1:nzt-2,k)>=1.e5)/1.e6
      meanoverstress(k)=sum(strengthexcess1(1:nxt,1:nzt-2,k),strengthexcess1(1:nxt,1:nzt-2,k)>=1.e5)*dh*dh/1.e6/nuclsize(k)/1.e6
      meanstrength(k)=sum(peak_xz(1:nxt,1:nzt-2)*slip1(1:nxt,1:nzt-2,k))/sum(slip1(1:nxt,1:nzt-2,k))

!dependence of slip on Dc
      if(mod(k,10)==0)then
        do j=1,nzt-2,10
          do i=1,nxt,10
            if(slip1(i,j,k)>0.05*slipmax)write(428,'(100E13.5)')slip1(i,j,k),Dc(i,j),strengthexcess1(i,j,k),peak_xz(i,j)
          enddo
        enddo
        write(428,*);write(428,*)
      endif
      
!Rupture velocity
      do j=2,nzt-3
        do i=2,nxt-1
          dum21(i,j)=sqrt((ruptime1(i+1,j,k)-ruptime1(i-1,j,k))**2+(ruptime1(i,j+1,k)-ruptime1(i,j-1,k))**2)/2./dh  !slowness
          rupvel1(i,j,k)=1.e-3/dum21(i,j)
        enddo
      enddo
      meanruptvel(k)=1./(sum(dum21(2:nxt-1,2:nzt-3)*slip1(2:nxt-1,2:nzt-3,k),ruptime1(2:nxt-1,2:nzt-3,k)>1.)/sum(slip1(2:nxt-1,2:nzt-3,k),ruptime1(2:nxt-1,2:nzt-3,k)>1.))/1.e3
      
!      EG(k)=sum(min(Dc(1:nxt,1:nzt-2),slip1(1:nxt,1:nzt-2,k))*peak_xz(1:nxt,1:nzt-2))/2.*dh*dh   ! Dissipated breakdown work (protoze je tam ten min)
      EG(k)=sum(peak_xz(1:nxt,1:nzt-2)*(Dc(1:nxt,1:nzt-2)-(Dc(1:nxt,1:nzt-2)-slip1(1:nxt,1:nzt-2,k))/Dc(1:nxt,1:nzt-2)*max(0.,Dc(1:nxt,1:nzt-2)-slip1(1:nxt,1:nzt-2,k))))/2.*dh*dh   ! Dissipated breakdown work (protoze je tam ten min)
      ER(k)=sum((strinix(1:nxt,1:nzt-2)+peak_xz(1:nxt,1:nzt-2)*max(0.,1.-slip1(1:nxt,1:nzt-2,k)/Dc(1:nxt,1:nzt-2)))*slip1(1:nxt,1:nzt-2,k)-peak_xz(1:nxt,1:nzt-2)*(Dc(1:nxt,1:nzt-2)-(max(0.,Dc(1:nxt,1:nzt-2)-slip1(1:nxt,1:nzt-2,k)))**2/Dc(1:nxt,1:nzt-2)))/2.*dh*dh
      RE(k)=ER(k)/(ER(k)+EG(k)) !Radiation efficiency
            
      meansd(k)=-sum(schange1(1:nxt,1:nzt-2,k)*slip1(1:nxt,1:nzt-2,k))/sum(slip1(1:nxt,1:nzt-2,k))/1.e6
      
      EGrate(k)=sum((peak_xz(1:nxt,1:nzt-2)*(Dc(1:nxt,1:nzt-2)-(Dc(1:nxt,1:nzt-2)-slip1(1:nxt,1:nzt-2,k))/Dc(1:nxt,1:nzt-2)*max(0.,Dc(1:nxt,1:nzt-2)-slip1(1:nxt,1:nzt-2,k))))*slip1(1:nxt,1:nzt-2,k))/2./sum(slip1(1:nxt,1:nzt-2,k))   ! Dissipated breakdown work rate
     
      meansl(k)=sum(strinix(1:nxt,1:nzt-2)/peak_xz(1:nxt,1:nzt-2)*slip1(1:nxt,1:nzt-2,k))/sum(slip1(1:nxt,1:nzt-2,k)) !Stress level (Eq. 4 in Ripperger et al., 2007)
      
      duration(k)=maxval(ruptime1(1:nxt,1:nzt-2,k),slip1(1:nxt,1:nzt-2,k)>0.05*slipmax)-minval(ruptime1(1:nxt,1:nzt-2,k),slip1(1:nxt,1:nzt-2,k)>0.05*slipmax)
      
      M0(k)=sum(slip1(1:nxt,1:nzt-2,k)*mu1(1:nxt,nyt-2,1:nzt-2))*dh*dh

      meanDc(k)=sum(Dc(1:nxt,1:nzt-2)*slip1(1:nxt,1:nzt-2,k))/sum(slip1(1:nxt,1:nzt-2,k))
      
      meanslip(k)=sum(slip1(1:nxt,1:nzt-2,k)*slip1(1:nxt,1:nzt-2,k))/sum(slip1(1:nxt,1:nzt-2,k))
      
      rupturearea(k)=count(slip1(1:nxt,1:nzt-2,k)>0.05*slipmax)*dh*dh
      
      meanStrengthExcess(k)=-sum(strengthexcess1(1:nxt,1:nzt-2,k)*slip1(1:nxt,1:nzt-2,k),strengthexcess1(1:nxt,1:nzt-2,k)<1.e5)/sum(slip1(1:nxt,1:nzt-2,k),strengthexcess1(1:nxt,1:nzt-2,k)<1.e5)/1.e6

!                               1         2       3         4            5        6    7      8       9             10          11      12             13               14       15 16      17             18            19             20
      write(201,'(100E13.5)')misfits(k),VRs(k),meansd(k),duration(k),nuclsize(k),EG(k),ER(k),RE(k),meansl(k),meanoverstress(k),M0(k),meanDc(k),meanStrengthExcess(k),meanslip(k),x0,z0,rupturearea(k),meanruptvel(k),meanstrength(k),Egrate(k)
    enddo
    close(201)

    write(*,*)'Mean stress drop: ',sum(meansd(:))/NM,'+-',sqrt(sum(meansd(:)**2)/NM-(sum(meansd(:))/NM)**2)
    write(*,*)'Fracture energy: ',sum(EG(:))/NM,'+-',sqrt(sum(EG(:)**2)/NM-(sum(EG(:))/NM)**2)
    write(*,*)'Radiated energy: ',sum(ER(:))/NM,'+-',sqrt(sum(ER(:)**2)/NM-(sum(ER(:))/NM)**2)
    write(*,*)'Mean Dc: ',sum(meanDc(:))/NM,'+-',sqrt(sum(meanDc(:)**2)/NM-(sum(meanDc(:))/NM)**2)
    write(*,*)'Mean slip: ',sum(meanslip(:))/NM,'+-',sqrt(sum(meanslip(:)**2)/NM-(sum(meanslip(:))/NM)**2)
    
    
    open(201,FILE='processBayes.slipmodels.dat')
    open(202,FILE='processBayes.strengthexcess.dat')
    do j=1,nzt-2
      do i=1,nxt
        write(201,'(10000E13.5)')slip1(i,j,1:NM)
        write(202,'(10000E13.5)')strengthexcess1(i,j,1:NM)/1.e6
      enddo
    enddo
    close(201)
    close(202)

    open(201,FILE='smaz.txt')
    do j=1,nzt-2
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
    CALL meansigma2(rupvel1(:,:,:),NM)
    close(201)
    
    END
    
    
    SUBROUTINE meansigma(arr,NM)
    USE fd3dparam_com
    USE interp_com
    USE interp_com
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
    do j=1,nzt-2
      write(201,'(10000E13.5)')(mean(i,j),i=1,nxt)
    enddo
    write(201,*);write(201,*)
    do j=1,nzt-2
      write(201,'(10000E13.5)')(sigma(i,j)*2.,i=1,nxt)
!      write(201,'(10000E13.5)')(sigma(i,j)/mean(i,j),i=1,nxt)  !Relative sigma
    enddo
    write(201,*);write(201,*)
    deallocate(meanA,sigmaA,mean,sigma)
    
    END
    
    
    SUBROUTINE meansigma2(arr,NM)
    USE fd3dparam_com
    IMPLICIT NONE
    INTEGER NM
    real arr(nxt,nzt,NM)
    REAL,ALLOCATABLE,DIMENSION(:,:):: mean,sigma
    INTEGER i,j
    
    allocate(mean(nxt,nzt),sigma(nxt,nzt))
    do j=1,nzt-2
      do i=1,nxt
        mean(i,j)=sum(arr(i,j,1:NM))/real(NM)
        sigma(i,j)=sqrt(sum(arr(i,j,:)**2)/real(NM)-mean(i,j)**2)
!        sigma(i,j)=(maxval(arr(i,j,:))-minval(arr(i,j,:)))    !Histogram maximum width
      enddo
    enddo
    do j=1,nzt-2
      write(201,'(10000E13.5)')(mean(i,j),i=1,nxt)
    enddo
    write(201,*);write(201,*)
    do j=1,nzt-2
      write(201,'(10000E13.5)')(sigma(i,j)*2.,i=1,nxt)
!      write(201,'(10000E13.5)')(sigma(i,j)/mean(i,j),i=1,nxt)  !Relative sigma
    enddo
    write(201,*);write(201,*)
    deallocate(mean,sigma)
    
    END
    
    
    SUBROUTINE interpolate(arrin,arrout)     ! Bilinear interpolation
    USE interp_com
    USE fd3dparam_com
    IMPLICIT NONE
    REAL arrin(NLI,NWI),arrout(nxt,nzt)
    REAL DL,DW,ZS,XS,t,u
    INTEGER i,k,ii,kk

    DL=dh*nxt/real(NLI-1)
    DW=dh*(nzt-2)/real(NWI-1)
    do k=1,nzt-2
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
    
