!-------------------------------------------------------
! Processing of Bayesian posterior samples from for fd3d_tsn_pt
!-------------------------------------------------------
! Authors: Frantisek Gallovic and Jan Premus (10/2023)
! Charles University in Prague, Faculty of Mathematics and Physics

! This code is published under the GNU General Public License. To any
! licensee is given permission to modify the work, as well as to copy
! and redistribute the work or any derivative version. Still we would
! like to kindly ask you to acknowledge the authors and don't remove
! their names from the code. This code is distributed in the hope
! that it will be useful, but WITHOUT ANY WARRANTY.
! ------------------------------------------------------

    MODULE interp_pb
      integer NLI,NWI
	  integer:: NLA, NWA
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
!#if defined DIPSLIP
 !     normstress=max(1.e5,8520.*dh*real(nzt-j)*sin(dip/180.*pi))
!#else
      normstress=max(1.e5,16200.*dh*real(nzt-j)*sin(dip/180.*pi))
!#endif
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
    REAL,ALLOCATABLE,DIMENSION(:):: normalstress,VRs,VRs2,misfits,meansd,meansl,duration,nuclsize,EG,ER,RE,meanoverstress,M0,meanDc,meanStrengthExcess,meanslip,rupturearea,meanruptvel,meanstrength,EGrate
    REAL,ALLOCATABLE,DIMENSION(:,:,:):: SEA,ruptime1,slip1,rise1,schange1,es1,strengthexcess1,rupvel1, slip2,schange2
	REAL,ALLOCATABLE,DIMENSION(:,:,:):: T0A,aA,baA,psiA,f0A,fwA,DcA,vwA,viniA, T0A2,aA2,baA2,f0A2,fwA2,DcA2
	REAL,ALLOCATABLE,DIMENSION(:,:):: T0int,aint,baint
	
    REAL,ALLOCATABLE,DIMENSION(:,:):: nucl
	REAL,ALLOCATABLE,DIMENSION(:,:):: dum11,dum12,dum13,dum14,dum15,dum16,dum17,dum18,dum19,dum21,dum22,dum23,dum24,dum25,ms1
	REAL,ALLOCATABLE,DIMENSION(:,:):: fr_interp, tslip_interp
	REAL,ALLOCATABLE,DIMENSION(:,:):: tslip 
	REAL,ALLOCATABLE,DIMENSION(:,:,:)::slipOUTA
    real, allocatable, dimension(:) :: MRate, dum10, dum20,meandepth
    real, allocatable, dimension(:,:) :: MomentRate, centre
	real, allocatable, dimension(:) :: M0A, EgA,ErA,TshiftA,VRgpsA, mfA
	real Nhisto, Mhisto
    REAL bestmisfit,misfitaccept,dum, temp, ntemp, sumx
    real bestvr, vraccept
    REAL vr,mf,x0,z0,x,z,slipmax
	real leng,widt
    INTEGER NM,NTOT,NM2,k2
    INTEGER i,j,k,ml(1),ncent, i_temp
    integer :: nsr,np
    real :: T,dtseis
	character(len=8) :: fm, fi, fj
	
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
    read(10,*) T!,dum,T1,T2
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
    read(10,*) leng,widt
    read(10,*)
    read(10,*) !epicL,epicW
    read(10,*)
    read(10,*) np
    close(10)

    dtseis=T/real(np)
    nSR=int(real(ntfd)/(dtseis/dt))+1
    allocate(MRate(nSR))

    allocate(lam1(nxt,nyt,nzt),mu1(nxt,nyt,nzt),d1(nxt,nyt,nzt))
    allocate(strinix(nxt,nzt),peak_xz(nxt,nzt),Dc(nxt,nzt))

    CALL readcrustalmodel(dip)

    open(10,FILE='inputinv.dat')
    read(10,*)
    read(10,*)NLI,NWI
	
    NLA=51 !Number of discrete points in quasidynamic calculation (horizontal)
    NWA=38 !Number of discrete points in quasidynamic calculation (vertical) 
	
    ALLOCATE(fr_interp(nxt,nzt), tslip_interp(nxt,nzt), tslip(NLA,NWA))
    ALLOCATE(misfits(NMAX),VRs(NMAX),VRs2(NMAX))
    ALLOCATE(dum10(5),dum11(NLI,NWI),dum12(NLI,NWI),dum13(NLI,NWI),dum14(NLI,NWI),dum15(NLI,NWI),dum16(NLI,NWI),dum17(NLI,NWI),dum18(NLI,NWI),dum19(NLI,NWI))
	allocate(dum20(5),dum21(nxt,nzt),dum22(nxt,nzt),dum23(nxt,nzt),dum24(nxt,nzt),dum25(2,NLA*NWA))
    ALLOCATE(normalstress(nzt))

    do i=1,nzt
       normalstress(i)=normstress(i)
    enddo

!------ Learn about misfits
    k=0
	k2=0
    open(101,FILE='sampls.dat',FORM='UNFORMATTED',ACCESS='STREAM')
10  read(101,END=11,ERR=11) mf,vr,dum10(:),dum11(:,:),dum12(:,:),dum13(:,:),dum14(:,:),dum15(:,:),dum16(:,:),dum17(:,:),dum18(:,:),dum19(:,:),& !dynamic model parameters
	  dum21(:,:),dum22(:,:),dum23(:,:),dum24(:,:),MRate(:),dum25(1,:),dum25(2,:),dum20(:)
    k=k+1
    misfits(k)=mf
	vrs2(k)=vr
    if(k==NMAX)stop 'Increase dimension!'
    goto 10
11  close(101)
    NTOT=k
    write(*,*)'Total number of models: ',NTOT
    bestmisfit=minval(abs(misfits(1:NTOT)))
	bestvr=maxval((vrs2(1:NTOT)))
    open(102,FILE='misfitss.dat')
	write(102,'(10000E13.5)')misfits(1:NTOT)
	close(102)
    write(*,*)'Best misfit: ',bestmisfit
	write(*,*)'Best Vr: ', bestvr
    
!    misfitaccept=maxval(misfits(1:ntot))
!    misfitaccept=bestmisfit-log(0.001) !Probability threashold (2% for real data)
!    misfitaccept=bestmisfit-log(0.001) !Probability threashold (1%% for inv1)
  ! misfitaccept=bestmisfit-log(0.0001) !Probability threashold (.1%% for pga)
       !misfitaccept=bestmisfit+150!log(0.0001) 
	   vraccept=0.43!bestvr-0.04
    NM=0
    do i=1,k
      !if(misfits(i)<=misfitaccept)NM=NM+1
	  if(vrs2(i)>=vraccept)NM=NM+1
    enddo
    write(*,*)'Number of accepted models: ',NM
    DEALLOCATE(misfits,VRs)

!------ Read accepted models
    ALLOCATE(misfits(NM),VRs(NM),meansd(NM),meansl(NM),duration(NM),nuclsize(NM),EG(NM),ER(NM),RE(NM),meanoverstress(NM),M0(NM),EGrate(NM),meandepth(NM))
    ALLOCATE(meanDc(NM),meanStrengthExcess(NM),meanslip(NM),rupturearea(NM),meanruptvel(NM),meanstrength(NM))
    allocate(nucl(5,NM),T0A(NLI,NWI,NM),aA(NLI,NWI,NM),baA(NLI,NWI,NM),psiA(NLI,NWI,NM),f0A(NLI,NWI,NM),fwA(NLI,NWI,NM),DcA(NLI,NWI,NM),vwA(NLI,NWI,NM),viniA(NLI,NWI,NM))
	allocate(T0A2(NLI,NWI,NM),aA2(NLI,NWI,NM),baA2(NLI,NWI,NM),f0A2(NLI,NWI,NM),fwA2(NLI,NWI,NM),DcA2(NLI,NWI,NM),slip2(nxt,nzt,NM),schange2(nxt,nzt,NM))
    allocate(ruptime1(nxt,nzt,NM),slip1(nxt,nzt,NM),rise1(nxt,nzt,NM),schange1(nxt,nzt,NM),es1(nxt,nzt,NM),ms1(nxt,nzt),strengthexcess1(nxt,nzt,NM),rupvel1(nxt,nzt,NM))
    allocate(MomentRate(nSr,NM),centre(2,NM))
	allocate(M0A(NM),EgA(NM),ErA(NM),TshiftA(NM),VRgpsA(NM), mfA(NTOT))
    allocate(slipOUTA(2,NLA*NWA,NM))
    open(101,FILE='sampls.dat',FORM='UNFORMATTED',ACCESS='STREAM')
    open(111,FILE='samples.dat',action='read')
    k=0

    do i=1,NTOT
      !read(101)mf,vr,dum11(:,:),dum12(:,:),dum13(:,:),dum21(:,:),dum22(:,:),dum23(:,:),dum24(:,:),Mrate(:)
	  read(101) mf,vr,dum10(:),dum11(:,:),dum12(:,:),dum13(:,:),dum14(:,:),dum15(:,:),dum16(:,:),dum17(:,:),dum18(:,:),dum19(:,:),& !dynamic model parameters
	  dum21(:,:),dum22(:,:),dum23(:,:),dum24(:,:),MRate(:),dum25(1,:),dum25(2,:),dum20(:)
      mfA(i)=mf
	  
      !if(mf<=misfitaccept)then
	  if(vr>=vraccept)then
        k=k+1
        misfits(k)=mf
        VRs(k)=vr
        !T0A(:,:,k)=dum11(:,:)
        !TsA(:,:,k)=dum12(:,:)
        !DcA(:,:,k)=dum13(:,:)
		
		
		nucl(:,k)=dum10(:)
		T0A(:,:,k)=dum11(:,:)
		aA(:,:,k)=dum12(:,:)
		baA(:,:,k)=dum13(:,:)
		psiA(:,:,k)=dum14(:,:)
		f0A(:,:,k)=dum15(:,:)
		fwA(:,:,k)=dum16(:,:)
		DcA(:,:,k)=dum17(:,:)
		vwA(:,:,k)=dum18(:,:)
		viniA(:,:,k)=dum19(:,:)
		
        
      !  do j=1,NWI
      !    SEA(:,j,k)=TsA(:,j,k)*normalstress((j-1)*((nzt)/(NWI-1))+1)-T0A(:,j,k)
      !  enddo
        
        ruptime1(:,:,k)=dum21(:,:)
        slip1(:,:,k)=dum22(:,:)
        rise1(:,:,k)=dum23(:,:)
        schange1(:,:,k)=dum24(:,:)
		slipOUTA(1,:,k)=dum25(1,:)
		slipOUTA(2,:,k)=dum25(2,:)
        MomentRate(:,k)=Mrate(:)
		M0A(k)=dum20(1)
		EgA(k)=dum20(2)
		ErA(k)=dum20(3)
		TshiftA(k)=dum20(4)
		VRgpsA(k)=dum20(5)
		

      endif
    enddo
    close(101)
	close(111)

		nucl(:,k)=dum10(:)
		T0A(:,:,k)=dum11(:,:)
		aA(:,:,k)=dum12(:,:)
		baA(:,:,k)=dum13(:,:)
		psiA(:,:,k)=dum14(:,:)
		f0A(:,:,k)=dum15(:,:)
		fwA(:,:,k)=dum16(:,:)
		DcA(:,:,k)=dum17(:,:)
		vwA(:,:,k)=dum18(:,:)
		viniA(:,:,k)=dum19(:,:)
		


! Save best and worst accepted models
    open(201,FILE='forwardmodel.best.dat')
    ml(:)=minloc(misfits(:))
    write(201,'(10000E13.5)')misfits(ml(1)),VRs(ml(1)),nucl(:,ml(1)),T0A(:,:,ml(1)),aA(:,:,ml(1)),baA(:,:,ml(1)),psiA(:,:,ml(1)), &
	f0A(:,:,ml(1)),fwA(:,:,ml(1)), DcA(:,:,ml(1)), vwA(:,:,ml(1)),viniA(:,:,ml(1))
    close(201)
    open(201,FILE='forwardmodel.worst.dat')
    ml(:)=maxloc(misfits(:))    
	write(201,'(10000E13.5)')misfits(ml(1)),VRs(ml(1)),nucl(:,ml(1)),T0A(:,:,ml(1)),aA(:,:,ml(1)),baA(:,:,ml(1)),psiA(:,:,ml(1)), &
	f0A(:,:,ml(1)),fwA(:,:,ml(1)), DcA(:,:,ml(1)), vwA(:,:,ml(1)),viniA(:,:,ml(1))
    close(201)
    
    open(201,FILE='processBayes.dat')
	!open(202,FILE='HistoValues.dat')
	
	fm='(I5.5)'
	
	open(202,FILE='baprofile.dat')
	open(2021,FILE='baprofile1.dat')
	open(2022,FILE='baprofile2.dat')
	open(2023,FILE='baprofile3.dat')
	open(203,FILE='aprofile.dat')
	open(204,FILE='t0profile.dat')
	open(205,FILE='dcprofile.dat')
	open(206,FILE='f0profile.dat')
	open(207,FILE='viniprofile.dat')
	open(208,FILE='vwprofile.dat')
	open(209,FILE='cslip.dat')
	open(2091,FILE='forwardmodel.picked.dat')
	i_temp=0
	
	
    do k=1,NM
      strinix=0.
      ms1=0.
      peak_xz=0.
      !CALL interpolate(T0A(:,:,k),strinix(:,:))
      !CALL interpolate(TSA(:,:,k),ms1(:,:))
      !CALL interpolate(DcA(:,:,k),Dc(:,:))
      !do i=1,nzt
      !  peak_xz(:,i)=ms1(:,i)*normalstress(i)
      !enddo
      !es1(:,:,k)=(peak_xz(:,:)-strinix(:,:))/max(1.,strinix(:,:))
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
      
	  centre(1,k)=0.
	  centre(2,k)=0.
	  sumx=0.
	  do j=1,nzt
        z=dh*(real(j)-0.5)
        do i=1,nxt
          x=dh*(real(i)-0.5)
			centre(1,k)=centre(1,k)+x*slip1(i,j,k)*mu1(i,nyt,j)
			centre(2,k)=centre(2,k)+z*slip1(i,j,k)*mu1(i,nyt,j)
			sumx=sumx+slip1(i,j,k)*mu1(i,nyt,j)
		enddo
	  enddo
 centre(1,k)=centre(1,k)/sumx
	  centre(2,k)=centre(2,k)/sumx
	  M0(k)=sum(slip1(:,:,k)*mu1(:,nyt,:))*dh*dh  !M0A(k)
	if ((centre(2,k)<15000-5430) .AND. (centre(2,k)>15000-5650) .AND. (M0(k)>1.85e18) .AND. (M0(k)<1.95e18)) then
		 k2=k2+1
!		if ((centre(2,k)<15000-0) .AND. (centre(2,k)>15000-15000) .AND. (M0(k)>0) .AND. (M0(k)<100e18)) then	 
	 !NM=NM+1
!Rupture velocity
      do j=2,nzt-1
        do i=2,nxt-1
          dum21(i,j)=sqrt((ruptime1(i+1,j,k)-ruptime1(i-1,j,k))**2+(ruptime1(i,j+1,k)-ruptime1(i,j-1,k))**2)/2./dh  !slowness
          rupvel1(i,j,k)=1.e-3/dum21(i,j)
        enddo
      enddo
      meanruptvel(k)=1./(sum(dum21(2:nxt-1,2:nzt-1)*slip1(2:nxt-1,2:nzt-1,k),ruptime1(2:nxt-1,2:nzt-1,k)>1.)/sum(slip1(2:nxt-1,2:nzt-1,k),ruptime1(2:nxt-1,2:nzt-1,k)>1.))/1.e3
      EG(k)=EgA(k)
	  ER(k)=ErA(k)
     ! EG(k)=sum(peak_xz(:,:)*(Dc(:,:)-(Dc(:,:)-slip1(:,:,k))/Dc(:,:)*max(0.,Dc(:,:)-slip1(:,:,k))))/2.*dh*dh   ! Dissipated breakdown work (protoze je tam ten min)
     ! ER(k)=sum((strinix(:,:)+peak_xz(:,:)*max(0.,1.-slip1(:,:,k)/Dc(:,:)))*slip1(:,:,k)-peak_xz(:,:)*(Dc(:,:)-(max(0.,Dc(:,:)-slip1(:,:,k)))**2/Dc(:,:)))/2.*dh*dh
      RE(k)=ER(k)/(ER(k)+EG(k)) !Radiation efficiency
            
      meansd(k)=-sum(schange1(:,:,k)*slip1(:,:,k))/sum(slip1(:,:,k))/1.e6

      !EGrate(k)=sum((peak_xz(1:nxt,1:nzt-2)*(Dc(1:nxt,1:nzt-2)-(Dc(1:nxt,1:nzt-2)-slip1(1:nxt,1:nzt-2,k))/Dc(1:nxt,1:nzt-2)*max(0.,Dc(1:nxt,1:nzt-2)-slip1(1:nxt,1:nzt-2,k))))*slip1(1:nxt,1:nzt-2,k))/2./sum(slip1(1:nxt,1:nzt-2,k))   ! Dissipated breakdown work rate
      
      meansl(k)=0.!sum(strinix(:,:)/peak_xz(:,:)*slip1(:,:,k))/sum(slip1(:,:,k)) !Stress level (Eq. 4 in Ripperger et al., 2007)
      
      duration(k)=maxval(ruptime1(:,:,k),slip1(:,:,k)>0.05*slipmax)-minval(ruptime1(:,:,k),slip1(:,:,k)>0.05*slipmax)

      meanDc(k)=sum(Dc(:,:)*slip1(:,:,k))/sum(slip1(:,:,k))
      
      meanslip(k)=sum(slip1(:,:,k)*slip1(:,:,k))/sum(slip1(:,:,k))
      
      rupturearea(k)=count(slip1(:,:,k)>0.05*slipmax)*dh*dh
	  EGrate(k)=EG(k)/rupturearea(k)

      
      meanStrengthExcess(k)=0.!-sum(strengthexcess1(:,:,k)*slip1(:,:,k),strengthexcess1(:,:,k)<1.e5)/sum(slip1(:,:,k),strengthexcess1(:,:,k)<1.e5)/1.e6
      meanoverstress(k)=0.
	  meandepth(k)=nucl(2,k)
	  
      slip2(:,:,k2)=dum22(:,:)
      schange2(:,:,k2)=dum24(:,:)

	  T0A2(:,:,k2)=dum11(:,:)
	  aA2(:,:,k2)=dum12(:,:)
	  baA2(:,:,k2)=dum13(:,:)
	  f0A2(:,:,k2)=dum15(:,:)
	  DcA2(:,:,k2)=dum17(:,:)  

!                               1         2       3         4            5        6    7      8       9             10          11      12             13               14       15 16      17             18            19            20        21           22            23			24			25
      write(201,'(100E13.5)')misfits(k),VRs(k), meansd(k),duration(k),nuclsize(k),EG(k),ER(k),RE(k),meansl(k),meanoverstress(k),M0(k),meanDc(k),meanStrengthExcess(k),meanslip(k),x0,z0,rupturearea(k),meanruptvel(k),meanstrength(k),Egrate(k), VRgpsA(k),  meandepth(k), nucl(3,k), centre(1,k), centre(2,k) 
	  write(209,'(100000E13.5)')slip1(:,:,k) 
	  i_temp=i_temp+1
	  if (i_temp>1) then
	  write(2091,'(1000000E13.5)')misfits(k),VRs(k),nucl(:,k),T0A(:,:,k),aA(:,:,k),baA(:,:,k),psiA(:,:,k), &
f0A(:,:,k),fwA(:,:,k), DcA(:,:,k), vwA(:,:,k),viniA(:,:,k)
	  i_temp=0
	  endif
	 do j=1,NWA
	   do i=1,NLA
		 tslip(i,j)=slipOUTA(2,(j-1)*NLA+i,k)
	   enddo
	 enddo
	 
	 call interpolate2(tslip(:,:), tslip_interp(:,:))
	 
	 call interpolate(baA(:,:,k), fr_interp(:,:))
	 do j=1,nzt
	 
	  temp=0.
	  ntemp=0
	   do i=1,nxt
		
		if (tslip_interp(i,j)>0.3) then

		  ntemp=ntemp+1
		 temp=temp+fr_interp(i,j)
		endif
		
	   enddo
      write(202,'(10000E13.5)',advance="no")temp/ntemp
	 enddo
	 write(202,*)	
	 
	 
	 do j=1,nzt
	 
	  temp=0.
	  ntemp=0
	   do i=1,75
		
		if (tslip_interp(i,j)>0.3) then

		  ntemp=ntemp+1
		 temp=temp+fr_interp(i,j)
		endif
		
	   enddo
      write(2021,'(10000E13.5)',advance="no")temp/ntemp
	 enddo
	 write(2021,*)	
	 
	 do j=1,nzt
	 
	  temp=0.
	  ntemp=0
	   do i=76,114
		
		if (tslip_interp(i,j)>0.3) then

		  ntemp=ntemp+1
		 temp=temp+fr_interp(i,j)
		endif
		
	   enddo
      write(2022,'(10000E13.5)',advance="no")temp/ntemp
	 enddo
	 write(2022,*)	
	 
	 do j=1,nzt
	 
	  temp=0.
	  ntemp=0
	   do i=115,151
		
		if (tslip_interp(i,j)>0.3) then

		 ntemp=ntemp+1
		 temp=temp+fr_interp(i,j)
		endif
		
	   enddo
      write(2023,'(10000E13.5)',advance="no")temp/ntemp
	 enddo
	 write(2023,*)	

	 
	 
	 
	 
	 call interpolate(aA(:,:,k), fr_interp(:,:))
	 do j=1,nzt
	 
	  temp=0.
	  ntemp=0
	   do i=1,nxt
		
		if (tslip_interp(i,j)>0.3) then

		  ntemp=ntemp+1
		 temp=temp+fr_interp(i,j)
		endif
		
	   enddo
      write(203,'(10000E13.5)',advance="no")temp/ntemp
	 enddo
	 write(203,*)	
	 
	 call interpolate(T0A(:,:,k), fr_interp(:,:))
	 do j=1,nzt
	 
	  temp=0.
	  ntemp=0
	   do i=1,nxt
		
		if (tslip_interp(i,j)>0.3) then

		  ntemp=ntemp+1
		 temp=temp+fr_interp(i,j)
		endif
		
	   enddo
      write(204,'(10000E13.5)',advance="no")temp/ntemp
	 enddo
	 write(204,*)	
	 
	 call interpolate(DcA(:,:,k), fr_interp(:,:))
	 do j=1,nzt
	 
	  temp=0.
	  ntemp=0
	   do i=1,nxt
		
		if (tslip_interp(i,j)>0.3) then

		  ntemp=ntemp+1
		 temp=temp+fr_interp(i,j)
		endif
		
	   enddo
      write(205,'(10000E13.5)',advance="no")temp/ntemp
	 enddo
	 write(205,*)	
	 
	 call interpolate(f0A(:,:,k), fr_interp(:,:))
	 do j=1,nzt
	 
	  temp=0.
	  ntemp=0
	   do i=1,nxt
		
		if (tslip_interp(i,j)>0.3) then

		  ntemp=ntemp+1
		 temp=temp+fr_interp(i,j)
		endif
		
	   enddo
      write(206,'(10000E13.5)',advance="no")temp/ntemp
	 enddo
	 write(206,*)	
	 
	 call interpolate(viniA(:,:,k), fr_interp(:,:))
	 do j=1,nzt
	 
	  temp=0.
	  ntemp=0
	   do i=1,nxt
		
		if (tslip_interp(i,j)>0.3) then

		  ntemp=ntemp+1
		 temp=temp+fr_interp(i,j)
		endif
		
	   enddo
      write(207,'(10000E13.5)',advance="no")temp/ntemp
	 enddo
	 write(207,*)	
	 
	 
	 call interpolate(vwA(:,:,k), fr_interp(:,:))
	 do j=1,nzt
	 
	  temp=0.
	  ntemp=0
	   do i=1,nxt
		
		if (tslip_interp(i,j)>0.3) then

		  ntemp=ntemp+1
		 temp=temp+fr_interp(i,j)
		endif
		
	   enddo
      write(208,'(10000E13.5)',advance="no")temp/ntemp
	 enddo
	 write(208,*)	
	!else
	!M0(k)=0
	!centre(1,k)=0
	!centre(2,k)=0
	
	
	
	
	
	endif
	
	
	enddo
	NM2=k2
	
    close(201)
    close(202)
    close(2021)
    close(2022)
	close(2023)
	close(203)
	close(204)	
	close(205)	
	close(206)	
	close(207)	
	close(208)	
	close(209)
	close(2091)
	
    open(202,FILE='bahisto.dat')
	do i=1,NLI
	  do j=1,NWI
	    write(202,'(10000E13.5)')baA(i,j,:)
	  enddo
	enddo
	close(202)
	
	open(202,FILE='ahisto.dat')
	do i=1,NLI
	  do j=1,NWI
	    write(202,'(10000E13.5)')aA(i,j,:)
	  enddo
	enddo
	close(202)

    open(202,FILE='T0histo.dat')
	do i=1,NLI
	  do j=1,NWI
	    write(202,'(10000E13.5)')T0A(i,j,:)
	  enddo
	enddo
	close(202)
	
	open(202,FILE='Lhisto.dat')
	do i=1,NLI
	  do j=1,NWI
	    write(202,'(10000E13.5)')DcA(i,j,:)
	  enddo
	enddo
	close(202)
	
	open(202,FILE='f0histo.dat')
	do i=1,NLI
	  do j=1,NWI
	    write(202,'(10000E13.5)')f0A(i,j,:)
	  enddo
	enddo
	close(202)
	
    open(202,FILE='vinihisto.dat')
	do i=1,NLI
	  do j=1,NWI
	    write(202,'(10000E13.5)')viniA(i,j,:)
	  enddo
	enddo
	close(202)
	
	
    open(202,FILE='slips.dat')
	do i=1,NM
	    write(202,'(10000E13.5)')slipOUTA(1,:,i)
	enddo
	close(202)
	
    open(202,FILE='afterslips.dat')
	do i=1,NM
	    write(202,'(10000E13.5)')slipOUTA(2,:,i)
	enddo
	close(202)
	
    open(202,FILE='rupvel.dat')
	do i=1,NM
	    write(202,'(100000E13.5)')rupvel1(:,:,i)
	enddo
	close(202)	
	
	open(202,FILE='schange.dat')
	do i=1,NM
	    write(202,'(100000E13.5)')schange1(:,:,i)
	enddo
	close(202)	


    write(*,*)'Mean stress drop: ',sum(meansd(:))/NM,'+-',sqrt(sum(meansd(:)**2)/NM-(sum(meansd(:))/NM)**2)
    write(*,*)'Fracture energy: ',(sum(EG(:))/NM)/(sum(rupturearea(:))/NM),'+-',(sum(EG(:))/NM)/(sum(rupturearea(:))/NM)*sqrt((sqrt(sum(EG(:)**2)/NM-(sum(EG(:))/NM)**2)/(sum(EG(:))/NM))**2 + (sqrt(sum(rupturearea(:)**2)/NM-(sum(rupturearea(:))/NM)**2)/(sum(rupturearea(:))/NM))**2)
    write(*,*)'Radiated energy: ',(sum(ER(:))/NM)/(sum(rupturearea(:))/NM),'+-',(sum(ER(:))/NM)/(sum(rupturearea(:))/NM)*sqrt((sqrt(sum(ER(:)**2)/NM-(sum(ER(:))/NM)**2)/(sum(ER(:))/NM))**2 + (sqrt(sum(rupturearea(:)**2)/NM-(sum(rupturearea(:))/NM)**2)/(sum(rupturearea(:))/NM))**2)
    write(*,*)'Mean Dc: ',sum(meanDc(:))/NM,'+-',sqrt(sum(meanDc(:)**2)/NM-(sum(meanDc(:))/NM)**2)
    write(*,*)'Mean slip: ',sum(meanslip(:))/NM,'+-',sqrt(sum(meanslip(:)**2)/NM-(sum(meanslip(:))/NM)**2)
    write(*,*)'Mean Ruptured Area: ',sum(rupturearea(:))/NM,'+-',sqrt(sum(rupturearea(:)**2)/NM-(sum(rupturearea(:))/NM)**2)
    write(*,*)'Mean Moment: ',sum(M0(:))/NM,'+-',sqrt(sum(M0(:)**2)/NM-(sum(M0(:))/NM)**2)
    write(*,*)'Mean Centroid Depth: ',sum(centre(2,:))/NM,'+-',sqrt(sum(centre(2,:)**2)/NM-(sum(centre(2,:))/NM)**2)

    
    
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

    open(201,FILE='processBayes.MomentRates.dat')
    do k=1,NM
      do j=1,nSR
        write(201,'(10000E13.5)')dtseis*(j-1),MomentRate(j,k)
      enddo
      write(201,*)
      write(201,*)
    enddo
    close(201)
    
    open(201,FILE='processBayes.sigma.dat')
    open(202,FILE='processBayes.mean.dat')
    
	CALL meansigma(T0A(:,:,:),NM)
	CALL meansigma(aA(:,:,:),NM)
	CALL meansigma(baA(:,:,:),NM)
	CALL meansigma(f0A(:,:,:),NM)
	CALL meansigma(fwA(:,:,:),NM)
	CALL meansigma(DcA(:,:,:),NM)
	CALL meansigma(vwA(:,:,:),NM)
	CALL meansigma(viniA(:,:,:),NM)
	
    close(201)
	close(202)
	
	
	open(201,FILE='processBayes.sigma2.dat')
    open(202,FILE='processBayes.mean2.dat')
    
	CALL meansigma2(slip1(:,:,:),NM)
	CALL meansigma2(schange1(:,:,:),NM)
	
    close(201)
	close(202)
	
	open(201,FILE='processBayes.sigma.subset.dat')
    open(202,FILE='processBayes.mean.subset.dat')
    
	CALL meansigma(T0A2(:,:,:),NM2)
	CALL meansigma(aA2(:,:,:),NM2)
	CALL meansigma(baA2(:,:,:),NM2)
	CALL meansigma(f0A2(:,:,:),NM2)
	CALL meansigma(DcA2(:,:,:),NM2)
	
    close(201)
	close(202)
    
	open(201,FILE='processBayes.sigma2.subset.dat')
    open(202,FILE='processBayes.mean2.subset.dat')
    
	CALL meansigma2(slip2(:,:,:),NM2)
	CALL meansigma2(schange2(:,:,:),NM2)
	
    close(201)
	close(202)	
	
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
      write(202,'(10000E13.5)')(mean(i,j),i=1,nxt)
    enddo
    write(201,*);write(201,*)
    do j=1,nzt
!      write(201,'(10000E13.5)')(sigma(i,j)*2.,i=1,nxt)
      write(201,'(10000E13.5)')(sigma(i,j)/abs(mean(i,j)),i=1,nxt)  !Relative sigma
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
      write(202,'(10000E13.5)')(mean(i,j),i=1,nxt)
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
    
    SUBROUTINE interpolate2(arrin,arrout)     ! Bilinear interpolation from the second grid
    USE interp_pb
    USE fd3dparam_pb
    IMPLICIT NONE
    REAL arrin(NLA,NWA),arrout(nxt,nzt)
    REAL DL,DW,ZS,XS,t,u
    INTEGER i,k,ii,kk

    DL=dh*nxt/real(NLA-1)
    DW=dh*nzt/real(NWA-1)
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

    SUBROUTINE interpolate3(arrin,arrout)     ! Bilinear interpolation from the first to second grid
    USE interp_pb
    USE fd3dparam_pb
    IMPLICIT NONE
    REAL arrin(NLI,NWI),arrout(NLA,NWA)
    REAL DL,DW,ZS,XS,t,u
    INTEGER i,k,ii,kk

    DL=dh*nxt/real(NLA-1)
    DW=dh*nzt/real(NWA-1)
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


