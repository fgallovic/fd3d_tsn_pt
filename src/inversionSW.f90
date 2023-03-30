! Parallel tempering inversion subroutines for a model with slip-weakening friction
!-------------------------------------------------------
! Authors: Frantisek Gallovic, Lubica Valentova (2019)
! Charles University in Prague, Faculty of Mathematics and Physics

! This code is published under the GNU General Public License. To any
! licensee is given permission to modify the work, as well as to copy
! and redistribute the work or any derivative version. Still we would
! like to kindly ask you to acknowledge the authors and don't remove
! their names from the code. This code is distributed in the hope
! that it will be useful, but WITHOUT ANY WARRANTY.
! ------------------------------------------------------

    module inversion_com

    real,allocatable,dimension(:,:):: DcI,T0I,TsI     !Test variables (for which misfit is calculated)
    real,allocatable,dimension(:,:,:):: DcA,T0A,TsA   !Array of variables in MC chains:
    real:: StepSizeT0,StepSizeTs,StepSizeD    
    integer:: RUNI,NLI,NWI  
    real,allocatable,dimension(:):: VRA, EgA, ErA, MisfitA,TshiftA, VRgpsA
    real,allocatable,dimension(:,:,:):: ruptimeA,riseA,slipA,schangeA
    real,allocatable :: pgaA(:,:,:),MwA(:),M0A(:),ruptdistA(:,:),MomentRateA(:,:)
    integer randseed,StepType
    
    end module



    module frictionconstraints_com
	
    real:: DcMin,DcMax,strinixMin,strinixMax,peak_xzMin,peak_xzMax
    integer:: ConstraintNucl
    real:: NuclConstraintL,NuclConstraintW,NuclConstraintR
    real :: OverstressConstraint
	
    end module


    subroutine inversion_init()
    use inversion_com
    use mod_ctrl
    use frictionconstraints_com
    use waveforms_com
    use fd3dparam_com
    use source_com
    use mod_pgamisf, only : GMPE_id
    use SlipRates_com
	use PostSeismic_com
    implicit none
	
    open(10,FILE='inputinv.dat')
    
    read(10,*)RUNI
    read(10,*)NLI,NWI
    read(10,*)randseed
    read(10,*)iwaveform, igps
    read(10,*)strinixMin,strinixMax
    read(10,*)peak_xzMin,peak_xzMax
    read(10,*)DcMin,DcMax
    read(10,*)ConstraintNucl,NuclConstraintL,NuclConstraintW,NuclConstraintR,OverstressConstraint	
    read(10,*)StepType,StepSizeT0,StepSizeTs,StepSizeD
	read(10,*)SigmaData,Mwsigma    !Mw constraint applies only when Mwsigma>0 (see evalmisfit())
    if (iwaveform==2) then !for gmpes read additional parameters:
       read(10,*) GMPE_id ! 1 for Zhao, 2 for Boore
       read(10,*) nper !here we define periods for which psa are calculated
       allocate(per(nper))
       read(10,*) per(:) !periods stored here
    endif
    if (iwaveform==4.or.iwaveform==5) then !for ASTFs read additional parameters:
      read(10,*)nper,VSt  !periods for smoothed spectra,S-wave velocity
      allocate(per(nper))
    endif
    close(10)

    allocate(DcI(NLI,NWI),T0I(NLI,NWI),TsI(NLI,NWI))
    allocate(DcA(NLI,NWI,nchains),T0A(NLI,NWI,nchains),TsA(NLI,NWI,nchains))
    allocate(ruptimeA(nxt,nzt,nchains),riseA(nxt,nzt,nchains),slipA(nxt,nzt,nchains),schangeA(nxt,nzt,nchains))
	allocate(VRA(nchains),VRgpsA(nchains),EgA(nchains),ErA(nchains),MisfitA(nchains),M0A(nchains),MwA(nchains),TshiftA(nchains))
    
    !Read GFs and seismograms
    call readGFs()
    if(iwaveform==1)call readwaveforms()
    if(iwaveform==4.or.iwaveform==5)call readastfs()

    !Read postseismic GFs and deformation 
    if (igps==1) call readSGFs()	

    allocate(MomentRateA(nSr,nchains))
    if (iwaveform==2) allocate(ruptdistA(NRseis,nchains),pgaA(NRseis,nper,nchains))

    end


    subroutine AdvanceChain(ichain,T,E,record_mcmc_now,iseed)   ! Old misfit is in E, it is updated if the proposed model is accepted
    use mod_ctrl, only : ifile
    use inversion_com
    use waveforms_com, only : misfit,VR,iwaveform,NRseis,ruptdist,pgaD,Tshift
    use SlipRates_com, only: Mw,M0,MomentRate
    use source_com
    use pml_com
    use fd3dparam_com
    use PostSeismic_com
    USE friction_com    !only for StepType=3
    USE frictionconstraints_com    !only for StepType=3
    implicit none
    
    real,parameter:: eps=1.e-6
    integer ichain,iseed
    real*8 T,E,prop12,prop21
    logical record_mcmc_now
    real*8 newmisfit
    real gasdev,gasdev1,SEold,SEnew
    logical  yn,modelinvalid
    integer i,j,jj,kk

    modelinvalid=.true.
    print *,'searching for model'
jj=0
    do while(modelinvalid)
jj=jj+1
      if(StepType==1)then   !Log-normal step
        prop12=0.
        prop21=0.
        do j=1,NWI
          do i=1,NLI
            T0I(i,j)=T0A(i,j,ichain)*exp(gasdev(iseed)*StepSizeT0)
            TsI(i,j)=TsA(i,j,ichain)*exp(gasdev(iseed)*StepSizeTs)
            DcI(i,j)=DcA(i,j,ichain)*exp(gasdev(iseed)*StepSizeD)
            prop12=prop12+log(T0A(i,j,ichain))+log(TsA(i,j,ichain))+log(DcA(i,j,ichain))
            prop21=prop21+log(T0I(i,j))+log(TsI(i,j))+log(DcI(i,j))
          enddo
        enddo
      elseif(StepType==2)then                 !Normal step
        prop12=log(1.)
        prop21=log(1.)
        do j=1,NWI
          do i=1,NLI
            T0I(i,j)=T0A(i,j,ichain)+gasdev(iseed)*StepSizeT0
            TsI(i,j)=TsA(i,j,ichain)+gasdev(iseed)*StepSizeTs
            DcI(i,j)=DcA(i,j,ichain)+gasdev(iseed)*StepSizeD
          enddo
        enddo
      else                                    !Testing new steps (log-normal + even periodic extension)

        prop12=0.
        prop21=0.
        do j=1,NWI
          do i=1,NLI
            T0I(i,j)=(T0A(i,j,ichain)-strinixMin)*exp(gasdev(iseed)*StepSizeT0)+strinixMin !Log-normal + even periodic extension
            if(T0I(i,j)>strinixMax)then
              T0I(i,j)=strinixMax**2/T0I(i,j)
              if(T0I(i,j)<strinixMin)T0I(i,j)=T0A(i,j,ichain)
            endif
            prop12=prop12+log(T0A(i,j,ichain)-strinixMin)
            prop21=prop21+log(T0I(i,j)-strinixMin)
            
            kk=int(dble(nzt-nfs-nabc-1)/dble(NWI-1)*(j-1)+nabc+1)
            SEold=TsA(i,j,ichain)*normstress(kk)+coh(1,kk)-T0A(i,j,ichain)                 !Assuming coh is constant along strike!
            SEnew=SEold*exp(gasdev(iseed)*StepSizeTs)
            prop12=prop12+log(abs(SEold))
            prop21=prop21+log(abs(SEnew))
            TsI(i,j)=(SEnew+T0I(i,j)-coh(1,kk))/normstress(kk)                             !Assuming coh is constant along strike!

            DcI(i,j)=(DcA(i,j,ichain)-Dcmin)*exp(gasdev(iseed)*StepSizeD)+Dcmin            !Log-normal + even periodic extension
            if(DcI(i,j)>Dcmax)then
                DcI(i,j)=Dcmax**2/DcI(i,j)
                if(DcI(i,j)<Dcmin)DcI(i,j)=DcA(i,j,ichain)
            endif
            prop12=prop12+log(DcA(i,j,ichain)-Dcmin)
            prop21=prop21+log(DcI(i,j)-Dcmin)

!            gasdev1=exp(gasdev(iseed)*StepSizeD)
!            DcI(i,j) = (Dcmin+Dcmax*(DcA(i,j,ichain)-Dcmin)/(Dcmax-DcA(i,j,ichain))*gasdev1)/(1.+(DcA(i,j,ichain)-Dcmin)/(Dcmax-DcA(i,j,ichain))*gasdev1);   !Logit-normal
!            prop12=prop12+log(DcA(i,j,ichain)-Dcmin)+log(Dcmax-DcA(i,j,ichain))
!            prop21=prop21+log(DcI(i,j)-Dcmin)+log(Dcmax-DcI(i,j))

          enddo
        enddo
      endif

      call inversion_modeltofd3d()
      call validatefd3dmodel(modelinvalid)
      continue

!if(jj>1e6)then
!  write(ifile+1000,'(100000E13.5)')MisfitA(ichain),VRA(ichain),T0A(:,:,ichain),TsA(:,:,ichain),DcA(:,:,ichain),EgA(ichain),ErA(ichain),TshiftA(ichain)
!  write(ifile+1000,'(100000E13.5)')0.,0.,T0I(:,:),TsI(:,:),DcI(:,:)
!  stop
!endif

    enddo
    print *,'done'
    call fd3d()
    call syntseis()
	if (igps==1) call CalcSyntGPS
    if (iwaveform==1) then
     call evalmisfit() 
    elseif (iwaveform==2) then
     call evalmisfit2()
    elseif (iwaveform==3) then
     call evalmisfitM()
    elseif (iwaveform==4) then
     call evalmisfitSspec()
    elseif (iwaveform==5) then
     call evalmisfitStime()
    endif
    newmisfit=misfit

    call PT_McMC_accept(T,E,prop21,newmisfit,prop12,yn,iseed)
    if (yn) then  !step accepted
      E=newmisfit
      MisfitA(ichain)=E
      T0A(:,:,ichain)=T0I(:,:)
      TsA(:,:,ichain)=TsI(:,:)
      DcA(:,:,ichain)=DcI(:,:)
      ruptimeA(:,:,ichain)=ruptime(:,:)
      riseA(:,:,ichain)=rise(:,:)
#if defined DIPSLIP
      slipA(:,:,ichain)=slipZ(:,:)
      schangeA(:,:,ichain)=schangeZ(:,:)
#else
      slipA(:,:,ichain)=slipX(:,:)
      schangeA(:,:,ichain)=schangeX(:,:)
#endif
      MomentRateA(:,ichain)=MomentRate(:)
      EgA(ichain)=Eg
      ErA(ichain)=Er
      M0A(ichain)=M0
      MwA(ichain)=Mw
      TshiftA(ichain)=Tshift
      if (iwaveform==2) then
        ruptdistA(:,ichain)=ruptdist(:)
        pgaA(:,:,ichain)=pgaD(:,:)
      endif
      VRA(ichain)=VR
      VRgpsA(ichain)=VRGPS
    endif

    !if (yn.and.(abs(T-1.0)<eps).and.record_mcmc_now) then  !write the accepted step
    if ((abs(T-1.0)<eps).and.record_mcmc_now) then  !write the present step whether accepted or not
      misfit=E
      write(ifile,'(1000000E13.5)')misfit,VRA(ichain),T0A(:,:,ichain),TsA(:,:,ichain),DcA(:,:,ichain),M0A(ichain),EgA(ichain),ErA(ichain),TshiftA(ichain),VRgpsA(ichain)
      flush(ifile)
      write(ifile+2)misfit,VRA(ichain),T0A(:,:,ichain),TsA(:,:,ichain),DcA(:,:,ichain),ruptimeA(nabc+1:nxt-nabc,nabc+1:nzt-nfs,ichain),slipA(nabc+1:nxt-nabc,nabc+1:nzt-nfs,ichain), &
          & riseA(nabc+1:nxt-nabc,nabc+1:nzt-nfs,ichain),schangeA(nabc+1:nxt-nabc,nabc+1:nzt-nfs,ichain),MomentRateA(:,ichain),M0A(ichain),EgA(ichain),ErA(ichain),TshiftA(ichain),VRgpsA(ichain)
      flush(ifile+2)
      if (iwaveform==2) then
        write(ifile*10) (misfit,mwA(ichain),ruptdistA(jj,ichain),pgaA(jj,:,ichain)/100., jj=1,nrseis)
        flush(ifile*10)
      endif
    endif
    
    end

	subroutine validatefd3dmodel(modelinvalid)
    use fd3dparam_com
    use friction_com
    use frictionconstraints_com
    use pml_com
    use source_com, only: ioutput
    use medium_com
    implicit none
    logical  modelinvalid
    real x,z,rr,x0,z0
    integer i,j,nuclOK,ncent,nuclsize,meanoverstress
    real, allocatable :: strengthexcess1(:,:)

    modelinvalid=.true.

!   write(*,*) minval(Dc(nabc+1:nzt-nfs,nabc+1:nxt-nabc)),maxval(Dc(nabc+1:nzt-nfs,nabc+1:nxt-nabc))
!   write(*,*) minval(striniX(nabc+1:nzt-nfs,nabc+1:nxt-nabc)),maxval(striniX(nabc+1:nzt-nfs,nabc+1:nxt-nabc))
!Constraints on Min/Max values
    do j=nabc+1,nzt-nfs
      do i=nabc+1,nxt-nabc
        if(Dc(i,j)<DcMin.or.Dc(i,j)>DcMax)then
!          write(*,*)'Dc',i,j,Dc(i,j),Dcmin,DcMax
          return
        endif
#if defined DIPSLIP
        if(striniZ(i,j)<strinixMin.or.striniZ(i,j)>strinixMax)then
!         write(*,*)'Strinix',i,j,striniZ(i,j),strinixMin,strinixMax
#else
        if(striniX(i,j)<strinixMin.or.striniX(i,j)>strinixMax)then
!          write(*,*)'Strinix',i,j,striniX(i,j),strinixMin,strinixMax
#endif
          return
        endif
        if((peak_xz(i,j)-dyn_xz(i,j))/normstress(j)<peak_xzMin.or.(peak_xz(i,j)-dyn_xz(i,j))/normstress(j)>peak_xzMax)then
!          write(*,*)'Peak_xz',i,j,peak_xz(i,j)
!          write(*,*)'Peak_xz',i,j,(peak_xz(i,j)-dyn_xz(i,j))/normstress(j)
          return
        endif
      enddo
    enddo

!Constraint on the nucleation zone
    if(ConstraintNucl==1)then
      nuclOK=0
      do j=nabc+1,nzt-nfs
        z=dh*real(j-1-nabc)
        do i=nabc+1,nxt-nabc
          x=dh*real(i-1-nabc)
          rr=(x-NuclConstraintL)**2+(z-NuclConstraintW)**2
#if defined DIPSLIP
          if(rr>NuclConstraintR**2.and.peak_xz(i,j)+coh(i,j)<=striniZ(i,j))then
!            write(*,*)'Nucl',x,z,peak_xz(i,j)-striniZ(i,j)
            return    !nucleation outside the nucleation zone
          endif
          if(rr<=NuclConstraintR**2.and.peak_xz(i,j)+coh(i,j)<=striniZ(i,j))nuclOK=1
#else
          if(rr>NuclConstraintR**2.and.peak_xz(i,j)+coh(i,j)<=striniX(i,j))then
!            write(*,*)'Nucl',x,z,peak_xz(i,j)-striniX(i,j)
            return    !nucleation outside the nucleation zone
          endif
          if(rr<=NuclConstraintR**2.and.peak_xz(i,j)+coh(i,j)<=striniX(i,j))nuclOK=1
#endif
        enddo
      enddo
      if(nuclOK==0)return
!Constraint on Mean Overstress (applied only if the constraint is larger than zero):
      if(overstressconstraint>0.)then
        allocate(strengthexcess1(nxt,nzt))
#if defined DIPSLIP
        strengthexcess1(:,:)=striniZ(:,:)-(peak_xz(:,:)+coh(:,:))
#else
        strengthexcess1(:,:)=striniX(:,:)-(peak_xz(:,:)+coh(:,:))
#endif
        nuclsize=dh*dh*COUNT(strengthexcess1(nabc+1:nxt-nabc,nabc+1:nzt-nfs)>=0.)
        meanoverstress=sum(strengthexcess1(nabc+1:nxt-nabc,nabc+1:nzt-nfs),strengthexcess1(nabc+1:nxt-nabc,nabc+1:nzt-nfs)>=0.)*dh*dh/nuclsize
        deallocate(strengthexcess1)    
        if (meanoverstress>overstressconstraint) then
!          print *,'meanoverstress:',meanoverstress
          return !mean overstress is too large
        endif
      endif
    elseif (ConstraintNucl==2) then
!      print *,'checking nucleation constraint 2'
!Constraint on Mean Overstress:   
      allocate(strengthexcess1(nxt,nzt))
#if defined DIPSLIP
      strengthexcess1(:,:)=striniZ(:,:)-(peak_xz(:,:)+coh(:,:))
#else
      strengthexcess1(:,:)=striniX(:,:)-(peak_xz(:,:)+coh(:,:))
#endif
      nuclsize=dh*dh*COUNT(strengthexcess1(nabc+1:nxt-nabc,nabc+1:nzt-nfs)>=0.)
      meanoverstress=sum(strengthexcess1(nabc+1:nxt-nabc,nabc+1:nzt-nfs),strengthexcess1(nabc+1:nxt-nabc,nabc+1:nzt-nfs)>=0.)*dh*dh/nuclsize
      deallocate(strengthexcess1)    
      if (meanoverstress>overstressconstraint) then
!        print *,'meanoverstress:',meanoverstress
        return !mean overstress is too large
      endif

      nuclOK=0
#if defined DIPSLIP
      if (minval(peak_xz(nabc+1:nxt-nabc,nabc+1:nzt-nfs)+coh(nabc+1:nxt-nabc,nabc+1:nzt-nfs)-striniZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs))<=0.) then !nucleation somewhere
#else
      if (minval(peak_xz(nabc+1:nxt-nabc,nabc+1:nzt-nfs)+coh(nabc+1:nxt-nabc,nabc+1:nzt-nfs)-striniX(nabc+1:nxt-nabc,nabc+1:nzt-nfs))<=0.) then !nucleation somewhere
#endif
!        print *,'model nucleating'
        x0=0.
        z0=0.
        ncent=0
        !find  center
        do j=nabc+1,nzt-nfs
          z=dh*real(j-1-nabc)
          do i=nabc+1,nxt-nabc
            x=dh*real(i-1-nabc)
#if defined DIPSLIP
            if (peak_xz(i,j)+coh(i,j)<=striniZ(i,j)) then 
#else
            if (peak_xz(i,j)+coh(i,j)<=striniX(i,j)) then 
#endif
              x0=x0+x
              z0=z0+z
              ncent=ncent+1
            endif
          enddo
        enddo
        x0=x0/ncent
        z0=z0/ncent
        do j=nabc+1,nzt-nfs
          z=dh*real(j-1-nabc)
          do i=nabc+1,nxt-nabc
            x=dh*real(i-1-nabc)
            rr=sqrt((x-x0)**2+(z-z0)**2)
#if defined DIPSLIP
            if (peak_xz(i,j)+coh(i,j)<=striniZ(i,j) .and. rr>NuclConstraintR) then
 !             print *,'nucleation zone:',x0,z0,x,z
              return !nucleations outside 
            endif
            if (peak_xz(i,j)+coh(i,j)<=striniZ(i,j) .and. rr<=NuclConstraintR) NuclOK=1 !nucleation is ok
#else
            if (peak_xz(i,j)+coh(i,j)<=striniX(i,j) .and. rr>NuclConstraintR) then
 !             print *,'nucleation zone:',x0,z0,x,z
              return !nucleations outside 
            endif
            if (peak_xz(i,j)+coh(i,j)<=striniX(i,j) .and. rr<=NuclConstraintR) NuclOK=1 !nucleation is ok
#endif
          enddo
        enddo
!        open(1122,file='nuclcenter',access='append',status='unknown')
!          write(1122,*) x0,z0,meanoverstress
!        close(1122)
      endif
      if (nuclOK==0) then 
!       write(*,*) 'no nucleation',nuclsize
       return
      endif
    endif
!    print *,'nucleation ok'

    if(ioutput==1)then   ! check process zone size
      open(592,FILE='processzone.dat')
      do j=nabc+1,nzt-nfs
        write(592,'(10000E13.5)')(mu1(i,nysc,j)*Dc(i,j)/(peak_xz(i,j)-dyn_xz(i,j))/dh,i=nabc+1,nxt-nabc)
      enddo
      close(592)
    endif
    
    modelinvalid=.false.   !All passed
    end


    subroutine InitiateChain(ichain,E,iseed) !Set all initial models
    use mod_ctrl, only: nchains,rname,ierr,ifile
    use inversion_com
    use waveforms_com, only : misfit,VR,Tshift,iwaveform,ruptdist,pgaD
    use source_com
    use SlipRates_com
    use PostSeismic_com    
    
	implicit none
    
	real*8 E
    real dum
    integer ichain,iseed
    logical modelinvalid

    if (RUNI==1) then
      open(244,FILE='initialmodel.dat')
      read(244,*)dum,dum,T0I(:,:),TsI(:,:),DcI(:,:)
      close(244)
    elseif (RUNI==2) then
      if(ichain==1)open(unit=ifile+1,file=trim(rname),status='old',iostat=ierr)
      read(ifile+1,*)dum,dum,T0I(:,:),TsI(:,:),DcI(:,:)
      if(ichain==nchains)close(ifile+1)
    endif
    
    call inversion_modeltofd3d()
    call validatefd3dmodel(modelinvalid)
    if(modelinvalid)write(*,*)'Initial model violates constraints!'

    call fd3d()
    call syntseis()
	if (igps==1) call CalcSyntGPS()
    if (iwaveform==1) then
      call evalmisfit()
      write(*,*)'Initial model VR: ',VR,' for shift',Tshift,'s'
    elseif (iwaveform==2) then
      call evalmisfit2()
    elseif (iwaveform==3) then
      call evalmisfitM()
    elseif (iwaveform==4) then
      call evalmisfitSspec()
    elseif (iwaveform==5) then
      call evalmisfitStime()
      write(*,*)'Initial model VR: ',VR,' for shift',Tshift,'s'
    endif

    E=misfit
    MisfitA(ichain)=E
    T0A(:,:,ichain)=T0I(:,:)
    TsA(:,:,ichain)=TsI(:,:)
    DcA(:,:,ichain)=DcI(:,:)
    ruptimeA(:,:,ichain)=ruptime(:,:)
    riseA(:,:,ichain)=rise(:,:)
#if defined DIPSLIP
    slipA(:,:,ichain)=slipZ(:,:)
    schangeA(:,:,ichain)=schangeZ(:,:)
#else
    slipA(:,:,ichain)=slipX(:,:)
    schangeA(:,:,ichain)=schangeX(:,:)
#endif
    VRA(ichain)=VR
    VRgpsA(ichain)=VRGPS
    MomentRateA(:,ichain)=MomentRate(:)
    EgA(ichain)=Eg
    ErA(ichain)=Er
    M0A(ichain)=M0
    MwA(ichain)=Mw
    TshiftA(ichain)=Tshift
    if (iwaveform==2) then
      ruptdistA(:,ichain)=ruptdist(:)
      pgaA(:,:,ichain)=pgaD(:,:)
    endif
    
    end


    subroutine saveforrestart()
    use mod_ctrl
    use inversion_com
    
	implicit none
    
	integer ichain
    
    open(unit=ifile+1,file=trim(rname),status='replace',iostat=ierr)
    do ichain=1,nchains
      write(ifile+1,'(100000E13.5)')MisfitA(ichain),VRA(ichain),T0A(:,:,ichain),TsA(:,:,ichain),DcA(:,:,ichain)
    enddo
    close(ifile+1)
    
    end	
	


    subroutine alloc_temp(iseed)
    use mod_ctrl
	
    implicit none
	
    double precision :: r1,aval,bval,dx
    integer :: i,iseed
    real ran3
    
    do i=1,nchains
       r1=ran3(iseed)
       modtemp(i)=thigh-(thigh-tlow)*r1
    enddo

    aval = log(tlow)
    bval = log(thigh)
    dx = (bval-aval)/(nchains-1)
    do i=1,nchains
      r1=ran3(iseed)
      modtemp(i) = exp(aval + r1*(bval-aval))
    end do
      
    i=nchains/4;modtemp(1:i) = tlow                   ! Force first i chains to be at tlow

    end
    
    
    subroutine PT_McMC_accept(T,logPPD1,logQ12,logPPD2,logQ21,yn,iseed)
    
	implicit none
    
	double precision :: logPPD1,logPPD2
    double precision :: logQ21,logQ12
    double precision :: delS
    double precision :: T
    logical          :: yn
    double precision :: a
    integer          :: iseed
    real             :: ran3

    yn = .false.
    delS = (logPPD1-logPPD2)/T
    delS = delS + logQ12 - logQ21
       
    a = ran3(iseed)
    if(log(a).le.delS)yn = .true.      ! swap successful
    
	end subroutine

    
    function gasdev(idum)
    
    integer idum, iset
    real gasdev
    real fac,gset,rsq,v1,v2,ran3
    save iset,gset
    data iset/0/
      
	if (idum.lt.0) iset=0
    if (iset.eq.0) then
1     v1=2.*ran3(idum)-1.
      v2=2.*ran3(idum)-1.
      rsq=v1**2+v2**2
      if(rsq.ge.1..or.rsq.eq.0.)goto 1
      fac=sqrt(-2.*log(rsq)/rsq)
      gset=v1*fac
      gasdev=v2*fac
      iset=1
    else
      gasdev=gset
      iset=0
    endif
    return
    
	end

            FUNCTION ran3(idum)
            INTEGER idum
            INTEGER MBIG,MSEED,MZ
!           REAL MBIG,MSEED,MZ
            REAL ran3,FAC
            PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!           PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
            INTEGER i,iff,ii,inext,inextp,k
            INTEGER mj,mk,ma(55)
!           REAL mj,mk,ma(55)
            SAVE iff,inext,inextp,ma
            DATA iff /0/
!           write(*,*)' idum ',idum
            if(idum.lt.0.or.iff.eq.0)then
               iff=1
               mj=MSEED-iabs(idum)
               mj=mod(mj,MBIG)
               ma(55)=mj
               mk=1
               do 11 i=1,54
                  ii=mod(21*i,55)
                  ma(ii)=mk
                  mk=mj-mk
                  if(mk.lt.MZ)mk=mk+MBIG
                  mj=ma(ii)
11             continue
               do 13 k=1,4
                  do 12 i=1,55
                     ma(i)=ma(i)-ma(1+mod(i+30,55))
                     if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12                continue
13             continue
               inext=0
               inextp=31
               idum=1
            endif
            inext=inext+1
            if(inext.eq.56)inext=1
            inextp=inextp+1
            if(inextp.eq.56)inextp=1
            mj=ma(inext)-ma(inextp)
            if(mj.lt.MZ)mj=mj+MBIG
            ma(inext)=mj
            ran3=mj*FAC
            return
            END FUNCTION
