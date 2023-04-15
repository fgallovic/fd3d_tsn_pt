! Parallel tempering inversion subroutines for a model with rate-and-state friction
!-------------------------------------------------------
! Authors: Frantisek Gallovic, Jan Premus (2020)
! Charles University in Prague, Faculty of Mathematics and Physics

! This code is published under the GNU General Public License. To any
! licensee is given permission to modify the work, as well as to copy
! and redistribute the work or any derivative version. Still we would
! like to kindly ask you to acknowledge the authors and don't remove
! their names from the code. This code is distributed in the hope
! that it will be useful, but WITHOUT ANY WARRANTY.
! ------------------------------------------------------

    module inversion_com
	
    real,allocatable,dimension(:,:):: T0I, aI, baI, psiI, f0I, fwI, DcI, vwI, viniI     !Test variables (for which misfit is calculated)
    real,allocatable,dimension(:,:,:):: T0A, aA, baA, psiA, f0A, fwA, DcA, vwA, viniA   !Array of variables in MC chains:
	real,allocatable,dimension(:)::nucl
	real,allocatable,dimension(:,:)::nuclA
    real::StepSizeT0, StepSizea, StepSizeba, StepSizepsi, StepSizef0, StepSizefw, StepSizeDc, StepSizevw, StepSizevini 
	real,allocatable,dimension(:):: StepSizenucl
    
	integer:: RUNI,NLI,NWI  
!    real,allocatable,dimension(:,:):: TsI   !Test variables (for which misfit is calculated)
!    real,allocatable,dimension(:,:,:):: TsA  !Array of variables in MC chains:
    real,allocatable,dimension(:):: VRA, EgA, ErA, MisfitA,TshiftA, VRgpsA
	real,allocatable,dimension(:,:,:):: ruptimeA,riseA,slipA,schangeA
	real,allocatable,dimension(:,:,:):: slipOUTA
	real,allocatable :: pgaA(:,:,:),MwA(:),M0A(:),ruptdistA(:,:),MomentRateA(:,:)
    integer randseed,StepType

	  
    end module

    module frictionconstraints_com
    
	real:: T0Min, aMin, baMin, psiMin, f0Min, fwMin, DcMin, vwMin, viniMax
    real:: T0Max, aMax, baMax, psiMax, f0Max, fwMax, DcMax, vwMax, viniMin
	real, allocatable, dimension(:):: nuclMin, nuclMax
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
    USE RATESTATE, only: MDIS, NDIS
    implicit none
	
    open(10,FILE='inputinv.dat')
    
    read(10,*)RUNI
    read(10,*)NLI,NWI
    read(10,*)randseed
    read(10,*)iwaveform, igps

	allocate(StepSizenucl(5), nuclMin(5), nuclMax(5))
	
    read(10,*)T0Min, aMin, baMin, psiMin, f0Min, fwMin, DcMin, vwMin, viniMin
    read(10,*)T0Max, aMax, baMax, psiMax, f0Max, fwMax, DcMax, vwMax, viniMax
	read(10,*)nuclMin(1:5)
	read(10,*)nuclMax(1:5)

    read(10,*)StepType
    read(10,*)StepSizeT0, StepSizea, StepSizeba, StepSizepsi, StepSizef0, StepSizefw, StepSizeDc, StepSizevw, StepSizevini
    read(10,*)StepSizenucl(1),StepSizenucl(2),StepSizenucl(3),StepSizenucl(4),StepSizenucl(5)
	
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

    allocate(T0I(NLI,NWI), aI(NLI,NWI), baI(NLI,NWI), psiI(NLI,NWI), f0I(NLI,NWI), fwI(NLI,NWI), DcI(NLI,NWI), vwI(NLI, NWI), viniI(NLI, NWI))
    allocate(T0A(NLI,NWI,nchains), aA(NLI,NWI,nchains), baA(NLI,NWI,nchains), psiA(NLI,NWI,nchains))
	allocate(f0A(NLI,NWI,nchains), fwA(NLI,NWI,nchains), DcA(NLI,NWI,nchains), vwA(NLI, NWI,nchains), viniA(NLI, NWI,nchains))
	allocate(slipOUTA(2,MDIS*NDIS,nchains))
	allocate(nucl(5), nuclA(5,nchains)) !hx0, hz0, RR2, TT2, perturb
    allocate(ruptimeA(nxt,nzt,nchains),riseA(nxt,nzt,nchains),slipA(nxt,nzt,nchains),schangeA(nxt,nzt,nchains))
	allocate(VRA(nchains),VRgpsA(nchains),EgA(nchains),ErA(nchains),MisfitA(nchains),M0A(nchains),MwA(nchains),TshiftA(nchains))
    
    !Read GFs and seismograms
    call readGFs()
    if(iwaveform==1)call readwaveforms()
	
	!Read postseismic GFs and deformation 
	if (igps==1) call readSGFs()
    if(iwaveform==4.or.iwaveform==5)call readastfs()
	
    allocate(MomentRateA(nSr,nchains))
    if (iwaveform==2) allocate(ruptdistA(NRseis,nchains),pgaA(NRseis,nper,nchains))

    end

    subroutine AdvanceChain(ichain,T,E,record_mcmc_now,iseed)   ! V E je stary misfit, nahradi se pripadne novym, pokud dojde k prijeti kroku
    use mod_ctrl, only : ifile
    use inversion_com
    use waveforms_com, only : misfit,VR,iwaveform,NRseis,ruptdist,pgaD,Tshift
    use SlipRates_com, only: Mw,M0,MomentRate
    use source_com
    use pml_com
    use fd3dparam_com
    use PostSeismic_com
	use frictionconstraints_com, only : viniMin, vwMin
	use RATESTATE, only:slipOUT
	
    implicit none
    
    real,parameter:: eps=1.e-6
    integer ichain,iseed
    real*8 T,E,prop12,prop21
    logical record_mcmc_now
    real*8 newmisfit
    real gasdev
    logical  yn,modelinvalid
    integer i,j,jj

    modelinvalid=.true.
    print *,'searching for model'
    do while(modelinvalid)
      if(StepType==1)then   !Log-normal step
        prop12=0.
        prop21=0.
		!DcI(1,1)=DcA(1,1,ichain)*exp(gasdev(iseed)*StepSizeDc)
		!prop12=prop12+log(DcA(1,1,ichain))
        !prop21=prop21+log(DcI(1,1))
		
        do j=1,NWI
          do i=1,NLI

		    T0I(i,j)=T0A(i,j,ichain)*exp(gasdev(iseed)*StepSizeT0)
            aI(i,j)=aA(i,j,ichain)*exp(gasdev(iseed)*StepSizea)
            f0I(i,j)=f0A(i,j,ichain)*exp(gasdev(iseed)*StepSizef0)
            DcI(i,j)=DcA(i,j,ichain)*exp(gasdev(iseed)*StepSizeDc)
            baI(i,j)=(baA(i,j,ichain)+aA(i,j,ichain))*exp(gasdev(iseed)*StepSizeba)-aA(i,j,ichain)
           ! DcI(i,j)=DcI(1,1)
			
			!normal step for ba
	        !baI(i,j)=baA(i,j,ichain)+gasdev(iseed)*StepSizeba
            prop12=prop12+log(T0A(i,j,ichain))+log(aA(i,j,ichain))+log(f0A(i,j,ichain))+log(DcA(i,j,ichain))+log(baA(i,j,ichain)+aA(i,j,ichain))
            prop21=prop21+log(T0I(i,j))+log(aI(i,j))+log(f0I(i,j))+log(DcI(i,j))+log(baI(i,j)+aI(i,j))			
			!vini and vw change only in velocity strenghtening zone
			if (baI(i,j)<0.) then
			  vwI(i,j)=vwA(i,j,ichain)*exp(gasdev(iseed)*StepSizevw)
			  viniI(i,j)=viniA(i,j,ichain)*exp(gasdev(iseed)*StepSizevini)
			  prop12=prop12+log(vwA(i,j,ichain))+log(viniA(i,j,ichain))
			  prop21=prop21+log(vwI(i,j))+log(viniI(i,j))
			else
			  vwI(i,j)=0.1
			  viniI(i,j)=1.e-12
			endif
			
		 enddo
        enddo
	  
        nucl(1)=nuclA(1,ichain)*exp(gasdev(iseed)*StepSizenucl(1))
        nucl(2)=nuclA(2,ichain)*exp(gasdev(iseed)*StepSizenucl(2))
        nucl(3)=nuclA(3,ichain)*exp(gasdev(iseed)*StepSizenucl(3))
        nucl(5)=nuclA(5,ichain)*exp(gasdev(iseed)*StepSizenucl(5))	  
	  
	    prop12=prop12+log(nuclA(1,ichain))+log(nuclA(2,ichain))+log(nuclA(3,ichain))+log(nuclA(5,ichain))
	    prop21=prop21+log(nucl(1))+log(nucl(2))+log(nucl(3))+log(nucl(5))
		
      else                  !Normal step	  
        prop12=log(1.)
        prop21=log(1.)
        do j=1,NWI
          do i=1,NLI
		  
            T0I(i,j)=T0A(i,j,ichain)+gasdev(iseed)*StepSizeT0
	        
			aI(i,j)=aA(i,j,ichain)+gasdev(iseed)*StepSizea
	        baI(i,j)=baA(i,j,ichain)+gasdev(iseed)*StepSizeba
			
	        !psiI(i,j)=psiA(i,j,ichain)+gasdev(iseed)*StepSizepsi
	        
			f0I(i,j)=f0A(i,j,ichain)+gasdev(iseed)*StepSizef0
	        
			!fwI(i,j)=fwA(i,j,ichain)+gasdev(iseed)*StepSizefw
	        
			DcI(i,j)=DcA(i,j,ichain)+gasdev(iseed)*StepSizeDc
	        
			!vwI(i,j)=vwA(i,j,ichain)+gasdev(iseed)*StepSizevw

          enddo
        enddo
		
		nucl(1)=nuclA(1,ichain)+gasdev(iseed)*StepSizenucl(1)	
		
		nucl(2)=nuclA(2,ichain)+gasdev(iseed)*StepSizenucl(2)	
		
		nucl(3)=nuclA(3,ichain)+gasdev(iseed)*StepSizenucl(3)	
	
		!nuclI(4)=nuclA(1,ichain)+gasdev(iseed)*StepSizenucl(4)	
		
		nucl(5)=nuclA(5,ichain)+gasdev(iseed)*StepSizenucl(5)	
		

		
      endif
	  
      call inversion_modeltofd3d()
      call validatefd3dmodel(modelinvalid)  
      continue
    enddo
    print *,'done'
    call fd3d()
    call syntseis()
	if (igps==1) then
	  call CalcSyntGPS()
	endif	 
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
	    aA(:,:,ichain)=aI(:,:)
	    baA(:,:,ichain)=baI(:,:)
	    psiA(:,:,ichain)=psiI(:,:)
	    f0A(:,:,ichain)=f0I(:,:)
	    fwA(:,:,ichain)=fwI(:,:)
	    DcA(:,:,ichain)=DcI(:,:)
	    vwA(:,:,ichain)=vwI(:,:)
	    viniA(:,:,ichain)=viniI(:,:)
	    nuclA(:,ichain)=nucl(:)
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
	    slipOUTA(:,:,ichain)=slipOUT(:,:)

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
      write(ifile,'(1000000E13.5)')misfit,VRA(ichain),nuclA(1:5, ichain),T0A(:,:,ichain),aA(:,:,ichain),baA(:,:,ichain),psiA(:,:,ichain),f0A(:,:,ichain),fwA(:,:,ichain),DcA(:,:,ichain),vwA(:,:,ichain),viniA(:,:,ichain),M0A(ichain),EgA(ichain),ErA(ichain),TshiftA(ichain),VRgpsA(ichain)
      flush(ifile)
      write(ifile+2)misfit,VRA(ichain),nuclA(1:5,ichain),T0A(:,:,ichain),aA(:,:,ichain),baA(:,:,ichain),psiA(:,:,ichain),f0A(:,:,ichain),fwA(:,:,ichain),DcA(:,:,ichain),vwA(:,:,ichain),viniA(:,:,ichain),ruptimeA(nabc+1:nxt-nabc,nabc+1:nzt-nfs,ichain),slipA(nabc+1:nxt-nabc,nabc+1:nzt-nfs,ichain), &
          & riseA(nabc+1:nxt-nabc,nabc+1:nzt-nfs,ichain),schangeA(nabc+1:nxt-nabc,nabc+1:nzt-nfs,ichain),MomentRateA(:,ichain),slipOUTA(1,:,ichain),slipOUTA(2,:,ichain),M0A(ichain),EgA(ichain),ErA(ichain),TshiftA(ichain),VRgpsA(ichain)
      flush(ifile+2)
      if (iwaveform==2) then
        write(ifile*10) (misfit,MwA(ichain),ruptdistA(jj,ichain),pgaA(jj,:,ichain)/100., jj=1,nrseis)
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
  !  read(10,*)T0Min, aMin, baMin, psiMin, f0Min, fwMin, DcMin, vwMin
  !  read(10,*)T0Max, aMax, baMax, psiMax, f0Max, fwMax, DcMax, vwMax
!    print *,minval(dc),maxval(dc)
!    write(*,*) minval(strinix),maxval(strinix)
!Constraints on Min/Max values
    do j=nabc+1,nzt-nfs
      do i=nabc+1,nxt-nabc

#if defined DIPSLIP
        if(striniZ(i,j)<T0Min.or.striniZ(i,j)>T0Max)then
        ! write(*,*)'Striniz',i,j,striniZ(i,j),T0Min,T0Max
#else
        if(striniX(i,j)<T0Min.or.striniX(i,j)>T0Max)then
        ! write(*,*)'Strinix',i,j,striniX(i,j),T0Min,T0Max
#endif
          return
        endif
		
        if(a(i,j)<aMin.or.a(i,j)>aMax)then
        ! write(*,*)'a',i,j,a(i,j)
          return
        endif
		
		if(ba(i,j)+a(i,j)<baMin.or.ba(i,j)+a(i,j)>baMax)then
          write(*,*)'ba',i,j,ba(i,j)
          return
        endif
		
	!	if(ba(i,j)+a(i,j)<0.)then
     !     write(*,*)'b<0.',i,j,ba(i,j)
      !    return
      !  endif
		
!		if(psi(i,j)<psiMin.or.psi(i,j)>psiMax)then	
!          write(*,*)'psi',i,j,psi(i,j)
 !         return
 !       endif
		
		if(f0(i,j)<f0Min.or.f0(i,j)>f0Max)then
    !      write(*,*)'f0',i,j,f0(i,j)
          return
        endif
		
 !       if(fw(i,j)<fwMin.or.fw(i,j)>fwMax)then
 !         write(*,*)'fw',i,j,fw(i,j)
 !         return
 !       endif
		
		if(Dc(i,j)<DcMin.or.Dc(i,j)>DcMax)then
     !     write(*,*)'Dc',i,j,Dc(i,j),Dcmin,DcMax
          return
        endif
        
		if (ba(i,j)<0.) then
		  if(vw(i,j)<vwMin.or.vw(i,j)>vwMax)then
      !      write(*,*)'vw',i,j,vw(i,j),vwmin,vwMax
            return
          endif
#if defined DIPSLIP
 		  if(wini(i,j)<viniMin.or.wini(i,j)>viniMax)then
       !     write(*,*)'vini',i,j,wini(i,j),vinimin,viniMax
            return
          endif
#else
 		  if(uini(i,j)<viniMin.or.uini(i,j)>viniMax)then
        !    write(*,*)'vini',i,j,uini(i,j),vinimin,viniMax
            return
          endif

#endif		  
		  
        endif 
      enddo
    enddo

    if(hx0-real(nabc)*dh<nuclMin(1).or.hx0-real(nabc)*dh>nuclMax(1))then	
 !      write(*,*)'nucl',i, nucl(i)
        return
    endif	
	
	if(hz0-real(nabc)*dh<nuclMin(2).or.hz0-real(nabc)*dh>nuclMax(2))then
!       write(*,*)'nucl',i, nucl(i)
        return
    endif	
	
	if(RR2<nuclMin(3).or.RR2>nuclMax(3))then	
    !   write(*,*)'nucl3', RR2, nuclMin(3), nuclMax(3)
       return
    endif	
	
!	if(TT2<nuclMin(4).or.TT2>nuclMax(4))then
!       write(*,*)'nucl',i, nucl(i)
!        return
!    endif		
	
	if(perturb<nuclMin(5).or.perturb>nuclMax(5))then
     !  write(*,*)'nucl5',i, perturb
       return
    endif	
    
    modelinvalid=.false.   !All passed
    end

    
    subroutine InitiateChain(ichain,E,iseed) !Set all initial models
    use mod_ctrl, only: nchains,rname,ierr,ifile
    use inversion_com
    use waveforms_com, only : misfit,VR,Tshift,iwaveform,ruptdist,pgaD,Tshift
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
      read(244,*)dum,dum,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:)
      close(244)
    elseif (RUNI==2) then
      if(ichain==1)open(unit=ifile+1,file=trim(rname),status='old',iostat=ierr)
      read(ifile+1,*)dum,dum,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:)
      if(ichain==nchains)close(ifile+1)
    endif
    
    call inversion_modeltofd3d()
    call validatefd3dmodel(modelinvalid)
    if(modelinvalid)write(*,*)'Initial model violates constraints!'

    call fd3d()
    call syntseis()
	if (igps==1)then
	  call CalcSyntGPS()
	endif	 
    if (iwaveform==1) then
      call evalmisfit()
      write(*,*)'Initial model VR: ',VR,' for shift',Tshift,'s', ', GPS VR: ', VRGPS
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
    aA(:,:,ichain)=aI(:,:)
    baA(:,:,ichain)=baI(:,:)
	psiA(:,:,ichain)=psiI(:,:)
    f0A(:,:,ichain)=f0I(:,:)
    fwA(:,:,ichain)=fwI(:,:)
	DCA(:,:,ichain)=DCI(:,:)
    vwA(:,:,ichain)=vwI(:,:)
    viniA(:,:,ichain)=viniI(:,:)
    nuclA(1:5,ichain)=nucl(1:5)
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
      write(ifile+1,'(100000E13.5)')MisfitA(ichain),VRA(ichain),nuclA(1:5,ichain),T0A(:,:,ichain),aA(:,:,ichain),baA(:,:,ichain),psiA(:,:,ichain),f0A(:,:,ichain),fwA(:,:,ichain),DcA(:,:,ichain),vwA(:,:,ichain),viniA(:,:,ichain),M0A(ichain),EgA(ichain),ErA(ichain),TshiftA(ichain),VRgpsA(ichain)
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
