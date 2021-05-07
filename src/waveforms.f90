    MODULE waveforms_com
!waveforms for misfit
      integer,allocatable,dimension(:,:):: stainfo
      integer:: NRseis,NSTAcomp,Nseis
      real,allocatable,dimension(:):: Dsynt,stasigma, ruptdist
      real,allocatable,dimension(:,:):: Dseis, StaDist
      logical, allocatable :: srfmask(:) !for estimation Rjb for nga-west2
      real SigmaData   ! Typically 0.05m
      real misfit,VR,Tshift
      integer iwaveform,nper
      real,allocatable :: pgaM(:,:),pgaD(:,:),per(:)
      real, allocatable :: PGAsigma(:,:),PGAtau(:,:) ! now tau and sigma of gmpe are on output..


!GFs for synthetic seismograms
      real:: T,T1,T2,artifDT,M0aprior,Mwsigma,leng,widt,elem,df
      real,allocatable,dimension(:,:):: H
      real,allocatable,dimension(:):: fc1,fc2
      integer,allocatable,dimension(:):: fcsta
      integer:: nfmax,np,nfc,iT1,iT2,nT
      character(100) :: GFfile
    END MODULE

    
    SUBROUTINE readGFs()
    USE waveforms_com
    USE SlipRates_com
    USE fd3dparam_com
    USE medium_com
    use mod_pgamisf
    use pml_com, only : nabc
    use mod_ctrl, only : mrank
    IMPLICIT NONE
#if defined MPI
    include 'mpif.h'
#endif
    REAL,PARAMETER:: PI=3.1415926535
    REAL,DIMENSION(:),ALLOCATABLE:: fltr4
    COMPLEX,DIMENSION(:,:),ALLOCATABLE:: cirN,cirE,cirZ
    COMPLEX,DIMENSION(:),ALLOCATABLE:: cseis
    complex :: cdtseis
    real,allocatable,dimension(:,:):: staweight
    integer i,j,k,jw,jl,jj,mm,dumi,ierr
    real dum,xgf,ygf,zgf,xst,yst,zst
    real*4 dumarr(6)
    logical fileex

    write(*,*)'Reading GFs...'

    open(10,file='input.dat',action='read')
    read(10,*)
    read(10,*) nfmax
    read(10,*)
    read(10,*) T,dum,T1,T2
    read(10,*)
    read(10,*) artifDT
    read(10,*)
    read(10,*) NRseis
    read(10,*)
    read(10,*) NL,NW
    read(10,*)
    read(10,*) M0aprior
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
    read(10,*)
    read(10,*) !vr
    read(10,*)
    read(10,*) nfc   !number of frequency bands
    

    allocate(fc1(nfc),fc2(nfc))
    do i=1,nfc
      read(10,*) fc1(i),fc2(i)
    enddo
    close(10)
    dL=leng/dble(NL);dW=widt/dble(NW)
    elem=dL*dW
    dtseis=T/real(np)
    df=1./T
    iT1=T1/dtseis+1
    iT2=T2/dtseis+1
    nT=iT2-iT1+1
    nSR=ceiling(real(ntfd)/(dtseis/dt))
    allocate(MSRX(NL*NW*nSR),MSRZ(NL*NW*nSR),MomentRate(nSR))

    if(iwaveform==0)return
    
    open(10,file='stainfo.dat',action='read')
    allocate(stainfo(3,NRseis),staweight(3,NRseis),fcsta(NRseis))
    do i=1,NRseis
      read(10,*)stainfo(:,i),staweight(:,i),fcsta(i)
    enddo
    close(10)
    NSTAcomp=sum(stainfo(:,:))
    allocate(stasigma(NSTAcomp))
    Nseis=nT*NSTAcomp
   
    fileex=.false.
    if (iwaveform==1) allocate(H(Nseis,NL*NW*nSR))

    if (iwaveform==2) then   
       write(GFfile,'(a)') 'GFspectr_st'
       inquire(file=trim(GFfile),exist=fileex)
       if (mrank==0) open(222,file=trim(GFfile),form='unformatted')
    endif
 
   if (.not.fileex) then
    open(20,form='unformatted',file='NEZsor.dat')  !if problems with NEZsor.dat file appears, check that resort.f90 writes unformatted file (and not binary)!
    allocate(cirN(min(nfmax,np),NL*NW),cirE(min(nfmax,np),NL*NW),cirZ(min(nfmax,np),NL*NW))
    j=0
    do jj=1,NRseis
      i=0
      read(20) dumi
      do jw=1,NW
        do jl=1,NL
          i=i+1
          read(20) dumi
          do k=1,nfmax
            read(20) dumarr
            if(k>np)cycle
            dum=sum(muSource(int(dL/dh*(jl-1))+nabc+1:nabc+int(dL/dh*jl),int(dW/dh*(jw-1))+nabc+1:nabc+int(dW/dh*jw)))/(dL/dh*dW/dh)
            cirN(k,i)=cmplx(dumarr(1),dumarr(4))*dum*elem
            cirE(k,i)=cmplx(dumarr(2),dumarr(5))*dum*elem
            cirZ(k,i)=cmplx(dumarr(3),dumarr(6))*dum*elem
          enddo
        enddo
      enddo

     do k=1,3
        if(stainfo(k,jj)==0)cycle
        j=j+1
        stasigma(j)=staweight(k,jj)
!$OMP parallel private(i,mm,cseis,fltr4) DEFAULT(SHARED)
        allocate(cseis(np),fltr4(np))
!$OMP do SCHEDULE(DYNAMIC,1)
        do i=1,NL*NW
          cseis=0.
          SELECT CASE (k)
          CASE(1)
            cseis(1:min(nfmax,np))=cirN(1:min(nfmax,np),i)*df
          CASE(2)
            cseis(1:min(nfmax,np))=cirE(1:min(nfmax,np),i)*df
          CASE(3)
            cseis(1:min(nfmax,np))=cirZ(1:min(nfmax,np),i)*df
          END SELECT
          cseis(np/2+2:np)=conjg(cseis(np/2:2:-1))
!          cseis(np/2+1)=real(cseis(np/2+1))
          call four1(cseis,np,1)
          fltr4=real(cseis)!*dtseis   !dt for time integration moved
          do mm=1,int(artifDT/dtseis)  ! REMOVING DWN ARTIFACT FROM THE SEISMOGRAM BEGINNING
            fltr4(mm)=fltr4(mm)*(cos((dtseis*real(mm-1)-artifDT)/artifDT*PI)/2.+.5);
          enddo
          if(fc1(fcsta(jj))>0.)then
            CALL XAPIIR(fltr4, np, 'BU', 0.0, 0.0, 4,'BP', fc1(fcsta(jj)), fc2(fcsta(jj)), dtseis, 1, np)
          else
            CALL XAPIIR(fltr4, np, 'BU', 0.0, 0.0, 4,'LP', fc1(fcsta(jj)), fc2(fcsta(jj)), dtseis, 1, np)
          endif
          if (iwaveform==1) then
            do mm=2,np   !time integration
              fltr4(mm)=fltr4(mm)+fltr4(mm-1)
            enddo
            fltr4=fltr4*dtseis
          elseif (iwaveform==2) then
!            do mm=2,np   !time integration
!              fltr4(mm)=fltr4(mm)+fltr4(mm-1)
!            enddo
!            fltr4=fltr4*dtseis
            do mm=1,np-1   !time derivative
              fltr4(mm)=fltr4(mm+1)-fltr4(mm)
            enddo
            fltr4=fltr4/dtseis
          endif
          
          if (iwaveform==2) then
           cseis=fltr4*dtseis
           call four1(cseis,np,-1)
           if (mrank==0)   write(222) cseis(1:np/2+1)
          elseif (iwaveform==1) then 
           do mm=1,nT
            H((j-1)*nT+mm,(i-1)*nSR+1:i*nSR)=dble(fltr4(iT1-1+mm:iT1+mm-nSR:-1))
           enddo
          endif
        enddo
!$omp end do
        deallocate(cseis,fltr4)
!$omp end parallel
      enddo
    enddo
    deallocate(cirN,cirE,cirZ,staweight)
    close(20)    
   endif
   if (iwaveform==2 .and. mrank==0) then
     close(222)
     print *,'file with GFspectr saved' 
   endif
#if defined MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif    
    allocate(Dsynt(Nseis))
   
    if (iwaveform==2) then
     allocate(stadist(nl*nw,nrseis),ruptdist(nrseis))
     !read sources.dat
     open(224,file='sources.dat',status='old')
     open(225,file='stations.dat',status='old')
      do j=1,nrseis
       read(225,*) xst,yst,zst
       do i=1,nl*nw
        read(224,*) dum,xgf,ygf,zgf
        stadist(i,j)=sqrt((xst-xgf)**2+(yst-ygf)**2+(zst-zgf)**2)
       enddo
       rewind(224)
      enddo
     close(224)
     close(225)
     call pga_init()
     allocate(pgaM(nrseis,nper),pgaD(nrseis,nper))
     allocate(PGAsigma(nrseis,nper),PGAtau(nrseis,nper))
    endif 
    
    END
	
    SUBROUTINE syntseis()
    USE waveforms_com
    USE source_com
    USE fd3dparam_com
    USE medium_com
    USE SlipRates_com
    use mod_pgamisf, only : gmpe_id !we need it for estimatio of Rjb
    IMPLICIT NONE
    real,allocatable,dimension(:):: slipGF
    complex, allocatable, dimension(:,:) :: sr,cseis
    integer i,j,k,m, jj,ii,ierr,kk
    integer ifrom,ito,jfrom,jto,kfrom,kto
    real dum,maxslip,dum1,dum2
    COMPLEX,DIMENSION(:),ALLOCATABLE:: seis1
    logical, allocatable :: slipmask(:,:)

    allocate(slipGF(nl*nw),sr(np,nl*nw))
    slipGF=0.;sr=0.

    if (iwaveform==2) then
      m=0
      do j=1,NW
        do i=1,NL
          jj=(j-1)*NL+i
          do k=1,nSR
            m=m+1
#if defined DIPSLIP 
            sr(k,jj)=MSRZ(m)
            slipGF(jj)=slipGF(jj)+MSRZ(m)*dtseis
#else
            sr(k,jj)=MSRX(m)
            slipGF(jj)=slipGF(jj)+MSRX(m)*dtseis			
#endif			
          enddo
        enddo
      enddo
    endif
 
    if (ioutput.eq.1) then
      open(297,FILE='mtildeX.dat',iostat=ierr)
      write(297,'(1E13.5)')MSRX
      close(297)
      open(297,FILE='mtildeZ.dat',iostat=ierr)
      write(297,'(1E13.5)')MSRZ
      close(297)
      open(297,FILE='mtildemomentrate.dat')
      do k=1,nSR
        write(297,*)dtseis*(k-1),Momentrate(k)
      enddo
      close(297)
    endif

    if(M0<1.e14)then
       Dsynt=0.
    else
      if(iwaveform==1)then
#if defined DIPSLIP 
        Dsynt=matmul(H,MSRZ)*dtseis
#else	
        Dsynt=matmul(H,MSRX)*dtseis
#endif
      elseif(iwaveform==2) then
        allocate(cseis(np,nl*nw),seis1(np))
        maxslip=maxval(slipGF)
        do jj=1,nl*nw
          call four1(sr(:,jj),np,-1)
        enddo
        sr=sr*dtseis 
        m=0
        write(GFfile,'(a)') 'GFspectr_st'
        open(243,file=trim(GFfile),form='unformatted',status='old',iostat=ierr)
        if (ierr/=0) print *,'error while openening file:',trim(GFfile),ierr
        do jj=1,NRseis
         select case (GMPE_id)
         case(1)
          ruptdist(jj)=minval(stadist(:,jj),mask=(slipGF>0.1*maxslip))
         case(2)
            allocate(slipmask(nl,nw),srfmask(nl*nw))
            srfmask=.false.
            slipmask=.false.
            do i=1,NL
             do j=1,NW
              ii=(j-1)*NL+i
              if (slipGF(ii)>0.1*maxslip) slipmask(i,j)=.true.
             enddo
              if (any(slipmask(i,:))) srfmask((nw-1)*nl+i)=.true.
             enddo
            ruptdist(jj)=minval(stadist(:,jj),mask=(srfmask))
            deallocate(slipmask,srfmask)
          end select
          do k=1,3
            if(stainfo(k,jj)==0)cycle
            do i=1,NL*NW
             read(243) cseis(1:np/2+1,i)
            enddo
            do ii=1,np/2+1
              seis1(ii)=sum(sr(ii,:)*cseis(ii,:))
            enddo
            seis1(np/2+2:np)=conjg(seis1(np/2:2:-1))
            call four1(seis1,np,1)
            seis1=seis1*df
            do ii=it1,it2
              m=m+1
              Dsynt(m)=real(seis1(ii))
            enddo 
          enddo
        enddo
        close(243)
        deallocate(seis1,cseis)
      endif
    endif

    if((iwaveform==1 .or. iwaveform==2) .and. ioutput==1)then
      open(297,FILE='svseisnez.dat')
      do i=iT1,iT2
        write(297,'(1000E13.5)')dtseis*(i-1)-artifDT,(Dsynt((j-1)*nT+i-iT1+1),j=1,NSTAcomp)
      enddo
      close(297)
    endif
   
    deallocate(slipGF,sr)

    END    
    
    
    SUBROUTINE readwaveforms()
    USE waveforms_com
    IMPLICIT NONE
    real dum
    integer i,j,jj,k,kk,mm
    
    allocate(Dseis(nT,NSTAcomp))
    write(*,*)'Reading data...'
    j=0
    do jj=1,NRseis
      do k=1,3
        if(stainfo(k,jj)==0)cycle
        j=j+1
        SELECT CASE (k)
        CASE(1)
          open(290,FILE='rvseisn.dat')
        CASE(2)
          open(290,FILE='rvseise.dat')
        CASE(3)
          open(290,FILE='rvseisz.dat')
        END SELECT
        i=0
        do kk=1,iT2
          if(kk<iT1)then
            read(290,*)
          else
             i=i+1
             read(290,*)(dum,mm=1,jj),Dseis(i,j)
          endif
        enddo
        close(290)
      enddo
    enddo

    END
    

    subroutine evalmisfit()
    use waveforms_com
    use SlipRates_com
	use PostSeismic_com
	use ieee_arithmetic
	use source_com, only: ioutput
    
	implicit none
	
    real,parameter:: maxTshift=2.
    real dumn,dump,normdatn,normdatp,dum2,norma2,dum3,norma3,temp2
    integer i,k,ims
    
    ims=int(maxTshift/dtseis)
    misfit=1.e30
    VR=-1.e30
    Tshift=0.
    
    do i=0,ims-1
      dumn=0.;dump=0.
      normdatn=0.;normdatp=0.
      do k=1,NSTAcomp
        dumn=dumn+.5*sum((Dsynt((k-1)*nT+1:k*nT-ims)-Dseis(i+1:nT+i-ims,k))**2)*stasigma(k)**2/SigmaData**2
        normdatn=normdatn+.5*sum((Dseis(i+1:nT+i-ims,k))**2)*stasigma(k)**2/SigmaData**2
        dump=dump+.5*sum((Dsynt((k-1)*nT+i+1:k*nT+i-ims)-Dseis(1:nT-ims,k))**2)*stasigma(k)**2/SigmaData**2
        normdatp=normdatp+.5*sum((Dseis(1:nT-ims,k))**2)*stasigma(k)**2/SigmaData**2
      enddo
      if(dumn<misfit)then
        Tshift=-dtseis*i !negative shifts
        misfit=dumn
        VR=1.-dumn/normdatn
      endif
      if(dump<misfit)then
        Tshift=dtseis*i !positive shifts
        misfit=dump
        VR=1.-dump/normdatp
      endif
    enddo
	
	print*,'seis misfit: ', misfit
	
	if (NGPS>0) then
	
	   dum2=0.
	   norma2=0.
	   dum3=0.
	   norma3=0.
	   do i=1,NGPS
	   
	    dum3=dum3+0.5*((gpssyntN(1,i)-gpsrealN(1,i))**2)*SigmaGPS**2/gpssigma(1,i)**2
	    dum3=dum3+0.5*((gpssyntE(1,i)-gpsrealE(1,i))**2)*SigmaGPS**2/gpssigma(2,i)**2
	    dum3=dum3+0.5*((gpssyntZ(1,i)-gpsrealZ(1,i))**2)*SigmaGPS**2/gpssigma(3,i)**2
		
	    norma3=norma3+0.5*((gpsrealN(1,i))**2)*SigmaGPS**2/gpssigma(1,i)**2
	    norma3=norma3+0.5*((gpsrealE(1,i))**2)*SigmaGPS**2/gpssigma(2,i)**2
	    norma3=norma3+0.5*((gpsrealZ(1,i))**2)*SigmaGPS**2/gpssigma(3,i)**2

	    if (gpsrealTN(i)>1) then
		   temp2=gpssyntN(1,i)-gpsrealN(1,i)
		   dum2=dum2+0.5*sum((gpssyntN(2:gpsrealTN(i),i)-gpsrealN(2:gpsrealTN(i),i)-temp2)**2)*SigmaGPS2**2/gpssigma(1,i)**2
	       temp2=gpssyntE(1,i)-gpsrealE(1,i)
		   dum2=dum2+0.5*sum((gpssyntE(2:gpsrealTN(i),i)-gpsrealE(2:gpsrealTN(i),i)-temp2)**2)*SigmaGPS2**2/gpssigma(2,i)**2
	       temp2=gpssyntZ(1,i)-gpsrealZ(1,i)
		   dum2=dum2+0.5*sum((gpssyntZ(2:gpsrealTN(i),i)-gpsrealZ(2:gpsrealTN(i),i)-temp2)**2)*SigmaGPS2**2/gpssigma(3,i)**2
	    norma2=norma2+0.5*sum((gpsrealN(1:gpsrealTN(i),i)-gpsrealN(1,i))**2)*SigmaGPS2**2/gpssigma(1,i)**2
		norma2=norma2+0.5*sum((gpsrealE(1:gpsrealTN(i),i)-gpsrealE(1,i))**2)*SigmaGPS2**2/gpssigma(2,i)**2
		norma2=norma2+0.5*sum((gpsrealZ(1:gpsrealTN(i),i)-gpsrealZ(1,i))**2)*SigmaGPS2**2/gpssigma(3,i)**2
		endif
	   enddo

!	   if (ieee_is_nan(dum2)) then
!	     dum2=1.e30
!		 dum3=1.e30
!		 norma2=1.e1
!	   endif
	   print*,'GPS misfit: ', dum2, dum3, norma2, norma3
	   misfit=misfit + dum2 + dum3
           if (norma2>0.) then
	     VRGPS=1.-(dum2/norma2+dum3/norma3)
           else
	     VRGPS=1.-(dum3/norma3)
	   endif

	   if (ioutput.EQ.1) then
	        open(3141,file='result/misfit.dat')  
		    write(3141,'(2E13.5)') misfit - (dum2+dum3), (dum2+dum3)
			close(3141)
	   endif
    endif
	
    if(Mwsigma>0.)misfit=misfit+0.5*(2./3.*log10(M0/M0aprior)/Mwsigma)**2
    
    END

    subroutine evalmisfit2() !old type of misfit calculation
    use mod_pgamisf
    use waveforms_com
    use SlipRates_com
    use source_com, only : ioutput,output_param
    use fd3dparam_com
    use pml_com, only: nabc,nfs
    implicit none
    integer :: jj,k,m,i
    real :: diff,mean,misf
    logical :: k1,k2
    logical, allocatable :: kmask(:)
    integer :: nstat,ierr
    real :: ctrl1,ctrl2
    real, allocatable :: hor1(:),hor2(:),rtd50(:)

    misfit=1.e30
    misf=0.
    diff=0.
    nstat=0 
    pgaM=0.
    pgaD=0.
    pgasigma=0.
    pgatau=0.

    !do not calculate and discard model if:
    if (Mw<5.) return !if moment too small - seismograms are zero
    if (abs(output_param(3)-(nxt-2*nabc)*(nzt-nabc-nfs)*dh*dh)<0.5*dh*dh) return !if whole fault ruptured -->  discard model : rupture may have continued on larger fault
    if (MomentRate(nSr)>0.) return ! if momentrate function not finished --> discard model, as the rupture may have continued
    allocate(kmask(nrseis))
    kmask=.false.

     if (ioutput==1) open(2223,file='dobliky.dat')
     m=0
     allocate(hor1(nT),hor2(nT),rtd50(nper))
     hor1=0.;hor2=0;rtd50=0.
     do jj=1,NRseis
        k1=.false.
        k2=.false.
      do k=1,3
        if(stainfo(k,jj)==0)cycle
        m=m+1
        if (k==1) then 
          hor1=Dsynt((m-1)*nT+1:m*nT)*100. !in cm
          k1=.true.
        elseif (k==2) then
          hor2=Dsynt((m-1)*nT+1:m*nT)*100.
          k2=.true.
        endif
        if (k1*k2 .and. k<3) then
          select case (GMPE_id)
          case(1) !ZHAO
           call  pcn05(nT,nT,nper,nper,dtseis,damp,per(:),hor1,sd,sv,sa1,psv,psa)
           call  pcn05(nT,nT,nper,nper,dtseis,damp,per(:),hor2,sd,sv,sa2,psv,psa)
           do i=1,nper
              pgaD(jj,i)=sqrt(sa1(i)*sa2(i))
           enddo
          case(2) !BOORE
            call get_rotd50(hor1,hor2,nT,nT,dtseis,per,nper,damp,rtd50)
            pgaD(jj,:)=rtd50(:)
          end select
          kmask(jj)=.true.
          nstat=nstat+1
        endif
       enddo
     enddo
     deallocate(hor1,hor2,rtd50)
     
     do i=1,nper
      do jj=1,nrseis
       if (kmask(jj)) then
          call pga_theor(ruptdist(jj),Mw,per(i),pgaM(jj,i),PGAsigma(jj,i),PGAtau(jj,i))
          misf=misf+(log(pgaD(jj,i))-pgam(jj,i))**2/(PGAsigma(jj,i)**2) !first part difference from the mean value for each station for particular event
          if (ioutput==1) write(2223,*) Mw,ruptdist(jj),per(i),exp(pgaM(jj,i))/100.,pgaD(jj,i)/100.,PGAsigma(jj,i),PGAtau(jj,i)
       endif
      enddo
      if (ioutput==1) write(2223,*)
      !control whether sigma and tau are not distance/magnitude dependent, otherwise how to misfit??!!
      ctrl1=sum(pgasigma(:,i))/nrseis
      ctrl2=sum(pgatau(:,i))/nrseis
      if ((any(abs(pgasigma(:,i)-ctrl1)>eps)) .or. (any(abs(pgatau(:,i)-ctrl2)>eps))) then
         print *,'sigma problem:', abs(pgasigma(:,i)-ctrl1) 
         print *,'tau problem:', abs(pgatau(:,i)-ctrl2) 
      else
        mean=sum(log(pgad(:,i)),mask=kmask)/nstat
        diff=sum(pgam(:,i),mask=kmask)/nstat 
        misf=misf-(mean-diff)**2*nstat**2*PGAtau(1,i)**2/(PGAsigma(1,i)**2+nstat*PGAtau(1,i)**2)/PGAsigma(1,i)**2
      endif
     enddo
    misf=misf/2./nper
    if (ioutput==1) close(2223)
!    deallocate(pgaD,pgaM,kmask)
    print *,'misfit:',misf
    misfit=misf 
    end subroutine 


    
    SUBROUTINE plotseis()
    USE waveforms_com
    USE SlipRates_com
    IMPLICIT NONE
    REAL, PARAMETER:: margin=0.05
    CHARACTER*5,ALLOCATABLE,DIMENSION(:):: staname
    REAL startx,starty,stept
    REAL,ALLOCATABLE,DIMENSION(:):: stepa,maxampl
    real dum
    INTEGER i,j,k,jj

    allocate(stepa(NRseis),maxampl(NRseis),staname(NRseis))
    open(233,FILE='stainfo.dat')
    do k=1,NRseis
      read(233,*)(dum,i=1,7),staname(k)
    enddo
    close(233)

    stept=1./dble(3)*(1.-margin)/(nT-1)
    jj=0
    do j=1,NRseis
      maxampl(j)=0.
      if(stainfo(1,j)==1)then
        jj=jj+1
        maxampl(j)=max(maxampl(j),maxval(abs(dseis(:,jj))))
      endif
      if(stainfo(2,j)==1)then
        jj=jj+1
        maxampl(j)=max(maxampl(j),maxval(abs(dseis(:,jj))))
      endif
      if(stainfo(3,j)==1)then
        jj=jj+1
        maxampl(j)=max(maxampl(j),maxval(abs(dseis(:,jj))))
      endif
      stepa(j)=.5/dble(NRseis)*(1.-margin)/maxampl(j)
    enddo

    open(205,FILE='seisplot.real.dat')
    jj=0
    do j=1,NRseis
      starty=((NRseis-j+1)+margin/2.)/dble(NRseis+1)
      if(stainfo(1,j)==1)then
        jj=jj+1
        startx=((1-1)+margin/2.)/3.
        do k=0,nT-1
          write(205,'(3E13.5)')startx+stept*(k),starty+stepa(j)*Dseis(k+1,jj)
        enddo
        write(205,*)
      endif
      if(stainfo(2,j)==1)then
        jj=jj+1
        startx=((2-1)+margin/2.)/3.
        do k=0,NT-1
          write(205,'(3E13.5)')startx+stept*(k),starty+stepa(j)*Dseis(k+1,jj)
        enddo
        write(205,*)
      endif
      if(stainfo(3,j)==1)then
        jj=jj+1
        startx=((3-1)+margin/2.)/3.
        do k=0,NT-1
          write(205,'(3E13.5)')startx+stept*(k),starty+stepa(j)*Dseis(k+1,jj)
        enddo
        write(205,*)
      endif
    enddo
    close(205)
    open(205,FILE='seisplot.synt.dat')
    jj=0
    do j=1,NRseis
      starty=((NRseis-j+1)+margin/2.)/dble(NRseis+1)
      if(stainfo(1,j)==1)then
        jj=jj+1
        startx=((1-1)+margin/2.)/3.
        do k=0,nT-1
          write(205,'(3E13.5)')startx+stept*(k-int(Tshift/dtseis)),starty+stepa(j)*Dsynt(k+1+(jj-1)*nT)
        enddo
        write(205,*)
      endif
      if(stainfo(2,j)==1)then
        jj=jj+1
        startx=((2-1)+margin/2.)/3.
        do k=0,NT-1
          write(205,'(3E13.5)')startx+stept*(k-int(Tshift/dtseis)),starty+stepa(j)*Dsynt(k+1+(jj-1)*nT)
        enddo
        write(205,*)
      endif
      if(stainfo(3,j)==1)then
        jj=jj+1
        startx=((3-1)+margin/2.)/3.
        do k=0,NT-1
          write(205,'(3E13.5)')startx+stept*(k-int(Tshift/dtseis)),starty+stepa(j)*Dsynt(k+1+(jj-1)*nT)
        enddo
        write(205,*)
      endif
    enddo
    close(205)

    open(201,FILE='seisplot.gp')
    write(201,*)'set term postscript portrait color solid enh'
    write(201,*)'set output "seisplot.ps"'
    write(201,*)'set xrange [0:1]'
    write(201,*)'set yrange [0:1]'
    write(201,*)'unset xtics'
    write(201,*)'unset ytics'
    write(201,*)'set border 0'
    do j=1,NRseis
      write(201,'(A11,F5.2,A23,G,A6)')'set label "',real(maxampl(j))*100.,'" at graph 1.09, graph ',((NRseis-j+1)+margin/2.)/dble(NRseis+1),' right'
      write(201,*)'set label "'//staname(j)//'" at graph -0.01, graph ',((NRseis-j+1)+margin/2.)/dble(NRseis+1),' right'
    enddo
    write(201,*)'set label "N-S" at .15,1'
    write(201,*)'set label "E-W" at .48,1'
    write(201,*)'set label "Z" at .86,1'
    write(201,*)'plot "seisplot.real.dat" u 1:2 notitle w l lt -1 lw 1,\'
    write(201,*)'"seisplot.synt.dat" u 1:2 notitle w l lt 1 lw 2'
    close(201)
    
    END
    
    
    
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
