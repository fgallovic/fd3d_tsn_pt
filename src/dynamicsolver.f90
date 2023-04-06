! Main body of the program - initialization, 
! reading input files and running the inversion. 
!-------------------------------------------------------
! Authors: Frantisek Gallovic, Jan Premus, Lubica Valentova (2019)
! Charles University in Prague, Faculty of Mathematics and Physics

! This code is published under the GNU General Public License. To any
! licensee is given permission to modify the work, as well as to copy
! and redistribute the work or any derivative version. Still we would
! like to kindly ask you to acknowledge the authors and don't remove
! their names from the code. This code is distributed in the hope
! that it will be useful, but WITHOUT ANY WARRANTY.
! ------------------------------------------------------   
    PROGRAM fd3d_tsn_pt
    USE source_com
    USE friction_com
    USE fd3dparam_com
    USE waveforms_com
    USE inversion_com
    USE mod_pt
    USE mod_ctrl
    USE pml_com
    USE SlipRates_com
   	USE PostSeismic_com
#if defined FVW
    USE RATESTATE, only:slipOUT
#endif
#if defined GPUMPI
    USE openacc
#endif
    IMPLICIT NONE
#if defined MPI
    include 'mpif.h'
#endif
    integer i,j,k,kk,it,itot,ncpu
    real dum
    real,allocatable:: msfts(:,:)
    real,allocatable,dimension (:,:):: sliprateoutX,sliprateoutZ
    integer ito,ifrom,jto,jfrom
    logical err

    allocate(modtemp(nchains))
    call pt (0,0,0,0,0,modtemp,0.d0,0.d0,0,0.0,0,dir,nproc,rank,iprint)   !Musi to byt prvni prikaz

#if defined MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD,mrank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#else
    mrank=0
#endif
    
    iseed = -(1000+mrank*mrank*1000) ! re-initialize random number generator
    call alloc_temp(iseed)

#if defined GPUMPI
    i=acc_get_num_devices(acc_device_nvidia)
    if(i>1.and.mrank>0)then    !assuming two GPUs at the node; selecting GPU according to even/odd mrank
      write(*,*)'  Node ',mrank,' has ',i,' GPUs; running on GPU #',mod(mrank,i),'.'
      call acc_set_device_num(mod(mrank,i), acc_device_nvidia)
    endif
#endif

    call fd3d_init()  !Reads FD parameters
    call inversion_init()  !Reads GFs, observed waveforms

    if(RUNI==1.or.RUNI==2)then
#if defined MPI
      if(mrank>0 .or. ncpu==1)then
        ifile=1111+mrank
        write(fname,'(a,i3.3)') 'samples',mrank
        open(unit=ifile,file=trim(fname),iostat=ierr,STATUS='REPLACE')
        write(fname,'(a,i3.3)') 'sampls',mrank
        open(unit=ifile+2,file=trim(fname),iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
        write(rname,'(a,i3.3)') 'restart',mrank
        if (iwaveform==2) then
         write(dname,'(a,i3.3)') 'SA', mrank
         open(unit=ifile*10,file=trim(dname),iostat=ierr,access='STREAM',STATUS='REPLACE')
        endif
      endif
#else
      ifile=1111+mrank
      write(fname,'(a,i3.3)') 'samples',mrank
      open(unit=ifile,file=trim(fname),iostat=ierr,STATUS='REPLACE')
      write(fname,'(a,i3.3)') 'sampls',mrank
      open(unit=ifile+2,file=trim(fname),iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
      write(rname,'(a,i3.3)') 'restart',mrank
      if (iwaveform==2) then
        write(dname,'(a,i3.3)') 'SA', mrank
        open(unit=ifile*10,file=trim(dname),iostat=ierr,access='STREAM',STATUS='REPLACE')
      endif
#endif
    endif
    
    if(RUNI==0.or.RUNI==10)then !    Call the forward modelling only
      write(*,*)'Running forward modeling:'

#if defined TPV5
		CALL forwardspecialTPV5()
#endif
#if defined TPV8
		CALL forwardspecialTPV8()
#endif
#if defined TPV9
		CALL forwardspecialTPV9()
#endif
#if defined TPV104
		CALL forwardspecialTPV104()
#endif

      CALL readinversionresult()
      ioutput=1
      if(iwaveform==0)write(*,*)'Note: No seismogram calculation.'
#if !defined FVW
      if(RUNI==10)then
        CALL randomdynmod(nxt,nzt,dh,striniZ,peak_xz,Dc)
      endif
#endif
      if(RUNI==0)then
        write(*,*)'Note: Saving normal stress profile.'
        open(719,FILE='normalstressprofile.dat')
        do i=nabc+1,nzt-nfs
          write(719,*)dh*real(nzt-nfs-i)/1.e3,normstress(i)/1.e6
        enddo
        close(719)
      endif
      write(*,*)'Running dynamic rupture simulation...'

      CALL validatefd3dmodel(err)  

      if(err)then
        write(*,*)'Forward model violates constraints!'
      else
        write(*,*)'Forward model satisfies constraints!'
      endif
      call fd3d()
      print *,'-----------------------------------------------------'
      print *,'Ruptured area:       ',output_param(3)
      print *,'Average speed:       ',output_param(1)
      print *,'Seismic moment:      ',output_param(2)
      print *,'Average stress drop: ',output_param(4)
      print *,'Fracture energy:     ',output_param(5)
      print *,'Radiated energy:     ',output_param(6)
      print *,'Energy efficiency:   ',output_param(6)/(output_param(5)+output_param(6))
      print *,'-----------------------------------------------------'
      call syntseis()
      if (igps.eq.1) call CalcSyntGPS()

#if defined FVW
      if(iwaveform==1)then
        call evalmisfit()
        write(*,*)'Variance reduction: ',VR,' for shift',Tshift,'s'
        if(igps==1)write(*,*)'GPS Variance reduction: ',VRgps
        call plotseis()
        open(594,FILE='forwardmodelsamples.dat',iostat=ierr)
        write(594,'(10000000E13.5)')misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),M0,Eg,Er,Tshift,VRgps
        close(594)
        open(594,FILE='forwardmodelsampls.dat',iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
#if defined DIPSLIP
        write(594)misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),slipOUT(1,:),slipOUT(2,:),M0,Eg,Er,Tshift,VRgps
#else
        write(594)misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeX(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),slipOUT(1,:),slipOUT(2,:),M0,Eg,Er,Tshift,VRgps
#endif
        close(594)
      elseif (iwaveform==2) then
        call evalmisfit2()
        open(594,FILE='forwardmodelsamples.dat',iostat=ierr)
        write(594,'(10000000E13.5)')misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),M0,Eg,Er,Tshift
        close(594)
        open(594,FILE='forwardmodelsampls.dat',iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
#if defined DIPSLIP
        write(594)misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),slipOUT(1,:),slipOUT(2,:),M0,Eg,Er,Tshift,VRgps
#else
        write(594)misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeX(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),slipOUT(1,:),slipOUT(2,:),M0,Eg,Er,Tshift,VRgps
#endif
        close(594)
      elseif (iwaveform==3) then
        call evalmisfitM()
        open(594,FILE='forwardmodelsamples.dat',iostat=ierr)
        write(594,'(10000000E13.5)')misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),M0,Eg,Er,-1000.,output_param(1),output_param(4),-1.,output_param(3)
        close(594)
        open(594,FILE='forwardmodelsampls.dat',iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
#ifdef DIPSLIP
        write(594)misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),slipOUT(1,:),slipOUT(2,:),M0,Eg,Er,-1000.,output_param(1),output_param(4),-1.,output_param(3)
#else
        write(594)misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeX(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),slipOUT(1,:),slipOUT(2,:),M0,Eg,Er,-1000.,output_param(1),output_param(4),-1.,output_param(3)
#endif
        close(594)
      elseif(iwaveform==4.or.iwaveform==5)then
        if(iwaveform==4)then
          call evalmisfitSspec()
          write(*,*)'Variance reduction: ',VR
        else  !iwaveform==5
          call evalmisfitStime()
          write(*,*)'Variance reduction: ',VR,' for shift',Tshift,'s'
        endif
        open(594,FILE='forwardmodelsamples.dat',iostat=ierr)
        write(594,'(10000000E13.5)')misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),M0,Eg,Er,Tshift,output_param(1),output_param(4),-1.,output_param(3)
        close(594)
        open(594,FILE='forwardmodelsampls.dat',iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
#if defined DIPSLIP
        write(594)misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),slipOUT(1,:),slipOUT(2,:),M0,Eg,Er,Tshift,output_param(1),output_param(4),-1.,output_param(3)
#else
        write(594)misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeX(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),slipOUT(1,:),slipOUT(2,:),M0,Eg,Er,Tshift,output_param(1),output_param(4),-1.,output_param(3)
#endif
        close(594)
      else
        open(594,FILE='forwardmodelsamples.dat',iostat=ierr)
        write(594,'(10000000E13.5)')misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),M0,Eg,Er,-1.,output_param(1),output_param(4),-1.,output_param(3)
        close(594)
        open(594,FILE='forwardmodelsampls.dat',iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
#if defined DIPSLIP
        write(594)misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),slipOUT(1,:),slipOUT(2,:),M0,Eg,Er,-1000.,output_param(1),output_param(4),-1.
#else
        write(594)misfit,VR,nucl(1:5),T0I(:,:),aI(:,:),baI(:,:),psiI(:,:),f0I(:,:),fwI(:,:),DcI(:,:),vwI(:,:),viniI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeX(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),slipOUT(1,:),slipOUT(2,:),M0,Eg,Er,-1000.,output_param(1),output_param(4),-1.
#endif
        close(594)
      endif
      
#else

      if(iwaveform==1)then
        call evalmisfit()
        write(*,*)'Variance reduction: ',VR,' for shift',Tshift,'s'
        if(igps==1)write(*,*)'GPS Variance reduction: ',VRgps
        call plotseis()
        open(594,FILE='forwardmodelsamples.dat',iostat=ierr)
        write(594,'(10000000E13.5)')misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),M0,Eg,Er,Tshift,output_param(1),output_param(4),VRgps,output_param(3)
        close(594)
        open(594,FILE='forwardmodelsampls.dat',iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
#if defined DIPSLIP
        write(594)misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),M0,Eg,Er,Tshift,VRgps,output_param(3)
#else
        write(594)misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeX(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),M0,Eg,Er,Tshift,VRgps,output_param(3)
#endif
        close(594)
      elseif (iwaveform==2) then
        call evalmisfit2()
        open(594,FILE='forwardmodelsamples.dat',iostat=ierr)
        write(594,'(10000000E13.5)')misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),M0,Eg,Er,Tshift,output_param(1),output_param(4),VRgps,output_param(3)
        close(594)
        open(594,FILE='forwardmodelsampls.dat',iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
#if defined DIPSLIP
        write(594)misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),M0,Eg,Er,Tshift,output_param(1),output_param(4),VRgps,output_param(3)
#else
        write(594)misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeX(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),M0,Eg,Er,Tshift,output_param(1),output_param(4),VRgps,output_param(3)
#endif
        close(594)
      elseif (iwaveform==3) then
        call evalmisfitM()
        open(594,FILE='forwardmodelsamples.dat',iostat=ierr)
        write(594,'(10000000E13.5)')misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),M0,Eg,Er,-1000.,output_param(1),output_param(4),-1.,output_param(3)
        close(594)
        open(594,FILE='forwardmodelsampls.dat',iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM')
#ifdef DIPSLIP
        write(594)misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),M0,Eg,Er,-1000.,output_param(1),output_param(4),-1.,output_param(3)
#else
        write(594)misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeX(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),M0,Eg,Er,-1000.,output_param(1),output_param(4),-1.,output_param(3)
#endif
        close(594)
      elseif(iwaveform==4.or.iwaveform==5)then
        if(iwaveform==4)then
          call evalmisfitSspec()
          write(*,*)'Variance reduction: ',VR
        else  !iwaveform==5
          call evalmisfitStime()
          write(*,*)'Variance reduction: ',VR,' for shift',Tshift,'s'
        endif
        open(594,FILE='forwardmodelsamples.dat',iostat=ierr)
        write(594,'(10000000E13.5)')misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),M0,Eg,Er,Tshift,output_param(1),output_param(4),-1.,output_param(3)
        close(594)
        open(594,FILE='forwardmodelsampls.dat',iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
#if defined DIPSLIP
        write(594)misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),M0,Eg,Er,Tshift,output_param(1),output_param(4),-1.,output_param(3)
#else
        write(594)misfit,VR,T0I(:,:),TsI(:,:),DcI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeX(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),M0,Eg,Er,Tshift,output_param(1),output_param(4),-1.,output_param(3)
#endif
        close(594)
      else
        open(594,FILE='forwardmodelsamples.dat',iostat=ierr)
        write(594,'(10000000E13.5)')-1.,-1.,T0I(:,:),TsI(:,:),DcI(:,:),M0,Eg,Er,-1.,output_param(1),output_param(4),-1.,output_param(3)
        close(594)
        open(594,FILE='forwardmodelsampls.dat',iostat=ierr,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE')
#if defined DIPSLIP
        write(594)-1.,-1.,T0I(:,:),TsI(:,:),DcI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),M0,Eg,Er,-1000.,output_param(1),output_param(4),-1.
#else
        write(594)-1.,-1.,T0I(:,:),TsI(:,:),DcI(:,:),ruptime(nabc+1:nxt-nabc,nabc+1:nzt-nfs),slipX(nabc+1:nxt-nabc,nabc+1:nzt-nfs), &
          & rise(nabc+1:nxt-nabc,nabc+1:nzt-nfs),schangeX(nabc+1:nxt-nabc,nabc+1:nzt-nfs),MomentRate(:),M0,Eg,Er,-1000.,output_param(1),output_param(4),-1.
#endif
        close(594)
      endif

#endif
      
    elseif(RUNI==1)then  ! Inversion from an initial model
      write(*,*)'Running inversion:'
      ioutput=0
      if (iwaveform==0) then 
       print *,'warning: need waveforms for misfit calculations'
       stop
      endif
      call pt (1,ialg,nchains,nsteps,iburn,modtemp,thigh,tlow,nbins,swaprate,randseed,dir,nproc,rank,iprint)
      close(ifile)
      close(ifile+2)
      call pt (99,0,0,0,0,modtemp,0.d0,0.d0,0,0.0,0,dir,nproc,rank,iprint)

    elseif(RUNI==2)then  ! Restart inversion
      write(*,*)'Restarting inversion:'
      ioutput=0
      if (iwaveform==0) then 
       print *,'warning: need waveforms for misfit calculations'
       stop
      endif
      call pt (1,ialg,nchains,nsteps,iburn,modtemp,thigh,tlow,nbins,swaprate,randseed,dir,nproc,rank,iprint)
      close(ifile)
      close(ifile+2)
      call pt (99,0,0,0,0,modtemp,0.d0,0.d0,0,0.0,0,dir,nproc,rank,iprint)
    
    elseif(RUNI==5)then  ! Create mtilde.dat (for plotting) and/or calculate seismograms from the latest result
      write(*,*)'Reading slip rates...'
      OPEN(25,file='result/sliprateX.res',FORM='UNFORMATTED',ACCESS='STREAM')
      OPEN(26,file='result/sliprateZ.res',FORM='UNFORMATTED',ACCESS='STREAM')
      MomentRate=0.
      MSRX=0.
      MSRZ=0.
      allocate(sliprateoutX(nxt,nzt),sliprateoutZ(nxt,nzt))
      do it = 1,ntfd
        read(25,IOSTAT=i)sliprateoutX(nabc+1:nxt-nabc,nabc+1:nzt-nfs)
        if(i==-1)exit
        read(26)sliprateoutZ(nabc+1:nxt-nabc,nabc+1:nzt-nfs)
        k=int(real(it-1)*dt/dtseis)+1
        if(k<1.or.k>nSR)write(*,*)'CHYBA!',k,nSR
        do j=1,NW
          jto=max(1,int(dW/dh*j))+1+nabc
          jfrom=min(jto,int(dW/dh*(j-1))+1)+1+nabc
          do i=1,NL
            ifrom=int(dL/dh*(i-1))+1+nabc
            ito=int(dL/dh*i)+nabc
            kk=((j-1)*NL+i-1)*nSR+k
            MSRX(kk)=MSRX(kk)+sum(sliprateoutX(ifrom:ito,jfrom:jto))/dble((ito-ifrom+1)*(jto-jfrom+1)*(dtseis/dt))
            MSRZ(kk)=MSRZ(kk)+sum(sliprateoutZ(ifrom:ito,jfrom:jto))/dble((ito-ifrom+1)*(jto-jfrom+1)*(dtseis/dt))
          enddo
        enddo
        do j = nabc+1,nzt-nfs
          do i = nabc+1,nxt-nabc
            MomentRate(k)=MomentRate(k)+sqrt(sliprateoutX(i,j)**2+sliprateoutZ(i,j)**2)*muSource(i,j)*dh*dh/(dtseis/dt)
          enddo
        enddo
      enddo
      close(25)
      close(26)
      deallocate(sliprateoutX,sliprateoutZ)
      M0=sum(MomentRate(:))*dtseis
      Mw=(log10(M0)-9.1)/1.5
      ioutput=1
      CALL syntseis()
      if(iwaveform==1)then
        call evalmisfit()
        write(*,*)'Variance reduction: ',VR,' for shift',Tshift,'s'
        call plotseis()
      endif
    
    elseif(RUNI==10)then ! Calculate misfits for varied parameter (hard-coded)
      write(*,*)'Calculating misfits for varied parameter (hard-coded):'
      CALL readinversionresult()
      ioutput=1
      iwaveform=1
      dum=DCI(6,4)
      write(*,*)'Best-model value of the perturbed parameter: ',dum
      itot=11   !better be odd number
      allocate(msfts(itot*2,3))
      do i=1,itot*2
        DCI(6,4)=dum*(1.+real(i-1-(itot-1)/2)/real(itot-1)*1.)
        msfts(i,1)=DCI(6,4)
        CALL inversion_modeltofd3d()
        print *,'Running dynamic rupture simulation #',i,'/',itot,'...'
        call fd3d()
        print *,'-----------------------------------------------------'
        print *,'Average speed:       ',output_param(1)
        print *,'Seismic moment:      ',output_param(2)
        print *,'Surface of rupture:  ',output_param(3)
        print *,'Average stress drop: ',output_param(4)
        print *,'Fracture energy: ',output_param(5)
        print *,'Radiated energy:    ',output_param(6)
        print *,'Energy efficiency:               ',output_param(6)/(output_param(5)+output_param(6))
        print *,'-----------------------------------------------------'
        call syntseis()
        call evalmisfit()
        write(*,*)'Variance reduction: ',VR,' for shift',Tshift,'s'
        msfts(i,2)=VR
        msfts(i,3)=misfit
      enddo
      dum=minval(msfts(:,3))
      do i=1,itot*2
        write(333,*)msfts(i,1),msfts(i,2),msfts(i,3)-dum
      enddo

    endif

    END PROGRAM
