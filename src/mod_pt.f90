module mod_ctrl !control parameters
!character(80) :: dir
!integer :: nproc, rank
!constrols:0
!integer :: ialg,nsteps,iburn,nbins,iseed0
!real :: swaprate,tlow,thigh

integer ::       ialg     = 0   ! Algorithm type: 
                                ! ialg=0, Parallel Tempering with T swap between all levels 
                                ! ialg=2, Parallel Tempering with T swap only allowed between 
                                ! neighbouring temperature levels 

real ::       swaprate = 1.0d0        ! Rate at which exchange swaps are proposed relative to within
                              ! chain steps. 1.0 = one exchage swap proposed for every McMC step.
                              ! Set this value to zero to turn off Parallel Tempering altogether.     
!real :: swaprate = 0.0d0       ! Turn off parallel tempering to allow performance comparison

integer ::      nsteps   = 500000         ! Number of chain steps per temperature
integer ::      iburn    = 0         ! Number of burn in samples
integer ::      nchains = 8            ! Define number of chains per processor
double precision ::       tlow  = 1.d0            ! Lowest temperature of chains
double precision ::      thigh = 100.d0           ! Highest temperature of chains
                               
!      iseed0 = 61254557       ! Random number seed
character(80) ::      dir = './'              ! Home directory for I/O files (some systems want full path names under MPI)
integer ::      nbins = 10              ! Number of temperature bins for diagnostics

integer :: mode,iseed0,nproc,rank
double precision, allocatable :: modtemp(:)   ! Temperatures of each chain
integer :: ifile,ierr,mrank
character(50) :: fname,rname,dname
integer :: iprint=10
integer :: ntemps

end module mod_ctrl


    
    module mod_pt

    implicit none
      integer                            :: iseed,nprocc,iprocc
      Double Precision, allocatable      :: Tbins(:)
      Real, allocatable                  :: sjd(:,:)
      Integer, allocatable               :: temptrans(:,:,:)
      Integer, allocatable               :: mcmctrans(:,:)
      Integer                            :: ntemps
      Integer                            :: iout
      Logical                            :: record_temp
      Logical                            :: record_mcmc
      Logical                            :: record_temp_now
      Logical                            :: record_mcmc_now
      Logical                            :: RestrictTemp
      Logical                            :: verbose
      Logical                            :: silent
      Real                               :: TimeStart,TimeEnd,time1,time2

contains
!-------------------------------------------------------------------------
!
!     Parallel Tempering routine 
!
!-------------------------------------------------------------------------
!
      Subroutine pt &
                 (mode,ialg,nchains,nsteps,iburn,modtemp,thigh,tlow,nbins,swaprate,&
                 &iseed0,basedir,nproc,iproc,iprint)
      IMPLICIT NONE
#if defined MPI
      include "mpif.h"
      Integer, dimension(MPI_STATUS_SIZE)           :: status
#endif

      Double precision, allocatable                 :: logPPD(:)
      Integer, allocatable                          :: finished(:)
      Double precision, dimension(4)                :: dmsg
      Double precision, dimension(3)                :: rmsg
      Double precision                              :: E1,E2,E
      Double precision                              :: t1,t2,T
      Double precision                              :: modtemp(*)
      Double precision                              :: thigh,tlow
      Integer, dimension(2)                         :: pair
      Integer                                       :: ipA,ipB,ipX,from
      Integer                                       :: tag,code,term
      Integer                                       :: totswapreq
      Integer                                       :: totswapacc
      Integer                                       :: totswapreqw
      Integer                                       :: totswapaccw
      Integer                                       :: nbins
      Logical, allocatable                          :: swapped(:)
      Logical                                       :: yes
      Character(len=100)                            :: filename
      Character(len=80)                             :: basedir
      Character*4                                   :: cproc

      integer :: mode,ialg,nchains,nsteps,iburn,iseed0,nproc,iproc
      real :: swaprate
      real :: a,cput,factormpi,factorT
      integer :: ic1,ic2,ichain,ierr,is,its,jchain,iprint
      real :: swapproc,swapratelocal,waittime,b,percentwait
      real:: ran3
      integer i
    
      iout=933
      
      if(mode.eq.0)then   ! Initialize routine
!				start MPI
#if defined MPI
         call MPI_Init(ierr )
         call MPI_Comm_size( MPI_COMM_WORLD, nproc, ierr )
         call MPI_Comm_rank( MPI_COMM_WORLD, iproc, ierr )
         write(*,*)'Hallo from proc. # ',iproc,'/',nproc
#else
         nproc = 1
         iproc = 0
#endif
         iprocc = iproc
         nprocc = nproc
!
! log file for each process -> easiest to see what's going on in MPI

         verbose = .false.
         silent = .false.
         write(cproc,'(I0)') iproc
         filename=trim(basedir)//'log/'//'log_'//trim(cproc)//'.txt'
         if(.not.silent)&
         open(unit=iout,file=trim(filename),action='write',status='replace') 
         ntemps = 0            ! initialize number of temperature bins to zero
         record_temp = .false. ! initialize diagnostic switches (to be reset by PT_diagnostics if called)
         record_mcmc = .false. ! initialize diagnostic switches (to be reset by PT_diagnostics if called)
         return

      else if(mode.eq.99)then   ! Finalize routine
         if(.not.silent)&
         close(iout)            ! close log files
#if defined MPI
         call MPI_BARRIER( MPI_COMM_WORLD, ierr)
         call MPI_Finalize(ierr)
         return
#endif
      else   ! Do the main work 

         waittime = 0.0
         CALL cpu_time(TimeStart)

         if(iprocc.eq.0)then                      ! write out pt log file
            filename=trim(basedir)//'pt.log'
            open(unit=123,file=filename,action='write',status='replace') 
            write(123,*)
            write(123,*)' Parallel Tempering log file'
            write(123,*)
            write(123,*)' Length of burn-in                   :',iburn
            write(123,*)' Length of each chain (post burn in) :',nsteps
            write(123,*)' Number of chains per processor      :',nchains
            write(123,*)' Number of processors                :',nprocc
            write(123,*)' Low temperature                     :',tlow
            write(123,*)' High temperature                    :',thigh
            write(123,*)' Number of temperature bins (diag)   :',nbins
            write(123,*)' Probability of Tempering swap       :',swaprate
            write(123,*)' Random number seed                  :',iseed0
            write(123,*)
            close(123)
         end if

         iseed = -abs(iseed0+iprocc*iprocc*1000) ! re-initialize random number generator
         a = ran3(iseed)
         record_temp_now = .false.
         record_mcmc_now = .false.

         allocate(logPPD(nchains))
         allocate(swapped(nchains))
         allocate(finished(nprocc)) 

         if(iproc>0.or.nproc==1)then
           do i=1,nchains
               call InitiateChain(i,logPPD(i),iseed)
           enddo         
         endif
         
!         write(*,*)logPPD(1:3)
         
         totswapreq=0
         totswapacc=0
         totswapreqw=0
         totswapaccw=0

         RestrictTemp = .false.
         if(ialg.eq.2)RestrictTemp = .true.

         ! The specified exchange swap rate is provided per chain and per step
         ! but exchange swaps occur after each chain has been advanced once and
         ! and hence swaprate need to be adjusted. 

                                            ! If we are restricting T-swaps to nearest neighbour Tbin 
                                            ! then temperature swap rate needs to be adjusted further
                                            ! (see notes)
         factorT = 1.0
         if(RestrictTemp.and.ntemps.gt.1)then  
           factorT = ntemps*ntemps/(2.0*(ntemps-1))    ! nearest neighbour Tswaps only
         else if(ntemps.gt.1)then
           factorT = ntemps/(ntemps-1)                 ! all neighbour Tswaps
         end if
         swapratelocal = factorT*swaprate

                                            ! Adjust temperature swap rate to account for
                                            ! parallel case (see notes)
         factormpi = 1.0
         if(iprocc.ne.0)factormpi = (2.0*(nprocc-1)-1)/(nprocc-1)
         swapratelocal = factormpi*swapratelocal
         its = 0

         if ( iprocc .eq. 0 .and. nprocc .gt. 1) then

                                            ! Master makes temperature swap decisions 
                                            ! for case where multiple processes are invoked
#if defined MPI
            finished=0
            finished(1:nprocc-1)=1
            !nbtotswap=ns*(nprocc-1)
            pair=0
            do while ( sum(finished) .ne. 0) 
               from = MPI_ANY_SOURCE
               tag = MPI_ANY_TAG
               call MPI_RECV(dmsg,4,MPI_DOUBLE_PRECISION,from,tag,MPI_COMM_WORLD,status,code)

               ! dmsg = double precision message sent by process who wants to swap
               ! dmsg(1) : iprocc
               ! dmsg(2) : 0 or 1 -> termination tag, 0 means finished
               ! dmsg(3) : temperature
               ! dmsg(4) : logPPD
               ! dmsg(2) equals 0 means that the process has finished its job

               if ( dmsg(2) .eq. 0.0 ) then
                    finished(int(dmsg(1)))=0
                    if(verbose)write(iout,*)&
                    'proc ',int(dmsg(1)),&
                    ' has sent its termination tag',int(dmsg(2))

               elseif ( (sum(pair) .eq. 0) .or. (pair(1) .eq. dmsg(1))) then

                    ! form the first process of a pair and prevent swapping 

                    pair(1)=dmsg(1)
                    t1=dmsg(3) ! store temperature for pair 1
                    E1=dmsg(4) ! store logPPD for pair 1
                    if(verbose)write(iout,*)&
                    'proc ',dmsg(1),' waiting for partnership...'
               else 

                    ! form the second process of a pair 

                    pair(2)=dmsg(1)
                    rmsg(2)=1 ! swap depending on result ...
                    rmsg(1)=pair(2)

                    t2=dmsg(3) !store temperature for pair 2
                    E2=dmsg(4) !store logPPD for pair 2

                    ! Master decides if the pair will swap or not based 
                    ! on temperature and logPPD of both processes
                    totswapreq=totswapreq+1

                    call tswap_accept(t1,t2,E1,E2,yes)  ! decide on T swap

                    if (yes) then
                       rmsg(2)=1
                       totswapacc=totswapacc+1
                    else 
                       rmsg(2)=0
                    endif

                    ! rmsg = double precision message sent by Master to both processes
                    ! rmsg(1) : process number to swap with
                    ! rmsg(2) : 0 or 1 -> swap tag, 1 means YES to swap
                    ! rmsg(3) : temperature value to swap

                    if(verbose)write(iout,*)&

                    'pair',pair(1),pair(2),'will swap at temperature',t1,t2
                    rmsg(3)=t2

                    call MPI_SEND &
                         (rmsg,3,MPI_DOUBLE_PRECISION,pair(1),&
                          2002,MPI_COMM_WORLD,code)

                    rmsg(1)=pair(1)
                    rmsg(3)=t1

                    call MPI_SEND &
                         (rmsg,3,MPI_DOUBLE_PRECISION,pair(2),&
                          2002,MPI_COMM_WORLD,code)

                    if(verbose)write(iout,*) ''
                    pair=0
               endif
                                                ! There is one process left 
                                                ! and no one can form a pair with it 

               if( (sum(finished) .eq. 1) .and. ( (sum(pair) .ne. 0)) ) then
                       rmsg(1)=pair(1)
                       rmsg(2)=0


                       call MPI_SEND &
                            (rmsg,3,MPI_DOUBLE_PRECISION,pair(1),&
                             2002,MPI_COMM_WORLD,code)

                       !print *,'MASTER says to ',dmsg(1),' not to swap and carry on'
               endif
            enddo
            if(.not.silent)then
            write(iout,*) '********************* MASTER HAS FINISHED ', &
                          '**************************'
            write(iout,*)'__________________________________________'
            write(iout,*)'Total number of swaps requested by slaves   :',&
                         totswapreq
            write(iout,*)'Total number of swaps accepted by master    :',&
                     totswapacc
            end if
#endif
         else                  ! Slave does its work 
!
!                              ! initialize default message 
            dmsg(1)=iprocc      ! dmsg(1) = process id to send to MASTER if swapping
            dmsg(2)=1.0        ! dmsg(2) = 1 :wants to swap or 0 :job finished 

            totswapreq=0
            totswapacc=0
            totswapreqw=0
            totswapaccw=0

            swapproc = 1.0     ! Set probability ratio for within and between node proposal
            if(nprocc.gt.1)swapproc = 1.0/(2.0*(nprocc-1)-1.0)
            if(nchains.eq.1)swapproc = -1.0

            do is = 1,iburn+nsteps  ! Main loop over chain steps        

!              if(is.gt.iburn)call printlog(is-iburn)

               if(is.gt.iburn)record_temp_now = record_temp   ! turn on recording of diagnostics after burnin
               if(is.gt.iburn .and. (mod(is,iprint)==0))record_mcmc_now = .true.   ! turn on recording of diagnostics after burnin

!		        Advance all chains over current step

               do ichain = 1,nchains    ! loop over chains for current step
!
                                        ! advance of current chain 

                  T = modtemp(ichain)   ! get temperature of current chain

                  ! Advance the chain using the user supplied algorithm
                  ! (This could be McMC or any model space sampling/optimization routine)
!			print *,ichain,is
                  E = logPPD(ichain)   ! F.G.
                  
                  call AdvanceChain(ichain,T,E,record_mcmc_now,iseed)  

                  logPPD(ichain) = E    ! update new energy for current chain
               end do

               if (record_mcmc_now) CALL saveforrestart() !Save for restart

               
               record_mcmc_now=.false.

!		        Temperature swapping
!
                                                     ! do nothing unless we have more than one chain

               !if(is.gt.iburn.and.(nchains.gt.1.or.nprocc.gt.2))then ! start PT only after burn in

               if(nchains.gt.1.or.nprocc.gt.2)then  ! Turn PT on from start

                  swapped = .false.

                  do jchain=1,nchains


                     a = ran3(iseed)
                     if(a.le.swapratelocal)then ! do we perform a temperature swap?
                        a = ran3(iseed) 
                                              ! do we swap within or between processors?
                        if(a.le.swapproc)then ! swap within current proessor

                           ichain = jchain
                           ic1 = 1 + nchains*ran3(iseed)
                           ic1 = min(nchains,ic1)
                           ic1 = max(1,ic1)
                           ic2 = 1 + nchains*ran3(iseed)
                           ic2 = min(nchains,ic2)
                           ic2 = max(1,ic2)

                           t1 = modtemp(ic1)
                           t2 = modtemp(ic2)
                           !if(.not.swapped(ic1).and..not.swapped(ic2))then
                              totswapreqw = totswapreqw + 1
                              E1 = logPPD(ic1)
                              E2 = logPPD(ic2)

                              call tswap_accept(t1,t2,E1,E2,yes)

!                                        accept the swap between chains

                              if(yes)then
                                 !swap models ic1 and ic2 in array model
                                 !call swapmodels(ic1,ic2)
                                 !swap logPPD values between ic1 and ic2

                                 !logPPD(ic1) = E2       !changed by MS 16/6/2014
                                 !logPPD(ic2) = E1       !changed by MS 16/6/2014

                                 modtemp(ic1) = t2
                                 modtemp(ic2) = t1
                                 swapped(ic1) = .true.
                                 swapped(ic2) = .true.
                                 if(verbose)write(iout,*)&
                                 'Swapmodels within node: chains',&
                                  ic1,ic2,'at temps',t1,' and ',t2
                                  totswapaccw = totswapaccw + 1
                              end if
                           !end if
 
                        else          ! swap with another processor
#if defined MPI
                                      ! The number of chain steps at each temperature
                                      ! is statistically identical between nodes. Exact 
                                      ! equivalence between nodes is not guaranteed. 

                           ichain = 1 + nchains*ran3(iseed)
                           ichain = min(nchains,ichain)
                           ichain = max(1,ichain)

                                      ! we send current temperature 
                                      ! and logPPD to master and request swap

                           dmsg(3)= modtemp(ichain)
                           dmsg(4)=logPPD(ichain) 

                                      ! send swap request to master 

                           CALL cpu_time(Time1)
                           call MPI_SEND &
                                (dmsg,4,MPI_DOUBLE_PRECISION,&
                                 0,2001,MPI_COMM_WORLD,code)

                                      ! receives answer from master

                           call MPI_RECV &
                                 (rmsg,3,MPI_DOUBLE_PRECISION,&
                                 0,2002,MPI_COMM_WORLD,status,code)
                           CALL cpu_time(Time2)
                           waittime = waittime + Time2-Time1

                           totswapreq=totswapreq+1

                           if (rmsg(2) .eq. 1) then

                                      ! master says we have a new temperature (=rmsg(3))
                                      ! obtained from processor rmsg(1)

                           if(verbose)&
                           write(iout,*)' Proc ',iprocc,' chain ',ichain,&
                           ' with temperature',modtemp(ichain),&
                           ' has new temperature',rmsg(3),'from ',rmsg(1)
                           modtemp(ichain)=rmsg(3)
                           totswapacc=totswapacc+1

                           else
                                      ! master says no to swap so carry on

                           if(verbose)&
                           write(iout,*) ' CARRY ON, no swap !'
                           endif
#endif
                        end if
                     end if
                  end do                  ! end loop over chains
!
              endif
!
!             end loop over chain steps

          end do

          if(.not.silent)then
          write(iout,*)'proc', iprocc,' has terminated'
          write(iout,*)'__________________________________________'
          write(iout,*)'Total number of swaps requested within node   :',&
                        totswapreqw
          write(iout,*)'Total number of swaps accepted  within node   :',&
                        totswapaccw
          write(iout,*)'Total number of swaps requests to master      :',&
                       totswapreq
          write(iout,*)'Total number of swaps accepted by master      :',&
                     totswapacc
          end if
          !write(iout,*)'Total number of nonswap        :',totswapreq-totswapacc
          dmsg(1)=iprocc
          dmsg(2)=0.0
          dmsg(3)=0.0
          dmsg(4)=0.0
          ! send to master the end of swapping phase (msg(2)=0)
#if defined MPI
          call MPI_SEND(dmsg,4,MPI_DOUBLE_PRECISION,0,2006,MPI_COMM_WORLD,code)
#endif
          !       end case processes
        endif
        CALL cpu_time(TimeEnd)

        cput = TimeEnd - TimeStart
        percentwait = 100.0*waittime/cput
        if(.not.silent)then
        write(iout,*)
        write(iout,*)' Time taken by the PT routine on this node  was ',cput,'seconds'
        write(iout,*)' Percentage time waiting for partner ',percentwait,'%'
        write(iout,*)
        end if

#if defined MPI
        call MPI_BARRIER( MPI_COMM_WORLD, ierr)
        !write(*,*)' iprocc',iprocc,' after barrier end of pt ierr',ierr
#endif


      end if

      return
      end subroutine 

                 
! -----------------------------------------------------------
!
!       tswap_accept -> decides whether to accept a proposed temperature swap 
!                       using the Parallel Tempering extension of 
!                       Metropolis-Hastings acceptance criterion.
!
! 	tf(4) 	= T1 -> temperature 1
!		= E1 -> Energy 1 
!		= T2 -> temperature 2
!		= E2 -> Energy 2
!
!	yn	= .true. or .false.
!
!       Notes:
!             In an optimization problem Energy, e.g. E1,or E2 should
!             be the property of the model to be minimized.
!
!             In a posterior PDF sampling problem E should be
!             the negative log of the posterior PDF plus the 
!             log of the proposal distribution from state 1 to state 2.
!
!             i.e. E2 = -log[p(m2|d)] + log[q(m2|m1)]
!
!             The proposal distribution must correspond to that specified
!             in used to perturb model 1 to model 2.
!
!             E need only be evaluated up to an additive constant
!	      because the decision to accept or reject this swap
!             depends only on differences in E.
!           
!             If the prior is a constant and the proposal 
!             distribution is symmetric then E may be set to 
!             to the negative log-Likelihood, or data misfit. 
!
!             This routine optionally records the history of 
!             successful temperature swap for diagnostic purposes.
!
! ---------------------------------------------------------------
!
       Subroutine tswap_accept(T1,T2,E1,E2,yn)
       IMPLICIT NONE
       Double precision              :: E1,E2,delE
       Double precision              :: T1,T2,delT,delS
       Logical                       :: yn
! from pt global
      Integer, allocatable               :: temptrans(:,:,:)
      Real, allocatable                  :: sjd(:,:)
      real :: a
      integer :: it,jt,k1,k2
      Logical                            :: record_temp_now
      Logical                            :: RestrictTemp
      real :: ran3

 
       it = 0
       jt = 0
       yn = .false.
       delE = E1-E2
       delT = 1.0/T1 - 1.0/T2
       delS = delE*delT
       a = ran3(iseed)

       if(log(a).le.delS)yn = .true.      ! swap successful

                                          ! find temperature bins

!       if(RestrictTemp.or.record_temp_now)call PT_findTbin(T1,k1)
!       if(RestrictTemp.or.record_temp_now)call PT_findTbin(T2,k2)

!       if(RestrictTemp.and.abs(k1-k2).ne.1)yn = .false. ! restrict temperature swaps to neighbours
!       if(record_temp_now)then  ! record successful temperature swaps for diagnostics
!          if(yn)then
!             temptrans(k1,k2,1) = temptrans(k1,k2,1) + 1  ! record success
!             temptrans(k2,k1,1) = temptrans(k1,k2,1)
!             sjd(k1,k2) = sjd(k1,k2) + real(delT*delT)
!             sjd(k2,k1) = sjd(k1,k2)
!          end if
!          temptrans(k1,k2,2) = temptrans(k1,k2,2) + 1  ! record attempt
!          temptrans(k2,k1,2) = temptrans(k1,k2,2)
!       endif

       end subroutine
!
! ---------------------------------------------------------------
!
!      PT_findTbin - A utility routine to locate the bin containing a 
!                    given temperature, T.
!
! ---------------------------------------------------------------
!
       Subroutine PT_findTbin(T,kbin)


       Double precision              :: T
       Integer                       :: kbin
      Double Precision, allocatable      :: Tbins(:)
       integer :: i,ntemps

       kbin = ntemps
       do i=ntemps-1,1,-1
          if(T.le.Tbins(i))then
            kbin = i
          end if
       end do

       return
       end subroutine 
       

end module 
