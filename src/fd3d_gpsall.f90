!-------------------------------------------------------
! Forward modeling of static (GNSS) displacements for fd3d_tsn_pt
!-------------------------------------------------------
! Authors: Jan Premus (10/2023)
! Charles University in Prague, Faculty of Mathematics and Physics

! This code is published under the GNU General Public License. To any
! licensee is given permission to modify the work, as well as to copy
! and redistribute the work or any derivative version. Still we would
! like to kindly ask you to acknowledge the authors and don't remove
! their names from the code. This code is distributed in the hope
! that it will be useful, but WITHOUT ANY WARRANTY.
! ------------------------------------------------------

    MODULE source_com
      REAL,ALLOCATABLE,DIMENSION(:,:):: ruptime,rise,slipZ,schangeX,schangeZ,sliptime,slipX
      real    :: output_param(6)
      integer :: ioutput
      integer:: Nstations
      integer,allocatable:: staX(:), staY(:), staZ(:)
      REAL,allocatable :: seisU(:), seisV(:), seisW(:)
	  real :: waveT
	  real :: Eg, Er
	 ! REAL,allocatable :: waveU(:,:,:),waveV(:,:,:), waveW(:,:,:)
	  
    END MODULE

   MODULE SlipRates_com
      INTEGER nSR,NL,NW
      REAL dL,dW,dtseis
      REAL M0,Mw
      REAL,allocatable,dimension(:):: MSRX,MSRZ,MomentRate
      REAL,allocatable,dimension(:,:):: muSource
   END MODULE
	
   MODULE PostSeismic_com
	
      REAL,allocatable,dimension(:):: MSX, MSZ, TS
      REAL,allocatable,dimension(:,:):: gpssyntN, gpssyntE, gpssyntZ
      REAL,allocatable,dimension(:,:):: gpsrealN, gpsrealE, gpsrealZ, gpsrealT
      REAL,allocatable,dimension(:,:):: gpsgfN, gpsgfE, gpsgfZ, gpssigma
      integer,allocatable,dimension(:)::gpsrealTN
      REAL :: T1S,T2S,SigmaGPS,SigmaGPS2, gpsweight, VRGPS
      integer :: NTS, NGPS, igps, NTSrv
	  
   END MODULE

PROGRAM fd3d_gpsall

    use source_com
    use SlipRates_com
    use PostSeismic_com

    implicit none
    
    double precision T, tw0, artifDT, dum, Tshift
    double precision dip, hypodepth, leng, widt, epicL,epicW
    double precision,ALLOCATABLE,DIMENSION(:):: fc1,fc2, MSXtemp, MSZtemp, sl
    double precision, DIMENSION(:,:),ALLOCATABLE:: sr, stalocGPS
    integer nfmax, NR, np, np_source, ifc
    integer i, jw, jl, j, k, kk

    double precision dt, df, elem

    integer npRIK


    write(*,*)'Reading parameters...'

    open(10,file='input.dat',action='read')
    read(10,*)
    read(10,*) nfmax
    read(10,*)
    read(10,*) T,tw0
    read(10,*)
    read(10,*) artifDT,dum,Tshift
    read(10,*)
    read(10,*) NR, NGPS
    read(10,*)
    read(10,*) NL,NW
    read(10,*)
    read(10,*) 
    read(10,*)
    read(10,*) dip,dip
    read(10,*)
    read(10,*) hypodepth
    read(10,*)
    read(10,*) leng,widt
    read(10,*)
    read(10,*) epicL,epicW
    read(10,*)
    read(10,*) np, np_source
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*) ifc   !number of frequency bands
    allocate (fc1(ifc),fc2(ifc))
    do i=1,ifc
      read(10,*) fc1(i),fc2(i)
    enddo
    close(10)
    
    dt=T/dble(np)
    npRIK=int(tw0/DT+.999d0)
    write(*,*)npRIK, NL, NW, npRIK*NL*NW

    df=1.d0/T
    dL=leng/dble(NL)
    dW=widt/dble(NW)
    elem=dL*dW

    write(*,*)'Reading slip rates...'
    ALLOCATE(sr(np,NL*NW))
    ALLOCATE(sl(NL*NW))
    
    !Read slip rates and calculate slip
    sr=0.d0
    sl=0.d0
    open(201,FILE='mtildeX.dat')
    i=0
    do jw=1,NW
      do jl=1,NL
        i=(jw-1)*NL+jl
        do j=1,np_source  !41
          read(201,*)sr(j,i)
	  sl(i)=sl(i)+sr(j,i)*dt
        enddo
       ! sr(:,i)=sr(:,i)*mu(jl,jw)*elem
      enddo
    enddo
    close(201)

    allocate(MSXtemp(NL*NW),MSZtemp(NL*NW))
    allocate(gpsgfN(NL*NW,NGPS),gpsgfE(NL*NW,NGPS),gpsgfZ(NL*NW,NGPS))
    allocate(gpssyntN(1,NGPS),gpssyntE(1,NGPS),gpssyntZ(1,NGPS))  
    allocate(stalocGPS(2,NGPS))

    open(10,file='NEZsorGPS.dat',action='read')
	
    do i=1,NGPS
        kk=0
        do j=1,NW
            do k=1,NL
                kk=kk+1
                read(10,*) gpsgfN(kk,i),gpsgfE(kk,i),gpsgfZ(kk,i)
            enddo
        enddo
    enddo
    close (10)

    !calculate GPS displacements

    open(10,file='GPSdist.dat',action='write')    	
    do i=1,NGPS	
	    
        MSXtemp(1:NL*NW)=sl(1:NL*NW)
       ! MSZtemp(1:NL*NW)=sl(1:NL*NW)
        gpssyntN(1,i)=dot_product(MSXtemp(1:NL*NW),gpsgfN(1:NL*NW,i))
	gpssyntE(1,i)=dot_product(MSXtemp(1:NL*NW),gpsgfE(1:NL*NW,i))
	gpssyntZ(1,i)=dot_product(MSXtemp(1:NL*NW),gpsgfZ(1:NL*NW,i))	
        write (10,'(3E13.5)') gpssyntN(1,i), gpssyntE(1,i), gpssyntZ(1,i)			
	print*, gpssyntN(1,i), gpssyntE(1,i), gpssyntZ(1,i)
    enddo
    
    close(10)	

    open(1,file='stations-GPS.dat')
    open(10,file='GPSarrows.dat')

    do i=1, NGPS
        read(1,*) stalocGPS(2,i), stalocGPS(1,i)
        print*, stalocGPS(2,i), stalocGPS(1,i), gpssyntN(1,i), gpssyntE(1,i)
        write (10,'(4E13.5)') stalocGPS(2,i), stalocGPS(1,i), gpssyntN(1,i)*100, gpssyntE(1,i)*100
    enddo
    close(1)
    close(10)


END PROGRAM 
