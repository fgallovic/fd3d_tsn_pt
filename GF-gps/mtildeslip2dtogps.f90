Program mtildeslip2dtogps

implicit none

real, parameter :: pi=3.141592654
real df,aw1,t0
complex rseis,ui,freq
real,allocatable,dimension(:):: leng,widt,strike,dip,rake,mtildeslip2d
real,allocatable:: offsets(:,:),sour(:,:)
real gleng,gwidt
real TM(3,3),ITM(3,3)
real NEZhypo(3),hypo2(3),xi(3)
real dx1,dx2
real,allocatable:: x1a(:),x2a(:)
real, allocatable:: stalocGPS(:,:)
integer,allocatable:: ng1(:),ng2(:)
integer NSeg,np
integer i,j,k,kk,l,IRET
integer nc,nfreq,nrGPS,ns,ikmax
real Stat(2)
real tl,aw,xl,uconv,fref
real lambda, mu, dL, dW, dumGPS, temp, ALPHA
real gpsgfN,gpsgfE,gpsgfZ,gpsN,gpsE,gpsZ
character*6 filename



  open(1,file='input.dat')

  read(1,*)
  read(1,*) nfreq
  read(1,*)
  read(1,*) tl
  read(1,*)
  read(1,*) t0,NSeg
  allocate(ng1(NSeg),ng2(NSeg),leng(NSeg),widt(NSeg))
  allocate(offsets(2,NSeg))
  read(1,*)
  read(1,*) temp, nrGPS
  read(1,*)
  read(1,*) (ng1(kk),ng2(kk),kk=1,NSeg)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*) 
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*) (leng(kk),widt(kk),kk=1,NSeg)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*) np
  
  close(1)
  
  lambda =1.
  mu=1.
  ALPHA=1.40047868354247!(lambda+mu)/(lambda+2.*mu)

  kk=1

  dL=leng(kk)/(ng1(kk))
  dW=widt(kk)/(ng2(kk))
  
!  allocate(stalocGPS(3,nrGPS))
!  open(1,file='stations-GPS.dat')
!  do k=1,nrGPS
!    read(1,*) stalocGPS(1,k),stalocGPS(2,k)
!    stalocGPS(3,k)=0.
!  enddo
!  close(1)

  nrGPS=50*50
  allocate(stalocGPS(3,nrGPS))
  l=0
  do j=1,50
    do i=1,50
      l=l+1
      stalocGPS(1,l)=(i-1)*1.-25.
      stalocGPS(2,l)=(j-1)*1.-25.
      stalocGPS(3,l)=0.
    enddo
  enddo

  stalocGPS=stalocGPS*1.e3
  
  allocate(mtildeslip2d(ng1(kk)*ng2(kk)))
  open(1,file='mtildeslip2D.dat')
  k=1
  do j=1,ng2(kk)
    read(1,*) mtildeslip2d(k:k+ng1(kk)-1)
    k=k+ng1(kk)
  enddo
  close(1)
  
  allocate(sour(3,ng1(kk)*ng2(kk)),strike(ng1(kk)*ng2(kk)),dip(ng1(kk)*ng2(kk)),rake(ng1(kk)*ng2(kk)))
  open(1,file='sources.dat')
  k=0
  do j=1,ng2(kk)
    do i=1,ng1(kk)
      k=k+1
      read(1,*) filename,sour(1:3,k),strike(k),dip(k),rake(k)
    enddo
  enddo
  close(1)
  sour=sour*1000

  open(2,file='mtildeslip2D-sGPS.out') 
  do k=1,NRgps
    gpsN=0.
    gpsE=0.
    gpsZ=0.
    l=0
    do j=1,ng2(kk)
      do i=1,ng1(kk)
        l=l+1
        dumGPS=1.
        CALL DC3Dmodif(ALPHA,sour(1,l),sour(2,l),sour(3,l),strike(l),dip(l),rake(l),dL,dW,dumGPS,stalocGPS(1,k),stalocGPS(2,k),stalocGPS(3,k),gpsgfN,gpsgfE,gpsgfZ,IRET)
        gpsN=gpsN+gpsgfN*mtildeslip2d(l)
        gpsE=gpsE+gpsgfE*mtildeslip2d(l)
        gpsZ=gpsZ+gpsgfZ*mtildeslip2d(l)
      enddo
    enddo
    write(2,*)stalocGPS(1,k)/1.e3,stalocGPS(2,k)/1.e3,gpsN,gpsE,gpsZ
  enddo
  close(2)
  
	
end program

