Program CrGPSGf  

implicit none

integer, parameter :: ngmax=10000
real, parameter :: pi=3.141592654
real df,aw1,t0
complex rseis,ui,freq
real,allocatable,dimension(:):: strike,dip,rake,leng,widt,hhypo
real,allocatable:: hypo(:,:),offsets(:,:)
real gleng,gwidt
real TM(3,3),ITM(3,3)
real NEZhypo(3),hypo2(3),xi(3),sour(3)
real dx1,dx2
real,allocatable:: x1a(:),x2a(:)
real, allocatable:: stalocGPS(:,:),dataGPS(:,:)
integer,allocatable:: ng1(:),ng2(:)
integer NSeg,np
integer i,j,k,kk, IRET
integer nc,nfreq,nrGPS,ns,ikmax
real Stat(2)
real tl,aw,xl,uconv,fref
real lambda, mu, dL, dW, dumGPS, temp, ALPHA
real gpsgfN,gpsgfE,gpsgfZ
character*6 filename



  open(1,file='input.dat')

  read(1,*)
  read(1,*) nfreq
  read(1,*)
  read(1,*) tl
  read(1,*)
  read(1,*) t0,NSeg
  allocate(ng1(NSeg),ng2(NSeg),strike(NSeg),dip(NSeg),rake(NSeg),hhypo(NSeg),leng(NSeg),widt(NSeg),hypo(3,NSeg))
  allocate(offsets(2,NSeg))
  read(1,*)
  read(1,*) temp, nrGPS
  read(1,*)
  read(1,*) (ng1(kk),ng2(kk),kk=1,NSeg)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*) (strike(kk),dip(kk),rake(kk),kk=1,NSeg)
  read(1,*)
  read(1,*) (hhypo(kk),kk=1,NSeg)
  read(1,*)
  read(1,*) (leng(kk),widt(kk),kk=1,NSeg)
  read(1,*)
  read(1,*) (hypo(2,kk),hypo(1,kk),kk=1,NSeg)   !WARNING! IN PREVIOUS VERSION THE ORDER WAS hypo(1),hypo(2)!
  read(1,*)
  read(1,*) np
  
  close(1)
  
  lambda =1.
  mu=1.
  ALPHA=1.40047868354247!(lambda+mu)/(lambda+2.*mu)
  dL=leng(1)/(ng1(1))
  dW=widt(1)/(ng2(1))
  
  allocate(stalocGPS(3,nrGPS),dataGPS(3,NRgps))
  open(1,file='stations-GPS.dat')
  
  
  do k=1, nrGPS
  
    read(1,*) stalocGPS(1,k),stalocGPS(2,k)
    stalocGPS(3,k)=0.
    print*, stalocGPS(1,k),stalocGPS(2,k),stalocGPS(3,k)
  enddo
  close(1)
  stalocGPS=stalocGPS*1.e3
  
  open(2,file='NEZsorGPS.dat') 
  kk=1
  print*, ng2(kk), ng1(kk), DL, DW
  do k=1,NRgps
    open(1,file='sources.dat')
    do j=1,ng2(kk)
      do i=1,ng1(kk)
        read(1,*) filename,sour(1),sour(2),sour(3),strike(1),dip(1),rake(1)
		sour=sour*1000
		!print*, sour
        dumGPS=1.
        

              CALL DC3Dmodif(ALPHA,sour(1),sour(2),sour(3),strike(1),dip(1),rake(1),dL,dW,dumGPS,stalocGPS(1,k),stalocGPS(2,k),stalocGPS(3,k),gpsgfN,gpsgfE,gpsgfZ,IRET)
			  !print*, ALPHA,sour(1),sour(2),sour(3),strike(1),dip(1),rake(1),dL,dW,dumGPS,stalocGPS(1,k),stalocGPS(2,k),stalocGPS(3,k)
			  !print*, gpsgfN,gpsgfE,gpsgfZ,IRET
			  !pause
              write(2,*) gpsgfN,gpsgfE,gpsgfZ
			  !G(Nseis+(k-1)*3+1,SegShift+(j-1)*Ssvd*NL(kk)+(i-1)*Ssvd+1:SegShift+(j-1)*Ssvd*NL(kk)+i*Ssvd)=gpsgfN*dt/sigmaGPS(1,k)
              !G(Nseis+(k-1)*3+2,SegShift+(j-1)*Ssvd*NL(kk)+(i-1)*Ssvd+1:SegShift+(j-1)*Ssvd*NL(kk)+i*Ssvd)=gpsgfE*dt/sigmaGPS(2,k)
              !G(Nseis+(k-1)*3+3,SegShift+(j-1)*Ssvd*NL(kk)+(i-1)*Ssvd+1:SegShift+(j-1)*Ssvd*NL(kk)+i*Ssvd)=gpsgfZ*dt/sigmaGPS(3,k)

      enddo
    enddo
	close(1)
  enddo

  close(1)
  close(2)
  
  
  
 deallocate(dataGPS)

	
	
	
	
end program

