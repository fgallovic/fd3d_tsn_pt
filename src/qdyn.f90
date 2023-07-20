#if defined FVW
   MODULE RATESTATE
      INTEGER MDIS,NDIS,MNDIS,MNDIS2,fft(2),MDISFT,NDISFT
      REAL*8 DX,DY,VS, G2B,Gtemp
      REAL*8 dtRATE,dtRUPT,dtSLIP
      REAL*8,ALLOCATABLE:: stresstep(:),tangstep(:),rupttime(:),timedcff(:)
      INTEGER,ALLOCATABLE:: ruptloc(:)
      REAL*8,ALLOCATABLE:: slip(:),ruptsnapshot(:,:),slipOUT(:,:),VPL
	  REAL*8,ALLOCATABLE:: f0PS(:,:),baPS(:,:),fwPS(:,:),vwPS(:,:),aPS(:,:),dcPS(:,:),snPS(:,:)
      COMPLEX*16,ALLOCATABLE:: Kij(:),varsc(:)
      integer gkoef
      double precision dhh
    END MODULE
	
    Subroutine QDYN3D
    
    use RATESTATE
    use fd3dparam_com
    use medium_com
    use friction_com
    use source_com
    use outputs_com
    use PostSeismic_com

    ! FFTW add-on
    !use mFFT_NR
    use mFFT_FFTW

    implicit none
    INTEGER,PARAMETER:: KMAXX=20000
    INTEGER kmax,kount
    REAL*8 dxsav,hstar,saveittime
    COMMON /path/ kmax,kount,dxsav
    real*8 time,time0,time1,time2,dum,LC,Z0,eps,h1,hmin,I1m,I1p,I2m,I2p,I1pI2p,I1mI2m,I1mI2p,I1pI2m,koeftemp
    integer nok,nbad,ruptongoing,scenarios,ifrom,ito,saveit
    real*8,allocatable:: vars(:),rupts(:)
	real*8,allocatable:: sliptot(:,:)
    integer,allocatable:: locs(:)
    integer :: ierr, rank, sizerank
	integer :: px, pz
    integer i,j,k,istatus,nsaves, i2, i3, j2, j3, ni, i4, j4

    ! FFTW add-on
    real cpu_time1,cpu_time2
    call cpu_time(cpu_time1)
    print *,'init_FFTW in QDYN3D'
    call init_FFTW

	!gkoef je v fd3d_init
	dhh=dh*gkoef  
	
	do i=1,100
	  if (MDIS<=2**i) then
	    MDISFT=2**i
	    exit
	  endif
	enddo
	
	do i=1,100
	  if (NDIS<=2**i) then
	    NDISFT=2**i
	    exit
	  endif
	enddo
	print*, MDISFT, NDISFT
    MNDIS=MDIS*NDIS;MNDIS2=MDIS*NDIS*2
    scenarios=1
    allocate(Kij(MDISFT*NDISFT*8))
    allocate(slip(MNDIS),varsc(MDISFT*NDISFT*8))
    allocate(vars(MNDIS2))
    allocate(stresstep(scenarios),tangstep(scenarios),timedcff(scenarios))
    allocate(ruptsnapshot(scenarios,MNDIS))
    timedcff(1)=0.
	stresstep(1)=0.
	tangstep(1)=0.
	slipOUT=0.
    Kij=0.
	
    fft(1)=NDISFT*4;fft(2)=MDISFT*2
    VPL=5.e-8

    DX=dhh
    DY=dhh
    LC=2000.d0*2.d0
    z0=10000.d0
    write(*,*)MDIS,'x',NDIS

	!open(102,FILE='g2b.txt')
  
    VS=v0
    k=0
    do j=1,NDIS
      do i=1,MDIS
        k=k+1
		i3=nabc+1+(i-1)*gkoef
		j3=nabc+1+(j-1)*gkoef
		slip(k)=slipX(i3,j3)
		slipOUT(1,k)=slip(k)
		f0PS(i,j)=f0Z(i3,j3)
		baPS(i,j)=baZ(i3,j3)
		fwPS(i,j)=fwZ(i3,j3)
		vwPS(i,j)=vwZ(i3,j3)
		aPS(i,j)=aZ(i3,j3)
		dcPS(i,j)=DCZ(i3,j3)
		snPS(i,j)=Sn(i3,j3)
      enddo
    enddo
	px=nxt/2
	pz=nzt/2
	Gtemp=mu1(px,nyt,pz)
	G2B=Gtemp/2.d0/sqrt(Gtemp/d1(px,nyt,pz))
	koeftemp=2*(lam1(px,nyt,pz)+Gtemp)/(lam1(px,nyt,pz)+2*Gtemp)
	!print*, sqrt(mu1(px,nyt,pz)/d1(px,nyt,pz)), d1(px,nyt,pz), SnZ(px,pz)
	!print*, koeftemp
	
  !close(102)
    hstar=2.d0*Gtemp/pi/SnZ(nxt/2,nzt/2)/maxval(baZ/DcZ(nxt/2,nzt/2))
	write(*,*)'G=  ',Gtemp
    write(*,*)'h*=  ',hstar
    write(*,*)'h/h*=',max(DX,DY)/hstar
    write(*,*)'hmax=',2.d0*Gtemp/pi/SnZ(nxt/2,nzt/2)/maxval(baZ)*maxval(DcZ)

    do i=1,MDISFT*2
      do j=1,NDISFT*4
        if(i>MDISFT)then
          I1m=(dble(MDISFT*2-i+1)-.5d0)*DX
          I1p=(dble(MDISFT*2-i+1)+.5d0)*DX
        else
          I1m=(dble(i-1)-.5d0)*DX
          I1p=(dble(i-1)+.5d0)*DX
        endif
        if(j>NDISFT*2)then
          I2m=(dble(NDISFT*4-j+1)-.5d0)*DY
          I2p=(dble(NDISFT*4-j+1)+.5d0)*DY
        else
          I2m=(dble(j-1)-.5d0)*DY
          I2p=(dble(j-1)+.5d0)*DY
        endif
        I1mI2p=sqrt(I1m**2+I2p**2)
        I1mI2m=sqrt(I1m**2+I2m**2)
        I1pI2p=sqrt(I1p**2+I2p**2)
        I1pI2m=sqrt(I1p**2+I2m**2)
		
		! if ((j<NDIS) .OR. (j>3*NDIS)) then
		  ! Gtemp=G(1)
		  ! koeftemp=koef(1)
		! else
		
		  ! if ((i<MDIS/2+1) .OR. (i>MDIS+MDIS/2+1)) then
		    ! Gtemp=G(1)  !?
		    ! koeftemp=koef(1)
		    
		  ! elseif (i>=MDIS/2+1) then
		    ! Gtemp=G((j-1-NDIS-1)*MDIS+i-MDIS/2)
		    ! koeftemp=koef((j-1-NDIS-1)*MDIS+i-MDIS/2)
		  ! else
			! Gtemp=G(MDIS*2*NDIS*4-((j-1-NDIS-1)*MDIS+i-MDIS/2))
			! koeftemp=koef(MDIS*2*NDIS*4-((j-1-NDIS-1)*MDIS+i-MDIS/2))
		  ! endif
		
		! endif
		Kij((i-1)*NDISFT*4+j)=Gtemp/4.d0/pi*(koeftemp*(1.d0/I1m*(I2p/I1mI2p-I2m/I1mI2m)-1.d0/I1p*(I2p/I1pI2p-I2m/I1pI2m))+1.d0/I2m*(I1p/I1pI2m-I1m/I1mI2m)-1.d0/I2p*(I1p/I1pI2p-I1m/I1mI2p))/dble(MDISFT*NDISFT*8)
      enddo
    enddo
	
    ! FFTW add-on
    call init_fourn(Kij,fft,2)
    call fourn(Kij,fft,2,-1)
    call destroy_fourn

    ifrom=1
    ito=1
    if (ioutput.eq.1) then	
	    open(114,FILE='schangep.txt')
        open(113,FILE='postmom.txt')
        open(112,FILE='slip.txt')
        open(111,FILE='state.txt')
        open(110,FILE='varsmax2.txt')
    endif
   ! open(104,FILE='friction.txt')
   ! open(114,FILE='slipslices.txt')
   ! open(115,form='binary',FILE='timefunc.dat')
   ! open(113,FILE='ruptslip.txt')
   ! open(117,FILE='ruptvelo.txt')
   ! open(120,FILE='epicsxy.dat')
   ! open(121,FILE='ruptsnapshot.dat')

    k=0
    do j=1,NDIS
      do i=1,MDIS
        k=k+1

        vars(k)=0.
		vars(k+MNDIS)=0.
		ni=0
		if (gkoef>1) then
          do j2=1,gkoef+1
          do i2=1,gkoef+1
		    i3=nabc+1+(i-1)*gkoef - gkoef/2 +i2
		    j3=nabc+1+(j-1)*gkoef - gkoef/2 +j2

		  !i3=nabc+i!+1
		  !j3=nabc+j!+1
            if ((i3>nabc+1).AND.(i3<nxt-nabc)) then	
            if ((j3>nabc+1).AND.(j3<nzt-nfs)) then		
		      ni=ni+1
		      vars(k)=vars(k)+abs(sqrt((sliprateoutX(i3,j3)+2*uini(i3,j3))**2))!+sliprateoutZ(i3,j3)**2))!vars(k)+sqrt(**2)!+sliprateoutZ(i3,j3)**2)
              vars(k+MNDIS)=vars(k+MNDIS)+psiout(i3,j3)
		    endif
		    endif
          enddo
          enddo
		  vars(k)=vars(k)/ni
		  vars(k+MNDIS)=vars(k+MNDIS)/ni
		else
		  vars(k)=abs(sqrt((sliprateoutX(nabc+i,nabc+j)+2*uini(nabc+i,nabc+j))**2+sliprateoutZ(nabc+i,nabc+j)**2))!vars(k)+sqrt(**2)!+sliprateoutZ(i3,j3)**2)
          vars(k+MNDIS)=psiout(nabc+i,nabc+j)		
		
		endif
      enddo
    !  dum=dy*(dble(j)-.5d0)
    !  if(dum>Z0-LC/2.d0.and.dum<Z0+LC/2.d0) vars(MDIS/4:MDIS/4*3+1,j)=10.d0*vars(MDIS/4:MDIS/4*3+1,j)
    enddo

    kmax=KMAXX
    eps=1.d-5
    h1=1.d-10
    hmin=1.d-16
    ruptongoing=0
    time0=T1S
    time2=T2S
    dxsav=(time2-time0)/dble(kmax-1)
    if (ioutput.eq.1) then
	open(102,FILE='vars.txt')
        write(102,'(E18.11)')vars
        close(102)	
    endif
    
    ! FFTW add-on
    call init_fourn(varsc,fft,2)
    CALL odeint(vars,time0,time2,time,eps,h1,hmin,nok,nbad,ruptongoing,1,1)
    call destroy_fourn

    do k=1,MNDIS
	    slipOUT(2,k)=slip(k)
	enddo
	
	!interpolating final (total) slip back to the FD grid
    if (ioutput.eq.1) then
	
	allocate(sliptot(nxt,nzt))
	sliptot=0.!slipX	!slightly dirty, we assign value of slip at the border
	k=0
    do j=1,NDIS
      do i=1,MDIS
        k=k+1
		i3=nabc+1+(i-1)*gkoef
		j3=nabc+1+(j-1)*gkoef
		!slip(k)=slipX(i3,j3)
		
		!if ((j.NE.1) .AND. (j.NE.NDIS) .AND. (i.NE.1) .AND. (i.NE.MDIS)) then
		
		    do i4=-int(gkoef/2),int(gkoef/2)
			    do j4=-int(gkoef/2),int(gkoef/2)
				
					sliptot(i3+i4,j3+j4)=slip(k) 
					
		        enddo
		    enddo
			
		!endif
      enddo
    enddo
	
	open(102,FILE='result/totalslip.txt')
	do k = nabc+1,nzt-nfs
        do i = nabc+1,nxt-nabc
            write(102,*) sliptot(i,k)
        enddo
    enddo
	close(102)

    ! FFTW add-on
    print *,'cleanup_FFTW in QDYN3D'
    call cleanup_FFTW
    call cpu_time(cpu_time2)
    print *,'CPU time of qdyn: ',cpu_time2-cpu_time1
    
	!open(102,FILE='vars.txt')
    !write(102,'(E18.11)')vars
    ! close(102)	
	    close(114)
        close(113)
        close(112)
		close(111)
        close(110)
    endif
    deallocate(Kij)
    deallocate(slip,varsc)
    deallocate(vars)
    deallocate(stresstep,tangstep,timedcff)
    deallocate(ruptsnapshot)

    END SUBROUTINE qdyn3d


    SUBROUTINE rates(t,vars,s,varsdx)
    USE RATESTATE	
	use fd3dparam_com
	use medium_com
	use friction_com

    ! FFTW add-on
    !use mFFT_NR
    use mFFT_FFTW
	
    IMPLICIT NONE
    REAL*8 t
    INTEGER scennum
    REAL*8 vars(MNDIS2),varsdx(MNDIS2)
    REAL*8 varsdum,psiss,fss,flv,dsdv,dsds,xx
    integer i,j,k,istatus,idum,s,i3,j3
    varsc=0.
!$OMP parallel do private(i,j,idum,varsdum) DEFAULT(SHARED) SCHEDULE(DYNAMIC,10)
    do i=1,MDIS
      idum=(i+MDISFT/2)*NDISFT*4-NDISFT*2
      do j=1,NDIS
        varsdum=vars(i+(j-1)*MDIS)-VPL
        varsc(idum-j+1)=varsdum
        varsc(idum+j)=varsdum
      enddo
    enddo
!$OMP end parallel do
    call fourn(varsc,fft,2,-1)
    varsc=varsc*Kij
    call fourn(varsc,fft,2,1)
  !  write(114,*) varsc(1:8*MNDIS)
	!write(114,*) ' '
!$OMP parallel do private(i,j,k,i3,j3,flv,fss,psiss,xx,dsdv,dsds) DEFAULT(SHARED) SCHEDULE(DYNAMIC,10)
    do j=1,NDIS
      do i=1,MDIS
            k=(j-1)*MDIS+i
	    !right side - evolution of state variable
	    i3=nabc+1+(i-1)*gkoef
        j3=nabc+1+(j-1)*gkoef	
	    flv = f0PS(i,j) - baPS(i,j)*log(abs(vars(k))/v0)
	    fss=fwPS(i,j) + (flv - fwPS(i,j))/((1. + (abs(vars(k))/vwPS(i,j))**8)**(1./8.))
        psiss = aPS(i,j)*(log(sinh(fss/aPS(i,j))) + log(2.*v0/abs(vars(k)))) 		
        varsdx(k+MNDIS)=-abs(vars(k))/dcPS(i,j)*(vars(MNDIS+k)-psiss)   !1.d0-vars(k)*vars(k+MNDIS)/DC(k)-ALPHAl(k)*vars(k+MNDIS)/Bl(k)*SIGD(i,j)/SIG(i,j)
	    !right side - velocity	
		!strsum=strsum+SnZ(i3,j3)*aII*asinh(xx)
		!frsum=frsum+f0Z(i3,j3)*SnZ(i3,j3)
		xx=abs(vars(k))/(2.*v0)*exp(vars(MNDIS+k)/aPS(i,j))
		if (vars(k)>0) then
	    !dsdv=aII*1./(abs(vel)/(2.*v0)*sqrt(1.+1./(xx*xx)))*(1.+xx/sqrt(xx*xx+1.))*1./(2.*v0)
		   dsdv=aPS(i,j)*1./sqrt(xx*xx+1)*xx/vars(k)
		else 
		   dsdv=-aPS(i,j)*1./sqrt(xx*xx+1)*xx/vars(k)
        endif
	    dsds=aPS(i,j)*1./sqrt(xx*xx+1)*xx/aPS(i,j)  !1./(abs(vel)/(2.*v0)*sqrt(1.+1./(xx*xx)))*(1.+xx/sqrt(xx*xx+1))*abs(vel)/(2.*v0)
 	    varsdx(k)=dble(varsc((i+MDISFT/2-1)*NDISFT*4+j+NDISFT*2)-SnPS(i,j)*dsds*varsdx(k+MNDIS))/(SnPS(i,j)*dsdv+G2B)!G2B(k))
      enddo
    enddo
!$OMP end parallel do
    END
	
	SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,scennum)
	IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION h,x,dydx(n),y(n),yerr(n),yout(n)
      INTEGER scennum
      DOUBLE PRECISION ak2(n),ak3(n),ak4(n),ak5(n),ak6(n),&
      ytemp(n),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,&
      B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,B31&
      =3.d0/40.d0,&
      B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,B52&
      =2.5d0,&
      B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0,B62=175.d0&
      /512.d0,&
      B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,C1&
      =37.d0/378.d0,&
      C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,DC1=C1-2825.d0&
      /27648.d0,&
      DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,DC5=-277.d0&
      /14336.d0,&
      DC6=C6-.25d0)
      ytemp=y+B21*h*dydx
      call rates(x+A2*h,ytemp,scennum,ak2)
      ytemp=y+h*(B31*dydx+B32*ak2)
      call rates(x+A3*h,ytemp,scennum,ak3)
      ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
      call rates(x+A4*h,ytemp,scennum,ak4)
      ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
      call rates(x+A5*h,ytemp,scennum,ak5)
      ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
      call rates(x+A6*h,ytemp,scennum,ak6)
      yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
      yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
      return
    END

    SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,scennum)
	IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      INTEGER scennum
      DOUBLE PRECISION errmax,h,htemp,xnew,yerr(n),ytemp(n),SAFETY
      DOUBLE PRECISION PGROW,PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89d-4)

      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,scennum)
      errmax=maxval(abs(yerr/yscal))
      errmax=errmax/eps
      if(errmax.gt.1.d0)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1d0*abs(h)),h)
        xnew=x+h
        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.d0*h
        endif
        hdid=h
        x=x+h
        y=ytemp
        return
      endif
    END

    
    SUBROUTINE odeint(ystart,x1,x2,x,eps,h1,hmin,nok,nbad,ruptongoing,scennum)
      USE RATESTATE
	  USE PostSeismic_com
	  USE SlipRates_com
	  use source_com, only: ioutput
	  use friction_com, only: v0
      IMPLICIT NONE
	  INTEGER nbad,nok,MAXSTP,ruptongoing,scennum
      DOUBLE PRECISION eps,h1,hmin,x1,x2,ystart(MNDIS2),TINY
      PARAMETER (MAXSTP=0,TINY=1.d-30)
      INTEGER kmax,kount,nstp
      DOUBLE PRECISION dxsav,h,hdid,hnext,x,xsav,dydx(MNDIS2)
	  double precision y(MNDIS2),dtout,yscal(MNDIS2)
      double precision slip0(MNDIS),strs(MNDIS),sumslip,sumvel,sumacc,maxvel,maxacc, ind, postmom, vini, postmom1, postmom2
	  integer i, j, jto, jfrom, ito, ifrom,k1,k2,kk, i2, j2, kout,iout, k
	  real tout,flv,fss,psiss
      COMMON /path/ kmax,kount,dxsav

      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
	  hdid=0.d0
	  sumslip=0.d0
      y=ystart
      if (kmax.gt.0) xsav=x-2.d0*dxsav
	  nstp=1
	  dtout=(x2-x1)/(NTS-1)
	  tout=0.
	  k1=0
	  k2=0
	  kout=0
	  iout=0
	  MSX=0.
	  strs=0.
	  
	  !open(111,FILE='mtildeS.txt')
16    continue
        !print*, x, tout
        call rates(x,y,scennum,dydx)
        yscal=abs(y)+abs(h*dydx)+TINY
	    slip(1:MNDIS)=slip(1:MNDIS)+y(1:MNDIS)*hdid
		postmom=sum(slip)*gtemp*dhh*dhh
		postmom1=0.
		postmom2=0.
		k=0
		 do j=1,NDIS
            do i=1,MDIS
		      k=k+1
		      if (j>27) then
				postmom1=postmom1+slip(k)
			  else
				postmom2=postmom2+slip(k)
			  endif
			    
	        enddo
	      enddo
		postmom1=postmom1*gtemp*dhh*dhh
		postmom2=postmom2*gtemp*dhh*dhh
		
        maxvel=maxval(y(1:MNDIS))
        maxacc=maxval(dydx(1:MNDIS))
	    sumvel=sum(y(1:MNDIS))*DX*DY
	    sumacc=sum(dydx(1:MNDIS))*DX*DY
        sumslip=sumslip+sumvel*hdid

	
		if (x>tout) then
		kout=kout+1
		TS(kout)=x
		  do j=1,NW
            jto=max(1,int(dW/dhh*j))+1
            jfrom=min(jto,int(dW/dhh*(j-1))+1)
            do i=1,NL
			  k2=k2+1
              ifrom=int(dL/dhh*(i-1))+1
              ito=int(dL/dhh*i)+1
              !kk=k2!((j-1)*NL+i-1)*NTS+k2
			  ind=0.
			  do j2=jfrom,jto
			  do i2=ifrom,ito
			    ind=ind+1.
                MSX(k2)=MSX(k2)+slip((j2-1)*MDIS+i2)
			  enddo
		      enddo
			  MSX(k2)=MSX(k2)/ind
			  !write(111,*)MSX(k2)
            enddo
          enddo
		  
		  k=0
          do j=1,NDIS
            do i=1,MDIS
		      k=k+1
		      strs(k)=SnPS(i,j)*aPS(i,j)*asinh(y(k)/v0*exp(y(k+MNDIS)/aPS(i,j)))	
	        enddo
	      enddo
	
		  tout=TS(kout)+dtout

		endif
	    
     !   write(110,'(10000E13.5)')x,y(1:MNDIS)
      !  write(115)x,hdid,sumacc,sumvel,sumslip,maxvel,maxacc

        if (ioutput.eq.1) then
          iout=iout+1
          !if (iout>50) then
            !print*,x, maxvel, hnext
			write(114,'(E18.11,1000000E13.5)') x, strs
			write(113,'(E18.11,1000000E13.5)') x, postmom, postmom1, postmom2
            write(112,'(E18.11,1000000E13.5)')x,slip(1:MNDIS)!, y(1+MNDIS:2*MNDIS)
            write(111,'(E18.11,1000000E13.5)')x,y(1+MNDIS:2*MNDIS)
            iout=0
		  !endif
          write(110,'(E18.11,10E13.5)') x,hdid,sumacc,sumvel,sumslip,maxvel,maxacc!,strsum,frsum
        endif


        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
#ifndef MPI
	     ! write(*,*)kount+1,':',x/31536000.d0
#endif
            if(kount.lt.kmax-1)then
              kount=kount+1
              xsav=x
            endif
          endif
        endif
		
        if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x

      call rkqs(y,dydx,MNDIS2,x,h,eps,yscal,hdid,hnext,scennum)
	  
          if(hdid.eq.h)then
            nok=nok+1
          else
            nbad=nbad+1
          endif
		  
          if((x-x2)*(x2-x1).ge.0.d0)then
            
	        ystart=y
            if(kmax.ne.0)then
              kount=kount+1
            endif
			if (kout<NTS) then
			print*, tout
			    do j=kout,NTS
			        TS(j)=tout
					MSX((j-1)*NL*NW+1:j*NL*NW)=MSX((kout-2)*NL*NW+1:(kout-1)*NL*NW)
					tout=tout+dtout
			    enddo
			endif
            return
          endif
		
       !   if(abs(hnext).lt.hmin) pause 'stepsize smaller than minimum in odeint'
            h=hnext
            if(MAXSTP>0.and.nstp>MAXSTP)then
              print*, 'too many steps in odeint'
			  if (kout<NTS) then
			    do j=kout,NTS
			        TS(j)=tout
					MSX((j-1)*NL*NW+1:j*NL*NW)=MSX((kout-2)*NL*NW+1:(kout-1)*NL*NW)
					tout=tout+dtout
			    enddo
			  endif
              return
	        endif
	      if (abs(maxvel)<1.e-12) then
		    print*,'no further movement'
			if (kout<NTS) then
			    do j=kout,NTS
			        TS(j)=tout
					MSX((j-1)*NL*NW+1:j*NL*NW)=MSX((kout-2)*NL*NW+1:(kout-1)*NL*NW)
					tout=tout+dtout
			    enddo
			endif
		    return
	      endif
	goto 16
	!close(111)
17  END

#endif

	SUBROUTINE readSGFs
	
	use SlipRates_com
	use PostSeismic_com
	use SlipRates_com
	use source_com, only:ioutput
	
	implicit none
	
	integer i,j,k,kk
	real dum
	
	open(10,file='inputgps.dat',action='read')
    read(10,*) T1S,T2S
	read(10,*) NTS, NTSrv
	read(10,*) NGPS
	read(10,*) SigmaGPS
	read(10,*) SigmaGPS2
	
	close (10)
    
	allocate(MSX(NL*NW*NTS),MSZ(NL*NW*NTS),TS(NTS))
	allocate(gpsgfN(NL*NW,NGPS),gpsgfE(NL*NW,NGPS),gpsgfZ(NL*NW,NGPS))
	allocate(gpssyntN(NTSrv,NGPS),gpssyntE(NTSrv,NGPS),gpssyntZ(NTSrv,NGPS))
	allocate(gpsrealN(NTSrv,NGPS), gpsrealE(NTSrv,NGPS), gpsrealZ(NTSrv,NGPS), gpsrealT(NTSrv,NGPS))
	allocate(gpsrealTN(NGPS))

	MSX=0.
	MSZ=0.
	TS=0.
	gpsgfN=0.
	gpsgfE=0.
	gpsgfZ=0.
	gpssyntN=0.
	gpssyntE=0.
	gpssyntZ=0.
	gpsrealN=0.
	gpsrealE=0.
	gpsrealZ=0.
	gpsrealT=0.
	
	!read Green function
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
    close(10)
	
	!read measured GPS movement
	open(10,file='rvgpsnez.dat',action='read')
	!open(11,file='rvgpse.dat',action='read')
	!open(12,file='rvgpsz.dat',action='read')
	

	!read(11,*)gpsrealTN(1:NGPS)
	!read(12,*)gpsrealTN(1:NGPS)
	
	do i=1,NGPS
	
	   read(10,*)gpsrealTN(i)
	   
	   do j=1,gpsrealTN(i)
		  read(10,*) gpsrealT(j,i),gpsrealN(j,i), gpsrealE(j,i), gpsrealZ(j,i)
	   enddo
	   
	   read(10,*)
	
	enddo
	close(10)
!	close(11)
!	close(12)
	
	!read sigma
	allocate(gpssigma(3,NGPS))
	open(10,file='gpsinfo.dat',action='read')
	do j=1,NGPS
    read(10,*) dum, dum, dum,gpssigma(1,j),gpssigma(2,j),gpssigma(3,j)
	
	enddo
	
	
	END 
	
	
	SUBROUTINE CalcSyntGPS()
	
	use SlipRates_com
	use PostSeismic_com
	use SlipRates_com
	use source_com, only:ioutput
	
	implicit none
	
	integer i,j,k,kk, jto2, jfrom2, jto, jfrom, i2, j2
	real, allocatable, dimension(:):: MSXtemp,MSZtemp 
	real stallocN, stallocE
	
	allocate(MSXtemp(NL*NW),MSZtemp(NL*NW))
	
	MSXtemp=0.
	MSZtemp=0.
	if (ioutput.eq.1) then
	   open(10,file='gpssyntnez.dat',action='write')
	endif
	gpssyntN=0.
	gpssyntE=0.
	gpssyntZ=0.	
		  
!	do j=1,NTS
!	    write (10,'(3E13.5)', advance="no") TS(j)
!		jfrom=(j-1)*NL*NW+1
!		jto=(j)*NL*NW
!		do i=1,NGPS
 !           gpssyntN(j,i)=dot_product(MSX(jfrom:jto),gpsgfN(1:NL*NW,i))
	!        gpssyntE(j,i)=dot_product(MSX(jfrom:jto),gpsgfE(1:NL*NW,i))
	 !       gpssyntZ(j,i)=dot_product(MSX(jfrom:jto),gpsgfZ(1:NL*NW,i))
		!  write (10,'(3E13.5)', advance="no") gpssyntN(j,i), gpssyntE(j,i), gpssyntZ(j,i)
!		enddo
!        write(10,*)
!	enddo
	
	do i=1,NGPS	
	  do j=1,gpsrealTN(i)
	    do j2=1,NTS
	    
		if (gpsrealT(j,i)<=TS(j2)) then
	
		  if (NTS.eq.1) then
		    MSXtemp(1:NL*NW)=MSX(1:NL*NW)
			MSZtemp(1:NL*NW)=MSZ(1:NL*NW)
		  else
		    jfrom=(j2-2)*NL*NW+1
		    jto=(j2-1)*NL*NW		
		    jfrom2=(j2-1)*NL*NW+1
		    jto2=(j2)*NL*NW	
		    if (j2==1) then
			  print*, 'error - gpsreaT(1) has to be > 0'
			endif
			MSXtemp(1:NL*NW)=(MSX(jfrom:jto)*(TS(j2)-gpsrealT(j,i))+MSX(jfrom2:jto2)*(gpsrealT(j,i)-TS(j2-1)))&
		    /(TS(j2)-TS(j2-1))	
			!no postseismic slip for dipslip earthquakes!
		  endif
		  if (ioutput.eq.1) then  
		     write (10,'(3E13.5)', advance="no") gpsrealT(j,i)
          endif
#if defined DIPSLIP
            gpssyntN(j,i)=dot_product(MSZtemp(1:NL*NW),gpsgfN(1:NL*NW,i))
	        gpssyntE(j,i)=dot_product(MSZtemp(1:NL*NW),gpsgfE(1:NL*NW,i))
	        gpssyntZ(j,i)=dot_product(MSZtemp(1:NL*NW),gpsgfZ(1:NL*NW,i))	
#else
            gpssyntN(j,i)=dot_product(MSXtemp(1:NL*NW),gpsgfN(1:NL*NW,i))
	        gpssyntE(j,i)=dot_product(MSXtemp(1:NL*NW),gpsgfE(1:NL*NW,i))
	        gpssyntZ(j,i)=dot_product(MSXtemp(1:NL*NW),gpsgfZ(1:NL*NW,i))	
#endif
            if (ioutput.eq.1) then
		       write (10,'(3E13.5)') gpssyntN(j,i), gpssyntE(j,i), gpssyntZ(j,i)			
			endif
		  exit
		endif
	    enddo
      enddo
	    if (ioutput.eq.1) then
	       write(10,*)
	    endif
	enddo
	if (ioutput.eq.1) then
      close(10)	
	
   	  open(10,file='mtildeslip2D-sGPS.out',action='write')
	  open(11,file='stations-GPS.dat')
	  open(12,file='mtildeslip2D-sGPSend.out',action='write')
	  do i=1,NGPS
		read (11,*) stallocN, stallocE 
	    write(10,'(5E13.5)') stallocN, stallocE, gpssyntN(1,i), gpssyntE(1,i), gpssyntZ(1,i)	
	    write(12,'(5E13.5)') stallocN, stallocE, gpssyntN(NTSrv,i)-gpssyntN(1,i), gpssyntE(NTSrv,i)-gpssyntE(1,i), gpssyntZ(NTSrv,i)-gpssyntZ(1,i)	
	  enddo
	  close(10)
	  close(11)
	  close(12)
	 call plotgps_old()
	endif
	deallocate(MSXtemp)


	END 
	
	SUBROUTINE plotgps_old()
    USE waveforms_com
    USE SlipRates_com
	USE PostSeismic_com
	
    IMPLICIT NONE
    REAL, PARAMETER:: margin=0.05
    CHARACTER*5,ALLOCATABLE,DIMENSION(:):: staname
    REAL startx,starty,stept
    REAL,ALLOCATABLE,DIMENSION(:):: stepa,maxampl
    real dum
    INTEGER i,j,k,jj

    allocate(stepa(NGPS),maxampl(NGPS),staname(NGPS))
    open(233,FILE='gpsinfo.dat')
    do k=1,NGPS
      read(233,*)(dum,i=1,7),staname(k)
    enddo
    close(233)

    !stept=1./dble(3)*(1.-margin)/(NTSrv)
    jj=0
    do j=1,NGPS
      if (gpsrealTN(j)>1) then
      maxampl(j)=0.
      maxampl(j)=max(maxampl(j),maxval(abs(gpssyntN(:,j)-gpssyntN(1,j))))
	  maxampl(j)=max(maxampl(j),maxval(abs(gpsrealN(:,j)-gpsrealN(1,j))))
      maxampl(j)=max(maxampl(j),maxval(abs(gpssyntE(:,j)-gpssyntE(1,j))))
	  maxampl(j)=max(maxampl(j),maxval(abs(gpsrealE(:,j)-gpsrealE(1,j))))
      maxampl(j)=max(maxampl(j),maxval(abs(gpssyntZ(:,j)-gpssyntZ(1,j))))
	  maxampl(j)=max(maxampl(j),maxval(abs(gpsrealZ(:,j)-gpsrealZ(1,j))))
      stepa(j)=.5/dble(NGPS)*(1.-margin)/maxampl(j)
      endif
    enddo

    open(205,FILE='gpsplot.real.dat')
    jj=0
    do j=1,NGPS
	 if (gpsrealTN(j)>1) then
      starty=((NGPS-j+1)+margin/2.)/dble(NGPS+1)
      startx=((1-1)+margin/2.)/2.
	  write(205,'(3E13.5)')startx,starty
	  stept=1./dble(3)*(1.-margin)/(gpsrealT(gpsrealTN(j),j))
      do k=1,gpsrealTN(j)
        write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(gpsrealN(k,j)-gpsrealN(1,j))
      enddo
      write(205,*)
      startx=((2-1)+margin/2.)/2.
	  write(205,'(3E13.5)')startx,starty
      do k=1,gpsrealTN(j)
        write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(gpsrealE(k,j)-gpsrealE(1,j))
      enddo
      write(205,*)
      !startx=((3-1)+margin/2.)/3.
	  !write(205,'(3E13.5)')startx,starty
      !do k=1,NTSrv
      !  write(205,'(3E13.5)')startx+stept*(k),starty+stepa(j)*(gpsrealZ(k,j)-gpsrealZ(1,j))
      !enddo
      !write(205,*)
	 endif
    enddo
	
    close(205)
    open(205,FILE='gpsplot.synt.dat')
    do j=1,NGPS
	 if (gpsrealTN(j)>1) then
      starty=((NGPS-j+1)+margin/2.)/dble(NGPS+1)
      startx=((1-1)+margin/2.)/2.
	  stept=1./dble(3)*(1.-margin)/(gpsrealT(gpsrealTN(j),j))
	  !write(205,'(3E13.5)')startx,starty
      do k=1,gpsrealTN(j)
        write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(gpssyntN(k,j)-gpssyntN(1,j))
      enddo
      write(205,*)
      startx=((2-1)+margin/2.)/2.
	  !write(205,'(3E13.5)')startx,starty
      do k=1,gpsrealTN(j)
        write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(gpssyntE(k,j)-gpssyntE(1,j))
      enddo
      write(205,*)
      !startx=((3-1)+margin/2.)/3.
	  !!write(205,'(3E13.5)')startx,starty
      !do k=1,NTSrv
      !  write(205,'(3E13.5)')startx+stept*k,starty+stepa(j)*(gpssyntZ(k,j)-gpssyntZ(1,j))
      !enddo
      !write(205,*)
	 endif
    enddo
    close(205)

    open(201,FILE='gpsplot.gp')
    write(201,*)'set term postscript portrait color solid enh'
    write(201,*)'set output "gpsplot.ps"'
    write(201,*)'set xrange [0:1]'
    write(201,*)'set yrange [0:1]'
    write(201,*)'unset xtics'
    write(201,*)'unset ytics'
    write(201,*)'set border 0'
    do j=1,NGPS
	 if (gpsrealTN(j)>1) then
      write(201,'(A11,F5.2,A23,G,A6)')'set label "',real(maxampl(j))*100.,'" at graph 1.09, graph ',((NGPS-j+1)+margin/2.)/dble(NGPS+1),' right'
      write(201,*)'set label "'//staname(j)//'" at graph -0.01, graph ',((NGPS-j+1)+margin/2.)/dble(NGPS+1),' right'
     endif
	enddo
    write(201,*)'set label "N-S" at .15,1'
    write(201,*)'set label "E-W" at .48,1'
    write(201,*)'set label "Z" at .86,1'
    write(201,*)'plot "gpsplot.real.dat" u 1:2 notitle w l lt -1 lw 1,\'
    write(201,*)'"gpsplot.synt.dat" u 1:2 notitle w l lt 1 lw 2'
    close(201)
    
    END

	
	SUBROUTINE plotgps()
    USE waveforms_com
    USE SlipRates_com
	USE PostSeismic_com
	
    IMPLICIT NONE
    REAL, PARAMETER:: margin=0.02
    CHARACTER*5,ALLOCATABLE,DIMENSION(:):: staname
    REAL startx,starty,stept
    REAL,ALLOCATABLE,DIMENSION(:):: stepa,maxampl
    real dum
    INTEGER i,j,k,jj, n_align, n_cos

	n_align=4
	n_cos=48
	
    allocate(stepa(NGPS),maxampl(NGPS),staname(NGPS))
    open(233,FILE='gpsinfo.dat')
    do k=1,NGPS-n_cos
      read(233,*)(dum,i=1,7),staname(k)
    enddo
    close(233)
	


    !stept=1./dble(3)*(1.-margin)/(NTSrv)
    jj=0
    do j=1,NGPS-n_align-n_cos
      if (gpsrealTN(j)>1) then
      maxampl(j)=0.
      maxampl(j)=max(maxampl(j),maxval(abs(gpssyntN(:,j)-gpssyntN(1,j))))
	  maxampl(j)=max(maxampl(j),maxval(abs(gpsrealN(:,j)-gpsrealN(1,j))))
      maxampl(j)=max(maxampl(j),maxval(abs(gpssyntE(:,j)-gpssyntE(1,j))))
	  maxampl(j)=max(maxampl(j),maxval(abs(gpsrealE(:,j)-gpsrealE(1,j))))
      maxampl(j)=max(maxampl(j),maxval(abs(gpssyntZ(:,j)-gpssyntZ(1,j))))
	  maxampl(j)=max(maxampl(j),maxval(abs(gpsrealZ(:,j)-gpsrealZ(1,j))))
	  maxampl(j)=maxampl(j)/2
      stepa(j)=.5/dble(NGPS-n_cos)*(1.-margin)/maxampl(j)
      endif
    enddo
	
	do j=NGPS-n_align+1-n_cos,NGPS-n_cos
      if (gpsrealTN(j)>1) then
      maxampl(j)=0.
      maxampl(j)=max(maxampl(j),maxval(abs(sqrt(gpssyntN(:,j)**2+gpssyntE(:,j)**2)-sqrt(gpssyntN(1,j)**2+gpssyntE(1,j)**2))))
	  maxampl(j)=max(maxampl(j),maxval(abs(sqrt(gpsrealN(:,j)**2+gpsrealE(:,j)**2)-sqrt(gpsrealN(1,j)**2+gpsrealE(1,j)**2))))
	  maxampl(j)=maxampl(j)/2
      stepa(j)=.5/dble(NGPS-n_cos)*(1.-margin)/maxampl(j)
      endif
    enddo

    open(205,FILE='gpsplot.real.dat')
    jj=0
    do j=1,NGPS-n_align-n_cos
	 if (gpsrealTN(j)>1) then
      starty=((NGPS-n_cos-j+1)+margin/2.)/dble(NGPS-n_cos+1)
      startx=((1-1)+margin/2.)/2.
	  !write(205,'(3E13.5)')startx,starty,0.01*gpssigma(1,j)/maxampl(j)
	  stept=1./dble(3)*(1.-margin)/(gpsrealT(gpsrealTN(j),j))
      do k=1,gpsrealTN(j)
        write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(gpsrealN(k,j)-gpsrealN(1,j)),0.01*gpssigma(1,j)/maxampl(j)
      enddo
      write(205,*)
      startx=((2-1)+margin/2.)/2.
	 ! write(205,'(3E13.5)')startx,starty,0.01*gpssigma(2,j)/maxampl(j)
      do k=1,gpsrealTN(j)
        write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(gpsrealE(k,j)-gpsrealE(1,j)),0.01*gpssigma(2,j)/maxampl(j)
      enddo
      write(205,*)
      !startx=((3-1)+margin/2.)/3.
	  !write(205,'(3E13.5)')startx,starty
      !do k=1,NTSrv
      !  write(205,'(3E13.5)')startx+stept*(k),starty+stepa(j)*(gpsrealZ(k,j)-gpsrealZ(1,j))
      !enddo
      !write(205,*)
	 endif
    enddo
	
	do j=NGPS-n_align+1-n_cos,NGPS-n_cos
	 if (gpsrealTN(j)>1) then
      starty=((NGPS-n_cos-j+1)+margin/2.)/dble(NGPS-n_cos+1)
      startx=((1.5-1)+margin/2.)/2.
	  !write(205,'(3E13.5)')startx,starty,0.01*gpssigma(1,j)/maxampl(j)
	  stept=1./dble(3)*(1.-margin)/(gpsrealT(gpsrealTN(j),j))
      do k=1,gpsrealTN(j)
        write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(sqrt(gpsrealN(k,j)**2+gpsrealE(k,j)**2)-sqrt(gpsrealN(1,j)**2+gpsrealE(1,j)**2)),0.01*gpssigma(1,j)/maxampl(j)
      enddo
      write(205,*)
    !  startx=((2-1)+margin/2.)/2.
	 ! write(205,'(3E13.5)')startx,starty,0.01*gpssigma(2,j)/maxampl(j)
   !   do k=1,gpsrealTN(j)
   !     write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(gpsrealE(k,j)-gpsrealE(1,j)),0.01*gpssigma(2,j)/maxampl(j)
   !   enddo
   !   write(205,*)
      !startx=((3-1)+margin/2.)/3.
	  !write(205,'(3E13.5)')startx,starty
      !do k=1,NTSrv
      !  write(205,'(3E13.5)')startx+stept*(k),starty+stepa(j)*(gpsrealZ(k,j)-gpsrealZ(1,j))
      !enddo
      !write(205,*)
	 endif
    enddo
	
    close(205)
	
	
    open(205,FILE='gpsplot.synt.dat')
    do j=1,NGPS-n_align-n_cos
	 if (gpsrealTN(j)>1) then
      starty=((NGPS-n_cos-j+1)+margin/2.)/dble(NGPS-n_cos+1)
      startx=((1-1)+margin/2.)/2.
	  stept=1./dble(3)*(1.-margin)/(gpsrealT(gpsrealTN(j),j))
	  !write(205,'(3E13.5)')startx,starty
      do k=1,gpsrealTN(j)
        write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(gpssyntN(k,j)-gpssyntN(1,j))
      enddo
      write(205,*)
      startx=((2-1)+margin/2.)/2.
	  !write(205,'(3E13.5)')startx,starty
      do k=1,gpsrealTN(j)
        write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(gpssyntE(k,j)-gpssyntE(1,j))
      enddo
      write(205,*)
      !startx=((3-1)+margin/2.)/3.
	  !!write(205,'(3E13.5)')startx,starty
      !do k=1,NTSrv
      !  write(205,'(3E13.5)')startx+stept*k,starty+stepa(j)*(gpssyntZ(k,j)-gpssyntZ(1,j))
      !enddo
      !write(205,*)
	 endif
    enddo
	
	do j=NGPS-n_align+1-n_cos,NGPS-n_cos
	 if (gpsrealTN(j)>1) then
      starty=((NGPS-n_cos-j+1)+margin/2.)/dble(NGPS-n_cos+1)
      startx=((1.5-1)+margin/2.)/2.
	  stept=1./dble(3)*(1.-margin)/(gpsrealT(gpsrealTN(j),j))
	  !write(205,'(3E13.5)')startx,starty
      do k=1,gpsrealTN(j)
        write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(sqrt(gpssyntN(k,j)**2+gpssyntE(k,j)**2)-sqrt(gpssyntN(1,j)**2+gpssyntE(1,j)**2))
      enddo
      write(205,*)
    !  startx=((2-1)+margin/2.)/2.
	  !write(205,'(3E13.5)')startx,starty
     ! do k=1,gpsrealTN(j)
     !   write(205,'(3E13.5)')startx+stept*gpsrealT(k,j),starty+stepa(j)*(gpssyntE(k,j)-gpssyntE(1,j))
    !  enddo
     ! write(205,*)
      !startx=((3-1)+margin/2.)/3.
	  !!write(205,'(3E13.5)')startx,starty
      !do k=1,NTSrv
      !  write(205,'(3E13.5)')startx+stept*k,starty+stepa(j)*(gpssyntZ(k,j)-gpssyntZ(1,j))
      !enddo
      !write(205,*)
	 endif
    enddo
	
    close(205)

    open(201,FILE='gpsplot.gp')
    write(201,*)'set term postscript portrait color solid enh'
    write(201,*)'set output "gpsplot.ps"'
    write(201,*)'set xrange [0:1]'
    write(201,*)'set yrange [0:1]'
    write(201,*)'unset xtics'
    write(201,*)'unset ytics'
    write(201,*)'set border 0'
    do j=1,NGPS-n_cos
	 if (gpsrealTN(j)>1) then
      write(201,'(A11,F5.2,A23,G,A6)')'set label "',real(maxampl(j))*100.,'" at graph 1.09, graph ',((NGPS-n_cos-j+1)+margin/2.)/dble(NGPS-n_cos+1),' right'
      write(201,*)'set label "'//staname(j)//'" at graph -0.01, graph ',((NGPS-n_cos-j+1)+margin/2.)/dble(NGPS-n_cos+1),' right'
     endif
	enddo
    write(201,*)'set label "N-S" at .15,1'
    write(201,*)'set label "E-W" at .48,1'
    write(201,*)'set label "Z" at .86,1'
    write(201,*)'plot "gpsplot.real.dat" notitle with errorbars lt rgb "black",\'
    write(201,*)'"gpsplot.synt.dat" u 1:2 notitle w l lt 1 lw 2'
    close(201)
    
    END

	
