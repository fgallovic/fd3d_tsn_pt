#if defined FVW

#else
    MODULE randomdynmod_com
    REAL*8 dkx,dky
    INTEGER idum
    INTEGER NLout,NWout
    END MODULE
    
    SUBROUTINE randomdynmod(NL,NW,dh,T0I,TsI,DcI)
    USE randomdynmod_com
    USE frictionconstraints_com
    IMPLICIT NONE
    INTEGER NL,NW
    REAL dh
    REAL,DIMENSION(NL,NW):: DcI,TsI,T0I                        !original parameters (input)
    REAL*8,ALLOCATABLE,DIMENSION(:,:):: t0,t1,Dc,t0R,t1R,DcR   !randomized parameters (output)
    REAL*8 EXPt0,KCt0,Strengtht0,Multiplyt0
    REAL*8 EXPt1,KCt1,Strengtht1,Multiplyt1
    REAL*8 EXPDc,KCDc,StrengthDc,MultiplyDc
    INTEGER APPLYt0,APPLYt1,APPLYDc
    REAL*8 dum
    INTEGER i,j

    write(*,*)'Randomizing dynamic parameters:'
!Read the input parameters
    open(100,FILE='randomdynmod.in')
    NLout=2**int(log(dble(NL))/log(2.d0)+1)
    NWout=2**int(log(dble(NW))/log(2.d0)+1)
    read(100,*)
    read(100,*)EXPt0,Strengtht0,KCt0,Multiplyt0,APPLYt0
    read(100,*)EXPt1,Strengtht1,KCt1,Multiplyt1,APPLYt1
    read(100,*)EXPDc,StrengthDc,KCDc,MultiplyDc,APPLYDc
    read(100,*)
    read(100,*)idum
    close(100)
    allocate(t0(NLout,NWout),t1(NLout,NWout),Dc(NLout,NWout))
    allocate(t0R(NLout,NWout),t1R(NLout,NWout),DcR(NLout,NWout))
    dkx=1.d0/dble(NL);dky=1.d0/dble(NW)

!Save the original model in Gnuplot format
    open(201,FILE='randomdynmod.gnuplot')
    do j=1,NW
      write(201,'(10000E13.5)')(T0I(i,j)/1.d6,i=1,NL)
    enddo
    write(201,*)
    write(201,*)
    do j=1,NW
      write(201,'(10000E13.5)')(TsI(i,j)/1.d6,i=1,NL)
    enddo
    write(201,*)
    write(201,*)
    do j=1,NW
      write(201,'(10000E13.5)')((TsI(i,j)-T0I(i,j))/1.d6,i=1,NL)
    enddo
    write(201,*)
    write(201,*)
    do j=1,NW
      write(201,'(10000E13.5)')(DcI(i,j),i=1,NL)
    enddo
    write(201,*)
    write(201,*)
    
    t0=0.d0;t1=0.d0;Dc=0.d0
    t0(1:NL,1:NW)=T0I(1:NL,1:NW)
    t1(1:NL,1:NW)=TsI(1:NL,1:NW)
    Dc(1:NL,1:NW)=DcI(1:NL,1:NW)
    
! Apply the randomization
    CALL randomize(t0,t0R,NL,NW,EXPt0,KCt0,Strengtht0,Multiplyt0,APPLYt0)
    CALL randomize(t1,t1R,NL,NW,EXPt1,KCt1,Strengtht1,Multiplyt1,APPLYt1)
    CALL randomize(Dc,DcR,NL,NW,EXPDc,KCDc,StrengthDc,MultiplyDc,APPLYDc)
    
    do j=1,NW
      do i=1,NL
        if(t1R(i,j)<peak_xzMin)then
          t1R(i,j)=0.d0
        endif
        if(DcR(i,j)>DcMin)then
          DcI(i,j)=DcR(i,j)
        endif
        if(T0I(i,j)<TSI(i,j).and.t1R(i,j)>t0R(i,j))then    !Overwrite outside the actual nucl. zone without creating new nucl.zone
          T0I(i,j)=t0R(i,j)
          TsI(i,j)=t1R(i,j)    
        endif
      enddo
    enddo
    
    write(*,*)'T0 min,max: ',minval(T0I),maxval(T0I)
    write(*,*)'Ts min,max: ',minval(TsI),maxval(TsI)
    write(*,*)'Dc min,max: ',minval(DcI),maxval(DcI)
    
!Save the randomized model in Gnuplot format
    do j=1,NW
      write(201,'(10000E13.5)')(T0I(i,j)/1.d6,i=1,NL)
    enddo
    write(201,*)
    write(201,*)
    do j=1,NW
      write(201,'(10000E13.5)')(TsI(i,j)/1.d6,i=1,NL)
    enddo
    write(201,*)
    write(201,*)
    do j=1,NW
      write(201,'(10000E13.5)')((TsI(i,j)-T0I(i,j))/1.d6,i=1,NL)
    enddo
    write(201,*)
    write(201,*)
    do j=1,NW
      write(201,'(10000E13.5)')(DcI(i,j),i=1,NL)
    enddo
    close(201)
    
    END
    
    SUBROUTINE randomize(A,B,NL,NW,EXPONENT,KC,Strength,Multiply,Apply)
    USE randomdynmod_com
    IMPLICIT NONE
    REAL*8 A(NLout,NWout),B(NLout,NWout)
    REAL*8 EXPONENT,KC,Strength,Multiply
    INTEGER Apply
    INTEGER NL,NW
    REAL*8 Bdum(NLout,NWout),gasdev2,kx,ky,krad,kcrad,KCx,KCy,rms
    COMPLEX*16 BC(NLout/2,NWout),speqBC(NWout)
    COMPLEX*16 BdumC(NLout/2,NWout),speqBdumC(NWout)
    INTEGER i,j

    if(Apply==0)return
!FT of the deterministic part
    B=A
    rms=sqrt(max(0.d0,sum(B(1:NL,1:NW)**2)/dble(NL*NW)-(sum(B(1:NL,1:NW))/dble(NL*NW))**2))

    CALL rlft3(B,speqBC,NLout,NWout,1,1)
    do i=1,NLout/2
      BC(i,:)=cmplx(B(2*i-1,:),B(2*i,:))/dble(NLout/2*NWout)
    enddo
    speqBC=speqBC/dble(NLout/2*NWout)

!FT of white noise
    do j=1,NWout
      do i=1,NLout
        Bdum(i,j)=gasdev2(idum)
      enddo
    enddo
    Bdum=Bdum*rms*Strength
    CALL rlft3(Bdum,speqBdumC,NLout,NWout,1,1)
    do i=1,NLout/2
      BdumC(i,:)=cmplx(Bdum(2*i-1,:),Bdum(2*i,:))/dble(NLout/2*NWout)
    enddo
    speqBdumC=speqBdumC/dble(NLout/2*NWout)
    
!Adding k^-2 by using the white noise spectrum
    kcrad=KC*sqrt(dkx**2+dky**2)                !Corner wave-number after which the random spectrum is applied
    KCx=sqrt(dkx**2+dky**2)              !Corner wave-number for along-strike direction
    KCy=sqrt(dkx**2+dky**2)              !Corner wave-number for down-dip direction

    do j=1,NWout
      if(j<=NWout/2+1)then
        ky=1./dble(NWout)*real(j-1)
      else
        ky=-1./dble(NWout)*real(NWout-j+1)
      endif
      do i=1,NLout/2+1
        kx=1./dble(NLout)*real(i-1)
        krad=sqrt(kx**2+ky**2)
        if(i<NLout/2+1)then
          if(krad>=kcrad)then
            BC(i,j)=BdumC(i,j)/sqrt(1.+((kx/KCx)**2+(ky/KCy)**2)**EXPONENT)
          endif
        elseif(krad>=kcrad)then
          speqBC(j)=speqBdumC(j)/sqrt(1.+((kx/KCx)**2+(ky/KCy)**2)**EXPONENT)
        endif
      enddo
    enddo

!Back FT
    CALL rlft3(BC,speqBC,NLout,NWout,1,-1)
    do i=1,NLout/2
      B(2*i-1,:)=real(BC(i,:))*Multiply
      B(2*i,:)=imag(BC(i,:))*Multiply
    enddo

    END
    
!NR 
    
      FUNCTION gasdev2(idum)
      INTEGER idum
      DOUBLE PRECISION gasdev2
      INTEGER iset
      DOUBLE PRECISION fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
1       v1=2.d0*ran1(idum)-1.d0
        v2=2.d0*ran1(idum)-1.d0
        rsq=v1**2+v2**2
        if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1
        fac=sqrt(-2.d0*log(rsq)/rsq)
        gset=v1*fac
        gasdev2=v2*fac
        iset=1
      else
        gasdev2=gset
        iset=0
      endif
      return
    END

     FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      DOUBLE PRECISION ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1.d0/IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=3.d-16,RNMX=1.d0-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END


      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
      INTEGER isign,nn1,nn2,nn3
      COMPLEX*16 data(nn1/2,nn2,nn3),speq(nn2,nn3)
      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      COMPLEX*16 c1,c2,h1,h2,w
      c1=dcmplx(0.5d0,0.0d0)
      c2=dcmplx(0.0d0,-0.5d0*isign)
      theta=6.28318530717959d0/dble(isign*nn1)
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      nn(1)=nn1/2
      nn(2)=nn2
      nn(3)=nn3
      if(isign.eq.1)then
        call fourn(data,nn,3,isign)
        do 12 i3=1,nn3
          do 11 i2=1,nn2
            speq(i2,i3)=data(1,i2,i3)
11        continue
12      continue
      endif
      do 15 i3=1,nn3
        j3=1
        if (i3.ne.1) j3=nn3-i3+2
        wr=1.0d0
        wi=0.0d0
        do 14 i1=1,nn1/4+1
          j1=nn1/2-i1+2
          do 13 i2=1,nn2
            j2=1
            if (i2.ne.1) j2=nn2-i2+2
            if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
            else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
            endif
13        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          w=dcmplx(dble(wr),dble(wi))
14      continue
15    continue
      if(isign.eq.-1)then
        call fourn(data,nn,3,isign)
      endif
      return
      END

    
      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      DOUBLE PRECISION data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=dble(wr)*data(k2)-dble(wi)*data(k2+1)
                tempi=dble(wr)*data(k2+1)+dble(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END
#endif
    