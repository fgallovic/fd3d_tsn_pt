    IMPLICIT NONE

    INTEGER,PARAMETER:: densify=1
    REAL*8,PARAMETER::  PI=3.1415926535d0, perturb=0.1
    REAL*8 L,W,rake,hypoW,hdepth
    REAL*4 dip
    INTEGER NLI,NWI,NLIs,NWIs
    REAL*8 DL,DW,DLs,DWs
    REAL*8 dum,u,t
    REAL*8 dhx,dhz
    REAL*8,ALLOCATABLE,DIMENSION(:,:):: T0I,TSI,DcI,T0Is,TSIs,DcIs
    REAL*8 XS,ZS,NucL,NucW,NucR
    REAL*8,ALLOCATABLE,DIMENSION(:,:):: fractals
    INTEGER Nlevels,D,N0,Nfractals,idum
    REAL*8 diam,dc,Dcmin,Dcmax,Rmax,ran2
    INTEGER n,inuc,jnuc
    INTEGER i,j,k,ii,kk

!--------------------      
! Read the input file from the inversions
!--------------------
    open(10,FILE='input.dat',status='old')
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)L,W
    close(10)

!--------------------      
! Read the input file for the dynamic simulations
!--------------------
    open(11, file='inputinv.dat', status='old')
    read(11,*)
    read(11,*) NLIs,NWIs
    read(11,*)
    read(11,*)
    read(11,*)
    read(11,*)
    read(11,*)
    read(11,*) dum,NucL,NucW,NucR
    read(11,*)
    read(11,*)
    read(11,*) idum
    close(11)
    allocate(T0Is(NLIs,NWIs),TSIs(NLIs,NWIs),DcIs(NLIs,NWIs))
    write(*,*)'Densified from ',NLIs,NWIs

    NLI=(NLIs-1)*densify+1
    NWI=(NWIs-1)*densify+1
    allocate(T0I(NLI,NWI),TSI(NLI,NWI),DcI(NLI,NWI))
    write(*,*)'Densified to   ',NLI,NWI
    dhx=L/real(NLI-1)
    dhz=W/real(NWI-1)

    open(102,FILE='originalmodel.dat')
10  read(102,*)dum,dum,T0Is(:,:),TsIs(:,:),DcIs(:,:)
    close(102)

    DL=1./real(NLI-1)
    DW=1./real(NWI-1)
    DLs=1./real(NLIs-1)
    DWs=1./real(NWIs-1)
    do k=1,NWI
      ZS=real(k-1)*DW
      kk=min(NWIs-1,int(ZS/DWs)+1)
      u=(ZS-DWs*(kk-1))/DWs
      do i=1,NLI
        XS=real(i-1)*DL
        ii=min(NLIs-1,int(XS/DLs)+1)
        t=(XS-DLs*(ii-1))/DLs
        T0I(i,k)     =(1.-t)*(1.-u)*T0Is(ii,kk)+t*(1.-u)*T0Is(ii+1,kk)+t*u*T0Is(ii+1,kk+1)+(1.-t)*u*T0Is(ii,kk+1)
        TSI(i,k)     =(1.-t)*(1.-u)*TSIs(ii,kk)+t*(1.-u)*TSIs(ii+1,kk)+t*u*TSIs(ii+1,kk+1)+(1.-t)*u*TSIs(ii,kk+1)
        DcI(i,k)     =(1.-t)*(1.-u)*DcIs(ii,kk)+t*(1.-u)*DcIs(ii+1,kk)+t*u*DcIs(ii+1,kk+1)+(1.-t)*u*DcIs(ii,kk+1)
      enddo
    enddo
    deallocate(T0Is,TSIs,DcIs)

! Prepare fractal distributions (all normalized to 1)
    D=2
    Nlevels=5
    Dcmin=1./2**(Nlevels)
    Rmax=min(L,W)/8.

    allocate(T0Is(NLI,NWI),TSIs(NLI,NWI),DcIs(NLI,NWI))
    N0=floor(max(L/W,W/L)*2**(D*(Nlevels))) ! Number of the smallest subsources (n=0)
    Dcmax=Dcmin*2**(Nlevels)
    TSIs(:,:)=1.
    T0Is(:,:)=1.
    DcIs(:,:)=1.

    i=0
    do n=0,Nlevels-1
      do k=1,N0/2**(D*n)
        i=i+1
      enddo
    enddo
    Nfractals=i
    write(*,*)Nfractals
    allocate(fractals(4,Nfractals)) !x-coord, y-coord, diameter, Dc
    i=0
    do n=0,Nlevels-1
      diam=Rmax/2**(Nlevels-n-1)
      dc=Dcmin*2**(n)
write(*,*)N0/2**(D*n),diam,dc
      do k=1,N0/2**(D*n)
        i=i+1
        fractals(3,i)=diam
        fractals(4,i)=Dc
        XS=ran2(idum)*(L-2*diam)+diam
        ZS=ran2(idum)*(W-2*diam)+diam
        fractals(1,i)=XS
        fractals(2,i)=ZS
      enddo
    enddo
    do j=1,NWI
      ZS=dhz*(j-1)
      do i=1,NLI
        XS=dhx*(i-1)
        if((NucL-XS)**2+(NucW-ZS)**2>NucR**2)then
          do n=Nfractals,1,-1
            if((XS-fractals(1,n))**2+(ZS-fractals(2,n))**2<fractals(3,n)**2)then
              if(TSI(i,j)>0.)DcIs(i,j)=fractals(4,n)
            endif
          enddo
          if(TSI(i,j)>0.)then
            T0Is(i,j)=T0Is(i,j)*(1.-.2*log(DcIs(i,j)/Dcmin)/log(Dcmax/Dcmin)) !Jen aby tam byla nejaka nahodnost, navic korelovana s Dc
            TSIs(i,j)=TSIs(i,j)*(1.+.4*log(DcIs(i,j)/Dcmin)/log(Dcmax/Dcmin)) ! To by melo efektivne jeste zvetsovat Gc
          else
            T0Is(i,j)=T0Is(i,j)!*(1.+.4*log(DcIs(i,j)/Dcmin)/log(Dcmax/Dcmin))
            TSIs(i,j)=TSIs(i,j)!*(1.-.2*log(DcIs(i,j)/Dcmin)/log(Dcmax/Dcmin))
          endif
        endif
      enddo
    enddo

! Apply fractals
    T0I(:,:)=T0Is(:,:)*T0I(:,:)
    TSI(:,:)=TSIs(:,:)*TSI(:,:)
    DcI(:,:)=DcIs(:,:)*DcI(:,:)
    deallocate(T0Is,TSIs,DcIs)

    open(201,FILE='forwardmodel.dat')   !perturbed model
    write(201,'(1000000E13.5)')-1.,-1.,T0I(:,:),TsI(:,:),DcI(:,:)
    close(201)

    END

    
    !NR 

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      DOUBLE PRECISION ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,&
IQ2=52774,IR1=12211,IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=3.d-16,RNMX=1.d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END

    
