    IMPLICIT NONE
    REAL*8,PARAMETER:: PI=3.1415926535d0
    REAL*8 L,W,hypoW,hdepth
    REAL*4 dip
    INTEGER nxt,nyt,nzt
    REAL*8 dhx,dhz
    REAL*8 dum
    REAL*8,ALLOCATABLE,DIMENSION(:,:):: T0I,TSI,DcI
! Ellipse parameters:
    REAL*8 ElL,ElW,ElNucL,ElNucW
    REAL*8 DcMin,DcRate
    REAL*8 InitStress,FricDrop

    INTEGER i,j,k,n,inuc,jnuc

!--------------------      
! Read the fault input file
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
    read(10,*)dum,dip
    read(10,*)
    read(10,*)hdepth
    read(10,*)
    read(10,*)L,W
    read(10,*)
    read(10,*)dum,hypoW
    close(10)

!--------------------      
! Read the input file for the dynamic simulations
!--------------------
    open(11, file='inputinv.dat', status='old')
    read(11,*)
    read(11,*) nxt,nzt
    close(11)
    dhx=L/real(nxt-1)
    dhz=W/real(nzt-1)
    allocate(T0I(nxt,nzt),TSI(nxt,nzt),DcI(nxt,nzt))

!----------------------------
    open(10,FILE='ellipticalmodelgener.inp')
    read(10,*) ElL,ElW !Dimensions of the Ellipse
    read(10,*) ElNucL,ElNucW !Position of the nucleation
    read(10,*) InitStress !Initial stress
    read(10,*) FricDrop !Friction drop - InitStress
    read(10,*) DcMin,DcRate ! Minimum Dc and its increase rate

    TSI(:,:)=InitStress+FricDrop
    do j=1,nzt
      T0I(:,j)=InitStress*16200.*dhz*(nzt-j) !*100.e6
    enddo
    inuc=int(ElNucL/L*(nxt-1))+1
    jnuc=int(ElNucW/W*(nzt-1))+1

    do j=1,nzt
      do i=1,nxt
        if(real(i-nxt/2)**2/(ElL/2./L*nxt)**2+real(j-nzt/2)**2/(ElW/2./W*nzt)**2>1)then
          TSI(i,j)=-InitStress
          T0I(i,j)=-(InitStress+FricDrop)*16200.*dhz*(nzt-j) !*100.e6
        endif
        DcI(i,j)=Dcmin+Dcrate*(sqrt((inuc-i)**2*dhx**2+(jnuc-j)**2*dhz**2)/1.e3)
      enddo
    enddo
    write(*,*)minval(DcI),maxval(DcI)

    do j=jnuc-1,jnuc+1
      T0I(inuc-1:inuc+1,j)=TSI(inuc-1:inuc+1,j)*1.1*16200.*dhz*(nzt-j) !*100.e6
    enddo

    open(201,FILE='forwardmodel.dat')   !model output
    write(201,'(1000000E13.5)')-1.,-1.,T0I(:,:),TsI(:,:),DcI(:,:)
    close(201)

    END


