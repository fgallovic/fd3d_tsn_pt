module mod_pgamisf
integer :: nper
real, allocatable :: per(:),a(:),b(:),c(:),d(:),c1(:),tau(:)
real, allocatable :: sd(:),sv(:),sa1(:),sa2(:),psv(:),psa(:)
!real :: pga_theor
real :: damp=0.05
real :: eps=1.e-6

!contains
end module

subroutine pga_theor(x,mw,per1,y) !result(y)
use mod_pgamisf
implicit none
integer :: i
real :: y,x,mw,per1
i=1
!print *,i,per1
do while (abs(per1-per(i))>eps .and. i<=nper)
 i=i+1
enddo
!print *,i
if (i>nper) then
 print *,'period inconsistent'
 y=0/0.
else
 y=a(i)*mw+b(i)*x-log(x+c(i)*exp(d(i)*mw))+c1(i)
endif
!print *,i,pga_theor
end subroutine

SUBROUTINE PGA_INIT()
use mod_pgamisf
implicit none
integer :: i,j,ierr
real :: dum,dum1

nper=0
!read coefficients
open(unit=1111,file='coeff1.dat',status='old')
 do while (ierr==0)
  read(1111,*,iostat=ierr) dum
  if (ierr==0) nper=nper+1
 enddo
!allocate
allocate(per(nper),a(nper),b(nper),c(nper),d(nper),c1(nper),tau(nper))
rewind(1111)
 do i=1,nper
  read(1111,*) per(i),a(i),b(i),c(i),d(i)
 enddo
close(1111)

open(unit=1112,file='coeff2.dat',status='old')
 do i=1,nper
  read(1112,*) dum1,dum,c1(i),(dum, j=1,4),tau(i)
  if (abs(dum1-per(i))>eps) then 
!   read(1112,*) 
!  else
   print *,'files with coefficients for PGA inconsistent',dum,per(i)
   stop
  endif
 enddo
close(1112)

allocate(sd(nper),sv(nper),sa1(nper),sa2(nper),psv(nper),psa(nper))
END SUBROUTINE



!-----------------------------------------------------------------------
      SUBROUTINE PCN05(NMX,N,NPMX,NP,DT,DAMP,P,ACC,XSD,XSV,XSA,XPSV,XPSA)
!-----------------------------------------------------------------------
!      ORIGINAL ONLY GAVE SA. MODIFIED BY JGA TO GIVE SD,SV,SA,PSV.PSA
!      new input:
!      NMX:   dimension of array containing acceleration trace
!      N:     number of points of acceleration trace
!      NPMX:  dimension of array containing frequencies of SDF oscillator
!      NP:    number of frequencies of SDF oscillator
!      DT:    delta t of acceleration trace
!      DAMP:  damping
!      P:     periods
!      ACC:   acceleration trace
!      output:
!      XSD:   spectral displacement
!      XSV:   spectral velocity
!      XSA:   spectral acceleration
!      XPSV:  pseudo velocity
!      XPSA:  pseudo acceleration
!-----------------------------------------------------------------------
      DIMENSION A(2,2),B(2,2),ACC(NMX),P(NPMX)
      DIMENSION XSD(NPMX),XSV(NPMX),XSA(NPMX),XPSV(NPMX)
      DIMENSION XPSA(NPMX),XFS(NPMX)
!-----------------------------------------------------------------------
      DO 10 IP=1,NP
         W=    6.2832/P(IP)
         DELT= P(IP)/10.
         L=    DT/DELT+1.0-1.E-05
         VERTL=1.0/FLOAT(L)
         DELT= DT*VERTL
         CALL PCN04(DAMP,W,DELT,A,B)
         XIP=    0.0
         XIPD=   0.0
         XSA(IP)=0.0
         XSV(IP)=0.0
         XSD(IP)=0.0
         NN=N-1
         DO 1 J=1,NN
            AI=ACC(J)
            SL=(ACC(J+1)-ACC(J))*VERTL
            DO 2 JJ=1,L
               AF=   AI+SL
               XIP1= XIP
               XIPD1=XIPD
               XIP=  A(1,1)*XIP1+A(1,2)*XIPD1+B(1,1)*AI+B(1,2)*AF
               XIPD= A(2,1)*XIP1+A(2,2)*XIPD1+B(2,1)*AI+B(2,2)*AF
               AI=   AF
               ABSOL=ABS(2.*DAMP*W*XIPD+W*W*XIP)
               IF(ABS(XIP).GT.XSD(IP))  XSD(IP)=ABS(XIP)
               IF(ABS(XIPD).GT.XSV(IP)) XSV(IP)=ABS(XIPD)
               IF(ABSOL.GT.XSA(IP))     XSA(IP)=ABSOL
2              CONTINUE
1           CONTINUE
         XPSV(IP)= XSD(IP)*W
         XPSA(IP)= XSV(IP)*W
         XFS(IP)=  SQRT((W*XIP)**2+XIPD**2)
10       CONTINUE
      RETURN
      END
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE PCN04(D,W,DELT,A,B)
!-----------------------------------------------------------------------
!     called by pcn05
!     D:     damping
!     W:     frequency omega
!     DELT:  something like integration step (?)
!     A:     2x2 matrix, returned
!     B:     2x2 matrix, returned
!-----------------------------------------------------------------------
      DIMENSION A(2,2),B(2,2)
!-----------------------------------------------------------------------
      DW=D*W
      D2=D*D
      A0=EXP(-DW*DELT)
      A1=W*SQRT(1.-D2)
      AD1=A1*DELT
      A2=SIN(AD1)
      A3=COS(AD1)
      A7=1.0/(W*W)
      A4=(2.0*D2-1.0)*A7
      A5=D/W
      A6= 2.0*A5*A7
      A8=1.0/A1
      A9=-(A1*A2+DW*A3)*A0
      A10=(A3-DW*A2*A8)*A0
      A11=A2*A8
      A12=A11*A0
      A13=A0*A3
      A14=A10*A4
      A15=A12*A4
      A16=A6*A13
      A17=A9*A6
      A(1,1)=A0*(DW*A11+A3)
      A(1,2)=A12
      A(2,1)=A10*DW+A9
      A(2,2)=A10
      DINV=1.0/DELT
      B(1,1)=(-A15-A16+A6)*DINV-A12*A5-A7*A13
      B(1,2)=(A15+A16-A6)*DINV+A7
      B(2,1)=(-A14-A17-A7)*DINV-A10*A5-A9*A7
      B(2,2)=(A14+A17+A7)*DINV
      RETURN
      END
!-----------------------------------------------------------------------
