! Forward (negative exponent) and backward (positive exponent) FFT by Numerical Recipes
! gfortran -O3 mFFT_NR.f90
! ifort -O3 -heap-arrays mFFT_NR.f90
! fastest: call four1(real,int,int); call fourn(real,int,int,int)
! slower:  call four1(complex,int,int); call fourn(complex,int,int,int)
! operators: .FFTn., .FFTp. (alternatives to four1)

MODULE mFFT_NR
IMPLICIT NONE

INTERFACE OPERATOR (.FFTn.)
MODULE PROCEDURE FFTn_r,FFTn_c
END INTERFACE

INTERFACE OPERATOR (.FFTp.)
MODULE PROCEDURE FFTp_r,FFTP_c
END INTERFACE

INTERFACE four1
MODULE PROCEDURE four1_r,four1_c
END INTERFACE

INTERFACE fourn
MODULE PROCEDURE fourn_r,fourn_c
END INTERFACE

CONTAINS

FUNCTION FFTn_r(data) RESULT (r)
REAL(8),INTENT(IN) :: data(:)
REAL(8) r(size(data))
r=data
call four1(r,size(r)/2,-1)
END FUNCTION

FUNCTION FFTn_c(data) RESULT (r)
COMPLEX(8),INTENT(IN) :: data(:)
COMPLEX(8) r(size(data))
r=data
call four1(r,size(r),-1)
END FUNCTION

FUNCTION FFTp_r(data) RESULT (r)
REAL(8),INTENT(IN) :: data(:)
REAL(8) r(size(data))
r=data
call four1(r,size(r)/2,1)
END FUNCTION

FUNCTION FFTp_c(data) RESULT (r)
COMPLEX(8),INTENT(IN) :: data(:)
COMPLEX(8) r(size(data))
r=data
call four1(r,size(r),1)
END FUNCTION

SUBROUTINE four1_r(data,nn,isign)
INTEGER isign,nn
REAL(8) data(2*nn)
INTEGER i,istep,j,m,mmax,n
REAL(8) tempi,tempr
REAL(8) theta,wi,wpi,wpr,wr,wtemp
n=2*nn
j=1
do i=1,n,2
  if (j.gt.i)then
    tempr=data(j)
    tempi=data(j+1)
    data(j)=data(i)
    data(j+1)=data(i+1)
    data(i)=tempr
    data(i+1)=tempi
  endif
  m=n/2
  do while ((m.ge.2).and.(j.gt.m))
    j=j-m
    m=m/2
  enddo
  j=j+m
enddo
mmax=2
do while (n.gt.mmax)
  istep=2*mmax
  theta=6.28318530717959d0/(isign*mmax)
  wpr=-2.d0*sin(0.5d0*theta)**2
  wpi=sin(theta)
  wr=1.d0
  wi=0.d0
  do m=1,mmax,2
    do i=m,n,istep
      j=i+mmax
      tempr=dble(wr)*data(j)-dble(wi)*data(j+1)
      tempi=dble(wr)*data(j+1)+dble(wi)*data(j)
      data(j)=data(i)-tempr
      data(j+1)=data(i+1)-tempi
      data(i)=data(i)+tempr
      data(i+1)=data(i+1)+tempi
    enddo
    wtemp=wr
    wr=wr*wpr-wi*wpi+wr
    wi=wi*wpr+wtemp*wpi+wi
  enddo
  mmax=istep
enddo
return
END SUBROUTINE

SUBROUTINE four1_c(data,nn,isign)
INTEGER isign,nn
COMPLEX(8) data(nn)
REAL(8) data_r(2*nn)
data_r=transfer(data,data_r)
call four1_r(data_r,nn,isign)
data=transfer(data_r,data)
END SUBROUTINE

SUBROUTINE fourn_r(data,nn,ndim,isign)
INTEGER isign,ndim,nn(ndim)
REAL(8) data(2*product(nn))
INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
REAL(8) tempi,tempr
REAL(8) theta,wi,wpi,wpr,wr,wtemp
ntot=1
do idim=1,ndim
  ntot=ntot*nn(idim)
enddo
nprev=1
do idim=1,ndim
  n=nn(idim)
  nrem=ntot/(n*nprev)
  ip1=2*nprev
  ip2=ip1*n
  ip3=ip2*nrem
  i2rev=1
  do i2=1,ip2,ip1
    if (i2.lt.i2rev)then
      do i1=i2,i2+ip1-2,2
        do i3=i1,ip3,ip2
          i3rev=i2rev+i3-i2
          tempr=data(i3)
          tempi=data(i3+1)
          data(i3)=data(i3rev)
          data(i3+1)=data(i3rev+1)
          data(i3rev)=tempr
          data(i3rev+1)=tempi
        enddo
      enddo
    endif
    ibit=ip2/2
    do while ((ibit.ge.ip1).and.(i2rev.gt.ibit))
      i2rev=i2rev-ibit
      ibit=ibit/2
    enddo
    i2rev=i2rev+ibit
  enddo
  ifp1=ip1
  do while (ifp1.lt.ip2)
    ifp2=2*ifp1
    theta=isign*6.28318530717959d0/(ifp2/ip1)
    wpr=-2.d0*sin(0.5d0*theta)**2
    wpi=sin(theta)
    wr=1.d0
    wi=0.d0
    do i3=1,ifp1,ip1
      do i1=i3,i3+ip1-2,2
        do i2=i1,ip3,ifp2
          k1=i2
          k2=k1+ifp1
          tempr=dble(wr)*data(k2)-dble(wi)*data(k2+1)
          tempi=dble(wr)*data(k2+1)+dble(wi)*data(k2)
          data(k2)=data(k1)-tempr
          data(k2+1)=data(k1+1)-tempi
          data(k1)=data(k1)+tempr
          data(k1+1)=data(k1+1)+tempi
        enddo
      enddo
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
    enddo
    ifp1=ifp2
  enddo
  nprev=n*nprev
enddo
return
END SUBROUTINE

SUBROUTINE fourn_c(data,nn,ndim,isign)
INTEGER isign,ndim,nn(ndim)
COMPLEX(8) data(product(nn))
REAL(8) data_r(2*product(nn))
data_r=transfer(data,data_r)
call fourn_r(data_r,nn,ndim,isign)
data=transfer(data_r,data)
END SUBROUTINE

SUBROUTINE init_fourn(data,nn,ndim)
INTEGER ndim,nn(ndim)
COMPLEX(8) data(product(nn))
END SUBROUTINE

SUBROUTINE destroy_fourn
END SUBROUTINE

SUBROUTINE init_FFTW
END SUBROUTINE

SUBROUTINE cleanup_FFTW
END SUBROUTINE

END MODULE mFFT_NR
