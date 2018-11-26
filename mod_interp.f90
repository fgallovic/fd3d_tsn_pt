SUBROUTINE interpolate(xintrpl,m,yintrpl,n,xnew,nx,ynew,ny,tintrpl,D)
implicit none
real :: yintrpl(n),xintrpl(m),xnew(nx),ynew(ny),tintrpl(m,n),D(nx,ny)
integer :: i,j,m,n,nx,ny
real, allocatable :: yspln(:), xspln(:),d1(:,:)
real :: dum

allocate(yspln(n),xspln(m),d1(m,ny))
yspln=0.
xspln=0.
d1=0.
do i=1,M
  CALL spline(yintrpl,tintrpl(i,:),N,1.e30,1.e30,yspln)
  do j=1,NY
    dum=ynew(j)
    call splint(yintrpl,tintrpl(i,:),yspln,N,dum,D1(i,j))
  enddo
enddo
do j=1,NY
  CALL spline(xintrpl,D1(:,j),M,1.e30,1.e30,xspln)
  do i=1,NX
    dum=xnew(i)
    call splint(xintrpl,D1(:,j),xspln,M,dum,D(i,j))
  enddo
enddo
deallocate(yspln,xspln,d1)
end subroutine interpolate


      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n
      REAL yp1,ypn,x(n),y(n),y2(n)
      INTEGER i,k
      REAL p,qn,sig,un
      REAL,ALLOCATABLE:: u(:)
      allocate(u(n))
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      deallocate(u)
      return
      END

      
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
!1     if (khi-klo.gt.1) then
      do while(khi-klo.gt.1)
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
!      goto 1
!      endif
      enddo
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END

