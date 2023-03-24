module mod_pgamisf
integer :: GMPE_id !1 for zhao, 2 for boore
integer :: nper_gmpe !number of periods with coefficients
!for Zhao GMPEs:
real, allocatable :: aZhao(:),bZhao(:),cZhao(:),dZhao(:),c0Zhao(:),stZhao(:),Zhao_sigma(:),Zhao_tau(:)
!for Boore GMPEs:
 integer, PARAMETER ::  ncoeff_stage1 = 3, ncoeff_stage2 = 7, nregions = 3
 real,allocatable :: e_gmpe(:,:), amh_gmpe(:), c_gmpe(:,:), amref_gmpe(:), rref_gmpe(:), h_gmpe(:)     
 real, allocatable ::  delta_c3(:,:), clin(:), vclin(:), vref(:), f1(:), f3(:), f4(:), f5(:), f6(:), f7(:) 
 real,allocatable ::  r1(:), r2(:), delta_phiR(:), delta_phiV(:), v1(:), v2(:), phi1(:), phi2(:), tau1(:), tau2(:)
 real, allocatable :: im_per(:),a_cita(:),b1_cita(:),b2_cita(:),c1_cita(:),c2_cita(:),c3_cita(:),Mref_cita(:),h_cita(:),sigma_cita(:),tau_cita(:),phi_s2s_cita(:),tau_l2l_cita(:),phi_p2p_cita(:),phi0_cita(:)
 real, allocatable :: a_ita18(:),b1_ita18(:),b2_ita18(:),c1_ita18(:),c2_ita18(:),c3_ita18(:),Mref_ita18(:),h_ita18(:),k_ita18(:),tau_ita18(:),phi_s2s_ita18(:),phi0_ita18(:),mh_ita18(:)
 real :: rref, mh,v0a
 integer :: indx_pga,indx_permax,indx_pgv
 real :: v30
!!!
real,allocatable :: per_gmpe(:) !periods where gmpes coeffs are defined
real, allocatable :: sd(:),sv(:),sa1(:),sa2(:),psv(:),psa(:)
!real :: pga_theor
real :: damp=0.05
real :: eps=1.e-6

contains
      FUNCTION select(k,n,arr)
      implicit none
      INTEGER k,n
      REAL select,arr(n)
      INTEGER i,ir,j,l,mid
      REAL a,temp
      l=1
      ir=n
1     if(ir-l.le.1)then
        if(ir-l.eq.1)then
          if(arr(ir).lt.arr(l))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
          endif
        endif
        select=arr(k)
        return
      else
        mid=(l+ir)/2
        temp=arr(mid)
        arr(mid)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        if(j.ge.k)ir=j-1
        if(j.le.k)l=i
      endif
      goto 1
      END function
end module

subroutine pga_theor(x,mw,perT,y,sigma,tau) !result(y)
use mod_pgamisf
implicit none
integer :: i
real, intent(in) :: x,mw,perT
real, intent(out) :: y,sigma,tau 
!parameters for boore gmpes:
real :: rjb, z1=-1, per_desired
integer :: iregion=0,irelation=1,mech=1
real :: r_pga4nl,per1,per2,r_t1,r_t2,slope_logy,slope_phi,slope_tau,slope_sigma,pga4nl,fe_pga4nl,fpb_pga4nl,fp_pga4nl,gspread_pga4nl
real :: y_t1,fe_t1,fpb_t1,fsb_t1,fs_t1,phi_t1,tau_t1,sigma_t1,gspread_t1,amp_lin_t1,amp_nl_t1,amp_total_t1,fp_t1
real :: y_t2,fe_t2,fpb_t2,fsb_t2,fs_t2,phi_t2,tau_t2,sigma_t2,gspread_t2,amp_lin_t2,amp_nl_t2,amp_total_t2,fp_t2
real :: y1,y2,p1,p2,t1,t2,s1,s2,pp,ll,l1,l2,peso1,peso2,tt,ss,ypred
integer :: ind1, ind2,iper,clust_id,indx
real :: yg,ycgs,phi
integer :: indx1,indx2,iflag_per

select case(GMPE_id)
case(1) !ZHAO
i=1
!print *,i,per1
do while (abs(perT-per_gmpe(i))>eps .and. i<=nper_gmpe)
 i=i+1
enddo
!print *,i
if (i>nper_gmpe) then
 print *,'period inconsistent'
 y=-99999.!1/0.
else !so far no interpolation between periods...
 y=aZhao(i)*mw+bZhao(i)*x-log(x+cZhao(i)*exp(dZhao(i)*mw))+c0Zhao(i)
 sigma=Zhao_sigma(i)
 tau=Zhao_tau(i)
endif
!print *,i,pga_theor

case(2) !BOORE
per_desired=perT; rjb=x; mech=1 !others taken from defs
!check for period consistency:
! Trap for per outside the range of the tabulated periods (but only if not pga or pgv):

            if (per_desired < per_gmpe(1) .or. per_desired > per_gmpe(indx_permax))   then
              print *,' ERROR: per = ', per_desired, ' is outside the tabulated range in per_gmpe table'
              print *,' per_gmpe '
              do i = 1, nper_gmpe
                print *, i, per_gmpe(i)
              end do
              print *, ' QUITTING!!!'
              stop
            end if
! Find lower index of the two tabulated periods that span the interval
! containing the input period
! assume that the periods are sorted from small to large, and that period value
! for index nper_gmpe corresponds to the maximum period (and thus 
! the indices 1 and 2 are for pgv and pgv (per_gmpe = -1 and 0)).


            iflag_per = 0
            IF (per_desired == per_gmpe(nper_gmpe)) THEN
              indx1 = nper_gmpe-1
              indx2 = nper_gmpe
              iflag_per = 1
            ELSE
              DO i = 1, nper_gmpe-1  
                IF (per_desired >= per_gmpe(i) .and. per_desired < per_gmpe(i+1) ) THEN
                  indx1 = i
                  indx2 = i+1
                  iflag_per = 1
                  EXIT
                END IF
              END DO
            END IF
            if (iflag_per == 0) then
              PRINT *,' ERROR: could not find per = ',per_desired,' in per_gmpe table'
              PRINT *,' per_gmpe '
              DO i = 1, nper_gmpe
                PRINT *, per_gmpe(i)
              END DO
              PRINT *, ' QUITTING!!!'
              STOP
            end if
  
            per1 = per_gmpe(indx1)
            per2 = per_gmpe(indx2)
! Now evaluate PSA and sigma for the periods on either side of the desired period 
! (the evaluations are so quick that there is no need to treat
! the case of the desired and tabulated periods being the same).

            r_t1 = sqrt(rjb**2+h_gmpe(indx1)**2)
      
            call bssa14_gm_sub4y(per_desired, mw, rjb, r_t1, v30, mech, iregion, z1, irelation, pga4nl,c_gmpe(:,indx1), amref_gmpe(indx1), rref_gmpe(indx1), delta_c3(:,indx1),&
                                 e_gmpe(:,indx1), amh_gmpe(indx1), clin(indx1), vclin(indx1), vref(indx1), f1(indx1), f3(indx1), f4(indx1), f5(indx1), f6(indx1), f7(indx1), &
                                 phi1(indx1), phi2(indx1),  tau1(indx1), tau2(indx1), R1(indx1), R2(indx1), delta_phiR(indx1), v1(indx1), v2(indx1), delta_phiV(indx1), y_t1, &
                                 fe_t1, fpb_t1, fp_t1, fsb_t1, fs_t1, phi_t1, tau_t1, sigma_t1, gspread_t1, amp_lin_t1, amp_nl_t1, amp_total_t1)    
     
     
            r_t2 = sqrt(rjb**2+h_gmpe(indx2)**2)
      
            call bssa14_gm_sub4y(per_desired, mw, rjb, r_t2, v30, mech, iregion, z1, irelation, pga4nl,c_gmpe(:,indx2), amref_gmpe(indx2), rref_gmpe(indx2), delta_c3(:,indx2), &
                                 e_gmpe(:,indx2), amh_gmpe(indx2), clin(indx2), vclin(indx2), vref(indx2), f1(indx2), f3(indx2), f4(indx2), f5(indx2), f6(indx2), f7(indx2), &
                                 phi1(indx2), phi2(indx2), tau1(indx2), tau2(indx2), R1(indx2), R2(indx2), delta_phiR(indx2), v1(indx2), v2(indx2), delta_phiV(indx2), y_t2,&
                                 fe_t2, fpb_t2, fp_t2, fsb_t2, fs_t2, phi_t2, tau_t2, sigma_t2, gspread_t2, amp_lin_t2, amp_nl_t2, amp_total_t2)     

            slope_logy =  (alog10(y_t2) - alog10(y_t1))/alog10(per2/per1)
            y = 10**(alog10(y_t1) +  slope_logy*alog10(per_desired/per1)) 
        
            slope_phi =  (phi_t2 - phi_t1)/alog10(per2/per1); phi = phi_t1 + slope_phi*alog10(per_desired/per1) 
            slope_tau =  (tau_t2 - tau_t1)/alog10(per2/per1); tau = tau_t1 + slope_tau*alog10(per_desired/per1)

!....we do not need the other parts by now...
 
!            slope_sigma =  (sigma_t2 - sigma_t1)/ alog10(per2/per1); sigma = sigma_t1 + slope_sigma*alog10(per_desired/per1)

 
!            slope_fe =  (fe_t2 - fe_t1)/alog10(per2/per1); fe =  fe_t1 + slope_fe*alog10(per_desired/per1) 
!            slope_fpb =  (fpb_t2 - fpb_t1)/alog10(per2/per1) ; fpb =  fpb_t1 + slope_fpb*alog10(per_desired/per1) 
!            slope_fp =  (fp_t2 - fp_t1)/alog10(per2/per1); fp =  fp_t1 + slope_fp*alog10(per_desired/per1) 
!            slope_fsb =  (fsb_t2 - fsb_t1)/alog10(per2/per1); fsb =  fsb_t1 + slope_fsb*alog10(per_desired/per1) 
!            slope_fs =  (fs_t2 - fs_t1)/alog10(per2/per1); fs =  fs_t1 + slope_fs*alog10(per_desired/per1) 
!            slope_gspread =  (gspread_t2 - gspread_t1)/alog10(per2/per1); gspread =  gspread_t1 + slope_gspread*alog10(per_desired/per1) 
!            slope_logamp_lin = (alog10(amp_lin_t2)-alog10(amp_lin_t1))/alog10(per2/per1); amp_lin = 10**(alog10(amp_lin_t1) + slope_logamp_lin*alog10(per_desired/per1)) 
!            slope_logamp_nl =  (alog10(amp_nl_t2) - alog10(amp_nl_t1))/alog10(per2/per1); amp_nl = 10**(alog10(amp_nl_t1) + slope_logamp_nl*alog10(per_desired/per1)) 
!            slope_logamp_total =  (alog10(amp_total_t2) - alog10(amp_total_t1))/alog10(per2/per1); amp_total = 10**(alog10(amp_total_t1) + slope_logamp_total*alog10(per_desired/per1)) 
!            slope_logr =  (alog10(r_t2) - alog10(r_t1))/alog10(per2/per1); r = 10**(alog10(r_t1) + slope_logr*alog10(per_desired/per1)) 

!CONVERT from g to cgs
        if (per_desired >= 0.0) then
          yg   =y
          ycgs = 981.0 * yg
        else
          yg   = y 
          ycgs = yg
        end if
        sigma=phi
        y=log(ycgs) !so that definition same with zhao: y in cm/s-2 and natural logarithm..

case(3) !ITALIAN GMPE
      !find if perT directly in table:
      if (perT>im_per(nper_gmpe)) then
          print *,'error! period ',perT, 'outside range. Max period:', im_per(nper_gmpe)
          stop
      endif
      indx=0
      do i=1,nper_gmpe
        if (abs(perT-im_per(i))<1.e-8) then
           indx=i
        endif
      enddo
      if (indx>0) then
       ypred=a_cita(indx)
       if (mw<=mh) then
           ypred=ypred+b1_cita(indx)*(mw-mh)
       else
           ypred=ypred+b2_cita(indx)*(mw-mh)
       endif
       ypred=ypred+(c1_cita(indx)*(mw-mref_cita(indx))+c2_cita(indx))*log10(sqrt(x**2+(h_cita(indx))**2)/Rref)+c3_cita(indx)*(sqrt(x**2+(h_cita(indx))**2)-Rref)
       sigma=phi0_cita(indx)**2!sigma_cita(indx)**2
       tau=tau_cita(indx)**2
       !phi_s2s not accounted
       !phi_p2p accounted - contains hanging wall and footwall dependency
       sigma=sigma+phi_p2p_cita(indx)**2
       !tau_l2l add to tau if not corrected for clusters:
       if (clust_id>0) then
          !!!CLUSTER CORRECTION ADD HERE
       else
         tau=tau+tau_l2l_cita(indx)**2
       endif
       !sqrt of sigmas and taus
       sigma=sqrt(sigma)
       tau=sqrt(tau)
 
      else
       do iper=1,nper_gmpe-1
        if (perT>im_per(iper) .and. perT<im_per(iper+1)) then
           ind1=iper
           ind2=iper+1
        endif
       enddo
       y1=a_cita(ind1)
       if (mw<=mh) then
           y1=y1+b1_cita(ind1)*(mw-mh)
       else
           y1=y1+b2_cita(ind1)*(mw-mh)
       endif
       y1=y1+(c1_cita(ind1)*(mw-mref_cita(ind1))+c2_cita(ind1))*log10(sqrt(x**2+(h_cita(ind1))**2)/Rref)+c3_cita(ind1)*(sqrt(x**2+(h_cita(ind1))**2)-Rref)
       s1=phi0_cita(ind1)!sigma_cita(ind1)
       t1=tau_cita(ind1)
       p1=phi_p2p_cita(ind1)
       l1=tau_l2l_cita(ind1)
       peso1=abs(im_per(ind1)-perT)
       
       y2=a_cita(ind2)
       if (mw<=mh) then
           y2=y2+b1_cita(ind2)*(mw-mh)
       else
           y2=y2+b2_cita(ind2)*(mw-mh)
       endif
       y2=y2+(c1_cita(ind2)*(mw-mref_cita(ind2))+c2_cita(ind2))*log10(sqrt(x**2+(h_cita(ind2))**2)/Rref)+c3_cita(ind2)*(sqrt(x**2+(h_cita(ind2))**2)-Rref)
       s2=phi0_cita(ind2)!sigma_cita(ind2)
       t2=tau_cita(ind2)
       p2=phi_p2p_cita(ind2)
       l2=tau_l2l_cita(ind2)
       peso2=abs(im_per(ind2)-perT)

       ypred=(y1*peso2+y2*peso1)/(peso1+peso2)
       ss=(s1*peso2+s2*peso1)/(peso1+peso2)
       tt=(t1*peso2+t2*peso1)/(peso1+peso2)
       pp=(p1*peso2+p2*peso1)/(peso1+peso2)
       ll=(l1*peso2+l2*peso1)/(peso1+peso2)
       sigma=sqrt(ss**2+pp**2)
       !
       if (clust_id>0) then
          !!!CLUSTER CORRECTION ADD HERE
       else !increase variance
         tau=sqrt(tt**2+ll**2)
       endif
      endif
      !CHECK JEDNOTKY, tu: ypred aj sigma je v log10, cm/s^2
      y=ypred*log(10.)
      tau=tau*log(10.)
      sigma=sigma*log(10.)


case(4,5) !ITALIAN ITA18 GMPE by Lanzano et al.
      !find if perT directly in table:
      if (perT>im_per(nper_gmpe)) then
          print *,'error! period ',perT, 'outside range. Max period:', im_per(nper_gmpe)
          stop
      endif
      indx=0
      do i=1,nper_gmpe
        if (abs(perT-im_per(i))<1.e-8) then
           indx=i
        endif
      enddo
      if (indx>0) then
       ypred=a_ita18(indx)
       if (mw<=mh_ita18(indx)) then
           ypred=ypred+b1_ita18(indx)*(mw-mh_ita18(indx))
       else
           ypred=ypred+b2_ita18(indx)*(mw-mh_ita18(indx))
       endif
       ypred=ypred+(c1_ita18(indx)*(mw-mref_ita18(indx))+c2_ita18(indx))*log10(sqrt(x**2+(h_ita18(indx))**2))+c3_ita18(indx)*(sqrt(x**2+(h_ita18(indx))**2))
       v0a=min(v30,1500.)
       ypred=ypred+k_ita18(indx)*log10(v0a/800.)
       sigma=phi0_ita18(indx)**2!sigma_cita(indx)**2
       tau=tau_ita18(indx)**2
       !phi_s2s not accounted
       !phi_p2p accounted - contains hanging wall and footwall dependency
       !sigma=sigma+phi_p2p_ita18(indx)**2
       !tau_l2l add to tau if not corrected for clusters:
       !sqrt of sigmas and taus
       sigma=sqrt(sigma)
       tau=sqrt(tau)
 
      else
       do iper=1,nper_gmpe-1
        if (perT>im_per(iper) .and. perT<im_per(iper+1)) then
           ind1=iper
           ind2=iper+1
        endif
       enddo
       y1=a_ita18(ind1)
       if (mw<=mh_ita18(ind1)) then
           y1=y1+b1_ita18(ind1)*(mw-mh_ita18(ind1))
       else
           y1=y1+b2_ita18(ind1)*(mw-mh_ita18(ind1))
       endif
       y1=y1+(c1_ita18(ind1)*(mw-mref_ita18(ind1))+c2_ita18(ind1))*log10(sqrt(x**2+(h_ita18(ind1))**2))+c3_ita18(ind1)*(sqrt(x**2+(h_ita18(ind1))**2))
       v0a=min(v30,1500.)
       y1=y1+k_ita18(indx)*log10(v0a/800.)
       s1=phi0_ita18(ind1)!sigma_cita(ind1)
       t1=tau_ita18(ind1)
       !p1=phi_p2p_ita18(ind1)
       peso1=abs(im_per(ind1)-perT)
       
       y2=a_ita18(ind2)
       if (mw<=mh_ita18(ind2)) then
           y2=y2+b1_ita18(ind2)*(mw-mh_ita18(ind2))
       else
           y2=y2+b2_ita18(ind2)*(mw-mh_ita18(ind2))
       endif
       y2=y2+(c1_ita18(ind2)*(mw-mref_ita18(ind2))+c2_ita18(ind2))*log10(sqrt(x**2+(h_ita18(ind2))**2))+c3_ita18(ind2)*(sqrt(x**2+(h_ita18(ind2))**2))
       y2=y2+k_ita18(ind2)*log10(v0a/800.)
       s2=phi0_ita18(ind2)!sigma_cita(ind2)
       t2=tau_ita18(ind2)
       !p2=phi_p2p_ita18(ind2)
       peso2=abs(im_per(ind2)-perT)

       ypred=(y1*peso2+y2*peso1)/(peso1+peso2)
       ss=(s1*peso2+s2*peso1)/(peso1+peso2)
       tt=(t1*peso2+t2*peso1)/(peso1+peso2)
       pp=0.!(p1*peso2+p2*peso1)/(peso1+peso2)
       sigma=sqrt(ss**2+pp**2)
       !
      endif
     !CHECK JEDNOTKY, tu: ypred aj sigma je v log10, cm/s^2
      y=ypred*log(10.)
      tau=tau*log(10.)
      sigma=sigma*log(10.)

end select 

end subroutine

SUBROUTINE PGA_INIT()
use mod_pgamisf
use waveforms_com, only: nper
implicit none
integer :: i,j,ierr
real :: dum,dum1,per_max
integer :: i_read_status
character(300) :: f_coeff
real, allocatable :: coeff(:,:)
integer :: nu_coeff

!read coefficients

select case(GMPE_id)
case(1) !ZHAO
ierr=0
nper_gmpe=0
open(unit=1111,file='coeff1.dat',status='old')
 do while (ierr==0)
  read(1111,*,iostat=ierr) dum
  if (ierr==0) nper_gmpe=nper_gmpe+1
 enddo
!allocate
allocate(per_gmpe(nper_gmpe),aZhao(nper_gmpe),bZhao(nper_gmpe),cZhao(nper_gmpe),dZhao(nper_gmpe),c0Zhao(nper_gmpe),stZhao(nper_gmpe),Zhao_sigma(nper_gmpe),Zhao_tau(nper_gmpe))
rewind(1111)
 do i=1,nper_gmpe
  read(1111,*) per_gmpe(i),aZhao(i),bZhao(i),cZhao(i),dZhao(i) 
 enddo
close(1111)

open(unit=1112,file='coeff2.dat',status='old')
 do i=1,nper_gmpe
  read(1112,*) dum1, c0Zhao(i),dum,(dum, j=1,3),Zhao_sigma(i),Zhao_tau(i),stZhao(i) 
  if (abs(dum1-per_gmpe(i))>eps) then 
!   read(1112,*) 
!  else
   print *,'files with coefficients for ZHAO gmpes inconsistent',dum,per_gmpe(i)
   stop
  endif
 enddo
close(1112)
allocate(sd(nper),sv(nper),sa1(nper),sa2(nper),psv(nper),psa(nper))

case(2) !BOORE
nu_coeff=1111
f_coeff='BSSA14_Coefficients_071314_Revisedf4_071514.csv'
!STOLEN:
! Open coefficient files and
! read coefficients for regression coefficients:
      open(unit=nu_coeff,file=trim(f_coeff), status='old')

      read(nu_coeff,*)
      i = 0
      
      read_coefficients: DO

!Coefficient file:
! 1	         2	3	4	5	6	7	8	9	10	11	12	13	14	          15	        16	         17	                 18	         19	     20	      21	22	23	24	25	26	        27	   28	29	30	31	32	33	34	35	36	37
!Period(sec)	e0	e1	e2	e3	e4	e5	e6	Mh	c1	c2	c3	Mref	Rref(km)	h(km)	Dc3(globalCATWNZ)	Dc3(ChinaTurkey)	Dc3(ItalyJapan)	  c	Vc(m/s)	Vref(m/s)	f1	f3	f4	f5	f6(1/km)	f7	R1(km)	R2(km)	DfR	DfV	V1(m/s)	V2(m/s)	f1	f2	t1	t2
!period	        e0	e1	e2	e3	e4	e5	e6	Mh	c1	c2	c3	Mref	Rref	        h	Dc3.globalCATWNZ	Dc3.ChinaTurkey	        Dc3.ItalyJapan	clin	Vc	Vref	        f1	f3	f4	f5	f6	        f7	R1	R2	DfR	DfV	V1	V2	f1	f2	t1	t2
!     :                      ncoeff_stage1 = 3,
!     :                      ncoeff_stage2 = 7,
!     :                      nregions = 3

!get per_gmpe
        read(nu_coeff,*,iostat=i_read_status) dum 

        if (i_read_status /= 0) exit
        
        i = i + 1
      
      END DO read_coefficients

      nper_gmpe = i

      allocate(per_gmpe(nper_gmpe),e_gmpe(0:ncoeff_stage2-1,nper_gmpe), amh_gmpe(nper_gmpe), &
      c_gmpe(1:ncoeff_stage1,nper_gmpe),amref_gmpe(nper_gmpe),rref_gmpe(nper_gmpe),h_gmpe(nper_gmpe),&
      delta_c3(0:nregions,nper_gmpe),clin(nper_gmpe),vclin(nper_gmpe), vref(nper_gmpe),&
      f1(nper_gmpe),f3(nper_gmpe),f4(nper_gmpe),f5(nper_gmpe),f6(nper_gmpe),f7(nper_gmpe),&
      r1(nper_gmpe),r2(nper_gmpe),delta_phiR(nper_gmpe),delta_phiV(nper_gmpe), &
      v1(nper_gmpe),v2(nper_gmpe), phi1(nper_gmpe),phi2(nper_gmpe),tau1(nper_gmpe),tau2(nper_gmpe))
      delta_c3=0.
     rewind(nu_coeff)
     read(nu_coeff,*) !header
     do i=1,nper_gmpe
      read(nu_coeff,*,iostat=ierr) per_gmpe(i), (e_gmpe(j,i), j = 0, ncoeff_stage2-1), amh_gmpe(i), &
                       (c_gmpe(j,i), j = 1, ncoeff_stage1), amref_gmpe(i), rref_gmpe(i), &
                        h_gmpe(i), (delta_c3(j,i), j = 1, nregions), clin(i), vclin(i), vref(i), &
                        f1(i), f3(i), f4(i), f5(i),f6(i), f7(i), r1(i), r2(i), delta_phiR(i), &
                        delta_phiV(i), v1(i), v2(i), phi1(i), phi2(i), tau1(i), tau2(i)
      if (ierr/=0) stop
     enddo          
      close(nu_coeff)

    indx_pgv = 0
      find_pgv_index: DO i = 1, nper_gmpe
        IF (per_gmpe(i) == -1.0) THEN
          indx_pgv = i
          EXIT find_pgv_index
        END IF
      END DO find_pgv_index
         indx_pga = 0
      find_pga_index: DO i = 1, nper_gmpe
        IF (per_gmpe(i) == 0.0) THEN
          indx_pga = i
          EXIT find_pga_index
        END IF
      END DO find_pga_index

      per_max = -1.0
      DO i = 1, nper_gmpe
        IF (per_gmpe(i) > per_max) THEN
          per_max = per_gmpe(i)
          indx_permax = i
        END IF
      END DO
case(3) !ITALIAN GMPE
     Rref=1
     Mh=5

     ierr=0
     nper_gmpe=0
     open(1111,file='neclust_coeff.dat',status='old')
     do while (ierr==0)
        read(1111,*,iostat=ierr) dum
        if (ierr==0) nper_gmpe=nper_gmpe+1
     enddo
     allocate(coeff(nper_gmpe,16))
     rewind(1111)
      do i=1,nper_gmpe
        read(1111,*) coeff(i,:)
      enddo
     close(1111)
     allocate(im_per(nper_gmpe),a_cita(nper_gmpe),b1_cita(nper_gmpe),b2_cita(nper_gmpe),c1_cita(nper_gmpe),c2_cita(nper_gmpe),c3_cita(nper_gmpe),Mref_cita(nper_gmpe),h_cita(nper_gmpe),sigma_cita(nper_gmpe),tau_cita(nper_gmpe),phi_s2s_cita(nper_gmpe),tau_l2l_cita(nper_gmpe),phi_p2p_cita(nper_gmpe),phi0_cita(nper_gmpe))
      im_per=coeff(:,1)
      a_cita=coeff(:,2)
      b1_cita=coeff(:,3)
      b2_cita=coeff(:,4)
      c1_cita=coeff(:,5)
      c2_cita=coeff(:,6)
      c3_cita=coeff(:,7)
      Mref_cita=coeff(:,9)
      h_cita=coeff(:,10)
      sigma_cita=coeff(:,16)
      tau_cita=coeff(:,11)
      phi_s2s_cita=coeff(:,12)
      tau_l2l_cita=coeff(:,13)
      phi_p2p_cita=coeff(:,14)
      phi0_cita=coeff(:,15)
      !GEOMEtriCAL MEAN oF HROIZONTAL COMP
      allocate(sd(nper),sv(nper),sa1(nper),sa2(nper),psv(nper),psa(nper))

case(4) !Italian GMPE ITA18 Lanzano 2019
     open(1111,file='coeff_ita18_JB.dat',status='old')
     do while (ierr==0)
        read(1111,*,iostat=ierr) dum
        if (ierr==0) nper_gmpe=nper_gmpe+1
     enddo
     allocate(coeff(nper_gmpe,16))
     rewind(1111)
      do i=1,nper_gmpe
        read(1111,*) coeff(i,:)
      enddo
     close(1111)
     allocate(im_per(nper_gmpe),a_ita18(nper_gmpe),b1_ita18(nper_gmpe),b2_ita18(nper_gmpe),c1_ita18(nper_gmpe),c2_ita18(nper_gmpe),c3_ita18(nper_gmpe),k_ita18(nper_gmpe),Mref_ita18(nper_gmpe),h_ita18(nper_gmpe),tau_ita18(nper_gmpe),phi_s2s_ita18(nper_gmpe),phi0_ita18(nper_gmpe),mh_ita18(nper_gmpe))
      im_per=coeff(:,1)
      a_ita18=coeff(:,2)
      b1_ita18=coeff(:,3)
      b2_ita18=coeff(:,4)
      c1_ita18=coeff(:,5)
      c2_ita18=coeff(:,6)
      c3_ita18=coeff(:,7)
      k_ita18=coeff(:,8)
      Mref_ita18=coeff(:,15)
      h_ita18=coeff(:,16)
      tau_ita18=coeff(:,11)
      phi_s2s_ita18=coeff(:,12)
      phi0_ita18=coeff(:,13)
       mh_ita18=coeff(:,14)
      
      !allocate(sd(nper),sv(nper),sa1(nper),sa2(nper),psv(nper),psa(nper))

case(5) !Italian GMPE ITA18 Lanzano 2019, rupture distance
     open(1111,file='coeff_ita18_rupt.dat',status='old')
     do while (ierr==0)
        read(1111,*,iostat=ierr) dum
        if (ierr==0) nper_gmpe=nper_gmpe+1
     enddo
     allocate(coeff(nper_gmpe,16))
     rewind(1111)
      do i=1,nper_gmpe
        read(1111,*) coeff(i,:)
      enddo
     close(1111)
     allocate(im_per(nper_gmpe),a_ita18(nper_gmpe),b1_ita18(nper_gmpe),b2_ita18(nper_gmpe),c1_ita18(nper_gmpe),c2_ita18(nper_gmpe),c3_ita18(nper_gmpe),k_ita18(nper_gmpe),Mref_ita18(nper_gmpe),h_ita18(nper_gmpe),tau_ita18(nper_gmpe),phi_s2s_ita18(nper_gmpe),phi0_ita18(nper_gmpe),mh_ita18(nper_gmpe))
      im_per=coeff(:,1)
      a_ita18=coeff(:,2)
      b1_ita18=coeff(:,3)
      b2_ita18=coeff(:,4)
      c1_ita18=coeff(:,5)
      c2_ita18=coeff(:,6)
      c3_ita18=coeff(:,7)
      k_ita18=coeff(:,8)
      Mref_ita18=coeff(:,15)
      h_ita18=coeff(:,16)
      tau_ita18=coeff(:,11)
      phi_s2s_ita18=coeff(:,12)
      phi0_ita18=coeff(:,13)
       mh_ita18=coeff(:,14)

end select

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
         VERTL=1.0/REAL(L)
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
! --------------------------------------- y_bssa14_no_site_amps --------------------------------
      subroutine y_bssa14_no_site_amps(m, r, mech, iregion,c, mref, rref, delta_c3,e, mh,y, fe, fpb, fp, gspread)     
     
! **** USE THIS VERSION WITH COEFFICIENTS IN PEER REPORT FORMAT ****

! NOTE: y in g, unless pgv!!!!  

! ncoeff_s1, ncoeff_s2 are not included as arguments because it is
! assumed that they are 3 and 7, respectively.

! Assume no site amp

! mech = 0, 1, 2, 3 for unspecified, SS, NS, RS, respectively

! Dates: 02/27/13 - Modified from y_ngaw2_no_site_amps
!        07/01/13 - Renamed from y_bssa13_no_site_amps
 
      IMPLICIT none
      
      real :: m, r,c(1:3), mref, rref, delta_c3(0:3), e(0:6), mh
      real ::  alny, y, fpb, fp, fe, fs, fault_type_term, gspread     
     
      integer :: mech, iregion
     

 
      if ( e(mech) == 0.0) then   ! no coefficient
        print *,' From y_bssa14_no_site_amps, '//' e(mech) == 0.0; QUITTING'
        stop
      end if
        
      if (mech /=0 .and. mech /= 1 .and. mech /= 2 .and. mech /= 3)then
        print *,' From y_bssa14_no_site_amps, mech = ', mech,' not a valid value; QUITTING'
        stop
      else
        fault_type_term = e(mech)
      end if
        

      if (m < mh ) then      
        fe = e(4)*(m-mh) + e(5)*(m-mh)**2 
      else
        fe = e(6)*(m-mh)
      end if
 

      fe = fault_type_term + fe
 
      gspread = c(1) + c(2)*(m-mref)

      fpb = gspread*alog(r/rref) + c(3)*(r-rref)

      fp = fpb + delta_c3(iregion)*(r-rref)
      
      fs = 0.0  ! No amps
         
      alny = fe + fp + fs 
      
      y = exp(alny)
      
      return
      end
! --------------------------------------- y_bssa14_no_site_amps --------------------------------
! --------------------------------------------------------- bssa14_gm_sub4y
      subroutine bssa14_gm_sub4y(per, m, rjb, r, v30, mech, iregion, z1, irelation,pga4nl,c, amref, rref, delta_c3,e, amh,clin, vclin, vref,f1, f3, f4, f5,f6, f7,phi1, phi2, tau1, tau2,r1, r2, delta_phiR,v1, v2, delta_phiV,y, fe, fpb, fp, fsb, fs, phi, tau, sigma, gspread, amp_lin, amp_nl, amp_total)     
      
!Input arguments:

!          per, m, rjb, 
!          mech, v30, 
!          e_gmpe, amh_gmpe, c_gmpe, amref_gmpe,
!          rref_gmpe, h_gmpe, 
!          v30ref, 
!          sigt_gmpe,
!          e_pga4nl, amh_pga4nl, c_pga4nl, amref_pga4nl,
!          rref_pga4nl, h_pga4nl,


!Output arguments:

!          y, expsiglny,
!          pga4nl, fm_pga4nl, fd_pga4nl, 
!          fs_lin_pga4nl, fs_nonlin_pga4nl,
!          bnl, amp_lin, amp_nl, amp_total,
!          fm_gmpe, fd_gmpe, fs_gmpe,
!          gspread


! Computes NGAW2 NGA motions for given input variables

! Note: For generality, include a magnitude-dependent anelastic
! term, although the coefficients are 0.0 for the 2007 BA NGA GMPEs

 
! Dates: 02/27/13 - Written by D.M. Boore, based on ngaw2_gm_sub4y
!        04/09/13 - Revise in light of the new report and coefficient table.
!                 - Sediment depth term in terms of delta_z1, which requires the relation
!                   between z1 and Vs30.  There are two relations in the report
!                   (equations (4.9a) and (4.9b)), one for Japan and one for California.
!                   I need to add an input argument to specify which relation to use, but for
!                   now, just skip the sediment term altogether (even if it is included, specifying
!                   Z1<0 will turn it off).
!                 - Return a large value of phi if Rjb> 300
!        04/11/13 - Ignore restriction on Rjb for phi computation
!        07/01/13 - Renamed from bssa13_gm_sub4y
!        01/09/14 - Adjust mu1_vs30 output from units of m to km.

      implicit none

      
      real,intent(in) ::   per, m, rjb, r, v30, z1,  c(1:3), amref, rref, delta_c3(0:3), e(0:6), amh
      real  pga4nl
      real  clin, vclin, vref, f1, f2, f3, f4, f5, f6, f7, r1, r2, delta_phiR, delta_phiV, v1, v2
      real  phi1, phi2, tau1, tau2,  fe, fpb, fp, fsb, fs, phiM, phiMR
      real,intent(out) ::  phi, tau, sigma, gspread,y
      real,intent(out) ::  amp_lin, amp_nl, amp_total
      real delta_z1, f_delta_z1,mu1_vs30
      real y_xamp,bnl
     
      integer,intent(in) :: mech, iregion, irelation

        
 !  print *,per, m,rjb,r,v30
!   print *,pga4nl

!GET Y FOR GMPE, WITHOUT SITE AMPS:

      call y_bssa14_no_site_amps(m, r, mech, iregion,c, amref, rref, delta_c3,e, amh,y_xamp, fe, fpb, fp, gspread)     
      
!     print *,pga4nl
!Compute site amps

      if (v30 <= vclin) then
        amp_lin  =  (v30/vref)**clin
      else
        amp_lin  =  (vclin/vref)**clin
      end if
        
      f2 = f4 * (exp(f5*(amin1(v30,760.0)-360.0)) - exp(f5*(760.0-360.0)           )   )
      
      bnl = f2
      
      amp_nl = exp(f1+f2*alog( (pga4nl+f3)/f3 ))

      amp_total = amp_lin * amp_nl
      
      fsb  = alog(amp_total)   ! Natural log!!
      
      if (z1 < 0.0) then
        delta_z1 = 0.0
      else
        if (irelation == 1 .or. irelation == 2) then
          delta_z1 = z1 - mu1_vs30(v30, irelation)/1000.0
        else
          delta_z1 = 0.0
        end if
      end if
      
      if (per < 0.65) then
        f_delta_z1 = 0.0
      else
        if (delta_z1 > f7/f6) then
          f_delta_z1 = f7
        else
          f_delta_z1 = f6*delta_z1
        end if
      end if
      
      fs = fsb + f_delta_z1
      
      y = exp(fs) * y_xamp
       
!Compute phi, tau, and sigma   

      if (m <= 4.5) then
        tau = tau1
      else if (m >= 5.5) then
        tau = tau2
      else 
        tau = tau1+(tau2-tau1)*(m-4.5)
      end if
      
      if (m <= 4.5) then
        phiM = phi1
      else if (m >= 5.5) then
        phiM = phi2
      else 
        phiM = phi1+(phi2-phi1)*(m-4.5)
      end if
      
      if (Rjb <= R1) then
        phiMR = phiM
      else if (Rjb > R2) then
        phiMR=phiM+delta_phiR
      else
        phiMR = phiM+delta_phiR*(alog(Rjb/R1)/alog(R2/R1))
      end if

      if (v30 <= v1) then
        phi = phiMR - delta_phiV
      else if (v30 >= v2) then
        phi = phiMR
      else
        phi = phiMR - delta_phiV*(alog(v2/v30)/alog(v2/v1))
      end if
      
      sigma = sqrt(phi**2 + tau**2)
         
      return
      end subroutine bssa14_gm_sub4y
! --------------------------------------------------------- bssa14_gm_sub4y

! --------------------------------------------------------- mu1_vs30
      function mu1_vs30(vs30, irelation)     

! irelation   relation
!       1     California
!       2     Japan

! NOTE: the units of mu1_vs30 are m. not km

! Dates: 05/10/13 - Written by D.M. Boore 

      implicit none
      
      real :: ln_mu1, mu1_vs30, vs30 
     
      integer :: irelation
      
      if (irelation == 1) then 
      
        ln_mu1 = (-7.15/4.0)*alog((vs30**4+570.94**4)/(1360.0**4+570.94**4))
        mu1_vs30 = exp(ln_mu1)
        
      else if (irelation == 2) then
      
        ln_mu1 = (-5.23/2.0)*alog((vs30**2+412.39**2)/(1360.0**2+412.39**2))
        mu1_vs30 = exp(ln_mu1)      
      
      else
      
        ! No relation, return negative value
        mu1_vs30 = -9.9
      
      end if
      

      return
      
      end function mu1_vs30


subroutine get_rotd50(a1,a2,nt1,nt2,dt,per1,nper1,damp,rotd50)
!input: a1(nt1),a2(nt2)
!       per1(nper1)
!       damp
use mod_pgamisf, only : select
implicit none
real :: a1(1:nt1),a2(1:nt2),per1(nper1),damp,rotd50(nper1)
!local:
real, allocatable :: arot(:,:)
real, allocatable :: ts_rd(:),ts_rv(:)
integer :: nrot,n50,nt,irot, iper
real :: deg2rad,angrot,omega,d0,v0,dt,rd,rv,aa
real, allocatable :: psarot(:,:),wrk(:)
integer :: nper1,nt1,nt2
real*8,parameter:: PI=3.1415926535
!check the size
!nt1=size(a1);nt2=size(a2)
!print *,nt1,nt2


if (nt1.ne.nt2) then
 print *,'different length in horizontal components, taking the smaller one'
 if (nt1<nt2 .and. nt1>0) then
  nt=nt1
 elseif (nt2>0) then 
  nt=nt2
 else
  print *,'problem with horizontal compononent lengths: ',nt1,nt2
  stop
 endif
else
 nt=nt1
endif

nrot=180
n50=nrot/2
deg2rad=pi/180.
allocate(arot(nt,nrot))

allocate(ts_rd(nt),ts_rv(nt),wrk(nrot))
allocate(psarot(nrot,nper1))

do irot=1,nrot
 angrot=(irot-1.)*deg2rad
 arot(:,irot)=a1(:)*cos(angrot)+a2(:)*sin(angrot)
 !get response spectra of angrot
  do iper=1,nper1
   omega=2.*pi/per1(iper)
   d0=0.;v0=0.
!  print *,damp,minval(arot(:,irot))
   call rdrvaa_rd_rv_ts(arot(:,irot),nt,omega,damp,dt,rd,rv,aa,d0,v0,ts_rd(:),ts_rv(:))
   psarot(irot,iper)=omega**2*rd 
!   print *,psarot(irot,iper),rd
  enddo
enddo


!sort and find n50
do iper=1,nper1
!get rotd50
 wrk=psarot(:,iper)
! print *,maxval(wrk),maxval(psarot(:,iper))
 rotd50(iper)=select(n50,nrot,wrk)
enddo

deallocate(arot,ts_rd,ts_rv,wrk,psarot)

end subroutine


     subroutine rdrvaa_rd_rv_ts(acc,na,omega,damp,dt,rd,rv,aa,d0,v0,ts_rd, ts_rv)
! This is a modified version of "Quake.For", originally
! written by J.M. Roesset in 1971 and modified by
! Stavros A. Anagnostopoulos, Oct. 1986.  The formulation is that of
! Nigam and Jennings (BSSA, v. 59, 909-922, 1969).  This modification
! returns the time series of the relative displacement, in addition to
! rd, rv, and aa

!   acc = acceleration time series
!    na = length of time series
! omega = 2*pi/per
!  damp = fractional damping (e.g., 0.05)
!    dt = time spacing of input
!    rd = relative displacement of oscillator
!    rv = relative velocity of oscillator
!    aa = absolute acceleration of oscillator
! d0,v0 = initial displacement and velocity (usually set to 0.0)

! Dates: 05/06/95 - Modified by David M. Boore
!        04/15/96 - Changed name to RD_CALC and added comment lines
!                   indicating changed for storing the oscillator time series
!                   and computing the relative velocity and absolute
!                   acceleration
!        04/16/96 - This is RD_CALC, with the time series of the relative
!                   displacement added and the name changed
!        03/14/01 - Made this double precision, but did not rename it as I
!                   did rc_calc and rdrvaa.  Also added initial displacement 
!                   and velocity.
!        01/31/03 - Moved implicit statement before the type declarations
!        10/10/07 - Initial variable assignments and iteration loop modified 
!                   to double-precision (Chris Stephens)
!        02/08/08 - Set ts(1)=0.0 and fill ts array up to na points.  The previous
!                   version of the program filled ts to na-1, but the program 
!                   calling this subroutine thought that the time series had
!                   a length of na.  What I do now is shift the ts by dt, assuming
!                   that ts(1) = 0.0
!        08/14/12 - Included 12/22/10 revision made to rscalc_ts: 
!                     Remove minus sign in front of the initialization of y, ydot. 
!                     The minus sign was a remnant of an earlier version where I did not
!                     understand the meaning of y and ydot.
!                 - Included 12/26/10 revision made to rscalc_ts: 
!                     Correct error that first vel and dis points = 0.0
!        09/22/12 - Renamed from rdcalcts, and now include rv and aa in output
!                   (following rdrvaa.for).
!        02/07/18 - Modified from rdrvaa_rd_ts.for, and nows writes out ts_rv, following rscalc_ts. 
!                   I am trying to decrease the number of programs that compute oscillator 
!                   response time series, as things have gotten quite confusing.  The
!                   program smc2psa_rot_gmrot_interp_acc_rot_osc_ts.for had used the
!                   subroutine rdrvaa_rd_ts, perhaps to avoid the extra time in computing
!                   ts_rv, but that extra cost is minor compared to the confusion of having
!                   too many programs that do almost the same thing.  So now 
!                   smc2psa_rot_gmrot_interp_acc_rot_osc_ts.for and smc2rs_ts.for will both
!                   call the subroutine rdrvaa_rd_rv_ts.  The subroutines rscalc_ts.for, rdcalcts.for, 
!                   and rdrvaa_rd_ts.for have all been saved with their latest dates as an extension,
!                   and the files with the actual names have been modified by deleting all
!                   content and adding the line
!                   ! Use rdrvaa_rd_rv_ts.for instead

!                   

      IMPLICIT NONE
      REAL acc(*), omega, damp, dt, rd, rv, aa, d0, v0, ts_rd(*), ts_rv(*)
!locals:
      real :: omt, d2,bom,d3,omd,om2,omdt,c1,c2,c3,c4,ss,cc,bomt,ee,s1,s2,s3,a11,a12,a21,a22,s4,s5
      real :: b11,b12,b21,b22,y,ydot,y1,z,z1,z2,ra
      integer :: i,na


      omt=dble(omega)*dble(dt)
      d2=1-dble(damp)*dble(damp)
      d2=sqrt(d2)
      bom=dble(damp)*dble(omega)
      d3 = 2.*bom                 ! for aa
      omd=dble(omega)*d2
      om2=dble(omega)*dble(omega)
      omdt=omd*dble(dt)
      c1=1.d0/om2
      c2=2.d0*dble(damp)/(om2*omt)
      c3=c1+c2
      c4=1.d0/(dble(omega)*omt)
      ss=sin(omdt)
      cc=cos(omdt)
      bomt=dble(damp)*omt
      ee=exp(-bomt)
      ss=ss*ee
      cc=cc*ee
      s1=ss/omd
      s2=s1*bom
      s3=s2+cc
      a11=s3
      a12=s1
      a21=-om2*s1
      a22=cc-s2
      s4=c4*(1.d0-s3)
      s5=s1*c4+c2
      b11=s3*c3-s5
      b12=-c2*s3+s5-c1
      b21=-s1+s4
      b22=-s4
      rd=0.
      rv = 0.                           ! for rv
      aa = 0.                           ! for aa

      y = dble(d0)
      ydot = dble(v0) ! These are initial values of the oscillator
!      y=0.
!      ydot=0.

      ts_rd(1) = d0  ! 14aug12
      ts_rv(1) = v0  ! 14aug12
!      ts(1) = 0.0  ! 08feb08

      DO i=1, na-1   ! 08feb08 (used new variable n1=na-1, but I replaced the upper limit with 
                     ! na-1 on 09/22/12, as n1 was only used here)

        y1=a11*y+a12*ydot+b11*dble(acc(i))+b12*dble(acc(i+1))
        ydot=a21*y+a22*ydot+b21*dble(acc(i))+b22*dble(acc(i+1))
        y=y1    ! y is the oscillator output at time corresponding to index i+1 
                ! changed "i" in this comment to "i+1" on 14aug12
        ts_rd(i+1) = y   !08feb08
        ts_rv(i+1) = ydot    
        z=abs(y)
        if (z > rd) rd=z
        z1 = abs(ydot)                   ! for rv
        if (z1 > rv) rv = z1            ! for rv
        ra = -d3*ydot -om2*y1            ! for aa
        z2 = abs(ra)                     ! for aa
        if (z2 > aa) aa = z2            ! for aa

      END DO

      RETURN
      END SUBROUTINE
!nFreq=nper
!!     Convert periods to freq in Radians
!    do iFreq=1,nFreq
!       w(iFreq) = 2.0*3.14159 / per(iFreq)
!    enddo
!!       Check that the two time series have the same number of points.  If not, reset to smaller value
!        if (npts1 .lt. npts2) then
!          npts0 = npts1
!        elseif (npts2 .lt. npts1) then
!          npts0 = npts2
!        elseif (npts1 .eq. npts2) then
!          npts0 = npts1
!        endif
!!       Check that the two time series have the same dt.
!        if (dt1 .ne. dt2) then
!          write (*,*) 'DT values are not equal!!!'
!!          write (*,*) 'DT1 = ', dt1
!!          write (*,*) 'DT2 = ', dt2
!        else
!          dt = dt1
!        endif
!        dt0 = dt 
!!        Copy to new array for interpolating ( from original acc )
!         do i=1,npts0
!           acc1(i) = acc1_0(i)
!           acc2(i) = acc2_0(i)
!         enddo    
!         npts = npts0
!         dt = dt0
!!        Loop over each oscilator frequency
!         do iFreq=1,nFreq 
!         Compute the oscillator time histoires for the two components. 
!          call CalcRspTH ( acc1, npts, dt, w(iFreq), damping, rspTH1 )
!          call CalcRspTH ( acc2, npts, dt, w(iFreq), damping, rspTH2 )
!!         Fill new array with points with amplitude on one component at least SaMin/1.5
!!         This sets the points for the rotation to speed up the calculation
!          call Calc_Sa ( rspTH1, sa1, npts)
!          call Calc_Sa ( rspTH2, sa2, npts )
!          test = amin1(sa1, sa2) / 1.5
!          j = 1
!          do i=1,npts
!            amp1 = abs(rspTH1(i))
!            amp2 = abs(rspTH2(i))
!            if ( amp2 .gt. amp1 ) amp1 = amp2
!            if ( amp1 .gt. test .and. iFlag .eq. 0 ) then
!              rsp1(j) = rspTH1(i)
!              rsp2(j) = rspTH2(i)
!              j = j + 1
!            endif 
!          enddo
!          npts1 = j -1
!!         Loop over different rotation angles and compute response spectra by rotating the Oscillator TH
!          do j=1,90
!            rotangle = real(((j-1)*3.14159)/180.0)
!            cos1 = cos(rotangle)
!            sin1 = sin(rotangle)
!            do i=1,npts1
!              x(i)=cos1*rsp1(i) - sin1*rsp2(i)
!              y(i)=sin1*rsp1(i) + cos1*rsp2(i)
!            enddo
!!          Find the maximum response for X and Y and load into a single Sa array
!            call Calc_Sa ( x, saX, npts1 )
!            call Calc_Sa ( y, saY, npts1 )
!            sa(j) = saX
!            sa(j+90) = SaY
!          enddo      
!         Get the as-recorded PSa
!          psa5E(iFreq) = sa(1)
!          psa5N(iFreq) = sa(91)

!         Sort the Sa array to find the median value.
!          n1 = 180
!          call SORT(Sa,WorkArray,N1)
!          rotD50(iFreq) = ( Sa(90) + Sa(91) ) /2.
!         enddo
!!        Find the Famp1.5 (assumes order of freq are high to low)
!         do iFreq=2,nFreq 
!           shape1 = rotD50(iFreq)/rotD50(1)
!           if ( shape1 .ge. 1.5 ) then
!             famp15 = 1./rsp_period(iFreq)
!!             goto 105
!           endif
!         enddo
!!  105    continue



