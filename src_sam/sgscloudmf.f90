
subroutine sgscloudmf(lenv)


!output:
! cfrac_tot ... sum of stratiform cover in non-convective region
!               and convective cloud cover
! qn        ... domain mean non-precip cloud condensate
! qv        ... domain mean watervapor
! qcl       ... domain mean liquid condensate
! qci       ... domain mean frozen condensate
! tabs      ... domain mean absolute temperature
! theta     ... domain mean potential temperature
! thetav    ... domain mean virtual potential temperature
! thetar    ... domain mean density potential temperature (no precip considered since 
!               precip is grid-scale thus no effect on plume buoyancy)
! c0,c1     ... coefficients needed for buoyancy flux

!input:   ... mean total water and mean moist static energy
!         ... plume total water and mse
!         ... grid-scale precip
!         ... plume condensate (ice and liquid)
!         ... plume area fraction


! strategy:
! iterate over the following steps
! 1) compute environmental mse and qt from the domain-mean and in-cloud plume values
! 2) get first guess for tabs in environment assuming no condensate in environment (only in condensed plumes)

! 3) get environmental thetali from temperature and condensate (which is zero in this first step only)
! 4) get vertical gradients of thetali and qt by using difference to level below
! 5) get variance of saturation deficit, sigmas, needed for PDF scheme
!    following Bechtold 1992,1995 and Sommerica and Deardorff (1977), extended here to include ice
! 6) get environmental cloud condensate, split into liquid and ice based on temperature
! 7) get new environmental temperature based on environemtnal mse and environmental condensate


! iterate steps 3-7

! 8) get domain-mean qcl and qci as area-weighted sum of plume and environmental qc and qi
! 9) get domain-mean water vapor content qv=qt-qn
! 10) get domain mean tabs, theta, thetav, thetar

use vars
use params
use grid
use micro_params


implicit none

logical, intent(in) :: lenv

!local variables
integer k, kb, kc, n,ks
real dtabs, an, bn, ap, bp, om, ag, omp, omn
real, dimension(nzm) :: qte
real fac1,fac2, fac_a
real fff,dfff,dqsat
real lstarn,dlstarn,lstarp,dlstarp, dtabsa, hlip, hlie, frac_mf2, qce, qie, qne, tabse, ds
integer,parameter :: niter=100
real :: lambdaf, alphaf, qsl, totheta, tl, qsw, qsi, cfrac_l, cfrac_t, tabsnoql, dqndt
real,dimension (2) :: tabs1, tabs2

real, parameter :: zsig_max = 1.0e-2

real, parameter :: q_crit=1.6

real :: cab,ckk,leps
real,dimension(nzm) :: qtqt, thlthl, qtthl
real, dimension(2:nzm-1) :: astar,bstar,cstar,dstar
real, dimension(nzm-2) :: atri,btri,ctri,dtri


if (dosgscloud) then

an = 1./(tbgmax-tbgmin) 
bn = tbgmin * an
ap = 1./(tprmax-tprmin) 
bp = tprmin * ap
fac1 = fac_cond+(1+bp)*fac_fus
fac2 = fac_fus*ap
ag = 1./(tgrmax-tgrmin) 

cfrac_t=0.5 * (1. + erf( 3.  ))
cfrac_l=0.5 * (1. + erf( -3. ))

hlip=0.

do k=1,nzm

kb=k-1
kc=k+1
if(k.eq.1) then
 kb=1
 kc=2
elseif (k.eq.nzm) then
 kb=nzm-1
 kc=nzm
end if

totheta=(pres(k)/p00)**(rgas/cp)

! environment properties
! interpolate plume area fraction to mass level
 frac_mf2 = 0.5*(frac_mf(k)+frac_mf(k+1))
if (frac_mf2.gt.0.0.and.lenv)  then
  if (frac_mf(k+1).gt.0.0) then
    ! get plume moist static energy
    hlip = 0.5*(cp*tabs_mf(k) + ggr*zi(k) - lcond * qcsgs_mf(k) - lsub * qisgs_mf(k) +&
            cp*tabs_mf(k+1) + ggr*zi(k+1) - lcond * qcsgs_mf(k+1) - lsub * qisgs_mf(k+1))
    ! total water in environment
    qte(k) = (qt(k) - frac_mf2*0.5*(qtsgs_mf(k)+qtsgs_mf(k+1)))/(1.-frac_mf2)
  else
    ! get plume moist static energy
    hlip = cp*tabs_mf(k) + ggr*zi(k) - lcond * qcsgs_mf(k) - lsub * qisgs_mf(k)
    ! total water in environment
    qte(k) = (qt(k) - frac_mf2*qtsgs_mf(k))/(1.-frac_mf2)
  end if
  hlie   = ((t(k)+lcond*qpl(k)+lsub*qpi(k)-frac_mf2*hlip))/(1.-frac_mf2)
else
  qte(k) = qt(k)
  hlie   = t(k)+lcond*qpl(k)+lsub*qpi(k)
end if


! thetali in environment
 thetali(k) = (hlie-ggr*z(k))  / cp / totheta 

end do

do k=1,nzm
! compute vertical gradients

kb=k-1
kc=k+1
if(k.eq.1) then
 kb=1
 kc=2
elseif (k.eq.nzm) then
 kb=nzm-1
 kc=nzm
end if

 thetaligrad(k) = (thetali(kc)-thetali(kb))/(z(kc)-z(kb))
 qtgrad(k)      = (qte(kc)-qte(kb))/(z(kc)-z(kb))
 tke(k)         = max(0.0,tke(k))

end do

if (dovartrans) then

  ! include transport of (co-)variance and numerically solve 2nd
  ! order ODE
  ! set lower BC: balance of diss and prod
  qtqt(1)  = lcld**2*qtgrad(1)**2.
  thlthl(1)= lcld**2*thetaligrad(1)**2.
  qtthl(1) = lcld**2*qtgrad(k)*thetaligrad(1) 
  ! set upper BC: zero variance
  qtqt(nzm)  = 0.0
  thlthl(nzm)= 0.0
  qtthl(nzm) = 0.0


  ! compute coefficients to solve for thlthl
  cab = 2.5   ! same as in Bechtold 1992, JAS
  ckk = 0.5   ! has to be same as in TKE scheme
  do k=2,nzm-1
    ks=k-1
    leps = (smix(k)+1.d-10)/2.5  ! Same as in TKE scheme
    astar(k) = cab / ckk / leps / (smix(k)+1.0d-10)
    dstar(k) = 2. * thetaligrad(k) * thetaligrad(k)
    if (tke(k).eq.0.0) then
      cstar(k)=0.
      bstar(k)=0.
    else  
      cstar(k)=1.
      bstar(k) = (log(smix(k+1)*sqrt(tke(k+1))+1.d-12) - log(smix(k-1)*sqrt(tke(k-1))+1.d-12))/(adz(k) + adzw(k+1))/dz
    end if
    
    if (k.gt.2.and.k.lt.nzm-1) then
      atri(ks) = cstar(k)/adz(k)/adzw(k)/dz/dz  - bstar(k)/(adz(k) + adzw(k+1) )/dz
      btri(ks) = astar(k) - (cstar(k)/adzw(k)/dz+1./adzw(k+1)/dz)/adz(k)/dz
      ctri(ks) = cstar(k)/adz(k)/adzw(k+1)/dz/dz  + bstar(k)/(adz(k) + adzw(k+1) )/dz
      dtri(ks) = dstar(k)
    elseif (k.eq.2) then
      atri(ks) = 0.0
      btri(ks) = 1. - (cstar(k)/adzw(k)/dz+1./adzw(k+1)/dz)/adz(k)/dz
      ctri(ks) = cstar(k)/adz(k)/adzw(k+1)/dz/dz  + bstar(k)/(adz(k) + adzw(k+1) )/dz
      dtri(ks) = dstar(k) - thlthl(1) * (1./adz(k)/adzw(k)/dz/dz  - bstar(k)/(adz(k) + adzw(k+1) )/dz)
    else
      atri(ks) = cstar(k)/adz(k)/adzw(k)/dz/dz  - bstar(k)/(adz(k) + adzw(k+1) )/dz
      btri(ks) = 1. - (cstar(k)/adzw(k)/dz+1./adzw(k+1)/dz)/adz(k)/dz
      ctri(ks) = 0.0
      dtri(ks) = dstar(k)
    end if 
  end do
  !solve tridiagonal matrix
  call tridiag(atri,btri,ctri,dtri,nzm-2)
  thlthl(2:nzm-1) = dtri 
  varwrt1=0.
  varwrt1(2:nzm-1) = thlthl(2:nzm-1)

  qtqt = 0.
  qtthl = 0.
  
   
else

   ! use balance between dissipation and production
   ! to diagnose (co-)variances, then apply limiters

   ! assuming leps=l/2.5 the variances are 2 l^2 dadz dbdz
   ! use constant length scale of 250 m from Cheineit, which equals l=250 * sqrt(2)  in the variance formulation
   do k=1,nzm
    qtqt(k)  = lcld**2*qtgrad(k)**2. 
    thlthl(k)= lcld**2*thetaligrad(k)**2.
    qtthl(k) = lcld**2*qtgrad(k)*thetaligrad(k)

    ! limit variances
    ! Kay's limiters
    !thlthl(k)=max(min(thlthl(k),10.),0.01)
    !qtqt(k)=max(min(qtqt(k),1.e-5),1.e-8)
    !qtthl(k)=sign(max(min(abs(qtthl(k)),1.e-2),1.e-5),qtthl(k))
    ! max/min taken from Heinze 2005, JAMES
    thlthl(k)=max(min(thlthl(k),3.),0.001)
    qtqt(k)=max(min(qtqt(k),2.e-6),1.e-8)
    qtthl(k)=max(min(qtthl(k),0.02d-3),-4.d-3)
   end do
end if

do k=1,nzm
 ! get liquid water temperature Tl
 totheta=(pres(k)/p00)**(rgas/cp)
 tl     = totheta * thetali(k)
 tabsnoql = tl 
 tabse=tabsnoql

 ! get saturation mixing ratio at Tl
 omn = max(0.,min(1.,(tl-tbgmin)*an))
 qsw=qsatw(tl,pres(k))
 qsi=qsati(tl,pres(k))
 qsl = qsw ! omn * qsw + (1.-omn) * qsi
 ! get dqsdT
 !alphaf = (omn * qsw * lcond/(tl**2*rv) + (1.-omn) * qsi * lsub/(tl**2*rv))
 alphaf = qsw * lcond/(tl**2*rv) 
 !lambdaf = (1. + alphaf * (omn*fac_cond + (1.-omn) * fac_sub) )**(-1.)
 lambdaf = (1. + alphaf *fac_cond )**(-1.)
 alphaf =alphaf * totheta


 ds=lambdaf *(qte(k)-qsl)
 sigmas(k) = lambdaf * (max(0.0,qtqt(k) + alphaf**2 * thlthl(k)**2. - 2. * alphaf * qtthl(k)))**(0.5)
 !sigmas(k) = MIN ( zsig_max, sigmas(k) )

 n=0
 dtabs=100.
 do while (abs(dtabs).gt.0.01.and.n.le.niter) 

 if (n.gt.0) then ! not first guess anymore
   tabse  = tabse+dtabs
   ! compute mean saturation deficit based on mean temperature
   ! get saturation mixing ratio at T
   omn = max(0.,min(1.,(tabse-tbgmin)*an))
   qsw=qsatw(tabse,pres(k))
   qsi=qsati(tabse,pres(k))
   qsl = qsw !omn * qsw + (1.-omn) * qsi
   ds=qte(k)-qsl
 end if

 IF (sigmas(k).le.0.0) then
   cfrac_pdf(k) =   ABS ( (SIGN(1.0,ds)+1.0)*0.5 )
   qne = cfrac_pdf(k) * ds
   q1(k)=-999.
   dqndt = -dtqsatw(tabse,pres(k))
 ELSE
   q1(k)= ds/sigmas(k)
   cfrac_pdf(k) = MIN ( 1.0, MAX ( 0.0, &
                                        0.5 * (1.0+q1(k)/q_crit) ) )
   IF ( q1(k) .le. - q_crit ) THEN
      qne = 0.0
      dqndt = 0.
   ELSEIF ( q1(k) .ge. q_crit ) THEN
      qne = sigmas(k) * q1(k) 
      dqndt = -dtqsatw(tabse,pres(k))
   ELSE
      qne = sigmas(k) * (q1(k)+q_crit) * (q1(k)+q_crit) / (2.*(q_crit+q_crit))
      dqndt = - (q1(k)+q_crit)/(2.*q_crit) * dtqsatw(tabse,pres(k))
   ENDIF
 END IF
 
 !if (n.gt.0) then
   !dtabs= (tl + an*qne*tbgmin*fac_fus + qne*fac_sub  )   &
   !       / (1. + an*qne*fac_fus ) - tabse
   dtabs= (tabsnoql - tabse + fac_cond * qne) / (1.- fac_cond * dqndt)
 !else
   !tabse  = (tl + an*qne*tbgmin*fac_fus + qne*fac_sub  )   &
   !       / (1. + an*qne*fac_fus )
 !end if

 n=n+1

 end do

 ! partition condensate into liquid and ice using temperature
 omn = max(0.,min(1.,(tabse-tbgmin)*an))
 qce = qne*omn
 qie = qne*(1.-omn)


 cthl(k) = cfrac_pdf(k) * (1. - (fac_cond/totheta - (1.+epsv) * tabse/totheta)*lambdaf * alphaf) & ! see bechtold 95 (a,b,alpha,beta) and my own notes
 !fac_a = (1.-qte + (1.+epsv)*qsl)*(1.+epsv*lcond/rgas/tabse)  & 
 !          / (1.+epsv*lcond/rgas/tl**2 *fac_cond *qsl)
 !cthl(k) = cfrac_pdf(k) * fac_a & ! see bechtold 95 (a,b,alpha,beta) and my own notes
      + (1.-cfrac_pdf(k)) * (1.+epsv*qte(k))
 cqt(k) =  cfrac_pdf(k) * (epsv*tabse/totheta + (fac_cond/totheta + (1.+epsv) * tabse/totheta) * lambdaf) &
 !cqt(k) =  cfrac_pdf(k) * (fac_cond/tabse*fac_a - 1.)*tabse/totheta &
      + (1.-cfrac_pdf(k)) * epsv*thetali(k)

 
 ! get final domain averages (convective and environment)
 qcl(k) = (1.-frac_mf2) * qce + 0.5 * frac_mf2 * (qcsgs_mf(k)+qcsgs_mf(k+1))
 qci(k) = (1.-frac_mf2) * qie + 0.5 * frac_mf2 * (qisgs_mf(k)+qisgs_mf(k+1))
 cfrac_tot(k) = min(frac_mf2,0.5*(cfrac_mf(k+1)+cfrac_mf(k))) + (1.-frac_mf2) * cfrac_pdf(k)
 qn(k) = max(qcl(k) + qci(k),0.)
 qv(k) = max(0.,qt(k) - qn(k))
 
 tabs(k)  = (t(k)-ggr*z(k))/cp + fac_cond*(qpl(k)+qcl(k)) + fac_sub*(qpi(k)+qci(k)) 
 theta(k)  = tabs(k)/totheta  
 thetav(k) = theta(k)  * (1.+epsv*qv(k))
 thetar(k) = theta(k)  * (1.+epsv*qv(k)-qn(k))
 thetali(k)= theta(k) - 1./totheta *(fac_cond*qcl(k) + fac_sub*qci(k)) 

end do

else  ! use all-or-nothing scheme

  do k=1,nzm
    call cloud(k)
  end do ! loop over k

end if

lwp=0.
do k=1,nzm
  lwp = lwp + rho(k) * qcl(k) * adz(k)*dz
end do
lwp = lwp * 1000. ! convert to g/m2


end subroutine sgscloudmf

