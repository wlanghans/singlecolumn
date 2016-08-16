
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
integer k, kb, kc
real dtabs, an, bn, ap, bp, om, ag, omp, omn
real fac1,fac2
real fff,dfff,dqsat
real lstarn,dlstarn,lstarp,dlstarp, dtabsa, hlip, hlie, frac_mf2, qce, qie, qte, qtel, qne, tabse
integer niter
real :: lambdaf, alphaf, qsl, totheta, tl, qsw, qsi, cfrac_l, cfrac_t
real,dimension (2) :: tabs1, tabs2

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


if (dosgscloud) then


if (k.eq.1.or.tke(k).le.0.0) then

! no PDF scheme needed if no variability (or if in first layer)
! if condensate is present, then it comes from all-or-nothing adjustment
! and will not be modified here
 tabs(k)=(t(k)-ggr*z(k))/cp + fac_cond*(qcl(k)+qpl(k)) + fac_sub*(qci(k)+qpi(k)) 
 thetali(k) = (tabs(k)- fac_cond*(qcl(k)) - fac_sub*(qci(k)) ) / totheta 
 theta(k)   = tabs(k)/totheta
 thetav(k)  = theta(k) * (1.+epsv*qv(k)) 
 ! density potential temperature will later be used in edmf to get buoyancy. 
 ! since qp is thought to be spread out evenly over grid-box it will not contribute to B
 thetar(k)  = thetav(k) - tabs(k)/totheta * qn(k)
 cfrac_pdf(k) = 0.

 omn = max(0.,min(1.,(tabs(k)-tbgmin)*an))
 qsw=qsatw(tabs(k),pres(k))
 qsi=qsati(tabs(k),pres(k))
 qsl = omn * qsw + (1.-omn) * qsi

 cthl(k) = 1.+epsv*qt(k)
 cqt(k) =  epsv*thetali(k)
 q1(k)  = qt(k) - qsl
 qn(k)  = max(0.0,q1(k))
 qcl(k) = omn*qn(k) 
 qci(k) = (1.-omn)*qn(k)
 qv(k)  = qt(k)-qn(k) 
 if (qn(k).gt.0.0) then
   cfrac_tot(k) = 1.
 else 
   cfrac_tot(k) = 0.0
 end if
 qte = qt(k)
 sigmas(k)  = 0.
 varwrt1(k) = 100.

else

! interpolate plume area fraction to mass level
 frac_mf2 = 0.5*(frac_mf(k)+frac_mf(k+1))
if (frac_mf2.gt.0.0.and.lenv)  then
  if (frac_mf(k+1).gt.0.0) then
    ! get plume moist static energy
    hlip = 0.5*(cp*tabs_mf(k) + ggr*zi(k) - lcond * qcsgs_mf(k) - lsub * qisgs_mf(k) +&
            cp*tabs_mf(k+1) + ggr*zi(k+1) - lcond * qcsgs_mf(k+1) - lsub * qisgs_mf(k+1))
    ! total water in environment
    qte = (qt(k) - frac_mf2*0.5*(qtsgs_mf(k)+qtsgs_mf(k+1)))/(1.-frac_mf2)
  else
    ! get plume moist static energy
    hlip = cp*tabs_mf(k) + ggr*zi(k) - lcond * qcsgs_mf(k) - lsub * qisgs_mf(k)
    ! total water in environment
    qte = (qt(k) - frac_mf2*qtsgs_mf(k))/(1.-frac_mf2)
  end if
  hlie   = ((t(k)+lcond*qpl(k)+lsub*qpi(k)-frac_mf2*hlip))/(1.-frac_mf2)
else
  qte = qt(k)
  hlie   = t(k)+lcond*qpl(k)+lsub*qpi(k)
end if

! first guess for tabs in environment outside of plumes: assume no condensate in environment
 tabse  = (hlie - ggr*z(k))/cp
! thetali in environment: assume no condensate present
 thetali(k) = tabse / totheta 
 
 niter=0
 dtabs = 100.
 do while(abs(dtabs).gt.0.01.and.niter.lt.20)

 if (niter.gt.0) then 
   tabse= tabse + 0.5 * dtabs
   thetali(k) = tabse / totheta - (lcond * qce + lsub*qie)/totheta
 end if

    ! environmental properties

    ! compute gradients of qt and thetali
    thetaligrad(k) = (thetali(k)-thetali(k-1))/(z(k)-z(k-1))
    qtgrad(k)      = (qte-qtel)/(z(k)-z(k-1))
    if (k.eq.2) then
      thetaligrad(1) = thetaligrad(2)
      qtgrad(1)      = qtgrad(2)
    end if

    ! get liquid water temperature Tl
    tl     = totheta * thetali(k)
    ! get saturation mixing ratio at Tl
    omn = max(0.,min(1.,(tl-tbgmin)*an))
    qsw=qsatw(tl,pres(k))
    qsi=qsati(tl,pres(k))
    qsl = omn * qsw + (1.-omn) * qsi
    ! get dqsdT
    alphaf = (omn * qsw * lcond/(tl**2*rv) + (1.-omn) * qsi * lsub/(tl**2*rv))
    lambdaf = (1. + alphaf * (omn*fac_cond + (1.-omn) * fac_sub) )**(-1.)
    alphaf =alphaf * totheta
    ! use constant length scale of 250 m
    sigmas(k) = lcld*(max(0.0,qtgrad(k)**2. + alphaf**2*thetaligrad(k)**2. - 2. * alphaf * thetaligrad(k)*qtgrad(k)))**(0.5)
    sigmas(k) = max(sigmas(k),1.d-12)
    varwrt1(k) = 0.

    ! compute saturation deficit
    q1(k)= qte-qsl

    ! get cloud fraction and mean liquid water content
    cfrac_pdf(k) = 0.5 * (1. + erf(   max(-3.,min(3.,q1(k)/(sigmas(k)*1.41)))    ))                           ! Bechtold Eq (20)
    if (cfrac_pdf(k).eq.cfrac_l) cfrac_pdf(k) = 0.
    if (cfrac_pdf(k).eq.cfrac_t) cfrac_pdf(k) = 1.
    if (cfrac_pdf(k).eq.0.0.or.cfrac_pdf(k).eq.1.) sigmas(k) = 0.0
    cthl(k) = cfrac_pdf(k) * (1. - (fac_cond/totheta - (1.+epsv) * tabse/totheta)*lambdaf * alphaf) & ! see bechtold 95 (a,b,alpha,beta) and my own notes
         + (1.-cfrac_pdf(k)) * (1.+epsv*qte)
    cqt(k) =  cfrac_pdf(k) * (epsv*tabse/totheta + (fac_cond/totheta + (1.+epsv) * tabse/totheta) * lambdaf) &
         + (1.-cfrac_pdf(k)) * epsv*thetali(k)

    if (sigmas(k).ne.0.0) then
      qne = max(0.0,lambdaf  * (cfrac_pdf(k) * q1(k) + sigmas(k)*exp(-((q1(k)/sigmas(k))**2)/2.)/2.51 ))  ! Bechtold Eq (21)
    else
      qne = max(0.0,lambdaf  * (cfrac_pdf(k) * q1(k)))
    end if
    if (qne.lt.qp_threshold) qne=0.
 
    ! re-define  sigmas and q1 to be in line with definition in Bechtold
    sigmas(k)  = sigmas (k) * alphaf/2.
    q1(k)      = q1(k) * alphaf/2.
 
    ! partition condensate into liquid and ice using temperature
    omn = max(0.,min(1.,(tabse-tbgmin)*an))
    qce = qne*omn
    qie = qne*(1.-omn)

    ! get new estimates
    niter = niter+1
    dtabs = (hlie - ggr*z(k) + lcond * qce + lsub*qie )/cp - tabse

    end do ! end of iteration

    ! get final domain averages
    qcl(k) = (1.-frac_mf2) * qce + 0.5 * frac_mf2 * (qcsgs_mf(k)+qcsgs_mf(k+1))
    qci(k) = (1.-frac_mf2) * qie + 0.5 * frac_mf2 * (qisgs_mf(k)+qisgs_mf(k+1))
    cfrac_tot(k) = min(frac_mf2,0.5*(cfrac_mf(k+1)+cfrac_mf(k))) + (1.-frac_mf2) * cfrac_pdf(k)
    qn(k) = max(qcl(k) + qci(k),0.)
    qv(k) = max(0.,qt(k) - qn(k))
    
    tabs(k)  = (t(k)-ggr*z(k))/cp + fac_cond*(qpl(k)+qcl(k)) + fac_sub*(qpi(k)+qci(k)) 
    theta(k)  = tabs(k)/totheta  
    thetav(k) = theta(k)  * (1.+epsv*qv(k))
    thetar(k) = theta(k)  * (1.+epsv*qv(k)-qn(k))

end if ! if tke is present

! store qte for computation of gradient at next level higher up
qtel=qte 


else ! dont use sgs cloud scheme

tabs(k)=(t(k)-ggr*z(k))/cp + fac_cond*(qcl(k)+qpl(k)) + fac_sub*(qci(k)+qpi(k))
thetali(k) = (tabs(k)- fac_cond*(qcl(k)) - fac_sub*(qci(k)) ) / totheta
theta(k)   = tabs(k)/totheta
thetav(k)  = theta(k) * (1.+epsv*qv(k))
thetar(k)  = thetav(k) - tabs(k)/totheta * qn(k)
cfrac_pdf(k) = 0.

omn = max(0.,min(1.,(tabs(k)-tbgmin)*an))
qsw=qsatw(tabs(k),pres(k))
qsi=qsati(tabs(k),pres(k))
qsl = omn * qsw + (1.-omn) * qsi

cthl(k) = 1.+epsv*qt(k)
cqt(k) =  epsv*thetali(k)
q1(k)  = qt(k) - qsl
qn(k)  = max(0.0,q1(k))
qcl(k) = omn*qn(k)
qci(k) = (1.-omn)*qn(k)
qv(k)  = qt(k)-qn(k)
if (qn(k).gt.0.0) then
  cfrac_tot(k) = 1.
else
  cfrac_tot(k) = 0.0
end if
sigmas(k)  = 0.

 
end if

end do ! loop over k

end subroutine sgscloudmf

