
subroutine cloud(k)

!  Condensation of cloud water/cloud ice.

use vars
use micro_params
use params

implicit none

integer,intent(in) :: k

integer i,j,kb, kc
real dtabs, tabs1, an, bn, ap, bp, om, ag, omp, omn, alphaf, lambdaf
real fac1,fac2, totheta  
real fff,dfff,qsatt,dqsat
real lstarn,dlstarn,lstarp,dlstarp
integer niter

an = 1./(tbgmax-tbgmin)	
bn = tbgmin * an
ap = 1./(tprmax-tprmin)	
bp = tprmin * ap
fac1 = fac_cond+(1+bp)*fac_fus
fac2 = fac_fus*ap
ag = 1./(tgrmax-tgrmin)	

    totheta=(pres(k)/p00)**(rgas/cp)


    qt(k)=max(0.,qt(k))


! Initail guess for temperature assuming no cloud water/ice:


    tabs(k) = (t(k)-ggr*z(k))/cp


    !tabs1=(tabs(k)+fac1*qp(k))/(1.+fac2*qp(k))
    !WL the correct form should be 
    tabs1= tabs(k)*(1.- fac2*qp(k)) + fac1*qp(k)

! Warm cloud:

    if(tabs1.ge.tbgmax) then

      tabs1=tabs(k)+fac_cond*qp(k)
      qsatt = qsatw(tabs1,pres(k))

! Ice cloud:

    elseif(tabs1.le.tbgmin) then

      tabs1=tabs(k)+fac_sub*qp(k)
      qsatt = qsati(tabs1,pres(k))

! Mixed-phase cloud:

    else

      om = an*tabs1-bn
      qsatt = om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k))

    endif


!  Test if condensation is possible:


    if(qt(k).gt.qsatt) then

      niter=0
      dtabs = 100.
      do while(abs(dtabs).gt.0.01.and.niter.lt.10)
	if(tabs1.ge.tbgmax) then
	   om=1.
	   lstarn=fac_cond
	   dlstarn=0.
	   qsatt=qsatw(tabs1,pres(k))
	   dqsat=dtqsatw(tabs1,pres(k))
        else if(tabs1.le.tbgmin) then
	   om=0.
	   lstarn=fac_sub
	   dlstarn=0.
	   qsatt=qsati(tabs1,pres(k))
	   dqsat=dtqsati(tabs1,pres(k))
	else
	   om=an*tabs1-bn
	   lstarn=om*fac_cond+(1.-om)*fac_fus
	   dlstarn=an*fac_fus
	   qsatt=om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k))
	   dqsat=om*dtqsatw(tabs1,pres(k))+(1.-om)*dtqsati(tabs1,pres(k))
	endif
	if(tabs1.ge.tprmax) then
	   omp=1.
	   lstarp=fac_cond
	   dlstarp=0.
        else if(tabs1.le.tprmin) then
	   omp=0.
	   lstarp=fac_sub
	   dlstarp=0.
	else
	   omp=ap*tabs1-bp
	   lstarp=fac_cond+(1.-omp)*fac_fus
	   dlstarp=ap*fac_fus
	endif
	fff = tabs(k)-tabs1+lstarn*(qt(k)-qsatt)+lstarp*qp(k)
	dfff=dlstarn*(qt(k)-qsatt)+dlstarp*qp(k)-lstarn*dqsat-1.
	dtabs=-fff/dfff
	niter=niter+1
	tabs1=tabs1+dtabs
      end do   

      qsatt = qsatt + dqsat * dtabs
      qn(k) = max(0.,qt(k)-qsatt)

    else

      qn(k) = 0.

    endif

    tabs(k) = tabs1
    qp(k) = max(0.,qp(k)) ! just in case
    omn = max(0.,min(1.,(tabs(k)-tbgmin)*an))
    qcl(k) = omn*qn(k)
    qci(k) = (1.-omn)*qn(k)
    qv(k)  = qt(k) - qn(k)
    thetali(k) = (tabs(k)- fac_cond*(qcl(k)) - fac_sub*(qci(k)) ) / totheta
    theta(k)   = tabs(k)/totheta
    thetav(k)  = theta(k) * (1.+epsv*qv(k))
    thetar(k)  = thetav(k) - tabs(k)/totheta * qn(k)
    cfrac_pdf(k) = 0.
    if (qn(k).gt.0.0) then
      cfrac_tot(k) = 1.
    else
      cfrac_tot(k) = 0.0
    end if
    sigmas(k)  = 0.

    ! liquid water temp Tl
    tabs1     = totheta * thetali(k)
    ! get saturation mixing ratio at Tl
    omn = max(0.,min(1.,(tabs1-tbgmin)*an))
    qsatt = omn * qsatw(tabs1,pres(k)) + (1.-omn) * qsati(tabs1,pres(k))
    ! get dqsdT
    alphaf = (omn * qsatw(tabs1,pres(k)) * lcond/(tabs1**2*rv) + (1.-omn) * qsati(tabs1,pres(k)) * lsub/(tabs1**2*rv))
    lambdaf = (1. + alphaf * (omn*fac_cond + (1.-omn) * fac_sub) )**(-1.)
    alphaf =alphaf * totheta

    cthl(k) = cfrac_tot(k) * (1. - (fac_cond/totheta - (1.+epsv) * tabs(k)/totheta)*lambdaf * alphaf) & ! see bechtold 95 (a,b,alpha,beta) and my own notes
         + (1.-cfrac_tot(k)) * (1.+epsv*qt(k))
    cqt(k) =  cfrac_tot(k) * (epsv*tabs(k)/totheta + (fac_cond/totheta + (1.+epsv) * tabs(k)/totheta) * lambdaf) &
         + (1.-cfrac_tot(k)) * epsv*thetali(k)


end subroutine cloud

