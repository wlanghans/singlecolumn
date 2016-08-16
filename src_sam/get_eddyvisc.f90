
subroutine get_eddyvisc()

! this subroutine computes the eddy viscosities and eventually tke tendencies

use vars
use params
use grid
implicit none

real grd,betdz,Ck,Ce,Ces,Ce1,Ce2,Pr,Cee,Cs
real ratio,a_prod_sh,a_prod_bu,a_diss
real lstarn, lstarp, bbb, omn, omp, tketau
real qsatt,dqsat, dtkedtsum, dtkedtmin, l23
real :: thetalt, thetalk, thetall, qtt, qtk, qtl, covarqtthetal, varqt, varthetal, wthl, wqt
integer i,j,k,kc,kb

! call t_startf('tke_full')

if (doteixpbl) then
  Ck=0.5
elseif (dowitekpbl.or.dolanghanspbl) then
  Ck=0.425
else
  Ck=0.1
end if
Cs = 0.15
Ce=Ck**3/Cs**4
Ces=Ce/0.7*3.0	


do k=1,nzm      
  kb=k-1
  kc=k+1

  grd=dz*adz(k)
  tke(k)=max(0.,tke(k))
 
  betdz=ggr/tabs0(k)/dz/(adzw(kc)+adzw(k))
  Ce1=Ce/0.7*0.19
  Ce2=Ce/0.7*0.51
  if(k.eq.1) then
    kb=1
    kc=2
    betdz=ggr/tabs0(k)/dz/adzw(kc)
    Ce1=Ces/0.7*0.19
    Ce2=Ces/0.7*0.51
  end if
  if(k.eq.nzm) then
    kb=nzm-1
    kc=nzm
    betdz=ggr/tabs0(k)/dz/adzw(k)
    Ce1=Ces/0.7*0.19
    Ce2=Ces/0.7*0.51
  end if

   !N**2 based on density potential temperature
   buoy_sgs(k)=ggr/thetav(k) * (thetav(kc)-thetav(kb))/ (z(kc)-z(kb))

   if (doteixpbl.or.dowitekpbl.or.dolanghanspbl) then 
   ! always use Teixeira's mixing length (Witek's is very similar)
      if (fixedtau) then
         if (doteixpbl)  tketau=600.
         if (dowitekpbl) tketau=600.
      else
         tketau= max(0.5 * pblh /  ((ggr/thetav(1)*sgs_thv_flux(1)*pblh)**(1./3.)),0.0)
      end if
      l23 =    (tketau*sqrt(tke(k)))**(-1) 
      if (buoy_sgs(k).gt.0.0.and.dowitekpbl) l23 = l23 + (max(0.7*sqrt(tke(k))/sqrt(buoy_sgs(k)),adz(k)*dz/2.))**(-1)
      l23 = l23**(-1)
      smix(k)=  l23 +(xkar*z(k)-l23)*exp(-z(k)/100.)
   else
     if(buoy_sgs(k).le.0.) then
       smix(k)=grd
     else
       ! WL replaced SAM's formulation with the equivalent formulation based on tke, which is l=0.76 * (e/N^2)^(1/2)
       ! Deardorff 1976
       ! WL this way we don't need to store tk for the next timestep, but only tke 
       ! WL also, TKE will have been advected and thus provides a better esimate than the local Km from the previous step
       ! smix=min(grd,max(0.1*grd, sqrt(0.76*tk(i,j,k)/Ck/sqrt(buoy_sgs+1.e-10))))
       smix(k)=min(grd,max(0.1*grd, 0.76*sqrt(tke(k)/ (buoy_sgs(k)+1.e-10)))) 
     end if
   end if
 
   ! Prantdl number
   if (dowitekpbl.or.dolanghanspbl) then
     Pr = 0.5882 
   else
     Pr=1. 
   end if

!   Pr=1. +2.*ratio
   if (doteixpbl) then
      Cee= 0.16 * 2.5
   elseif (dowitekpbl.or.dolanghanspbl) then
      Cee=0.304 
   else
      ratio=smix(k)/grd
      Cee=Ce1+Ce2*ratio
   end if


   if (doconsttk) then
   ! ==================================
   ! set eddy viscosities to a constant value
   ! ==================================

     tk(k)=tkconst
     a_prod_sh=(tk(k)+0.001)*def2(k)
     a_prod_bu=-(tk(k)+0.001)*Pr*buoy_sgs(k)
     a_diss=a_prod_sh+a_prod_bu

     ! use explicit estimates for hil, qt, and qp fluxes to get thli and qt fluxes
     wthl = -Pr * tk(k) * (t(kc)-t(kb))  &
         * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     if (qp(k).gt.0.0) then
       wthl = wthl + ((t(kc)-t(kb)) + (lcond*qpl(k)+lsub*qpi(k))/qp(k) * (qp(kc)-qp(kb)) ) &
         * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     end if
       
     wqt  = -Pr * tk(k) * (qt(kc)-qt(kb)) / (z(kc)-z(kb))
     ! get buoyancy flux from PDF scheme
     tke_thvflx(k) = cthl(k) * wthl + cqt(k) * wqt

   elseif (dosmagor) then
   ! ==================================
   ! Use Smagorinski closure; TKE is diagnosic
   ! ==================================

     tk(k)=sqrt(Ck**3/Cee*max(0.,def2(k)-Pr*buoy_sgs(k)))*smix(k)**2
     a_prod_sh=(tk(k)+0.001)*def2(k)
     a_prod_bu=-(tk(k)+0.001)*Pr*buoy_sgs(k)
     a_diss=a_prod_sh+a_prod_bu

     ! use explicit estimates for hil, qt, and qp fluxes to get thli and qt fluxes
     wthl = -Pr * tk(k) * (t(kc)-t(kb))  &
         * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     if (qp(k).gt.0.0) then
       wthl = wthl + ((t(kc)-t(kb)) + (lcond*qpl(k)+lsub*qpi(k))/qp(k) * (qp(kc)-qp(kb)) ) &
         * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     end if
     wqt  = -Pr * tk(k) * (qt(kc)-qt(kb)) / (z(kc)-z(kb))
     ! get buoyancy flux from PDF scheme
     tke_thvflx(k) = cthl(k) * wthl + cqt(k) * wqt

   elseif (doteixpbl) then
   ! ==================================
   ! Use Teixeira CPBL closure; TKE is prognostic 
   ! ==================================

     !if (dosgscloud.and. k.gt.1.and.k.lt.nzm) then
     !thetalt = theta(k+1) - fac_cond * qcl(k+1)/(1.+qv(k+1)) * (pres(k+1)/p00)**(-rgas/cp)
     !thetalk = theta(k) - fac_cond * qcl(k)/(1.+qv(k)) * (pres(k)/p00)**(-rgas/cp)
     !thetall = theta(k-1) - fac_cond * qcl(k-1)/(1.+qv(k-1)) * (pres(k-1)/p00)**(-rgas/cp)
     !thetalgrad(k) = 0.5 * ((thetalt-thetalk)/(z(k+1)-z(k)) + (thetalk-thetall)/(z(k)-z(k-1))  )
     !qtgrad(k) = 0.5 * ((qtt-qtk)/(z(k+1)-z(k)) + (qtk-qtl)/(z(k)-z(k-1))  )
     !qtt = (qv(k+1) + qcl(k+1))/(1.+qv(k+1))
     !qtk = (qv(k)   + qcl(k))/(1.+qv(k))
     !qtl = (qv(k-1) + qcl(k-1))/(1.+qv(k-1))
     !thetalflux = -(tk(k)+0.001) * Pr * thetalgrad(k)
     !qtflux = -(tk(k)+0.001) * Pr * qtgrad(k)


     ! taken from Bechtold et al. (1992)
     !varqt = 2. *Ck/1.8*smix(k)**2 * qtgrad(k)**2
     !varthetal = 2. *Ck/2.6*smix(k)**2 * thetalgrad(k)**2
     !covarqtthetal = 2. *Ck/1.2*smix(k)**2 * thetalgrad(k) * qtgrad(k)
     

     ! call sgs cloud scheme
     !call sgscloud(k,qtflux,thetalflux,varqt,varthetal,covarqtthetal,thetalk,qtk,pres(k), & ! input
     !                  cfrac_ed(k),qlsgs_ed(k),tke_thvflx(k))                               ! output
     !call sgscloud2(k,qtflux,thetalflux,thetalk,qtk,pres(k), & ! input
     !                  cfrac_ed(k),qlsgs_ed(k),tke_thvflx(k))                               ! output

     !else
     !  cfrac_pdf(k) = 0.
     !  qlsgs_ed(k) = qcl(k)
     !  tk(k) = Ck*smix(k)*sqrt(tke(k))
     !  tke_thvflx(k) = -(tk(k)+0.001)*Pr*buoy_sgs(k) * thetav(k)/ggr
     !end if


     ! get Km 
     tk(k) = Ck*smix(k)*sqrt(tke(k))

     ! use explicit estimates for hil, qt, and qp fluxes to get thli and qt fluxes
     wthl = -Pr * tk(k) * (t(kc)-t(kb))  &
         * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     if (qp(k).gt.0.0) then
       wthl = wthl + ((t(kc)-t(kb)) + (lcond*qpl(k)+lsub*qpi(k))/qp(k) * (qp(kc)-qp(kb)) ) &
         * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     end if
     wqt  = -Pr * tk(k) * (qt(kc)-qt(kb)) / (z(kc)-z(kb)) 
     ! get buoyancy flux from PDF scheme
     tke_thvflx(k) = cthl(k) * wthl + cqt(k) * wqt

     a_prod_sh=(tk(k)+0.001)*def2(k)
     a_prod_bu=ggr/thetar(k) * tke_thvflx(k)
     a_diss=Cee / smix(k)*tke(k)**1.5 ! cap the diss rate (useful for large time steps
     dtkedtsum = a_prod_sh+a_prod_bu-a_diss
     dtkedtmin = -tke(k)/dt
     if (dtkedtsum.lt.dtkedtmin) then
       dtkedtsum = dtkedtmin / dtkedtsum
       a_prod_bu = dtkedtsum * a_prod_bu
       a_prod_sh = dtkedtsum * a_prod_sh
       a_diss    = dtkedtsum * a_diss
     end if
     if (dosequential) tke(k)=tke(k)+dt*(a_prod_sh+a_prod_bu-a_diss)
     wstar=max(0.,(ggr/thetav(1)*sgs_thv_flux(1)*pblh)**(1./3.))

   elseif (dowitekpbl) then
   ! ==================================
   ! Use Witek CPBL closure; TKE is prognostic 
   ! ==================================

     ! get Km 
     tk(k) = Ck*smix(k)*sqrt(tke(k))
     a_prod_sh=(tk(k)+0.001)*def2(k)

     ! use explicit estimates for hil, qt, and qp fluxes to get thli and qt fluxes
     wthl = -Pr * tk(k) * (t(kc)-t(kb))  &
         * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     if (qp(k).gt.0.0) then
       wthl = wthl + ((t(kc)-t(kb)) + (lcond*qpl(k)+lsub*qpi(k))/qp(k) * (qp(kc)-qp(kb)) ) &
         * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     end if
     wqt  = -Pr * tk(k) * (qt(kc)-qt(kb)) / (z(kc)-z(kb))
     ! get buoyancy flux from PDF scheme
     tke_thvflx(k) = cthl(k) * wthl + cqt(k) * wqt

     a_prod_bu=  ggr/thetar(k) * tke_thvflx(k) 
     a_diss=Cee / smix(k)*tke(k)**1.5 ! cap the diss rate (useful for large time steps
     dtkedtsum = a_prod_sh+a_prod_bu-a_diss
     dtkedtmin = -tke(k)/dt
     if (dtkedtsum.lt.dtkedtmin) then
       dtkedtsum = dtkedtmin / dtkedtsum
       a_prod_bu = dtkedtsum * a_prod_bu
       a_prod_sh = dtkedtsum * a_prod_sh
       a_diss    = dtkedtsum * a_diss
     end if
     if (dosequential) tke(k)=tke(k)+dt*(a_prod_sh+a_prod_bu-a_diss)
     wstar=max(0.,(ggr/thetav(1)*sgs_thv_flux(1)*pblh)**(1./3.))

    elseif (dolanghanspbl) then
   ! ==================================
   ! Use Langhans CPBL closure with sgs clouds; TKE is prognostic 
   ! ==================================

   elseif (dotkeles) then
   ! ==================================
   ! Use 1.5 TKE closure; TKE is prognostic 
   ! ==================================
   ! same as above but different length scale 
   ! and different coefficients
     tk(k) = Ck*smix(k)*sqrt(tke(k))

     ! use explicit estimates for hil, qt, and qp fluxes to get thli and qt fluxes
     wthl = -Pr * tk(k) * (t(kc)-t(kb))  &
         * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     if (qp(k).gt.0.0) then
       wthl = wthl + ((t(kc)-t(kb)) + (lcond*qpl(k)+lsub*qpi(k))/qp(k) * (qp(kc)-qp(kb)) ) &
         * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     end if
     wqt  = -Pr * tk(k) * (qt(kc)-qt(kb)) / (z(kc)-z(kb))
     ! get buoyancy flux from PDF scheme
     tke_thvflx(k) = cthl(k) * wthl + cqt(k) * wqt

     a_prod_sh=(tk(k)+0.001)*def2(k)
     a_prod_bu=-(tk(k)+0.001)*Pr*buoy_sgs(k)
     a_diss=min(tke(k)/(4.*dt),Cee/smix(k)*tke(k)**1.5) ! cap the diss rate (useful for large time steps
     dtkedtsum = a_prod_sh+a_prod_bu-a_diss
     dtkedtmin = -tke(k)/dt
     if (dtkedtsum.lt.dtkedtmin) then
       dtkedtsum = dtkedtmin / dtkedtsum
       a_prod_bu = dtkedtsum * a_prod_bu
       a_prod_sh = dtkedtsum * a_prod_sh
       a_diss    = dtkedtsum * a_diss
     end if
     if (dosequential) tke(k)=tke(k)+dt*(a_prod_sh+a_prod_bu-a_diss)

   else
     write(*,*) 'ERROR in get_eddyvisc: Closure not implemented' 
     stop
   end if

   tkh(k)=Pr*tk(k)

   ! ===================================================
   ! save tendencies for tke
   ! ===================================================
   if (progtke) then
     tend_tke_buoy (k)   =     a_prod_bu
     tend_tke_shear(k)   =     a_prod_sh
     tend_tke_diss (k)   =  -  a_diss 
   else
     ! closure is diagnostic, but tke will still be advected, since the
     ! mixing length at the next step will be computed based on the old
     ! (but advected) tke 
     tend_tke_buoy (k)   =    ((tk(k)/(Ck*smix(k)))**2 &
                                   - tke0(k)) / dt
     tend_tke_shear(k)   =    0.
     tend_tke_diss (k)   =    0.
   end if
 
   ! apply stability limiter if fully explicit
   if (betam.eq.1.) then
     tk(k)  =min(0.5 * grd**2/dt, tk(k)  )
     tkh(k) =min(0.5 * grd**2/dt, tkh(k)  )
   end if

end do ! k

end subroutine get_eddyvisc

