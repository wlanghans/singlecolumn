
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


if (doteixpbl.or.dolanghanspbl) then
  Ck=0.5
elseif (dowitekpbl) then
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
         if (doteixpbl)  tketau=400.
         if (dowitekpbl) tketau=600.
      else
         tketau= max(ctketau * pblh /  ((ggr/thetav(1)*sgs_thv_flux(1)*pblh)**(1./3.)),0.0)
      end if
      l23 =    (tketau*sqrt(tke(k)+1.0d-10))**(-1) 
      if (buoy_sgs(k).gt.0.0) l23 = l23 +&
        (max(0.7*sqrt(tke(k)+1.0d-10)/sqrt(buoy_sgs(k)),adz(k)*dz/2.))**(-1)
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
   if (dowitekpbl) then
     Pr = 0.5882 
   else
     Pr=1. 
   end if

!   Pr=1. +2.*ratio
   if (doteixpbl.or.dolanghanspbl) then
      Cee= 0.16 * 2.5
   elseif (dowitekpbl) then
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
     !wthl = -Pr * tk(k) * (t(kc)-t(kb))  &
     !    * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     !if (qp(k).gt.0.0) then
     !  wthl = wthl + ((t(kc)-t(kb)) + (lcond*qpl(k)+lsub*qpi(k))/qp(k) * (qp(kc)-qp(kb)) ) &
     !    * (p00/pres(k))**(rgas/cp) / (cp*(z(kc)-z(kb)))
     !end if
     !  
     !wqt  = -Pr * tk(k) * (qt(kc)-qt(kb)) / (z(kc)-z(kb))

     ! use previously evaluated fluxes
     wthl = (0.5*(t_flux_ed(k) + t_flux_ed(k+1))) / cp 
     if (qp(k).gt.0.0) then
       wthl = wthl +  ((lcond*qpl(k)+lsub*qpi(k))/qp(k) * 0.5*(qp_flux_ed(k) + qp_flux_ed(k+1)))  &
              * (p00/pres(k))**(rgas/cp) / cp
     end if
     wqt  = 0.5*(qt_flux_ed(k) + qt_flux_ed(k+1))

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
 
     ! use previously evaluated fluxes
     wthl = (0.5*(t_flux_ed(k) + t_flux_ed(k+1))) / cp 
     if (qp(k).gt.0.0) then
       wthl = wthl +  ((fac_cond*qpl(k)+fac_sub*qpi(k))/qp(k) * 0.5*(qp_flux_ed(k) + qp_flux_ed(k+1))) * (p00/pres(k))**(rgas/cp)
     end if
     wqt  = 0.5*(qt_flux_ed(k) + qt_flux_ed(k+1))

     ! get buoyancy flux from PDF scheme
     tke_thvflx(k) = cthl(k) * wthl + cqt(k) * wqt

   elseif (doteixpbl) then
   ! ==================================
   ! Use Teixeira CPBL closure; TKE is prognostic 
   ! ==================================

     ! use previously evaluated fluxes
     wthl = (0.5*(t_flux_ed(k) + t_flux_ed(k+1))) / cp 
     if (qp(k).gt.0.0) then
       wthl = wthl +  ((fac_cond*qpl(k)+fac_sub*qpi(k))/qp(k) * 0.5*(qp_flux_ed(k) + qp_flux_ed(k+1))) * (p00/pres(k))**(rgas/cp)
     end if
     wqt  = 0.5*(qt_flux_ed(k) + qt_flux_ed(k+1))

     ! get buoyancy flux from PDF scheme
     tke_thvflx(k) = cthl(k) * wthl + cqt(k) * wqt

     ! get Km 
     tk(k) = Ck*smix(k)*sqrt(tke(k))

     a_prod_sh=(tk(k)+0.001)*def2(k)
     a_prod_bu= ggr/thetar(k) * tke_thvflx(k) !-(tk(k)+0.001) * Pr * buoy_sgs(k) !        ggr/thetar(k) * tke_thvflx(k)
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

     ! use previously evaluated fluxes
     wthl = (0.5*(t_flux_ed(k) + t_flux_ed(k+1)+t_flux_mf(k) + t_flux_mf(k+1))) / cp 
     if (qp(k).gt.0.0) then
       wthl = wthl +  ((lcond*qpl(k)+lsub*qpi(k))/qp(k) * 0.5*(qp_flux_ed(k) + qp_flux_ed(k+1)+qp_flux_mf(k) + qp_flux_mf(k+1)))  &
              * (p00/pres(k))**(rgas/cp) / cp
     end if
     wqt  = 0.5*(qt_flux_ed(k) + qt_flux_ed(k+1)+qt_flux_mf(k) + qt_flux_mf(k+1))

     ! get buoyancy flux from PDF scheme
     tke_thvflx(k) = cthl(k) * wthl + cqt(k) * wqt

     ! get Km 
     tk(k) = Ck*smix(k)*sqrt(tke(k))
     a_prod_sh=(tk(k)+0.001)*def2(k)


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

    ! use previously evaluated fluxes
     wthl = (0.5*(t_flux_ed(k) + t_flux_ed(k+1))) / cp
     if (qp(k).gt.0.0) then
       wthl = wthl +  ((fac_cond*qpl(k)+fac_sub*qpi(k))/qp(k) * 0.5*(qp_flux_ed(k) + qp_flux_ed(k+1))) * (p00/pres(k))**(rgas/cp)
     end if
     wqt  = 0.5*(qt_flux_ed(k) + qt_flux_ed(k+1))

     ! get buoyancy flux from PDF scheme for environment 
     tke_thvflx(k) = (1.-0.5*(frac_mf(k)+frac_mf(k+1) )) * (cthl(k) * wthl + cqt(k) * wqt)
     ! now add buoyancy flux in plumes
     if (k.eq.1) then
       tke_thvflx(k) = tke_thvflx(k) + 0.5 * (0.5*(frac_mf(k)+frac_mf(k+1))*  t_flux_ed(1)/cp+sumMthv(2))
     else
       tke_thvflx(k) = tke_thvflx(k) + 0.5* (sumMthv(k)+sumMthv(k+1))
     end if

     ! get Km 
     tk(k) = Ck*smix(k)*sqrt(tke(k))

     ! production from mean flow
     a_prod_sh=(tk(k)+0.001)*def2(k)
     ! add production from plumes
     if ((frac_mf(k)+frac_mf(k+1) ).gt.0.0) then
       a_prod_sh=a_prod_sh + &
      ( Ck*pblh*sqrt( (tke_mf(k)+tke_mf(k+1))/(frac_mf(k)+frac_mf(k+1) ) ) +0.001)*0.5*(sumDEF2(k)+sumDEF2(k+1))
     end if
     a_prod_bu= ggr/thetar(k) * tke_thvflx(k) !-(tk(k)+0.001) * Pr * buoy_sgs(k) !        ggr/thetar(k) * tke_thvflx(k)
     ! dissipation in environment
     a_diss= Cee / smix(k)*tke(k)**1.5
     !a_diss=(1.-0.5*(frac_mf(k)+frac_mf(k+1) )) * Cee / smix(k)*tke(k)**1.5
     ! add dissipation in plumes (which is prop to sum of tke's)
     !if ((frac_mf(k)+frac_mf(k+1) ).gt.0.0  ) then
     !  a_diss=a_diss + 0.5*(frac_mf(k)+frac_mf(k+1) ) *&
     !  Cee / smix(k)*(tke(k)+0.5*(tke_mf(k)+tke_mf(k+1))/(0.5*(frac_mf(k)+frac_mf(k+1) )))**1.5 
     !end if
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


   elseif (dotkeles) then
   ! ==================================
   ! Use 1.5 TKE closure; TKE is prognostic 
   ! ==================================
   ! same as above but different length scale 
   ! and different coefficients
     tk(k) = Ck*smix(k)*sqrt(tke(k))

     ! use previously evaluated fluxes
     wthl = (0.5*(t_flux_ed(k) + t_flux_ed(k+1)))/ cp 
     if (qp(k).gt.0.0) then
       wthl = wthl +  ((lcond*qpl(k)+lsub*qpi(k))/qp(k) * 0.5*(qp_flux_ed(k) + qp_flux_ed(k+1)))  &
              * (p00/pres(k))**(rgas/cp) / cp
     end if
     wqt  = 0.5*(qt_flux_ed(k) + qt_flux_ed(k+1))
  
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

