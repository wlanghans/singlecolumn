
subroutine tke_full()

!	this subroutine solves the TKE equation

use vars
use params
use grid
implicit none

real grd,betdz,Ck,Ce,Ces,Ce1,Ce2,Pr,Cee,Cs
real ratio,a_prod_sh,a_prod_bu,a_diss
real lstarn, lstarp, bbb, omn, omp, tketau
real qsatt,dqsat, dtkedtsum, dtkedtmin
integer i,j,k,kc,kb
real, dimension(nz) :: tkelediss, tkesbdiss, tkesbshear, tkesbbuoy

! call t_startf('tke_full')

!Cs = 0.1944
Cs = 0.15
Ck=0.1
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
  tkelediss(k) = 0.
  tkesbdiss(k) = 0.
  tkesbshear(k)= 0.
  tkesbbuoy(k) = 0.

!  SGS buoyancy flux

!bloss: removed temperature diagnostics for omn.
!         - use mass weighted qsat, dqsat and latent heat for cloud 
!         - separate buoyancy contributions for precipitating water and ice.
   
   
     if(qcl(k)+qci(k) .gt. 0.) then
      
       omn = qcl(k)/(qcl(k)+qci(k)+1.e-20)
       lstarn = fac_cond+(1.-omn)*fac_fus
      
       dqsat = omn*dtqsatw(tabs(k),pres(k))+ &
                             (1.-omn)*dtqsati(tabs(k),pres(k))
       qsatt = omn*qsatw(tabs(k),pres(k))+(1.-omn)*qsati(tabs(k),pres(k))
       bbb = 1. + epsv*qsatt-qcl(k)-qci(k) -qpl(k)-qpi(k)+1.61*tabs(k)*dqsat
       bbb = bbb / (1.+lstarn*dqsat)
       buoy_sgs(k)=betdz*(bbb*(t(kc)-t(kb)) &
         +(bbb*lstarn - (1.+lstarn*dqsat)*tabs(k))* &
             (qv(kc)+qcl(kc)+qci(kc)-qv(kb)-qcl(kb)-qci(kb)) & 
         + (bbb*fac_cond - (1.+fac_cond*dqsat)*tabs(k))*(qpl(kc)-qpl(kb))  &
         + (bbb*fac_sub  - (1.+fac_sub *dqsat)*tabs(k))*(qpi(kc)-qpi(kb)) )
!bloss   +(bbb*lstarp - (1.+lstarp*dqsat)*tabs(i,j,k))* &
!bloss            (qpl(i,j,kc)+qpi(i,j,kc)-qpl(i,j,kb)-qpi(i,j,kb)) )
     else
      
        bbb = 1.+epsv*qv(k)-qpl(k)-qpi(k)
        buoy_sgs(k)=betdz*( bbb*(t(kc)-t(kb)) &
          +epsv*tabs(k)* &
           (qv(kc)+qcl(kc)+qci(kc)-qv(kb)-qcl(kb)-qci(kb)) &
          +(bbb*fac_cond-tabs(k))*(qpl(kc)-qpl(kb)) &
          +(bbb*fac_sub -tabs(k))*(qpi(kc)-qpi(kb)) )
!bloss    +(bbb*lstarp-tabs(i,j,k))* &
!bloss         (qpl(i,j,kc)+qpi(i,j,kc)-qpl(i,j,kb)-qpi(i,j,kb)) )
      end if
       !WL use simpler formulation for N**2
       buoy_sgs(k)=ggr/thetav(k) * (thetav(kc)-thetav(kb))/ (z(kc)-z(kb))

   if(buoy_sgs(k).le.0.) then
     smix(k)=grd
   else
     ! WL replaced SAM's formulation with the equivalent formulation based on tke, which is l=0.76 * (e/Nm^2)^(1/2)
     ! WL this way we don't need to store tk for the next timestep, but only tke 
     ! WL also, TKE will have been advected and thus provides a better esimate than the local Km from the previous step
     ! smix=min(grd,max(0.1*grd, sqrt(0.76*tk(i,j,k)/Ck/sqrt(buoy_sgs+1.e-10))))
     smix(k)=min(grd,max(0.1*grd, 0.76*tke(k)/sqrt(buoy_sgs(k)+1.e-10))) 
   end if

   if (dolteix) then 
      if (fixedtau) then
         tketau=600.
      else
         tketau= max(0.5 * pblh /  ((ggr/300.*sgs_sens_heat_flux(1)*pblh)**(1./3.)),0.0)
      end if
      smix(k)=  tketau*sqrt(tke(k)) +(xkar*z(k)-tketau*sqrt(tke(k)))*exp(-z(k)/100.)
   end if
 

   ratio=smix(k)/grd
   Pr=1. 
!   Pr=1. +2.*ratio
   Cee=Ce1+Ce2*ratio

   if (doconsttk) then

       tk(k)=tkconst
       tke(k) = (tk(k)/(Ck*smix(k)))**2
       a_prod_sh=(tk(k)+0.001)*def2(k)
       a_prod_bu=-(tk(k)+0.001)*Pr*buoy_sgs(k)
       a_diss=a_prod_sh+a_prod_bu

   else
     if(dosmagor) then

       tk(k)=sqrt(Ck**3/Cee*max(0.,def2(k)-Pr*buoy_sgs(k)))*smix(k)**2
       tke(k) = (tk(k)/(Ck*smix(k)))**2
       a_prod_sh=(tk(k)+0.001)*def2(k)
       a_prod_bu=-(tk(k)+0.001)*Pr*buoy_sgs(k)
       a_diss=a_prod_sh+a_prod_bu

     else

       ! get Km 
       tk(k) = 0.5*smix(k)*sqrt(tke(k))
       a_prod_sh=(tk(k)+0.001)*def2(k)
       a_prod_bu=-(tk(k)+0.001)*Pr*buoy_sgs(k)
       !a_diss=min(tke(k)/(4.*dt),Cee/smix(k)*tke(k)**1.5) ! cap the diss rate (useful for large time steps
       !a_diss=min(tke(k)/(4.*dt),0.16*2.5 / smix(k)*tke(k)**1.5) ! cap the diss rate (useful for large time steps
       a_diss=0.16*2.5 / smix(k)*tke(k)**1.5 ! cap the diss rate (useful for large time steps
       !a_prod_bu = max(a_prod_bu,-a_prod_sh)
       dtkedtsum = a_prod_sh+a_prod_bu-a_diss
       dtkedtmin = -tke(k)/dt
       if (dtkedtsum.lt.dtkedtmin) then
         dtkedtsum = dtkedtmin / dtkedtsum
         a_prod_bu = dtkedtsum * a_prod_bu
         a_prod_sh = dtkedtsum * a_prod_sh
         a_diss    = dtkedtsum * a_diss
       end if
       ! this is a preliminary euler step to get K(n+1), which is used for mixing 
       ! we will later advance rhotke if progtke=true
       tke(k)=max(0.0,tke(k)+dt*(a_prod_sh+a_prod_bu-a_diss))

       tk(k)=Ck*smix(k)*sqrt(tke(k))
       !tk(k)=0.5*smix(k)*sqrt(tke(k))

     end if
   end if
	
   tkh(k)=Pr*tk(k)

   if (progtke) then
     tend_rho_tke_buoy (k)   =    rho(k) * a_prod_bu
     tend_rho_tke_shear(k)   =    rho(k) * a_prod_sh
     tend_rho_tke_diss (k)   =  - rho(k) * a_diss 
   end if

   ! apply stability limiter if fully explicit
   if (betam.eq.1.) then
     tk(k)  =min(0.5 * grd**2/dt, tk(k)  )
     tkh(k) =min(0.5 * grd**2/dt, tkh(k)  )
   end if

end do ! k

! call t_stopf('tke_full')

end

