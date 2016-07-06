subroutine add_tendencies()

use vars
use params

implicit none

  rho          = rho       + dt*(tend_sgs_rho + tend_sub_rho+tend_force_rho)
  rhoqv        = rhoqv     + dt*(tend_sgs_rho_qv + tend_sub_rho_qv+tend_force_rho_qv)
  rhothv       = rhothv    + dt*(tend_sgs_rho_thetav + tend_sub_rho_thetav+ &
                           tend_force_rho_thetav+tend_rad_rho_thetav)
  rhou         = rhou      + dt*(tend_sgs_rho_u + tend_sub_rho_u+tend_force_rho_u)
  rhov         = rhov      + dt*(tend_sgs_rho_v + tend_sub_rho_v+tend_force_rho_v)

  if (progtke) then
    rhotke       = rhotke    + dt*(tend_sgs_rho_tke &
                 + tend_rho_tke_buoy + tend_rho_tke_shear + tend_rho_tke_diss)
  else
    ! tke has been diagnosed
    rhotke       = rho * tke 
  end if
  rhotke = rhotke + dt * (tend_sub_rho_tke + tend_force_rho_tke)

end subroutine add_tendencies
