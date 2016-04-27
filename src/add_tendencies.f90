subroutine add_tendencies()

use vars
use params

implicit none

  rho          = rho       + dt*tend_sgs_rho
  rhoqv        = rhoqv     + dt*tend_sgs_rho_qv
  rhothv       = rhothv    + dt*tend_sgs_rho_thetav
  rhou         = rhou      + dt*tend_sgs_rho_u
  rhov         = rhov      + dt*tend_sgs_rho_v
  if (progtke) then
    rhotke       = rhotke    + dt*(tend_sgs_rho_tke &
                 + tend_rho_tke_buoy + tend_rho_tke_shear + tend_rho_tke_diss)
  else
    ! tke has been diagnosed
    rhotke       = rho * tke
  end if

end subroutine add_tendencies
