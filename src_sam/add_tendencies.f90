subroutine add_tendencies()

use vars
use params

implicit none

  t         = t0      + dt*(tend_sgs_t  + tend_sub_t +tend_force_t+tend_rad_t )
  qt        = qt0     + dt*(tend_sgs_qt + tend_sub_qt+tend_force_qt)
  qp        = qp0     + dt*(tend_sgs_qp + tend_sub_qp)
  u         = u0      + dt*(tend_sgs_u + tend_sub_u+tend_force_u)
  v         = v0      + dt*(tend_sgs_v + tend_sub_v+tend_force_v)

  if (progtke) then
    tke     = tke0    + dt*(tend_sgs_tke &
                 + tend_tke_buoy + tend_tke_shear + tend_tke_diss)
  else
    ! tke has been diagnosed, nothing to do
    tke       = tke0 + dt*tend_tke_buoy  
  end if
  ! finally advect/force tke if necessary
  tke = tke + dt * (tend_sub_tke + tend_force_tke)

end subroutine add_tendencies
