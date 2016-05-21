program singlecolumn

use grid
use params
use vars
use snapshot

implicit none

!set all parameters and constants
call setparm()

!allocate memory
call alloc()

!initialize vertical grid
call setgrid()

!initialize variables
call setdata()
call get_state()


!initialize netCDF output files
call snapshot_init()

!beginning of timestepping
do while(nstep.lt.nstop)

nstep = nstep + 1
time = time + dt

write(*,*) 'Working on timestep ', nstep

! ======================================= 
! set all tendecies to zero and clip non-negative scalars
! ======================================= 
 call zero_stuff()

! ======================================= 
  ! get drag/transfer coefficients, explicit fluxes, surface values, and wind speed on 1st model level
! ======================================= 
  if (dosurface) then
       !expects rho to be air density and qv to be mass fraction
      call surface()
  end if

! ======================================= 
  ! convert mass fractions to mixing ratios
  ! air density to dry air density
  ! and get h/cp
! ======================================= 
 call convert_to_sam_state()

! ======================================= 
  ! get eddy-diffusivities and tke
! ======================================= 
  if (dosgs) call get_def2()
  !requires sam state: rho=rhod and mixing ratios qv=rv
  if (dosgs) call get_eddyvisc()

! ======================================= 
  ! get PBL height
! ======================================= 
  if ((dosgs.and.doedmf).or.dopblh) call get_pblh() ! if true is passed than over sea, over land otherwise

! ======================================= 
  ! call edmf 
! ======================================= 
  if (dosgs.and.doedmf) then
      call edmf()
  end if

! ======================================= 
  ! get sgs/diffusion tendencies (tridiagonal solver)
! ======================================= 
  if (dosurface.or.dosgs) then
       !requires sam state rho=rhod and mixing ratios qv=rv
       ! note that in case of tke, mixing will act on the preliminary tke from tke_full; not sure if this is matters
       call sgs_tendencies()
       ! optional diagnostic call to get dimensionless K's
       !if (dosgs) call get_dimlessdiff()
  end if

! ======================================= 
  ! convert mixing ratios back to mass fractions
  ! dry air density to air density
  ! Note: this reconversion is only needed here in case tke is diagnosed (rhtke=rho*tke in add_tendencies)
! ======================================= 
 call reconvert_sam_to_state()

! ======================================= 
  ! apply tendencies to prog. variables
! ======================================= 
  call add_tendencies()

  call get_state()

  call snapshot_update()

end do ! 
! end of timestepping

call snapshot_close()

!deallocate memory
call dealloc()

write(*,*) 'Done'

end program singlecolumn
