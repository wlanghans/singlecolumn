program singlecolumn

use grid
use params
use vars
use snapshot
use advect_mod

implicit none

!set all parameters and constants
call setparm()

!allocate memory
call alloc()

!initialize vertical grid
call setgrid()

!initialize variables
call setdata()

! ======================================= 
! call PDF condensation scheme to get sgs
! cloud cover and condensed water in stratiform part
! iterative solve also provides tabs, theta, thetav, etc.
! in case of dosgscloud=false, no sgs cloud is considered
! ======================================= 
call sgscloudmf(.false.)  ! if false, then the domain mean qt and t are used to get cfrac, qcl, ...
                          ! if true, then the environmental qt and t are used to get cfrac, qcl, ...

! ======================================= 
! recompute hydrostatic pressure, since thetav might have changed 
! under presence of sgs clouds
! input needed: presi(1),thetav 
! ======================================= 
call get_pressure()


!initialize netCDF output files
call snapshot_init()

!beginning of timestepping
do while(nstep.lt.nstop)

nstep = nstep + 1
time = time + dt

write(*,*) 'Working on timestep ', nstep

! ======================================= 
! set all tendecies to zero and clip non-negative scalars
! save prog variables at beginning of time step
! ======================================= 
  call zero_stuff()

! ======================================= 
! compute tendencies from large scale forcing
! coriolis and dqtdt, dthetadt
! ======================================= 
  call forcing()

! ======================================= 
! get drag/transfer coefficients, explicit fluxes, surface values, and wind speed on 1st model level
! ======================================= 
  if (dosurface) then
      call surface()
  end if

! ======================================= 
  ! get PBL height
! ======================================= 
  ! needs thetav as input
  if (dopblh) call get_pblh() 

! ======================================= 
  ! call edmf 
! ======================================= 
  ! needs pblh, surface fluxes, and 
  ! domain-mean density pot temp (for buoyancy) as input
  if (dosgs.and.doedmf) then
      call edmf()
  end if

! ======================================= 
  ! get eddy-diffusivities and tke
  ! needs domain-mean virtual pot temperature to get static stability
! ======================================= 
  if (dosgs) call get_def2()
  if (dosgs) call get_eddyvisc()

! ======================================= 
! compute tendencies from subsidence 
! ======================================= 
  if (dosubsidence) then
   call subsidence()
  end if

! ======================================= 
  ! get sgs/diffusion tendencies (tridiagonal solver)
! ======================================= 
  if (dosurface.or.dosgs) then
       call sgs_tendencies()
  end if

! ======================================= 
  ! call radiative transfer scheme
! ======================================= 
  if(dolongwave.or.doshortwave) then
     call radiation()
  end if

! ======================================= 
  ! apply tendencies to prog. variables
  ! note that if dosequential=t, then tendencies are based 
  ! on atm. state that resulted from previously called
  ! modules
! ======================================= 
  call add_tendencies()

! ======================================= 
! call PDF condensation scheme to get sgs
! cloud cover and condensed water in stratiform part
! iterative solve also provides tabs, theta, thetav, etc.
! in case of dosgscloud=false, no sgs cloud is considered
! ======================================= 
  call sgscloudmf(.false.)  ! if false, then the domain mean qt and t are used to get cfrac, qcl, ...
                              ! if true, then the environmental qt and t are used to get cfrac, qcl, ...

! ======================================= 
! get hydrostatic pressure
! input needed: presi(1),thetav 
! ======================================= 
  call get_pressure()

  call snapshot_update()

end do ! 
! end of timestepping

call snapshot_close()

!deallocate memory
call dealloc()

write(*,*) 'Done'

end program singlecolumn
