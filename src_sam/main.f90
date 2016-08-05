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


!initialize netCDF output files
call snapshot_init()

!beginning of timestepping
do while(nstep.lt.nstop)

nstep = nstep + 1
time = time + dt

write(*,*) 'Working on timestep ', nstep

! ======================================= 
! set all tendecies to zero and clip non-negative scalars
! copy prog varaibles at beginning of time step
! ======================================= 
 call zero_stuff()


! ======================================= 
! compute tendencies from subsidence 
! (subsidence vel has to be available on faces)
! ======================================= 
 if (dosubsidence) then
   call subsidence()
 end if

! ======================================= 
! compute tendencies from large scale forcing
! coriolis and dqvdt, dthetadt
! ======================================= 
   call forcing()

! ======================================= 
! get drag/transfer coefficients, explicit fluxes, surface values, and wind speed on 1st model level
! ======================================= 
  if (dosurface) then
       !expects rho to be air density and qv to be mass fraction
      call surface()
  end if

! ======================================= 
  ! get PBL height
! ======================================= 
  if (dopblh) call get_pblh() 

! ======================================= 
  ! call edmf 
! ======================================= 
  ! needs pblh and surface fluxes as input
  if (dosgs.and.doedmf) then
      call edmf()
  end if

! ======================================= 
  ! get eddy-diffusivities and tke
! ======================================= 
  if (dosgs) call get_def2()
  if (dosgs) call get_eddyvisc()

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
! ======================================= 
  call add_tendencies()

! ======================================= 
! call PDF condensation scheme to get sgs
! cloud cover and condensed water
! iterative solve also provides tabs, theta, thetav, etc.
! ======================================= 
  if (dosgscloud) then
    call sgscloud2()
  end if

! ======================================= 
! get pressure from eos
! input needed: tabs, qv, and thetav 
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
