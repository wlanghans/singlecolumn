module grid


implicit none

integer, parameter :: real_bytes = 4
integer, parameter :: double_bytes = 8

integer :: nzm = 60 
integer :: nz  = 61 
integer,parameter :: nx  = 1 
integer,parameter :: ny  = 1 

real,allocatable,dimension(:) :: z, zi, adz, adzw

integer :: nstep = 0
real :: time = 0.
integer:: nrad =1
real :: dt 
real :: dz
integer :: nstop =0
logical :: doconstdz  = .true.
character(80) case


logical :: compute_reffc = .false.
logical :: compute_reffi = .false.

real :: day0=0.  ! starting day (including fraction)
real :: day        ! current day (including fraction)

character(80) :: path   ! The relative path from  main directory to the run directory


contains

      !=======================!
      ! Function dir_writable !
      !=======================!
      function dir_writable(path)

         implicit none

         ! In
         character(*), intent(in) :: path

         ! Out
         logical :: dir_writable

         integer :: ios

         dir_writable = .true.   ! All processors other than the masterproc get .true.
         open(10,file=trim(path)//'/dam_temporary_file_please_delete.temp',iostat=ios)
         if (ios .eq. 0) then
            dir_writable = .true.
         else
            dir_writable = .false.
         end if
         close(10,status='delete')

      end function dir_writable


end module grid
