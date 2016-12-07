module chunking_interface

   integer, parameter :: nf90_hdf5 = 0

   contains

      function nf90_def_var_chunking(ncid, varid, contiguous, chunksizes)
         integer, intent(in) :: ncid
         integer, intent(in) :: varid
         integer, intent(in) :: contiguous
         integer, dimension(:), intent(in) :: chunksizes
         integer :: nf90_def_var_chunking
         nf90_def_var_chunking = 0
      end function nf90_def_var_chunking

end module chunking_interface


module netcdf_mod

use chunking_interface   ! This will be overridden by the netcdf module if it is version 4

implicit none

type netcdf_field
   integer :: varID
   character(50) :: short_name
   character(100) :: long_name
   character(20) :: units
   integer, dimension(:), allocatable :: dimIDs
   integer, dimension(:), allocatable :: location   ! The size of location indicates which of data0 or data1 associated
   logical :: time   ! True if time is one of the field's dimensions, false otherwise
   integer, dimension(:), allocatable :: normalize_dim   ! The dimension along which to normalize a histogram, none if zero
   real,                     pointer :: data0
   real, dimension(:),       pointer :: data1
   real, dimension(:,:),     pointer :: data2
end type netcdf_field

type netcdf_data
   integer :: ncid
   character(50) :: title
   character(200) :: file_path
   character(20), dimension(:), allocatable :: dimNames
   character(20), dimension(:), allocatable :: dimUnits
   integer, dimension(:), allocatable :: dimLengths
   integer, dimension(:), allocatable :: dimIDs
   integer, dimension(:), allocatable :: dimVarIDs
   real, dimension(:,:), allocatable :: dimValues   ! First index ranges over 1 to the number of spatial dimensions
                                                    ! Second index ranges along the dimension's axis
   integer :: size_as_float   ! Size of this netcdf_data as a float
   logical :: time   ! True if time is one of the file's dimensions, false otherwise
   integer :: timeCounter   ! Tracks the number of time slices written
   integer :: writePeriod   ! The number of steps between writing time slices to file
   integer :: numWrites   ! Number of times already written to current output file
   logical :: zero_out   ! If true, set all the data to zero after writing to the NetCDF file
   logical :: as_double   ! If true, data is saved in double precision
   real :: factor   ! Divide the data by this factor before writing to the NetCDF file; should be set to 1 by default
   type(netcdf_field), dimension(:), allocatable :: field
end type netcdf_data


contains      

   !========================!
   ! Subroutine netcdf_open !
   !========================!
   subroutine netcdf_open(nc)

      use netcdf
      use grid

      implicit none
   
      ! In
      type(netcdf_data) :: nc

      integer :: ios, i, data_type
      character(8) :: date
      character(10) :: timenow
      character(2000) :: single_line
      character(20000) :: all_lines
      integer :: netcdf_version

      if (nc%as_double) then
         data_type = nf90_double
      else
         data_type = nf90_real
      end if

      !====================!
      ! Create NetCDF file !
      !====================!
      if (.not.dir_writable(nc%file_path(1:(scan(nc%file_path,'/',back=.true.)-1)))) then
         write(*,*)'Error in netcdf_open: Cannot write to directory ',   &
            nc%file_path(1:(scan(nc%file_path,'/',back=.true.)-1))
            stop
      end if
      netcdf_version = first_integer(nf90_inq_libvers())
      if (netcdf_version .eq. 3) then
         call handle_err( nf90_create(path = nc%file_path, cmode = nf90_clobber, ncid = nc%ncid) )
      else if (netcdf_version .eq. 4) then
         call handle_err( nf90_create(path = nc%file_path, cmode = nf90_hdf5, ncid = nc%ncid) )
      else
         write (*,*) 'Error in netcdf_open: netCDF version number must be either 3 or 4'
         stop
      end if
      call handle_err( nf90_put_att(nc%ncid, nf90_global, "title", nc%title) )


      !=======================!
      ! Add the prm attribute !
      !=======================!
      
      ! Add the header to all_lines
      call date_and_time(date,timenow)
      all_lines = '!'
      do i = 1,max(21,len_trim(path)+2)
         all_lines = trim(all_lines)//'='
      end do
      all_lines = trim(all_lines)//'!\n! '//trim(path)
      i = len_trim(all_lines)-len_trim(path)+max(21,len_trim(path)+2)
      all_lines(i:i) = '!'
      all_lines = trim(all_lines)//'\n! '//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '&
              //timenow(1:2)//':'//timenow(3:4)//'.'//timenow(5:6)
      i = len_trim(all_lines)-19+max(21,len_trim(path)+2)
      all_lines(i:i) = '!'
      all_lines = trim(all_lines)//'\n!'
      do i = 1,max(21,len_trim(path)+2)
         all_lines = trim(all_lines)//'='
      end do
      all_lines = trim(all_lines)//'!\n\n'

      ! Read the nml file
      open(8, file=trim(path)//'/prm', status='old', form='formatted', iostat=ios)
      if (ios .gt. 0) then
          write(*,*) 'Error in netcdf_open: Unable to open nml file'
          stop
      end if
      read(8,"(a)",iostat=ios) single_line
      do while (ios .eq. 0)
         all_lines = trim(all_lines)//trim(single_line)//'\n'
         read(8,"(a)",iostat=ios) single_line
      end do
      close(8)
      
      ! Write to the nml attribute
      call handle_err( nf90_put_att(nc%ncid, nf90_global, "prm", trim(all_lines)) )
      

      !===================!
      ! Define dimensions !
      !===================!
      do i = 1,size(nc%dimNames)
         call handle_err( nf90_def_dim(nc%ncid, nc%dimNames(i), nc%dimLengths(i), nc%dimIDs(i)) )
         call handle_err( nf90_def_var(nc%ncid, nc%dimNames(i), nf90_double, (/ nc%dimIDs(i) /), nc%dimVarIDs(i)) )
         call handle_err( nf90_put_att(nc%ncid, nc%dimVarIDs(i), "units", nc%dimUnits(i)) )
      end do


      !==================!
      ! Define variables !
      !==================!
      do i = 1,size(nc%field)
         if (allocated(nc%field(i)%dimIDs)) then
            call handle_err( nf90_def_var(nc%ncid, nc%field(i)%short_name, data_type, &
               nc%field(i)%dimIDs, nc%field(i)%varID) )
            if (associated(nc%field(i)%data2)) then
               call handle_err( nf90_def_var_chunking(ncid=nc%ncid,varid=nc%field(i)%varID,contiguous=0, &
                  chunksizes=(/shape(nc%field(i)%data2),1/)) )
            else if (associated(nc%field(i)%data1)) then
               call handle_err( nf90_def_var_chunking(ncid=nc%ncid,varid=nc%field(i)%varID,contiguous=0, &
                  chunksizes=(/shape(nc%field(i)%data1),1/)) )
            end if
         else
            call handle_err( nf90_def_var(nc%ncid, nc%field(i)%short_name, data_type, &
               nc%field(i)%varID) )
         end if
         call handle_err( nf90_put_att(nc%ncid, nc%field(i)%varID, "units", nc%field(i)%units) )
         call handle_err( nf90_put_att(nc%ncid, nc%field(i)%varID, "long_name", nc%field(i)%long_name) )
      end do


      !=================!
      ! End define mode !
      !=================!
      call handle_err( nf90_enddef(nc%ncid) )
    

      !=====================!
      ! Define spatial axes !
      !=====================!
      if (nc%time) then
         do i = 1,size(nc%dimNames)-1
          !  if (nc%as_double) then
               call handle_err( nf90_put_var(nc%ncid, nc%dimIDs(i), real(nc%dimValues(i,1:nc%dimLengths(i)),double_bytes)) )
          !  else
          !     call handle_err( nf90_put_var(nc%ncid, nc%dimIDs(i), real(nc%dimValues(i,1:nc%dimLengths(i)),real_bytes)) )
          !  end if
        end do
      else
         do i = 1,size(nc%dimNames)
          !  if (nc%as_double) then
               call handle_err( nf90_put_var(nc%ncid, nc%dimIDs(i), real(nc%dimValues(i,1:nc%dimLengths(i)),double_bytes)) )
          !  else
          !     call handle_err( nf90_put_var(nc%ncid, nc%dimIDs(i), real(nc%dimValues(i,1:nc%dimLengths(i)),real_bytes)) )
          !  end if
         end do
      end if   

      ! Sync the file
      call netcdf_sync(nc)

   end subroutine netcdf_open


   !========================!
   ! Subroutine netcdf_sync !
   !========================!
   subroutine netcdf_sync(nc)
     
      use netcdf
      use grid

      implicit none

      ! In
      type(netcdf_data) :: nc

      call handle_err( nf90_sync(nc%ncid) )

   end subroutine netcdf_sync


   !=================================!
   ! Subroutine netcdf_write_allproc !
   !=================================!
   subroutine netcdf_write_allproc(nc)

      use grid, only : real_bytes, double_bytes

      implicit none

      ! In/Out
      type(netcdf_data) :: nc

      integer :: i, j, request, dummy_rank, dummy_tag
      integer, dimension(4) :: s
      real, dimension(:,:,:), allocatable :: summation


      ! Update timeCounter
      if (nc%time) nc%timeCounter = nc%timeCounter + 1

      ! Divide by factor
      if (nc%factor .ne. 1.) then
         do i = 1,size(nc%field) 
            if (associated(nc%field(i)%data2)) then
               nc%field(i)%data2 = nc%field(i)%data2/nc%factor 
            else if (associated(nc%field(i)%data1)) then
               nc%field(i)%data1 = nc%field(i)%data1/nc%factor
            else if (associated(nc%field(i)%data0)) then
               nc%field(i)%data0 = nc%field(i)%data0/nc%factor
            end if
         end do
      end if
    

      ! Write masterproc's data
      call netcdf_write_singleproc(nc)

      ! Zero out the data
      if (nc%zero_out) then
         do i = 1,size(nc%field) 
            if (associated(nc%field(i)%data2)) then
               nc%field(i)%data2 = 0.
            else if (associated(nc%field(i)%data1)) then
               nc%field(i)%data1 = 0.
            else if (associated(nc%field(i)%data0)) then
               nc%field(i)%data0 = 0.
            end if
         end do
      end if

      ! Increment numWrites
      !nc%numWrites = nc%numWrites + 1

      ! Sync the file
      call netcdf_sync(nc)

   end subroutine netcdf_write_allproc


   !====================================!
   ! Subroutine netcdf_write_singleproc !
   !====================================!
   subroutine netcdf_write_singleproc(nc)
   
      use netcdf
      use grid, only: real_bytes, double_bytes, time
   
      implicit none
      
      ! In/Out
      type(netcdf_data) :: nc
   
      ! Dummy
      integer :: i
      integer, dimension(:), allocatable :: start


      ! Set the corresponding time variable to time
      if (nc%time) then
         if (nc%as_double) then
            call handle_err( nf90_put_var(nc%ncid, nc%dimVarIDs(size(nc%dimVarIDs)), &
                 real(time,double_bytes), start=(/ nc%timeCounter /)) )
         else
            call handle_err( nf90_put_var(nc%ncid, nc%dimVarIDs(size(nc%dimVarIDs)), &
                 real(time,real_bytes), start=(/ nc%timeCounter /)) )
         end if
      end if


      ! Write the time-independent fields   !================================================!    
      if (nc%timeCounter .eq. 0) then       ! If timeCounter equals 0, that means that there !    
                                            ! is no time dimension, so write all fields      ! 
         do i = 1,size(nc%field)            !================================================!            
            if (allocated(nc%field(i)%location)) then   ! Either data1, data2, data3, or data4 ...
               if (size(nc%field(i)%location) .eq. 2) then
                  if (nc%as_double) then
                     call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data2,double_bytes), &   ! data2
                        start=nc%field(i)%location) )
                  else
                     call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data2,real_bytes), &   ! data2
                        start=nc%field(i)%location) )
                  end if
               else if (size(nc%field(i)%location) .eq. 1) then 
                  if (nc%as_double) then
                     call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data1,double_bytes), &   ! data1
                        start=nc%field(i)%location) )
                  else
                     call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data1,real_bytes), &   ! data1
                        start=nc%field(i)%location) )                  
                  end if
               else
                  print *,'Error: In netcdf_write, size of location is ',size(nc%field(i)%location)
                  stop
               end if
            else                                        ! ... or else data0
               if (nc%as_double) then
                  call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data0,double_bytes)) )       ! data0
               else
                  call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data0,real_bytes)) )       ! data0               
               end if
            end if
         end do  
                                               !==========================================!
      else if (nc%timeCounter .eq. 1) then     ! If timeCounter equals 1, that means that !
                                               ! this is the first call, so write all     !
         do i = 1,size(nc%field)               ! fields that do not depend on time        !
            if (.not. nc%field(i)%time) then   !==========================================!

               if (allocated(nc%field(i)%location)) then   ! Either data4, data3, data2, or data1 ...    

                  if (size(nc%field(i)%location) .eq. 2) then
                     if (nc%as_double) then
                        call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data2,double_bytes), &   ! data2
                           start=nc%field(i)%location) )
                     else
                        call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data2,real_bytes), &   ! data2
                           start=nc%field(i)%location) )
                     end if
                  else if (size(nc%field(i)%location) .eq. 1) then 
                     if (nc%as_double) then             
                        call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data1,double_bytes), &   ! data1
                           start=nc%field(i)%location) )
                     else
                        call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data1,real_bytes), &   ! data1
                           start=nc%field(i)%location) )                     
                     end if
                  else
                     print *,'Error: In netcdf_write, size of location is ',size(nc%field(i)%location)
                     stop
                  end if

               else                                        ! ... or else data0

                  if (nc%as_double) then
                     call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data0,double_bytes)) )
                  else
                     call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data0,real_bytes)) )
                  end if
                  
               end if
            end if
         end do  
      
      end if


      ! Write the time-dependent fields
      do i = 1,size(nc%field)   ! Cycle through the fields

         if (nc%field(i)%time) then   ! Only write time-dependent fields

            ! Allocate start and set equal to location + time
            if (allocated(nc%field(i)%location)) then   ! If there are spatial dimensions, add them plus time to start ...
               allocate(start(size(nc%field(i)%location)+1))
               start(1:(size(start)-1)) = nc%field(i)%location
               start(size(start)) = nc%timeCounter
            else if (nc%field(i)%time .and. .not. allocated(nc%field(i)%location)) then   ! ... otherwise, add just time to start
               allocate(start(1))
               start(1) = nc%timeCounter
            end if     
            
            ! Write the time-dependent field to the NetCDF file      
            if (size(start) .eq. 3) then
               if (nc%as_double) then
                  call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data2,double_bytes), &   ! data2
                     start=start) )
               else
                  call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data2,real_bytes), &   ! data2
                     start=start) )
               end if
            else if (size(start) .eq. 2) then 
               if (nc%as_double) then            
                  call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data1,double_bytes), &   ! data1
                     start=start) )
               else
                  call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data1,real_bytes), &   ! data1
                     start=start) )               
               end if
            else if (size(start) .eq. 1) then
               if (nc%as_double) then
                  call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data0,double_bytes), &   ! data0
                     start=start) )
               else
                  call handle_err( nf90_put_var(nc%ncid, nc%field(i)%varID, real(nc%field(i)%data0,real_bytes), &   ! data0
                     start=start) )               
               end if
            else
               print *,'Error: In netcdf_write, size of location is ',size(nc%field(i)%location)
               stop
            end if
            
            deallocate(start)

         end if

      end do

   end subroutine netcdf_write_singleproc
   

   !=========================!
   ! Subroutine netcdf_close !
   !=========================!
   subroutine netcdf_close(nc)
   
      use netcdf

      implicit none
      
      ! In
      type(netcdf_data) :: nc
      
      call handle_err( nf90_close(nc%ncid) )
   
   end subroutine netcdf_close

   !=====================================!
   ! Subroutine netcdf_set_size_as_float !
   !=====================================!
   subroutine netcdf_set_size_as_float(snap)

      implicit none

      ! In
      type(netcdf_data) :: snap

      integer :: i, j

      ! First find the required size of buffer
      j = 0
      do i = 1,size(snap%field)
         j = j + 4   ! To hold location
         if (allocated(snap%field(i)%location)) then
            if (size(snap%field(i)%location) .eq. 2) then
               j = j + size(snap%field(i)%data2)
            else if (size(snap%field(i)%location) .eq. 1) then
               j = j + size(snap%field(i)%data1)
            end if
         else
            j = j + 1
         end if
      end do

      snap%size_as_float = j

   end subroutine netcdf_set_size_as_float

   
   
   !=======================!
   ! Subroutine handle_err !
   !=======================!
   subroutine handle_err(status)

      use netcdf

      implicit none
      
      ! In
      integer :: status
      
      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop
      end if
      
   end subroutine handle_err


   !==============================================!
   ! Function first_integer                       !
   !                                              !
   ! This function of a character string returns  !
   ! the first integer that appear in the string. !
   !==============================================!
   function first_integer(string)

      implicit none

      ! In
      character(*), intent(in) :: string

      ! Out
      integer :: first_integer

      integer :: i

      i = 1
      do while (i .le. len(trim(string)) .and. index('0123456789',string(i:i)) .eq. 0)
         i = i + 1
      end do

      if (i .le. len(trim(string))) then
         first_integer = index('0123456789',string(i:i)) - 1
      else
         first_integer = 0
      end if

   end function first_integer

end module netcdf_mod
