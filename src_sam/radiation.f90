
subroutine radiation()

!	Radiation interface

use params, only: doradsimple
implicit none

	
if(doradsimple) then

!  A simple predefined radiation (longwave only)

       call rad_simple()
	 
else


! Call full radiation package:
 
    write(*,*)'ERROR: Full radiation not implemented yet.'
    stop
    !call rad_full()	
 
endif

end subroutine


