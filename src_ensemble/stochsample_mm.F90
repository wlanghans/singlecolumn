Program Parallel_Run_Ensemble
  use mpi
  use netcdf
  USE IFPORT
  implicit none

  character(300) :: param, path_out,pathout,path_ncl
  character(500) :: fileout
  character(700) ::command
  character(3) :: char_id
  character(4) :: charnruns

  integer :: my_id, Nproc, ierr, N_local_runs, N_max
  integer, allocatable, dimension(:) :: mpi_id

  integer :: i,j,k
  integer,parameter :: nruns=1000
  integer :: itmp(10)
  integer :: run_start, run_end, status

  logical :: fileexist, filenowexist
 
  real(8) :: t1,t2



  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nproc,ierr)

  allocate(mpi_id(Nproc-1))

  itmp(1) = 0
  do i = 0, Nproc-1
    if(i /= my_id) then
      itmp(1) = itmp(1) + 1
      mpi_id( itmp(1) ) = i
    endif
  enddo

  call getarg(1,param)
  path_out = trim(adjustl(param))
  pathout='"'//trim(path_out)//'"'

  call getarg(2,param)
  path_ncl='"'//trim(adjustl(param))//'"'



  write(char_id,'(i3)') my_id
  write(charnruns,'(i4)') nruns



   fileout='"'//trim(path_out)//"/RMSE"//"_"//trim(adjustl(char_id))//'"'
   write(*,*) 'Writing file ', fileout

!  char_id='"'//trim(char_id)//'"'
!  charnruns='"'//trim(charnruns)//'"'

   command="cd "//trim(path_out)//"; ncl -Q "//trim(path_ncl)//"/getrmse_mm_rmseqtdthvd_stochastic.ncl benchdir="//"'"//trim(path_ncl)//"'"//" fileout="//"'"//trim(fileout)//"'"//" mpiid="//"'"//trim(char_id)//"'"//" nruns="//"'"//trim(charnruns)//"'"

   status = system( trim(command))

   
   fileout=trim(path_out)//"/RMSE"//"_"//trim(adjustl(char_id))

   INQUIRE(FILE=trim(fileout),EXIST=fileexist)
   if (.not.fileexist) then

      t1  = MPI_WTIME() 
 
      do while (.not.fileexist)
         INQUIRE(FILE=trim(fileout),EXIST=fileexist)
         t2 = MPI_WTIME() 
         if (t2-t1.gt.300.) then
            write(*,*) 'ERROR: Waited for 8 min to get postprocessing done but something is wrong. Stopping.'
            write(*,*) 'File ',fileout,' does not exist'
            call MPI_FINALIZE(ierr)
         end if  
      end do

   end if


  call MPI_barrier(MPI_COMM_WORLD,ierr)
  call MPI_Finalize(ierr)

  if(my_id == 0) write(*,*) 'MPI finalized ...'

  
contains

  function char2int(c)
    implicit none

    character(*), intent(in) :: c

    integer :: char2int
    integer :: i
    character(100) :: c2

    c2 = trim(c)   ! Eliminate prefixed blanks
    char2int = 0
    do i = 1,len(trim(c2))
      char2int = 10*char2int + index('0123456789',c2(i:i)) - 1
    end do

  end function char2int


end Program Parallel_Run_Ensemble

