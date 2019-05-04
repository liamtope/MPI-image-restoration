!-----------------------------------------------------------------------
!  ____  ____    ____                        _                         |
! |___ \|  _ \  |  _ \  ___  _ __ ___   __ _(_)_ __                    |
!   __) | | | | | | | |/ _ \| '_ ` _ \ / _` | | '_ \                   |
!  / __/| |_| | | |_| | (_) | | | | | | (_| | | | | |                  |
! |_____|____/  |____/ \___/|_| |_| |_|\__,_|_|_| |_|   _              |
! |  _ \  ___  ___ ___  _ __ ___  _ __   ___  ___(_) |_(_) ___  _ __   |
! | | | |/ _ \/ __/ _ \| '_ ` _ \| '_ \ / _ \/ __| | __| |/ _ \| '_ \  |
! | |_| |  __/ (_| (_) | | | | | | |_) | (_) \__ \ | |_| | (_) | | | | |
! |____/ \___|\___\___/|_| |_| |_| .__/ \___/|___/_|\__|_|\___/|_| |_| |
!                                |_|                                   |
!                                                                      |
!    -->   2D decomposition of a 2D domain                             |
!                                                                      |
! ----------------------------------------------------------------------

module decomp_2d
  use mpi 
  use decomp_1d
  use mp_image

  implicit none

  contains

  ! ----------------------------------------------
  ! ____  _     _        _ _           _       
  ! |  _ \(_)___| |_ _ __(_) |__  _   _| |_ ___ 
  ! | | | | / __| __| '__| | '_ \| | | | __/ _ \
  ! | |_| | \__ \ |_| |  | | |_) | |_| | ||  __/
  ! |____/|_|___/\__|_|  |_|_.__/ \__,_|\__\___|
  ! 
  ! ---------------------------------------------                                              
  subroutine DistributeData2D(sourcedata, recvdata, procid, communicator, domainsize)
    implicit none
    ! Arguments
    ! ---------
    real(kind=8), intent(in),  allocatable, dimension(:,:)  :: sourcedata    ! Initial Data
    real(kind=8), intent(out), allocatable, dimension(:,:)  :: recvdata      ! Process specific data
    integer,      intent(in)                                :: procid        ! Process ID
    integer,      intent(in)                                :: communicator  ! Cartesian communicator
    integer,      intent(in),               dimension(2)    :: domainsize    ! Golbal domain size
    ! Other declerations
    ! ------------------
    ! Information about/from global communicator
    integer                                                 :: globalcoords(2), globaldims(2)
    integer                                                 :: rootsize(2), dummycoords(2), globalsize(2)
    logical                                                 :: globalperiods(2)
    integer                                                 :: row, col
    ! Information about/from row/col communicator
    integer                                                 :: rowcomm_id, colcomm_id
    integer                                                 :: rowcomm_size, colcomm_size
    integer                                                 :: rcomm, ccomm, rowcomm, colcomm
    integer                                                 :: rcomm_coords(1), rcomm_dims(1)
    logical                                                 :: rcomm_periods(1), ccomm_periods(1)
    integer                                                 :: ccomm_coords(1), ccomm_dims(1)
    ! Data to be sent/recieved
    real(kind=8),              allocatable, dimension(:,:)  :: dummycolData, dummyrecv, colData
    ! MPI subroutine arguments
    integer                                                 :: ierr, decompdims = 2
    ! Sub-domain size 
    integer                                                 :: nrow, ncol, nx, ny
 
    globalsize = domainsize
    nx = globalsize(1);  ny = globalsize(2)

    ! ------------------------------------------------------------------
    ! Create communicators for processes in the same cartesian rows/cols
    ! ------------------------------------------------------------------
    call mpi_cart_get(communicator, 2, globaldims, globalperiods, globalcoords, ierr)

    row = globalcoords(1);  col = globalcoords(2)

    rcomm_periods(:) = .false.
    ! 'Same row' communicator
    call mpi_comm_split(communicator, row, procid, rcomm, ierr)
    rcomm_dims(:) = globaldims(2)
    call mpi_cart_create(rcomm, 1, rcomm_dims, rcomm_periods, .false., rowcomm, ierr)

    ccomm_periods(:) = .false.
    ! 'Same column' communicator
    call mpi_comm_split(communicator, col, procid, ccomm, ierr)
    ccomm_dims = globaldims(1)
    call mpi_cart_create(ccomm, 1, ccomm_dims, ccomm_periods, .false., colcomm, ierr)

    ! ----------------------
    ! Distributing the data --> 2D Scatter performed with two 1D scatters 
    ! ----------------------
    ! First send columns to those in same row as root process
    if (row.eq.0) then

      call mpi_comm_rank(rowcomm, rowcomm_id, ierr)
      call DistributeData1D(sourcedata, colData, rowcomm_id, rowcomm, globalsize)

      ! Take the transpose so the next send can be contiguous
      allocate(dummycolData(size(colData,2), size(colData,1)))
      dummycolData = transpose(colData)
      deallocate(colData)

    end if

    globalsize = domainsize
    ! All processes need to know size of source for next distribution
    dummycoords(1) = 0;  dummycoords(2) = col
    call mpSubDomainSize(dummycoords, globaldims, globalsize)
    
    globalsize(1) = globalsize(2)
    globalsize(2) = domainsize(1)

    call mpi_comm_rank(colcomm, colcomm_id, ierr)
    call DistributeData1D(dummycolData, dummyrecv, colcomm_id, colcomm, globalsize)
      
    if(row.eq.0) deallocate(dummycolData)

    ! Transpose back to normal shape
    allocate(recvdata(size(dummyrecv,2), size(dummyrecv,1)))
    recvdata = transpose(dummyrecv)
      
    deallocate(dummyrecv)

    
  
  end subroutine

  
  ! ------------------------------------------------
  !   ____                                        
  !  |  _ \ ___        __ _ _ __ ___  _   _ _ __  
  !  | |_) / _ \_____ / _` | '__/ _ \| | | | '_ \ 
  !  |  _ <  __/_____| (_| | | | (_) | |_| | |_) |
  !  |_| \_\___|      \__, |_|  \___/ \__,_| .__/ 
  !                   |___/                |_|    
  !
  ! ------------------------------------------------  
  subroutine GroupData2D(procdata, groupdata, procid, communicator)
    implicit none
    ! Arguments
    ! ---------
    real(kind=8), intent(in),   dimension(:,:)    :: procdata     ! Data specific to a given thread
    real(kind=8), intent(out),  dimension(:,:)    :: groupdata    ! Grouped data array
    integer,      intent(in)                      :: procid       ! Process ID
    integer,      intent(in)                      :: communicator ! Cartesian communicator
    ! Other declerations
    ! ------------------
    ! Data to be sent/received
    real(kind=8),  allocatable, dimension(:,:)    :: nextsend, dummyrecv, dummysend
    ! MPI subroutine args
    integer                                       :: ierr, decompdims = 2
    ! Information about/from global domain
    integer                                       :: row, col, n, nx, ny
    integer                                       :: globalcoords(2), globaldims(2)
    integer                                       :: rootsize(2), dummycoords(2), globalsize(2)
    logical                                       :: globalperiods(2)
    ! Row/Column communicator information
    integer                                       :: rowcomm_id, colcomm_id
    integer                                       :: rowcomm_size, colcomm_size
    integer                                       :: rowcomm, colcomm, rcomm, ccomm
    integer                                       :: rcomm_dims(1),  ccomm_dims(1)
    logical                                       :: rcomm_periods(1), ccomm_periods(1)
    ! Sub-domain information
    integer,                    dimension(1)      :: nblocks
    logical,                    dimension(1)      :: periods

    nx = size(groupdata,1);   ny = size(groupdata,2)

    ! ------------------------------------------------------------------
    ! Create communicators for processes in the same cartesian rows/cols
    ! ------------------------------------------------------------------
    call mpi_cart_get(communicator, 2, globaldims, globalperiods, globalcoords, ierr)

    if (globalcoords(2).lt.0) then 
      ! 1D Distribution if negative number of columns
      call GroupData1D(procdata, groupdata, procid, communicator)
    else

      row = globalcoords(1);  col = globalcoords(2)

      rcomm_periods(:) = .false.
      ! 'Same row' communicator
      call mpi_comm_split(communicator, row, procid, rcomm, ierr)
      rcomm_dims(:) = globaldims(2)
      call mpi_cart_create(rcomm, 1, rcomm_dims, rcomm_periods, .false., rowcomm, ierr)

      ccomm_periods(:) = .false.
      ! 'Same column' communicator
      call mpi_comm_split(communicator, col, procid, ccomm, ierr)
      ccomm_dims = globaldims(1)
      call mpi_cart_create(ccomm, 1, ccomm_dims, ccomm_periods, .false., colcomm, ierr)


      ! ----------------------
      !  Re-groupong the data --> 2D gather performed with two 1D gathers
      ! ----------------------

      ! Transpose the process data for a contiguous send
      allocate(dummysend(size(procdata,2), size(procdata,1)))
      dummysend = transpose(procdata)

      ! Alocate the receive of all the rows
      allocate(dummyrecv(size(procdata,2), nx))

      ! Get all the rows data to the zeroth row
      call mpi_comm_rank(colcomm, colcomm_id, ierr)
      call GroupData1D(dummysend, dummyrecv, colcomm_id, colcomm)

      deallocate(dummysend)

      if (row.eq.0) then
        ! Transpose what we received to be back in correct shape
        allocate(nextsend(nx, size(procdata,2)))
        nextsend = transpose(dummyrecv)

        ! Send via the 'same row' communicator
        call mpi_comm_rank(rowcomm, rowcomm_id, ierr)
        call GroupData1D(nextsend, groupdata, rowcomm_id, rowcomm)

        deallocate(nextsend, dummyrecv)

      end if
    end if

  end subroutine

end module