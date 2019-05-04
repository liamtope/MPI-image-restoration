!------------------------------------------------------------------------
!  _ ____    ____                        _                              |
! / |  _ \  |  _ \  ___  _ __ ___   __ _(_)_ __                         |
! | | | | | | | | |/ _ \| '_ ` _ \ / _` | | '_ \                        |
! | | |_| | | |_| | (_) | | | | | | (_| | | | | |                       |
! |_|____/  |____/ \___/|_| |_| |_|\__,_|_|_| |_| _ _   _               |
! |  _ \  ___  ___ ___  _ __ ___  _ __   ___  ___(_) |_(_) ___  _ __    |
! | | | |/ _ \/ __/ _ \| '_ ` _ \| '_ \ / _ \/ __| | __| |/ _ \| '_ \   |
! | |_| |  __/ (_| (_) | | | | | | |_) | (_) \__ \ | |_| | (_) | | | |  |
! |____/ \___|\___\___/|_| |_| |_| .__/ \___/|___/_|\__|_|\___/|_| |_|  |
!                                |_|                                    |
!                                                                       |
!     -->  1D Decomposition of a 2D domain                              |
!                                                                       |
!------------------------------------------------------------------------

module decomp_1d
  use mpi
  use mp_image
  implicit none

  contains

    ! Distribute information to sub-domains
    ! --------------------------------------
    subroutine DistributeData1D(sourcedata, recvdata, procid, communicator, domainsize)
      implicit none

      ! Arguments
      real(kind=8), intent(in),               dimension(:,:)    :: sourcedata     ! Initial Data
      real(kind=8), intent(out), allocatable, dimension(:,:)    :: recvdata       ! Process specific data
      integer,      intent(in)                                  :: procid         ! Process ID
      integer,      intent(in)                                  :: communicator   ! Cartesian communicator
      integer,      intent(in),  optional,     dimension(:)     :: domainsize     ! Shape of initial data
      ! Other declerations
      integer,                   allocatable,  dimension(:)     :: senddisp, sendcounts
      integer                                                   :: ii, nx, ny, inx, iny, mynx, myny
      integer                                                   :: dtype = mpi_double_precision
      integer                                                   :: ierr, root = 0
      integer                                                   :: startx, starty, endx, endy
      integer                                                   :: globalsize(2)
      integer,                                  dimension(1)    :: coords, cartdims, dummycoords
      logical,                                  dimension(1)    :: cartperiods

      if (.not.present(domainsize)) then
        nx = size(sourcedata,1)
        ny = size(sourcedata,2)
      else
        nx = domainsize(1)
        ny = domainsize(2)
      end if

      ! Cartesian domain information
      call mpi_cart_get(communicator, 1, cartdims, cartperiods, coords, ierr)

      ! Initialise information for scatterv
      allocate(senddisp(cartdims(1)), sendcounts(cartdims(1)))
      senddisp(1) = 0

      ! Find all sizes and start points of sub-domains
      do ii=0,cartdims(1)-1
      
        ! Get size of sub-domains
        dummycoords(1) = ii
        globalsize(1) = nx;  globalsize(2) = ny

        call mpSubDomainSize(dummycoords, cartdims, globalsize)
        inx = globalsize(1);  iny = globalsize(2)
        
        ! Bounds for this process
        if (ii.eq.coords(1)) then
          mynx = inx
          myny = iny
        end if
        
        ! Size of this send
        sendcounts(ii+1) = inx * iny
        ! Displacement of this send
        if (ii>0) senddisp(ii+1) = senddisp(ii) + sendcounts(ii) 
      end do
      
      allocate(recvdata(mynx, myny))
      
      call mpi_scatterv(sourcedata, sendcounts, senddisp, dtype, recvdata, size(recvdata), &
                        dtype, root, communicator, ierr)

      deallocate(senddisp, sendcounts)

    end subroutine


    ! Group back all data from sub-domains
    ! ------------------------------------
    subroutine GroupData1d(procdata, groupdata, procid, communicator)
      implicit none
      ! Arguments
      real(kind=8), intent(in),        dimension(:,:)    :: procdata     ! Data specific to a given thread
      real(kind=8), intent(out),       dimension(:,:)    :: groupdata    ! Grouped data array
      integer                                            :: procid       ! Process ID in communicator
      integer,      intent(in)                           :: communicator ! Cartesian communicator
      ! Other declerations 
      integer,            allocatable, dimension(:)      :: recvdisp, recvcounts
      integer                                            :: dtype=mpi_double_precision 
      integer                                            :: ierr, nx, ny, ii, inx, iny
      integer                                            :: root = 0
      integer                                            :: dummycoords(1), coords(1), cartdims(1)
      integer                                            :: globalsize(2), decompdims
      logical                                            :: cartperiods(1)

      ! Size of global domain
      nx = size(groupdata,1)
      ny = size(groupdata,2)

      ! Cartesian domain dimensions
      call mpi_cart_get(communicator, 1, cartdims, cartperiods, coords, ierr)
      dummycoords(1) = coords(1)

      ! Initialise information for gatherv
      allocate(recvcounts(cartdims(1)), recvdisp(cartdims(1)))
      recvdisp(1) = 0

      ! Find all sizes and start points of sub-domains
      do ii=0,cartdims(1)-1
        
        dummycoords(1) = ii
        ! Get bounds of sub-domains 
        globalsize(1) = nx;  globalsize(2) = ny
        call mpSubDomainSize(dummycoords, cartdims, globalsize)
        inx = globalsize(1);  iny = globalsize(2)
    
        ! Size of receive
        recvcounts(ii+1) = inx * iny
        ! Displacement of receive
        if (ii>0) recvdisp(ii+1) = recvdisp(ii) + recvcounts(ii)
    
      end do

      call mpi_gatherv(procdata, size(procdata), dtype, groupdata, recvcounts, recvdisp, &
                      dtype, root, communicator, ierr)

      deallocate(recvcounts, recvdisp)

    end subroutine

end module