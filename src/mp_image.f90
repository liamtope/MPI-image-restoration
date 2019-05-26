module mp_image
  use mpi
  implicit none

  ! ------------------------------------------------------------------  
  !  mp_image : parallel routines for the image reconstruction
  !  		program. 
  ! 		No MPI function calls required in the main program,
  !		all are contained in this module
  ! ------------------------------------------------------------------  

  contains

  subroutine mpStart(procid, nprocs, communicator)
    implicit none

    ! Parameters
    integer, intent(out)                :: procid        ! Process ID
    integer, intent(out)                :: nprocs        ! Number of processes
    integer, intent(out)                :: communicator  ! Communicator containing processes
    ! Other declerations
    integer                             :: ierr, world = mpi_comm_world   ! Initial communicator
    integer                             :: decompdims = 2   ! 2 by default (ammended later)
    logical                             :: reorder = .false.
    integer, allocatable, dimension(:)  :: dims
    logical, allocatable, dimension(:)  :: periods

    call mpi_init(ierr)
    call mpi_comm_size(world, nprocs, ierr)

    if (is_prime(nprocs)) decompdims = 1
    allocate(dims(decompdims), periods(decompdims))

    ! Creating cartesian communicator
    dims(:) = 0
    ! Use appropriate dimensions
    if (decompdims.gt.1) then
      call mpi_dims_create(nprocs, decompdims, dims, ierr)
    else
      dims(1) = nprocs
    end if
    ! Create the cartesian topology
    periods(:) = .false.
    call mpi_cart_create(world, decompdims, dims, periods, reorder, communicator, ierr)
    ! Find rank
    call mpi_comm_rank(communicator, procid, ierr)
    
  end subroutine


  subroutine mpStop()
    implicit none
  
    ! No arguments - only error variable
    integer :: ierr

    call mpi_finalize(ierr)

  end subroutine


  subroutine mpSubDomainSize(coords, cartdims, sizes)
    implicit none

    ! Arguments
    integer, intent(in), dimension(:)      :: coords       ! Cartesian coordinates of this process
    integer, intent(in), dimension(:)      :: cartdims     ! Dimensions of cartesian grid communicator
    integer, intent(inout), dimension(:)   :: sizes        ! Size of domain
    ! Other declerations
    integer                                  :: ierr, procits, ii, jj, amax, diff
    integer                                  :: ndims_domain, decompdims

    ! Dimension of the domain
    ndims_domain = size(sizes)
    decompdims = size(coords)

    jj = decompdims
    ! Decompose image
    do ii=ndims_domain, ndims_domain-decompdims+1, -1

      ! Iterations per process
      procits = sizes(ii) / cartdims(jj)
      ! Remainder from division
      diff = procits * cartdims(jj) - sizes(ii)
      ! Offset to maximum bound
      amax = logic2int(coords(jj).lt.diff) 
    
      ! Iteration limit
      sizes(ii) = procits + amax
    
      jj = jj - 1

    end do
    
  end subroutine
 

  subroutine mpSwapHalos(arr, communicator)
    implicit none
    ! Arguments
    ! ---------
    real(kind=8),     dimension(:,:)  :: arr            ! Array to update halos
    integer                           :: communicator   ! Cartesian communicator
    ! Other declerations
    ! -------------------
    ! Neighbouring processes
    integer                           :: left, right, up, down
    ! Array size
    integer                           :: imin, imax, jmin, jmax
    ! MPI messaging variables
    integer                           :: dtype = mpi_double_precision
    integer                           :: up_send_request, down_send_request
    integer                           :: left_send_request, right_send_request
    integer                           :: up_recv_request, down_recv_request
    integer                           :: left_recv_request, right_recv_request
    integer                           :: ierr, status(mpi_status_size)
    ! Vector-type variables
    integer                           :: blocklen, stride, count
    integer                           :: vector
    ! Communicator information
    integer                           :: nprocs

    ! Communicator information
    call mpi_comm_size(communicator, nprocs, ierr)

    imin = lbound(arr,1);   imax = ubound(arr,1)
    jmin = lbound(arr,2);   jmax = ubound(arr,2)

    ! If we have been passed a 1D decomposed domain call 1D routine
    if (is_prime(nprocs)) then

      ! Vertical neighbours
      call mpi_cart_shift(communicator, 0, 1, down, up, ierr)
      ! Horizontal neighbours
      left = mpi_proc_null
      right = mpi_proc_null

    else  ! 2D Decomposition
    
      ! Vertical neighbours
      call mpi_cart_shift(communicator, 1, 1, down, up, ierr)
      ! Horizintal neighbours
      call mpi_cart_shift(communicator, 0, 1, left, right, ierr)
    
    end if

    ! Utilise non-blocking synchronous sends
    ! --------------------------------------
    ! Send up
    call mpi_issend(arr(imin, jmax-1), size(arr,1), dtype, up, 0, &
                    communicator, up_send_request, ierr)
    ! Send down
    call mpi_issend(arr(imin, jmin+1), size(arr,1), dtype, down, 0, &
                    communicator, down_send_request, ierr)

    ! Receive bottom halo
    call mpi_recv(arr(imin,jmin), size(arr,1), dtype, down, 0, &
                   communicator, status, ierr)
    ! Recieve top halo
    call mpi_recv(arr(imin,jmax), size(arr,1), dtype, up, 0, &
                   communicator, status, ierr)

    ! Wait for inital sends to finish
    call mpi_wait(up_send_request, status, ierr)
    call mpi_wait(down_send_request, status, ierr)

    ! Construct vectors for left and right sends 
    blocklen = 1
    stride = size(arr,1)
    count = size(arr,2) - 2

    call mpi_type_vector(count, blocklen, stride, dtype, vector, ierr)
    call mpi_type_commit(vector, ierr)

    ! Send right
    call mpi_issend(arr(imax-1,jmin+1), 1, vector, right, 0, &
                    communicator, right_send_request, ierr)
    ! Send left
    call mpi_issend(arr(imin+1,jmin+1), 1, vector, left, 0, &
                    communicator, left_send_request, ierr)

    ! Receive right halo
    call mpi_recv(arr(imax,jmin+1), 1, vector, right, 0, &
                   communicator, status, ierr)
    ! Receive left halo
    call mpi_recv(arr(imin,jmin+1), 1, vector, left, 0, &
                   communicator, status, ierr)

    ! Now wait for sends to go through
    call mpi_wait(right_send_request, status, ierr)
    call mpi_wait(left_send_request, status, ierr)

  end subroutine

  function within_precision(old, new, eps, communicator) result(finished)

    ! Arguments
    real(kind=8), dimension(:,:)               :: old, new      ! Arrays before and after algorithm iteration
    real(kind=8)                               :: eps           ! Degree of precision
    integer                                    :: communicator  ! Communicator
    ! Output
    logical                                    :: finished

    ! Other declerations
    real(kind=8)                               :: maxdiff, maxdiff_sum
    integer                                    :: ierr, dtype = mpi_double_precision

    ! Maximum difference on this process
    maxdiff = maxval(old - new)

     ! Check other processes also 
    call mpi_allreduce(maxdiff, maxdiff_sum, 1, dtype, mpi_sum, communicator, ierr)

    finished = .false.
    if (maxdiff_sum .lt. eps)  finished = .true.

  end function


  function logic2int(logic) result(integ)
  
    logical, intent(in) :: logic
    integer             :: integ

    integ = 0
    if(logic) integ = 1
    
  end function


  function is_prime(n) result(yn)
    
    ! Arguments
    integer :: n 
    logical :: yn
    ! Other declerations
    integer :: d, num_divisors

    yn = .false.
    num_divisors=0

    do d=1,n
      if (mod(n,d)==0) num_divisors = num_divisors + 1
    end do

    if (num_divisors.eq.2) yn = .true.

  end function

end module
