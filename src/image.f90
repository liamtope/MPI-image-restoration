program image
  use pgmio
  use mp_image
  use decomp_2d

  implicit none

  ! Declerations
  integer                                     :: nx, ny, nxloc, nyloc, imagesize(2)
  integer                                     :: nn, ii, jj, rank, nprocs, comm
  logical                                     :: converged
  real(kind=8)                                :: old_ave, new_ave
  real(kind=8), allocatable, dimension(:,:)   :: newimg, oldimg, edges, localedges
  character(len=50)                           :: filename, newfile
  integer                                     :: root = 0             ! Master process
  real(kind=8)                                :: eps =  1.0d-3        ! Degree of precision
  real(kind=8)                                :: start, stopp         ! Timing variables

  ! Initialise message passing
  call mpStart(rank, nprocs, comm)
  
  ! ---------------------------------------
  !        EXTRACTING DATA FROM FILE      |
  ! ---------------------------------------
 
  ! File we want to use  (CL arg)
  call getarg(1, filename)
  
  ! Image dimensions
  call pgmsize(trim(filename), nx, ny)

  ! Read the contents of the image file to edge array
  allocate(edges(1:nx, 1:ny))
  if (rank.eq.root) call pgmread(trim(filename), edges)

  ! Start timer
  start = mpi_wtime()

  ! --------------------------------------
  !         DECOMPOSITION OF IMAGE        |
  ! --------------------------------------

  imagesize(1) = nx;   imagesize(2) = ny
  call DistributeData2D(edges, localedges, rank, comm, imagesize)

  nxloc = size(localedges,1)
  nyloc = size(localedges,2)

  ! Allocate the image and edge arrays with halos
  allocate(newimg(0:nxloc+1, 0:nyloc+1))
  allocate(oldimg(0:nxloc+1, 0:nyloc+1))

  ! --------------------------------------
  !       IMAGE RESTORATION ALGORITHM    |
  ! --------------------------------------

  ! Set Initial guess to all white
  oldimg(:,:) = 255

  ! Iterate algorithm to converge
  nn = 1
  do

    ! Condition for printing which iteration we are on
    !if ((mod(nn,500) .eq. 0).and.(rank.eq.0)) write(*, "('Iteration ', i5)") nn 

    ! Image restoration algorithm
    do jj=1,nyloc
      do ii=1,nxloc

        ! Update the image
        newimg(ii,jj) = .25 * (oldimg(ii-1,jj) + oldimg(ii+1,jj) + oldimg(ii,jj-1) &
                        + oldimg(ii,jj+1) - localedges(ii,jj))    
      end do
    end do

    ! Is newimg within specified degree of precision?
    converged = within_precision(oldimg(1:nxloc, 1:nyloc), newimg(1:nxloc, 1:nyloc), eps, comm)
    if (converged) exit

    ! Otherwise, reset oldimg to newimg before next iteration
    oldimg(1:nxloc, 1:nyloc) = newimg(1:nxloc, 1:nyloc)  
    
    ! Update halos
    call mpSwapHalos(oldimg, comm)

    ! Increment iterable
    nn = nn + 1

  end do
  
  ! Algorithm done  - clean up used data
  deallocate(oldimg, localedges)
  
  ! --------------------------------------
  !      GROUP ALL DATA BACK TOGETHER    |
  ! --------------------------------------

  ! Combine all process' data
  call GroupData2D(newimg(1:nxloc, 1:nyloc), edges, rank, comm) 
   
  ! Stop timer
  stopp = mpi_wtime()

  if(rank.eq.root) write(*, "(a20, ',  ', i3, ',  ', f14.9)") trim(filename), nprocs, stopp-start

  ! Finished with newimg now
  deallocate(newimg)
  
  ! Write to file
  if (rank.eq.0) then
    write(newfile, '("image", i3, "x", i3, ".pgm")') nx, ny
    call pgmwrite(trim(newfile), edges)
  end if

  ! Final clean up
  deallocate(edges)

  ! Stop message passing
  call mpStop()

end program
