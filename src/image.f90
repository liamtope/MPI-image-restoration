program image
  use pgmio
  use mp_image
  use decomp_2d

  implicit none

  ! Declerations
  integer :: nx, ny, niter, nxloc, nyloc, imagesize(2)
  integer :: nn, ii, jj, rank, nprocs, comm
  real(kind=8), allocatable, dimension(:,:) :: newimg, oldimg, edges, localedges
  character(len=50) :: filename, ch_niter, newfile
  integer :: root = 0     ! Master process
  
  ! Initialise message passing
  call mpStart(rank, nprocs, comm)
  
  ! ---------------------------------------
  !        EXTRACTING DATA FROM FILE      |
  ! ---------------------------------------
 
  ! File we want to use  (CL arg)
  call getarg(1, filename)

  ! Number of algorithm iterations (CL arg)
  call getarg(2, ch_niter)
  read(ch_niter, *) niter
  
  ! Image dimensions
  call pgmsize(trim(filename), nx, ny)

  ! Read the contents of the image file to edge array
  allocate(edges(1:nx, 1:ny))
  if (rank.eq.root) call pgmread(trim(filename), edges)


  ! --------------------------------------
  !         DECOMPOSITION OF IMAGE        |
  ! --------------------------------------

  !write(*,*) 'Starting to distribute data on rank ', rank
  ! Decompose and distribute the image
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
  do nn=1,niter

    if ((mod(nn,100) .eq. 0).and.(rank.eq.0)) write(*, "('Iteration ', i4)") nn

    do jj=1,nyloc
      do ii=1,nxloc
      
        ! Update the image
        newimg(ii,jj) = .25 * (oldimg(ii-1,jj) + oldimg(ii+1,jj) + oldimg(ii,jj-1) &
                        + oldimg(ii,jj+1) - localedges(ii,jj))    
      end do
    end do

    ! Reset oldimg to newimg before next iteration
    oldimg(1:nxloc, 1:nyloc) = newimg(1:nxloc, 1:nyloc)  
    
    ! Update halos
    call mpSwapHalos(oldimg, comm)

  end do  

  ! Clean up used data
  deallocate(oldimg, localedges)
  
  ! --------------------------------------
  !      GROUP ALL DATA BACK TOGETHER    |
  ! --------------------------------------

  ! Combine all process' data
  call GroupData2D(newimg(1:nxloc, 1:nyloc), edges, rank, comm) 
 
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
