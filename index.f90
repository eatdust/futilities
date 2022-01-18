module index
  use precision, only : fdp
  implicit none

  private

  integer, dimension(:), allocatable, save :: shuffleIndex

  logical, parameter :: useSimpleMapping = .false.

  integer, parameter :: shuffleseed = 1

  public create_shuffle_index, index_mapping, free_shuffle_index

contains


  subroutine free_shuffle_index()
    implicit none

    if (allocated(shuffleIndex)) then
       deallocate(shuffleIndex)
    else
       stop 'free_schuffle_index: not allocated!'
    endif

  end subroutine free_shuffle_index

  
  subroutine create_shuffle_index(nsize)
    implicit none
    integer, intent(in) :: nsize
    integer, dimension(:), allocatable :: inseeds

    integer :: n, i,j,shufj
    real(fdp) :: harv

    if (allocated(shuffleIndex)) then
       stop 'create_shuffle_index: already allocated!'
    endif

    !initialize random generator
    call random_seed(size=n)
    allocate(inseeds(n))
    inseeds = shuffleseed* (/(i - 1, i = 1, n)/)
    call random_seed(put=inseeds)
    deallocate(inseeds)

    write(*,*)
    write(*,*)'random generator seeded with: ', shuffleseed   
    write(*,*)

    allocate(shuffleIndex(1:nsize))

    forall (i=1:nsize)
       shuffleIndex(i) = i
    end forall

    do i=1,nsize
       call random_number(harv)
!in [i,nsize]
       j = i + int((nsize-i+1)*harv)
       shufj = shuffleIndex(j)
       shuffleIndex(j) = shuffleIndex(i)
       shuffleIndex(i) = shufj
    enddo


  end subroutine create_shuffle_index




  function index_mapping(i)
    implicit none
    integer, intent(in) :: i
    integer :: index_mapping

    if (useSimpleMapping) then
       index_mapping = i
    else
       index_mapping = shuffleIndex(i)
    endif
  end function index_mapping





end module index
