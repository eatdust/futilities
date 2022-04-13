module iopara
  use mpi

  implicit none

  integer, parameter :: dp = kind(1.0_8)
  
  private



  public open_parallel_file, close_parallel_file
  public ascii_write_parallel_vectors
  

contains


  subroutine open_parallel_file(filename,ides)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ides
    integer :: icode

    character(len=MPI_MAX_ERROR_STRING) :: errorMessage
    integer :: lenErrorMessage
    integer :: icodeBis

    call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_RDWR + MPI_MODE_CREATE, &
         MPI_INFO_NULL,ides,icode)

    if (icode .ne. MPI_SUCCESS) then
       call MPI_ERROR_STRING(icode,errorMessage,lenErrorMessage,icodeBis)
       write(*,*)'open_parallel_file:'
       write(*,*)errorMessage(1:lenErrorMessage)
       call MPI_ABORT(MPI_COMM_WORLD,42,icode)
    endif

  end subroutine open_parallel_file



  subroutine ascii_write_parallel_vectors(ides,data)
    implicit none
    integer :: ides
    real(dp), dimension(:), intent(in) :: data

!they should all match!
    integer, parameter :: lencbuffer = 16
    character(len=*), parameter :: cformat = '(ES16.8)'
    character(len=lencbuffer) :: cbuffer

!variable size string
    character(:), allocatable :: cvec
    
    integer :: i, nsize

    integer :: icode
    integer, dimension(MPI_STATUS_SIZE) :: status
    
    nsize = size(data,1)

    cvec = ''
    do i=1,nsize
       write(cbuffer,fmt=cformat)data(i)
       cvec = cvec//cbuffer
    enddo

    if (len(cvec).ne.nsize*lencbuffer) then
       write(*,*)'asciiwrite_parallel_vectors:'
       write(*,*)'len(cvec)= ',len(cvec),nsize*lencbuffer
       stop 'character lengtb error!'
    endif

    call MPI_FILE_WRITE_ORDERED(ides,cvec//new_line(cbuffer),nsize*lencbuffer+1,MPI_CHARACTER &
         ,status,icode)
    
    deallocate(cvec)

  end subroutine ascii_write_parallel_vectors
  

  subroutine close_parallel_file(ides)
    implicit none
    integer :: ides
    integer :: icode
    
    call MPI_FILE_CLOSE(ides,icode)

    if (icode.ne.MPI_SUCCESS) then
       write(*,*)'close_parallel_file: error!'
       call MPI_ABORT(MPI_COMM_WORLD,2,icode)
    endif
    
  end subroutine close_parallel_file
  


end module iopara
