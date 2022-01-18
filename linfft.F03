module linfft
  use precision, only : fdp
  use, intrinsic :: iso_c_binding
  
  implicit none

  include 'fftw3.f03'
  
  private

  type lingrid
     integer :: n
     real(C_DOUBLE), dimension(:), pointer :: val => null()
     complex(C_DOUBLE_COMPLEX), dimension(:), pointer :: cxval => null()
     type(C_PTR) :: sto
     type(C_PTR), pointer :: plan => null()
  end type lingrid


  real(fdp), parameter :: pi = 3.1415926535897932384626433832795_fdp


  public lingrid, pi
  public check_lingrid, allocate_lingrid, free_lingrid, reset_lingrid
  public fft_1D_plan_forw, fft_1D_plan_back
  public fft_1D_forw, fft_1D_back, free_plan
  public index_conjugate, unit_wavenumber, unit_position

contains

  function check_lingrid(ag)
    implicit none
    logical :: check_lingrid
    type(lingrid), intent(in) :: ag

    check_lingrid = associated(ag%val)

  end function check_lingrid


  subroutine allocate_lingrid(ngrid,ag)
    implicit none
    integer, intent(in) :: ngrid
    type(lingrid), intent(inout) :: ag
    
    integer(C_SIZE_T) :: ndata

    if (check_lingrid(ag)) stop 'allocate_lingrid: already associated!'

!    allocate(ag%val(2*(ngrid/2+1)))
    
    ndata = int((ngrid/2+1), C_SIZE_T)

    ag%n=ngrid
    ag%sto = fftw_alloc_complex(ndata)
    call c_f_pointer(ag%sto, ag%val, [2*(ngrid/2+1)] )
    call c_f_pointer(ag%sto, ag%cxval, [ngrid/2+1] )

  end subroutine allocate_lingrid



  subroutine reset_lingrid(ag)
    implicit none
    integer :: nsize
    type(lingrid), intent(inout) :: ag

    integer :: i

    if (.not.check_lingrid(ag)) stop 'reset_lingrid: not associated!'

    nsize = size(ag%val,1)

    forall (i=1:nsize)
       ag%val(i) = 0._C_DOUBLE
    end forall

  end subroutine reset_lingrid



  subroutine free_lingrid(ag)
    implicit none
    type(lingrid), intent(inout) :: ag

    if (.not.check_lingrid(ag)) stop 'free_lingrid: not associated!'

!    deallocate(ag%val)

    ag%val => null()
    ag%cxval => null()
    ag%n=0
    call fftw_free(ag%sto)

    
  end subroutine free_lingrid



  subroutine fft_1D_plan_forw(ag)
    implicit none
    type(lingrid), intent(inout) :: ag
    integer(C_INT) :: n

    if (.not.check_lingrid(ag)) stop 'ftt_1D_plan_forw: not associated!'
    
    n = ag%n

    if (associated(ag%plan)) stop 'fft_1D_plan_forw: plan already exists!'

    allocate(ag%plan)

#ifdef PATIENT
    write(*,*)'fft_1D_plan_forw: PATIENT plan...'
    ag%plan = fftw_plan_dft_r2c_1d(n,ag%val,ag%cxval,FFTW_PATIENT)
#else
    write(*,*)'fft_1D_plan_forw: MEASURE plan...'
    ag%plan = fftw_plan_dft_r2c_1d(n,ag%val,ag%cxval,FFTW_MEASURE)
#endif

  end subroutine fft_1D_plan_forw


  subroutine fft_1D_forw(ag)
    implicit none
    type(lingrid), intent(inout) :: ag
    integer(C_INT) :: n

    if (.not.check_lingrid(ag)) stop 'ftt_1D_forw: not associated!'

    n = ag%n
    
    if (.not.associated(ag%plan)) then
       write(*,*)'fft_1D_forw: plan not found, ESTIMATE one...'
       ag%plan = fftw_plan_dft_r2c_1d(n,ag%val,ag%cxval,FFTW_ESTIMATE)
    endif
    
    call fftw_execute_dft_r2c(ag%plan,ag%val,ag%cxval)

  end subroutine fft_1D_forw



  subroutine fft_1D_plan_back(ag)
    implicit none
    type(lingrid), intent(inout) :: ag
    integer(C_INT) :: n

    if (.not.check_lingrid(ag)) stop 'ftt_1D_plan_back: not associated!'
    
    n = ag%n

    if (associated(ag%plan)) stop 'fft_1D_plan_back: plan already exists!'

    allocate(ag%plan)

#ifdef PATIENT
    write(*,*)'fft_1D_plan_back: PATIENT plan...'
    ag%plan = fftw_plan_dft_c2r_1d(n,ag%cxval,ag%val,FFTW_PATIENT)
#else
    write(*,*)'fft_1D_plan_back: MEASURE plan...'
    ag%plan = fftw_plan_dft_c2r_1d(n,ag%cxval,ag%val,FFTW_MEASURE)
#endif

  end subroutine fft_1D_plan_back
  


  subroutine fft_1D_back(ag)
    implicit none
    type(lingrid), intent(inout) :: ag
    integer(C_INT) :: n

    if (.not.check_lingrid(ag)) stop 'ftt_1D_back: not associated!'

    n = ag%n
    
    if (.not.associated(ag%plan)) then
       write(*,*)'fft_1D_back: plan not found, ESTIMATE one...'
       write(*,*)
       ag%plan = fftw_plan_dft_c2r_1d(n,ag%cxval,ag%val,FFTW_ESTIMATE)
    endif
    
    call fftw_execute_dft_c2r(ag%plan,ag%cxval,ag%val)

  end subroutine fft_1D_back
  


  subroutine free_plan(ag)
    implicit none
    type(lingrid), intent(inout) :: ag

    if (.not.associated(ag%plan)) stop 'free_plan: not associated!'

    call fftw_destroy_plan(ag%plan)

    deallocate(ag%plan)

    ag%plan => null()

  end subroutine free_plan


  function index_conjugate(ngrid,q)
    implicit none
    integer, intent(in) :: ngrid,q
    integer :: index_conjugate

    if (q.eq.1) then
!this is the index of k=0
       index_conjugate = 1
    else
!and of -k
       index_conjugate = ngrid + 2 - q
    endif

  end function index_conjugate


  function unit_wavenumber(ngrid,q)
    implicit none
    integer, intent(in) :: ngrid,q
    real(fdp) :: unit_wavenumber
    
    if (q.le.1+ngrid/2) then
       unit_wavenumber = 2.*pi*real(q-1,fdp)
    else
       unit_wavenumber = 2.*pi*real(q-1-ngrid,fdp)
    endif

  end function unit_wavenumber


  function unit_position(ngrid,q)
    implicit none
    integer, intent(in) :: ngrid,q
    real(fdp) :: unit_position
    
    if (q.le.1+ngrid/2) then
       unit_position = real(q-1,fdp)/real(ngrid,fdp)
    else
       unit_position = real(q-1-ngrid,fdp)/real(ngrid,fdp)
    endif

  end function unit_position


  
end module linfft
