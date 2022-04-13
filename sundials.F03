module sundials
!currently support only serial cvode in linear and non-linear solver,
!to be expanded...
  use, intrinsic :: iso_c_binding
  use precision, only : fdp

  use fcvode_mod

  use fsundials_nvector_mod
  use fnvector_serial_mod

  use fsunlinsol_dense_mod
  use fsundials_linearsolver_mod
  
  use fsundials_nonlinearsolver_mod
  use fsunnonlinsol_fixedpoint_mod

  use fsunmatrix_dense_mod
  use fsundials_matrix_mod
  
  implicit none

  private

  include 'sundials.h'

  
  type fsun
     type(c_ptr) :: mem
     type(N_Vector), pointer :: vecY
     type(SUNMatrix), pointer :: matA
     type(SUNLinearSolver), pointer :: ls 
     type(SUNNonlinearSolver), pointer :: nls
  end type fsun


  type tcvode
     integer :: ndim
     procedure(cvoderhs), pointer, nopass :: ptr_rhs => null()
     procedure(cvodejac), pointer, nopass :: ptr_jac => null()
  end type tcvode

  
  type(tcvode) :: cvodedata
!$omp threadprivate(cvodedata)
  
  integer(c_long), parameter :: mxsteps = 100000
   
  public easycvode, stepscvode_nls, stepscvode_ls
  
contains


  function easycvode(n,cvoderhs,tini,yini,tend,tol)
    implicit none
    include 'sundials.h'
    integer, intent(in) :: n
    real(fdp), intent(in) :: tini
    real(fdp), intent(in) :: tend
    real(fdp), intent(in), optional :: tol
    real(fdp), dimension(n), intent(in) :: yini
    real(fdp), dimension(n) :: easycvode

    real(fdp), dimension(n,1) :: output
       
    output = stepscvode_nls(n,cvoderhs,tini,yini,(/tend/),tol)

    easycvode = output(1:n,1) 
    
  end function easycvode


  

  function stepscvode_nls(n,cvoderhs,tini,yini,tend,tol)
    implicit none

    include 'sundials.h'

    integer, intent(in) :: n
    real(fdp), intent(in) :: tini
    real(fdp), dimension(:), intent(in) :: tend
    real(fdp), dimension(n), intent(in) :: yini
    real(fdp), dimension(n,size(tend,1)) :: stepscvode_nls
    real(fdp), intent(in) , optional :: tol
    
    integer(c_long) :: neq
    real(c_double) :: tstart, tout
    real(c_double), dimension(1) :: tcur
    real(c_double), dimension(n) :: yvec

    type(fsun) :: sun

    integer(c_int) :: ierr

    
    
    real(c_double) :: rtol
    real(c_double), parameter :: atol = 100._c_double * epsilon(1._c_double)
  
    integer :: i,nsteps

    if (present(tol)) then
       rtol = real(tol,c_double)
    else
       rtol = 1d-10
    endif
       
    
    neq = int(n,c_long)
    tstart = real(tini,c_double)
    yvec = real(yini,c_double)

    nsteps = size(tend,1)
    
!initialize module procedure pointer
    cvodedata%ndim = neq
    cvodedata%ptr_rhs => cvoderhs
    
!create sundial vect (yvec is target)
    sun%vecY => FN_VMake_Serial(neq,yvec)


!allocate CVode solver (ADAMS method or BDF for stiff)
    sun%mem = FCVodeCreate(CV_ADAMS)
    
!initialization
    ierr = FCVodeInit(sun%mem, c_funloc(cvode_rhs), tstart, sun%vecY)

    if (ierr.ne.0) then
       stop ' stepscvode_nls: initialization error in FCVodeInit!'
    endif

    ierr = FCVodeSStolerances(sun%mem, rtol, atol)

    if (ierr.ne.0) then
       stop ' stepscvode_nls: error in setting tolerances!'
    endif

!default is 500    
    ierr = FCVodeSetMaxNumSteps(sun%mem,mxsteps)
    
!create a non-linear solver    
    sun%nls => FSUNNonlinSol_FixedPoint(sun%vecY, 0)


!attach nonlinear solver object to CVode
    ierr = FCVodeSetNonlinearSolver(sun%mem, sun%nls)
    
    if (ierr.ne.0) then
       write(*,*)'ierr = ',ierr
       stop ' stepscvode_nls: cannot attach non-linear solver'
    endif

!solving, tcur is overwritten by tout, possibly allowing for multiple
!steps
    tcur = tstart

    do i=1,nsteps
       tout = real(tend(i),c_double)
!       if (tout.le.tcur(1)) then
!         write(*,*) 'stepscvode_nls: tend vector is not steadily increasing!'
!       endif
       ierr = FCVode(sun%mem, tout, sun%vecY, tcur, CV_NORMAL)

       if (ierr.ne.0) then
          stop 'stepscvode_nls: error in FCVode, integration failed!'
       endif

       stepscvode_nls(1:n,i) = yvec(1:n)
    enddo
    
!freeing        
    cvodedata%ndim = 0
    cvodedata%ptr_rhs => null()
    
    call FCVodeFree(sun%mem)
    ierr = FSUNNonLinSolFree(sun%nls)
    call FN_VDestroy(sun%vecY)
    
  end function stepscvode_nls
  




  function stepscvode_ls(n,cvoderhs,cvodejac,tini,yini,tend)
    implicit none

    include 'sundials.h'

    integer, intent(in) :: n
    real(fdp), intent(in) :: tini
    real(fdp), dimension(:), intent(in) :: tend
    real(fdp), dimension(n), intent(in) :: yini
    real(fdp), dimension(n,size(tend,1)) :: stepscvode_ls
    
    integer(c_long) :: neq
    real(c_double) :: tstart, tout
    real(c_double), dimension(1) :: tcur
    real(c_double), dimension(n) :: yvec

    type(fsun) :: sun

    integer(c_int) :: ierr
    
    
    real(c_double), parameter :: rtol = 1d-10
    real(c_double), parameter :: atol = 100._c_double * epsilon(1._c_double)
  
    integer :: i,nsteps
    
    
    neq = int(n,c_long)
    tstart = real(tini,c_double)
    yvec = real(yini,c_double)

    nsteps = size(tend,1)
    
!initialize module procedure pointer
    cvodedata%ndim = neq
    cvodedata%ptr_rhs => cvoderhs
    cvodedata%ptr_jac => cvodejac
    
!create sundial vect (yvec is target)
    sun%vecY => FN_VMake_Serial(neq,yvec)

!create a dense matrix
    sun%matA => FSUNDenseMatrix(neq, neq)

!create a linear solver
    sun%ls => FSUNDenseLinearSolver(sun%vecY, sun%matA)
    
!allocate CVode solver (ADAMS method or BDF for stiff)
    sun%mem = FCVodeCreate(CV_ADAMS)
    
!initialization
    ierr = FCVodeInit(sun%mem, c_funloc(cvode_rhs), tstart, sun%vecY)

    if (ierr.ne.0) then
       stop ' stepscvode_ls: initialization error in FCVodeInit!'
    endif

    ierr = FCVodeSStolerances(sun%mem, rtol, atol)

    if (ierr.ne.0) then
       stop ' stepscvode_ls: error in setting tolerances!'
    endif

!default is 500    
    ierr = FCVodeSetMaxNumSteps(sun%mem,mxsteps)
    
!attach the liner solver
    ierr = FCVodeSetLinearSolver(sun%mem, sun%ls, sun%matA);

    if (ierr.ne.0) then
       write(*,*)'ierr = ',ierr
       stop ' stepscvode_ls: cannot attach linear solver'
    endif

!set the Jacobian routine
    ierr = FCVodeSetJacFn(sun%mem, c_funloc(cvode_jac))
  

!solving, tcur is overwritten by tout, possibly allowing for multiple
!steps
    tcur = tstart

    do i=1,nsteps
       tout = real(tend(i),c_double)
       if (tout.le.tcur(1)) then
          stop 'stepscvode_ls: tend vector is not steadily increasing!'
       endif
       ierr = FCVode(sun%mem, tout, sun%vecY, tcur, CV_NORMAL)

       if (ierr.ne.0) then
          stop 'stepscvode_ls: error in FCVode, integration failed!'
       endif

       stepscvode_ls(1:n,i) = yvec(1:n)
    enddo
    
!freeing        
    cvodedata%ndim = 0
    cvodedata%ptr_rhs => null()
    
    call FCVodeFree(sun%mem)
    ierr = FSUNLinSolFree(sun%ls)
    call FSUNMatDestroy(sun%matA)
    call FN_VDestroy(sun%vecY)
    
  end function stepscvode_ls
  





  
  function cvode_rhs(t,y,f,userdata) bind(C,name='RhsFn')
    implicit none
    integer(c_int) :: cvode_rhs

    real(c_double), value :: t
    type(N_Vector) :: y, f
    type(c_ptr), value :: userdata

    real(c_double), dimension(:), pointer :: yvec
    real(c_double), dimension(:), pointer :: fvec
        
    yvec => FN_VGetArrayPointer(y)
    fvec => FN_VGetArrayPointer(f)
    
    fvec = cvodedata%ptr_rhs(cvodedata%ndim,t,yvec)

!success return value    
    cvode_rhs = 0    

  end function cvode_rhs


  
  function cvode_jac(t,y,f,dfody,userdata,tmp1,tmp2,tmp3) bind(C,name='JacFn')
    implicit none
    integer(c_int) :: cvode_jac

    real(c_double), value :: t
    type(N_Vector) :: y, f
    type(SUNMatrix) :: dfody
    type(c_ptr), value :: userdata
    type(N_Vector) :: tmp1, tmp2, tmp3
    
    real(c_double), dimension(:), pointer :: yvec
    real(c_double), dimension(:), pointer :: fvec
    real(c_double), dimension(:), pointer :: jmat
        
    yvec => FN_VGetArrayPointer(y)
    fvec => FN_VGetArrayPointer(f)
    jmat =>  FSUNDenseMatrix_Data(dfody)

!to be checked and debugged, should be column-major ordering (fortran standard)
    jmat = cvodedata%ptr_jac(cvodedata%ndim,t,yvec,fvec)

    cvode_jac = 0    

  end function cvode_jac
  
  
end module sundials
