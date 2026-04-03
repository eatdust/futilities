module sundials
!currently support only serial cvode in linear and non-linear solver,
!to be expanded...
  use, intrinsic :: iso_c_binding
  use precision, only : fdp

  use fcvodes_mod
  use fsundials_core_mod

  use fnvector_serial_mod
#ifdef SUNOPENMP
  use fnvector_openmp_mod
#endif  
  
  use fsunlinsol_dense_mod
  use fsunnonlinsol_fixedpoint_mod
  use fsunnonlinsol_newton_mod

  use fsunmatrix_dense_mod

  
  implicit none

  private

  include 'sundials.h'

  
  type fsun
     type(c_ptr) :: mem
     type(c_ptr) :: ctx
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

  
  integer, save :: inst = 1
!$omp threadprivate(inst)
  
  integer, parameter :: ninstmax = 3
  type(tcvode), dimension(ninstmax) :: cvodedata
!$omp threadprivate(cvodedata)
  
  integer(c_long), parameter :: mxsteps = 10000000

#ifdef SUNOPENMP  
  integer, parameter :: nthreads = 2
#endif
  
  public easycvode, easycvquad
  public cvode_set_instance, cvode_get_instance
  public stepscvode_nls, stepscvode_ls
  
contains


  subroutine cvode_set_instance(setinst)
    implicit none
    integer, intent(in) :: setinst
    
    if (setinst.gt.ninstmax) stop 'set_instance: setinst > ninstmax!'

    inst = setinst

  end subroutine cvode_set_instance


  function cvode_get_instance()
    implicit none
    integer :: cvode_get_instance

    cvode_get_instance = inst

  end function  cvode_get_instance
  
  
!non-linear solver
  function easycvode(n,cvoderhs,tini,yini,tend,tol,tola)
    implicit none
    include 'sundials.h'
    integer, intent(in) :: n
    real(fdp), intent(in) :: tini
    real(fdp), intent(in) :: tend
    real(fdp), intent(in), optional :: tol,tola
    real(fdp), dimension(n), intent(in) :: yini
    real(fdp), dimension(n) :: easycvode

    real(fdp), dimension(n,1) :: output
       
    output = stepscvode_nls(n,cvoderhs,tini,yini,(/tend/),tol,tola)
    
    easycvode = output(1:n,1) 
    
  end function easycvode


!linear solver with vanishing jacobian  
  function easycvquad(n,cvoderhs,tini,yini,tend,tol,tola)
    implicit none
    include 'sundials.h'
    integer, intent(in) :: n
    real(fdp), intent(in) :: tini
    real(fdp), intent(in) :: tend
    real(fdp), intent(in), optional :: tol,tola
    real(fdp), dimension(n), intent(in) :: yini
    real(fdp), dimension(n) :: easycvquad

    real(fdp), dimension(n,1) :: output

    output = stepscvode_ls(n,cvoderhs,nulljac,tini,yini,(/tend/),tol,tola)

    easycvquad = output(1:n,1) 

  end function easycvquad


  

  function stepscvode_nls(n,cvoderhs,tini,yini,tend,tolr,tola)
#ifdef MPI
    use mpi
#endif    
    implicit none

    include 'sundials.h'

    integer, intent(in) :: n
    real(fdp), intent(in) :: tini
    real(fdp), dimension(:), intent(in) :: tend
    real(fdp), dimension(n), intent(in) :: yini
    real(fdp), dimension(n,size(tend,1)) :: stepscvode_nls
    real(fdp), intent(in) , optional :: tolr, tola
    
    integer(c_long) :: neq
    real(c_double) :: tstart, tout
    real(c_double), dimension(1) :: tcur
    real(c_double), dimension(n) :: yvec

    type(fsun) :: sun

    integer(c_int) :: ierr
        
    real(c_double) :: rtol, atol
    
    integer :: i,nsteps

    if (present(tolr)) then
       rtol = real(tolr,c_double)
    else
       rtol = 100000._c_double * epsilon(1._c_double)
    endif

    if (present(tola)) then
       atol = real(tola,c_double)
    else
       atol = 100._c_double * epsilon(1._c_double)
    endif
    
    neq = int(n,c_long)
    tstart = real(tini,c_double)
    yvec = real(yini,c_double)

    nsteps = size(tend,1)

!context creation (sundials >= 7.0)
#ifndef MPI
    ierr = FSUNContext_Create(SUN_COMM_NULL,sun%ctx)
#else
    ierr = FSUNContext_Create(MPI_COMM_NULL,sun%ctx)
#endif    

    if (ierr.ne.0) stop 'error in sunctx creation!'
    
!initialize module procedure pointer
    cvodedata(inst)%ndim = neq
    cvodedata(inst)%ptr_rhs => cvoderhs
    
!create sundial vect (yvec is target)
#ifndef SUNOPENMP
    sun%vecY => FN_VMake_Serial(neq, yvec, sun%ctx)
#else
    sun%vecY => FN_VMake_OpenMP(neq, yvec, nthreads, sun%ctx)
#endif
    
    

!allocate CVode solver (ADAMS method or BDF for stiff)
    sun%mem = FCVodeCreate(CV_ADAMS ,sun%ctx)

!default is 12 for adams, 5 for bdf
!    ierr = FCVodeSetMaxOrd(sun%mem,8)
    
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
    sun%nls => FSUNNonlinSol_FixedPoint(sun%vecY, 0, sun%ctx)


!attach nonlinear solver object to CVode
    ierr = FCVodeSetNonlinearSolver(sun%mem, sun%nls)
    
    if (ierr.ne.0) then
       write(*,*)'ierr = ',ierr
       stop ' stepscvode_nls: cannot attach non-linear solver'
    endif

!default is 3
!    ierr = FCVodeSetMaxNonlinIters(sun%mem,3)
    
!null pointer to userdata as we don't use it
    ierr = FCVodeSetUserData(sun%mem, c_null_ptr)
    
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
    cvodedata(inst)%ndim = 0
    cvodedata(inst)%ptr_rhs => null()
    
    call FCVodeFree(sun%mem)
    sun%mem = c_null_ptr
    ierr = FSUNNonLinSolFree(sun%nls)
    sun%nls => null()
    call FN_VDestroy(sun%vecY)
    sun%vecY => null()
    
    
!delete context    
    ierr = FSUNContext_Free(sun%ctx)
    sun%ctx = c_null_ptr
    
  end function stepscvode_nls
  




  function stepscvode_ls(n,cvoderhs,cvodejac,tini,yini,tend,tolr,tola)
#ifdef MPI
    use mpi
#endinf    
    implicit none

    include 'sundials.h'

    integer, intent(in) :: n
    real(fdp), intent(in) :: tini
    real(fdp), dimension(:), intent(in) :: tend
    real(fdp), dimension(n), intent(in) :: yini
    real(fdp), intent(in) , optional :: tolr, tola
    real(fdp), dimension(n,size(tend,1)) :: stepscvode_ls
    
    integer(c_long) :: neq
    real(c_double) :: tstart, tout
    real(c_double), dimension(1) :: tcur
    real(c_double), dimension(n) :: yvec
    
    type(fsun) :: sun

    integer(c_int) :: ierr
    
    
    real(c_double) :: rtol
    real(c_double) :: atol
  
    integer :: i,nsteps, nthreads
        
    if (present(tolr)) then
       rtol = real(tolr,c_double)
    else
       rtol = 100000._c_double * epsilon(1._c_double)
    endif

    if (present(tola)) then
       atol = real(tola,c_double)
    else
       atol = 100._c_double * epsilon(1._c_double)
    endif
    
    neq = int(n,c_long)
    tstart = real(tini,c_double)
    yvec = real(yini,c_double)

    nsteps = size(tend,1)
    
!initialize module procedure pointer
    cvodedata(inst)%ndim = neq
    cvodedata(inst)%ptr_rhs => cvoderhs
    cvodedata(inst)%ptr_jac => cvodejac

!context creation (sundials >= 7.0)
#ifndef MPI
    ierr = FSUNContext_Create(SUN_COMM_NULL,sun%ctx)
#else
    ierr = FSUNContext_Create(MPI_COMM_NULL,sun%ctx)
#endif    
    if (ierr.ne.0) stop 'error in sunctx creation!'
    
!create sundial vect (yvec is target)
#ifndef SUNOPENMP
    sun%vecY => FN_VMake_Serial(neq, yvec, sun%ctx)
#else
    sun%vecY => FN_VMake_OpenMP(neq, yvec, nthreads, sun%ctx)
#endif
    
!create a dense matrix
    sun%matA => FSUNDenseMatrix(neq, neq, sun%ctx)

!create a linear solver
    sun%ls =>  FSUNLinSol_Dense(sun%vecY, sun%matA, sun%ctx)
    
!allocate CVode solver (ADAMS method or BDF for stiff)
    sun%mem = FCVodeCreate(CV_ADAMS, sun%ctx)

!default is 12 for adams, 5 for bdf
!    ierr = FCVodeSetMaxOrd(sun%mem,8)
    
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

!null pointer to userdata as we don't use it
    ierr = FCVodeSetUserData(sun%mem, c_null_ptr)
    
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
    cvodedata(inst)%ndim = 0
    cvodedata(inst)%ptr_rhs => null()
    
    call FCVodeFree(sun%mem)
    sun%mem = c_null_ptr
    ierr = FSUNLinSolFree(sun%ls)
    sun%ls => null()
    call FSUNMatDestroy(sun%matA)
    sun%matA => null()
    call FN_VDestroy(sun%vecY)
    sun%vecY => null()
    
!delete context    
    ierr = FSUNContext_Free(sun%ctx)
    sun%ctx = c_null_ptr
    
  end function stepscvode_ls
  





  
  function cvode_rhs(t,y,f,userdata) bind(C,name='CVRhsFn')
    implicit none
    integer(c_int) :: cvode_rhs

    real(c_double), value :: t
    type(N_Vector) :: y, f
    type(c_ptr) :: userdata
    
    real(c_double), dimension(:), pointer :: yvec
    real(c_double), dimension(:), pointer :: fvec

    yvec => FN_VGetArrayPointer(y)
    fvec => FN_VGetArrayPointer(f)
    
    fvec = cvodedata(inst)%ptr_rhs(cvodedata(inst)%ndim,t,yvec)

!success return value    
    cvode_rhs = 0    

  end function cvode_rhs

  
  function cvode_jac(t,y,f,dfody,userdata,tmp1,tmp2,tmp3) bind(C,name='CVJacFn')
    implicit none
    integer(c_int) :: cvode_jac

    real(c_double), value :: t
    type(N_Vector) :: y, f
    type(SUNMatrix) :: dfody
    type(c_ptr), pointer :: userdata
    type(N_Vector) :: tmp1, tmp2, tmp3
    
    real(c_double), dimension(:), pointer :: yvec
    real(c_double), dimension(:), pointer :: fvec
    real(c_double), dimension(:), pointer :: jmat

    yvec => FN_VGetArrayPointer(y)
    fvec => FN_VGetArrayPointer(f)
    jmat =>  FSUNDenseMatrix_Data(dfody)

!to be checked and debugged, should be column-major ordering (fortran standard)
    jmat = cvodedata(inst)%ptr_jac(cvodedata(inst)%ndim,t,yvec,fvec)

    cvode_jac = 0    

  end function cvode_jac



!for test and easy quadrature  
  function nulljac(n,t,y,f)
    use precision, only : fdp
    implicit none
    integer, intent(in) :: n
    real(fdp), intent(in) :: t
    real(fdp), dimension(n), intent(in) :: y,f
    real(fdp), dimension(n*n) :: nulljac

    nulljac = 0._fdp

  end function nulljac


  
end module sundials
