module cuba
  implicit none

  private    

  integer, parameter :: dp = kind(1._8)

  integer, parameter :: isp = kind(1_4)
  integer, parameter :: idp = kind(1_8)
  
  logical, parameter :: display = .false.

  
!parameter common to all integrators

!default minimum and maximum numbers of iteration
  integer(idp), parameter :: minevaldef = 1024
  integer(idp), parameter :: maxevaldef = 1048576

  
!statefile name (empty is no statefile)
  character(len=*), parameter :: statefile=""

!spin=-1 is automatic allocation, and most importantly, deallocation
!of spinning cores
  integer(kind(1_8)), save :: spin = -1
  
!absolute accuracy (pushed to machine precision)  
  real(dp), parameter :: epsabs = epsilon(1._dp)
 

  interface simplevegas
     module procedure sp_simplevegas, dp_simplevegas
  end interface simplevegas

  interface simplesuave
     module procedure sp_simplesuave, dp_simplesuave
  end interface simplesuave

  interface simplecuhre
     module procedure sp_simplecuhre, dp_simplecuhre
  end interface simplecuhre
  
  interface simpledivonne
     module procedure sp_simpledivonne, dp_simpledivonne
  end interface simpledivonne
  
  public cuba_set_spin, cuba_close_spin, cuba_set_cores

  public simplevegas, simplesuave, simplecuhre, simpledivonne

  
contains

!also accessible with environment variables CUBACORES, CUBACORESMAX  
  subroutine cuba_set_cores(ncores,pcores)
    implicit none
    integer, intent(in) :: ncores,pcores

    call cubacores(ncores,pcores)
    
  end subroutine cuba_set_cores
  

  

  subroutine cuba_set_spin(inispin)
    integer, intent(in) :: inispin

    spin = inispin

  end subroutine cuba_set_spin

  subroutine cuba_close_spin()

    call cubawait(spin)
    
  end subroutine cuba_close_spin

  
  
!settings bit per bit
  subroutine set_seedflags(iseed,iflag)
    implicit none
    integer, intent(out) :: iseed
    integer, intent(out) :: iflag

    iseed = 0
    iflag = 0

!control verbosity 0 to 3 on 2 bits
    iflag = ibclr(iflag,0)
    iflag = ibclr(iflag,1)

!keep all samples 0, or only the last 1
    iflag = ibclr(iflag,2)

!Vegas and Suave: Apply additional smoothing 0 or no smoothing 1
!(do not use smoothing on sharp edges)
    iflag = ibclr(iflag,3)

!Delete the state file at the end of the integration 0; keep it 1.
    iflag = ibclr(iflag,4)

!Vegas only, take the integrator' state from the state file if one is
!present 0; or reset the integrator' state while keeping the grid if
!one is present 1 (with bit4=1, you can reuse a grid)
    iflag = ibclr(iflag,5)    

!Bits 8-31 are for the random number generator (level=4 for highest
!possible luxury when iseed is non-zero):
    if (iseed.ne.0) then
       iflag = ibclr(iflag,8)
       iflag = ibclr(iflag,9)
       iflag = ibclr(iflag,10)
    endif

    if (display) then
       write(*,*)'iflag= iseed= ',iflag, iseed
    endif
    
  end subroutine set_seedflags



  
  subroutine diagnostics(nregions,neval,ifail,integral,error,prob,name,tolasked)
    implicit none
    integer, intent(in) :: nregions,neval,ifail
    real(dp), dimension(:), intent(in) :: integral, error, prob
    real(dp), intent(in), optional :: tolasked
    character(len=*), intent(in) :: name

    real(dp) :: toldiag

    if (present(tolasked)) then
       toldiag = epsilon(1._dp)/tolasked
    else
       toldiag = 0._dp
    endif
    
    if (ifail.ne.0) then
       write(*,*)'cuba: accuracy goal failed!'
       write(*,*)'integrator is: ',name
       write(*,*)'neval= nregions=',neval,nregions
       write(*,*)'ifail= error= ',ifail, error
    endif

    if (any(prob.eq.1._dp).and.(maxval(abs(integral)).gt.toldiag)) then
       write(*,*)'cuba: error most likely underestimated!'
       write(*,*)'integrator is: ',name
       write(*,*)'integral= error= ',integral, error
       write(*,*)'prob= ',prob
    endif
    
    if (display) then
       write(*,*)'cuba: neval= ',neval
       write(*,*)'integrator is: ',name
       write(*,*)'integral= ',integral(:)
       write(*,*)'error   = ',error(:)
       write(*,*)'prob    = ',prob(:)
       write(*,*)
    end if
  end subroutine diagnostics




  function sp_simplevegas(ndim,ncomp,integrand,tol,nvsize,maxiter)
    implicit none
    integer(isp), intent(in) :: ndim,ncomp
    real(dp), intent(in) :: tol
    integer(isp), intent(in), optional :: nvsize, maxiter
    real(dp), dimension(ncomp) :: sp_simplevegas
    
    include 'cuba.h'
    
    integer, parameter :: nregions = 1
    integer :: neval, ifail
    real(dp), dimension(ncomp) :: integral, error, prob

    integer, parameter :: userdata = 0
    integer :: nvec, iflag, iseed, maxeval, mineval

!specific to vegas
    integer, parameter :: nstart = 1000
    integer, parameter :: nincrease = 500
    integer, parameter :: nbatch = 1000
    integer, parameter :: gridno = 0


    call set_seedflags(iseed,iflag)

    if (present(nvsize)) then
       nvec = nvsize
    else
       nvec = 1
    endif
   
    if (present(maxiter)) then
       maxeval = maxiter
    else
       maxeval = maxevaldef
    endif

    mineval = minevaldef
    
    call vegas(ndim,ncomp,integrand,userdata,nvec,tol,epsabs,iflag &
         ,iseed,mineval,maxeval,nstart,nincrease,nbatch,gridno &
         ,statefile,spin,neval,ifail,integral,error,prob)


    call diagnostics(nregions,neval,ifail,integral,error,prob,'vegas')        

    sp_simplevegas = integral

    
  end function sp_simplevegas


  function dp_simplevegas(ndim,ncomp,integrand,tol,nvsize,maxiter)
    implicit none
    integer(idp), intent(in) :: ndim,ncomp
    real(dp), intent(in) :: tol
    integer(idp), intent(in), optional :: nvsize, maxiter
    real(dp), dimension(ncomp) :: dp_simplevegas
    
    include 'cuba.h'
    
    integer, parameter :: nregions = 1
    integer :: ifail
    real(dp), dimension(ncomp) :: integral, error, prob

    integer, parameter :: userdata = 0
    integer :: iflag, iseed

    integer(idp) :: neval, nvec, maxeval, mineval
    
!specific to vegas
    integer(idp), parameter :: nstart = 1000
    integer(idp), parameter :: nincrease = 500
    integer(idp), parameter :: nbatch = 10000
    integer, parameter :: gridno = 0


    call set_seedflags(iseed,iflag)

    if (present(nvsize)) then
       nvec = nvsize
    else
       nvec = 1
    endif
   
    if (present(maxiter)) then
       maxeval = maxiter
    else
       maxeval = maxevaldef
    endif

    mineval = minevaldef
    
    call llvegas(ndim,ncomp,integrand,userdata,nvec,tol,epsabs,iflag &
         ,iseed,mineval,maxeval,nstart,nincrease,nbatch,gridno &
         ,statefile,spin,neval,ifail,integral,error,prob)


    call diagnostics(nregions,int(neval,isp),ifail,integral,error,prob,'vegas')        

    dp_simplevegas = integral

    
  end function dp_simplevegas


  


  function sp_simplesuave(ndim,ncomp,integrand,tol,nvsize,maxiter,flat)
    implicit none    
    integer(isp), intent(in) :: ndim,ncomp
    real(dp), intent(in) :: tol
    integer(isp), intent(in), optional :: nvsize, maxiter, flat
    real(dp), dimension(ncomp) :: sp_simplesuave

    include 'cuba.h'
    
    integer :: nregions, neval, ifail
    real(dp), dimension(ncomp) :: integral, error, prob

    integer, parameter :: userdata = 0
    integer :: nvec, iflag, iseed, maxeval, mineval
    
!suave specific    
!number of new integrand evaluations in each subdivision    
    integer, parameter :: nmin = 128
    integer, parameter :: nnew = 16384
!should be large for flat integrand, small for volatile
    real(dp):: flatness

    
    call set_seedflags(iseed,iflag)

    if (present(nvsize)) then
       nvec = nvsize
    else
       nvec = 1
    endif
   
    if (present(maxiter)) then
       maxeval = maxiter
    else
       maxeval = maxevaldef
    endif

    mineval = minevaldef

    if (present(flat)) then
       flatness = flat
    else
       flatness = 10
    endif
    
    call suave(ndim,ncomp,integrand,userdata,nvec,tol,epsabs,iflag &
         ,iseed,mineval,maxeval,nnew,nmin,flatness,statefile,spin &
         ,nregions,neval,ifail,integral,error,prob)

    call diagnostics(nregions,neval,ifail,integral,error,prob,'suave')        

    sp_simplesuave = integral
        
  end function sp_simplesuave


  function dp_simplesuave(ndim,ncomp,integrand,tol,nvsize,maxiter,flat)
    implicit none    
    integer(idp), intent(in) :: ndim,ncomp
    real(dp), intent(in) :: tol
    integer(idp), intent(in), optional :: nvsize, maxiter, flat
    real(dp), dimension(ncomp) :: dp_simplesuave

    include 'cuba.h'
    
    integer :: nregions, ifail
    real(dp), dimension(ncomp) :: integral, error, prob

    integer, parameter :: userdata = 0
    integer :: iflag, iseed

    integer(idp) :: neval, nvec, maxeval, mineval
    
!suave specific    
!number of new integrand evaluations in each subdivision    
    integer(idp), parameter :: nmin = 128
     integer(idp), parameter :: nnew = 16384
!should be large for flat integrand, small for volatile
    real(dp):: flatness

    
    call set_seedflags(iseed,iflag)

    if (present(nvsize)) then
       nvec = nvsize
    else
       nvec = 1
    endif
   
    if (present(maxiter)) then
       maxeval = maxiter
    else
       maxeval = maxevaldef
    endif

    mineval = minevaldef

    if (present(flat)) then
       flatness = real(flat,dp)
    else
       flatness = 10._dp
    endif
    
    call llsuave(int(ndim,isp),int(ncomp,isp),integrand,userdata,nvec,tol,epsabs,iflag &
         ,iseed,mineval,maxeval,nnew,nmin,flatness,statefile,spin &
         ,nregions,neval,ifail,integral,error,prob)

    call diagnostics(nregions,int(neval,isp),ifail,integral,error,prob,'suave')        

    dp_simplesuave = integral
        
  end function dp_simplesuave
  


  function sp_simpledivonne(ndim,ncomp,integrand,tol,nvsize,maxiter)
    implicit none        
    integer(isp), intent(in) :: ndim,ncomp
    real(dp), intent(in) :: tol
    integer(isp), intent(in), optional :: nvsize, maxiter
    real(dp), dimension(ncomp) :: sp_simpledivonne

    include 'cuba.h'
    
    integer :: nregions, neval, ifail
    real(dp), dimension(ncomp) :: integral, error, prob

    integer, parameter :: userdata = 0
    integer :: nvec, iflag, iseed, maxeval, mineval

!divonne specific (cubature rule of degree key12 :7,9,11,13)
    integer, parameter :: key1 = 47
    integer, parameter :: key2 = 1
    integer, parameter :: key3 = 1
    integer, parameter :: maxpass = 5

    real(dp), parameter :: border = 0
    real(dp), parameter :: maxchisq = 10._dp
    real(dp), parameter :: mindeviation = 0.25_dp

    integer, parameter :: ngiven = 0
    real(dp), dimension(ndim,ngiven) :: xgiven
    integer, parameter :: nextra = 0
    integer, pointer :: peakfinder 

    integer :: ldxgiven
    ldxgiven = ndim
    
    call set_seedflags(iseed,iflag)

    if (present(nvsize)) then
       nvec = nvsize
    else
       nvec = 1
    endif

    if (present(maxiter)) then
       maxeval = maxiter
    else
       maxeval = maxevaldef
    endif

    mineval = minevaldef
    
    call divonne(ndim,ncomp,integrand,userdata,nvec,tol,epsabs,iflag,iseed &
         ,mineval,maxeval,key1,key2,key3,maxpass,border,maxchisq,mindeviation &
         ,ngiven,ldxgiven,xgiven,nextra,peakfinder &
         ,statefile,spin,nregions,neval,ifail &
         ,integral,error,prob)
    
    call diagnostics(nregions,neval,ifail,integral,error,prob,'divonne')

    sp_simpledivonne = integral
    
  end function sp_simpledivonne
  

  function dp_simpledivonne(ndim,ncomp,integrand,tol,nvsize,maxiter)
    implicit none        
    integer(idp), intent(in) :: ndim,ncomp
    real(dp), intent(in) :: tol
    integer(idp), intent(in), optional :: nvsize, maxiter
    real(dp), dimension(ncomp) :: dp_simpledivonne

    include 'cuba.h'
    
    integer :: nregions, ifail
    real(dp), dimension(ncomp) :: integral, error, prob

    integer, parameter :: userdata = 0
    integer :: iflag, iseed

    integer(idp) :: neval, nvec, maxeval, mineval

    
!divonne specific (cubature rule of degree key12 :7,9,11,13)
    integer, parameter :: key1 = 47
    integer, parameter :: key2 = 1
    integer, parameter :: key3 = 1
    integer, parameter :: maxpass = 5

    
    real(dp), parameter :: border = 0
    real(dp), parameter :: maxchisq = 10._dp
    real(dp), parameter :: mindeviation = 0.25_dp

    integer(idp), parameter :: ngiven = 0
    real(dp), dimension(ndim,ngiven) :: xgiven
    integer(idp), parameter :: nextra = 0
    integer, pointer :: peakfinder 

    integer :: ldxgiven
    ldxgiven = ndim
    
    call set_seedflags(iseed,iflag)

    if (present(nvsize)) then
       nvec = nvsize
    else
       nvec = 1
    endif

    if (present(maxiter)) then
       maxeval = maxiter
    else
       maxeval = maxevaldef
    endif

    mineval = minevaldef
    
    call lldivonne(int(ndim,isp),int(ncomp,isp),integrand,userdata,nvec,tol,epsabs,iflag,iseed &
         ,mineval,maxeval,key1,key2,key3,maxpass,border,maxchisq,mindeviation &
         ,ngiven,ldxgiven,xgiven,nextra,peakfinder &
         ,statefile,spin,nregions,neval,ifail &
         ,integral,error,prob)
    
    call diagnostics(nregions,int(neval,isp),ifail,integral,error,prob,'divonne')

    dp_simpledivonne = integral
    
  end function dp_simpledivonne
  


  function sp_simplecuhre(ndim,ncomp,integrand,tol,nvsize,maxiter)
    implicit none       
    integer(isp), intent(in) :: ndim,ncomp
    real(dp), intent(in) :: tol
    integer(isp), intent(in), optional :: nvsize, maxiter
    real(dp), dimension(ncomp) :: sp_simplecuhre
    
    include 'cuba.h'
    
    integer :: nregions, neval, ifail
    real(dp), dimension(ncomp) :: integral, error, prob

    integer, parameter :: userdata = 0
    integer :: nvec, iflag, iseed, maxeval, mineval

!cuhre specific (cubature rule of degree key, 7,9,11,13, other value
!triggers default)
!    integer, parameter :: key = 9
    integer, parameter :: key = 13

    
    call set_seedflags(iseed,iflag)

    if (present(nvsize)) then
       nvec = nvsize
    else
       nvec = 1
    endif

    if (present(maxiter)) then
       maxeval = maxiter
    else
       maxeval = maxevaldef
    endif

    mineval = minevaldef
    
    call cuhre(ndim,ncomp,integrand,userdata,nvec,tol,epsabs,iflag &
         ,mineval,maxeval,key,statefile,spin,nregions,neval,ifail &
         ,integral,error,prob)
    
    call diagnostics(nregions,neval,ifail,integral,error,prob,'cuhre',tol)

    sp_simplecuhre = integral
    
  end function sp_simplecuhre
    

  function dp_simplecuhre(ndim,ncomp,integrand,tol,nvsize,maxiter)
    implicit none       
    integer(idp), intent(in) :: ndim,ncomp
    real(dp), intent(in) :: tol
    integer(idp), intent(in), optional :: nvsize, maxiter
    real(dp), dimension(ncomp) :: dp_simplecuhre
    
    include 'cuba.h'
    
    integer :: nregions, ifail
    real(dp), dimension(ncomp) :: integral, error, prob

    integer, parameter :: userdata = 0
    integer :: iflag, iseed

    integer(idp) :: nvec, neval, maxeval, mineval
    
!cuhre specific (cubature rule of degree key, 7,9,11,13, other value
!triggers default)
!    integer, parameter :: key = 9
    integer, parameter :: key = 13

    
    call set_seedflags(iseed,iflag)

    if (present(nvsize)) then
       nvec = nvsize
    else
       nvec = 1
    endif

    if (present(maxiter)) then
       maxeval = maxiter
    else
       maxeval = maxevaldef
    endif

    mineval = minevaldef
    
    call llcuhre(int(ndim,isp),int(ncomp,isp),integrand,userdata,nvec,tol,epsabs,iflag &
         ,mineval,maxeval,key,statefile,spin,nregions,neval,ifail &
         ,integral,error,prob)
    
    call diagnostics(nregions,int(neval,isp),ifail,integral,error,prob,'cuhre',tol)

    dp_simplecuhre = integral
    
  end function dp_simplecuhre
    
  

end module cuba
