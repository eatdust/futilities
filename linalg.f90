module linalg
  implicit none

  private

  integer, parameter :: sp=kind(1.0_4)
  integer, parameter :: dp=kind(1.0_8)


  interface sqrt
     module procedure dp_square_root
  end interface sqrt

  interface rschur
     module procedure dp_real_schur_form
  end interface rschur

  interface rsylvester
     module procedure dp_real_sylvester
  end interface rsylvester

  interface svd
     module procedure dp_svd
  end interface svd

  interface matdiv
     module procedure dp_matdiv
  end interface matdiv

  interface inverse
     module procedure dp_pseudo_inverse
  end interface inverse
  
  interface linsyssym
     module procedure dp_diagpiv_linsys_symmetric
  end interface linsyssym


 interface eigensym
    module procedure dp_divconq_eigen_symmetric
    module procedure dp_iram_eigen_symmetric
    module procedure dp_iram_eigen_band_symmetric
  end interface eigensym

  interface divconq_eigensym
     module procedure dp_divconq_eigen_symmetric
  end interface divconq_eigensym

  interface iram_eigensym
     module procedure dp_iram_eigen_symmetric
     module procedure dp_iram_eigen_band_symmetric
  end interface

 

  public iram_eigensym, divconq_eigensym

  public matdiv, inverse, eigensym
  public sqrt, svd, rschur, rsylvester



contains


  function dp_matdiv(a,b)
    implicit none
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(:,:) :: b
    real(dp), dimension(size(a,1),size(b,1)) :: dp_matdiv

    if (size(a,2).ne.size(b,1)) stop 'dp_matdiv: size mismatch'

    dp_matdiv = matmul(a,dp_pseudo_inverse(b))
    
  end function dp_matdiv



  function dp_pseudo_inverse(a)
    implicit none
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(size(a,2),size(a,1)) :: dp_pseudo_inverse

    real(dp), dimension(size(a,1),size(a,1)) :: u
    real(dp), dimension(min(size(a,1),size(a,2))) :: sdiag
    real(dp), dimension(size(a,2),size(a,2)) :: vt

    real(dp), dimension(size(a,2),size(a,1)) :: spinv
    real(dp), parameter :: zero = 0._dp

    integer :: n,m,i

    m = size(a,1)
    n = size(a,2)

    call dp_svd(a,u,sdiag,vt)

    spinv = 0._dp

    do i=1,min(n,m)
       if (sdiag(i).gt.zero) spinv(i,i)=1._dp/sdiag(i)
    enddo

    dp_pseudo_inverse = matmul(transpose(vt),matmul(spinv, transpose(u)))
    

  end function dp_pseudo_inverse



  function dp_square_root(a)
    implicit none
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(size(a,1),size(a,1)) :: dp_square_root
   
    real(dp), dimension(size(a,1),size(a,1)) :: t, z, r

    call dp_real_schur_form(a,z,t)

    call uptriang_square_root(t,r)
    
    dp_square_root = matmul(z,matmul(r,transpose(z)))
   
  end function dp_square_root



!computes for an N-by-N real nonsymmetric matrix A, the eigenvalues,
!the real Schur form T, and, optionally, the matrix of Schur
!vectors Z
  subroutine dp_real_schur_form(a,z,t)
!m = z.t.z^T
    implicit none
   
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(:,:), intent(out) :: z,t
    
!force computation of schur vectors
    character, parameter  :: JOBVS = 'V'
!do not sort eigenvalue
    character, parameter :: SORT = 'N'
    integer :: SDIM


    integer :: INFO, LDA, LDVS, LWORK, N
    logical, dimension(:), allocatable :: BWORK

    real(dp), dimension(:,:), allocatable :: atemp
    real(dp), dimension(:,:), allocatable :: VS
    real(dp), dimension(:), allocatable :: WI, WR, WORK

    N = size(a,1)
    LDA = N
    LDVS = N
    LWORK = max(1,3*N)

    allocate(atemp(LDA,N))
    allocate(WR(N),WI(N))
    allocate(VS(LDVS,N))

    allocate(WORK(LWORK))
    allocate(BWORK(N))

    atemp = a

    call dgees(JOBVS, SORT, eigensort , N, atemp, LDA,  SDIM,  WR,  WI,  VS,  LDVS &
         ,WORK, LWORK, BWORK, INFO )

    if (INFO.ne.0) stop 'dp_schur_form: error in dgees!'

    t = atemp
    z = VS

    deallocate(atemp, WR, WI, VS)
    deallocate(WORK, BWORK)

  end subroutine dp_real_schur_form

  logical function eigensort(a,b)
    implicit none   
    real(dp) :: a, b

    eigensort = .true.

  end function eigensort
  

!Get a square root for a up-triangular matrix
  subroutine uptriang_square_root(t,r,blocksize)
    implicit none
    real(dp), dimension(:,:), intent(in) :: t
    real(dp), dimension(size(t,1),size(t,1)), intent(out) :: r

    integer, intent(in), optional :: blocksize
    integer, parameter :: bsdefault = 64

    real(dp), dimension(size(t,1)) :: tdiag

    real(dp), dimension(:,:), allocatable :: s,x
    real(dp) :: ss, denom, scale

    integer :: bs
    integer :: n,i,j,k,l
    integer :: nblocks, ns
    integer :: nlarge, nsmall, blarge, bsmall   

    integer, parameter :: nmaxpairs = 8
    integer, dimension(nmaxpairs,2) :: startStopPairs

    integer :: nstart, ncount, nsize
    integer :: istart, istop, jstart, jstop, lstart, lstop
     
    if (present(blocksize)) then
       bs = blocksize
    else
       bs= bsdefault
    end if


    n = size(t,1)

    do j=1,n
       do i=j+1,n
          if (t(i,j).ne.0._dp) stop 't is not upper triangular!'
       enddo
    enddo


    forall (i=1:n)
       tdiag(i) = t(i,i)
    end forall

    r = 0._dp

    if (any(tdiag.lt.0._dp)) stop 'uptriang_square_root: eigenvalue < 0'

    forall (i=1:n)
       r(i,i) = sqrt(tdiag(i))
    end forall

!thanks to numpy
    nblocks = max(1,n/bs)
   
    bsmall = n/nblocks
    blarge = bsmall + 1

    nlarge = modulo(n,nblocks)
    nsmall = nblocks - nlarge

!    print *,'ns',nblocks,nlarge,blarge,nsmall,bsmall

    if (nsmall*bsmall + nlarge*blarge.ne.n) then
       stop 'uptriang_square_root: block screwed!'
    endif

!Define the index range covered by each block.
    nstart = 0
    k=0
    
    do i=1,nsmall
       k=k+1
       startStopPairs(k,:) = (/nstart,nstart+bsmall/)
       nstart = nstart+bsmall
    enddo

    do i=1,nlarge
       k=k+1
       startStopPairs(k,:) = (/nstart,nstart+blarge/)
       nstart = nstart+blarge
    enddo

    startstopPairs = startstopPairs + 1

!    print *,'start',startStopPairs(1:k,1)
!    print *,'stopp',startStopPairs(1:k,2)

!Within-block interactions.

    do l=1,k
       lstart = startStopPairs(l,1)
       lstop = startStopPairs(l,2)
   
       do j=lstart,lstop-1
          do i=j-1,lstart,-1
             ss = 0
             if ((j-i)>1) ss = dot_product(r(i,i+1:j),r(i+1:j,j))
             denom = r(i,i)+r(j,j)
             if (denom.eq.0._dp) stop 'denom is null!'
             r(i,j) = (t(i,j) - ss)/denom
          enddo
       enddo
    enddo

!Between-block interactions.
    do j=1,nblocks-1
       jstart = startStopPairs(j,1)
       jstop = startStopPairs(j,2)

!fix: was i=j-1,0,-1
       
       do i=j-1,1,-1
          istart = startStopPairs(i,1)
          istop = startStopPairs(i,1)
          allocate(s(istop-istart+1,jstop-jstart+1))
          allocate(x(istop-istart+1,jstop-jstart+1))
          s = t(istart:istop,jstart:jstop)
          if (j-i>1) s = s - matmul(r(istart:istop,istop:jstart) &
               , r(istop:jstart,jstart:jstop))
         
          call dp_real_sylvester(r(istart:istop, istart:istop) &
               ,r(jstart:jstop, jstart:jstop),s,x,scale)

          r(istart:istop,jstart:jstop) = x * scale

          deallocate(s,x)

       enddo
    enddo

  end subroutine uptriang_square_root



! Solves the real Sylvester matrix equation:
!
! op(A)*X + X*op(B) = scale*C
!
! or
!
! op(A)*X - X*op(B) = scale*C
!
! where op(A) = A or A**T, and A and B are both upper quasi-
! triangular. A is M-by-M and B is N-by-N; the right hand side C and
! the solution X are M-by-N; and scale is an output scale factor, set
! <= 1 to avoid overflow in X.  A and B must be in Schur canonical
! form (as returned by DHSEQR), that is, block upper triangular with
! 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has
! its diagonal elements equal and its off-diagonal elements of
! opposite sign.

  subroutine dp_real_sylvester(a,b,c,x,scale)
    implicit none
    real(dp), dimension(:,:), intent(in) :: a,b,c
    real(dp), dimension(size(a,1),size(b,1)), intent(out) :: x

    real(dp), dimension(:,:), allocatable :: ctemp
    real(dp) :: scale

    character, parameter :: TRANA = 'N'
    character, parameter :: TRANB = 'N'
    integer, parameter :: ISGN = 1

    integer :: INFO, LDA, LDB, LDC, M, N
  
    LDA = size(a,1)
    M = size(a,2)
    LDB = size(b,1)
    N = size(b,2)
    LDC = size(c,1)

    allocate(ctemp(LDC,N))

    ctemp = c

    call dtrsyl(TRANA,TRANB,ISGN,M,N,A,LDA,B,LDB,ctemp,LDC,SCALE,INFO)
    
    x = ctemp

    deallocate(ctemp)
    
  end subroutine dp_real_sylvester


! computes the singular value decomposition (SVD) of a real M-by-N
! matrix A, optionally computing the left and/or right singular
! vectors.
!
! A = U * SIGMA * transpose(V) where SIGMA is an M-by-N matrix which
! is zero except for its min(m,n) diagonal elements, U is an
! M-by-M orthogonal matrix, and V is an N-by-N orthogonal
! matrix.  The diagonal elements of SIGMA are the singular
! values of A; they are real and non-negative, and are returned
! in descending order.  The first min(m,n) columns of U and V
! are the left and right singular vectors of A.  Note that the
! routine returns V**T, not V.

  subroutine dp_svd(a,u,sdiag,vt)
    implicit none
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(size(a,1),size(a,1)), intent(out) :: u
    real(dp), dimension(size(a,2),size(a,2)), intent(out) :: vt
    real(dp), dimension(min(size(a,1),size(a,2))), intent(out) :: sdiag
    
    character, parameter :: JOBU = 'A'
    character, parameter :: JOBVT = 'A'
    integer :: M,N
    integer :: INFO, LDA, LDU, LDVT, LWORK

    real(dp), dimension(:), allocatable :: WORK
    real(dp), dimension(:,:), allocatable :: atemp

    M = size(a,1)
    N = size(a,2)

    LDA = M
    LDU = M
    LDVT = N

    LWORK=max( 3 * min ( m, n ) + max ( m, n ), 5 * min ( m, n ) )

    allocate(WORK(LWORK))
    allocate(atemp(m,n))

    atemp = a
    
    call dgesvd(JOBU, JOBVT, M, N, atemp, LDA, sdiag,  U,  LDU,  VT,  LDVT &
         ,WORK, LWORK, INFO )

    deallocate(WORK,atemp)

    if (INFO.ne.0) then
       write(*,*) 'INFO= ',INFO
       stop 'dp_svd: dgesvd failed!'
    endif
    
  end subroutine dp_svd



  function dp_diagpiv_linsys_symmetric(a,b)
    implicit none
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(size(b,1)) :: dp_diagpiv_linsys_symmetric

    character, parameter :: UPLO='U'
    integer, parameter :: NRHS = 1   

    integer :: n
    integer :: info, lwork
    integer, dimension(size(b,1)) :: ipiv
    integer, dimension(size(a,1)) :: iwork

    integer, parameter :: factor = 3
    real(dp) :: rcond
    real(dp), dimension(size(a,1),size(a,2)) :: abuf, af
    real(dp), dimension(size(b,1)) :: x, berr, ferr
    real(dp), dimension(factor*size(b,1)) :: work

    n = size(a,1)   
    if (size(a,2).ne.n) stop 'dp_diagpiv_linsys_symmetric: matrix not symmetric!'
    if (size(b,1).ne.n) stop 'dp_diagpiv_linsys_symmetric: b is not of size a!'
    
    x = b
    abuf = a
    lwork = factor*n
    
    call dsysv(UPLO, n, NRHS, abuf, n, ipiv, x, n, work, lwork, info)
    
    if (info.ne.0) then
       write(*,*)'dp_diagpiv_linsys_symmetric:'
       write(*,*)'info= ',info
       if (info.le.n) stop
    endif

    dp_diagpiv_linsys_symmetric = x
    
  end function dp_diagpiv_linsys_symmetric


!Divide and Conquer Algorithm for Eigenproblem with Symmetric matrices
  subroutine dp_divconq_eigen_symmetric(a,w,v)

    implicit none
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(size(a,1)), intent(out) :: w
    real(dp), dimension(size(a,1),size(a,2)), intent(out) :: v

    character, parameter :: JOBZ='V'
    character, parameter :: UPLO='U'
    integer :: INFO, LDA, LIWORK, LWORK, N

    integer, dimension(:), allocatable :: IWORK
    real(dp), dimension(:), allocatable :: WORK
    real(dp), dimension(:,:), allocatable :: atemp

    N = size(a,1)
    if (size(a,2).ne.N) stop 'dp_eigen_symmetric: a not squared!'

    LDA = N
    LWORK = 1 + 6*N + 2*N**2
    LIWORK = 3  +  5*N

    allocate(WORK(LWORK), IWORK(LIWORK))
    allocate(atemp(N,N))

    atemp = a
    
    call dsyevd(JOBZ, UPLO, N, atemp, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO)

    if (INFO.ne.0) stop 'dp_eigen_symmetric: dsyevd failed!'

    v = atemp

    deallocate(atemp)
    deallocate(WORK,IWORK)

  end subroutine dp_divconq_eigen_symmetric


!Implicitely Restarted Arnoldi Method for Eigenproblem with Symmetric matrices
  subroutine dp_iram_eigen_symmetric(a,evalues,evects,which,sig)
    implicit none
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(:), intent(out) :: evalues
    real(dp), dimension(size(a,1),size(evalues,1)), intent(out) :: evects
    character(len=2), intent(in) :: which
    real(dp), intent(in), optional :: sig

    integer :: ldv
    integer, parameter :: factor = 2

    real(dp), dimension(size(a,1),factor*size(evalues,1)) :: v
    real(dp), dimension(factor*size(evalues,1)*(factor*size(evalues,1)+8)) :: workl
    real(dp), dimension(3*size(a,1)) :: workd
    real(dp), dimension(factor*size(evalues,1),2) :: d
    real(dp), dimension(size(a,1)) :: resid

    real(dp), dimension(size(a,1),size(a,2)) :: sid, amsid

    logical, dimension(factor*size(evalues,1)) :: selection
    integer, dimension(11) :: iparam, ipntr
    
    character :: bmat
    integer :: ido, nev, ncv, n, lworkl, info, ierr
    integer :: i,j, nconv, maxitr, mode, ishfts

    logical, parameter :: rvec = .true.

    real(dp) :: sigma
    real(dp), parameter :: tol = epsilon(1._dp)


    logical, parameter :: display = .true.
    logical, parameter :: shownorm = .true.


!dimension
    n = size(a,1)
    if (size(a,2).ne.n) stop 'dp_iram_eigensym: a is not symmetric!'    
    ldv = n

!NEV  asks for the number of eigenvalues to be  computed.                             
    nev = size(evalues,1)
   
    if (nev.ge.n) then
       write(*,*)'dp_iram_eigen_symmetric:'
       write(*,*)'n= nev= ',n,nev
       stop 'restarted Arnoldi methods cannot compute all eigenvalues!'
    endif

!NCV sets the length of the Arnoldi factorization 
    ncv = factor*nev

    if (ncv.gt.n) stop 'dp_iram_eigensym: ncv > n!'



!This is a standard problem
    bmat ='I'


!Ask for the NEV eigenvalues of largest magnitude (indicated by which
!= 'LM') See documentation in DSAUPD for the other options SM, LA, SA,
!LI, SI.
    
    select case (which)
    case ('LM','SM','LA','SA','BE')

    case default
       stop 'dp_iram_eigensym: parameter WHICH incorrect!'
    end select

   
    lworkl = ncv*(ncv+8)
    info = 0
    ido = 0


! This program uses the exact shift strategy (indicated by setting
! PARAM(1) = 1).  IPARAM(3) specifies the maximum number of Arnoldi
! iterations allowed.  Mode 1 of DSAUPD is used (IPARAM(7) = 1). All
! these options can be changed by the user. For details see the
! documentation in DSAUPD.

    ishfts = 1
    maxitr = 1000

    if (present(sig)) then
!shit inverted
       sigma = sig
       mode = 3
    else
!regular
       mode = 1
    endif

    iparam(1) = ishfts 
    iparam(3) = maxitr 
    iparam(7) = mode 


    select case (mode)

!regular
    case (1)

       do 
          call dsaupd ( ido, bmat, n, which, nev, tol, resid &
            , ncv, v, ldv, iparam, ipntr, workd, workl &
            , lworkl, info )

          if ((ido.eq.1).or.(ido.eq.-1)) then
!             A.x -> y
             workd(ipntr(2):ipntr(2)+n-1) = matmul(a,workd(ipntr(1):ipntr(1)+n-1))
          else
             if (ido.ne.99) stop 'dp_iram_eigen_symmetric: ido >< 99'
             exit
          endif
          
       enddo

    case (3)
!sigma * id
       forall(i = 1:n, j = 1:n) sid(i,j) = sigma*((i/j)*(j/i))
     
       amsid = a - sid

       do 
          call dsaupd ( ido, bmat, n, which, nev, tol, resid &
            , ncv, v, ldv, iparam, ipntr, workd, workl &
            , lworkl, info )

          if ((ido.eq.1).or.(ido.eq.-1)) then
![(A-sigma I)^-1].x -> y
!nobody sane wants to calculate an inverse, let's solve [A-sigma I].y = x             
             workd(ipntr(2):ipntr(2)+n-1) = linsyssym(amsid,workd(ipntr(1):ipntr(1)+n-1))
          else
             if (ido.ne.99) stop 'dp_iram_eigen_symmetric: ido >< 99'
             exit
          endif
          
       enddo

    case default
       stop 'dp_iram_eigensym: mode not implemented!'

    end select


    if (info.lt.0) then
       write(*,*)'dp_iram_eigen_symmetric:'
       write(*,*)'error in dsaupd with info= ',info
       stop
    elseif (info.eq.1) then
       write(*,*)'dp_iram_eigen_symmetric:'
       write(*,*)'maximum number of iteration reached!',maxitr
       stop
    elseif (info.eq.3) then
       write(*,*)'dp_iram_eigen_symmetric:'
       write(*,*)'no shifts could be applied during implicit Arnoldi update!'
       write(*,*)'try increasing ncv: ncv= ',ncv
       stop
    endif


!Eigenvalues are returned in the first column of the two dimensional
!array D and the corresponding eigenvectors are returned in the first
!NEV columns of the two dimensional array V if requested.  Otherwise,
!an orthogonal basis for the invariant subspace corresponding to the
!eigenvalues in D is returned in V.

    call dseupd ( rvec, 'All', selection, d, v, ldv, sigma &
         , bmat, n, which, nev, tol, resid, ncv, v, ldv &
         , iparam, ipntr, workd, workl, lworkl, ierr )


    if (ierr.ne.0) then
       write(*,*)'dp_iram_eigen_symmetric:'
       write(*,*)'error in dseupd with ierr= ',ierr
       stop
    else
       if (display) then
          write(*,*)'dp_iram_eigen_symmetric:'
          write(*,*)'implicit Arnoldi iter=    ',iparam(3)
          write(*,*)'number of OP.x operation= ',iparam(9)
          write(*,*)'number of converged Ritz= ',iparam(5)
          write(*,*)'which= ncv= tol=          ',which,ncv,tol
       endif
    end if

    if (shownorm) then
       do i=1,nev
          write(*,*)'||A.x - l.x||= ',norm2(matmul(a,v(1:n,i))-d(i,1)*v(1:n,i))
       enddo
    endif


    evalues(1:nev)=d(1:nev,1)
    evects(1:n,1:nev)=v(1:n,1:nev)


  end subroutine dp_iram_eigen_symmetric


 
!Implicitely Restarted Arnoldi Method for Eigenproblem with BAND STORAGE Symmetric matrices
  subroutine dp_iram_eigen_band_symmetric(ab,kl,ku,evalues,evects,which,sig)
    use arpex, only : dsband
    implicit none
    real(dp), dimension(:,:), intent(in) :: ab
    integer, intent(in) :: kl,ku
    real(dp), dimension(:), intent(out) :: evalues
    real(dp), dimension(size(ab,2),size(evalues,1)), intent(out) :: evects

    character(len=2), intent(in) :: which
    real(dp), intent(in), optional :: sig

    integer :: lda, ldv, ldz
    integer, parameter :: factor = 2

    real(dp), dimension(size(evalues,1),2) :: d
    real(dp), dimension(size(ab,2),factor*size(evalues,1)) :: z
    real(dp), dimension(size(ab,1),size(ab,2)) :: mb,rfac

    real(dp), dimension(size(ab,2)) :: resid
    real(dp), dimension(size(ab,2),factor*size(evalues,1)) :: v

    real(dp), dimension(3*size(ab,2)) :: workd
    real(dp), dimension(:), allocatable :: workl, ax


    logical, dimension(factor*size(evalues,1)) :: selection

    integer, dimension(size(ab,2)) :: iwork
    integer, dimension(11) :: iparam, ipntr
    
    character :: bmat
    integer :: ido, nev, ncv, n, lworkl, info, ierr
    integer :: i,j, nconv, maxitr, mode, ishfts

    logical, parameter :: rvec = .true.

    real(dp) :: sigma
    real(dp), parameter :: tol = epsilon(1._dp)


    logical, parameter :: display = .true.
    logical, parameter :: shownorm = .true.


!dimension
    lda = size(ab,1)
    n = size(ab,2)
    ldv = n
    ldz = n

!NEV  asks for the number of eigenvalues to be  computed.                             
    nev = size(evalues,1)
   
    if (nev.ge.n) then
       write(*,*)'dp_iram_eigen_band_symmetric:'
       write(*,*)'n= nev= ',n,nev
       stop 'restarted Arnoldi methods cannot compute all eigenvalues!'
    endif

!NCV sets the length of the Arnoldi factorization 
    ncv = factor*nev

    if (ncv.gt.n) stop 'dp_iram_band_eigensym: ncv > n!'



!This is a standard problem
    bmat ='I'


!Ask for the NEV eigenvalues of largest magnitude (indicated by which
!= 'LM') See documentation in DSAUPD for the other options SM, LA, SA,
!LI, SI.
    
    select case (which)
    case ('LM','SM','LA','SA','BE')

    case default
       stop 'dp_iram_eigensym: parameter WHICH incorrect!'
    end select

   
    info = 0
    ido = 0


! This program uses the exact shift strategy (indicated by setting
! PARAM(1) = 1).  IPARAM(3) specifies the maximum number of Arnoldi
! iterations allowed.  Mode 1 of DSAUPD is used (IPARAM(7) = 1). All
! these options can be changed by the user. For details see the
! documentation in DSAUPD.

    ishfts = 1
    maxitr = 1000

    if (present(sig)) then
!shit inverted
       sigma = sig
       mode = 3
       lworkl = 3*ncv**2+6*ncv
    else
!regular
       mode = 1
       lworkl = ncv*(ncv+8)
    endif

    allocate(workl(lworkl))

    iparam(1) = ishfts 
    iparam(3) = maxitr 
    iparam(7) = mode 


!for shifted-inverse method, dsband call a LU factorization of A-sigma
!I; which is less robust for singular matrices than a diagonal
!pivoting method. So, you might get in these situations errors larger
!than using the shifted-inverse method directely on non-band storage
!matrices.
   
    call dsband ( rvec, 'A', selection, d, v, ldv, sigma, n, ab, mb, lda &
         , rfac, kl, ku, which, bmat, nev, tol &
         , resid, ncv, v, ldv, iparam, workd, workl, lworkl &
         , iwork, info)
      

    if (info.lt.0) then
       write(*,*)'dp_iram_eigen_band_symmetric:'
       write(*,*)'error in dsaupd with info= ',info
       stop
    elseif (info.eq.1) then
       write(*,*)'dp_iram_eigen_band_symmetric:'
       write(*,*)'maximum number of iteration reached!',maxitr
       stop
    elseif (info.eq.3) then
       write(*,*)'dp_iram_eigen_band_symmetric:'
       write(*,*)'no shifts could be applied during implicit Arnoldi update!'
       write(*,*)'try increasing ncv: ncv= ',ncv
       stop
    endif

    if (display) then
       write(*,*)'dp_iram_eigen_band_symmetric:'
       write(*,*)'implicit Arnoldi iter=    ',iparam(3)
       write(*,*)'number of OP.x operation= ',iparam(9)
       write(*,*)'number of converged Ritz= ',iparam(5)
       write(*,*)'which= ncv= tol=          ',which,ncv,tol
    endif

    if (shownorm) then
       allocate(ax(n))
       do i=1,nev
          call dgbmv ('N',n, n, kl, ku, 1._dp,ab(kl+1,1), lda, v(1,i) &
               ,1,0._dp,ax,1)
     
          write(*,*)'||A.x - l.x||= ',norm2(ax(1:n)-d(i,1)*v(1:n,i))
       enddo
       deallocate(ax)
    endif


    evalues(1:nev)=d(1:nev,1)
    evects(1:n,1:nev)=v(1:n,1:nev)

    deallocate(workl)

  end subroutine dp_iram_eigen_band_symmetric



end module linalg
