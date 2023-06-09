module algprs
  use bvpprec
  implicit none
  public
  integer,save :: nminit =8
  integer,save :: iprint =0
  integer,save :: idum=0
  real(dp), save :: uval0 = 0._dp
  logical,save :: pdebug =.false.
  logical,save :: use_c =.true.
  logical,save :: comp_c =.true.
end module algprs


module sepbvp
  use bvpprec
  implicit none
  
  private

  public easybvpc



contains


  subroutine easybvpc(nEq,nBc,r,y,alltol)
    use algprs, only : uval0
    implicit none
    integer, intent(in) :: nEq,nBc
    real(dp), dimension(:), pointer :: r
    real(dp), dimension(:,:), pointer ::y
    real(dp), intent(in) :: alltol

!internal working arrays of twpbvpbc      
    integer, parameter :: nxxdim = 10000000
    integer, save :: lwrkfl= 10000000,lwrkin=10000000

    integer, parameter :: lwrkMax=100000000

!input of twpbvpbc
    logical :: linear, givmsh, giveu
    integer :: ncomp,nlbc
    integer :: nfxpnt
    integer :: ntol
    integer :: nmsh,nmax
    integer :: nudim, iflbvp

    real(dp) :: aleft, aright
    real(dp), dimension(1) :: fixpnt

    integer, dimension(:), allocatable :: ltol
    real(dp), dimension(:), allocatable :: tol
    
    real(dp), dimension(:,:), allocatable :: u
    real(dp), dimension(:), allocatable :: xx


    real(dp), dimension(:), allocatable :: wrk
    integer, dimension(:), allocatable :: iwrk

    real(dp) :: ckappa,ckappa1,gamma1
    real(dp), dimension(1) :: rpar
    integer, dimension(1) :: ipar
    
!local scope
    logical, parameter :: checkInput=.true.
    integer :: i,j,n
    
!interface to the user supplied functions
    include 'sepbvp.h'

    if ((.not.associated(r)).or.(.not.associated(y))) then
       stop 'easybvp: trial solutions not found!'
    else
       n = size(r,1)
    endif

    if (checkInput) then
       if (n.gt.nxxdim) then
          write(*,*)'nmax= ',nxxdim
          stop 'r-vector too large'
       endif
       if ((size(y,1).ne.neq).or.(size(y,2).ne.n)) then
          stop 'y-vector size incorrect!'
       endif
       if ((nBc.gt.neq).or.(nBc.lt.0)) then
          write(*,*)'Improper number of boundary conditions!'
          write(*,*)'nEq= nBc= ',nEq,nBc
          stop
       endif
       if (r(1).ge.r(n)) stop 'rmin > rmax!'
    endif



!boundary points
    aleft=r(1)
    aright=r(n)

!no fix point added on the mesh
    nfxpnt=0
    fixpnt(1)=(aleft+aright)/2.

!tolerance vector
    ntol=neq
    allocate(ltol(ntol),tol(ntol))

    do i=1,ntol
       ltol(i)=i
       tol(i)=alltol
    enddo

!mesh and array sizes        
    nlbc=nBc
    ncomp=neq
    nudim=ncomp
    nmsh=n
    
!user given mesh
    allocate(xx(1:nxxdim))
    xx(1:nmsh)=r(1:nmsh)

!the solution vector
    allocate(u(nudim,nxxdim))
    u(1:ncomp,1:nmsh)=y(1:neq,1:n)


    linear = .false.
    givmsh=.true.
    giveu=.true.

    do

       allocate(wrk(lwrkfl))
       allocate(iwrk(lwrkin))

       call twpbvpc(ncomp,nlbc,aleft,aright,nfxpnt,fixpnt,ntol,ltol,tol &
            ,linear,givmsh,giveu,nmsh,nxxdim,xx,nudim,u,nmax,lwrkfl,wrk &
            ,lwrkin,iwrk,fsub,dfsub,gsub,dgsub &
            ,ckappa1,gamma1,ckappa,rpar,ipar,iflbvp)
    
       if (nmax.gt.nxxdim) stop 'nmax > nxxdim'
       write(*,*)'twpbvpc on exit: iflbvp= ',iflbvp
       write(*,*)'nmsh= ',nmsh
       write(*,*)'nmax= ',nmax

       select case (iflbvp)

       case (0)
          if (nmsh.ne.n) then
             deallocate(r,y)
             allocate(r(nmsh),y(ncomp,nmsh))
             r = xx(1:nmsh)
             y(1:ncomp,1:nmsh)=u(1:ncomp,1:nmsh)
          endif

          exit

       case (-1)
          stop 'easybvp: invalid input parameters!'
          
       case (1)
          write(*,*) 'easybvp: maximal numbers of mesh points exceeded!'

!start again with bigger working arrays          
          write(*,*) '...increasing working array sizes...'
          lwrkfl = lwrkfl*5
          lwrkin = lwrkin*5
          write(*,*) 'lwrkfl= lwrkin= ',lwrkfl,lwrkin
          deallocate(wrk,iwrk)
!reset to the trial solution
          nlbc=nBc
          ncomp=neq
          nudim=ncomp
          nmsh=n
          xx(1:nmsh)=r(1:nmsh)
          u(1:ncomp,1:nmsh)=y(1:neq,1:n)
          call random_number(uval0)
!check if we are not killing all memory
          if (max(lwrkin,lwrkfl).gt.lwrkMax) exit

          cycle

       case default
          stop 'easybvp: internal error!'
          
       end select

    enddo

    deallocate(ltol,tol,u,xx)

  end subroutine easybvpc
    


  

!-------------------------------------------------------------------------
!encapsulated into f90 module.
!block data removed + algprs common block converted to algprs module.
!added implicit declaration of integer in between i-o
!external statements replaced by interfaces + include statements.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!credit: http://www.ma.ic.ac.uk/~jcash/BVP_software/readme.php#docs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------------------

  

  subroutine twpbvpc(ncomp, nlbc, aleft, aright,                        &
     &       nfxpnt, fixpnt, ntol, ltol, tol,                           &
     &       linear, givmsh, giveu, nmsh,                               &
     &       nxxdim, xx, nudim, u, nmax,                                &
     &       lwrkfl, wrk, lwrkin, iwrk,                                 &
     &       fsub, dfsub, gsub, dgsub,                                  &
     &       ckappa1,gamma1,ckappa,rpar,ipar,iflbvp)                    
!                                                                       
!  The subroutine twpbvpc is intended to solve two-point boundary       
!  value problems.                                                      
!                                                                       
!                                                                       
! References:                                                           
!  Cash, J. R.; Mazzia, F.                                              
!   A new mesh selection algorithm, based on conditioning,              
!   for two-point boundary value codes.                                 
!   J. Comput. Appl. Math.  184  (2005), no. 2, 362--381.               
!                                                                       
!   Revision History                                                    
!                                                                       
! revision  July 3,  2006                                               
!   added rpar and ipar in the functions                                
!   DFSUB, DGSUB, FSUB, GSUB                                            
!   changed the name of the variable double in ddouble                  
!                                                                       
! revision   September 20, 2004                                         
!   This is a modified version of twpbvp that uses the conditioning     
!   in the mesh selection.                                              
!                                                                       
!     New subroutines not included in the old version:                  
!           condestim                                                   
!           estimkappa                                                  
!           moncond                                                     
!           selcondmsh                                                  
!           selconderrmsh                                               
!           smpselcondmsh                                               
!                                                                       
!     Updates subroutines:                                              
!           bvpsol                                                      
!           conv4                                                       
!           fail4                                                       
!           fail6                                                       
!           conv8                                                       
!           fail8                                                       
!           newteq                                                      
!           mshref                                                      
!           wtchdg                                                      
!                                                                       
!     Auxiliary function not used by the old version                    
!           donest                                                      
!           abdnrm                                                      
!                                                                       
!     The common block algprs contains two more variable                
!     logical pdebug, use_c, comp_c                                     
!     common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
!                                                                       
! The starting subroutine is twpbvp by J.R. Cash and M.H. Wright        
!                                                                       
                                                                        
 !module
      use algprs
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension rpar(*),ipar(*) 
      dimension fixpnt(*), ltol(*), tol(*) 
      dimension xx(*), u(nudim,*) 
      dimension wrk(lwrkfl), iwrk(lwrkin) 
                                                                        
      logical linear, givmsh, giveu 
!      external fsub, dfsub, gsub, dgsub 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      intrinsic abs, min 
                                                                        
      parameter ( zero = 0.0_dp ) 
 
!interface
      include 'sepbvp.h'
      
                                                                       
!  Check for invalid input parameters.  If any parameters are           
!  invalid, exit with the flag iflbvp set to -1.                        
                                                                        
      iflbvp = -1 
      if (ncomp .le. 0)  return 
      if (nlbc .lt. 0 .or. nlbc .gt. ncomp) return 
      if (aleft .ge. aright) return 
                                                                        
      if (nfxpnt .lt. 0)  return 
      if (givmsh .and. nmsh .lt. nfxpnt+2) return 
      if (givmsh .and. xx(1) .ne. aleft) return 
!     SCMODIFIED add an extra condition to avoid accessing xx(0)        
      if (nmsh .gt. 0) then 
        if (givmsh .and. xx(nmsh) .ne. aright) return 
      end if 
      if (nfxpnt .gt. 0) then 
         if (fixpnt(1) .le. aleft) return 
         if (fixpnt(nfxpnt) .ge. aright) return 
         do 50 i = 1, nfxpnt-1 
            if (fixpnt(i+1) .le. fixpnt(i)) return 
   50    continue 
      endif 
                                                                        
      if (ntol .lt. 1) return 
      do 60 i = 1, ntol 
         if (ltol(i) .lt. 0 .or. ltol(i) .gt. ncomp) return 
         if (tol(i) .le. zero) return 
   60 continue 
                                                                        
      if (giveu .and. .not. givmsh) return 
                                                                        
      if (use_c .and. .not. comp_c) return 
                                                                        
      if (nudim .le. 0) return 
      if (lwrkfl .le. 0 .or. lwrkin .le. 0) return 
                                                                        
                                                                        
!  Calculate maximum number of mesh points possible with the            
!  given floating-point and integer workspace.                          
                                                                        
      isp = lwrkfl - 2*ntol - 22*ncomp - 6*ncomp*ncomp 
      iden = 6*ncomp*ncomp + 22*ncomp + 3 
      nmax1 = isp/iden 
                                                                        
      isp = lwrkin - 2*ncomp 
      nmax2 = isp/(2*ncomp+3) 
                                                                        
      nmax = min(nmax1,nmax2) 
! nmax from workspace                                                   
      nmax = min(nmax, nxxdim) 
! nmax from size of u and xx                                            
                                                                        
      if (iprint .ge. 0) write(6,901) nmax 
  901 format(1h ,'nmax from workspace =',i8) 
                                                                        
      if (nmax .le. 1) return 
                                                                        
                                                                        
!  Partition floating point workspace.                                  
                                                                        
      irhs = 1 
      lrhs = ncomp*nmax 
                                                                        
      itpblk = irhs + lrhs 
      ltpblk = ncomp*nlbc 
                                                                        
      ibtblk = itpblk + ltpblk 
      lbtblk = ncomp*(ncomp - Nlbc) 
                                                                        
      iajac = ibtblk + lbtblk 
      lajac = 2*ncomp*ncomp*nmax 
                                                                        
      ibhold = iajac + lajac 
      lbhold = ncomp*ncomp*nmax 
                                                                        
      ichold = ibhold + lbhold 
      lchold = ncomp*ncomp*nmax 
                                                                        
      ifval = ichold + lchold 
      lfval = ncomp*nmax 
                                                                        
      idef = ifval + lfval 
      ldef = ncomp*(nmax-1) 
                                                                        
      idefex = idef + ldef 
      ldefex = ncomp*(nmax-1) 
                                                                        
!  def6 uses the same space as defexp                                   
                                                                        
      idef6 = idefex 
      ldef6 = ncomp*(nmax-1) 
                                                                        
      idefim = idef6 + ldef6 
      ldefim = ncomp*(nmax-1) 
                                                                        
!  def8 uses the same space as defimp                                   
                                                                        
      idef8 = idefim 
      ldef8 = ncomp*(nmax-1) 
                                                                        
      iusve = idef8 + ldef8 
      lusve = ncomp*nmax 
                                                                        
      iuold = iusve + lusve 
      luold = ncomp*nmax 
                                                                        
      itmrhs = iuold + luold 
      ltmrhs = ncomp*nmax 
                                                                        
      irhtri = itmrhs + ltmrhs 
      lrhtri = ncomp*nmax 
                                                                        
      idelu = irhtri + lrhtri 
      ldelu = ncomp*nmax 
                                                                        
      ixmer = idelu + ldelu 
      lxmer = ncomp*nmax 
                                                                        
!  rerr occupies the same space as xmerit                               
      irerr = ixmer 
      lrerr = ncomp*nmax 
                                                                        
      iutri = irerr + lrerr 
      lutri = ncomp*nmax 
                                                                        
      iermx = iutri + lutri 
      lermx = nmax 
                                                                        
      irtdc = iermx + lermx 
      lrtdc = nmax 
                                                                        
      ixxold = irtdc + lrtdc 
      lxxold = nmax 
                                                                        
      iuint = ixxold + lxxold 
      luint = ncomp 
                                                                        
      iftmp = iuint + luint 
      lftmp = ncomp 
                                                                        
      idgtm = iftmp + lftmp 
      ldgtm = ncomp 
                                                                        
      idftm1 = idgtm + ldgtm 
      ldftm1 = ncomp*ncomp 
                                                                        
      idftm2 = idftm1 + ldftm1 
      ldftm2 = ncomp*ncomp 
                                                                        
      itmp = idftm2 + ldftm2 
      ltmp = ncomp*8 
                                                                        
      idsq = itmp + ltmp 
      ldsq = ncomp*ncomp 
                                                                        
      idexr = idsq + ldsq 
      ldexr = ncomp 
                                                                        
      ietst6 = idexr + ldexr 
      letst6 = ntol 
                                                                        
      ietst8 = ietst6 + letst6 
      letst8 = ntol 
                                                                        
      iamg = ietst8 + letst8 
      lamg = ncomp*nmax 
                                                                        
      ic1 = iamg + lamg 
      lc1 = ncomp*ncomp*nmax 
                                                                        
                                                                        
      iwrkrhs = ic1+lc1 
      lwrkrhs = ncomp*nmax 
                                                                        
      ilast = iwrkrhs +  lwrkrhs 
                                                                        
                                                                        
      if (iprint .eq. 1) write(6,903) ilast 
  903  format(1h ,'ilast',i10) 
                                                                        
                                                                        
!  Partition integer workspace.                                         
                                                                        
      iiref = 1 
      liref = nmax 
                                                                        
      iihcom = iiref + liref 
      lihcom = nmax 
                                                                        
      iipvbk = iihcom + lihcom 
      lipvbk = ncomp*nmax 
                                                                        
      iipvlu = iipvbk + lipvbk 
      lipvlu = ncomp 
                                                                        
      iisign = iipvlu + lipvlu 
      lisign = ncomp*nmax 
                                                                        
      call bvpsol(ncomp, nmsh, nlbc, aleft, aright,                     &
     &   nfxpnt, fixpnt, ntol, ltol, tol, nmax, linear,                 &
     &   giveu, givmsh, xx, nudim, u,                                   &
     &   wrk(idefex), wrk(idefim), wrk(idef), wrk(idelu),               &
     &   wrk(irhs), wrk(ifval),                                         &
     &   wrk(itpblk), wrk(ibtblk), wrk(iajac), wrk(ibhold),             &
     &   wrk(ichold), iwrk(iipvbk), iwrk(iipvlu),iwrk(iisign),          &
     &   wrk(iuint), wrk(iftmp), wrk(itmrhs),                           &
     &   wrk(idftm1), wrk(idftm2), wrk(idgtm),                          &
     &   wrk(iutri), wrk(irhtri), wrk(ixmer),                           &
     &   wrk(ixxold), wrk(iuold), wrk(iusve),                           &
     &   wrk(itmp), wrk(idsq), wrk(idexr), wrk(irtdc),                  &
     &   wrk(irerr), wrk(ietst6), wrk(ietst8), wrk(iermx),              &
     &   iwrk(iihcom), iwrk(iiref), wrk(idef6), wrk(idef8),             &
     &   fsub,dfsub,gsub,dgsub,iflbvp,                                  &
     &   wrk(iamg),wrk(ic1),wrk(iwrkrhs),                               &
     &   ckappa1,gamma1,ckappa,rpar, ipar)                              
      return 
      end subroutine                                           
                                                                        
                                                                        
                                                                        
      subroutine bvpsol(ncomp, nmsh, nlbc, aleft, aright,               &
     &   nfxpnt, fixpnt,                                                &
     &   ntol, ltol, tol, nmax, linear, giveu, givmsh,                  &
     &   xx, nudim, u, defexp, defimp, def, delu, rhs, fval,            &
     &   topblk, botblk, ajac, bhold, chold, ipvblk, ipivlu,isign,      &
     &   uint, ftmp, tmprhs, dftmp1, dftmp2, dgtm,                      &
     &   utrial, rhstri, xmerit, xxold, uold, usave,                    &
     &   tmp, dsq, dexr, ratdc, rerr,                                   &
     &   etest6, etest8, ermx, ihcomp, irefin,                          &
     &   def6, def8, fsub, dfsub, gsub, dgsub, iflbvp,                  &
     &   amg, c1, wrkrhs,ckappa1,gamma1,ckappa,rpar,ipar)               
!module
      use algprs
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension rpar(*), ipar(*) 
      dimension  fixpnt(*), ltol(ntol), tol(ntol) 
      dimension  xx(*), u(nudim, *) 
      dimension  defexp(ncomp,*), defimp(ncomp,*), def(ncomp,*) 
      dimension  delu(ncomp, *), rhs(*), fval(ncomp,*) 
      dimension  topblk(nlbc,*), botblk(ncomp-nlbc,*) 
      dimension  ajac(ncomp, 2*ncomp, *) 
      dimension  bhold(ncomp, ncomp, *), chold(ncomp, ncomp, *) 
      dimension  ipivlu(*), ipvblk(*), isign(*) 
      dimension  uint(ncomp), ftmp(ncomp) 
      dimension  dgtm(ncomp), tmprhs(*) 
      dimension  dftmp1(ncomp, ncomp), dftmp2(ncomp, ncomp) 
      dimension  utrial(ncomp,*), rhstri(*) 
      dimension  xmerit(ncomp, *) 
      dimension  xxold(*), uold(ncomp,*), usave(ncomp,*) 
      dimension  tmp(ncomp,8) 
      dimension  dsq(ncomp,ncomp), dexr(ncomp) 
      dimension  ratdc(*), rerr(ncomp,*) 
      dimension  etest6(*), etest8(*), ermx(*) 
      dimension  ihcomp(*), irefin(*) 
      dimension  def6(ncomp,*), def8(ncomp,*) 
      dimension  amg(*), c1(ncomp,*), wrkrhs(*) 
                                                                        
      logical linear, giveu, givmsh, ddouble 
                                                                        
!      external fsub, dfsub, gsub, dgsub 
                                                                        
      common/mchprs/flmin, flmax, epsmch 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      intrinsic max 
                                                                        
      logical smooth, succes, strctr, trst6, reaft6 
      logical onto6, onto8, ludone, rhsgiv, maxmsh 
      logical first4, first8 
                                                                        
      logical mchset 
      save mchset 
                                                                        
      logical frscal 
!      save frscal                                                      
                                                                        
                                                                        
      logical stab_kappa, stab_gamma, stab_cond, stiff_cond 
                                                                        
      parameter (zero = 0.0_dp, one = 1.0_dp) 
      parameter (third = 0.33_dp, fourth = 0.25_dp) 
      parameter (quan6 = 0.1_dp ) 
      integer, parameter ::itcondmax = 10
                                                                        
!  blas: dload                                                          
!  real(dp) d1mach                                              
!interface
      include 'sepbvp.h'

      data mchset/.true./ 
      data fxfct/10.0_dp/ 
      data maxmsh/.false./ 
                                                                        
                                                                        
      frscal = .true. 
      if (mchset) then 
         flmin = d1mach(1) 
         flmax = d1mach(2) 
         epsmch = d1mach(3) 
         if (pdebug) write(6,901) epsmch 
         mchset = .false. 
      endif 
                                                                        
!  The routine stcons calculates integration constants stored in        
!  labeled common consts.                                               
                                                                        
      call stcons 
                                                                        
!  Set up arrays for the error tests.                                   
                                                                        
                                                                        
      if (.not. linear) then 
         call dload(ntol, one, etest6, 1) 
      else 
         do 10 i = 1, ntol 
            etest6(i) = one/max(quan6, tol(i)**third) 
   10    continue 
      endif 
      nmold = 1 
      smooth = .false. 
      strctr = .false. 
      trst6 = .true. 
      reaft6 = .false. 
      numbig = 0 
      nummed = 0 
      first4 = .true. 
      first8 = .true. 
      onto6  = .false. 
                                                                        
!     initialize parameter for the conditioning estimation              
      itcond  = 0 
      gamma1old  = flmax 
      gamma1     = flmax 
      ckappa1old  = flmax 
      ckappa1     = flmax 
      ckappa      = flmax 
      stiff_cond = .false. 
      stab_cond  = .false. 
      tolmin = flmax 
      do i=1,ntol 
         tolmin = min(tol(i),tolmin) 
      end do 
!  If givmsh is .true., the initial number of mesh points must be       
!  provided by the user in nmsh, and the mesh points must be            
!  contained in the array xx (of dimension nmsh).                       
!  Otherwise, nmsh is set to its default value, and a                   
!  uniform initial mesh is created.                                     
                                                                        
      if (.not. giveu .and. .not. givmsh) then 
         nmsh = nminit 
         if (nmsh .lt. nfxpnt+2) nmsh = nfxpnt + 2 
         call unimsh(nmsh, aleft, aright, nfxpnt, fixpnt, xx) 
      endif 
      if (pdebug) then 
         write(6,902) 
         call sprt(nmsh, xx) 
      endif 
                                                                        
      if (.not. giveu) call initu(ncomp, nmsh, xx, nudim, u,            &
     &     rpar,ipar  )                                                 
                                                                        
                                                                        
!**** top of logic for 4th order solution ****                          
                                                                        
  400 continue 
      if (iprint .eq. 1) write(6,903) nmsh 
                                                                        
!  Set the def (deferred correction) array to zero.                     
                                                                        
      call mtload(ncomp, nmsh-1, zero, ncomp, def) 
      iorder = 4 
                                                                        
!  The routine fneval calls fsub at the mesh xx and the                 
!  solution u, and saves the values in the array fval.                  
                                                                        
      call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,rpar,ipar) 
                                                                        
                                                                        
!  Try to compute a 4th order solution by solving a system of nonlinear 
!  equations.                                                           
                                                                        
      if (linear) then 
         ludone = .false. 
                                                                        
          call lineq( ncomp, nmsh, nlbc,                                &
     &    ludone, xx, nudim, u, def,                                    &
     &    delu, rhs, fval,                                              &
     &    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,                     &
     &    ajac, topblk, botblk, bhold, chold, ipvblk,                   &
     &    fsub, dfsub, gsub, dgsub, iflnwt,rpar,ipar)                   
                                                                        
!  Call fneval to evaluate the fval array at the new solution u.        
!  (Such a call is not necessary for the nonlinear case because         
!  fval is called within newteq for the new u.)                         
                                                                        
         call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,rpar,ipar) 
                                                                        
      else 
                                                                        
         rhsgiv = .false. 
         call newteq(ncomp, nmsh, nlbc,                                 &
     &    rhsgiv, ntol, ltol, tol,                                      &
     &    xx, nudim, u, def,                                            &
     &    delu, rhs, fval,                                              &
     &    utrial, rhstri,                                               &
     &    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,             &
     &    ajac, topblk, botblk, bhold, chold, ipvblk,                   &
     &    fsub, dfsub, gsub, dgsub, itnwt, iflnwt,rpar,ipar,frscal)     
                                                                        
      endif 
!                                                                       
!  these flags are used in the mesh selection strategy                  
!                                                                       
                                                                        
      if (iflnwt .eq. 0) then 
                                                                        
                                                                        
!       COMPUTE ESTIMATIONS OF CONDITION NUMBERS KAPPA1 and GAMMA1      
!       related to perturbation to the boundary conditions              
!                                                                       
                                                                        
        N =nmsh*ncomp 
        ninter=nmsh-1 
        if (comp_c) then 
          gamma1old = gamma1 
          ckappa1old = ckappa1 
          call CONDESTIM(aleft,aright,nmsh,ncomp,N,xx,topblk,nlbc,      &
     &      ncomp, ajac, ncomp,2*ncomp,ninter,botblk,ncomp-nlbc,        &
     &      ipvblk,isign,amg,c1,wrkrhs,ckappa1,gamma1)                  
                                                                        
          if (iprint .ge. 0) then 
!                                                                       
!     COMPUTE ESTIMATION OF THE CONDITION NUMBER KAPPA                  
!                                                                       
            call ESTIMKAPPA(nmsh,ncomp,N,xx,topblk,                     &
     &            nlbc,ncomp,ajac, ncomp,2*ncomp,                       &
     &          ninter,botblk,ncomp-nlbc,ipvblk,isign,c1,wrkrhs,ckappa) 
          end if 
                                                                        
          if (iprint .ge. 0) then 
            write(6,1001) ckappa1/gamma1 
            write(6,1002) gamma1 
            write(6,1003) ckappa1 
            write(6,1004) ckappa 
          end if 
                                                                        
          stab_kappa = abs(ckappa1old-ckappa1)/(1d0+ckappa1).lt.5d-2    &
     &      .and. ckappa1 .lt. flmax .and. gamma1.lt.flmax              
                                                                        
          stab_gamma = abs(gamma1old-gamma1)/(1d0+gamma1).lt.5d-2       &
     &     .and. gamma1.lt.flmax .and. ckappa1 .lt. flmax               
                                                                        
          stab_cond = stab_kappa .and. stab_gamma 
                                                                        
          stiff_cond = (( (ckappa1/gamma1 .ge. 10d0  ))) 
                                                                        
                                                                        
           if (iprint .eq. 1) then 
             write(6,*) 'stab_kappa = ', stab_kappa 
             write(6,*) 'stab_gamma = ', stab_gamma 
             write(6,*) 'stiff_cond = ', stiff_cond 
           end if 
                                                                        
         end if 
! endif if (comp_c)                                                     
                                                                        
                                                                        
         call conv4( ncomp, nmsh, ntol, ltol, tol, linear, nmax,        &
     &           xx, nudim, u, defexp, defimp, def, fval,               &
     &           tmp, bhold, chold, dsq, dexr, usave,                   &
     &           ratdc, rerr, ipivlu, nmold, xxold,                     &
     &           smooth, reaft6, onto6, strctr, trst6, ddouble ,        &
     &           fsub, maxmsh, succes, first4,                          &
     &           amg,stab_cond,ckappa,stiff_cond,wrkrhs,                &
     &              nfxpnt, fixpnt, irefin,rpar,ipar)                   
                                                                        
                                                                        
        if (pdebug .and. .not. onto6) write (6,904) 
                                                                        
      else 
                                                                        
       if (comp_c) then 
         if (iflnwt .ne. -1) then 
           gamma1old = gamma1 
           ckappa1old = ckappa1 
           N =nmsh*ncomp 
           ninter=nmsh-1 
           call CONDESTIM(aleft,aright,nmsh,ncomp,N,xx,topblk,nlbc,     &
     &       ncomp, ajac, ncomp,2*ncomp,ninter,botblk,ncomp-nlbc,       &
     &       ipvblk,isign,amg,c1,wrkrhs,ckappa1,gamma1)                 
           if (iprint .ge.0) then 
            call ESTIMKAPPA(nmsh,ncomp,N,xx,topblk,                     &
     &            nlbc,ncomp,ajac, ncomp,2*ncomp,                       &
     &          ninter,botblk,ncomp-nlbc,ipvblk,isign,c1,wrkrhs,ckappa) 
           end if 
         end if 
                                                                        
         if (iprint .ge. 0) then 
            write(6,1001) ckappa1/gamma1 
            write(6,1002) gamma1 
            write(6,1003) ckappa1 
            write(6,1004) ckappa 
         end if 
                                                                        
         stab_kappa = abs(ckappa1old-ckappa1)/(1d0+ckappa1).lt.5d-2     &
     &      .and. ckappa1 .lt. flmax .and. gamma1.lt.flmax              
                                                                        
         stab_gamma = abs(gamma1old-gamma1)/(1d0+gamma1).lt.5d-2        &
     &     .and. gamma1.lt.flmax .and. ckappa1 .lt. flmax               
                                                                        
         stab_cond = stab_kappa .and. stab_gamma 
         stiff_cond = (( (ckappa1/gamma1 .ge. 5d0  ))) 
                                                                        
         if (iprint .eq. 1) then 
           write(6,*) 'stab_kappa = ',stab_kappa 
           write(6,*) 'stab_gamma = ', stab_gamma 
           write(6,*) 'stiff_cond = ', stiff_cond 
         end if 
       end if 
! end if if(comp_c)                                                     
                                                                        
         succes = .false. 
         onto6 = .false. 
         reaft6 = .false. 
                                                                        
                                                                        
         call fail4( ncomp, nmsh, nlbc, ntol, ltol,                     &
     &           xx, nudim, u, rhs, linear, nmax,                       &
     &           nmold, xxold, uold, ratdc,                             &
     &           iorder, iflnwt, itnwt, ddouble , maxmsh,               &
     &           numbig, nummed,wrkrhs,amg,stab_cond,stiff_cond,        &
     &              nfxpnt, fixpnt, irefin,itcond,itcondmax,rpar,ipar)  
                                                                        
                                                                        
      endif 
                                                                        
      if (succes) then 
         if (comp_c) then 
           call ESTIMKAPPA(nmsh,ncomp,N,xx,topblk,                      &
     &            nlbc,ncomp,ajac, ncomp,2*ncomp,                       &
     &          ninter,botblk,ncomp-nlbc,ipvblk,isign,c1,wrkrhs,ckappa) 
         end if 
          if (iprint .ne. -1 .and. comp_c ) then 
            if (ckappa .ge. tolmin/epsmch) write(6,1005) 
          end if 
          iflbvp = 0 
          return 
      elseif (maxmsh) then 
          go to 900 
      elseif (.not. onto6)  then 
          go to 400 
      endif 
                                                                        
!  To reach here, onto6 must be .true.                                  
                                                                        
!*** logic for 6th order ****                                           
                                                                        
      if (iprint .eq. 1) write(6,905) 
                                                                        
!  Save the 4th order solution on this mesh in uold.                    
      call matcop(nudim, ncomp, ncomp, nmsh, u, uold) 
!  Copy the current mesh into the xxold array.                          
      nmold = nmsh 
      call dcopy(nmold, xx, 1, xxold, 1) 
                                                                        
                                                                        
      iorder = 6 
                                                                        
      if (linear) then 
         call lineq( ncomp, nmsh, nlbc,                                 &
     &    ludone, xx, nudim, u, def,                                    &
     &    delu, rhs, fval,                                              &
     &    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,                     &
     &    ajac, topblk, botblk, bhold, chold, ipvblk,                   &
     &    fsub, dfsub, gsub, dgsub, iflnwt,rpar,ipar)                   
                                                                        
      else 
                                                                        
         call fixjac(ncomp, nmsh, nlbc,                                 &
     &     iorder, ntol, ltol, tol,                                     &
     &     xx, nudim, u, def, def, delu,                                &
     &     rhs, fval, utrial, rhstri,                                   &
     &     rnsq, uint, ftmp, tmprhs,                                    &
     &     ajac, topblk, botblk, ipvblk,                                &
     &     fsub, gsub, iflnwt,rpar,ipar)                                
                                                                        
!  If the fixed Jacobian iterations fail but rnsq is small,             
!  try a Newton procedure.  Set rhsgiv to indicate that                 
!  the right-hand side and fval have already been evaluated             
!  at the given u.                                                      
                                                                        
         if (iflnwt .eq. -3 .and. rnsq .lt. fxfct*epsmch) then 
            rhsgiv = .true. 
                                                                        
            call newteq(ncomp, nmsh, nlbc,                              &
     &           rhsgiv, ntol, ltol, tol,                               &
     &           xx, nudim, u, def,                                     &
     &           delu, rhs, fval,                                       &
     &           utrial, rhstri,                                        &
     &           uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,      &
     &           ajac, topblk, botblk, bhold, chold, ipvblk,            &
     &          fsub, dfsub, gsub, dgsub, iter, iflnwt,rpar,ipar,frscal)
                                                                        
         endif 
      endif 
                                                                        
      if (iflnwt .eq. 0) then 
                                                                        
         call conv6(ncomp, nmsh, ntol, ltol, tol,                       &
     &             nudim, u, uold, etest6, err6,                        &
     &             trst6, onto8, reaft6, succes)                        
      else 
                                                                        
         onto8 = .false. 
                                                                        
         call fail6( ncomp, nmsh, nlbc, ntol, ltol, tol,                &
     &              nfxpnt, fixpnt, iorder, nmax,                       &
     &              xx, nudim, u, rhs, usave, xxold, uold, nmold,       &
     &              ihcomp, irefin,                                     &
     &              rerr, ermx, ratdc,                                  &
     &              reaft6, ddouble , succes, maxmsh,                   &
     &              numbig, nummed,                                     &
     &              wrkrhs,amg, stab_cond,ckappa1,gamma1,ckappa,        &
     &    stiff_cond,      itcond, itcondmax)                           
                                                                        
                                                                        
      endif 
                                                                        
      if (succes) then 
         if (comp_c) then 
           call ESTIMKAPPA(nmsh,ncomp,N,xx,topblk,                      &
     &            nlbc,ncomp,ajac, ncomp,2*ncomp,                       &
     &          ninter,botblk,ncomp-nlbc,ipvblk,isign,c1,wrkrhs,ckappa) 
         end if 
         if (iprint .ne. -1 .and. comp_c ) then 
           if ( ckappa .ge. tolmin/epsmch) write(6,1005) 
         end if 
         iflbvp = 0 
         return 
       elseif (maxmsh) then 
         go to 900 
       elseif (.not. onto8) then 
         go to 400 
       endif 
                                                                        
!**** logic for trying to calculate 8th order solution *****            
                                                                        
      if (iprint .eq. 1) write(6,906) 
                                                                        
      call matcop(nudim, ncomp, ncomp, nmsh, u, uold) 
!  Copy the current mesh into the xxold array.                          
      nmold = nmsh 
      call dcopy(nmold, xx, 1, xxold, 1) 
                                                                        
                                                                        
!  Save the old deferred correction vector def in def6.                 
                                                                        
      call matcop(ncomp, ncomp, ncomp, nmsh-1, def, def6) 
                                                                        
!  For linear problems, calculate the fval array for the                
!  new solution u.                                                      
                                                                        
      if (linear) call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,    &
     &       rpar,ipar)                                                 
                                                                        
!  Calculate 8th order deferred corrections (the array def8).           
                                                                        
      call df8cal (ncomp, nmsh, xx, nudim, u, fval, def8,               &
     &      tmp, fsub,rpar,ipar)                                        
                                                                        
!  For linear problems, the def array is the def8 array.                
!  For nonlinear problems, add the def8 array to the                    
!  already-calculated def array.                                        
                                                                        
      if (linear) then 
         call matcop(ncomp, ncomp, ncomp, nmsh-1, def8, def) 
      else 
         call maxpy(ncomp, nmsh-1, one, def8, ncomp, def) 
      endif 
                                                                        
      iorder = 8 
                                                                        
      if (linear) then 
         call lineq( ncomp, nmsh, nlbc,                                 &
     &    ludone, xx, nudim, u, def,                                    &
     &    delu, rhs, fval,                                              &
     &    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,                     &
     &    ajac, topblk, botblk, bhold, chold, ipvblk,                   &
     &    fsub, dfsub, gsub, dgsub, iflnwt,rpar,ipar)                   
                                                                        
      else 
                                                                        
         call fixjac(ncomp, nmsh, nlbc,                                 &
     &     iorder, ntol, ltol, tol,                                     &
     &     xx, nudim, u, def, def8, delu,                               &
     &     rhs, fval, utrial, rhstri,                                   &
     &     rnsq, uint, ftmp, tmprhs,                                    &
     &     ajac, topblk, botblk, ipvblk,                                &
     &     fsub, gsub, iflnwt,rpar,ipar)                                
                                                                        
!  If the fixed Jacobian iterations fail but rnsq is small,             
!  try a Newton procedure.  Set rhsgiv to indicate that                 
!  the right-hand side and fval have already been evaluated             
!  at the given u.                                                      
                                                                        
         if (iflnwt .eq. -3 .and. rnsq .lt. fxfct*epsmch) then 
            rhsgiv = .true. 
                                                                        
            call newteq(ncomp, nmsh, nlbc,                              &
     &           rhsgiv, ntol, ltol, tol,                               &
     &           xx, nudim, u, def,                                     &
     &           delu, rhs, fval,                                       &
     &           utrial, rhstri,                                        &
     &           uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,      &
     &           ajac, topblk, botblk, bhold, chold, ipvblk,            &
     &          fsub, dfsub, gsub, dgsub, iter, iflnwt,rpar,ipar,frscal)
                                                                        
                                                                        
         endif 
      endif 
                                                                        
      if (iflnwt .eq. 0) then 
                                                                        
         call conv8( ncomp, nmsh, ntol, ltol, tol,                      &
     &              nfxpnt, fixpnt, linear, nmax,                       &
     &              xx, nudim, u, def, def6, def8, uold,                &
     &              ihcomp, irefin, ermx, err6,                         &
     &              etest8, strctr,                                     &
     &              ddouble , nmold, xxold, maxmsh, succes, first8,     &
     &              wrkrhs,amg, stab_cond,ckappa1,gamma1,ckappa,        &
     &     stiff_cond,rpar,ipar)                                        
                                                                        
      else 
                                                                        
         succes = .false. 
         call  fail8(ncomp, nmsh, nfxpnt, fixpnt, nmax,                 &
     &      ntol, ltol, tol, nmold,                                     &
     &      xx, nudim, u, def6, xxold, uold,                            &
     &      ihcomp, irefin, ermx, ddouble , maxmsh,                     &
     &       wrkrhs,amg, stiff_cond, stab_cond)                         
                                                                        
      endif 
      if (maxmsh) then 
         go to 900 
      elseif (.not. succes) then 
         go to 400 
      endif 
                                                                        
!  Successful termination.                                              
       if (comp_c) then 
           call ESTIMKAPPA(nmsh,ncomp,N,xx,topblk,                      &
     &            nlbc,ncomp,ajac, ncomp,2*ncomp,                       &
     &          ninter,botblk,ncomp-nlbc,ipvblk,isign,c1,wrkrhs,ckappa) 
        end if 
      if (iprint .ne. -1 .and. comp_c ) then 
        if ( ckappa .ge. tolmin/epsmch) write(6,1005) 
      end if 
      iflbvp = 0 
                                                                        
      return 
                                                                        
  900 continue 
                                                                        
! Error exit---too many mesh points.                                    
                                                                        
      iflbvp = 1 
      write(6,*) 'Terminated, too many mesh points' 
      if (iprint .ne. -1 .and. comp_c ) then 
        if (linear .and. ckappa .ge. tolmin/epsmch) write(6,1006) 
        if (.not.linear .and. ckappa .ge. tolmin/epsmch) write(6,1007) 
      end if 
      return 
                                                                        
  901 format(1h ,'epsmch',1pe10.3) 
  902 format(1h ,'initial mesh') 
  903 format(1h ,'start 4th order, nmsh',i5) 
  904 format(1h ,'do not go on to 6th') 
  905 format(1h ,'start 6th order') 
  906 format(1h ,'start 8th order') 
 1001 format(1h ,'stiffness = ',1pe11.3) 
 1002 format(1h ,'gamma1    = ',1pe11.3) 
 1003 format(1h ,'kappa1    = ',1pe11.3) 
 1004 format(1h ,'kappa     = ',1pe11.3) 
 1005 format(1h ,'The problem is ill-conditioned, ',                    &
     &     ' the solution could be inaccurate')                         
 1006 format(1h ,'The problem is ill-conditioned,',                     &
     &     ' try with a less stringent tolerance')                      
 1007 format(1h ,'The problem is ill-conditioned,',                     &
     &     ' try with a less stringent tolerance',                      &
     &     ' or with a different initial guess' )                       
      END subroutine                                          
                                                                        
                                                                        
                                   
                                                                        
        SUBROUTINE CONDESTIM(ALEFT,ARIGHT,NMSH,NCOMP,N,XX,TOPBLK,       &
     &            NRWTOP,NOVRLP,ARRAY,                                  &
     &          NRWBLK,NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,ISIGN,OMG,     &
     &          C1,WORK,KPPA,GAMMA)                                     
                                                                        
!     **************************************************************    
!                                                                       
!     COMPUTES THE FIRST AND LAST BLOCK COLUMN OF INVERSE MATRIX AND    
!     ESTIMATE THE CONDITION NUMBERS KPPA, GAMMA                        
!                                                                       
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
        IMPLICIT NONE 
        REAL(DP) TOPBLK,ARRAY,BOTBLK,WORK,C1,                   &
     &          OMG,GAMMA,gamma1,KPPA,MINMG                             
        REAL(DP) ALEFT,ARIGHT,XX,  CSUM 
        INTEGER ISIGN(*) 
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,IPVCD 
        INTEGER NCOMP,NMSH,idmx,idmn,idomg,job 
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),OMG(*),       &
     &          BOTBLK(NRWBOT,*),WORK(*),C1(ncomp,*),XX(*),             &
     &          IPVCD(*)                                                
        INTEGER k,i,j,l 
                                                                        
                                                                        
                                                                        
        MINMG=0.0_dp 
        KPPA=0.0_dp 
        do 300 k=1,ncomp 
           do 330 l=1,N 
              WORK(l)=0.0d0 
              if (k.le.NRWTOP) then 
                   WORK(k)=1.0d0 
              else 
                   WORK(N-ncomp+k)=1.0d0 
              endif 
  330   continue 
                                                                        
        CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,                  &
     &                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,               &
     &                 WORK,0)                                          
                                                                        
                                                                        
           do 340 l=1,N 
             C1(k,l)=WORK(l) 
  340      continue 
                                                                        
  300           continue 
                                                                        
!                                                                       
!     ESTIMATION OF CONDITION NUMBER BY BRUGNANO,TRIGIANTE & MAZZIA     
!                                                                       
!     infinity-norm                                                     
!                                                                       
        l=1 
        do 400 j=1,N-ncomp+1,ncomp 
           OMG(l) = 0d0 
           do  i = j,j+ncomp-1 
               CSUM = 0d0 
               do k=1,ncomp 
                  CSUM = CSUM + ABS(C1(k,i)) 
               end do 
               OMG(l) = MAX( OMG(l),CSUM) 
           end do 
                                                                        
                                                                        
           l=l+1 
  400                 continue 
!                                                                       
!     1-norm                                                            
!                                                                       
!          l=1                                                          
!        do 400 j=1,N-ncomp+1,ncomp                                     
!           OMG(l) = 0d0                                                
!           do  k = 1,ncomp                                             
!               CSUM = 0d0                                              
!               do i=j,j+ncomp-1                                        
!                  CSUM = CSUM + ABS(C1(k,i))                          
!               end do                                                  
!               OMG(l) = MAX( OMG(l),CSUM)                            
!           end do                                                      
!                                                                       
!                                                                       
!           l=l+1                                                       
! 400                  continue                                         
                                                                        
                                                                        
        GAMMA=0.0D0 
        do 450 i=2,NBLOKS+1 
           if (OMG(I).GT.OMG(I-1)) then 
               gamma1=OMG(I)* (XX(I)-XX(I-1)) 
           else 
               gamma1=OMG(I-1)* (XX(I)-XX(I-1)) 
           end if 
           GAMMA=GAMMA + gamma1 
                                                                        
  450                 continue 
        IDOMG=IDAMAX(NBLOKS+1,OMG,1) 
                                                                        
        GAMMA=GAMMA/(ARIGHT-ALEFT) 
                                                                        
        idmx=1 
        idmn=1 
        MINMG = OMG(1) 
        KPPA = OMG(1) 
                                                                        
        do  500 j=2,NBLOKS+1 
           if (OMG(J).GT.KPPA) then 
               KPPA=OMG(J) 
               idmx=j 
           endif 
           if (OMG(J).LT.MINMG) then 
               MINMG=OMG(J) 
               idmn=j 
           endif 
  500          continue 
                                                                        
                                                                        
                                                                        
         RETURN 
      END subroutine
                                                                        
        SUBROUTINE ESTIMKAPPA(NMSH,NCOMP,N,XX,TOPBLK,                   &
     &            NRWTOP,NOVRLP,ARRAY, NRWBLK,                          &
     &          NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,ISIGN,C1,WORK,ckappa) 
                                                                        
!     **************************************************************    
!                                                                       
!                                                                       
!     ESTIMATE THE CONDITION NUMBER  kappa                              
!                                                                       
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
        IMPLICIT NONE 
        REAL(DP) TOPBLK,ARRAY,BOTBLK,WORK,C1 
        REAL(DP) XX 
        INTEGER ISIGN(*) 
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,IPVCD 
        INTEGER NCOMP,NMSH,job 
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),              &
     &          BOTBLK(NRWBOT,*),WORK(*),XX(*),C1(NCOMP,*),             &
     &          IPVCD(*)                                                
                                                                        
                                                                        
        real(dp) ckappa 
        INTEGER KASE,ISOLVE 
        INTEGER k,i,j,l,jj 
                                                                        
                                                                        
!                                                                       
!       DO THE CONDITION NUMBER ESTIMATION BY HIGHAM:                   
!       (DONE IN COLROW) infinity-norm                                  
!                                                                       
                                                                        
        ISOLVE = 0 
        KASE = 0 
   55             CALL DONEST(N,C1,WORK,ISIGN,ckappa,KASE) 
                                                                        
        IF (KASE .NE. 0) THEN 
           ISOLVE = ISOLVE+1 
                                                                        
           IF (KASE .eq. 1) THEN 
               JOB = 1 
           ELSE 
               JOB = 0 
           END IF 
                                                                        
           IF (JOB .eq. 0) THEN 
             DO i=1,NBLOKS 
               DO j=1,NRWBLK 
                  WORK( (i-1)*NRWBLK+NRWTOP+j ) =                       &
     &             (XX(i+1)-XX(i))*WORK( (i-1)*NRWBLK+NRWTOP+j)         
               ENDDO 
             END DO 
           END IF 
           CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,               &
     &                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,               &
     &                 WORK,JOB)                                        
           IF (JOB .eq. 1) THEN 
             DO i=1,NBLOKS 
               DO j=1,NRWBLK 
                  WORK( (i-1)*NRWBLK+NRWTOP+j ) =                       &
     &             (XX(i+1)-XX(i))*WORK( (i-1)*NRWBLK+NRWTOP+j)         
               ENDDO 
             END DO 
           END IF 
           GOTO 55 
        END IF 
                                                                        
                                                                        
      RETURN 
      end subroutine                                        
                                                                        
                                                                        
      subroutine conv4( ncomp, nmsh, ntol, ltol, tol, linear, nmax,     &
     &             xx, nudim, u, defexp, defimp, def, fval,             &
     &             tmp, bhold, chold, dsq, dexr, usave,                 &
     &             ratdc, rerr, ipivot, nmold, xxold,                   &
     &             smooth, reaft6, onto6, strctr, trst6, ddouble ,      &
     &             fsub, maxmsh, succes, first4,                        &
     &             amg,stab_cond,ckappa,stiff_cond,r4,                  &
     &              nfxpnt, fixpnt, irefin,rpar,ipar)                   
                                                                        
!module
      use algprs
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension rpar(*), ipar(*) 
      dimension ltol(ntol), ipivot(*) 
      dimension xx(*), u(nudim,*), tol(ntol) 
      dimension defexp(ncomp,*), defimp(ncomp,*), def(ncomp,*) 
      dimension fval(ncomp,*), tmp(ncomp,4) 
      dimension bhold(ncomp,ncomp,*), chold(ncomp, ncomp,*) 
      dimension dsq(ncomp,ncomp), dexr(ncomp), usave(ncomp,*) 
      dimension xxold(*) 
      dimension ratdc(*), rerr(ncomp,*), amg(*) 
      dimension fixpnt(*), irefin(*), r4(*) 
                                                                        
      logical linear, smooth, reaft6, onto6, strctr, trst6 
      logical ddouble , succes, maxmsh, first4 
      logical stab_cond, stiff_cond 
!      external fsub 
                                                                        
      logical callrt, oscchk, reposs, savedu 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
      common/mchprs/flmin, flmax, epsmch 
                                                                        
      parameter (hundth = 1.0d+5, rerfct = 1.5_dp) 
      parameter (power = 1.0_dp/6.0_dp, one = 1.0_dp) 
      parameter (zero = 0.0_dp, huge = 1._dp/epsilon(1._dp)) 
                                                                        
      intrinsic abs, max 
                                                                        
      save  dfold, oldrt1, savedu, reposs 
                                                                        
      logical adjrer 
!interface
      include 'sepbvp.h'

!  The Newton iteration converged for a 4th order solution.             

                                                                       
      if (iprint .eq. 1) write(6,901) 
                                                                        
      if (first4) then 
         dfold = zero 
         oldrt1 = huge 
         savedu = .false. 
         reposs = .false. 
         first4 = .false. 
      endif 
                                                                        
      succes = .false. 
      maxmsh = .false. 
                                                                        
                                                                        
!  Compute the explicit deferred correction array defexp.               
!  The parameter fval is an ncomp by nmsh array, calculated             
!  by a previous call to fneval for this mesh and solution.             
                                                                        
      call dfexcl( ncomp, nmsh, xx, nudim, u, defexp, fval,             &
     &                tmp, fsub,rpar,ipar )                             
                                                                        
      call matcop(ncomp, ncomp, ncomp, nmsh-1, defexp, def) 
                                                                        
!  If the problem has been declared as smooth, or we might wish         
!  to perform Richardson extrapolation after trying to calculate a      
!  6th order solution, just go on to the 6th order logic (don't         
!  bother to calculate implicit deferred corrections).                  
                                                                        
      if ( smooth .or. reaft6) then 
         if (smooth .and. pdebug) write(6,902) 
         if (reaft6 .and. pdebug) write(6,903) 
         onto6 = .true. 
         return 
      endif 
                                                                        
!  Compute the cheap implicit deferred correction array defimp.         
!  The array chold must be unchanged from a call to jaccal              
!  for this mesh and solution.                                          
!  The temporary arrays dsq and dexr are calculated inside dfimcl.      
                                                                        
      call dfimcl( ncomp, nmsh, defexp, chold, dsq, dexr,               &
     &                 ipivot, defimp)                                  
                                                                        
                                                                        
!  Call dccal to calculate: dfexmx, the maximum-magnitude element       
!  of defexp in components for which tolerances are specified;          
!  incmp, inmsh, and intol, the indices of the component,               
!  mesh interval, and tolerance of dfexmx; derivm, the                  
!  maximum-magnitude element of fval(incmp,*) for all mesh              
!  points; dfimmx, the maximum-magnitude element of defimp in           
!  component incmp; the ratios rat1 and rat2; and the array             
!  ratdc (used in osc).                                                 
      dfctol = hundth*epsmch 
      call dccal( ncomp, nmsh, ntol, ltol,                              &
     &                  defexp, defimp, dfctol, fval,                   &
     &                  ratdc, dfexmx, incmp, inmsh, intol,             &
     &                  derivm, dfimmx, rat1, rat2)                     
                                                                        
!  decid4 sets logical flags that determine the next step of the        
!  algorithm.                                                           
                                                                        
      call decid4(linear, rat1, rat2, dfexmx, dfimmx,                   &
     &     derivm, dfold, tol(intol), oldrt1,                           &
     &     onto6, smooth, callrt, strctr, oscchk, ddouble , reposs)     
                                                                        
                                                                        
      if (pdebug) then 
         if (smooth) write(6,904) 
         if (callrt) write(6,905) 
         if (oscchk) write(6,906) 
         if (strctr) write(6,907) 
         if (reposs) write(6,908) 
         if (savedu) write(6,909) 
         if (onto6) write(6,910) 
         write(6,911) rat1, oldrt1 
      endif 
      oldrt1 = rat1 
      dfold = dfexmx 
      if (callrt) then 
                                                                        
!  ratcor calculates a more expensive implicit deferred correction.     
!  The matrix bhold must be unchanged from the last call to jaccal.     
!  If callrt is true, onto6 is always true also.                        
                                                                        
         call ratcor( ncomp, nmsh, xx, defimp, bhold, def) 
                                                                        
      elseif (linear) then 
         if (oscchk) then 
            call osc (ncomp, nmsh, dfexmx, incmp,                       &
     &            defexp, ratdc, ddouble , inmsh, onto6, trst6,         &
     &            smooth)                                               
                                                                        
         elseif (reposs)  then 
                                                                        
!  If reposs (`Richardson extrapolation possible') is true              
!  for two successive 4th order solutions, a special termination        
!  test may be used based on Richardson extrapolation.                  
                                                                        
!  If reposs is true for the first time, savedu will                    
!  be false; savedu can be true only when reposs is true for a          
!  second consecutive mesh (a doubled version of the first).            
                                                                        
            if (savedu) then 
                                                                        
!  The flag savedu is .true. when the immediately preceding             
!  converged 4th order solution was saved in the array usave,           
!  and the mesh was then doubled.                                       
!  In this case, the routine rerrvl is called to compute a              
!  Richardson extrapolation (RE) error estimate remax.                  
!  The rerr array does not need to be adjusted, so adjrer is false.     
                                                                        
               adjrer = .false. 
               call rerrvl( ncomp, nmsh, nudim, u, usave, ntol,         &
     &              ltol, rerr, remax, itlmx, adjrer )                  
                                                                        
               if (remax .lt. rerfct*tol(itlmx)) then 
                  succes = .true. 
                  return 
               endif 
!           end of logic for savedu = .true.                            
            endif 
                                                                        
!  Here, reposs is .true., but either we hadn't saved the previous      
!  solution, or else the RE error estimate is not small.                
!  Save u in usave, and set savedu to .true.                            
!  Set ddouble  to .true. to indicate that the mesh should be doubled.  
!  Set .onto6. to .false. to return to the top of the 4th order         
!  iteration loop.                                                      
                                                                        
      call matcop(nudim, ncomp, ncomp, nmsh, u, usave) 
            ddouble= .true. 
            savedu = .true. 
            onto6 = .false. 
!        end of logic for reposs = .true.                               
         endif 
!     end of logic for linear                                           
      endif 
                                                                        
      if (pdebug .and. reposs .and. .not. onto6) write(6,912) 
                                                                        
!      if (stiff_cond .and. stab_cond) onto6 = .true.                   
                                                                        
      if (.not. onto6) then 
                                                                        
!  NB: onto6 can be false only for linear problems                      
                                                                        
                                                                        
         if (ddouble ) then 
           if ((use_c) ) then 
             if (stiff_cond .and. .not. stab_cond) then 
               call selcondmsh(ncomp, nmsh,                             &
     &            nfxpnt, fixpnt,  nmax, xx,  irefin,                   &
     &            nmold, xxold, ddouble , maxmsh,r4,amg)                
                  ddouble= .false. 
              else 
                call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
              end if 
            else 
               call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
           end if 
         else 
                                                                        
                                                                        
                                                                        
!  Refine the mesh near interval inmsh; the number of points            
!  to be added depends on the relative size of dfexmx compared          
!  to u and tol.                                                        
                                                                        
!  If the conditining parameter are not stabilized and the problem      
!  is considered stiff  the points are added and removed using          
!  the monitor function given by the conditioning parameters            
!  otherwise  the old technique and the conditioning parameters         
!  are used together.                                                   
!  new functions: selcondmsh, smpselcondmsh                             
!                                                                       
        if (use_c) then 
         if ((stiff_cond)) then 
             drat = dfexmx/                                             &
     &               (max(one, abs(u(incmp,inmsh)))*tol(intol))         
                                                                        
            if (pdebug) write(6,913) drat, u(incmp,inmsh), tol(intol) 
            numadd = drat**power 
            call smpselcondmsh(ncomp, nmsh,                             &
     &        nfxpnt, fixpnt,  nmax, xx,  irefin,inmsh,numadd,          &
     &        nmold, xxold, ddouble , maxmsh,r4,amg)                    
         else 
             drat = dfexmx/                                             &
     &               (max(one, abs(u(incmp,inmsh)))*tol(intol))         
           if (pdebug) write(6,913) drat, u(incmp,inmsh), tol(intol) 
           numadd = drat**power 
             call smpmsh (nmsh, nmax, xx, inmsh, numadd,                &
     &             nmold, xxold, maxmsh)                                
                                                                        
         endif 
        else 
            drat = dfexmx/                                              &
     &               (max(one, abs(u(incmp,inmsh)))*tol(intol))         
            if (pdebug) write(6,913) drat, u(incmp,inmsh), tol(intol) 
            numadd = drat**power 
             call smpmsh (nmsh, nmax, xx, inmsh, numadd,                &
     &             nmold, xxold, maxmsh)                                
                                                                        
         end if 
      end if 
      if (.not. maxmsh)  call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar) 
                                                                        
!     end of logic for .not. onto6                                      
      endif 
                                                                        
                                                                        
                                                                        
      return 
                                                                        
  901 format(1h ,'conv4') 
  902 format(1h ,'smooth') 
  903 format(1h ,'reaft6') 
  904 format(1h ,'smooth') 
  905 format(1h ,'callrt') 
  906 format(1h ,'oscchk') 
  907 format(1h ,'strctr') 
  908 format(1h ,'reposs') 
  909 format(1h ,'savedu') 
  910 format(1h ,'onto6 after decid4') 
  911 format(1h ,'rat1,oldrt1',2(1pe11.3)) 
  912 format(1h ,'reposs and not onto6') 
  913 format(1h ,'drat,u,tol',3(1pe11.3)) 
                                                                        
      end subroutine                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      subroutine fail4( ncomp, nmsh, nlbc, ntol, ltol,                  &
     &             xx, nudim, u, rhs, linear, nmax,                     &
     &             nmold, xxold, uold, tmwork,                          &
     &             iorder, iflnwt, itnwt, ddouble , maxmsh,             &
     &             numbig, nummed,r4,amg,stab_cond,stiff_cond,          &
     &              nfxpnt, fixpnt, irefin,itcond,itcondmax,            &
     &                                                      rpar,ipar)  
                                                                        
 !module
      use algprs                  
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension rpar(*),ipar(*) 
      dimension ltol(ntol) 
      dimension xx(*), u(nudim, *), rhs(*) 
      dimension xxold(*), uold(ncomp, *), tmwork(*) 
      logical linear, ddouble , maxmsh 
      logical stab_cond,stiff_cond 
      dimension amg(*),r4(*), fixpnt(*), irefin(*) 
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
!  The Newton procedure failed to obtain a 4th order solution.          
                                                                        
      if (iprint .eq. 1) write(6,901) 
                                                                        
      maxmsh = .false. 
                                                                        
      if (iflnwt .eq. -1) then 
                                                                        
!  iflnwt = -1 means that the Jacobian was considered singular.         
!  (This is the only possible failure for a linear problem.)            
         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
         call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar) 
      else 
!  The routine mshref decides how to refine the mesh and then           
!  performs the refinement, either by doubling or based on              
!  the rhs vector at the best point.                                    
                                                                        
                                                                        
         call mshref(ncomp, nmsh, nlbc, ntol, ltol,                     &
     &             iorder, rhs, tmwork,                                 &
     &             nmax, xx, nmold, xxold, ddouble , maxmsh,            &
     &             numbig, nummed,amg,stab_cond,stiff_cond,             &
     &             r4, nfxpnt,fixpnt, irefin,itcond,itcondmax)          
                                                                        
                                                                        
         if (.not. maxmsh) then 
              if (linear  .or. itnwt .eq.0) then 
! .or.                                                                  
!     *                            (itnwt.le.2 .and. iflnwt.ne.0)) then 
                call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar) 
              else 
!  Interpolate the partially converged solution.                        
               call matcop(nudim, ncomp, ncomp, nmold, u, uold) 
               call interp(ncomp, nmsh, xx, nudim, u, nmold,            &
     &                 xxold, uold)                                     
              endif 
         endif 
                                                                        
!     End of logic for failure because of some reason other than        
!     singularity.                                                      
      endif 
                                                                        
      return 
  901 format(1h ,'fail4') 
      end subroutine                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      subroutine conv6(ncomp, nmsh, ntol, ltol, tol,                    &
     &             nudim, u, uold, etest6, err6,                        &
     &             trst6, onto8, reaft6, succes)                        

 !module
      use algprs                               
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension ltol(ntol) 
      dimension u(nudim,*), tol(ntol) 
      dimension uold(ncomp,*), etest6(*) 
      logical trst6, onto8, reaft6, succes 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      logical errok 
                                                                        
!  The Newton iteration converged for a 6th-order solution.             
                                                                        
      if (iprint .eq. 1) write(6,901) 
                                                                        
      succes = .false. 
                                                                        
!  The logical flag reaft6 is true only when the 6th order solution     
!  failed.  Since the 6th order solution converged, reaft6 is false.    
                                                                        
      reaft6 = .false. 
      onto8 = .true. 
                                                                        
! Calculate the error estimates for the 4th order solution.             
                                                                        
      call errest (ncomp, nmsh, ntol, ltol, tol,                        &
     &          nudim, u, uold, etest6, err6, errok)                    
                                                                        
      if (trst6 .and. errok) then 
         succes = .true. 
         return 
      endif 
                                                                        
      return 
  901 format(1h ,'conv6') 
      end subroutine                                        
                                                                        
                                                                        
                                                                        
      subroutine fail6( ncomp, nmsh, nlbc, ntol, ltol, tol,             &
     &              nfxpnt, fixpnt, iorder, nmax,                       &
     &              xx, nudim, u, rhs, usave, xxold, uold, nmold,       &
     &              ihcomp, irefin, rerr, ermx, tmwork,                 &
     &              reaft6, ddouble , succes, maxmsh,                   &
     &              numbig, nummed,                                     &
     &             r4,amg, stab_cond,ckappa1,gamma1,ckappa,stiff_cond,  &
     &            itcond, itcondmax)                                    
!module
      use algprs 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension ltol(ntol) 
      dimension fixpnt(*), tol(*), rhs(*) 
      dimension xx(*), u(nudim,*), xxold(*) 
      dimension uold(ncomp,*), usave(ncomp,*) 
      dimension ihcomp(*), irefin(*) 
      dimension rerr(ncomp,*), ermx(*), tmwork(*) 
      logical reaft6, ddouble , succes, maxmsh 
      logical adjrer 
      logical stab_cond, stiff_cond 
      dimension amg(*), r4(*) 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
!  blas: dcopy                                                          
                                                                        
      parameter (eight = 8.0_dp) 
                                                                        
!  nmpt is the standard number of mesh points to be added to selected   
!  intervals when the mesh is altered based on the distribution         
!  in size of the rhs vector.                                           
      integer, parameter :: nmpt = 15
                                                                        
!  Non-convergence of 6th order.                                        
                                                                        
      if (iprint .eq. 1) write(6,901) 
      succes = .false. 
      maxmsh = .false. 
                                                                        
!  NB: the problem must be nonlinear.  Linear problems will either      
!  fail to converge for 4th order, or else, once they've converged      
!  for 4th order, must converge for 6th order.                          
                                                                        
!  Restore the u array to the previous 4th order solution.              
                                                                        
      call matcop(ncomp, nudim, ncomp, nmold, uold, u) 
                                                                        
      if (reaft6 .and. iprint .ge. 0) write(6,9999) 
 9999 format(1h ,'in fail6, reaft6is true') 
      if (.not.reaft6 .and. iprint .ge. 0) write(6,9998) 
 9998 format(1h ,'in fail6, not reaft6') 
      if (ddouble.and. iprint .ge. 0) write(6,9997) 
 9997 format(1h ,'in fail6, ddouble  is true') 
      if (.not.ddouble.and. iprint.ge.0 ) write(6,9996) 
 9996 format(1h ,'in fail6, not double') 
                                                                        
      if (.not. reaft6 .or. .not. ddouble ) then 
                                                                        
!  Here, either                                                         
!  (1) the mesh for which this 6th order solution failed                
!  is not a doubled version of the immediately preceding mesh, or       
!  (2) for the immediately preceding mesh, it is not true               
!  that the 4th order converged and the 6th order failed.               
                                                                        
!  Setting reaft6 to .true. signals that Richardson extrapolation       
!  may be possible if the next 6th order solution fails.  When          
!  reaft6 is true, the routine conv4 immediately sets onto6 to true.    
                                                                        
         reaft6 = .true. 
                                                                        
!  Save the current 4th order solution in usave.                        
                                                                        
         call matcop(nudim, ncomp, ncomp, nmsh, u, usave) 
                                                                        
!  Use the distribution of sizes in the rhs vector to refine the mesh.  
                                                                        
                                                                        
                                                                        
         call mshref(ncomp, nmsh, nlbc, ntol, ltol,                     &
     &              iorder, rhs, tmwork,                                &
     &              nmax, xx, nmold, xxold, ddouble , maxmsh,           &
     &              numbig, nummed,amg,stab_cond,stiff_cond,            &
     &              r4,nfxpnt, fixpnt,irefin,itcond,itcondmax)          
         if (.not. maxmsh) call interp(ncomp, nmsh, xx, nudim, u,       &
     &                        nmold, xxold, uold)                       
                                                                        
      else 
                                                                        
!  Here, reaft6 and ddoubleare both true.  So for two consecutive       
!  meshes, the 4th order converged and the 6th order failed,            
!  and the second mesh is the ddoubleof the first.                      
                                                                        
!  Calculate an error estimate from Richardson extrapolation            
!  with the current and previous 4th order solutions.                   
!  (usave is the 4th order solution saved from the previous (halved)    
!  mesh.)                                                               
!  Set addrer to .true. to signal that the rerr array should            
!  be adjusted.                                                         
                                                                        
         adjrer = .true. 
                                                                        
         call rerrvl( ncomp, nmsh, nudim, u, usave, ntol, ltol,         &
     &          rerr, remax, itlmx, adjrer )                            
         write(6,9994) 
 9994    format(1h ,'***in fail6') 
         write(6,9993) remax, eight*tol(itlmx) 
 9993    format(1h ,'remax',1pe14.4,5x,'8*tol',1pe14.4) 
         if (remax .lt. eight*tol(itlmx)) then 
            succes = .true. 
         else 
                                                                        
!  Richardson extrapolation did not give sufficient accuracy.           
!  Perform selective mesh refinement on the OLD (NB: old!) mesh         
!  and the old (saved) solution, using the error estimate from          
!  Richardson extrapolation to guide where the mesh points are placed.  
!f controllare se posso mettere 2                                       
            nmsh = 1 + (nmsh-1)/2 
            call dcopy(nmsh, xxold, 2, xx, 1) 
            call dcopy(nmsh, amg, 2, tmwork, 1) 
            ipow = 4 
                                                                        
!  The rerr array is overwritten by selmsh.                             
        if (use_c) then 
          if (.not. stiff_cond ) then 
             call selmsh(ncomp, nmsh, ntol, ltol, tol,                  &
     &             nfxpnt, fixpnt, ipow, nmax,                          &
     &             xx, ncomp, usave, rerr, irefin, ihcomp,              &
     &             nmold, xxold, ermx, ddouble , maxmsh)                
                                                                        
           else 
!  The rerr array is overwritten by selconderrmsh.                      
           call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,             &
     &             nfxpnt, fixpnt, ipow, nmax,                          &
     &             xx, ncomp, usave, rerr, irefin, ihcomp,              &
     &             nmold, xxold, ermx, ddouble , maxmsh,                &
     &             r4,tmwork,stab_cond)                                 
                                                                        
          end if 
        else 
              call selmsh(ncomp, nmsh, ntol, ltol, tol,                 &
     &             nfxpnt, fixpnt, ipow, nmax,                          &
     &             xx, ncomp, usave, rerr, irefin, ihcomp,              &
     &             nmold, xxold, ermx, ddouble , maxmsh)                
        end if 
                                                                        
!  If ddoubleis false on exit from selmsh, the call to selmsh has       
!  produced a different (non-doubled) mesh.   Interpolate the           
!  saved solution (from the old mesh) onto the mesh newly created       
!  by selmsh.                                                           
!  NB: Because ddouble is false, we won't try Richardson extrapolation  
!  if the next 4th order converges and 6th order fails.                 
                                                                        
            if (.not. maxmsh) then 
               if (.not. ddouble ) then 
                          call interp(ncomp, nmsh, xx, nudim, u, nmold, &
     &                    xxold, usave)                                 
                else 
                                                                        
!  Selective mesh refinement based on the old mesh simply               
!  produced the same mesh that we started with.  So now refine          
!  starting with the doubled mesh (from before) and the solution.       
                                                                        
                  reaft6 = .true. 
!  Save the solution in usave in case we can carry out Richardson       
!  extrapolation in the same circumstances next time.                   
                  call matcop(nudim, ncomp, ncomp, nmsh, u, usave) 
                                                                        
!  Use the distribution of sizes in the rhs vector to refine the mesh.  
                                                                        
                  call mshref(ncomp, nmsh, nlbc, ntol, ltol,            &
     &                   iorder, rhs, tmwork,                           &
     &                   nmax, xx, nmold, xxold, ddouble , maxmsh,      &
     &                   numbig, nummed,amg,stab_cond,stiff_cond,       &
     &                   r4,nfxpnt, fixpnt,irefin,itcond,itcondmax)     
                                                                        
                  if (.not. maxmsh)                                     &
     &                call interp(ncomp, nmsh, xx, nudim, u,            &
     &                        nmold, xxold, usave)                      
                                                                        
                                                                        
!              end of logic for needing to refine (again) based on the  
!              current mesh                                             
               endif 
                                                                        
!           end of logic for not too many mesh points                   
            endif 
                                                                        
!        end of logic for failure of Richardson extrapolation           
!        to produce a converged solution                                
         endif 
                                                                        
!     end of logic for both reaft6 and ddouble being true               
      endif 
                                                                        
      return 
  901 format(1h ,'fail6') 
      end subroutine                                        
                                                                        
                                                                        
                                                                        
                                                                        
      subroutine conv8( ncomp, nmsh, ntol, ltol, tol,                   &
     &              nfxpnt, fixpnt, linear, nmax,                       &
     &              xx, nudim, u, def, def6, def8, uold,                &
     &              ihcomp, irefin, ermx, err6,                         &
     &              etest8, strctr,                                     &
     &              ddouble , nmold, xxold, maxmsh, succes, first8,     &
     &              r4, amg,stab_cond,ckappa1,gamma1,ckappa,stiff_cond, &
     &    rpar,ipar)                                                    

 !module
      use algprs                                           
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension rpar(*), ipar(*) 
      dimension ltol(ntol), tol(ntol) 
      dimension fixpnt(*) 
      dimension etest8(ntol) 
      dimension xx(*), u(nudim,*), def(ncomp,*) 
      dimension def6(ncomp,*), def8(ncomp,*), uold(ncomp,*) 
      dimension ihcomp(*), irefin(*) 
      dimension ermx(*), xxold(*), amg(*), r4(*) 
      logical linear, strctr, ddouble , maxmsh, succes, first8 
      logical stab_cond, stiff_cond 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      intrinsic max 
                                                                        
      logical errok 
                                                                        
      parameter (one = 1.0_dp, fourth = 0.25_dp, quan8 = 0.025_dp) 
      parameter ( efact  = 100.0_dp, huge = 1._dp/epsilon(1._dp) ) 
                                                                        
!  blas: dload                                                          
                                                                        
      save er6old, er8old 
                                                                        
!  The Newton iteration converged for the 8th order solution.           
                                                                        
                                                                        
      if (iprint .eq. 1) write(6,901) 
                                                                        
      if (first8) then 
         er6old = huge 
         er8old = huge 
         first8 = .false. 
      endif 
                                                                        
      if (.not. linear) then 
         call dload(ntol, one, etest8, 1) 
      else 
         do 10 i = 1, ntol 
            etest8(i) = one/max(quan8, tol(i)**fourth) 
   10    continue 
      endif 
      succes = .false. 
      maxmsh = .false. 
                                                                        
!  Check estimated error.  For a nonlinear problem, all components      
!  of etest8 (the ratios used in testing the error) are set to one.     
!  For a linear problem, the components of etest8 are in general        
!  larger than one.  But if strctr is .true. and the number of mesh     
!  points decreased, we set the elements of etest8 to one (which        
!  makes a stricter test).                                              
                                                                        
      if (linear .and. strctr .and. nmsh .lt. nmold)                    &
     &   call dload(ntol, one, etest8, 1)                               
                                                                        
      call errest (ncomp, nmsh, ntol, ltol, tol,                        &
     &          nudim, u, uold, etest8, err8, errok)                    
      if (errok)  then 
         succes = .true. 
         return 
      endif 
                                                                        
!  At this point, the 8th order solution converged, but did not         
!  satisfy the test for termination.                                    
      if (pdebug) write(6,902) err6, err8, er6old, er8old 
      if (nmsh .lt. nmold .and.                                         &
     &         err6 .gt. efact*er6old .and.                             &
     &         err8 .gt. efact*er8old) then                             
                                                                        
!  If the number of mesh points decreased and the errors in the         
!  6th and 8th order solutions did not decrease sufficiently compared   
!  to the previous mesh, double the mesh and go back to try to          
!  calculate a 4th order solution.                                      
                                                                        
                                                                        
         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
         if (.not. maxmsh) then 
            er6old = err6 
            er8old = err8 
! If the problem is not linear we use the old solution                  
! instead the the first value                                           
! old code                                                              
!            call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar)            
! new code                                                              
            if (linear) then 
               call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar) 
            else 
! we do not use u but uold                                              
               call matcop(nudim, ncomp, ncomp, nmold, u, uold) 
               call interp(ncomp, nmsh, xx, nudim, u,                   &
     &                       nmold, xxold, uold)                        
            endif 
         endif 
         return 
      endif 
                                                                        
!   Here, we know that                                                  
!   (1) the number of mesh points exceeds that for the previous mesh; or
!   (2) the number of mesh points decreased, the 6th and 8th order      
!       errors did not satisfy the termination tests, but they did not  
!       increase so much that the mesh needed to be doubled.            
                                                                        
      er6old = err6 
      er8old = err8 
      if (err8 .le. err6) then 
                                                                        
!  Perform selective mesh refinement based on the 8th order deferred    
!  corrections.  The value of ipow indicates that the error estimate    
!  is of order 6.  Then, for a nonlinear problem, interpolate the       
!  latest solution onto the new mesh.                                   
                                                                        
         ipow = 6 
                                                                        
!  NB: The array def8 will be overwritten by selmsh.                    
                                                                        
        if (use_c) then 
          if (.not. stiff_cond ) then 
            call selmsh(ncomp, nmsh, ntol, ltol, tol,                   &
     &            nfxpnt, fixpnt, ipow, nmax,                           &
     &            xx, ncomp, uold, def8, irefin, ihcomp,                &
     &            nmold, xxold, ermx, ddouble , maxmsh)                 
          else 
!  NB: The array def8 will be overwritten by selconderrmsh.             
            call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,            &
     &            nfxpnt, fixpnt, ipow, nmax,                           &
     &            xx, ncomp, uold, def8, irefin, ihcomp,                &
     &            nmold, xxold, ermx, ddouble , maxmsh,r4,amg,stab_cond)
         end if 
        else 
            call selmsh(ncomp, nmsh, ntol, ltol, tol,                   &
     &            nfxpnt, fixpnt, ipow, nmax,                           &
     &            xx, ncomp, uold, def8, irefin, ihcomp,                &
     &            nmold, xxold, ermx, ddouble , maxmsh)                 
        end if 
                                                                        
                                                                        
                                                                        
         if (.not. maxmsh) then 
            if (linear) then 
               call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar) 
            else 
               call matcop(nudim, ncomp, ncomp, nmold, u, uold) 
               call interp(ncomp, nmsh, xx, nudim, u,                   &
     &                       nmold, xxold, uold)                        
            endif 
         endif 
      else 
                                                                        
!  err8 is greater than err6                                            
                                                                        
!  For a linear problem, set all elements of etest8 to one,             
!  which makes the error test stricter.  (The elements of etest8        
!  may have already been set to one earlier in this routine.)           
                                                                        
         if (linear) call dload(ntol, one, etest8, 1) 
                                                                        
!  Selectively refine the mesh using the old solution and the           
!  6th order deferred correction.  Then, for a nonlinear prpblem,       
!  interpolate the old solution onto the new mesh.                      
                                                                        
         ipow = 4 
                                                                        
!  The array def6 will be overwritten by selmsh.                        
      if (use_c) then 
         if (.not. stiff_cond ) then 
          call selmsh(ncomp, nmsh, ntol, ltol, tol,                     &
     &            nfxpnt, fixpnt, ipow, nmax,                           &
     &            xx, ncomp, uold, def6, irefin, ihcomp,                &
     &            nmold, xxold, ermx, ddouble , maxmsh)                 
           else 
            call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,            &
     &            nfxpnt, fixpnt, ipow, nmax,                           &
     &            xx, ncomp, uold, def6, irefin, ihcomp,                &
     &        nmold, xxold, ermx, ddouble , maxmsh,r4,amg,stab_cond)    
                                                                        
         end if 
        else 
          call selmsh(ncomp, nmsh, ntol, ltol, tol,                     &
     &            nfxpnt, fixpnt, ipow, nmax,                           &
     &            xx, ncomp, uold, def6, irefin, ihcomp,                &
     &            nmold, xxold, ermx, ddouble , maxmsh)                 
        end if 
         if (.not. maxmsh) then 
            if (linear) then 
               call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar) 
            else 
               call interp(ncomp, nmsh, xx, nudim, u,                   &
     &                       nmold, xxold, uold)                        
            endif 
         endif 
!     end of logic for err8 greater than err6                           
      endif 
                                                                        
      if (pdebug .and. .not.succes) write(6,903) 
      return 
                                                                        
  901 format(1h ,'conv8') 
  902 format(1h ,'err6, err8, er6old, er8old',4(1pe11.3)) 
  903 format(1h ,'8th order fails error tests.') 
      end subroutine                                        
                                                                        
                                                                        
                                                                        
                                                                        
      subroutine fail8(ncomp, nmsh, nfxpnt, fixpnt, nmax,               &
     &      ntol, ltol, tol, nmold,                                     &
     &      xx, nudim, u, def6, xxold, uold,                            &
     &      ihcomp, irefin, ermx, ddouble , maxmsh,                     &
     &       r4,amg, stiff_cond, stab_cond)                             
                                                                        
 !module
      use algprs
      implicit integer (i-n)
      implicit real(dp)(a-h,o-z) 
      dimension fixpnt(*), ltol(ntol), tol(ntol) 
      dimension xx(*), u(nudim,*), def6(ncomp,*) 
      dimension xxold(*), uold(ncomp,*) 
      dimension ihcomp(*), irefin(*) 
      dimension ermx(*) 
      dimension amg(*), r4(*) 
      logical ddouble , maxmsh 
      logical stiff_cond, stab_cond 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      if (pdebug) write(6,901) 
                                                                        
!  8th order solution did not converge (the problem must be nonlinear)  
                                                                        
      ipow = 4 
                                                                        
!  Selectively refine the mesh based on the 6th order deferred          
!  correction and the old solution.                                     
                                                                        
!  The def6 array is overwritten by selmsh.                             
                                                                        
        if ( use_c .and. stiff_cond) then 
         call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,               &
     &            nfxpnt, fixpnt, ipow, nmax,                           &
     &            xx, ncomp, uold, def6, irefin, ihcomp,                &
     &      nmold, xxold, ermx, ddouble , maxmsh,r4,amg,stab_cond)      
        else 
         call selmsh(ncomp, nmsh, ntol, ltol, tol,                      &
     &        nfxpnt, fixpnt, ipow, nmax,                               &
     &        xx, ncomp, uold, def6, irefin, ihcomp,                    &
     &        nmold, xxold, ermx, ddouble , maxmsh)                     
        end if 
                                                                        
!  Interpolate to obtain the new initial solution.                      
                                                                        
      if (.not. maxmsh) then 
         call interp(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold) 
      endif 
                                                                        
      return 
  901 format(1h ,'fail8') 
      end subroutine                                        
                                                                        
      subroutine dccal( ncomp, nmsh, ntol, ltol,                        &
     &                     defexp, defimp, dfctol, fval,                &
     &                     ratdc, dfexmx, incmp, inmsh, intol,          &
     &                     derivm, dfimmx, rat1, rat2)                  
 !module
      use algprs               
      implicit integer (i-n)
      implicit real(dp)  (a-h,o-z) 
                                                                        
      dimension  ltol(ntol) 
      dimension  defexp(ncomp,nmsh-1), defimp(ncomp,nmsh-1) 
      dimension  fval(ncomp,nmsh), ratdc(nmsh-1) 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      parameter ( zero = 0.0_dp, one = 1.0_dp ) 
      parameter ( rtst = 50.0_dp, tstrat = 0.1_dp ) 
                                                                        
      intrinsic abs, max 
                                                                        
!  blas: idamax                                                         
                                                                        
!  Find dfexmx, the maximum-magnitude element of defexp                 
!  in components for which a tolerance is specified.                    
!  The component index of dfexmx (incmp), its mesh                      
!  interval index (inmsh), and its tolerance index (intol),             
!  are output parameters.                                               
                                                                        
      dfexmx = zero 
      do 10 it = 1, ntol 
         icmp = ltol(it) 
         idmx = idamax(nmsh-1, defexp(icmp, 1), ncomp) 
         dval = abs(defexp(icmp, idmx)) 
         if (dval .ge. dfexmx) then 
            dfexmx = dval 
            incmp = icmp 
            inmsh = idmx 
            intol = it 
         endif 
   10 continue 
                                                                        
      if (pdebug) then 
         write(6,901) 
         write(6,902) dfexmx, incmp, inmsh, intol 
      endif 
                                                                        
!  Find derivm (maximum-magnitude element of fval(incmp,*))             
!  for all mesh points.                                                 
                                                                        
      idmx = idamax(nmsh, fval(incmp, 1), ncomp) 
      derivm = abs(fval(incmp, idmx)) 
      if (pdebug) write(6,903) derivm 
                                                                        
!  For component incmp, go through the mesh intervals to calculate      
!  (1) dfimmx, the maximum implicit deferred correction;                
!  (2) two crucial ratios, rat1 and rat2, used in deciding whether      
!      to refine the mesh;                                              
!  (3) the array ratdc of deferred-correction ratios (explicit to       
!      implicit).                                                       
                                                                        
!  In defining rat1 and rat2, we consider only intervals for            
!  which the explicit deferred correction (defexp) exceeds the          
!  tolerance dfctol in magnitude.  If it does not, the associated       
!  interval does not affect rat1 or rat2, and the value of ratdc        
!  is taken as 1.                                                       
                                                                        
!  rat2 is the maximum-magnitude ratio of sufficiently large            
!  defexp to the larger of defimp and dfctol.                           
                                                                        
!  rat1 is the maximum-magnitude ratio of sufficiently large            
!  defexp to the larger in magnitude of defimp and dfctol, but only     
!  for those values of defexp greater than tstrat*dfexmx.               
!  Thus by construction rat1 is less than or equal to rat2.             
                                                                        
      rat1 = zero 
      rat2 = zero 
      dfimmx = zero 
      smtest = tstrat*dfexmx 
                                                                        
      do 100 im = 1, nmsh-1 
         texp = defexp(incmp, im) 
         timp = defimp(incmp, im) 
         dfimmx = max(dfimmx, abs(timp)) 
         abtexp = abs(texp) 
         if (abtexp .le. dfctol) then 
            ratdc(im) = one 
         else 
            if (abs(timp) .lt. dfctol) timp = dfctol 
            ratdc(im) = texp/timp 
            abrat = abs(ratdc(im)) 
            rat2 = max(rat2, abrat) 
            if (abs(texp) .ge. smtest                                   &
     &                .and. abrat .ge. rat1)  rat1 = abrat              
         endif 
!         if (pdebug) write(6,904) im, texp, timp, ratdc(im), dfctol    
  100 continue 
                                                                        
      if (pdebug) write(6,905) rat1, rat2, dfimmx 
                                                                        
      return 
                                                                        
  901 format(1h ,'dccal') 
  902 format(1h ,'dfexmx, incmp, inmsh, intol',1pe11.3,3i5) 
  903 format(1h ,'derivm',1pe11.3) 
  904 format(1h ,'im, texp, timp, ratdc, dfctol',i5,4(1pe11.3)) 
  905 format(1h ,'rat1, rat2, dfimmx', 3(1pe11.3)) 
      end subroutine                                        
                                                                        
                                                                        
      subroutine decid4(linear, rat1, rat2, dfexmx, dfimmx,             &
     &     derivm, dfold, tolval, oldrt1,                               &
     &     onto6, smooth, callrt, strctr, oscchk, ddouble , reposs)     
                                                                        
 !module
      use algprs                                                           
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
                                                                        
      logical linear, onto6, smooth, callrt, strctr,                    &
     &     oscchk, ddouble , reposs                                     
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      logical stest 
                                                                        
      parameter ( tenth = 0.1_dp, one = 1.0_dp, two = 2.0_dp ) 
      parameter ( thrtwo = 32.0_dp ) 
      parameter ( rtst = 50.0_dp, derval = 50.0_dp ) 
                                                                        
!  decid4 evaluates information about the deferred corrections          
!  and the nature of the problem, and sets various logical              
!  flags that control the subsequent logic of the algorithm.            
                                                                        
!  This code has been written for clarity, and NOT for efficiency.      
                                                                        
      onto6 = .true. 
      callrt = .false. 
      smooth = .false. 
      oscchk = .false. 
      strctr = .false. 
      reposs = .false. 
      ddouble = .false. 
                                                                        
!  rat2 is always greater than or equal to rat1.                        
                                                                        
      if (pdebug) write(6,901) 
      if (pdebug) write(6,902) tolval, rtst 
                                                                        
      stest = .true. 
      if (linear) stest = dfexmx .lt. tenth*dfold 
                                                                        
                                                                        
                                                                        
      if (rat2 .lt. rtst) then 
         if (stest) then 
            smooth = .true. 
         else 
            oscchk = .true. 
         endif 
                                                                        
         return 
      endif 
                                                                        
!  We know now that rat2 .ge. rtst.                                     
                                                                        
      thttol = thrtwo*tolval 
      if (pdebug) write(6,903) thttol 
                                                                        
      if (rat1 .lt. rtst .and. dfexmx .lt. thttol) then 
          if (stest) then 
             smooth = .true. 
          else 
             oscchk = .true. 
          endif 
                                                                        
          return 
      endif 
                                                                        
      if (rat1 .lt. rtst .and. dfexmx .ge. thttol) then 
         callrt = .true. 
                                                                        
         return 
      endif 
                                                                        
!  We know now that rat1 .ge. rtst (and that rat2 .ge. rtst).           
                                                                        
      if (derivm .gt. derval .and. dfexmx .lt. thttol) then 
         if (stest) then 
            smooth = .true. 
         else 
            oscchk = .true. 
         endif 
                                                                        
         return 
      endif 
                                                                        
      if (derivm .gt. derval .and. dfexmx .gt. thttol) then 
         if (dfimmx .lt. one) then 
            callrt = .true. 
                                                                        
         else 
            strctr = .true. 
            if (linear) then 
               onto6 = .false. 
               if (two*rat1 .ge. oldrt1) ddouble = .true. 
                                                                        
!           end of logic for linear                                     
            endif 
!        end of logic for dfimmx .ge. one                               
         endif 
                                                                        
         return 
!     end of logic for derivm .gt. derval .and dfexmx .gt. thttol       
      endif 
                                                                        
!  To reach this point in the code, both of the following must hold:    
!    rat1 .ge. rtst (which means that rat2 .ge. rtst)                   
!    derivm .le. derval                                                 
                                                                        
!  On linear problems, a special termination rule is tried if two       
!  conditions hold:                                                     
!    (1) the 4th order solution has been computed on two consecutive    
!        meshes, one of which is the double of the other                
!        (this holds when ddouble  is .true.), and                      
!    (2) on both meshes, rat1 .ge. rtst, and derivm is small.  When     
!        the conditions in (2) hold for a particular mesh, decid4       
!        sets reposs to .true. to indicate that Richardson              
!        extrapolation may be possible.                                 
!  This set of tests is to take care of transients kept out by          
!  initial conditions.                                                  
                                                                        
      if (linear) reposs = .true. 
      return 
                                                                        
  901 format(1h ,'decid4') 
  902 format(1h ,'tolval, rtst',2(1pe11.3)) 
  903 format(1h ,'thttol',1pe11.3) 
      end subroutine                                        
                                                                        
                                                                        
      subroutine dfexcl (ncomp, nmsh, xx, nudim, u, defexp, fval,       &
     &               tmp, fsub,rpar,ipar)                               
 !module
      use algprs
      implicit integer (i-n)
      implicit real(dp) (a-h, o-z) 
      integer ncomp, nmsh 
      dimension rpar(*), ipar(*) 
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp, nmsh) 
      dimension  defexp(ncomp,nmsh-1), tmp(ncomp,4) 
!      external fsub 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      parameter ( half = 0.5_dp, fourth = 0.25_dp, thfrth= 0.75_dp ) 
                                                                        
      common /consts/ alp1, alp2, alp3, bet0, bet2, bet3, bet4,         &
     &    a1, b1, c1, d1, e2, f2, c2, d2, b2, a2,                       &
     &    p3, q3, e3, f3, c3, d3, a3, b3, a4, p4, x4, e4, c4,           &
     &    a5, b5, c5, d5, e5, f5, a6, b6, c6                            

!interface 
      include 'sepbvp.h'
                                                                        
!  Given the nmsh mesh points xx, the estimated solution                
!  u and the array fval of function values at (xx(im), u(*,im)),        
!  im = 1,...,nmsh, dfexcl calculates sixth-order explicit              
!  deferred correction, stored in the array defexp, indexed             
!  over the components and mesh intervals.                              
                                                                        
!  The array tmp is workspace for 4 intermediate vectors of             
!  dimension ncomp.                                                     
                                                                        
      do 50 im = 1, nmsh-1 
         hmsh = xx(im+1) - xx(im) 
         do 10 ic = 1, ncomp 
            tmp(ic,1) = (a5*u(ic, im+1) + b5*u(ic, im))                 &
     &         + hmsh*(c5*fval(ic,im)- d5*fval(ic,im+1))                
            tmp(ic,2) = (b5*u(ic,im+1) + a5*u(ic,im))                   &
     &         + hmsh*(-c5*fval(ic,im+1) + d5*fval(ic,im))              
   10    continue 
                                                                        
         call fsub (ncomp, xx(im) + fourth*hmsh, tmp(1,1),              &
     &          tmp(1,3),rpar,ipar)                                     
         call fsub (ncomp, xx(im) + thfrth*hmsh, tmp(1,2),              &
     &          tmp(1,4),rpar,ipar)                                     
                                                                        
         do 20 ic = 1, ncomp 
            tmp(ic,1) = half*(u(ic,im+1) + u(ic,im))                    &
     &          + e5*hmsh*(fval(ic,im+1) - fval(ic,im))                 &
     &          - f5*hmsh*(tmp(ic,4) - tmp(ic,3))                       
   20    continue 
                                                                        
         call fsub(ncomp, half*(xx(im) + xx(im+1)), tmp(1,1),           &
     &          tmp(1,2),rpar,ipar)                                     
         do 30 ic = 1, ncomp 
            defexp(ic,im) = hmsh*(a6*(fval(ic,im+1) + fval(ic,im))      &
     &           + b6*(tmp(ic,3) + tmp(ic,4)) + c6*tmp(ic,2))           &
     &           - u(ic,im+1) + u(ic,im)                                
   30    continue 
                                                                        
   50 continue 
      return 
      end subroutine                                        
                                                                        
                                                                        
      subroutine df8cal (ncomp, nmsh, xx, nudim, u, fval, def8,         &
     &      tmp, fsub,rpar,ipar)                                        
                                                                        
!   Given the mesh points xx, the solution u, and the function          
!   values fval, df8cal computes eighth-order deferred corrections,     
!   which are stored in def8.                                           
!   The array tmp is workspace for 8 intermediate vectors.              
 !module
      use algprs
      implicit integer (i-n)
      implicit real(dp) (a-h, o-z) 
      integer ncomp, nmsh 
      dimension rpar(*),ipar(*) 
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp,nmsh) 
      dimension def8(ncomp, nmsh-1),  tmp(ncomp,8) 
!      external fsub 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      parameter ( half = 0.5_dp, two = 2.0_dp ) 
      parameter ( fc1 = 0.625_dp, fc2= 0.375_dp ) 
                                                                        
      common /consts/ alp1, alp2, alp3, bet0, bet2, bet3, bet4,         &
     &    a1, b1, c1, d1, e2, f2, c2, d2, b2, a2,                       &
     &    p3, q3, e3, f3, c3, d3, a3, b3, a4, p4, x4, e4, c4,           &
     &    a5, b5, c5, d5, e5, f5, a6, b6, c6                            
 
!interface
      include 'sepbvp.h'
                                                                       
      do 100 im = 1, nmsh-1 
                                                                        
         hmsh = xx(im+1) - xx(im) 
                                                                        
         do 10 ic = 1, ncomp 
            tmp(ic,1) = a1*u(ic,im+1) + b1*u(ic,im)                     &
     &         + hmsh*(c1*fval(ic,im+1) + d1*fval(ic,im))               
            tmp(ic,2) = b1*u(ic,im+1) + a1*u(ic,im)                     &
     &         - hmsh*(c1*fval(ic,im) + d1*fval(ic,im+1))               
   10    continue 
                                                                        
         call fsub(ncomp, xx(im) + fc1*hmsh, tmp(1,1), tmp(1,3),        &
     &    rpar,ipar)                                                    
         call fsub(ncomp, xx(im) + fc2*hmsh, tmp(1,2), tmp(1,4),        &
     &    rpar,ipar)                                                    
         do 20 ic = 1, ncomp 
            tmp(ic,1) = a2*u(ic,im+1) + b2*u(ic,im)                     &
     &         + hmsh*(c2*fval(ic,im+1) + d2*fval(ic,im)                &
     &         + e2*tmp(ic,3) + f2*tmp(ic,4))                           
            tmp(ic,2) = b2*u(ic,im+1) + a2*u(ic,im)                     &
     &         - hmsh*(d2*fval(ic,im+1) + c2*fval(ic,im)                &
     &         + f2*tmp(ic,3) + e2*tmp(ic,4))                           
   20    continue 
                                                                        
         call fsub(ncomp, xx(im) + (half + alp2)*hmsh, tmp(1,1),        &
     &               tmp(1,5),rpar,ipar)                                
         call fsub(ncomp, xx(im) + (half - alp2)*hmsh, tmp(1,2),        &
     &               tmp(1,6),rpar,ipar)                                
         do 30 ic = 1, ncomp 
            tmp(ic,1) = a3*u(ic,im+1) + b3*u(ic,im)                     &
     &         + hmsh*(c3*fval(ic,im+1) + d3*fval(ic,im)                &
     &         + e3*tmp(ic,3) + f3*tmp(ic,4)                            &
     &         + p3*tmp(ic,5) + q3*tmp(ic,6))                           
            tmp(ic,2) = b3*u(ic,im+1) + a3*u(ic,im)                     &
     &         - hmsh*(d3*fval(ic,im+1) + c3*fval(ic,im)                &
     &         + f3*tmp(ic,3) + e3*tmp(ic,4)                            &
     &         + q3*tmp(ic,5) + p3*tmp(ic,6))                           
   30    continue 
                                                                        
         call fsub (ncomp, xx(im) + (half + alp3)*hmsh, tmp(1,1),       &
     &                 tmp(1,7),rpar,ipar)                              
         call fsub (ncomp, xx(im) + (half - alp3)*hmsh, tmp(1,2),       &
     &                 tmp(1,8),rpar,ipar)                              
         do 40 ic = 1, ncomp 
            tmp(ic,1) = a4*(u(ic,im+1) + u(ic,im))                      &
     &         + hmsh*(c4*(fval(ic,im+1) - fval(ic,im))                 &
     &         + e4*(tmp(ic,3) - tmp(ic,4))                             &
     &         + x4*(tmp(ic,7) - tmp(ic,8)))                            
   40    continue 
                                                                        
         call fsub (ncomp, xx(im) + half*hmsh, tmp(1,1), tmp(1,2),      &
     &       rpar, ipar)                                                
         do 50 ic = 1, ncomp 
            def8(ic,im) =                                               &
     &          hmsh*(bet0*(fval(ic,im) + fval(ic,im+1))                &
     &          + bet2*(tmp(ic,5) + tmp(ic,6))                          &
     &          + bet3*(tmp(ic,7) + tmp(ic,8))                          &
     &          + two*bet4*tmp(ic,2))                                   &
     &          - u(ic,im+1) + u(ic,im)                                 
   50    continue 
                                                                        
  100 continue 
      return 
      end subroutine                                        
                                                                        
                                                                        
      subroutine dfimcl( ncomp, nmsh, defexp, chold, dsq, dexr,         &
     &                 ipivot, defimp)                                  
!module
      use algprs 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension defexp(ncomp, nmsh-1), chold(ncomp, ncomp, nmsh-1) 
      dimension dsq(ncomp, ncomp), dexr(ncomp) 
      dimension ipivot(ncomp), defimp(ncomp, nmsh-1) 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
!  blas: dcopy                                                          
                                                                        
      parameter (zero = 0.0_dp) 
                                                                        
!  dfimcl calculates the rational deferred correction array,            
!  which is indexed over the components and mesh intervals.             
                                                                        
      call mtload(ncomp, nmsh-1, zero, ncomp, defimp) 
                                                                        
      do 100 im = 1, nmsh-1 
                                                                        
         call dcopy(ncomp, defexp(1,im), 1, dexr(1), 1 ) 
         do 50 ic = 1, ncomp 
            call dcopy(ncomp, chold(1,ic,im), 1, dsq(1,ic), 1) 
   50    continue 
         call lufac (ncomp, ncomp, dsq, ipivot, ierlu) 
         if (ierlu .eq. 0) then 
            call lusol (ncomp, ncomp, dsq, ipivot, dexr,                &
     &          defimp(1, im))                                          
         endif 
  100 continue 
      return 
      end subroutine                                        
                                                                        
                                                                        
      subroutine osc (ncomp, nmsh, dfexmx, incmp,                       &
     &     defcor, ratdc, ddouble , inmsh, onto6, trst6, smooth)        
 !module
      use algprs
      implicit integer (i-n)
      implicit real(dp) (a-h, o-z) 
      dimension defcor(ncomp, *), ratdc(*) 
                                                                        
      logical ddouble , onto6, trst6, smooth 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      intrinsic abs 
                                                                        
      parameter ( zero = 0.0_dp, half = 0.5_dp ) 
      parameter ( frac1 = 0.1_dp, frac2 = 1.0d-2 ) 
                                                                        
!  For linear problems, subroutine osc performs heuristic tests         
!  to detect an oscillating solution.  the tests check whether          
!  (significant) explicit and implicit deferred corrections have        
!  different (componentwise) signs.                                     
                                                                        
!  dfexmx is the maximum-magnitude explicit deferred correction,        
!  and is known to occur in component incmp.                            
                                                                        
!  The array defcor contains the explicit deferred corrections.         
!  The array ratdc contains the ratios of explicit to implicit          
!  deferred corrections.                                                
                                                                        
!  jsndif counts the number of differences in sign.                     
                                                                        
      jsndif = 0 
      rmax = zero 
                                                                        
!  allsum is the sum of the magnitudes of all deferred corrections,     
!  smlsum is the sum of the magnitudes of the small deferred            
!  corrections, and bigsum is the sum of the magnitudes of the          
!  large deferred corrections.  Here, small is defined as less          
!  than half of the maximum.                                            
                                                                        
      ninter = nmsh - 1 
                                                                        
      if (pdebug) write(6,901) 
                                                                        
      allsum = zero 
      smlsum = zero 
      bigsum = zero 
      ibig = 0 
      ism = 0 
                                                                        
      do 30 im = 1, ninter 
         abdef = abs(defcor(incmp,im)) 
         allsum = allsum + abdef 
         if (abdef .lt. half*dfexmx) then 
              ism = ism + 1 
              smlsum = smlsum + abdef 
         else 
              ibig = ibig + 1 
              bigsum = bigsum + abdef 
         endif 
                                                                        
!  The counter of sign differences is incremented if (1) ratdc is negati
!  (which means that the two deferred corrections have opposite         
!  sign) and (2) the explicit deferred correction is not too small      
!  relative to the maximum.                                             
                                                                        
         if (pdebug) write(6,902) im, ratdc(im), abdef, frac2*dfexmx 
                                                                        
         if (ratdc(im).lt.zero .and. abdef.ge.frac2*dfexmx) then 
            jsndif = jsndif + 1 
                                                                        
!  If more than 4 sign differences have occurred, exit after setting    
!  ddouble to .true., which signals that the mesh                       
!  should be doubled (i.e., twice as many intervals).                   
                                                                        
            if (jsndif.gt.4) then 
               onto6 = .false. 
               ddouble = .true. 
               return 
            endif 
            if (abs(ratdc(im)).ge.rmax) then 
               rmax = abs(ratdc(im)) 
               inmsh = im 
            endif 
         endif 
   30 continue 
                                                                        
      if (pdebug) write(6,903) rmax, jsndif 
                                                                        
      avsm = zero 
      if (ism.gt.0) avsm = smlsum/ism 
      avbg = zero 
      if (ibig.gt.0) avbg = bigsum/ibig 
      ave = allsum/ninter 
                                                                        
      if (pdebug) write(6,904) ave, avsm, avbg 
                                                                        
      if (avsm.gt.frac1*avbg .or. ave.gt.half*avbg) then 
                                                                        
!  The error appears to be uniformly large.                             
!  Signal that the 6th order solution should be calculated.             
         onto6 = .true. 
                                                                        
      elseif (jsndif.eq.0) then 
!  If there were no sign changes, the problem appears to be smooth.     
         smooth = .true. 
         onto6 = .true. 
      else 
                                                                        
!  If the sign changed at between 1 and 4 points, don't go on to        
!  6th order, and don't ever accept a 6th order solution even if the    
!  error estimate at a later stage indicates that it is OK to do so.    
!  Set ddouble to .false., to signal that the mesh will not necessarily 
!  be doubled.                                                          
                                                                        
          ddouble = .false. 
          onto6 = .false. 
          trst6 = .false. 
      endif 
      return 
  901 format(1h ,'osc') 
  902 format(1h ,'im, ratdc, abdef, val',i5,3(1pe11.3)) 
  903 format(1h ,'rmax, jsndif', 1pe11.3,i5) 
  904 format(1h ,'ave, avsm, avbg', 3(1pe11.3)) 
      end subroutine                                        
                                                                        
                                                                        
      subroutine ratcor ( ncomp, nmsh, xx, defimp, bhold, dfrat) 
      implicit integer (i-n)
      implicit real(dp) (a-h, o-z) 
      dimension xx(nmsh), defimp(ncomp,nmsh-1) 
      dimension dfrat(ncomp,nmsh-1), bhold(ncomp,ncomp,nmsh-1) 
                                                                        
!  blas: ddot, dscal                                                    
                                                                        
      parameter (zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp) 
                                                                        
      ninter = nmsh - 1 
      do 10 im = 1, ninter 
         hmsh = xx(im+1) - xx(im) 
         call dscal(ncomp*ncomp, (-half*hmsh), bhold(1,1,im), 1) 
   10 continue 
                                                                        
      do 20 im = 1, ninter 
      do 20 ic = 1, ncomp 
         bhold(ic,ic,im) = bhold(ic,ic,im) + one 
   20 continue 
                                                                        
      do 30 im = 1, ninter 
      do 30 ic = 1, ncomp 
         dfrat(ic,im) = ddot(ncomp, bhold(ic,1,im), ncomp,              &
     &                        defimp(1,im), 1)                          
   30 continue 
      return 
      end subroutine                                        
                                                                        
                                                                        
      subroutine stcons 
                                                                        
        implicit integer (i-n)
      implicit real(dp)  (a-h,o-z) 
                                                                        
!  stcons computes constants needed in integration formulae             
!  and stores them in a labeled common area.                            
                                                                        
      parameter ( one = 1.0_dp,  four = 4.0_dp, two = 2.0_dp ) 
      parameter ( five = 5.0_dp, three = 3.0_dp ) 
      parameter ( half = 0.5_dp, fourth = 0.25_dp ) 
                                                                        
      common /consts/ alp1, alp2, alp3, bet0, bet2, bet3, bet4,         &
     &    a1, b1, c1, d1, e2, f2, c2, d2, b2, a2,                       &
     &    p3, q3, e3, f3, c3, d3, a3, b3, a4, p4, x4, e4, c4,           &
     &    a5, b5, c5, d5, e5, f5, a6, b6, c6                            
                                                                        
                                                                        
                                                                        
      alp1 = one/8.0_dp 
      alp1sq = one/64.0_dp 
      alp2sq = five/28.0_dp 
      alp3sq = five/84.0_dp 
      alp2 = sqrt(alp2sq) 
      alp3 = sqrt(alp3sq) 
      bet0 = one/6.0_dp - (10.0_dp - 28.0_dp*(alp2sq + alp3sq))/        &
     &         (105.0_dp*(one - four*alp2sq)*                           &
     &         (one - four*alp3sq))                                     
      bet2 = -(28.0_dp*alp3sq - three)/                                 &
     &         (1680.0_dp*alp2sq*                                       &
     &         (one - four*alp2sq)*(alp2sq - alp3sq))                   
      bet3 = (28.0_dp*alp2sq - three)/                                  &
     &          (1680.0_dp*alp3sq*                                      &
     &          (one - four*alp3sq)*(alp2sq - alp3sq))                  
      bet4 = half - bet0 - bet2-bet3 
      a1 = half*(one - alp1)*(two*alp1 + one)*                          &
     &          (two*alp1 + one)                                        
      b1 = half*(one + alp1)*                                           &
     &          (two*alp1 - one)*(two*alp1 - one)                       
      c1 = half*(alp1sq - fourth)*(two*alp1 + one) 
      d1 = half*(alp1sq - fourth)*(two*alp1 - one) 
      uu = alp2*((four*alp2sq - one)**2)/                               &
     &       ((four*alp1sq - one)*(20.0_dp*alp1*alp1 - one))            
      vv = ((four*alp2sq - one)**2)/                                    &
     &      (16.0_dp*alp1*(four*alp1sq - one))                          
      e2 = half*(uu + vv) 
      f2 = half*(uu - vv) 
      rr = half*(alp2*(four*alp2sq - one) +                             &
     &      (one - 12.0_dp*alp1sq)*(e2 + f2))                           
      ss = fourth*(four*alp2sq - one) - two*alp1*(e2 - f2) 
      c2 = half*(rr + ss) 
      d2 = half*(rr - ss) 
      ww = two*(alp2 - (c2 + d2 + e2 + f2)) 
      b2 = half*(one - ww) 
      a2 = half*(one + ww) 
      z1 = (three - 28.0_dp*alp3sq)/                                    &
     &         (1680.0_dp*alp2*(four*alp2sq - one)*                     &
     &         (alp2*alp2 - alp3sq)*bet3)                               
      z2 = one/(105.0_dp*alp3*bet3*                                     &
     &    (20.0_dp*alp2sq - one)*(four*alp2sq - one))                   
      p3 = half*(z1 + z2) 
      q3 = half*(z2 - z1) 
      u1 = (alp3*((four*alp3sq - one)**2)-                              &
     &       (p3 + q3)*(20.0_dp*alp2sq - one)                           &
     &       *(four*alp2sq - one))/                                     &
     &        ((four*alp1sq - one)*(20.0_dp*alp1sq - one))              
      v1 = (alp3sq*(one - two*alp3sq)-                                  &
     &         two*alp2*(one - four*alp2sq)*(p3 - q3)                   &
     &        -one/8.0_dp)/(two*alp1*(one - four*alp1sq))               
      e3 = half*(u1 + v1) 
      f3 = half*(u1 - v1) 
      r1 = half*(alp3*(four*alp3sq - one) +                             &
     &         (e3 + f3)*(one - 12.0_dp*alp1sq) +                       &
     &         (p3 + q3)*(one - 12.0_dp*alp2sq))                        
      s1 = alp3sq - fourth - two*alp1*(e3 - f3)                         &
     &         - two*alp2*(p3 - q3)                                     
      c3 = half*(r1 + s1) 
      d3 = half*(r1 - s1) 
      w1 = two*(alp3 - (c3 + d3 + e3 + f3 + p3 + q3)) 
      a3 = half*(one + w1) 
      b3 = half*(one - w1) 
      a4 = half 
      p4 = 0.0_dp 
      x4 = (three - 28.0_dp*alp2sq)/                                    &
     &        (3360.0_dp*alp3*bet4*(four*alp3sq - one)                  &
     &        *(alp3sq - alp2sq))                                       
      e4 = (0.125_dp + four*alp2*p4*                                    &
     &         (one - four*alp2sq) +                                    &
     &         four*alp3*x4*(one - four*alp3sq))/                       &
     &         (four*alp1*(four*alp1sq - one))                          
      c4 = -(0.125_dp + two*alp1*e4 + two*alp2*p4 +                     &
     &         two*alp3*x4)                                             
      a5 = five/32.0_dp 
      b5 = 27.0_dp/32.0_dp 
      c5 = 9.0_dp/64.0_dp 
      d5 = three/64.0_dp 
      e5 = five/24.0_dp 
      f5 = two/three 
      a6 = 7.0_dp/90.0_dp 
      b6 = 16.0_dp/45.0_dp 
      c6 = two/15.0_dp 
      return 
      end subroutine                                        
                                                                        
      subroutine fixjac(ncomp, nmsh, nlbc,                              &
     &    iorder, ntol, ltol, tol,                                      &
     &    xx, nudim, u, defcor, defnew, delu,                           &
     &    rhs, fval, utrial, rhstri,                                    &
     &    rnsq, uint, ftmp, tmprhs,                                     &
     &    ajac, topblk, botblk, ipivot,                                 &
     &    fsub, gsub, iflag,rpar,ipar)                                  
                                                                        
! Fixed Jacobian iterations.                                            
 !module
      use algprs                                                           
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension rpar(*), ipar(*) 
      dimension  ltol(ntol), tol(ntol) 
      dimension  xx(nmsh), u(nudim,nmsh), defcor(ncomp,nmsh-1) 
      dimension  defnew(ncomp,nmsh-1), delu(ncomp,nmsh) 
      dimension  rhs(ncomp*nmsh), fval(ncomp,nmsh) 
      dimension  utrial(ncomp,nmsh), rhstri(ncomp*nmsh) 
      dimension  uint(ncomp), ftmp(ncomp), tmprhs(ncomp*nmsh) 
      dimension  ajac(ncomp, 2*ncomp, nmsh-1) 
      dimension  topblk(nlbc,ncomp), botblk(ncomp-nlbc,ncomp) 
      dimension  ipivot(ncomp*nmsh) 
      logical    better 
                                                                        
!      external   fsub, gsub 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
      common/mchprs/flmin, flmax, epsmch 
                                                                        
      intrinsic  abs, max 
                                                                        
!  blas: dcopy, dssq                                                    
                                                                        
      parameter ( one    = 1.0_dp ) 
      parameter ( xlarge = 1.0d+6, huge = 1._dp/epsilon(1._dp) ) 
      parameter ( rngrow = 16.0_dp, rfact = 100.0_dp ) 
      parameter ( tolfct = 0.1_dp ) 
      integer, parameter :: lmtfrz = 8

!interface
      include 'sepbvp.h'
                                                                       
!  The iteration scheme uses a fixed Jacobian matrix to solve for       
!  correction vectors, once there has been convergence of the Newton    
!  iterations on this mesh.   It is assumed that the LU                 
!  factors of the Jacobian have been left unaltered since               
!  their calculation.                                                   
                                                                        
      if (iprint .eq. 1) write(6,901) 
      ninter = nmsh - 1 
      rnold = flmax 
      isize=nmsh*ncomp 
                                                                        
!  Evaluate the right-hand side rhstri at the initial solution u by     
!  adding the new deferred corrections to the already-calculated        
!  rhs vector.                                                          
                                                                        
      call dcopy(nlbc, rhs, 1, rhstri, 1) 
      ind = nlbc 
      do 10 im = 1, ninter 
      do 10 ic = 1, ncomp 
         ind = ind + 1 
         rhstri(ind) = rhs(ind) + defnew(ic, im) 
   10 continue 
      ind = ninter*nmsh + nlbc + 1 
      call dcopy(ncomp-nlbc, rhs, 1, rhstri, 1) 
                                                                        
      call dssq  ( nmsh*ncomp, rhstri, 1, scale, sumsq ) 
      rnsq = (scale**2)*sumsq 
                                                                        
      iter = 0 
                                                                        
!  If the initial right-hand side is too large, do not even attempt to  
!  solve the nonlinear equations.                                       
                                                                        
      if (rnsq.gt.huge .or.                                             &
     &      (iorder.eq. 8 .and. rnsq.gt.xlarge)) then                   
         if (iprint .eq. 1) write (6,902) rnsq 
         iflag = -2 
         return 
      end if 
      call dcopy(ncomp*nmsh, rhstri, 1, rhs, 1) 
                                                                        
!  Statement 100 is the top of the iteration loop.                      
                                                                        
  100 continue 
                                                                        
!  If rnsq is sufficiently small, terminate immediately.                
                                                                        
      if (iprint .eq. 1) write(6,903) iter, rnsq 
      if (rnsq .le. epsmch) then 
         iflag = 0 
         return 
      endif 
                                                                        
      iter = iter + 1 
                                                                        
!  Solve for the step delu by solving a system involving the fixed      
!  Jacobian (whose LU factors are saved).  Copy the rhs array into      
!  tmprhs, which is then overwritten by blkslv.                         
                                                                        
      call dcopy(ncomp*nmsh, rhs, 1, tmprhs, 1) 
      call dcopy(ncomp*nmsh,tmprhs,1,delu,1) 
      job = 0 
      call crslve(Topblk, Nlbc, Ncomp, Ajac, Ncomp, 2*Ncomp,            &
     & Ninter, Botblk, Ncomp-Nlbc, Ipivot, Delu, job)                   
                                                                        
!  Compute the trial point utrial by adding delu to u.                  
                                                                        
      call matcop( nudim, ncomp, ncomp, nmsh, u, utrial ) 
      call maxpy( ncomp, nmsh, one, delu, ncomp, utrial ) 
                                                                        
!  compute the right-hand side vector rhstri and its squared            
!  two-norm at the trial point.                                         
                                                                        
      rnold = rnsq 
      call fneval(ncomp, nmsh, xx, ncomp, utrial, fval, fsub,rpar,ipar) 
      call rhscal (ncomp, nmsh, nlbc, xx, ncomp, utrial, defcor,        &
     &   fsub, gsub, rhstri, rnsq, fval, ftmp, uint,rpar,ipar)          
                                                                        
!  If rnsq strictly decreased, update the solution vector u             
!  and the right-hand side rhs.                                         
                                                                        
      better = .false. 
      if (rnsq .lt. rnold) then 
         better = .true. 
         call matcop( ncomp, nudim, ncomp, nmsh, utrial, u ) 
         call dcopy( ncomp*nmsh, rhstri, 1, rhs, 1 ) 
      endif 
                                                                        
!  Stop the fixed Jacobian iterations if there have been too            
!  many iterations, or if rnsq has not decreased by a factor            
!  of at least rngrow.                                                  
                                                                        
      if (iter .ge. lmtfrz .or. rnsq .gt. (rnold/rngrow)) then 
         if (better) then 
                                                                        
!  Setting iflag to -3 signals that, although the fixed Jacobian        
!  iterations did not succeed, the current point was an improvement     
!  on the previous one.  Hence, if we switch to a Newton procedure,     
!  the right-hand side does not need to be recalculated.                
                                                                        
            iflag = -3 
         else 
            iflag = -2 
         endif 
         if (iprint .eq. 1) write(6,904) iflag 
         return 
      endif 
                                                                        
!  Test for convergence using the ratio abs((change in u)/max(u,1)).    
                                                                        
      do 150 im = 1, nmsh 
      do 150 it = 1, ntol 
         itol = ltol(it) 
         er = abs(delu(itol,im))/max(abs(u(itol,im)), one) 
         if (er .gt. tolfct*tol(it)) go to 100 
  150 continue 
                                                                        
!  To exit from the loop here, the convergence tests have               
!  been passed.                                                         
                                                                        
      if (iprint .ge. 0) write(6,905) iter, rnsq 
                                                                        
      iflag = 0 
      return 
  901 format(1h ,'fixed Jacobian iterations') 
  902 format(1h ,'Large residual, rnsq =',1pe12.4) 
  903 format(1h ,'iter, rnsq',i5,1pe11.3) 
  904 format(1h ,'failure of fixed Jacobian, iflag =',i5) 
  905 format(1h ,'fixed Jacobian convergence',i5,1pe11.3) 
      end subroutine                                        
                                                                        
!***********************************************************************
                                                                        
      subroutine lineq( ncomp, nmsh, nlbc,                              &
     &    ludone, xx, nudim, u, defcor,                                 &
     &    delu, rhs, fval,                                              &
     &    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,                     &
     &    ajac, topblk, botblk, bhold, chold, ipivot,                   &
     &    fsub, dfsub, gsub, dgsub, iflag,rpar,ipar)                    
!********************************************************************** 
!     call by: bvpsol                                                   
!     calls to: lnrhs, jaccal, dcopy, colrow, dload, dload, clrslve     
!           maxpy                                                       
!*********************************************************************  
        !module
      use algprs
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension rpar(*), ipar(*) 
      dimension  xx(nmsh), u(nudim, nmsh), defcor(ncomp, nmsh-1) 
      dimension  delu(ncomp, nmsh), rhs(ncomp*nmsh) 
      dimension  fval(ncomp,nmsh), uint(ncomp), ftmp(ncomp) 
      dimension  topblk(nlbc,*), botblk(ncomp-nlbc,*) 
      dimension  dftmp1(ncomp, ncomp), dftmp2(ncomp, ncomp) 
      dimension  dgtm(ncomp), tmprhs(ncomp*nmsh) 
      dimension  ajac(ncomp, 2*ncomp, nmsh-1),                          &
     &               bhold(ncomp, ncomp, nmsh-1),                       &
     &               chold(ncomp, ncomp, nmsh-1)                        
      dimension  ipivot(ncomp*nmsh) 
      integer    job 
      logical    ludone 
!      external   fsub, dfsub, gsub, dgsub 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
!  blas: dcopy, dload                                                   
                                                                        
      parameter  ( one = 1.0_dp, zero = 0.0_dp ) 
                                                                        
!interface
      include 'sepbvp.h'

!  The routine lineq calculates the Newton step for a linear            
!  problem.  The Newton step is exact unless the Jacobian               
!  matrix is singular.                                                  
                                                                        
      isize=nmsh*ncomp 
      ninter = nmsh - 1 
      iflag = 0 
      if (.not. ludone) then 
                                                                        
!  Compute the right-hand side vector rhs.                              
                                                                        
         call lnrhs (ncomp, nmsh, nlbc, xx, nudim, u,                   &
     &          fsub, gsub, rhs, rnsq, fval, ftmp, uint,rpar,ipar)      
                                                                        
                                                                        
!  If the Jacobian for this mesh has not previously been                
!  calulated and factorized successfully, call jaccal.                  
!  The block-structured Jacobian matrix is stored in three              
!  matrices (topblk, ajac, and botblk).                                 
!  The matrices bhold and chold are also calculated in jaccal,          
!  and are saved for later use in outer routines.                       
                                                                        
         call jaccal (ncomp, nmsh, nlbc,                                &
     &      xx, nudim, u, fval, dgtm, dftmp1, dftmp2, uint,             &
     &      ajac, topblk, botblk, bhold, chold,                         &
     &      dfsub, dgsub,rpar,ipar)                                     
                                                                        
!  Call blkdcm to calculate the LU factors of the Jacobian.             
!  The factors are overwritten on the matrices topblk, ajac and botblk. 
!  Interchanges are represented in the integer array ipivot.            
                                                                        
      call dcopy(ncomp*nmsh,rhs,1,tmprhs,1) 
      call dcopy(ncomp*nmsh,tmprhs,1,delu,1) 
                                                                        
      job = 0 
      call colrow(isize, topblk,nlbc, ncomp, ajac, ncomp, 2*ncomp,      &
     &   ninter, botblk, ncomp-nlbc, ipivot, delu, iflag, job)          
                                                                        
                                                                        
                                                                        
         ludone = .true. 
                                                                        
!  Copy the rhs into the temporary vector tmprhs, which will be         
!  overwritten by blkslv.                                               
                                                                        
                                                                        
      else 
                                                                        
!  The right-hand side is the deferred correction array,                
!  padded with zeros at the boundary conditions.                        
                                                                        
         call dload(nlbc, zero, tmprhs(1), 1) 
         do 100 im = 1, ninter 
            loc = (im-1)*ncomp + nlbc + 1 
            call dcopy(ncomp, defcor(1,im), 1, tmprhs(loc), 1) 
  100    continue 
         nrhs = ninter*ncomp + nlbc + 1 
         call dload(ncomp-nlbc, zero, tmprhs(nrhs), 1) 
         call dcopy(ncomp*nmsh,tmprhs,1,delu,1) 
         job = 0 
                                                                        
      call crslve(topblk,nlbc,ncomp,ajac,ncomp,2*ncomp,ninter,          &
     &   botblk,ncomp-nlbc,ipivot,delu,job)                             
                                                                        
      endif 
                                                                        
!  Since the problem is linear, the Newton step  is exact.  The         
!  new u array is obtained by adding delu to u.                         
                                                                        
      call maxpy ( ncomp, nmsh, one, delu, nudim, u ) 
                                                                        
!      iflag = 0                                                        
      return 
                                                                        
  901 format(1h ,'Singular matrix') 
      end subroutine                                        
                                                                        
!*****************************************************************      
      subroutine newteq(ncomp, nmsh, nlbc,                              &
     &    rhsgiv, ntol, ltol, tol,                                      &
     &    xx, nudim, u, defcor,                                         &
     &    delu, rhs, fval,                                              &
     &    utrial, rhstri,                                               &
     &    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,             &
     &    ajac, topblk, botblk, bhold, chold, ipivot,                   &
     &    fsub, dfsub, gsub, dgsub, iter, iflag,rpar,ipar,frscal)       
!*****************************************************************      
!     call by: bvpsol                                                   
!     calls to:                                                         
!*****************************************************************      
 !module
      use algprs
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension rpar(*), ipar(*) 
      dimension  ltol(*), tol(*), xx(*) 
      dimension  fval(ncomp,*) 
      dimension  u(nudim, *), delu(ncomp, *), utrial(ncomp,nmsh) 
      dimension  rhs(ncomp*nmsh),  defcor(ncomp,*) 
      dimension  ajac(ncomp, 2*ncomp, *) 
      dimension  topblk(nlbc,*), botblk(ncomp-nlbc,*) 
      dimension  ftmp(*), uint(*), dgtm(ncomp) 
      dimension  dftmp1(ncomp, ncomp), dftmp2(ncomp, ncomp) 
      dimension  bhold(ncomp, ncomp, *), chold(ncomp, ncomp, *) 
      dimension  ipivot(*) 
      dimension  rhstri(ncomp*nmsh) 
      dimension  tmprhs(ncomp*nmsh), xmerit(ncomp, nmsh) 
      logical rhsgiv 
                                                                        
                                                                        
!      external   fsub, dfsub, gsub, dgsub 
                                                                        
      parameter  ( zero   = 0.0_dp, one    = 1.0_dp ) 
      parameter ( two = 2.0_dp, half = 0.5_dp, fourth = 0.25_dp ) 
      parameter ( tenth = 0.1_dp, ten = 10.0_dp, hund = 100.0_dp ) 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
      common/mchprs/flmin, flmax, epsmch 
                                                                        
      logical gtpdeb, imprvd, braktd, crampd, extrap, vset, wset 
      save  gtpdeb, mfsrch, epsaf, epsag, eta, rmu, tolabs, alfmax 
      save  tolrel, toltny 
                                                                        
      intrinsic  abs, max 

!interface
      include 'sepbvp.h'

                                                                       
!  blas: dcopy                                                          
                                                                        
      logical frscal 
!      save frscal                                                      
      parameter (cnvfct = 0.1_dp ) 
!      data  frscal/.true./                                             
! frscal has been updated in bvpsol                                     
!                                                                       
      data  alfsml/1.0d-4/,  alfmax/1.1_dp/ 
      data  imerit/1/, lmtnwt/70/ 
      data  shrfct/100.0_dp/, stpfct/2.0_dp/ 
      data  gtpdeb/.false./, mfsrch/5/ 
      data  eta/.999999_dp/, rmu/1.0d-6/ 

                                                                        
!  The routine newteq performs Newton iterations with a line            
!  search, to solve the nonlinear equations.                            
                                                                        
                                                                        
!  Set up constants if this is the first call to newteq.                
                                                                        
      if (frscal) then 
         frscal = .false. 
         epsaf = epsmch 
         epsag = epsmch 
         tolabs = epsmch 
         tolrel = epsmch 
         toltny = epsmch 
      endif 
      ninter = nmsh - 1 
                                                                        
      if (iprint .eq. 1) write(6,901) 
                                                                        
!  A Newton method with line search and watchdog safeguarding           
!  is used to try to solve the nonlinear equations.                     
                                                                        
!  Initialize iter (the counter of Newton iterations) and alfold        
!  (the step taken at the previous iteration).                          
                                                                        
      iter = -1 
      alfold = one 
      alfa = zero 
                                                                        
                                                                        
                                                                        
      if (.not. rhsgiv) then 
!  If necessary, evaluate the right-hand side at the initial u.         
                                                                        
         call rhscal (ncomp, nmsh, nlbc, xx, nudim, u, defcor,          &
     &      fsub, gsub, rhs, rnsq, fval, ftmp, uint,rpar,ipar)          
      endif 
                                                                        
                                                                        
!  At any given Newton iteration, rnprev is the value of rnsq at        
!  the immediately preceding Newton iteration.                          
                                                                        
      rnprev = flmax 
      rnbest = flmax 
      if (.not. pdebug .and. iprint .ge. 0) write (6,902) 
                                                                        
!  Initialize counter of watchdog iterations.                           
                                                                        
      itwtch = 0 
                                                                        
!  Statement 100 is the top of the Newton iteration loop.               
                                                                        
  100 continue 
                                                                        
      iter = iter + 1 
                                                                        
      if (iprint .eq. 1) write(6,910) iter 
                                                                        
!  If there have been too many Newton iterations, terminate.            
                                                                        
      if (iter .ge. lmtnwt) then 
         if (iprint .ge. 0) write(6,903) 
         iflag = -2 
         return 
      endif 
                                                                        
!  The vector rhs is the right-hand side at the current iterate,        
!  and rnsq is its squared two-norm.                                    
!  Perform watchdog tests, using the unscaled merit function (rnsq)     
!  as the watchdog function.  The routine wtchdg updates rnbest         
!  and itwtch.  If iflwat is not zero, this sequence of Newton          
!  iterations is terminated.                                            
                                                                        
      iflwat = 0 
                                                                        
      call wtchdg ( iter, rnsq, rnbest, rnprev, itwtch,                 &
     &                alfold, iflwat )                                  
                                                                        
      if (iflwat .ne. 0) then 
         if (iprint .ge. 0) write(6,904) iter 
         iflag = -3 
         return 
      endif 
                                                                        
!  Watchdog tests are passed.  Proceed with the Newton iteration.       
!  Call jaccal to evaluate the block-structured Jacobian matrix,        
!  which is stored in three matrices (topblk, ajac, and botblk).        
!  The matrices bhold and chold are saved for use in later              
!  calculations in the outer routine.                                   
                                                                        
                                                                        
!  If rnsq is sufficiently small, terminate immediately.                
!  Note that the stored Jacobian does not correspond exactly            
!  to the final point.                                                  
                                                                        
      if (rnsq .le. epsmch) then 
         if (iprint .ge. 0)  write(6,906) iter, rnsq 
         iflag = 0 
         return 
      endif 
                                                                        
      call jaccal (ncomp, nmsh, nlbc,                                   &
     &    xx, nudim, u, fval, dgtm, dftmp1, dftmp2, uint,               &
     &    ajac, topblk, botblk, bhold, chold,                           &
     &    dfsub, dgsub,rpar,ipar)                                       
                                                                        
!  blkdcm is called to calculate the LU factors of the Jacobian,        
!  which are overwritten on topblk, ajac and botblk.                    
!  Interchanges are represented in the integer array ipivot.            
                                                                        
!   Solve for the Newton step delu.  Copy the rhs array into tmprhs,    
!   which is then overwritten by blkslv.                                
                                                                        
      call dcopy(ncomp*nmsh, rhs, 1, tmprhs, 1) 
      call dcopy(ncomp*nmsh,tmprhs,1,delu,1) 
      job = 0 
      call colrow(nmsh*ncomp,topblk,nlbc,ncomp,ajac,ncomp,2*ncomp,      &
     &   ninter,botblk,ncomp-nlbc,ipivot,delu,iflag,job)                
                                                                        
      if (iprint .ge. 0 .and. iflag.ne.0) write(6,905) iter 
      if (iflag .ne. 0) return 
!        the jacobian is singular                                       
                                                                        
                                                                        
!  If imerit = 1, the line search is based on the scaled merit function,
!  the squared two-norm of the solution xmerit of the linear system     
!      (Jacobian)*xmerit = rhs,                                         
!  where (Jacobian) is the Jacobian at the current Newton iterate.      
!  Thus the initial value of the scaled merit function is simply        
!  the squared two-norm of the Newton step delu itself.                 
                                                                        
      if (imerit.eq.1) then 
         call mssq( ncomp, nmsh, delu, xmscal, xmsq ) 
         fmtry = (xmscal**2)*xmsq 
      else 
!  The unscaled merit function is simply the squared two-norm of rhs.   
         fmtry = rnsq 
      end if 
                                                                        
!  fa and oldg represent the merit function and its gradient            
!  at the initial point of the line search.                             
                                                                        
      fa = fmtry 
      oldg = -two*fa 
      alfa = zero 
      if (iprint .eq. 1) write (6,908) alfa, fmtry, rnsq 
                                                                        
!  On the first Newton iteration, the initial trial step is unity.      
!  On subsequent iterations, the initial step is not allowed to         
!  be more than the factor stpfct larger than the final step at         
!  the immediately preceding iteration.                                 
                                                                        
                                                                        
      alfa = one 
      if(stpfct*alfold .lt. one) alfa = stpfct*alfold 
                                                                        
      if (alfa .lt. alfsml) alfa = alfsml 
                                                                        
      fmold = fa 
      inform = -1 
                                                                        
!  Statement 150 is the top of the inner line search iteration.         
!  The line search routine getptq has been altered so that it           
!  terminates with an indication of success as soon as a                
!  strictly lower value of the merit function is found.  Note that      
!  this is a much less strict requirement than the usual sufficient     
!  decrease conditions.                                                 
                                                                        
                                                                        
  150 continue 
      iwr = 6 
      call getptq (gtpdeb, mfsrch, iwr, alfmax, alfsml, alfuzz,         &
     &      epsaf, epsag,                                               &
     &      eta, fmtry, fmold, oldg, rmu, tolabs, tolrel, toltny,       &
     &      imprvd, inform, nfsrch, alfa, alfbst, fbest,                &
     &      braktd, crampd, extrap, vset, wset, nsamea, nsameb,         &
     &      alin, blin, fa, factor, fv, fw, xtry, xv, xw)               
                                                                        
                                                                        
                                                                        
                                                                        
!  inform = 1, 2 or 3 indicates success in finding an acceptable point. 
!  inform = 4 means alfmax is too small (this should never happen here, 
!  since alfmax is set always to 1.1).                                  
!  inform = 5 means that a decrease was not achieved for any step       
!  greater than alfsml.                                                 
!  inform = 6 means a better point could not be found (the minimum      
!  probably lies too close to alfa=0).                                  
!  inform = 7 means that the gradient at alfa=0 (oldg) is positive      
!  (this cannot happen here, since oldg=-two*fa, and fa is a non-negativ
!  number)                                                              
                                                                        
      if (pdebug) write(6,907) inform, alfa 
                                                                        
      if (inform .eq. 5) then 
          iflag = -5 
          return 
      elseif (inform .eq. 4 .or. inform .eq. 7) then 
         iflag = -4 
         return 
      elseif (inform .eq. 0) then 
                                                                        
!  inform = 0 means that a new function value should be obtained        
!  with the step alfa.                                                  
!  We may override alfa from getptq by requiring that the step is not   
!  allowed to decrease by more than a factor of shrfct during           
!  a line search iteration.                                             
!                                                                       
         if (alfa .lt. alfold/shrfct) alfa = alfold/shrfct 
         alfold = alfa 
                                                                        
!  Define the next iterate utrial = u + alfa*delu.                      
!  Call fneval and rhscal to evaluate the right-hand side               
!  rhstri at utrial.                                                    
!  The vector rhstri is stored separately, and rhs is overwritten       
!  only when an improved point is found.                                
                                                                        
         call matcop ( nudim, ncomp, ncomp, nmsh, u, utrial) 
         call maxpy ( ncomp, nmsh, alfa, delu, ncomp, utrial ) 
         call fneval(ncomp, nmsh, xx, ncomp, utrial, fval, fsub,        &
     &    rpar,ipar)                                                    
         call rhscal (ncomp, nmsh, nlbc, xx, ncomp, utrial, defcor,     &
     &      fsub, gsub, rhstri, rnsqtr, fval, ftmp, uint,rpar,ipar)     
                                                                        
         fmold = fmtry 
         if (imerit .eq. 1) then 
                                                                        
!  Solve a linear system to obtain the 2-d array xmerit whose squared   
!  norm is the scaled merit function.   The LU factors of the Jacobian  
!  have already been calculated by blkdcm.                              
!  Copy rhstri into tmprhs, which is overwritten by blkslv.             
                                                                        
      call dcopy(ncomp*nmsh, rhstri, 1, tmprhs, 1) 
      call dcopy(ncomp*nmsh,tmprhs,1,xmerit,1) 
      job = 0 
      call crslve(topblk, nlbc, ncomp, ajac, ncomp, 2*ncomp, ninter,    &
     &   botblk, ncomp-nlbc, ipivot, xmerit,job)                        
      call mssq( ncomp, nmsh, xmerit, xscale, xsolsq ) 
      fmtry = (xscale**2)*xsolsq 
         else 
                                                                        
!  The unscaled merit function is the squared two-norm of the right-hand
!  side.                                                                
            fmtry = rnsqtr 
         end if 
         if (iprint .eq. 1) write (6,908) alfa, fmtry, rnsqtr 
         go to 150 
      endif 
                                                                        
!  To reach here, inform must be 1, 2, 3, or 6, and the line search     
!  has found a strictly lower value of the merit function.              
!  Store the new Newton iterate in u, and the corresponding rhs         
!  vector in rhs.                                                       
                                                                        
      rnprev = rnsq 
      rnsq = rnsqtr 
      call matcop (ncomp, nudim, ncomp, nmsh, utrial, u) 
      call dcopy(ncomp*nmsh, rhstri, 1, rhs, 1) 
      if (iprint .ge. 0) write(6,909) iter, alfa, fmtry, rnsq 
                                                                        
!  Now test for convergence using the ratio of the Newton step          
!  for each component with max(1, abs(current solution estimate)).      
!  If the test fails for any element of u, branch back to the           
!  top of the Newton iteration.                                         
                                                                        
      do 160 im = 1, nmsh 
      do 160 it = 1, ntol 
         icmp = ltol(it) 
         er = abs(delu(icmp,im))/max(abs(u(icmp,im)), one) 
         if (er .gt. cnvfct*tol(it)) go to 100 
  160 continue 
                                                                        
                                                                        
      if (iprint .ge. 0) write(6, 906) iter+1, rnsq 
      iflag = 0 
                                                                        
!  To fall through the above loop, the termination test for a           
!  sufficiently small delu is satisfied.                                
!  Note that the stored Jacobian and its factorization do not           
!  correspond to the final solution.                                    
                                                                        
                                                                        
      return 
                                                                        
  901 format(1h ,'start Newton iterations') 
  902 format(1h ,' iter',                                               &
     &              7x,'alfa',6x,'merit',7x,'rnsq')                     
  903 format(1h ,'Too many Newton iterations') 
  904 format(1h ,'Watchdog tests fail, iter =', i5) 
  905 format(1h ,'Singular Jacobian, iter=',i5) 
  906 format(1h ,'Convergence, iter =',i5,4x,'rnsq =',1pe12.3) 
  907 format(1h ,'inform, alfa after getptq',i5,3x, 1pe11.3) 
  908 format(1h ,'alfa, merit, rnsq',3(1pe11.3)) 
  909 format(1h ,i5,3(1pe11.3)) 
  910 format(1h ,'Newton iteration',i5) 
      end subroutine                                        
                                                                        
                                                                        
      subroutine wtchdg ( iter, wmerit, wmbest, wmprev,                 &
     &      itwtch, alfold, iflag )                                     
                                                                        
!  Logic for watchdog tests.                                            
                                                                        
        implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      integer, parameter :: itonew = 5, itwtmx = 8
      parameter ( grfct = 80.0_dp ) 
      parameter ( half = 0.5_dp ) 
                                                                        
!  Perform watchdog tests in two forms:                                 
!  (1) to determine whether a sufficient decrease in the                
!  watchdog merit function has occurred within the most recent          
!  sequence of itwtmx iterations;                                       
!  (2) to determine whether the watchdog merit function has increased   
!  too much in a single iteration after itonew Newton iterations        
!  have been performed.  This allows the merit function to increase     
!  wildly only during the first itonew iterations.                      
                                                                        
!  wmbest is the smallest watchdog merit function achieved in this      
!  sequence of Newton iterations.                                       
!  wmprev is the watchdog merit function from the immediately           
!  preceding Newton iteration.                                          
                                                                        
!  itwtch counts the number of iterations without an improvement        
!  in the unscaled merit function.                                      
                                                                        
!      write(6,99) iter, wmerit, wmbest, wmprev                         
!      write(6,98) itwtch, alfold                                       
!   99 format(1h ,'iter,wmer,wbest,wprev',i5,3(1pe15.5))                
!   98 format(1h ,'itwtch,alfold',i5,1pe15.5)                           
      iflag = 0 
      if (wmerit .le. wmbest) then 
                                                                        
!  The current watchdog merit function is the best.                     
                                                                        
         wmbest = wmerit 
         itwtch = 0 
         return 
      endif 
                                                                        
!  The current merit function is not the best.                          
                                                                        
      itwtch = itwtch + 1 
                                                                        
!  Do not apply watchdog tests if (1) the previous step alfold          
!  exceeds 1/2, or (2) the watchdog merit function decreased in         
!  the immediately preceding iteration and itwtch does not              
!  exceed twice its maximum.                                            
                                                                        
      if (alfold .ge. half) return 
      if (wmerit .le. wmprev .and. itwtch .le. 2*itwtmx) return 
                                                                        
                                                                        
!  If more than itwtmx iterations have occurred without                 
!  an overall improvement in the watchdog merit function,               
!  signal for termination.                                              
                                                                        
      if (itwtch .ge. itwtmx) then 
         iflag = -1 
                                                                        
!  If a too-large increase in the watchdog merit function               
!  compared to the best value occurred, and iter .ge. itonew,           
!  signal for termination.                                              
                                                                        
      elseif (iter .ge. itonew .and.                                    &
     &          wmerit .gt. grfct*wmbest) then                          
          iflag = -1 
      endif 
      return 
      end subroutine                                        
                                                                        
      subroutine fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,rpar,ipar) 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension rpar(*),ipar(*) 
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp,nmsh) 
!      external fsub 

!interface
      include 'sepbvp.h'
                                                                        
!  fneval evaluates the function values (from fsub) for                 
!  a given mesh xx and array u, and stores the values                   
!  in the array fval.                                                   
                                                                        
      call fsub (ncomp, xx(1), u(1,1), fval(1,1),rpar,ipar) 
                                                                        
      do 50 im = 1, nmsh-1 
         hmsh = xx(im+1) - xx(im) 
         call fsub (ncomp, xx(im+1), u(1,im+1), fval(1,im+1),rpar,ipar) 
   50 continue 
      return 
      end subroutine                                        
                                                                        
                                                                        
      subroutine jaccal (ncomp, nmsh, nlbc,                             &
     &   xx, nudim, u, fval,                                            &
     &   dgtm, dftm1, dftm2, uint,                                      &
     &   ajac, topblk, botblk, bhold, chold,                            &
     &   dfsub, dgsub,rpar,ipar)                                        
 
!module
      use algprs                                                   
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension rpar(*), ipar(*) 
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp,nmsh) 
      dimension dgtm(ncomp) 
      dimension dftm1(ncomp, ncomp), dftm2(ncomp, ncomp),               &
     &             uint(ncomp)                                          
      dimension ajac(ncomp, 2*ncomp, nmsh-1) 
      dimension topblk(nlbc, ncomp), botblk(ncomp-nlbc,ncomp) 
      dimension bhold(ncomp, ncomp, nmsh-1),                            &
     &             chold(ncomp, ncomp, nmsh-1)                          
                                                                        
!      external  dfsub, dgsub 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
!  blas: dcopy, ddot                                                    
                                                                        
      parameter ( zero = 0.0_dp, half = 0.5_dp, eighth = 0.125_dp ) 
      parameter ( four = 4.0_dp, six = 6.0_dp ) 
      parameter ( one = 1.0_dp, three = 3.0_dp, twelve = 12.0_dp ) 
!interface
      include 'sepbvp.h'
                                                                        
!      if (pdebug) write(6,901)                                         
                                                                        
      ninter = nmsh - 1 
                                                                        
!      if (pdebug) write(6,902)                                         
      do 110 i = 1, nlbc 
         call dgsub (i, ncomp, u(1,1), dgtm,rpar,ipar) 
         call dcopy(ncomp, dgtm(1), 1, topblk(i,1), nlbc) 
!         if (pdebug) write(6,903) i, (topblk(i,j),j=1,ncomp)           
  110 continue 
                                                                        
      call dfsub (ncomp, xx(1), u(1,1), dftm1(1,1),rpar,ipar) 
                                                                        
!  on entry to jaccal, the array fval contains the function values      
!  at (xx(im), u(ic,im)), ic=1,...,ncomp and im = 1,...,nmsh,           
!  calculated by a preceding call of rhscal with the same xx and u      
!  arrays.                                                              
                                                                        
!      if (pdebug) write(6,904)                                         
                                                                        
      do 200 im = 1, ninter 
                                                                        
         hmsh = xx(im+1) - xx(im) 
                                                                        
         do 120 ic = 1, ncomp 
            uint(ic) = half*(u(ic,im) + u(ic,im+1))                     &
     &         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))              
  120    continue 
         xhalf = half*(xx(im+1) + xx(im)) 
         call dfsub (ncomp, xhalf, uint, dftm2(1,1),rpar,ipar) 
         do 140 ic = 1, ncomp 
            do 130 jc = 1, ncomp 
               dsq = ddot(ncomp, dftm2(ic,1), ncomp,                    &
     &                       dftm1(1,jc), 1)                            
               ajac(ic,jc,im) = -hmsh*(dftm1(ic,jc)/six                 &
     &             + dftm2(ic,jc)/three + hmsh*dsq/twelve)              
  130       continue 
            ajac(ic,ic,im) = ajac(ic,ic,im) - one 
!            if (pdebug) write(6,905) im, ic,                           
!     *            (ajac(ic,jc,im), jc=1,ncomp)                         
  140    continue 
                                                                        
         call dfsub (ncomp, xx(im+1), u(1,im+1), dftm1(1,1),rpar,ipar) 
         do 170 ic = 1, ncomp 
            do 160 jc = 1, ncomp 
               dsq = ddot(ncomp, dftm2(ic,1), ncomp,                    &
     &                         dftm1(1,jc), 1)                          
               ajac(ic,jc+ncomp,im) = -hmsh*(dftm1(ic,jc)/six           &
     &               + dftm2(ic,jc)/three - hmsh*dsq/twelve)            
  160       continue 
            call dcopy(ncomp, ajac(ic,ncomp+1,im), ncomp,               &
     &                   chold(ic,1,im), ncomp)                         
            call dcopy(ncomp, dftm1(ic,1), ncomp,                       &
     &                   bhold(ic,1,im), ncomp)                         
            ajac(ic,ic+ncomp,im) = ajac(ic,ic+ncomp,im) + one 
            chold(ic,ic,im) = ajac(ic,ic+ncomp,im) 
!            if (pdebug) write(6,905) im, ic,                           
!     *                    (ajac(ic,jc+ncomp,im),jc=1,ncomp)            
  170    continue 
                                                                        
                                                                        
  200 continue 
!      if (pdebug) write(6,906)                                         
      do 220 i = nlbc+1, ncomp 
         call dgsub (i, ncomp, u(1, nmsh), dgtm,rpar,ipar) 
         call dcopy(ncomp, dgtm(1), 1, botblk(i-nlbc,1), ncomp-nlbc) 
!         if (pdebug) write(6,903) i,(botblk(i-nlbc,j), j=1,ncomp)      
  220 continue 
!      write(6,991)                                                     
!991   format(1x,'topblk')                                              
!      write(6,992) topblk(1,1),topblk(1,2)                             
!992   format(1x,2g22.10)                                               
!      write(6,993)                                                     
!993   format(1x,'main jacobian')                                       
!      do 994 iu=1,ninter                                               
!      do 994 iv=1,ncomp                                                
!      write(6,996) ajac(iv,1,iu),ajac(iv,2,iu),ajac(iv,3,iu),          
!     +ajac(iv,4,iu)                                                    
!996   format(1x,4f12.7)                                                
!994   continue                                                         
!      write(6,995)                                                     
!995   format(1x,'botblk')                                              
!      write(6,992) botblk(1,1),botblk(1,2)                             
                                                                        
                                                                        
      return 
                                                                        
  901 format(1h ,'jaccal') 
  902 format(1h ,'topblk') 
  903 format(1h ,i5,6(1pe11.3)) 
  904 format(1h ,'ajac') 
  905 format(1h ,2i5,5(1pe11.3)) 
  906 format(1h ,'botblk') 
      end subroutine                                        
                                                                        
                                                                        
                                                                        
      subroutine lnrhs (ncomp, nmsh, nlbc,                              &
     &   xx, nudim, u, fsub, gsub,                                      &
     &   rhs, rnsq, fval, ftmp, uint,rpar,ipar)                         
                                                                        
        implicit integer (i-n)
       implicit real(dp)(a-h,o-z) 
                                                                        
!  This subroutine is designed to calculate the right-hand              
!  side for linear problems.                                            
      dimension rpar(*),ipar(*) 
      dimension xx(*), u(nudim,*) 
      dimension rhs(*), fval(ncomp,*), ftmp(*), uint(*) 
!      external fsub, gsub 
                                                                        
      parameter ( zero = 0.0_dp, half = 0.5_dp, eighth = 0.125_dp ) 
      parameter ( one = 1.0_dp, four = 4.0_dp, six = 6.0_dp ) 
                                                                        
      common/mchprs/flmin, flmax, epsmch 
      intrinsic abs 
 
!interface
      include 'sepbvp.h'
                                                                       
!  blas: dssq                                                           
                                                                        
      ninter = nmsh - 1 
      rnsq = zero 
                                                                        
!  first, process the left-hand boundary conditions.                    
                                                                        
      do 20 i = 1, nlbc 
         call gsub (i, ncomp, u(1,1), wg,rpar,ipar) 
         rhs(i) = -wg 
   20 continue 
                                                                        
!  Next, process the interior mesh points.  The fval array              
!  contains the function values from fsub at xx and u.                  
                                                                        
      do 50 im = 1, ninter 
         hmsh = xx(im+1) - xx(im) 
         do 30 ic = 1, ncomp 
            uint(ic) = half*(u(ic,im) + u(ic,im+1))                     &
     &         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))              
   30    continue 
         xhalf = half*(xx(im) + xx(im+1)) 
         call fsub (ncomp, xhalf, uint, ftmp,rpar,ipar) 
         loc = (im-1)*ncomp + nlbc 
         do 40 ic = 1, ncomp 
            rhs(loc+ic) = -u(ic,im+1) + u(ic,im)                        &
     &        + hmsh*                                                   &
     &            (fval(ic,im) + fval(ic,im+1) + four*ftmp(ic))/six     
   40    continue 
   50 continue 
                                                                        
      nrhs = ninter*ncomp 
      do 60 ii = nlbc+1, ncomp 
         call gsub (ii, ncomp, u(1,nmsh), wg,rpar,ipar) 
         rhs(nrhs+ii) = -wg 
   60 continue 
                                                                        
      call dssq  ( nmsh*ncomp, rhs, 1, scale, sumsq ) 
      rnsq = (scale**2)*sumsq 
                                                                        
      return 
      end subroutine                                        
                                                                        
                                                                        
      subroutine rhscal (ncomp, nmsh, nlbc,                             &
     &   xx, nudim, u, defcor,                                          &
     &   fsub, gsub,                                                    &
     &   rhs, rnsq, fval, ftmp, uint,rpar,ipar)                         
 !module
      use algprs                                                          
      implicit integer (i-n)
      implicit real(dp)(a-h,o-z) 
                                                                        
!  This subroutine constructs the (ncomp*nmsh)-dimensional              
!  vector rhs, which is the right-hand side of the Newton equations.    
!  The ncomp by nmsh array fval is assumed to have been calculated      
!  elsewhere by routine fneval.                                         
              dimension rpar(*),ipar(*) 
      dimension  xx(nmsh), u(nudim,nmsh), defcor(ncomp,nmsh-1) 
      dimension  rhs(ncomp*nmsh), fval(ncomp,nmsh) 
      dimension  ftmp(ncomp), uint(ncomp) 
!      external   fsub, gsub 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
      common/mchprs/flmin, flmax, epsmch 
      intrinsic abs 
 
!interface
      include 'sepbvp.h'
                                                                       
!  blas: dssq                                                           
                                                                        
      parameter ( zero = 0.0_dp, half = 0.5_dp, eighth = 0.125_dp ) 
      parameter ( one = 1.0_dp, four = 4.0_dp, six = 6.0_dp ) 
                                                                        
!      if (pdebug) write(6,901)                                         
                                                                        
!  ninter is the number of intervals in the mesh (one less than the     
!  number of mesh points)                                               
                                                                        
      ninter = nmsh - 1 
      rnsq = zero 
                                                                        
!  First, process the left-hand boundary conditions.                    
                                                                        
      do 20 i = 1, nlbc 
         call gsub (i, ncomp, u(1,1), wg,rpar,ipar) 
         rhs(i) = -wg 
   20 continue 
                                                                        
!  Next, process the interior mesh points.  The fval array              
!  contains the function values from fsub at xx and u.                  
                                                                        
      do 50 im = 1, ninter 
         hmsh = xx(im+1) - xx(im) 
         do 30 ic = 1, ncomp 
            uint(ic) = half*(u(ic,im) + u(ic,im+1))                     &
     &         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))              
   30    continue 
         xhalf = half*(xx(im) + xx(im+1)) 
         call fsub (ncomp, xhalf, uint, ftmp,rpar,ipar) 
         loc = (im-1)*ncomp + nlbc 
         do 40 ic = 1, ncomp 
            rhs(loc+ic) = -u(ic,im+1) + u(ic,im) + defcor(ic,im)        &
     &        + hmsh*                                                   &
     &            (fval(ic,im) + fval(ic,im+1) + four*ftmp(ic))/six     
                                                                        
   40    continue 
   50 continue 
                                                                        
      nrhs = ninter*ncomp 
      do 60 ii = nlbc+1, ncomp 
         call gsub (ii, ncomp, u(1,nmsh), wg,rpar,ipar) 
         rhs(nrhs+ii) = -wg 
   60 continue 
                                                                        
                                                                        
      call dssq  ( nmsh*ncomp, rhs, 1, scale, sumsq ) 
      rnsq = (scale**2)*sumsq 
                                                                        
      if (pdebug) then 
         write (6,902) rnsq 
         write(6,903) 
         write(6,904) (rhs(i), i=1,ncomp*nmsh) 
      endif 
      return 
                                                                        
  901 format(1h ,'rhscal') 
  902 format(1h ,'rnsq',1pe11.3) 
  903 format(1h ,'rhs vector') 
  904 format(1h ,(7(1pe11.3))) 
      end subroutine                                        
                                                                        
      subroutine dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
!module
      use algprs 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension xx(*), xxold(*) 
      logical maxmsh 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
!  blas: dcopy                                                          
                                                                        
      parameter (half = 0.5_dp) 
                                                                        
                                                                        
!  This routine is used to double the mesh, i.e., produce a mesh        
!  with twice as many intervals in which each new interval is           
!  half the corresponding old interval.                                 
                                                                        
!  On entry to dblmsh, the integer nmsh and the array xx                
!  specify a set of mesh points xx(1),..., xx(nmsh) (assumed            
!  to be in ascending order).                                           
                                                                        
!  If the number of mesh points in the doubled mesh would               
!  exceed the maximum allowed number nmax, the flag maxmsh is           
!  set to true, and we exit without changing any other parameters.      
                                                                        
!  Otherwise, nmold is set to the old number of mesh points,            
!  xxold is set to the old mesh, nmsh is the new number of mesh         
!  points, and xx contains the new mesh points.                         
                                                                        
      nmold = nmsh 
      call dcopy(nmold, xx, 1, xxold, 1) 
                                                                        
      ninnew = 2*(nmsh-1) 
      nmnew = ninnew + 1 
      if(nmnew .ge. nmax) then 
         if (iprint .ge. 0)  write(6,901) nmnew 
         nmsh = nmold 
         call dcopy(nmold, xxold, 1, xx, 1) 
         maxmsh = .true. 
         return 
      endif 
      maxmsh = .false. 
                                                                        
!  Loop backwards through the old mesh points to create the new ones.   
                                                                        
      xx(nmnew) = xx(nmsh) 
      do 100 i = ninnew, 4, -2 
         id2 = i/2 
         xx(i) = half*(xx(i+1) + xx(id2)) 
         xx(i-1) = xx(id2) 
  100 continue 
                                                                        
!  Calculate the new xx(2). xx(1) remains unchanged.                    
      xx(2) = half*(xx(3) + xx(1)) 
      nmsh = nmnew 
      if(iprint .ge. 0)  write(6,902) nmsh 
      return 
  901 format (1h , ' dblmsh.  maximum mesh exceeded, nmnew =', i8) 
  902 format (1h , ' dblmsh.  the doubled mesh has ', i8,' points.') 
      end subroutine                                        
                                                                        
                                                                        
      subroutine selmsh(ncomp, nmsh, ntol, ltol, tol,                   &
     &     nfxpnt, fixpnt, ipow, nmax,                                  &
     &     xx, nudim, u, ermeas, irefin, ihcomp,                        &
     &     nmold, xxold, ermx, ddouble , maxmsh)                        
 !module
      use algprs                                                         
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
                                                                        
      dimension  ltol(ntol), tol(ntol), fixpnt(*) 
      dimension  xx(*), u(nudim, *), ermeas(ncomp,*) 
      dimension  irefin(nmsh-1), ihcomp(nmsh-1) 
      dimension  xxold(*), ermx(*) 
      logical    ddouble , maxmsh 
                                                                        
      intrinsic abs, max, int 
                                                                        
!  blas: dcopy                                                          
!  real(dp) log                                                
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      parameter  ( zero = 0.0_dp, one = 1.0_dp, onep1 = 1.1_dp ) 
      parameter  ( erdcid = 5.0_dp ) 
      parameter  ( phitst = 0.1_dp ) 
                                                                        
      logical first 
      save    first, rlndec 
      data    first / .true. / 
                                                                        
!  The routine selmsh performs selective mesh refinement, depending     
!  on the error measure ermeas.                                         
                                                                        
      if (first) then 
         first = .false. 
         rlndec = log(erdcid) 
      endif 
                                                                        
      maxmsh = .false. 
                                                                        
      if (pdebug) write(6,901) nmsh, ipow 
                                                                        
      frcpow = one/ipow 
      ddouble = .false. 
      nmold = nmsh 
      ninter = nmsh - 1 
                                                                        
!  Copy the current mesh into the xxold array.                          
                                                                        
                                                                        
      call dcopy(nmold, xx, 1, xxold, 1) 
      ithres = 0 
      thres = one 
                                                                        
!  On input, the array ermeas represents some error measure defined     
!  over the components and mesh intervals (not mesh points).            
!  It is normalized in the following loop with respect to the           
!  tolerance array and the current solution.                            
!  The value errmax gives the maximum normalized error.                 
                                                                        
      errmax = zero 
      do 120 im = 1, ninter 
         ermx(im) = zero 
         do 110 it = 1, ntol 
            jcomp = ltol(it) 
            denom = tol(it)*max(one, abs(u(jcomp,im))) 
            ems = ermeas(jcomp,im) 
            ermeas(jcomp,im) = abs(ems)/denom 
!            if (pdebug .and. ermeas(jcomp,im) .ge. thres)              
!     *             write(6,902) im,jcomp,ems,ermeas(jcomp,im)          
            err = ermeas(jcomp, im) 
            if (err .ge. ermx(im)) then 
                ermx(im) = err 
                ihcomp(im) = jcomp 
            endif 
  110    continue 
         errmax = max(ermx(im), errmax) 
  120 continue 
                                                                        
                                                                        
      if (pdebug) write(6,903) errmax 
                                                                        
      if (errmax .gt. zero .and. errmax .le. erdcid) then 
                                                                        
!  If errmax > 0 and .le. erdcid, find the smallest integer exponent ii 
!  such that (erdcid**ii)*errmax > erdcid.                              
                                                                        
         if(errmax .gt. one) then 
            ii = 1 
            decii = erdcid 
         else 
            ilg = -log(errmax)/rlndec 
            ii = 2 + ilg 
            decii = erdcid**ii 
         endif 
                                                                        
!  Multiply error measures by erdcid**ii.                               
                                                                        
         errmax = decii*errmax 
         do 140 im = 1, ninter 
            ermx(im) = decii*ermx(im) 
            do 140 it = 1, ntol 
               jcomp = ltol(it) 
               ermeas(jcomp,im) = decii*ermeas(jcomp,im) 
  140    continue 
      endif 
                                                                        
  200 continue 
                                                                        
!  For each interval im,  the integer irefin(im) is calculated          
!  based on an equidistrbution procedure involving the                  
!  threshold value thres.  If irefin(im) > 1, we add                    
!  points to interval im.  If irefin(im) = 1, we check whether          
!  point im can be removed from the mesh.                               
                                                                        
!  nmest is a lower bound on the number of points in the new mesh.      
!  We do not know in advance how many points will be removed,           
!  so nmest is computed by assuming all eligible points are removed.    
                                                                        
      nmest = nmsh 
      do 220 im = 1, ninter 
         if (ermx(im).ge.thres) then 
            irefin(im) = int(ermx(im)**frcpow) + 1 
            nmest = nmest + irefin(im) - 1 
         else 
            irefin(im) = 1 
            nmest = nmest - 1 
         endif 
  220 continue 
!      if (pdebug) write(6,904) nmest, (irefin(i), i=1,ninter)          
                                                                        
      if (nmest .gt. nmax) then 
                                                                        
         go to 360 
                                                                        
      elseif (nmest-1 .gt. 3*ninter) then 
                                                                        
         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
         ddouble = .true. 
         return 
      endif 
                                                                        
!  It appears that we can perform the desired selective mesh            
!  refinement.                                                          
                                                                        
!  Now begin running through the mesh, adding and possibly deleting     
!  points as indicated by the irefin array.                             
                                                                        
!  The integer new is a count of the number of intervals in             
!  the tentative mesh being generated by the refinement strategy.       
                                                                        
      new = 1 
                                                                        
!  The first interval is treated as a special case, since xx(1)         
!  always remains in the mesh, and cannot be a fixed point.             
                                                                        
      rlen = xxold(2) - xx(1) 
      slen = rlen 
      if (irefin(1).gt.1) then 
         dx = rlen/irefin(1) 
         do 230 j = 2, irefin(1) 
            new = new + 1 
            xx(new) = xx(1) + (j-1)*dx 
  230    continue 
      endif 
                                                                        
!  The fixed points specified by the fixpnt array cannot be             
!  removed from the mesh.  The value fxnext indicates the 'next'        
!  fixed point.  When no further fixed points remain to be processed    
!  (or if nfxpnt = 0), fxnext is set to a value strictly larger than    
!  the last mesh point, so that no mesh point can equal fxnext.         
!  This way we need to compare only the values of xxold(i)              
!  and fxnext.                                                          
                                                                        
      ifxcnt = 1 
      if (nfxpnt .eq. 0) then 
         fxnext = onep1*abs(xxold(nmsh)) 
      else 
         fxnext = fixpnt(ifxcnt) 
      endif 
                                                                        
!  jtkout is a counter of the number of consecutive points that         
!  have been removed from the mesh.                                     
                                                                        
      jtkout = 0 
      do 330 im = 2, ninter 
         rlold = rlen 
         rlen = xxold(im+1) - xxold(im) 
                                                                        
!  If xxold(im) is the next fixed point, it cannot be removed           
!  and so we don't test its error estimates.                            
                                                                        
         if(xxold(im) .eq. fxnext)  then 
                                                                        
            ifxcnt = ifxcnt + 1 
            if(ifxcnt .gt. nfxpnt) then 
               fxnext = onep1*abs(xxold(nmsh)) 
            else 
               fxnext = fixpnt(ifxcnt) 
            endif 
                                                                        
         elseif (irefin(im).eq.1) then 
                                                                        
!  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish    
!  to remove point im from the mesh.                                    
                                                                        
!  If we are considering removing points and jtkout = 0, this           
!  is the first point in a possible consecutive set to be removed,      
!  and we initialize phihat, which represents a maximum of              
!  certain estimates.                                                   
!  If jtkout is not zero, previous points contiguous to this            
!  point have been removed, and phihat retains its previous value.      
                                                                        
            slen = slen + rlen 
                                                                        
            if (jtkout .eq. 0) then 
                ind1 = ihcomp(im-1) 
                phihat = ermeas(ind1,im-1)/(rlold**ipow) 
            endif 
            phihat = max(phihat,                                        &
     &                 ermeas(ihcomp(im),im)/(rlen**ipow))              
            val1 = phihat*(slen**ipow) 
            if (val1 .le. phitst                                        &
     &             .and. jtkout .lt. 4) then                            
                                                                        
!  Increment the counter of removed points.                             
!  'Remove' the mesh point xxold(im) by not including it.               
                                                                        
               jtkout = jtkout+1 
               go to 330 
            endif 
!        end of logic for irefin(im) = 1.                               
         endif 
                                                                        
         jtkout = 0 
         new = new + 1 
         xx(new) = xxold(im) 
                                                                        
         if (new + irefin(im)  .gt. nmax) goto 360 
                                                                        
         if (irefin(im) .gt. 1) then 
            dx = rlen/irefin(im) 
            do 300 j = 2, irefin(im) 
              new = new + 1 
              xx(new) = xxold(im) + (j-1)*dx 
  300       continue 
         endif 
         slen = rlen 
                                                                        
         if (new  .gt. nmax) then 
                                                                        
!  If the new mesh contains too many points, branch out of the          
!  loop to try alternative strategies.                                  
                                                                        
            go to 360 
                                                                        
         elseif (new .gt. 3*ninter)  then 
                                                                        
!  Here, the new mesh does not exceed the specified maximum,            
!  but has more than 3 times as many intervals as the old mesh.         
!  Try doubling the mesh if possible.                                   
                                                                        
            call dcopy(nmsh, xxold, 1, xx, 1) 
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
            ddouble = .true. 
            return 
                                                                        
          end if 
  330 continue 
                                                                        
!  To end up here, we have processed the entire interval,               
!  and have neither exceeded the specified maximum nor                  
!  exceeded three times the number of intervals in the old              
!  mesh.  The last mesh point remains unchanged.                        
                                                                        
      new = new + 1 
      xx(new) = xxold(nmsh) 
      nmsh = new 
      maxmsh = .false. 
      if (iprint .ge. 0) write(6,905) nmsh 
      return 
                                                                        
  360 continue 
!  To reach here, the number of mesh points created at some stage       
!  of the refinement process was larger than the maximum permitted      
!  value nmax.                                                          
                                                                        
!  Check whether the mesh can safely be doubled.                        
                                                                        
      if ((2*nmsh-1) .lt. nmax) then 
                                                                        
!  Double the mesh.                                                     
         call  dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
         ddouble = .true. 
                                                                        
!  If the number of intervals is too large and the mesh cannot be       
!  doubled, increase the threshold thres by a factor of erdcid and      
!  try the selective refinement again.                                  
!  If this happens three times without success or if thres exceeds      
!  or is equal to errmax, stop.  (In this case, we know already         
!  that doubling the mesh produces too many points.)                    
                                                                        
      elseif (thres .lt. errmax .and. ithres .lt. 3) then 
         ithres = ithres + 1 
         thres = erdcid*thres 
         if(thres .gt. errmax) thres = errmax 
         call dcopy(nmsh, xxold, 1, xx, 1) 
         go to 200 
      else 
!         nmsh = 2*nmsh - 1                                             
         nmsh = nmold 
         call dcopy(nmold, xxold, 1, xx, 1) 
         maxmsh = .true. 
      endif 
      return 
                                                                        
  901 format(1h ,'selmsh.  nmsh, ipow =',2i5) 
  902 format(1h ,'im, jcomp, ermeas, normalized er',2i5,2(1pe11.3)) 
  903 format(1h ,'errmax',1pe11.3) 
  904 format(1h ,'nmest, irefin',(10i5)) 
  905 format(1h ,'selmsh.  new nmsh =',i8) 
  910 format(1h ,'ihcomp',(10i5)) 
      end subroutine                                        
                                                                        
                                                                        
       subroutine smpmsh (nmsh, nmax, xx, intref, numadd,               &
     &      nmold, xxold, maxmsh)                                       
 !module
      use algprs                                                        
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      logical maxmsh 
      dimension xx(*), xxold(*) 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
!  blas: dcopy                                                          
                                                                        
!  The routine smpmsh performs simple mesh refinement by adding         
!  points to one or three interval(s) in the region indicated           
!  by the integer intref.                                               
!  numadd gives the trial number of points to be added in each          
!  interval.                                                            
                                                                        
      if (pdebug) write(6,901) 
                                                                        
      nmold = nmsh 
      call dcopy(nmold, xx, 1, xxold, 1) 
                                                                        
!  numadd is altered if necessary so that it lies between 4 and 49      
                                                                        
      if(numadd .gt. 49) then 
         numadd = 49 
      elseif (numadd .lt. 4) then 
         numadd = 4 
      endif 
      if (pdebug) write (6,902) nmsh, intref, numadd 
                                                                        
      maxmsh = .false. 
      if (intref .eq. 1) then 
                                                                        
!  Add numadd points to the first interval if intref = 1.               
                                                                        
         nmnew = nmsh + numadd 
         if (nmnew .gt. nmax) then 
            if (iprint .ge. 0)  write(6,903) nmnew 
            maxmsh = .true. 
            return 
         endif 
                                                                        
!  Renumber the later points in reverse order.                          
                                                                        
         nint = numadd + 1 
         do 100 i = nmnew, numadd+2, -1 
            xx(i) = xx(i-numadd) 
  100    continue 
         dx = (xx(2) - xx(1))/nint 
         do 110 i = 2, nint 
            xx(i) = xx(1) + (i-1)*dx 
  110    continue 
                                                                        
      elseif (intref .eq. nmsh-1) then 
                                                                        
!  Add numadd points to the last interval if intref = nmsh-1.           
                                                                        
         nmnew = nmsh + numadd 
         if (nmnew .gt. nmax) then 
            if (iprint .ge. 0)  write(6,903) nmnew 
            maxmsh = .true. 
            return 
         endif 
         nint = numadd + 1 
         dx = (xx(nmsh) - xx(nmsh-1))/nint 
         xx(nmnew) = xx(nmsh) 
         do 200 i = nmsh, nmnew-1 
            xx(i) = xx(nmsh-1) + (i-nmsh+1)*dx 
  200    continue 
                                                                        
      else 
                                                                        
         if (numadd .gt. 9) numadd = 9 
                                                                        
!  Here, intref lies between 2 and nmsh-2.  Add numadd points to        
!  each of the three intervals intref-1, intref and intref+1.           
                                                                        
         nmnew = nmsh + 3*numadd 
         if (nmnew .gt. nmax) then 
            if (iprint .ge. 0)  write(6,903) nmnew 
            maxmsh = .true. 
            return 
         endif 
                                                                        
!  noalt is the number of points at the right end of the interval       
!  whose numerical values remain the same, but whose indices change.    
!  nochsm is the smallest index in the new ordering of one of these     
!  points.                                                              
                                                                        
         noalt = nmsh - intref - 1 
         nochsm = nmnew - noalt + 1 
                                                                        
!  Renumber the noalt unchanged points at the right end of the interval 
!  (in reverse order).                                                  
                                                                        
         j = 0 
         do 300 i = nmnew, nochsm, -1 
            xx(i) = xx(nmsh-j) 
            j = j + 1 
  300    continue 
                                                                        
!  Add numadd points to the three contiguous intervals.                 
!  The remaining points at the left end of the interval retain          
!  their original indices, and are left unchanged.                      
                                                                        
         nint = numadd + 1 
         innew = nochsm - nint 
         do 320 i = intref+1, intref-1, -1 
            xx(innew) = xx(i) 
            dx = (xx(innew + nint) - xx(innew))/nint 
            do 310 j = 1, numadd 
               xx(innew + j) = xx(innew) + j*dx 
  310       continue 
            innew = innew - nint 
  320    continue 
                                                                        
      endif 
      nmsh = nmnew 
                                                                        
      if(iprint .ge. 0)  write(6,904) nmsh 
      return 
  901 format(1h , ' smpmsh') 
  902 format(1h ,'nmsh, intref, numadd',3i6) 
  903 format(1h , ' smpmsh.  maximum points exceeded, nmnew =',i6) 
  904 format(1h ,'smpmsh, new nmsh =',i7) 
      end subroutine                                        
                                                                        
                                                                        
                                                                        
      subroutine unimsh(nmsh, aleft, aright, nfxpnt, fixpnt, xx) 
!module
      use algprs 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      integer  nmsh, nfxpnt 
      dimension fixpnt(*), xx(nmsh) 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      intrinsic max 
                                                                        
                                                                        
!  Given a left endpoint aleft, a right endpoint aright,                
!  a set of nfxpnt fixed points fixpnt(i), i = 1,...,nfxpnt,            
!  (where fixpnt(i) is different from aleft and aright for all i),      
!  and an initial target number nmsh of mesh points,                    
!  the subroutine unimsh generates a piecewise uniform mesh             
!  beginning at aleft, ending at aright, and with equally               
!  spaced points between aleft and fixpnt(1), then between              
!  fixpnt(1) and fixpnt(2), ..., and finally between                    
!  fixpnt(nfxpnt) and aright.  The final number of intervals            
!  is the maximum of nfxpnt+2 and the initial value of nmsh.            
                                                                        
!  In the simplest case when nfxpnt = 0, unimsh generates a             
!  uniform mesh with nmsh intervals in the closed interval              
!  (aleft, aright).                                                     
                                                                        
!  On exit, the integer nmsh contains the number of mesh points         
!  (which is the maximum of the initial nmsh and nfxpnt).               
!  The array xx (of dimension nmsh) contains the mesh points.           
                                                                        
      if (iprint .ge. 0) write(6,901) nmsh 
                                                                        
      if (nfxpnt .eq. 0) then 
                                                                        
!  If there are no interior fixed points, the spacing is uniform        
!  throughout the interval.  Calculate the spacing dx                   
!  and set up the xx array.                                             
                                                                        
        ninter = nmsh - 1 
                                                                        
         dx = (aright - aleft)/ninter 
         do 10 i = 1, ninter 
            xx(i) = aleft + (i-1)*dx 
   10    continue 
         xx(nmsh) = aright 
         return 
      endif 
                                                                        
!  We know that there is at least one fixed point strictly between      
!  the endpoints.                                                       
                                                                        
      if (nmsh .lt. nfxpnt+2)  nmsh = nfxpnt + 2 
      ninter = nmsh - 1 
      xx(1) = aleft 
      ileft = 1 
      xleft = aleft 
      totint = aright - aleft 
      ndif = ninter - nfxpnt 
      do 50 j = 1, nfxpnt + 1 
                                                                        
!  Deal in turn with the subintervals defined by the interval           
!  boundaries and the fixed  points.                                    
                                                                        
         if (j .lt. nfxpnt+1) then 
                                                                        
!  The j-th fixed point is xright.  Calculate where it should           
!  fall in the mesh.                                                    
                                                                        
            xright = fixpnt(j) 
            nmin = ninter*(xright-aleft)/totint + 1.5_dp 
            if (nmin .gt. ndif+j) nmin = ndif + j 
            iright = max(ileft+1, nmin) 
         else 
            xright = aright 
            iright = nmsh 
         endif 
                                                                        
!  npt is the number of equally spaced points that should               
!  lie strictly between the (j-1)-th and j-th fixed points.             
                                                                        
         xx(iright) = xright 
         npt = iright - ileft - 1 
         dx = (xright - xleft)/(npt + 1) 
         do 30 i = 1, npt 
            xx(ileft+i) = xleft + i*dx 
   30    continue 
         ileft = iright 
         xleft = xright 
   50 continue 
                                                                        
      return 
!                                                                       
  901 format (1h ,'unimsh.  nmsh =',i5) 
      end subroutine                                        
                                                                        
                                                                        
      subroutine stats(len, elem, ebigst, esecnd, summod, index) 
      implicit integer (i-n)
      implicit real(dp) (a-h, o-z) 
      integer index 
      dimension elem(len) 
                                                                        
      intrinsic abs 
                                                                        
      parameter  ( zero = 0.0_dp ) 
                                                                        
!  Given the real array elem of length len, stats calculates            
!  the following:                                                       
!      - summod, the sum of the magnitudes of the elements of elem;     
!      - ebigst (the largest element in magnitude);                     
!      - index (the index in elem of ebigst); and                       
!      - esecnd (the second largest element in magnitude, strictly      
!          less than ebigst unless both are zero).                      
                                                                        
      index = 1 
      ebigst = zero 
      esecnd = zero 
      summod = zero 
                                                                        
      do 100 i = 1, len 
         elmod = abs(elem(i)) 
         summod = summod + elmod 
         if (elmod .gt. ebigst) then 
            esecnd = ebigst 
            ebigst = elmod 
            index = i 
         elseif (elmod .gt. esecnd) then 
            esecnd = elmod 
         endif 
  100 continue 
      return 
      end subroutine                                        
                                                                        
                                                                        
      subroutine mshref(ncomp, nmsh, nlbc, ntol, ltol,                  &
     &              iorder, rhs, tmwork,                                &
     &              nmax, xx, nmold, xxold, ddouble , maxmsh,           &
     &              numbig, nummed,amg,stab_cond,stiff_cond,r4,         &
     &              nfxpnt, fixpnt,irefin,itcond,itcondmax)             
                                                                        
!module
      use algprs
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 

      dimension  ltol(ntol) 
      dimension  rhs(ncomp*nmsh), tmwork(nmsh-1) 
      dimension  xx(nmsh), xxold(nmold) 
      dimension  fixpnt(*), irefin(*),r4(*) 
      logical    ddouble , maxmsh 
      dimension  amg(*) 
      logical stab_cond, stiff_cond, forcedouble 
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
!  blas: dcopy                                                          
      logical nodouble 
      parameter (two = 2.0_dp) 
      parameter (bigfac = 10.0_dp, small = 1.0d-2 )
      integer, parameter :: numpt = 14
                                                                        
                                                                        
!  This routine performs calculations leading to a decision             
!  about how the mesh will be refined, and then refines the mesh.       
                                                                        
!  The choices for mesh refinement in this routine are either to        
!  double the mesh (i.e., divide each existing interval in half),       
!  or to add points to just a few intervals.                            
                                                                        
!  The decision is made based on two criteria:  (1) the distribution    
!  of magnitudes of components of rhs (broadly speaking, if             
!  the maximum component of rhs is much larger than the                 
!  average, points are added near the corresponding interval),          
!  and (2) the history of previous mesh refinements (if points          
!  have been added to only a few intervals already, this strategy       
!  is abandoned for the moment and the mesh is doubled).                
                                                                        
!  The decision is indicated by setting the logical flag double         
!  and (if ddouble is .false.) the integer intref.                      
!  Setting ddouble to .true. means that the mesh should be doubled      
!  (i.e., the new mesh should contain twice as many intervals).         
!  The integer intref indicates the region of the mesh where the        
!  points should be added (see the routine smpmsh), and numadd          
!  indicates how many points are to be added.                           
!  The integers nummed and numbig represent running totals,             
!  used in deciding on the mesh refinement strategy.                    
                                                                        
!  If iorder = 4, meaning that we were just performing a Newton         
!  iteration for a 4th order solution, check all elements of rhs.       
!  If iorder .gt. 4, signalling that we were trying Newton              
!  iterations for order 6 or 8, check only elements of rhs              
!  corresponding to components for which a tolerance is specified.      
                                                                        
      if (pdebug) write(6,901) nummed, numbig 
                                                                        
      nodouble = ((iorder.eq.4) .and.                                   &
     &       (stiff_cond .and. .not. stab_cond) .and. (use_c))          
!      nodouble = nodouble                                              
!     *  .or.((iorder.gt.4) .and. ( stiff_cond .and. .not. stab_cond)   
!     *   .and. (use_c))                                                
                                                                        
                                                                        
                                                                        
      forcedouble = .false. 
      if (use_c) then 
         if ( itcond .eq. itcondmax) then 
            itcond = 0 
            forcedouble = .true. 
         else 
            forcedouble = .false. 
         endif 
      endif 
                                                                        
                                                                        
      ninter = nmsh-1 
      nup = ncomp 
      if (iorder .gt. 4) nup = ntol 
                                                                        
!  Check the vector rhs for a non-negligible component                  
!  whose magnitude is significantly larger than the average.            
!  (small defines negligible, and bigfac defines significantly larger.) 
                                                                        
      do 50 ic = 1, nup 
         icmp = ic 
         if (iorder .gt. 4) icmp = ltol(ic) 
                                                                        
!  For component icmp, examine the ninter elements of rhs not           
!  corresponding to boundary conditions.                                
                                                                        
!  subroutine stats calculates rbigst and rsecnd (the first- and        
!  second-largest elements of rhs in magnitude), and the index          
!  intref of the interval in which the largest value occurs.            
!  The value sumrhs is the sum of the magnitudes of the components      
!  of rhs.                                                              
                                                                        
         indrhs = nlbc + ic 
!  Copy the elements of rhs corresponding to interior mesh              
!  points for component icmp into a single vector tmwork.               
                                                                        
         call dcopy(ninter, rhs(indrhs), ncomp, tmwork, 1) 
                                                                        
         call stats(ninter, tmwork, rbigst, rsecnd,                     &
     &                 sumrhs, intref)                                  
         tstval = bigfac*(sumrhs-rbigst)/ninter 
         if (pdebug) write(6,902) ic, tstval, rbigst, rsecnd 
         if (rbigst .ge. small .and. rbigst .ge. tstval) go to 100 
   50 continue 
                                                                        
!  If we reach this point, no interval has a significantly larger       
!  element than average.  Set counters and double the mesh.             
                                                                        
      numbig = 0 
      nummed = 0 
      ddouble = .true. 
      if (pdebug) write(6,903) 
!f the mesh if not doubled if the problem is stiff and the order is 4   
       if (nodouble .and. .not. forcedouble) then 
          call selcondmsh(ncomp, nmsh,                                  &
     &     nfxpnt, fixpnt,  nmax, xx,  irefin,                          &
     &     nmold, xxold, ddouble , maxmsh,r4,amg)                       
           ddouble = .false. 
           itcond = itcond + 1 
        else 
          call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
            itcond = 0 
        end if 
!      call dblmsh(nmsh, nmax, xx, nmold, xxold, maxmsh)                
      return 
                                                                        
!  To reach statement 100, it must be true that for some component      
!  icmp,  rbigst is non-negligible and large relative to the            
!  average.  intref indicates the region to which points may be added.  
                                                                        
!  If too many specialized refinements (adding a few points to          
!  a small number of intervals) have been made, signal that             
!  the mesh should be doubled.                                          
!  Otherwise, increment counters and add numadd points as indicated.    
                                                                        
  100 continue 
      if (rbigst .lt. two*rsecnd) nummed = nummed + 1 
      numadd = numpt 
      numbig = numbig + 1 
      ddouble = .false. 
                                                                        
      if (rbigst .le. bigfac*rsecnd .or. numbig .gt. 8) then 
         numbig = 0 
         nummed = nummed + 1 
         if (nummed .ge. 4 .and. iorder .eq. 4) then 
            ddouble = .true. 
            nummed = 0 
         elseif (nummed .ge. 8 .and. iorder .gt. 4) then 
            ddouble = .true. 
            nummed = 0 
         endif 
      endif 
                                                                        
!  Refine the mesh.                                                     
      if (pdebug) write(6,904) numbig, nummed 
                                                                        
                                                                        
      if (ddouble) then 
!f the mesh if not doubled if the problem is stiff                      
         if  (nodouble .and. .not. forcedouble)  then 
             numadd = numpt 
            call smpselcondmsh(ncomp, nmsh,                             &
     &       nfxpnt, fixpnt,  nmax, xx,  irefin,intref,numadd,          &
     &       nmold, xxold, ddouble , maxmsh,r4,amg)                     
!               itcond = itcond + 1                                     
!           call selcondmsh(ncomp, nmsh,                                
!     *        nfxpnt, fixpnt,  nmax, xx,  irefin,                      
!     *        nmold, xxold, ddouble , maxmsh,r4,amg)                   
           ddouble = .false. 
           itcond = itcond + 1 
         else 
           call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
             itcond  = 0 
         end if 
      else 
!f if the problem is stiff  we use both the old technique               
!f and the conditioning                                                 
!f                                                                      
         if (nodouble .and. .not. forcedouble)  then 
           numadd = numpt 
           call smpselcondmsh(ncomp, nmsh,                              &
     &       nfxpnt, fixpnt,  nmax, xx,  irefin,intref,numadd,          &
     &       nmold, xxold, ddouble , maxmsh,r4,amg)                     
               itcond = itcond + 1 
        elseif (forcedouble .and. use_c) then 
            ddouble = .true. 
              itcond = 0 
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
        else 
                numadd = numpt 
                call smpmsh (nmsh, nmax, xx, intref, numadd,            &
     &          nmold, xxold, maxmsh)                                   
        endif 
                                                                        
      endif 
      if (ddouble  .and. pdebug) write(6,905) 
      return 
                                                                        
  901 format(1h ,'mshref. nummed, numbig =',2i5) 
  902 format(1h ,'ic, tst, bigst, second',i5, 3(1pe11.3)) 
  903 format(1h ,'No significantly large value') 
  904 format(1h ,'numbig, nummed =',2i5) 
  905 format(1h ,'double the mesh') 
      end subroutine                                        
                                                                        
      subroutine errest (ncomp, nmsh, ntol, ltol, tol,                  &
     &   nudim, u, uold, etest, errsum, errok)                          
 !module
      use algprs                                                          
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension ltol(ntol), tol(ntol), u(nudim,nmsh),                   &
     &              uold(ncomp,nmsh), etest(ntol)                       
      logical errok 
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
      intrinsic abs, max 
                                                                        
      parameter ( one = 1.0_dp, zero = 0.0_dp ) 
                                                                        
                                                                        
!  Given current and previous solutions u and uold on the same          
!  mesh, errest calculates an error measure for each                    
!  component for which a tolerance is specified.                        
!  The error measure is the usual relative error normalized             
!  by dividing by the tolerance.  On exit, errsum is the                
!  sum of these error measures.                                         
                                                                        
!  The array etest specifies the error test to be applied to each       
!  error measure.                                                       
                                                                        
!  On exit, the logical flag errok                                      
!   -- is .false. if any of the error measures exceeds the              
!      corresponding value in the array etest                           
!   -- is .true. if all error measures are less than the                
!      corresponding values of etest.                                   
                                                                        
      if (pdebug) write(6,900) 
      errsum = zero 
      errok = .true. 
                                                                        
      do 10 im = 1, nmsh 
      do 10 it = 1, ntol 
         icmp = ltol(it) 
         er = u(icmp,im) - uold(icmp,im) 
         denom = max(one, abs(uold(icmp,im))) 
         errel = abs(er/(tol(it)*denom)) 
         if (pdebug)  write(6,901) im, it, errel, etest(it) 
         errsum = errsum + errel 
         if (errel .gt. etest(it)) errok = .false. 
   10 continue 
                                                                        
      if (pdebug) write(6,902) errsum 
      return 
  900 format(1h ,'errest') 
  901 format(1h ,2i5,2(1pe11.3)) 
  902 format(1h ,'errsum',1pe11.3) 
      end subroutine                                        
                                                                        
                                                                        
      subroutine getptq( debug, mfsrch, nout, alfmax, alfsml, alfuzz,   &
     &                   epsaf, epsag, eta, ftry, oldf, oldg,           &
     &                   rmu, tolabs, tolrel, toltny, imprvd,           &
     &                   inform, nfsrch, alfa, alfbst, fbest,           &
     &                   braktd, crampd, extrap,vset,wset,nsamea,nsameb,&
     &                   a, b, fa, factor, fv, fw, xtry, xv, xw )       
                                                                        
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      logical            debug, imprvd 
      logical            braktd, crampd, extrap, vset, wset 
      integer            mfsrch, nout, inform, nfsrch, nsamea, nsameb 
!                                                                       
!  *********************************************************************
!  getptq  is a step-length algorithm for minimizing a function of one  
!  variable.  it will be called repeatedly by a search routine whose    
!  purpose is to estimate a point  alfa = alfbst  that minimizes some   
!  function  f(alfa)  over the closed interval (0, alfmax).             
!                                                                       
!  getptq  requires the function  f(alfa)  (but not its gradient)       
!  to be evaluated at various points within the interval.  new          
!  step-length estimates are computed using quadratic interpolation with
!  safeguards.                                                          
!                                                                       
!  reverse communication is used to allow the calling program to        
!  evaluate  f.  some of the parameters must be set or tested           
!  by the calling program.  the remainder would ordinarily be local     
!  variables.                                                           
!                                                                       
!                                                                       
!  input parameters (relevant to the calling program)                   
!  --------------------------------------------------                   
!                                                                       
!  debug         specifies whether detailed output is wanted.           
!                                                                       
!  inform        must be nonzero on the first entry (e.g., -1).         
!                it will be altered by  getptq  for later entries.      
!                                                                       
!  mfsrch        is an upper limit on the number of times  getptq  is   
!                to be entered consecutively with  inform = 0           
!                (following an initial entry with  inform lt 0).        
!                                                                       
!  nout          is the file number to be used for printed output       
!                if debug is true.                                      
!                                                                       
!  alfa          is the first estimate of the step length.  alfa  is    
!                subsequently altered by  getptq  (see below).          
!                                                                       
!  alfmax        is the upper limit of the interval to be searched.     
!                                                                       
!  alfsml        is intended to prevent inefficiency when the optimum   
!                step is very small, for cases where the calling        
!                program would prefer to re-define  f(alfa).  alfsml is 
!                allowed to be zero. early termination will occur if    
!                getptq  determines that the optimum step lies          
!                somewhere in the interval  (0, alfsml)  (but not if    
!                alfmax .le. alfsml).                                   
!                                                                       
!  epsaf         is an estimate of the absolute precision in the        
!                computed values of  f.                                 
!                                                                       
!  eta           controls the accuracy of the search.  it must lie      
!                in the range   0.0  le  eta  lt  1.0.  decreasing      
!                eta  tends to increase the accuracy of the search.     
!                                                                       
!  oldf          is the value of  f(0).                                 
!                                                                       
!  oldg          is an estimate of the gradient of  f  at  alfa = 0.    
!                it should be non-positive.                             
!                                                                       
!  rmu           controls what is meant by a significant decrease in  f.
!                the final  f(alfbst)  should lie on or below the line  
!                      l(alfa)  =  oldf + alfa*rmu*oldg.                
!                rmu  should be in the open interval (0, 0.5).          
!                the value  rmu = 1.0d-4  is good for most purposes.    
!                                                                       
!  tolabs,tolrel define a function  tol(alfa) = tolrel*alfa + tolabs    
!                such that if  f  has already been evaluated at step    
!                alfa,  then it will not be evaluated at any point      
!                closer than  tol(alfa).                                
!                these values may be reduced by  getptq  if they seem   
!                to be too large.                                       
!                                                                       
!  toltny        is the smallest value that  tolabs  is allowed to be   
!                reduced to.                                            
!                                                                       
!                                                                       
!  output parameters (relevant to the calling program)                  
!  ---------------------------------------------------                  
!                                                                       
!  imprvd        is true if the previous step  alfa  was the best       
!                point so far.  any related quantities (e.g., arrays)   
!                should be saved by the calling program before paying   
!                attention to  inform.                                  
!                                                                       
!  inform = 0    means the calling program should evaluate              
!                           ftry = f(alfa)                              
!                for the new trial step  alfa,  and then re-enter       
!                getptq.                                                
!                                                                       
!  inform = 1    means the search has terminated successfully           
!                with a step  alfbst  that is less than the             
!                upper bound  alfmax.                                   
!                                                                       
!  inform = 2    means the search has terminated successfully           
!                with a step  alfbst  that is equal to the              
!                upper bound  alfmax.                                   
!                                                                       
!  inform = 3    means that the search failed to find a point of        
!                sufficient decrease in  mfsrch  functions, but an      
!                improved point was found.                              
!                                                                       
!  inform = 4    means  alfmax  is so small that a search should        
!                not have been done.                                    
!                                                                       
!  inform = 5    means that the search was terminated prematurely       
!                because of the value of  alfsml  (see above).          
!                                                                       
!  inform = 6    means the search has failed to find a useful step.  if 
!                the subroutine for the function and gradient has been  
!                programmed correctly, this will usually occur if the   
!                minimum lies very close to  alfa = 0  or the gradient  
!                is not sufficiently accurate.                          
!                                                                       
!  inform = 7    means that the value of  g(0) was positive on entry.   
!                                                                       
!  alfa          is the step at which the next function value must be   
!                computed.                                              
!                                                                       
!  alfbst        should be accepted by the calling program as the       
!                required step-length estimate, whenever  getptq        
!                returns  inform = 1,  2  or  3.                        
!                                                                       
!  fbest         will be the corresponding value of  f.                 
!                                                                       
!                                                                       
!  the following parameters retain information between entries          
!  -----------------------------------------------------------          
!                                                                       
!  alfuzz        is such that, if the final  alfa  lies in the interval 
!                (0,alfuzz)  and  abs( f(alfa)-oldf ) le epsaf,  alfa   
!                cannot be guaranteed to be a point of sufficient       
!                decrease.                                              
!                                                                       
!  braktd        is false if  f  has not been evaluated at the far end  
!                of the interval of uncertainty.  in this case, the     
!                point  b  will be at  alfmax + tol(alfmax).            
!                                                                       
!  crampd        is true if  alfmax  is very small (le tolabs).         
!                if the search fails, this indicates that a zero        
!                step should be taken.                                  
!                                                                       
!  extrap        is true if alfbst has moved at least once and  xv      
!                lies outside the interval of uncertainty.  in this     
!                case, extra safeguards are applied to allow for        
!                instability in the polynomial fit.                     
!                                                                       
!  vset          records whether a third-best point has been            
!                determined.                                            
!                                                                       
!  wset          records whether a second-best point has been           
!                determined.  it will always be true by the             
!                time the convergence test is applied (label 300).      
!                                                                       
!  nsamea        is the number of consecutive times that the left-hand  
!                end of the interval of uncertainty has remained the    
!                same.                                                  
!                                                                       
!  nsameb        similarly for the right-hand end.                      
!                                                                       
!  a, b, alfbst  define the current interval of uncertainty.            
!                the required minimum lies somewhere within the         
!                closed interval  (alfbst + a, alfbst + b).             
!                                                                       
!  alfbst        is the best point so far.  it is strictly within the   
!                the interval of uncertainty except when it lies at the 
!                left-hand end when  alfbst  has not been moved.        
!                hence we have    a le 0,   b gt 0.                     
!                                                                       
!  fbest         is the value of  f  at the point  alfbst.              
!                                                                       
!  fa            is the value of  f  at the point  alfbst + a.          
!                                                                       
!  factor        controls the rate at which extrapolated estimates of   
!                alfa  may expand into the interval of uncertainty.     
!                factor is not used if the minimum has been bracketed   
!                (i.e., when the variable  braktd  is true).            
!                                                                       
!  fv, fw        are the values of  f  at the points  alfbst + xv,      
!                alfbst + xw.  they are not defined until  vset         
!                or  wset  (respectively) is true.                      
!                                                                       
!  ftry          is the value of  f  at the new point  alfbst + xtry.   
!                                                                       
!  xtry          is the trial point within the shifted interval (a, b). 
!                the new trial function value must be computed at the   
!                point  alfa  =  alfbst + xtry.                         
!                                                                       
!  xv            is such that  alfbst + xv  is the third-best point.    
!                it is not defined until  vset  is true.                
!                                                                       
!  xw            is such that  alfbst + xw  is the second-best point.   
!                it is not defined until  wset  is true.                
!                in some cases,  xw  will replace a previous  xw  that  
!                has a lower function but has just been excluded from   
!                the interval of uncertainty.                           
!                                                                       
!                                                                       
!  systems optimization laboratory, stanford university, california.    
!  original version february 1982.  rev. may 1983.                      
!  *********************************************************************
!                                                                       
      logical            closef, conv1, conv2, conv3, convrg 
      logical            moved, sigdec, xinxw 
      data               zero, point1, half/ 0.0_dp,  0.1_dp, 0.5_dp/ 
      data                one,  two,  five   / 1.0_dp, 2.0_dp,  5.0_dp/ 
      data                ten, eleven      /10.0_dp, 11.0_dp        / 
!                                                                       
!                                                                       
!  local variables                                                      
!  ---------------                                                      
!                                                                       
!  closef        is true if the worst function  fv  is within  epsaf    
!                of  fbest  (up or down).                               
!                                                                       
!  convrg        will be set to true if at least one of the convergence 
!                conditions holds at  alfbst.                           
!                                                                       
!  moved         is true if a better point has been found (alfbst gt 0).
!                                                                       
!  sigdec        says whether  fbest  represents a significant decrease 
!                in the function, compared to the initial value  oldf.  
!                                                                       
!  xinxw         is true if  xtry  is in  (xw,0)  or  (0,xw).           
!  ---------------------------------------------------------------------
!                                                                       
      imprvd = .false. 
      if (inform .ne. -1) go to 100 
!                                                                       
!  ---------------------------------------------------------------------
!  first entry.  initialize various quantities, check input data and    
!  prepare to evaluate the function at the initial step  alfa.          
!  ---------------------------------------------------------------------
      nfsrch = 0 
      alfbst = zero 
      fbest  = oldf 
      if (oldg   .gt.      zero) go to 970 
      if (oldg   .ge. (- epsag)) go to 960 
      if (alfmax .le.    toltny) go to 940 
!                                                                       
      braktd = .false. 
      crampd = alfmax .le. tolabs 
      extrap = .false. 
      vset   = .false. 
      wset   = .false. 
      nsamea = 0 
      nsameb = 0 
      alfuzz = two*epsaf/(rmu*abs( oldg )) 
      a      = zero 
      b      = alfmax + (tolrel*alfmax + tolabs) 
      fa     = oldf 
      factor = five 
      tol    = tolabs 
      xtry   = alfa 
      if (debug) write (nout, 1000) alfmax, oldf, oldg, tolabs,         &
     &   alfuzz, epsaf, epsag, tolrel, crampd                           
      go to 800 
!                                                                       
!  ---------------------------------------------------------------------
!  subsequent entries.                                                  
!  the function has just been evaluated at  alfa = alfbst + xtry,       
!  giving  ftry.                                                        
!  ---------------------------------------------------------------------
  100 nsamea = nsamea + 1 
      nsameb = nsameb + 1 
      xtry   = alfa - alfbst 
      moved  = alfbst .gt. zero 
!                                                                       
!  check if  xtry  is in the interval  (xw,0)  or  (0,xw).              
!                                                                       
      xinxw  = .false. 
      if (wset) xinxw =       zero .lt. xtry  .and.  xtry .le. xw       &
     &                  .or.    xw .le. xtry  .and.  xtry .lt. zero     
!                                                                       
!  see if the new step is better.                                       
!                                                                       
      deltaf = ftry   - oldf 
      ctry   = deltaf - alfa*rmu*oldg 
      if (alfa .le. alfuzz) sigdec = deltaf .le. (- epsaf) 
      if (alfa .gt. alfuzz) sigdec = ctry   .le.    epsaf 
      imprvd = sigdec  .and.  ( ftry - fbest ) .le. (- epsaf) 
!                                                                       
      if (debug) write (nout, 1100) alfa, ftry, ctry 
      if (.not. imprvd) go to 130 
!                                                                       
!  we seem to have an improvement.  the new point becomes the           
!  origin and other points are shifted accordingly.                     
!                                                                       
      if (.not. wset) go to 110 
      xv     = xw - xtry 
      fv     = fw 
      vset   = .true. 
  110 xw     = zero - xtry 
      fw     = fbest 
      wset   = .true. 
      fbest  = ftry 
      alfbst = alfa 
      a      =    a - xtry 
      b      =    b - xtry 
      moved  = .true. 
      extrap = .not. xinxw 
!                                                                       
!  decrease the length of the interval of uncertainty.                  
!                                                                       
      if (xtry .lt. zero) go to 120 
      a      = xw 
      fa     = fw 
      nsamea = 0 
      go to 300 
  120 b      = xw 
      nsameb = 0 
      braktd = .true. 
      go to 300 
!                                                                       
!  the new function value is no better than the best point found so far.
!  the point  xtry  must be a new end point of the interval of          
!  uncertainty.                                                         
!                                                                       
  130 if (xtry .ge. zero) go to 140 
      a      = xtry 
      fa     = ftry 
      nsamea = 0 
      go to 150 
  140 b      = xtry 
      nsameb = 0 
      braktd = .true. 
!                                                                       
!  the origin remains unchanged but  xtry  may qualify as  xw.          
!                                                                       
  150 if (.not. wset)   go to 160 
      if ((ftry - fw) .gt. epsaf) go to 170 
      xv     = xw 
      fv     = fw 
      vset   = .true. 
  160 xw     = xtry 
      fw     = ftry 
      wset   = .true. 
      if (moved) extrap = xinxw 
      go to 300 
!                                                                       
!  ftry  is no better than  fbest  or  fw.  if the best point has not   
!  been moved, there must be more than one minimum.                     
!                                                                       
  170 if (moved) go to 175 
      xw     = xtry 
      fw     = ftry 
      go to 300 
!                                                                       
!  ftry  is no better than  fbest  or  fw,  but  xtry  may become  xv.  
!  extrap  has the value set in the previous entry.                     
!                                                                       
  175 if (.not. vset) go to 180 
      if ((ftry - fv) .gt. epsaf  .and.  extrap) go to 300 
  180 if (xinxw) go to 190 
      xv     = xtry 
      fv     = ftry 
      vset   = .true. 
      go to 300 
  190 if (vset) xw = xv 
      if (vset) fw = fv 
      xv     = xtry 
      fv     = ftry 
      vset   = .true. 
!                                                                       
!  ---------------------------------------------------------------------
!  check the termination criteria.                                      
!  ---------------------------------------------------------------------
  300 tol    = tolrel*alfbst + tolabs 
      deltaf = fbest  - oldf 
                                                                        
      cbest  = deltaf - alfbst*rmu*oldg 
      if (alfbst .le. alfuzz) sigdec = deltaf .le. (- epsaf) 
      if (alfbst .gt. alfuzz) sigdec = cbest  .le.    epsaf 
      closef = .false. 
      if (vset) closef = abs( fbest - fv ) .le. epsaf 
!                                                                       
      conv1  = max( abs( a ), b )  .le.  (tol + tol) 
                                                                        
!  conv2 changed by mhw, 20 sept 1992, to allow it to be                
!  satified for any significant decrease in f                           
!      conv2  =  moved  .and.  sigdec                                   
!     *                 .and.  abs( fa - fbest )  .le.  a*eta*oldg      
      conv2 = moved .and. sigdec 
      conv3  = closef  .and.  (sigdec  .or.                             &
     &                        (.not. moved)  .and.  (b .le. alfuzz))    
      convrg = conv1  .or.  conv2  .or.  conv3 
!                                                                       
      atrue  = alfbst + a 
      btrue  = alfbst + b 
      alfaw  = alfbst + xw 
      gap    = b - a 
      if (debug) write (nout, 1200) atrue, btrue, gap, tol,             &
     &   nsamea, nsameb, braktd, closef, imprvd, conv1, conv2, conv3,   &
     &   extrap, alfbst, fbest, cbest, alfaw, fw                        
      if (vset) alfav  = alfbst + xv 
      if (debug  .and.  vset) write (nout, 1300) alfav, fv 
      if (convrg  .and.  moved) go to 910 
!                                                                       
!  exit if the step is too small.                                       
!                                                                       
      if (btrue   .lt.  alfsml) go to 950 
                                                                        
      if (nfsrch  .ge.  mfsrch) go to 930 
      if (.not. convrg) go to 400 
!                                                                       
!  a better point has not yet been found (the step  xw  is no better    
!  than step  zero).  check that the change in  f  is consistent with a 
!  perturbation in  x  of  tol, the estimate of the minimum spacing     
!  constant.  if the change in  f  is larger than  epsaf,  the value    
!  of  tol  is reduced.                                                 
!                                                                       
      tol    = xw/ten 
      tolabs = tol 
      if (abs(fw - oldf) .gt. epsaf  .and.  tol .gt. toltny) go to 400 
      if (crampd) go to 940 
      go to 960 
!                                                                       
!  ---------------------------------------------------------------------
!  proceed with the computation of a trial step length.                 
!  the choices are...                                                   
!  1. parabolic fit using function values only.                         
!  2. damped parabolic fit if the regular fit appears to be             
!     consistently over-estimating the distance to the minimum.         
!  3. bisection, geometric bisection, or a step of  tol  if the         
!     parabolic fit is unsatisfactory.                                  
!  ---------------------------------------------------------------------
  400 xmidpt = half*(a + b) 
      q      = zero 
      s      = zero 
!                                                                       
!  ---------------------------------------------------------------------
!  fit a parabola.                                                      
!  ---------------------------------------------------------------------
!                                                                       
!  check if there are two or three points for the parabolic fit.        
!                                                                       
      gw = (fw - fbest)/xw 
      if (vset  .and.  moved) go to 450 
!                                                                       
!  only two points available.  use  fbest,  fw  and the derivative      
!  oldg.                                                                
!                                                                       
      if (.not. moved) s = oldg 
      if (      moved) s = oldg - two*gw 
      q = two*(oldg - gw) 
      if (debug) write (nout, 2100) 
      go to 600 
!                                                                       
!  three points available.  use  fbest,  fw  and  fv.                   
!                                                                       
  450 gv = (fv - fbest)/xv 
      s  = gv - (xv/xw)*gw 
      q  = two*(gv - gw) 
      if (debug) write (nout, 2200) 
!                                                                       
!  ---------------------------------------------------------------------
!  construct an artificial interval  (artifa, artifb)  in which the     
!  new estimate of the step length must lie.  set a default value of    
!  xtry  that will be used if the polynomial fit is rejected.  in the   
!  following, the interval  (a,b)  is considered the sum of two         
!  intervals of lengths  dtry  and  daux, with common end point at the  
!  best point (zero).  dtry  is the length of the interval into which   
!  the default  xtry  will be placed and  endpnt  denotes its non-zero  
!  end point.  the magnitude of  xtry  is computed so that the exponents
!  of  dtry  and  daux  are approximately bisected.                     
!  ---------------------------------------------------------------------
  600 artifa = a 
      artifb = b 
      if (braktd) go to 610 
!                                                                       
!  the minimum has not been bracketed.  set an artificial upper bound   
!  by expanding the interval  xw  by a suitable factor.                 
!                                                                       
      xtry   = - factor*xw 
      artifb =   xtry 
      if (alfbst + xtry .lt. alfmax) factor = five*factor 
      go to 700 
!                                                                       
!  the minimum has been bracketed.                                      
!  if the gradient at the origin is being used for the                  
!  polynomial fit, the default  xtry  is one tenth of  xw.              
!                                                                       
  610 if (vset  .and.  moved) go to 620 
      xtry   = xw/ten 
      if (debug) write (nout, 2400) xtry 
      go to 700 
!                                                                       
!  three points exist in the interval of uncertainty.  check whether    
!  the points are configured for an extrapolation or interpolation.     
!                                                                       
  620 if (extrap) go to 660 
!                                                                       
!  if the interpolation appears to be consistently over-estimating the  
!  distance to the minimum,  damp the interpolation step.               
!                                                                       
      if (nsamea .lt. 3  .and.  nsameb .lt. 3) go to 630 
      factor = factor / five 
      s      = factor * s 
      go to 640 
  630 factor = one 
!                                                                       
!  the points are configured for an interpolation.  the artificial      
!  interval will be just  (a,b).   set  endpnt  so that  xtry           
!  lies in the larger of the intervals  (a,0)  and  (0,b).              
!                                                                       
  640 if (xmidpt .lt. zero) endpnt = a 
      if (xmidpt .gt. zero) endpnt = b 
!                                                                       
!  if a bound has remained the same for three iterations, set  endpnt   
!  so that  xtry  is likely to replace the offending bound.             
!                                                                       
      if (nsamea .ge. 3) endpnt = a 
      if (nsameb .ge. 3) endpnt = b 
      go to 680 
!                                                                       
!  the points are configured for an extrapolation.                      
!                                                                       
  660 if (xw .lt. zero) endpnt = b 
      if (xw .gt. zero) endpnt = a 
!                                                                       
!  compute the default value of  xtry.                                  
!                                                                       
  680 dtry = abs( endpnt ) 
      daux = gap - dtry 
      if (daux .ge. dtry)   xtry = five*dtry*(point1 + dtry/daux)/eleven 
      if (daux .lt. dtry)   xtry = half*sqrt( daux )*sqrt( dtry ) 
      if (endpnt .lt. zero) xtry = - xtry 
      if (debug) write (nout, 2500) xtry, daux, dtry 
!                                                                       
!  if the points are configured for an extrapolation set the artificial 
!  bounds so that the artificial interval lies strictly within  (a,b).  
!  if the polynomial fit is rejected,  xtry  will remain at the relevant
!  artificial bound.                                                    
!                                                                       
      if (extrap  .and.  xtry .le. zero) artifa = xtry 
      if (extrap  .and.  xtry .gt. zero) artifb = xtry 
!                                                                       
!  ---------------------------------------------------------------------
!  the polynomial fits give  (s/q)*xw  as the new step.                 
!  reject this step if it lies outside  (artifa, artifb).               
!  ---------------------------------------------------------------------
  700 if (q .eq. zero) go to 800 
      if (q .lt. zero) s = - s 
      if (q .lt. zero) q = - q 
      if (s*xw .lt. q*artifa   .or.   s*xw .gt. q*artifb) go to 800 
!                                                                       
!  accept the polynomial fit.                                           
!                                                                       
      xtry = zero 
      if (abs( s*xw ) .ge. q*tol) xtry = (s/q)*xw 
      if (debug) write (nout, 2600) xtry 
!                                                                       
!  ---------------------------------------------------------------------
!  test for  xtry  being larger than  alfmax  or too close to  a  or  b.
!  ---------------------------------------------------------------------
  800 if (braktd) go to 810 
!                                                                       
!  if the step is close to or larger than  alfmax,  replace it by       
!  alfmax  (to force evaluation of the function at the boundary).       
!                                                                       
      alfa   = alfbst + xtry 
      if (alfmax - alfa .gt. (tolrel*alfmax + tolabs)) go to 810 
      braktd = .true. 
      xtry   = alfmax - alfbst 
      alfa   = alfmax 
      go to 900 
!                                                                       
!  otherwise, the function must not be evaluated at a point too close   
!  to  a  or  b.  (it has already been evaluated at both those points.) 
!                                                                       
  810 xmidpt = half*(a + b) 
      if (xtry .gt. a + tol  .and.  xtry .lt. b - tol) go to 820 
      if (xmidpt .gt. zero) xtry =   tol 
      if (xmidpt .le. zero) xtry = - tol 
!                                                                       
!                                                                       
!  f  must not be calculated too close to  alfbst.                      
!                                                                       
  820 if (abs( xtry ) .lt. tol  .and.  xmidpt .lt. zero) xtry = - tol 
      if (abs( xtry ) .lt. tol  .and.  xmidpt .ge. zero) xtry =   tol 
      alfa   = alfbst + xtry 
!                                                                       
!  ---------------------------------------------------------------------
!  exit.                                                                
!  ---------------------------------------------------------------------
!                                                                       
!  new function value required.                                         
!                                                                       
  900 inform = 0 
      go to 990 
!                                                                       
!  convergence test satisfied.                                          
!                                                                       
  910 inform = 1 
      if (alfa .eq. alfmax) inform = 2 
      go to 990 
!                                                                       
!  mfsrch  function evaluations without sufficient decrease, but an     
!  improved point was found.                                            
!                                                                       
  930 if (.not. moved) go to 960 
      inform = 3 
      go to 990 
!                                                                       
!  zero step (alfmax too small).                                        
!                                                                       
  940 inform = 4 
      go to 990 
!                                                                       
!  premature termination.  the step is smaller than  alfsml.            
!                                                                       
  950 inform = 5 
      go to 990 
!                                                                       
!  zero step (a sufficiently better point could not be found).          
!                                                                       
  960 inform = 6 
      go to 990 
!                                                                       
!  zero step (positive gradient at the starting point).                 
!                                                                       
  970 inform = 7 
!                                                                       
!  exit.                                                                
!                                                                       
  990 if (debug) write (nout, 3000) 
      return 
!                                                                       
 1000 format(/ 31h alfmax  oldf    oldg    tolabs, 1p2e22.14, 1p2e16.8  &
     &       / 31h alfuzz  epsaf   epsag   tolrel, 1p2e22.14, 1p2e16.8  &
     &       / 31h crampd                        ,  l6)                 
 1100 format(/ 31h alfa    ftry    ctry          , 1p2e22.14, 1pe16.8) 
 1200 format(/ 31h a       b       b - a   tol   , 1p2e22.14, 1p2e16.8  &
     &       / 31h nsamea  nsameb  braktd  closef, 2i3, 2l6             &
     &       / 31h imprvd  convrg  extrap        ,  l6, 3x, 3l1, l6     &
     &       / 31h alfbst  fbest   cbest         , 1p2e22.14, 1pe16.8   &
     &       / 31h alfaw   fw                    , 1p2e22.14)           
 1300 format(  31h alfav   fv                    , 1p2e22.14 /) 
                                                                        
 2100 format(30h parabolic fit,    two points.) 
 2200 format(30h parabolic fit,  three points.) 
 2400 format(31h exponent reduced.  trial point, 1p1e22.14) 
 2500 format(31h geo. bisection. xtry,daux,dtry, 1p3e22.14) 
 2600 format(31h polynomial fit accepted.  xtry, 1p1e22.14) 
 3000 format(53h ---------------------------------------------------- /) 
!                                                                       
!  end of getptq                                                        
      end subroutine                                        
                                                                        
                                                                        
      subroutine interp(ncomp, nmsh, xx, nudim, u,                      &
     &                    nmold, xxold, uold)                           
 !module
      use algprs                                                             
      implicit integer (i-n)
      implicit real(dp) (a-h, o-z) 
      dimension xx(*), u(nudim,*), xxold(*), uold(ncomp,*) 
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
! blas: dcopy                                                           
                                                                        
      parameter (zero = 0.0_dp) 
                                                                        
!  interp performs piecewise linear interpolation of the old            
!  solution uold at the nmold old mesh points xxold onto the nmsh       
!  new mesh points xx, producing an interpolated solution u.            
!  Note that no assumption is made that the new mesh has                
!  more points than the old, nor that the new and old mesh              
!  points are related in a specific way (except that their first        
!  and last points are identical).                                      
                                                                        
      if (pdebug) write(6,900) 
                                                                        
!  By construction, xx(1) = xxold(1).  Copy the first ncomp             
!  components of uold into those of u.                                  
                                                                        
      call dcopy(ncomp, uold(1,1), 1, u(1,1), 1) 
                                                                        
      i = 2 
      do 100 im = 2, nmsh-1 
                                                                        
   50    continue 
         if (i .gt. nmold) return 
                                                                        
!  Check whether the im-th point in the new mesh lies strictly          
!  to the right of, or to the left of (or exactly on) the               
!  i-th point in the old mesh.                                          
                                                                        
                                                                        
         if (xx(im) .gt. xxold(i)) then 
            i = i + 1 
            go to 50 
         else 
            xdif = xxold(i) - xx(im) 
            if (xdif .eq. zero) then 
                                                                        
!  xx(im) and xxold(i) are identical.                                   
                                                                        
               call dcopy(ncomp, uold(1,i), 1, u(1,im), 1) 
               i = i + 1 
            else 
               xint = xxold(i) - xxold(i-1) 
               xrat = xdif/xint 
               do 70 k = 1, ncomp 
                  u(k,im) = uold(k,i) + xrat*(uold(k,i-1)-uold(k,i)) 
   70          continue 
            endif 
         endif 
                                                                        
  100 continue 
      call dcopy(ncomp, uold(1,nmold), 1, u(1,nmsh), 1) 
      return 
  900 format(1h ,'interp') 
      end subroutine                                        
                                                                        
                                                                        
      subroutine rerrvl( ncomp, nmsh, nudim, u, usvrex, ntol, ltol,     &
     &          rerr, remax, itlmx, adjrer )                            
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension ltol(*) 
      dimension u(nudim, *), usvrex(ncomp, *) 
      dimension rerr(ncomp, *) 
      logical adjrer 
                                                                        
      intrinsic abs, max 
                                                                        
      parameter ( one = 1.0_dp, zero = 0.0_dp ) 
                                                                        
!  rerrvl is used in considering Richardson extrapolation.              
!  The two solutions u and usvrex have a special relationship:          
!  u corresponds to a doubled mesh, with twice as many                  
!  intervals, where each interval is half the size of that for          
!  usvrex's mesh.   nmsh is the number of mesh points in the            
!  mesh for u.                                                          
                                                                        
!  remax is the maximum relative error, and itlmx is the                
!  index of the tolerance for which the maximum occurs.                 
                                                                        
!  The array rerr contains the absolute value of the difference         
!  between u and usvrex at the common mesh points, but is defined       
!  only for components for which an error tolerance is specified.       
                                                                        
!  The logical variable adjrer is true on entry if the values in        
!  rerr are to be adjusted for later use in selective mesh refinement.  
                                                                        
                                                                        
      itlmx = 1 
      remax = zero 
      nmold = 1 + (nmsh-1)/2 
      do 100 it = 1, ntol 
         icmp = ltol(it) 
         imnew = 1 
         do 50 im = 1, nmold 
            rerr(icmp, im) = abs(usvrex(icmp,im) - u(icmp,imnew)) 
            denom = max(one, abs(usvrex(icmp,im))) 
            rerel = rerr(icmp,im)/denom 
            if (rerel .gt. remax) then 
               remax = rerel 
               itlmx = it 
            endif 
            imnew = imnew + 2 
   50    continue 
  100 continue 
                                                                        
      if (adjrer) then 
                                                                        
!  Adjust the rerr array if it may be used later in selective           
!  mesh refinement.                                                     
                                                                        
         do 150 it = 1, ntol 
            icmp = ltol(it) 
            do 150 im = 1, nmold - 1 
               rerr(icmp,im) = max(rerr(icmp,im),                       &
     &                              rerr(icmp, im+1))                   
  150    continue 
      endif 
                                                                        
      return 
      end subroutine                                        
                                                                        
      function dasum(n,dx,incx) 
!#                                                                      
!#     takes the sum of the absolute values.                            
!#     jack dongarra, linpack, 3/11/78.                                 
!#     modified 3/93 to return if incx .le. 0.                          
!#     modified 12/3/93, array(1) declarations changed to array(*)      
!#                                                                      
      real(dp) dx(*),dtemp ,dasum
      integer i,incx,m,mp1,n,nincx 
!#                                                                      
      dasum = 0.0d0 
      dtemp = 0.0d0 
      if( n .le. 0 .or. incx .le. 0 )return 
      if(incx.eq.1)go to 20 
!#                                                                      
!#        code for increment not equal to 1                             
!#                                                                      
      nincx = n*incx 
      do 10 i = 1,nincx,incx 
        dtemp = dtemp + abs(dx(i)) 
   10 continue 
      dasum = dtemp 
      return 
!#                                                                      
!#        code for increment equal to 1                                 
!#                                                                      
!#                                                                      
!#        clean-up loop                                                 
!#                                                                      
   20 m = mod(n,6) 
      if( m .eq. 0 ) go to 40 
      do 30 i = 1,m 
        dtemp = dtemp + abs(dx(i)) 
   30 continue 
      if( n .lt. 6 ) go to 60 
   40 mp1 = m + 1 
      do 50 i = mp1,n,6 
        dtemp = dtemp + abs(dx(i)) + abs(dx(i + 1)) + abs(dx(i + 2)) &
     &  + abs(dx(i + 3)) + abs(dx(i + 4)) + abs(dx(i + 5))           
   50 continue 
   60 dasum = dtemp 
      return 
      end function dasum
                                                                        
                                                                        
                                                                        
                                                                        
      subroutine selconderrmsh(ncomp, nmsh, ntol, ltol, tol,            &
     &     nfxpnt, fixpnt, ipow, nmax,                                  &
     &     xx, nudim, u, ermeas, irefin, ihcomp, nmold, xxold,          &
     &     ermx, ddouble , maxmsh,r4,amg,stab_cond)                     
 
!module
      use algprs                                                             
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
                                                                        
      dimension  ltol(ntol), tol(ntol), fixpnt(*) 
      dimension  xx(*), u(nudim, *), ermeas(ncomp,*) 
      dimension  irefin(nmsh-1), ihcomp(nmsh-1) 
      dimension  xxold(*), ermx(*), amg(*), r4(nmsh) 
      logical    ddouble , maxmsh, stab_cond 
                                                                        
      intrinsic abs, max, int 
                                                                        
!  blas: dcopy                                                          
!  real(dp) log                                                
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      parameter  ( zero = 0.0_dp, one = 1.0_dp, onep1 = 1.1_dp ) 
      parameter  ( erdcid = 5.0_dp ) 
      parameter  ( phitst = 0.1_dp ) 
                                                                        
      logical first, add 
      save    first, rlndec 
      data    first / .true. / 
                                                                        
!  The routine selconderrmsh performs selective mesh refinement, dependi
!  on the error measure ermeas and the monitor function based on the    
!  conditioning.                                                        
                                                                        
      if (first) then 
         first = .false. 
         rlndec = log(erdcid) 
      endif 
                                                                        
      maxmsh = .false. 
                                                                        
      if (pdebug) write(6,901) nmsh, ipow 
                                                                        
      frcpow = one/ipow 
      ddouble = .false. 
      nmold = nmsh 
      ninter = nmsh - 1 
                                                                        
!  Copy the current mesh into the xxold array.                          
                                                                        
                                                                        
      call dcopy(nmold, xx, 1, xxold, 1) 
      ithres = 0 
      thres = one 
                                                                        
!  On input, the array ermeas represents some error measure defined     
!  over the components and mesh intervals (not mesh points).            
!  It is normalized in the following loop with respect to the           
!  tolerance array and the current solution.                            
!  The value errmax gives the maximum normalized error.                 
                                                                        
      errmax = zero 
      do 120 im = 1, ninter 
         ermx(im) = zero 
         do 110 it = 1, ntol 
            jcomp = ltol(it) 
            denom = tol(it)*max(one, abs(u(jcomp,im))) 
            ems = ermeas(jcomp,im) 
            ermeas(jcomp,im) = abs(ems)/denom 
!            if (pdebug .and. ermeas(jcomp,im) .ge. thres)              
!     *             write(6,902) im,jcomp,ems,ermeas(jcomp,im)          
            err = ermeas(jcomp, im) 
            if (err .ge. ermx(im)) then 
                ermx(im) = err 
                ihcomp(im) = jcomp 
            endif 
  110    continue 
         errmax = max(ermx(im), errmax) 
  120 continue 
                                                                        
!       if ((.not. stab_cond) .and. errmax .ge. 1.0d6  ) then           
!f  only the conditioning                                               
!            call  selcondmsh(ncomp, nmsh,                              
!     *         nfxpnt, fixpnt,  nmax, xx,  irefin,                     
!     *         nmold, xxold, ddouble , maxmsh,r4,amg)                  
!                                                                       
!      else                                                             
!f   the conditining and the error                                      
      call moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg) 
                                                                        
                                                                        
      if (pdebug) write(6,903) errmax 
                                                                        
      if (errmax .gt. zero .and. errmax .le. erdcid) then 
                                                                        
!  If errmax > 0 and .le. erdcid, find the smallest integer exponent ii 
!  such that (erdcid**ii)*errmax > erdcid.                              
                                                                        
         if(errmax .gt. one) then 
            ii = 1 
            decii = erdcid 
         else 
            ilg = -log(errmax)/rlndec 
            ii = 2 + ilg 
            decii = erdcid**ii 
         endif 
                                                                        
!  Multiply error measures by erdcid**ii.                               
                                                                        
         errmax = decii*errmax 
         do 140 im = 1, ninter 
            ermx(im) = decii*ermx(im) 
            do 140 it = 1, ntol 
               jcomp = ltol(it) 
               ermeas(jcomp,im) = decii*ermeas(jcomp,im) 
  140    continue 
      endif 
                                                                        
  200 continue 
                                                                        
!  For each interval im,  the integer irefin(im) is calculated          
!  based on an equidistrbution procedure involving the                  
!  threshold value thres.  If irefin(im) > 1, we add                    
!  points to interval im.  If irefin(im) = 1, we check whether          
!  point im can be removed from the mesh.                               
                                                                        
!  nmest is a lower bound on the number of points in the new mesh.      
!  We do not know in advance how many points will be removed,           
!  so nmest is computed by assuming all eligible points are removed.    
                                                                        
      nmest = nmsh 
      do 220 im = 1, ninter 
         if (ermx(im).ge.thres) then 
            irefin(im) = int(ermx(im)**frcpow) + 1 
            nmest = nmest + irefin(im) - 1 
         else 
            irefin(im) = 1 
            nmest = nmest - 1 
         endif 
  220 continue 
                                                                        
       if(nptcond.ge.4 ) then 
         add = .false. 
       nptcond = nptcond/2 
       do 221 im = 1, ninter-1 
         if  (max(r4(im), r4(im+1)) .gt. fatt_r1r3) then 
            if (.not. add) then 
              irefin(im) = max(nptcond, irefin(im)) 
!              nmest = nmest + nptcond - 1                              
            endif 
              irefin(im+1) = max(nptcond, irefin(im+1)) 
!                nmest = nmest + nptcond - 1                            
              add = .true. 
         else 
            irefin(im) = max(1, irefin(im)) 
            irefin(im+1) = max(1, irefin(im+1)) 
!            nmest = nmest - 1                                          
            add = .false. 
         endif 
  221    continue 
                                                                        
!        do 221 im = 1, ninter                                          
!         if ( r4(im) .gt. fatt_r1r3) then                              
!              irefin(im) = max(nptcond, irefin(im))                    
!c              nmest = nmest + nptcond - 1                             
!         else                                                          
!              irefin(im) = max(1, irefin(im))                          
!c              nmest = nmest - 1                                       
!         endif                                                         
! 221     continue                                                      
        end if 
                                                                        
      if (pdebug) write(6,904) nmest, (irefin(i), i=1,ninter) 
                                                                        
      if (nmest .gt. nmax) then 
                                                                        
         go to 360 
                                                                        
      elseif (nmest-1 .gt. 3*ninter) then 
                                                                        
         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
         ddouble = .true. 
         return 
      endif 
                                                                        
!  It appears that we can perform the desired selective mesh            
!  refinement.                                                          
                                                                        
!  Now begin running through the mesh, adding and possibly deleting     
!  points as indicated by the irefin array.                             
                                                                        
!  The integer new is a count of the number of intervals in             
!  the tentative mesh being generated by the refinement strategy.       
                                                                        
      new = 1 
                                                                        
!  The first interval is treated as a special case, since xx(1)         
!  always remains in the mesh, and cannot be a fixed point.             
                                                                        
      rlen = xxold(2) - xx(1) 
      slen = rlen 
      if (irefin(1).gt.1) then 
         dx = rlen/irefin(1) 
         do 230 j = 2, irefin(1) 
            new = new + 1 
            xx(new) = xx(1) + (j-1)*dx 
  230    continue 
      endif 
                                                                        
!  The fixed points specified by the fixpnt array cannot be             
!  removed from the mesh.  The value fxnext indicates the 'next'        
!  fixed point.  When no further fixed points remain to be processed    
!  (or if nfxpnt = 0), fxnext is set to a value strictly larger than    
!  the last mesh point, so that no mesh point can equal fxnext.         
!  This way we need to compare only the values of xxold(i)              
!  and fxnext.                                                          
                                                                        
      ifxcnt = 1 
      if (nfxpnt .eq. 0) then 
         fxnext = onep1*abs(xxold(nmsh)) 
      else 
         fxnext = fixpnt(ifxcnt) 
      endif 
                                                                        
!  jtkout is a counter of the number of consecutive points that         
!  have been removed from the mesh.                                     
                                                                        
      jtkout = 0 
      do 330 im = 2, ninter 
         rlold = rlen 
         rlen = xxold(im+1) - xxold(im) 
                                                                        
!  If xxold(im) is the next fixed point, it cannot be removed           
!  and so we don't test its error estimates.                            
                                                                        
         if(xxold(im) .eq. fxnext)  then 
                                                                        
            ifxcnt = ifxcnt + 1 
            if(ifxcnt .gt. nfxpnt) then 
               fxnext = onep1*abs(xxold(nmsh)) 
            else 
               fxnext = fixpnt(ifxcnt) 
            endif 
                                                                        
         elseif (irefin(im).eq.1) then 
                                                                        
!  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish    
!  to remove point im from the mesh.                                    
                                                                        
!  If we are considering removing points and jtkout = 0, this           
!  is the first point in a possible consecutive set to be removed,      
!  and we initialize phihat, which represents a maximum of              
!  certain estimates.                                                   
!  If jtkout is not zero, previous points contiguous to this            
!  point have been removed, and phihat retains its previous value.      
                                                                        
            slen = slen + rlen 
                                                                        
            if (jtkout .eq. 0) then 
                ind1 = ihcomp(im-1) 
                phihat = ermeas(ind1,im-1)/(rlold**ipow) 
            endif 
            phihat = max(phihat,                                        &
     &                 ermeas(ihcomp(im),im)/(rlen**ipow))              
            val1 = phihat*(slen**ipow) 
            if (val1 .le. phitst                                        &
     &             .and. jtkout .lt. 4                                  &
     &             .and. r4(im) .le. r3) then                           
                                                                        
!  Increment the counter of removed points.                             
!  'Remove' the mesh point xxold(im) by not including it.               
                                                                        
               jtkout = jtkout+1 
               go to 330 
            endif 
!        end of logic for irefin(im) = 1.                               
         endif 
                                                                        
         jtkout = 0 
         new = new + 1 
         xx(new) = xxold(im) 
         if (irefin(im) .gt. 1) then 
            dx = rlen/irefin(im) 
            do 300 j = 2, irefin(im) 
              new = new + 1 
              xx(new) = xxold(im) + (j-1)*dx 
  300       continue 
         endif 
         slen = rlen 
                                                                        
         if (new .gt. nmax) then 
                                                                        
!  If the new mesh contains too many points, branch out of the          
!  loop to try alternative strategies.                                  
                                                                        
            go to 360 
                                                                        
         elseif (new .gt. 3*ninter)  then 
                                                                        
!  Here, the new mesh does not exceed the specified maximum,            
!  but has more than 3 times as many intervals as the old mesh.         
!  Try doubling the mesh if possible.                                   
            nmsh = nmold 
            call dcopy(nmsh, xxold, 1, xx, 1) 
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
            ddouble = .true. 
            return 
                                                                        
         endif 
                                                                        
  330 continue 
                                                                        
!  To end up here, we have processed the entire interval,               
!  and have neither exceeded the specified maximum nor                  
!  exceeded three times the number of intervals in the old              
!  mesh.  The last mesh point remains unchanged.                        
                                                                        
      new = new + 1 
      xx(new) = xxold(nmsh) 
      nmsh = new 
      maxmsh = .false. 
      if (iprint .ge. 0) write(6,905) nmsh 
      return 
                                                                        
  360 continue 
!  To reach here, the number of mesh points created at some stage       
!  of the refinement process was larger than the maximum permitted      
!  value nmax.                                                          
                                                                        
!  Check whether the mesh can safely be doubled.                        
                                                                        
      if ((2*nmsh-1) .lt. nmax) then 
                                                                        
!  Double the mesh.                                                     
         call  dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
         ddouble = .true. 
                                                                        
!  If the number of intervals is too large and the mesh cannot be       
!  doubled, increase the threshold thres by a factor of erdcid and      
!  try the selective refinement again.                                  
!  If this happens three times without success or if thres exceeds      
!  or is equal to errmax, stop.  (In this case, we know already         
!  that doubling the mesh produces too many points.)                    
                                                                        
      elseif (thres .lt. errmax .and. ithres .lt. 3) then 
         ithres = ithres + 1 
         thres = erdcid*thres 
         if(thres .gt. errmax) thres = errmax 
         call dcopy(nmsh, xxold, 1, xx, 1) 
         go to 200 
      else 
!         nmsh = 2*nmsh - 1                                             
         nmsh = nmold 
         call dcopy(nmsh, xxold, 1, xx, 1) 
         maxmsh = .true. 
      endif 
                                                                        
!      end if                                                           
!   endif use the conditioning and the error                            
      return 
                                                                        
  901 format(1h ,'selconderrmsh.  nmsh, ipow =',2i5) 
  902 format(1h ,'im, jcomp, ermeas, normalized er',2i5,2(1pe11.3)) 
  903 format(1h ,'errmax',1pe11.3) 
  904 format(1h ,'nmest, irefin',(10i5)) 
  905 format(1h ,'selconderrmsh.  new nmsh =',i8) 
  910 format(1h ,'ihcomp',(10i5)) 
      end subroutine                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      subroutine selcondmsh(ncomp, nmsh,                                &
     &     nfxpnt, fixpnt,  nmax, xx,  irefin,                          &
     &     nmold, xxold, ddouble , maxmsh, r4, amg)                     
 !module
      use algprs                                                             
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
                                                                        
      dimension  fixpnt(*) 
      dimension  xx(*) 
      dimension  irefin(*) 
      dimension  xxold(*),  amg(*), r4(nmsh) 
      logical    ddouble , maxmsh, add 
                                                                        
      intrinsic abs, max, int 
                                                                        
!  blas: dcopy                                                          
!  real(dp) log                                                
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      parameter  ( zero = 0.0_dp, one = 1.0_dp, onep1 = 1.1_dp ) 
      parameter  ( erdcid = 5.0_dp ) 
                                                                        
                                                                        
                                                                        
                                                                        
!  The routine selcondmsh performs selective mesh refinement, depending 
!  on the monitor function based on the conditioning parameters.        
                                                                        
                                                                        
      maxmsh = .false. 
                                                                        
                                                                        
      if (pdebug) write(6,901) nmsh 
                                                                        
      ddouble = .false. 
      nmold = nmsh 
      ninter = nmsh - 1 
                                                                        
!  Copy the current mesh into the xxold array.                          
                                                                        
                                                                        
      call dcopy(nmold, xx, 1, xxold, 1) 
                                                                        
      ithres = 0 
      thres = one 
                                                                        
!  On input, the array amg represents the conditioning vector defined   
!  over the components and mesh intervals (not mesh points).            
!  We compute the monitor function and the related parameters in the    
!  followinf function                                                   
                                                                        
      call moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg) 
                                                                        
                                                                        
                                                                        
                                                                        
!  For each interval im,  the integer irefin(im) is calculated          
!  based on an equidistrbution procedure involving the                  
!  threshold value thres.  If irefin(im) > 1, we add                    
!  points to interval im.  If irefin(im) = 1, we check whether          
!  point im can be removed from the mesh.                               
                                                                        
!  nmest is a lower bound on the number of points in the new mesh.      
!  We do not know in advance how many points will be removed,           
!  so nmest is computed by assuming all eligible points are removed.    
                                                                        
      add = .false. 
      nmest = nmsh 
      do 220 im = 1, ninter-1 
         if  (max(r4(im), r4(im+1)) .gt. fatt_r1r3) then 
            if (.not. add) then 
              irefin(im) = nptcond 
              nmest = nmest + nptcond - 1 
            endif 
              irefin(im+1) = nptcond 
                nmest = nmest + nptcond - 1 
              add = .true. 
         else 
            irefin(im) = 1 
            irefin(im+1) = 1 
            nmest = nmest - 1 
            add = .false. 
         endif 
                                                                        
  220 continue 
                                                                        
                                                                        
      if (nmest .gt. nmax) then 
                                                                        
         go to 360 
                                                                        
      endif 
                                                                        
!  It appears that we can perform the desired selective mesh            
!  refinement.                                                          
                                                                        
!  Now begin running through the mesh, adding and possibly deleting     
!  points as indicated by the irefin array.                             
                                                                        
!  The integer new is a count of the number of intervals in             
!  the tentative mesh being generated by the refinement strategy.       
                                                                        
      new = 1 
                                                                        
!  The first interval is treated as a special case, since xx(1)         
!  always remains in the mesh, and cannot be a fixed point.             
                                                                        
      rlen = xxold(2) - xx(1) 
      slen = rlen 
      if (irefin(1).gt.1) then 
         dx = rlen/irefin(1) 
         do 230 j = 2, irefin(1) 
            new = new + 1 
            xx(new) = xx(1) + (j-1)*dx 
  230    continue 
      endif 
                                                                        
!  The fixed points specified by the fixpnt array cannot be             
!  removed from the mesh.  The value fxnext indicates the 'next'        
!  fixed point.  When no further fixed points remain to be processed    
!  (or if nfxpnt = 0), fxnext is set to a value strictly larger than    
!  the last mesh point, so that no mesh point can equal fxnext.         
!  This way we need to compare only the values of xxold(i)              
!  and fxnext.                                                          
                                                                        
      ifxcnt = 1 
      if (nfxpnt .eq. 0) then 
         fxnext = onep1*abs(xxold(nmsh)) 
      else 
         fxnext = fixpnt(ifxcnt) 
      endif 
                                                                        
!  jtkout is a counter of the number of consecutive points that         
!  have been removed from the mesh.                                     
                                                                        
      jtkout = 0 
      do 330 im = 2, ninter 
                                                                        
         rlold = rlen 
         rlen = xxold(im+1) - xxold(im) 
                                                                        
!  If xxold(im) is the next fixed point, it cannot be removed           
!  and so we don't test its error estimates.                            
                                                                        
         if(xxold(im) .eq. fxnext)  then 
                                                                        
            ifxcnt = ifxcnt + 1 
            if(ifxcnt .gt. nfxpnt) then 
               fxnext = onep1*abs(xxold(nmsh)) 
            else 
               fxnext = fixpnt(ifxcnt) 
            endif 
                                                                        
         elseif (irefin(im).eq.1) then 
                                                                        
!  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish    
!  to remove point im from the mesh.                                    
                                                                        
!  If we are considering removing points and jtkout = 0, this           
!  is the first point in a possible consecutive set to be removed,      
!  and we initialize phihat, which represents a maximum of              
!  certain estimates.                                                   
!  If jtkout is not zero, previous points contiguous to this            
!  point have been removed, and phihat retains its previous value.      
                                                                        
            slen = slen + rlen 
                                                                        
            if ( jtkout .lt. 1                                          &
     &             .and. r4(im) .le. fatt_r3 ) then                     
                                                                        
!  Increment the counter of removed points.                             
!  'Remove' the mesh point xxold(im) by not including it.               
                                                                        
               jtkout = jtkout+1 
               go to 330 
            endif 
!        end of logic for irefin(im) = 1.                               
         endif 
         if (new+irefin(im) .gt. nmax) then 
                                                                        
!  If the new mesh contains too many points, branch out of the          
!  loop to try alternative strategies.                                  
                                                                        
            go to 360 
         end if 
                                                                        
         jtkout = 0 
         new = new + 1 
         xx(new) = xxold(im) 
         if (irefin(im) .gt. 1) then 
            dx = rlen/irefin(im) 
            do 300 j = 2, irefin(im) 
              new = new + 1 
              xx(new) = xxold(im) + (j-1)*dx 
  300       continue 
         endif 
         slen = rlen 
                                                                        
         if (new .gt. nmax) then 
                                                                        
!  If the new mesh contains too many points, branch out of the          
!  loop to try alternative strategies.                                  
                                                                        
            go to 360 
                                                                        
         elseif (new .gt. 3*ninter)  then 
                                                                        
!  Here, the new mesh does not exceed the specified maximum,            
!  but has more than 3 times as many intervals as the old mesh.         
!  Try doubling the mesh if possible.                                   
            nmsh = nmold 
            call dcopy(nmsh, xxold, 1, xx, 1) 
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
            ddouble = .true. 
            return 
                                                                        
                                                                        
         endif 
                                                                        
                                                                        
                                                                        
  330 continue 
                                                                        
!  To end up here, we have processed the entire interval,               
!  and have neither exceeded the specified maximum nor                  
!  exceeded three times the number of intervals in the old              
!  mesh.  The last mesh point remains unchanged.                        
                                                                        
      new = new + 1 
      xx(new) = xxold(nmsh) 
      nmsh = new 
      maxmsh = .false. 
      if (iprint .ge. 0) write(6,905) nmsh 
      return 
                                                                        
  360 continue 
!  To reach here, the number of mesh points created at some stage       
!  of the refinement process was larger than the maximum permitted      
!  value nmax.                                                          
                                                                        
      nmsh = nmold 
      call dcopy(nmsh, xxold, 1, xx, 1) 
      maxmsh = .true. 
                                                                        
                                                                        
      return 
                                                                        
                                                                        
                                                                        
  901 format(1h ,'selcondmsh.  nmsh, ipow =',2i5) 
  902 format(1h ,'im, jcomp, ermeas, normalized er',2i5,2(1pe11.3)) 
  903 format(1h ,'errmax',1pe11.3) 
  904 format(1h ,'nmest, irefin',(10i5)) 
  905 format(1h ,'selcondmsh.  new nmsh =',i8) 
  910 format(1h ,'ihcomp',(10i5)) 
      end subroutine                                        
                                                                        
                                                                        
                                                                        
      subroutine smpselcondmsh(ncomp, nmsh,                             &
     &     nfxpnt, fixpnt,  nmax, xx,  irefin, intref, numadd,          &
     &     nmold, xxold, ddouble , maxmsh, r4, amg)                     
 !module
      use algprs                                                           
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
                                                                        
      dimension  fixpnt(*) 
      dimension  xx(*) 
      dimension  irefin(*) 
      dimension  xxold(*),  amg(*), r4(*) 
      logical    ddouble , maxmsh, add 
                                                                        
      intrinsic abs, max, int 
                                                                        
!  blas: dcopy                                                          
!  real(dp) log                                                
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      parameter  ( zero = 0.0_dp, one = 1.0_dp,onep1 = 1.1_dp) 
      parameter  ( erdcid = 5.0_dp ) 
                                                                        
                                                                        
                                                                        
                                                                        
!  The routine smpselcondmsh performs selective mesh refinement, by addi
!  points to one or three interval(s) in the region indicated           
!  by the integer intref (numadd gives the trial number of points to    
!  added in each  interval)  and by adding and removing point           
!  using  the monitor function based on the conditioning parameters.    
                                                                        
                                                                        
                                                                        
      maxmsh = .false. 
                                                                        
                                                                        
      if (pdebug) write(6,901) nmsh 
                                                                        
      ddouble = .false. 
      nmold = nmsh 
      ninter = nmsh - 1 
                                                                        
!  Copy the current mesh into the xxold array.                          
                                                                        
                                                                        
         call dcopy(nmold, xx, 1, xxold, 1) 
                                                                        
      ithres = 0 
      thres = one 
                                                                        
!  On input, the array amg represents the conditioning vector defined   
!  over the components and mesh intervals (not mesh points).            
                                                                        
      call moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg) 
                                                                        
                                                                        
                                                                        
!  For each interval im,  the integer irefin(im) is calculated          
!  based on an equidistrbution procedure involving the                  
!  threshold value thres.  If irefin(im) > 1, we add                    
!  points to interval im.  If irefin(im) = 1, we check whether          
!  point im can be removed from the mesh.                               
                                                                        
!  nmest is a lower bound on the number of points in the new mesh.      
!  We do not know in advance how many points will be removed,           
!  so nmest is computed by assuming all eligible points are removed.    
                                                                        
                                                                        
      nmest = nmsh 
                                                                        
                                                                        
      do 220 im = 1, ninter 
         if ( r4(im) .gt. fatt_r1r3) then 
              irefin(im) = nptcond 
              nmest = nmest + nptcond - 1 
         else 
              irefin(im) = 1 
              nmest = nmest - 1 
         endif 
                                                                        
  220 continue 
                                                                        
      if(numadd .gt. 49) then 
         numadd = 49 
      elseif (numadd .lt. 4) then 
         numadd = 4 
      endif 
                                                                        
        if ( intref .eq. 1) then 
            irefin(1) = max(numadd+1,irefin(1)) 
            nmest = nmest + numadd - 1 
        elseif (intref .eq. ninter) then 
            irefin(ninter) = max(numadd+1,irefin(ninter)) 
            nmest = nmest + numadd - 1 
        else 
          if(numadd .gt. 9) then 
             numadd = 9 
           elseif (numadd .lt. 4) then 
              numadd = 4 
            endif 
            irefin(intref-1) = max(numadd+1 , irefin(intref-1)) 
            irefin(intref)   = max(numadd+1 , irefin(intref)) 
            irefin(intref+1) = max(numadd+1 , irefin(intref+1)) 
            nmest = nmest + 3*numadd - 1 
        end if 
                                                                        
                                                                        
                                                                        
      if (nmest .gt. nmax) then 
                                                                        
         go to 360 
                                                                        
      endif 
                                                                        
!  It appears that we can perform the desired selective mesh            
!  refinement.                                                          
                                                                        
!  Now begin running through the mesh, adding and possibly deleting     
!  points as indicated by the irefin array.                             
                                                                        
!  The integer new is a count of the number of intervals in             
!  the tentative mesh being generated by the refinement strategy.       
                                                                        
      new = 1 
                                                                        
!  The first interval is treated as a special case, since xx(1)         
!  always remains in the mesh, and cannot be a fixed point.             
                                                                        
      rlen = xxold(2) - xx(1) 
      slen = rlen 
      if (irefin(1).gt.1) then 
         dx = rlen/irefin(1) 
         do 230 j = 2, irefin(1) 
            new = new + 1 
            xx(new) = xx(1) + (j-1)*dx 
  230    continue 
      endif 
                                                                        
!  The fixed points specified by the fixpnt array cannot be             
!  removed from the mesh.  The value fxnext indicates the 'next'        
!  fixed point.  When no further fixed points remain to be processed    
!  (or if nfxpnt = 0), fxnext is set to a value strictly larger than    
!  the last mesh point, so that no mesh point can equal fxnext.         
!  This way we need to compare only the values of xxold(i)              
!  and fxnext.                                                          
                                                                        
      ifxcnt = 1 
      if (nfxpnt .eq. 0) then 
         fxnext = onep1*abs(xxold(nmsh)) 
      else 
         fxnext = fixpnt(ifxcnt) 
      endif 
                                                                        
!  jtkout is a counter of the number of consecutive points that         
!  have been removed from the mesh.                                     
                                                                        
      jtkout = 0 
      do 330 im = 2, ninter 
                                                                        
         rlold = rlen 
         rlen = xxold(im+1) - xxold(im) 
                                                                        
!  If xxold(im) is the next fixed point, it cannot be removed           
!  and so we don't test its error estimates.                            
                                                                        
         if(xxold(im) .eq. fxnext)  then 
                                                                        
            ifxcnt = ifxcnt + 1 
            if(ifxcnt .gt. nfxpnt) then 
               fxnext = onep1*abs(xxold(nmsh)) 
            else 
               fxnext = fixpnt(ifxcnt) 
            endif 
                                                                        
         elseif (irefin(im).eq.1) then 
                                                                        
!  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish    
!  to remove point im from the mesh.                                    
                                                                        
!  If we are considering removing points and jtkout = 0, this           
!  is the first point in a possible consecutive set to be removed,      
!  and we initialize phihat, which represents a maximum of              
!  certain estimates.                                                   
!  If jtkout is not zero, previous points contiguous to this            
!  point have been removed, and phihat retains its previous value.      
                                                                        
            slen = slen + rlen 
                                                                        
            if ( jtkout .lt. 1                                          &
     &             .and. r4(im) .le. fatt_r3) then                      
                                                                        
!  Increment the counter of removed points.                             
!  'Remove' the mesh point xxold(im) by not including it.               
                                                                        
               jtkout = jtkout+1 
               go to 330 
            endif 
!        end of logic for irefin(im) = 1                                
         endif 
          if (new + irefin(im) .gt. nmax) then 
                                                                        
!  If the new mesh contains too many points, branch out of the          
!  loop to try alternative strategies.                                  
                                                                        
         go to 360 
         end if 
         jtkout = 0 
         new = new + 1 
         xx(new) = xxold(im) 
         if (irefin(im) .gt. 1) then 
            dx = rlen/irefin(im) 
            do 300 j = 2, irefin(im) 
              new = new + 1 
              xx(new) = xxold(im) + (j-1)*dx 
  300       continue 
         endif 
         slen = rlen 
                                                                        
         if (new .gt. nmax) then 
                                                                        
!  If the new mesh contains too many points, branch out of the          
!  loop to try alternative strategies.                                  
                                                                        
            go to 360 
                                                                        
         elseif (new .gt. 3*ninter)  then 
                                                                        
!  Here, the new mesh does not exceed the specified maximum,            
!  but has more than 3 times as many intervals as the old mesh.         
!  Try doubling the mesh if possible.                                   
            if (iprint .eq. 1) write(6,*) 'smpselcondmsh' 
            nmsh = nmold 
            call dcopy(nmsh, xxold, 1, xx, 1) 
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
            ddouble = .true. 
            return 
                                                                        
         endif 
                                                                        
                                                                        
                                                                        
  330 continue 
                                                                        
!  To end up here, we have processed the entire interval,               
!  and have neither exceeded the specified maximum nor                  
!  exceeded three times the number of intervals in the old              
!  mesh.  The last mesh point remains unchanged.                        
                                                                        
      new = new + 1 
      xx(new) = xxold(nmsh) 
      nmsh = new 
      maxmsh = .false. 
      if (iprint .ge. 0) write(6,905) nmsh 
      return 
                                                                        
  360 continue 
!  To reach here, the number of mesh points created at some stage       
!  of the refinement process was larger than the maximum permitted      
!  value nmax.                                                          
                                                                        
      nmsh = nmold 
      call dcopy(nmsh, xxold, 1, xx, 1) 
      maxmsh = .true. 
                                                                        
                                                                        
      return 
                                                                        
                                                                        
                                                                        
  901 format(1h ,'smpselcondmsh.  nmsh, ipow =',2i5) 
  902 format(1h ,'im, jcomp, ermeas, normalized er',2i5,2(1pe11.3)) 
  903 format(1h ,'errmax',1pe11.3) 
  904 format(1h ,'nmest, irefin',(10i5)) 
  905 format(1h ,'smpselcondmsh.  new nmsh =',i8) 
  910 format(1h ,'ihcomp',(10i5)) 
      end subroutine                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      subroutine moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,         &
     &                        nptcond,r4,amg)                           
 !module
      use algprs                                                          
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
                                                                        
      dimension  xx(*) 
      dimension  amg(*), r4(*) 
                                                                        
      intrinsic abs, max, int 
                                                                        
!  blas: dcopy                                                          
!  real(dp) log                                                
                                                                        
!      logical pdebug, use_c, comp_c 
!      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c 
                                                                        
      parameter  ( zero = 0.0_dp, one = 1.0_dp ) 
                                                                        
! the function moncond compute the monitor function based on the        
! conditioning parameters and the factor used to perform the mesh select
                                                                        
       do i=1,nmsh-1 
          r4(i) = (xx(i+1)-xx(i))*abs(amg(i+1)- amg(i)) 
       end do 
                                                                        
                                                                        
       r2 = r4(1) 
       do i=2,nmsh-1 
          r2 = r2 + r4(i) 
       end do 
       do i=1,nmsh-1 
          r4(i) = r4(i)+(xx(i+1)-xx(i))*(r2/(xx(nmsh)-xx(1)))*0.08d0 
       end do 
       r1 = r4(1) 
       do i=2,nmsh-1 
          r1 = max(r1,r4(i)) 
       end do 
       do i=1,nmsh-1 
         r4(i) = r4(i)/r1 
       end do 
       r1 = one 
                                                                        
        r2 = r4(1) 
        r1m = r4(1) 
       do i=2,nmsh-1 
          r1m = min(r1m,r4(i)) 
          r2 = r2 + r4(i) 
       end do 
                                                                        
       r3 = r2/(nmsh-1) 
       fatt_r3  = r3*0.00001d0 
       fatt_r1r3= max(r3,r1*0.5d0) 
                                                                        
       nptm = 0 
       nptr = 0 
       do i=1,nmsh-1 
           if (r4(i) .ge. fatt_r1r3)  nptm = nptm + 1 
           if (r4(i) .le. fatt_r3)  nptr = nptr + 1 
        enddo 
                                                                        
        if (nptm .le. 1) then 
           nptcond =  14 
        elseif (nptm .le. 2) then 
           nptcond =  10 
        elseif (nptm .le. 4) then 
           nptcond =  8 
        elseif (nptm .le. 8) then 
            nptcond = 6 
        elseif (nptm .le. nmsh/20 ) then 
            nptcond = 4 
        else 
            nptcond = 2 
        endif 
                                                                        
                                                                        
       if (iprint .eq. 1) write(6,901)r1,r3,fatt_r1r3,nptcond,nptm,nptr 
                                                                        
                                                                        
  901 format(1h ,'moncondmsh.', (1pe11.3), 2(1pe11.3), 3i10) 
                                                                        
      end subroutine                                        
                                                                        
!==================================================================     
!                                                                       
        REAL(DP) FUNCTION ABDNRM(NBLOKS,NTOP,NBOT,NOVRLP,       &
     &                            NRWBLK,NCLBLK,TOP,A,BOT)              
                                                                        
!******************************************************************     
!                                                                       
!       ABDNRM IS USED IN CONJUNCTION WITH DONEST TO COMPUTE THE        
!       CONDITION NUMBER OF AN ALMOST BLOCK DIAGONAL MATRIX LIKE        
!       THE ONES HANDLED BY COLROW. [SEE COMMENTS IN COLROW, CRDCMP]    
                                                                        
!               *****  AUXILIARY PROGRAMS  *****                        
                                                                        
!       THE BLAS ARE REQUIRED BY ABDNRM.                                
!       SPECIFICALLY, DASUM IS USED.                                    
!                                                                       
                                                                        
!******************************************************************     
                                                                        
        INTEGER NBLOKS,NTOP,NBOT,NOVRLP,NRWBLK,NCLBLK 
        REAL(DP) TOP(NTOP,*),A(NRWBLK,NCLBLK,*),BOT(NBOT,*) 
        INTEGER J,K 
        REAL(DP) MAX,TEMP !,DASUM
        TEMP = 0.0D0 
!                                                                       
!       FIRST, GO OVER THE COLUMNS OF TOP AND THE FIRST BLOCK:          
!                                                                       
        DO 10 J = 1,NOVRLP 
           TEMP = MAX(TEMP, DASUM(NTOP,TOP(1,J),1) +                    &
     &                      DASUM(NRWBLK,A(1,J,1),1))                   
   10         CONTINUE 
        DO 40 K = 1,NBLOKS-1 
!                                                                       
!          IN EACH BLOCK:                                               
!                                                                       
!          FIRST, THE COLUMNS FROM THE KTH BLOCK ALONE:                 
!                                                                       
           DO 20 J = NOVRLP+1,NRWBLK 
              TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,K),1)) 
   20               CONTINUE 
!                                                                       
!          NOW, TH COLUMNS WHICH INTERSECT BOTH THE KTH AND             
!          (K+1)ST BLOCKS.                                              
!                                                                       
           DO 30 J = NRWBLK+1,NCLBLK 
              TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,K),1) +               &
     &                         DASUM(NRWBLK,A(1,J-NRWBLK,K+1),1))       
   30               CONTINUE 
   40                  CONTINUE 
!                                                                       
!       FINALLY, THE COLUMNS OF THE LAST BLOCK WHICH DO NOT OVERLAP     
!       WITH ANYTHING.                                                  
!                                                                       
        DO 50 J = NOVRLP+1,NRWBLK 
           TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,NBLOKS),1)) 
   50         CONTINUE 
!                                                                       
!       AND THOSE COLUMNS OVERLAPPING WITH BOTH THE LAST BLOCK AND THE  
!       BOTTOM BLOCK.                                                   
!                                                                       
        DO 60 J = NRWBLK+1,NCLBLK 
           TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,NBLOKS),1) +             &
     &                      DASUM(NBOT,BOT(1,J-NRWBLK),1))              
   60         CONTINUE 
        ABDNRM = TEMP 
        RETURN 
      end function
                                                                        
      SUBROUTINE DONEST (N, V, X, ISGN, EST, KASE) 
      INTEGER N, ISGN(N), KASE 
      REAL(DP) V(N), X(N), EST 
                                                                        
!                                                                       
!     DONEST ESTIMATES THE 1-NORM OF A SQUARE, REAL(DP) MATRIX  
!     REVERSE COMMUNICATION IS USED FOR EVALUATING                      
!     MATRIX-VECTOR PRODUCTS.                                           
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX.  N .GE. 1.                    
!                                                                       
!        ISGN    INTEGER(N)                                             
!                USED AS WORKSPACE.                                     
!                                                                       
!        KASE    INTEGER                                                
!                = 0.                                                   
!                                                                       
!     ON INTERMEDIATE RETURNS                                           
!                                                                       
!        KASE    = 1 OR 2.                                              
!                                                                       
!        X       REAL(DP)(N)                                    
!                MUST BE OVERWRITTEN BY                                 
!                                                                       
!                     A*X,             IF KASE=1,                       
!                     TRANSPOSE(A)*X,  IF KASE=2,                       
!                                                                       
!                AND DONEST MUST BE RE-CALLED, WITH ALL THE OTHER       
!                PARAMETERS UNCHANGED.                                  
!                                                                       
!     ON FINAL RETURN                                                   
!                                                                       
!        KASE    = 0.                                                   
!                                                                       
!        EST     REAL(DP)                                       
!                CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).      
!                                                                       
!        V       REAL(DP)(N)                                    
!                = A*W,   WHERE  EST = NORM(V)/NORM(W)                  
!                         (W  IS NOT RETURNED).                         
!                                                                       
!     THIS VERSION DATED MARCH 16, 1988.                                
!     NICK HIGHAM, UNIVERSITY OF MANCHESTER.                            
!                                                                       
!     MODIFIED FOR REAL(DP) ON JUNE 11, 1996.                   
!                                                                       
!     REFERENCE                                                         
!     N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING                   
!     THE ONE-NORM OF A REAL OR COMPLEX MATRIX, WITH APPLICATIONS       
!     TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,      
!     UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.            
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!     BLAS     IDAMAX, DASUM, DCOPY                                     
!     GENERIC  ABS, NINT, FLOAT, SIGN                                   
!                                                                       
        INTRINSIC FLOAT 
!        REAL(DP) FLOAT                                         
                                                                        
        INTRINSIC ABS 
!        REAL(DP) ABS 
                                                                        
        INTRINSIC SIGN 
!        REAL(DP) SIGN 
                                                                        
                                                                        
!      REAL(DP) DASUM 
!      INTEGER IDAMAX 
                                                                        
      INTEGER ITMAX 
      PARAMETER (ITMAX = 5) 
      REAL(DP) ZERO,ONE,TWO 
      PARAMETER (ZERO = 0.0D0) 
      PARAMETER (ONE = 1.0D0) 
      PARAMETER (TWO = 2.0D0) 
!                                                                       
!     INTERNAL VARIABLES                                                
      INTEGER I, ITER, J, JLAST, JUMP 
      REAL(DP) ALTSGN, ESTOLD, TEMP 
!                                                                       
      SAVE 
!                                                                       
      IF (KASE .EQ. 0) THEN 
         DO 10,I = 1,N 
            X(I) = ONE/FLOAT(N) 
   10           CONTINUE 
         KASE = 1 
         JUMP = 1 
         RETURN 
      ENDIF 
!                                                                       
      GOTO (100, 200, 300, 400, 500) JUMP 
!                                                                       
!     ................ ENTRY   (JUMP = 1)                               
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.                  
!                                                                       
  100  CONTINUE 
      IF (N .EQ. 1) THEN 
         V(1) = X(1) 
         EST = ABS(V(1)) 
!        ... QUIT                                                       
         GOTO 510 
      ENDIF 
      EST = DASUM(N,X,1) 
!                                                                       
      DO 110,I = 1,N 
         X(I) = SIGN(ONE,X(I)) 
         ISGN(I) = NINT(X(I)) 
  110     CONTINUE 
      KASE = 2 
      JUMP = 2 
      RETURN 
!                                                                       
!     ................ ENTRY   (JUMP = 2)                               
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.       
!                                                                       
  200  CONTINUE 
      J = IDAMAX(N,X,1) 
      ITER = 2 
!                                                                       
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.                             
!                                                                       
  220  CONTINUE 
      DO 230,I = 1,N 
         X(I) = ZERO 
  230     CONTINUE 
      X(J) = ONE 
      KASE = 1 
      JUMP = 3 
      RETURN 
!                                                                       
!     ................ ENTRY   (JUMP = 3)                               
!     X HAS BEEN OVERWRITTEN BY A*X.                                    
!                                                                       
  300  CONTINUE 
      CALL DCOPY(N,X,1,V,1) 
      ESTOLD = EST 
      EST = DASUM(N,V,1) 
      DO 310,I = 1,N 
         IF ( NINT( SIGN(ONE,X(I)) ) .NE. ISGN(I) ) GOTO 320 
  310     CONTINUE 
!     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.     
      GOTO 410 
!                                                                       
  320  CONTINUE 
!     TEST FOR CYCLING.                                                 
      IF (EST .LE. ESTOLD) GOTO 410 
!                                                                       
      DO 330,I = 1,N 
         X(I) = SIGN(ONE,X(I)) 
         ISGN(I) = NINT(X(I)) 
  330     CONTINUE 
      KASE = 2 
      JUMP = 4 
      RETURN 
!                                                                       
!     ................ ENTRY   (JUMP = 4)                               
!     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.                         
!                                                                       
  400  CONTINUE 
      JLAST = J 
      J = IDAMAX(N,X,1) 
      IF (   (  X(JLAST) .NE. ABS(X(J))  ) .AND.                        &
     &       (ITER .LT. ITMAX)   ) THEN                                 
         ITER = ITER + 1 
         GOTO 220 
      ENDIF 
!                                                                       
!     ITERATION COMPLETE.  FINAL STAGE.                                 
!                                                                       
  410  CONTINUE 
      ALTSGN = ONE 
      DO 420,I = 1,N 
         X(I) = ALTSGN * (ONE + FLOAT(I-1)/FLOAT(N-1)) 
         ALTSGN = -ALTSGN 
  420     CONTINUE 
      KASE = 1 
      JUMP = 5 
      RETURN 
!                                                                       
!     ................ ENTRY   (JUMP = 5)                               
!     X HAS BEEN OVERWRITTEN BY A*X.                                    
!                                                                       
  500  CONTINUE 
      TEMP = TWO*DASUM(N,X,1)/FLOAT(3*N) 
      IF (TEMP .GT. EST) THEN 
         CALL DCOPY(N,X,1,V,1) 
         EST = TEMP 
      ENDIF 
!                                                                       
  510  KASE = 0 
      RETURN 
!                                                                       
      end subroutine                                        
                                                                        
! This is the code that has gone into netlib as a replacement.          
                                                                        
      REAL(DP) FUNCTION D1MACH(I) 
      INTEGER I 
!                                                                       
!  DOUBLE-PRECISION MACHINE CONSTANTS                                   
!                                                                       
!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.           
!                                                                       
!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.           
!                                                                       
!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.                 
!                                                                       
!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.                 
!                                                                       
!  D1MACH( 5) = LOG10(B)                                                
!                                                                       
!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,                 
!  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY            
!  REMOVING THE C FROM COLUMN 1.                                        
!  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.            
!  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)             
!                                                                       
!  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST     
!  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.  IF YOU DO NOT    
!  KNOW WHICH SET TO USE, TRY BOTH AND SEE WHICH GIVES PLAUSIBLE        
!  VALUES.                                                              
!                                                                       
!  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED     
!  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING     
!  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD          
!  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO    
!  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER      
!  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.                
!                                                                       
!  COMMENTS JUST BEFORE THE END STATEMENT (LINES STARTING WITH *)       
!  GIVE C SOURCE FOR D1MACH.                                            
!                                                                       
      INTEGER SMALL(2) 
      INTEGER LARGE(2) 
      INTEGER RIGHT(2) 
      INTEGER DIVER(2) 
      INTEGER LOG10(2) 
      INTEGER SC, CRAY1(38), J 
      COMMON /D9MACH/ CRAY1 
!/6S                                                                    
!/7S                                                                    
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC 
!/                                                                      
      REAL(DP) DMACH(5) 
!                                                                       
      EQUIVALENCE (DMACH(1),SMALL(1)) 
      EQUIVALENCE (DMACH(2),LARGE(1)) 
      EQUIVALENCE (DMACH(3),RIGHT(1)) 
      EQUIVALENCE (DMACH(4),DIVER(1)) 
      EQUIVALENCE (DMACH(5),LOG10(1)) 
!                                                                       
!     MACHINE CONSTANTS FOR BIG-ENDIAN IEEE ARITHMETIC (BINARY FORMAT)  
!     MACHINES IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST,      
!     SUCH AS THE AT&T 3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G.   
!     SUN 3), AND MACHINES THAT USE SPARC, HP, OR IBM RISC CHIPS.       
!                                                                       
!      DATA SMALL(1),SMALL(2) /    1048576,          0 /                
!      DATA LARGE(1),LARGE(2) / 2146435071,         -1 /                
!      DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /                
!      DATA DIVER(1),DIVER(2) / 1018167296,          0 /                
!      DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /, SC/987/       
!                                                                       
!     MACHINE CONSTANTS FOR LITTLE-ENDIAN (BINARY) IEEE ARITHMETIC      
!     MACHINES IN WHICH THE LEAST SIGNIFICANT BYTE IS STORED FIRST,     
!     E.G. IBM PCS AND OTHER MACHINES THAT USE INTEL 80X87 OR DEC       
!     ALPHA CHIPS.                                                      
!                                                                       
      DATA SMALL(1),SMALL(2) /          0,    1048576 / 
      DATA LARGE(1),LARGE(2) /         -1, 2146435071 / 
      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 / 
      DATA DIVER(1),DIVER(2) /          0, 1018167296 / 
!                                                                       
!     MACHINE CONSTANTS FOR AMDAHL MACHINES.                            
!                                                                       
!      DATA SMALL(1),SMALL(2) /    1048576,          0 /                
!      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /                
!      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /                
!      DATA DIVER(1),DIVER(2) /  873463808,          0 /                
!      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /, SC/987/       
!                                                                       
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.                  
!                                                                       
!      DATA SMALL(1) / ZC00800000 /                                     
!      DATA SMALL(2) / Z000000000 /                                     
!                                                                       
!      DATA LARGE(1) / ZDFFFFFFFF /                                     
!      DATA LARGE(2) / ZFFFFFFFFF /                                     
!                                                                       
!      DATA RIGHT(1) / ZCC5800000 /                                     
!      DATA RIGHT(2) / Z000000000 /                                     
!                                                                       
!      DATA DIVER(1) / ZCC6800000 /                                     
!      DATA DIVER(2) / Z000000000 /                                     
!                                                                       
!      DATA LOG10(1) / ZD00E730E7 /                                     
!      DATA LOG10(2) / ZC77800DC0 /, SC/987/                            
!                                                                       
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.                  
!                                                                       
!      DATA SMALL(1) / O1771000000000000 /                              
!      DATA SMALL(2) / O0000000000000000 /                              
!                                                                       
!      DATA LARGE(1) / O0777777777777777 /                              
!      DATA LARGE(2) / O0007777777777777 /                              
!                                                                       
!      DATA RIGHT(1) / O1461000000000000 /                              
!      DATA RIGHT(2) / O0000000000000000 /                              
!                                                                       
!      DATA DIVER(1) / O1451000000000000 /                              
!      DATA DIVER(2) / O0000000000000000 /                              
!                                                                       
!      DATA LOG10(1) / O1157163034761674 /                              
!      DATA LOG10(2) / O0006677466732724 /, SC/987/                     
!                                                                       
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.            
!                                                                       
!      DATA SMALL(1) / O1771000000000000 /                              
!      DATA SMALL(2) / O7770000000000000 /                              
!                                                                       
!      DATA LARGE(1) / O0777777777777777 /                              
!      DATA LARGE(2) / O7777777777777777 /                              
!                                                                       
!      DATA RIGHT(1) / O1461000000000000 /                              
!      DATA RIGHT(2) / O0000000000000000 /                              
!                                                                       
!      DATA DIVER(1) / O1451000000000000 /                              
!      DATA DIVER(2) / O0000000000000000 /                              
!                                                                       
!      DATA LOG10(1) / O1157163034761674 /                              
!      DATA LOG10(2) / O0006677466732724 /, SC/987/                     
!                                                                       
!     MACHINE CONSTANTS FOR FTN4 ON THE CDC 6000/7000 SERIES.           
!                                                                       
!      DATA SMALL(1) / 00564000000000000000B /                          
!      DATA SMALL(2) / 00000000000000000000B /                          
!                                                                       
!      DATA LARGE(1) / 37757777777777777777B /                          
!      DATA LARGE(2) / 37157777777777777774B /                          
!                                                                       
!      DATA RIGHT(1) / 15624000000000000000B /                          
!      DATA RIGHT(2) / 00000000000000000000B /                          
!                                                                       
!      DATA DIVER(1) / 15634000000000000000B /                          
!      DATA DIVER(2) / 00000000000000000000B /                          
!                                                                       
!      DATA LOG10(1) / 17164642023241175717B /                          
!      DATA LOG10(2) / 16367571421742254654B /, SC/987/                 
!                                                                       
!     MACHINE CONSTANTS FOR FTN5 ON THE CDC 6000/7000 SERIES.           
!                                                                       
!      DATA SMALL(1) / O"00564000000000000000" /                        
!      DATA SMALL(2) / O"00000000000000000000" /                        
!                                                                       
!      DATA LARGE(1) / O"37757777777777777777" /                        
!      DATA LARGE(2) / O"37157777777777777774" /                        
!                                                                       
!      DATA RIGHT(1) / O"15624000000000000000" /                        
!      DATA RIGHT(2) / O"00000000000000000000" /                        
!                                                                       
!      DATA DIVER(1) / O"15634000000000000000" /                        
!      DATA DIVER(2) / O"00000000000000000000" /                        
!                                                                       
!      DATA LOG10(1) / O"17164642023241175717" /                        
!      DATA LOG10(2) / O"16367571421742254654" /, SC/987/               
!                                                                       
!     MACHINE CONSTANTS FOR CONVEX C-1                                  
!                                                                       
!      DATA SMALL(1),SMALL(2) / '00100000'X, '00000000'X /              
!      DATA LARGE(1),LARGE(2) / '7FFFFFFF'X, 'FFFFFFFF'X /              
!      DATA RIGHT(1),RIGHT(2) / '3CC00000'X, '00000000'X /              
!      DATA DIVER(1),DIVER(2) / '3CD00000'X, '00000000'X /              
!      DATA LOG10(1),LOG10(2) / '3FF34413'X, '509F79FF'X /, SC/987/     
!                                                                       
!     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.                  
!                                                                       
!      DATA SMALL(1) / 201354000000000000000B /                         
!      DATA SMALL(2) / 000000000000000000000B /                         
!                                                                       
!      DATA LARGE(1) / 577767777777777777777B /                         
!      DATA LARGE(2) / 000007777777777777776B /                         
!                                                                       
!      DATA RIGHT(1) / 376434000000000000000B /                         
!      DATA RIGHT(2) / 000000000000000000000B /                         
!                                                                       
!      DATA DIVER(1) / 376444000000000000000B /                         
!      DATA DIVER(2) / 000000000000000000000B /                         
!                                                                       
!      DATA LOG10(1) / 377774642023241175717B /                         
!      DATA LOG10(2) / 000007571421742254654B /, SC/987/                
!                                                                       
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200              
!                                                                       
!     SMALL, LARGE, RIGHT, DIVER, LOG10 SHOULD BE DECLARED              
!     INTEGER SMALL(4), LARGE(4), RIGHT(4), DIVER(4), LOG10(4)          
!                                                                       
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -      
!     STATIC DMACH(5)                                                   
!                                                                       
!      DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/                      
!      DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/                         
!      DATA LOG10/40423K,42023K,50237K,74776K/, SC/987/                 
!                                                                       
!     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7              
!                                                                       
!      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /                  
!      DATA LARGE(1),LARGE(2) / '37777777, '37777577 /                  
!      DATA RIGHT(1),RIGHT(2) / '20000000, '00000333 /                  
!      DATA DIVER(1),DIVER(2) / '20000000, '00000334 /                  
!      DATA LOG10(1),LOG10(2) / '23210115, '10237777 /, SC/987/         
!                                                                       
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.              
!                                                                       
!      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /          
!      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /          
!      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /          
!      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /          
!      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/ 
!                                                                       
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,                     
!     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.                  
!                                                                       
!      DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /                  
!      DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /                  
!      DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /                  
!      DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /                  
!      DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /, SC/987/         
!                                                                       
!     MACHINE CONSTANTS FOR THE INTERDATA 8/32                          
!     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.                         
!                                                                       
!     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE                    
!     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.                        
!                                                                       
!      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /              
!      DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /              
!      DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /              
!      DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /              
!      DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /, SC/987/     
!                                                                       
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).                  
!                                                                       
!      DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /          
!      DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /          
!      DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /          
!      DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /          
!      DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /, SC/987/ 
!                                                                       
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).                  
!                                                                       
!      DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /          
!      DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /          
!      DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /          
!      DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /          
!      DATA LOG10(1),LOG10(2) / "177464202324, "047674776746 /, SC/987/ 
!                                                                       
!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING                  
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).                 
!                                                                       
!      DATA SMALL(1),SMALL(2) /    8388608,           0 /               
!      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /               
!      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /               
!      DATA DIVER(1),DIVER(2) /  620756992,           0 /               
!      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/      
!                                                                       
!      DATA SMALL(1),SMALL(2) / O00040000000, O00000000000 /            
!      DATA LARGE(1),LARGE(2) / O17777777777, O37777777777 /            
!      DATA RIGHT(1),RIGHT(2) / O04440000000, O00000000000 /            
!      DATA DIVER(1),DIVER(2) / O04500000000, O00000000000 /            
!      DATA LOG10(1),LOG10(2) / O07746420232, O20476747770 /, SC/987/   
!                                                                       
!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING                  
!     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).                 
      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /, SC/987/ 
!                                                                       
!     SMALL, LARGE, RIGHT, DIVER, LOG10 SHOULD BE DECLARED              
!     INTEGER SMALL(4), LARGE(4), RIGHT(4), DIVER(4), LOG10(4)          
!                                                                       
!      DATA SMALL(1),SMALL(2) /    128,      0 /                        
!      DATA SMALL(3),SMALL(4) /      0,      0 /                        
!                                                                       
!      DATA LARGE(1),LARGE(2) /  32767,     -1 /                        
!      DATA LARGE(3),LARGE(4) /     -1,     -1 /                        
!                                                                       
!      DATA RIGHT(1),RIGHT(2) /   9344,      0 /                        
!      DATA RIGHT(3),RIGHT(4) /      0,      0 /                        
!                                                                       
!      DATA DIVER(1),DIVER(2) /   9472,      0 /                        
!      DATA DIVER(3),DIVER(4) /      0,      0 /                        
!                                                                       
!      DATA LOG10(1),LOG10(2) /  16282,   8346 /                        
!      DATA LOG10(3),LOG10(4) / -31493, -12296 /, SC/987/               
!                                                                       
!      DATA SMALL(1),SMALL(2) / O000200, O000000 /                      
!      DATA SMALL(3),SMALL(4) / O000000, O000000 /                      
!                                                                       
!      DATA LARGE(1),LARGE(2) / O077777, O177777 /                      
!      DATA LARGE(3),LARGE(4) / O177777, O177777 /                      
!                                                                       
!      DATA RIGHT(1),RIGHT(2) / O022200, O000000 /                      
!      DATA RIGHT(3),RIGHT(4) / O000000, O000000 /                      
!                                                                       
!      DATA DIVER(1),DIVER(2) / O022400, O000000 /                      
!      DATA DIVER(3),DIVER(4) / O000000, O000000 /                      
!                                                                       
!      DATA LOG10(1),LOG10(2) / O037632, O020232 /                      
!      DATA LOG10(3),LOG10(4) / O102373, O147770 /, SC/987/             
!                                                                       
!     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS                 
!     WITH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,                   
!     SUPPLIED BY IGOR BRAY.                                            
!                                                                       
!      DATA SMALL(1),SMALL(2) / :10000000000, :00000100001 /            
!      DATA LARGE(1),LARGE(2) / :17777777777, :37777677775 /            
!      DATA RIGHT(1),RIGHT(2) / :10000000000, :00000000122 /            
!      DATA DIVER(1),DIVER(2) / :10000000000, :00000000123 /            
!      DATA LOG10(1),LOG10(2) / :11504046501, :07674600177 /, SC/987/   
!                                                                       
!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000                    
!                                                                       
!      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /                 
!      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /                 
!      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /                 
!      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /                 
!      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /, SC/987/        
!                                                                       
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.                     
!                                                                       
!      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /          
!      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /          
!      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /          
!      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /          
!      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/ 
!                                                                       
!     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER                   
!                                                                       
!      DATA SMALL(1),SMALL(2) /        128,           0 /               
!      DATA LARGE(1),LARGE(2) /     -32769,          -1 /               
!      DATA RIGHT(1),RIGHT(2) /       9344,           0 /               
!      DATA DIVER(1),DIVER(2) /       9472,           0 /               
!      DATA LOG10(1),LOG10(2) /  546979738,  -805796613 /, SC/987/      
!                                                                       
!     MACHINE CONSTANTS FOR THE VAX-11 WITH                             
!     FORTRAN IV-PLUS COMPILER                                          
!                                                                       
!      DATA SMALL(1),SMALL(2) / Z00000080, Z00000000 /                  
!      DATA LARGE(1),LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /                  
!      DATA RIGHT(1),RIGHT(2) / Z00002480, Z00000000 /                  
!      DATA DIVER(1),DIVER(2) / Z00002500, Z00000000 /                  
!      DATA LOG10(1),LOG10(2) / Z209A3F9A, ZCFF884FB /, SC/987/         
!                                                                       
!     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2                         
!                                                                       
!      DATA SMALL(1),SMALL(2) /       '80'X,        '0'X /              
!      DATA LARGE(1),LARGE(2) / 'FFFF7FFF'X, 'FFFFFFFF'X /              
!      DATA RIGHT(1),RIGHT(2) /     '2480'X,        '0'X /              
!      DATA DIVER(1),DIVER(2) /     '2500'X,        '0'X /              
!      DATA LOG10(1),LOG10(2) / '209A3F9A'X, 'CFF884FB'X /, SC/987/     
!                                                                       
!  ***  ISSUE STOP 779 IF ALL DATA STATEMENTS ARE COMMENTED...          
      IF (SC .NE. 987) THEN 
         DMACH(1) = 1.D13 
         IF (      SMALL(1) .EQ. 1117925532                             &
     &       .AND. SMALL(2) .EQ. -448790528) THEN                       
!           *** IEEE BIG ENDIAN ***                                     
            SMALL(1) = 1048576 
            SMALL(2) = 0 
            LARGE(1) = 2146435071 
            LARGE(2) = -1 
            RIGHT(1) = 1017118720 
            RIGHT(2) = 0 
            DIVER(1) = 1018167296 
            DIVER(2) = 0 
            LOG10(1) = 1070810131 
            LOG10(2) = 1352628735 
         ELSE IF ( SMALL(2) .EQ. 1117925532                             &
     &       .AND. SMALL(1) .EQ. -448790528) THEN                       
!           *** IEEE LITTLE ENDIAN ***                                  
            SMALL(2) = 1048576 
            SMALL(1) = 0 
            LARGE(2) = 2146435071 
            LARGE(1) = -1 
            RIGHT(2) = 1017118720 
            RIGHT(1) = 0 
            DIVER(2) = 1018167296 
            DIVER(1) = 0 
            LOG10(2) = 1070810131 
            LOG10(1) = 1352628735 
         ELSE IF ( SMALL(1) .EQ. -2065213935                            &
     &       .AND. SMALL(2) .EQ. 10752) THEN                            
!               *** VAX WITH D_FLOATING ***                             
            SMALL(1) = 128 
            SMALL(2) = 0 
            LARGE(1) = -32769 
            LARGE(2) = -1 
            RIGHT(1) = 9344 
            RIGHT(2) = 0 
            DIVER(1) = 9472 
            DIVER(2) = 0 
            LOG10(1) = 546979738 
            LOG10(2) = -805796613 
         ELSE IF ( SMALL(1) .EQ. 1267827943                             &
     &       .AND. SMALL(2) .EQ. 704643072) THEN                        
!               *** IBM MAINFRAME ***                                   
            SMALL(1) = 1048576 
            SMALL(2) = 0 
            LARGE(1) = 2147483647 
            LARGE(2) = -1 
            RIGHT(1) = 856686592 
            RIGHT(2) = 0 
            DIVER(1) = 873463808 
            DIVER(2) = 0 
            LOG10(1) = 1091781651 
            LOG10(2) = 1352628735 
         ELSE IF ( SMALL(1) .EQ. 1120022684                             &
     &       .AND. SMALL(2) .EQ. -448790528) THEN                       
!           *** CONVEX C-1 ***                                          
            SMALL(1) = 1048576 
            SMALL(2) = 0 
            LARGE(1) = 2147483647 
            LARGE(2) = -1 
            RIGHT(1) = 1019215872 
            RIGHT(2) = 0 
            DIVER(1) = 1020264448 
            DIVER(2) = 0 
            LOG10(1) = 1072907283 
            LOG10(2) = 1352628735 
         ELSE IF ( SMALL(1) .EQ. 815547074                              &
     &       .AND. SMALL(2) .EQ. 58688) THEN                            
!           *** VAX G-FLOATING ***                                      
            SMALL(1) = 16 
            SMALL(2) = 0 
            LARGE(1) = -32769 
            LARGE(2) = -1 
            RIGHT(1) = 15552 
            RIGHT(2) = 0 
            DIVER(1) = 15568 
            DIVER(2) = 0 
            LOG10(1) = 1142112243 
            LOG10(2) = 2046775455 
         ELSE 
            DMACH(2) = 1.D27 + 1 
            DMACH(3) = 1.D27 
            LARGE(2) = LARGE(2) - RIGHT(2) 
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN 
               CRAY1(1) = 67291416 
               DO 10 J = 1, 20 
   10             CRAY1(J+1) = CRAY1(J) + CRAY1(J) 
               CRAY1(22) = CRAY1(21) + 321322 
               DO 20 J = 22, 37 
   20             CRAY1(J+1) = CRAY1(J) + CRAY1(J) 
               IF (CRAY1(38) .EQ. SMALL(1)) THEN 
!                  *** CRAY ***                                         
!                 SMALL(1) = 2332160919536140288                        
                  SMALL(1) = 2332160 
                  SMALL(1) = 1000000*SMALL(1) + 919536 
                  SMALL(1) = 1000000*SMALL(1) + 140288 
                  SMALL(2) = 0 
!                 LARGE(1) = 6917247552664371199                        
                  LARGE(1) = 6917247 
                  LARGE(1) = 1000000*LARGE(1) + 552664 
                  LARGE(1) = 1000000*LARGE(1) + 371199 
!                 LARGE(2) = 281474976710654                            
                  LARGE(2) = 28147497 
                  LARGE(2) = 10000000*LARGE(2) + 6710654 
!                 RIGHT(1) = 4585649583081652224                        
                  RIGHT(1) = 4585649 
                  RIGHT(1) = 1000000*RIGHT(1) + 583081 
                  RIGHT(1) = 1000000*RIGHT(1) + 652224 
                  RIGHT(2) = 0 
!                 DIVER(1) = 4585931058058362880                        
                  DIVER(1) = 4585931 
                  DIVER(1) = 1000000*DIVER(1) + 058058 
                  DIVER(1) = 1000000*DIVER(1) + 362880 
                  DIVER(2) = 0 
!                 LOG10(1) = 4611574008272714703                        
                  LOG10(1) = 4611574 
                  LOG10(1) = 1000000*LOG10(1) +   8272 
                  LOG10(1) = 1000000*LOG10(1) + 714703 
!                 LOG10(2) = 272234615232940                            
                  LOG10(2) = 27223461 
                  LOG10(2) = 10000000*LOG10(2) + 5232940 
               ELSE 
                  WRITE(*,9000) 
                  STOP 779 
                  END IF 
            ELSE 
               WRITE(*,9000) 
               STOP 779 
               END IF 
            END IF 
         SC = 987 
         END IF 
!                                                                       
!  ***  ISSUE STOP 778 IF ALL DATA STATEMENTS ARE OBVIOUSLY WRONG...    
      IF (DMACH(4) .GE. 1.0D0) STOP 778 
!C/6S                                                                   
!C     IF (I .LT. 1  .OR.  I .GT. 5)                                    
!C    1   CALL SETERR(24HD1MACH - I OUT OF BOUNDS,24,1,2)               
!C/7S                                                                   
!      IF (I .LT. 1  .OR.  I .GT. 5)                                    
!     1   CALL SETERR('D1MACH - I OUT OF BOUNDS',24,1,2)                
!C/                                                                     
      IF (I .LT. 1 .OR. I .GT. 5) THEN 
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.' 
         STOP 
         END IF 
      D1MACH = DMACH(I) 
      RETURN 
!/6S                                                                    
!9000 FORMAT(/46H Adjust D1MACH by uncommenting data statements/        
!    *30H appropriate for your machine.)                                
!/7S                                                                    
 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/         &
     &' appropriate for your machine.')                                 
!/                                                                      
!                                                                       
! /* C source for D1MACH -- remove the * in column 1 */                 
!#include <stdio.h>                                                     
!#include <float.h>                                                     
!#include <math.h>                                                      
!                                                                       
!double d1mach_(long *i)                                                
!{                                                                      
!     switch(*i){                                                       
!       case 1: return DBL_MIN;                                         
!       case 2: return DBL_MAX;                                         
!       case 3: return DBL_EPSILON/FLT_RADIX;                           
!       case 4: return DBL_EPSILON;                                     
!       case 5: return log10(FLT_RADIX);                                
!       }                                                               
!                                                                       
!     fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);           
!     exit(1);                                                          
!     return 0; /* for compilers that complain of missing return values 
!     }                                                                 
      end function
                                                                        
                                                                        
      subroutine sprt(n, array) 
! sprt prints an array.                                                 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension array(n) 
                                                                        
      write(6,900) (array(i), i=1,n) 
  900 format(1h ,(7(1pe11.3))) 
      return 
      end subroutine                                        
                                                                        
      subroutine mprt(nrowd, nrow, ncol, array) 
! mprt prints a matrix.                                                 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension array(nrowd, ncol) 
                                                                        
      do 400 i = 1, nrow 
         write(6,900) i,(array(i,j), j=1,ncol) 
  400 continue 
  900 format(1h ,i5,(6(1pe11.3))) 
      return 
      end subroutine                                        
                                                                        
                                                                        
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
!     THE AUGUST 27 1992 VERSION OF COLROW IN WHICH X IS NO LONGER      
!     REQUIRED, WITH THE SOLUTION BEING RETURNED IN B, THE RIGHT        
!     HAND SIDE.  IN ADDITION, ALL VARIABLES ARE EXPLICITLY DECLARED.   
!     A PARAMETER "JOB" IS INCLUDED, TO SPECIFY WHICH OF A.X = B OR     
!     TRANSPOSE(A).X = B IS TO BE SOLVED.                               
                                                                        
      SUBROUTINE COLROW(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,            &
     &             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,IFLAG,JOB)       
!                                                                       
!***************************************************************        
!                                                                       
!  THIS PROGRAM SOLVES ONE OF THE LINEAR SYSTEMS  A*X = B OR            
!  TRANSPOSE(A)*X = B, WHERE  A IS AN ALMOST BLOCK DIAGONAL             
!  MATRIX OF THE FORM                                                   
!                                                                       
!               TOPBLK                                                  
!               ARRAY(1)                                                
!                     ARRAY(2)                                          
!                          .                                            
!                             .                                         
!                                .                                      
!                                   .                                   
!                                    ARRAY(NBLOKS)                      
!                                           BOTBLK                      
!                                                                       
!  WHERE                                                                
!           TOPBLK IS  NRWTOP  BY NOVRLP                                
!           ARRAY(K), K=1,NBLOKS, ARE NRWBLK BY NRWBLK+NOVRLP           
!           BOTBLK IS NRWBOT BY NOVRLP,                                 
!  AND                                                                  
!           NOVRLP = NRWTOP + NRWBOT                                    
!  WITH                                                                 
!           NOVRLP.LE.NRWBLK .                                          
!                                                                       
!  THE LINEAR SYSTEM IS OF ORDER  N = NBLOKS*NRWBLK + NOVRLP.           
!                                                                       
!  THE METHOD IMPLEMENTED IS BASED ON GAUSS ELIMINATION WITH            
!  ALTERNATE ROW AND COLUMN ELIMINATION WITH PARTIAL PIVOTING,          
!  WHICH PRODUCES A STABLE DECOMPOSITION OF THE MATRIX  A               
!  WITHOUT INTRODUCING FILL-IN.                                         
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!  TO OBTAIN A SINGLE PRECISION VERSION OF THIS PACKAGE, REMOVE         
!  ALL REAL(DP) STATEMENTS.  THERE IS ONE SUCH STATEMENT        
!  IN C O L R O W, THREE IN C R D C M P, AND TWO IN C R S O L V.        
!  IN ADDITION, REFERENCES TO BUILT-IN FUNCTIONS ABS AND MAX         
!  MUST BE REPLACED BY ABS AND AMAX1, RESPECTIVELY.  ABS OCCURS        
!  NINE TIMES, IN C R D C M P.  MAX OCCURS FOUR TIMES, IN             
!  C R D C M P.  FINALLY, ZERO IS INITIALISED TO 0.D0 IN A              
!  DATA STATEMENT IN C R D C M P.  THIS MUST BE REPLACED BY:            
!               DATA ZERO/0.0/                                          
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               *****  PARAMETERS  *****                                
!                                                                       
!       *** ON ENTRY ...                                                
!                                                                       
!               N      - INTEGER                                        
!                         THE ORDER OF THE LINEAR SYSTEM,               
!                         GIVEN BY NBLOKS*NRWBLK + NOVRLP               
!                                                                       
!               TOPBLK - REAL(DP)(NRWTOP,NOVRLP)                
!                         THE FIRST BLOCK OF THE ALMOST BLOCK           
!                         DIAGONAL MATRIX A                             
!                                                                       
!               NRWTOP - INTEGER                                        
!                         NUMBER OF ROWS IN THE BLOCK TOPBLK            
!                                                                       
!               NOVRLP - INTEGER                                        
!                         THE NUMBER OF COLUMNS IN WHICH SUCC-          
!                         ESSIVE BLOCKS OVERLAP, WHERE                  
!                                NOVRLP = NRWTOP + NRWBOT               
!                                                                       
!               ARRAY  - REAL(DP)(NRWBLK,NCLBLK,NBLOKS)         
!                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK           
!                         BY NCLBLK BLOCK OF THE MATRIX A               
!                                                                       
!               NRWBLK - INTEGER                                        
!                         NUMBER OF ROWS IN K-TH BLOCK                  
!                                                                       
!               NCLBLK - INTEGER                                        
!                         NUMBER OF COLUMNS IN K-TH BLOCK               
!                                                                       
!               NBLOKS - INTEGER                                        
!                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN          
!                         THE MATRIX A                                  
!                                                                       
!               BOTBLK - REAL(DP)(NRWBOT,NOVRLP)                
!                         THE LAST BLOCK OF THE MATRIX A                
!                                                                       
!               NRWBOT - INTEGER                                        
!                         NUMBER OF ROWS IN THE BLOCK BOTBLK            
!                                                                       
!                PIVOT - INTEGER(N)                                     
!                         WORK SPACE                                    
!                                                                       
!                    B - REAL(DP)(N)                            
!                         THE RIGHT HAND SIDE VECTOR                    
!                                                                       
!               JOB    - INTEGER, INDICATING:                           
!                      = 0: SOLVE A*X = B;                              
!                      NON-ZERO: SOLVE TRANSPOSE(A)*X = B.              
!                                                                       
!       *** ON RETURN  ...                                              
!                                                                       
!               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE             
!                        DESIRED DECOMPOSITION OF THE MATRIX A          
!                        (IF IFLAG = 0)                                 
!                                                                       
!                PIVOT - INTEGER(N)                                     
!                         RECORDS THE PIVOTING INDICES DETER-           
!                         MINED IN THE DECOMPOSITION                    
!                                                                       
!                    B - REAL(DP)(N)                            
!                         THE SOLUTION VECTOR (IF IFLAG = 0)            
!                                                                       
!               IFLAG  - INTEGER                                        
!                         =  1, IF INPUT PARAMETERS ARE INVALID         
!                         = -1, IF MATRIX IS SINGULAR                   
!                         =  0, OTHERWISE                               
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               *****  AUXILIARY PROGRAMS  *****                        
!                                                                       
!       CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,       
!    *     BOTBLK,NRWBOT,PIVOT,IFLAG)                                   
!            - DECOMPOSES THE MATRIX  A  USING MODIFIED                 
!              ALTERNATE ROW AND COLUMN ELIMINATON WITH                 
!              PARTIAL PIVOTING, AND IS USED FOR THIS                   
!              PURPOSE IN C O L R O W.                                  
!              THE ARGUMENTS ARE AS IN C O L R O W.                     
!                                                                       
!       CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,         
!    *     BOTBLK,NRWBOT,PIVOT,B,JOB)                                   
!            - SOLVES THE SYSTEM A*X = B ONCE A IS DECOMPOSED.          
!              THE ARGUMENTS ARE ALL AS IN C O L R O W.                 
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!       THE SUBROUTINE  C O L R O W  AUTOMATICALLY SOLVES THE           
!  INPUT SYSTEM WHEN IFLAG=0.  C O L R O W  IS CALLED ONLY ONCE         
!  FOR A GIVEN SYSTEM. THE SOLUTION FOR A SEQUENCE OF P RIGHT           
!  HAND SIDES CAN BE OBTAINED BY ONE CALL TO  C O L R O W  AND          
!  P-1 CALLS TO CRSLVE ONLY. SINCE THE ARRAYS TOPBLK,ARRAY,             
!  BOTBLK AND PIVOT CONTAIN THE DECOMPOSITION OF THE GIVEN              
!  COEFFICIENT MATRIX AND PIVOTING INFORMATION ON RETURN FROM           
!  C O L R O W , THEY MUST NOT BE ALTERED BETWEEN SUCCESSIVE            
!  CALLS TO CRSLVE WITH THE SAME LEFT HAND SIDES. FOR THE               
!  SAME REASON, IF THE USER WISHES TO SAVE THE COEFFICIENT              
!  MATRIX, THE ARRAYS TOPBLK,ARRAY,BOTBLK MUST BE COPIED                
!  BEFORE A CALL TO  C O L R O W .                                      
!                                                                       
!***********************************************************************
        IMPLICIT NONE 
        REAL(DP) TOPBLK,ARRAY,BOTBLK,B 
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*),   &
     &          IFLAG,JOB, IDAMAX,i,IFAIL                               
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),              &
     &          BOTBLK(NRWBOT,*),B(*)                                   
                                                                        
!        NOUT = 6                                                       
                                                                        
!       DO THE FACTORIZATION USING CRDCMP:                              
!                                                                       
        CALL CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,  &
     &          BOTBLK,NRWBOT,PIVOT,IFLAG)                              
        IF(IFLAG.NE.0) RETURN 
!                                                                       
!     *****************solving the linear system********************    
        job=0 
                                                                        
        CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,    &
     &          BOTBLK,NRWBOT,PIVOT,B,JOB)                              
        RETURN 
      end subroutine                                        
                                                                        
                                                                        
      SUBROUTINE INVERSE(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,    &
     &          NBLOKS,BOTBLK,NRWBOT,PIVOT,INMAT)                       
!******************************************************************     
!     INVERSE COMPUTES THE INVERSE OF A MATRIX                          
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
        IMPLICIT NONE 
        REAL(DP) TOPBLK,ARRAY,BOTBLK,WORK,                      &
     &          INMAT                                                   
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*) 
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),              &
     &          BOTBLK(NRWBOT,*),WORK(N),                               &
     &          INMAT(N,N)                                              
        INTEGER K,L,J 
                                                                        
        do 300 k=1,N 
           do 330 l=1,N 
              WORK(l)=0.0d0 
              if (l.eq.k) WORK(l)=1.0d0 
  330      continue 
                                                                        
          CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,                &
     &                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,               &
     &                 WORK,0)                                          
                                                                        
          do 440 l=1,N 
             INMAT(l,k)=WORK(l) 
  440        continue 
  300     continue 
                                                                        
          RETURN 
      end subroutine                                        
                                                                        
!******************************************************************     
                                                                        
                                                                        
                                                                        
      SUBROUTINE CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,     &
     &   NBLOKS,BOTBLK,NRWBOT,PIVOT,IFLAG)                              
!                                                                       
!***************************************************************        
!                                                                       
!  C R D C M P DECOMPOSES THE ALMOST BLOCK DIAGONAL MATRIX A            
!  USING MODIFIED ALTERNATE ROW AND COLUMN ELIMINATION WITH             
!  PARTIAL PIVOTING.  THE MATRIX  A  IS STORED IN THE ARRAYS            
!  TOPBLK, ARRAY, AND BOTBLK.                                           
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               *****  PARAMETERS  *****                                
!                                                                       
!       *** ON ENTRY ...                                                
!                                                                       
!               N      - INTEGER                                        
!                         THE ORDER OF THE LINEAR SYSTEM,               
!                         GIVEN BY NBLOKS*NRWBLK + NOVRLP               
!                                                                       
!               TOPBLK - REAL(DP)(NRWTOP,NOVRLP)                
!                         THE FIRST BLOCK OF THE ALMOST BLOCK           
!                         DIAGONAL MATRIX A TO BE DECOMPOSED            
!                                                                       
!               NRWTOP - INTEGER                                        
!                         NUMBER OF ROWS IN THE BLOCK TOPBLK            
!                                                                       
!               NOVRLP - INTEGER                                        
!                         THE NUMBER OF COLUMNS IN WHICH SUCC-          
!                         ESSIVE BLOCKS OVERLAP, WHERE                  
!                                NOVRLP = NRWTOP + NRWBOT               
!                                                                       
!               ARRAY  - REAL(DP)(NRWBLK,NCLBLK,NBLOKS)         
!                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK           
!                         BY NCLBLK BLOCK OF THE MATRIX A               
!                                                                       
!               NRWBLK - INTEGER                                        
!                         NUMBER OF ROWS IN K-TH BLOCK                  
!                                                                       
!               NCLBLK - INTEGER                                        
!                         NUMBER OF COLUMNS IN K-TH BLOCK               
!                                                                       
!               NBLOKS - INTEGER                                        
!                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN          
!                         THE MATRIX A                                  
!                                                                       
!               BOTBLK - REAL(DP)(NRWBOT,NOVRLP)                
!                         THE LAST BLOCK OF THE MATRIX A                
!                                                                       
!               NRWBOT - INTEGER                                        
!                         NUMBER OF ROWS IN THE BLOCK BOTBLK            
!                                                                       
!                PIVOT - INTEGER(N)                                     
!                         WORK SPACE                                    
!                                                                       
!       *** ON RETURN  ...                                              
!                                                                       
!               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE             
!                        DESIRED DECOMPOSITION OF THE MATRIX A          
!                        (IF IFLAG = 0)                                 
!                                                                       
!                PIVOT - INTEGER(N)                                     
!                         RECORDS THE PIVOTING INDICES DETER-           
!                         MINED IN THE DECOMPOSITION                    
!                                                                       
!               IFLAG  - INTEGER                                        
!                         =  1, IF INPUT PARAMETERS ARE INVALID         
!                         = -1, IF MATRIX IS SINGULAR                   
!                         =  0, OTHERWISE                               
!                                                                       
!***************************************************************        
!                                                                       
      implicit integer (i-n)
      IMPLICIT REAL(DP)(A-H,O-Z) 
      INTEGER PIVOT(*) 
      DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),BOTBLK(NRWBOT,*) 
      DATA ZERO / 0.0_dp / 
!                                                                       
!***************************************************************        
!                                                                       
!          ****  DEFINE THE CONSDTANTS USED THROUGHOUT  ****            
!                                                                       
!***************************************************************        
!                                                                       
      IFLAG = 0 
      PIVMAX = ZERO 
      NRWTP1 = NRWTOP+1 
      NROWEL = NRWBLK-NRWTOP 
      NRWEL1 = NROWEL+1 
      NVRLP0 = NOVRLP-1 
!                                                                       
!***************************************************************        
!                                                                       
!          ****  CHECK VALIDITY OF THE INPUT PARAMETERS....             
!                                                                       
!               IF PARAMETERS ARE INVALID THEN TERMINATE AT 10;         
!                                         ELSE CONTINUE AT 100.         
!                                                                       
!***************************************************************        
!                                                                       
      IF (N.NE.NBLOKS*NRWBLK+NOVRLP) GO TO 10 
      IF (NOVRLP.NE.NRWTOP+NRWBOT) GO TO 10 
      IF (NCLBLK.NE.NOVRLP+NRWBLK) GO TO 10 
      IF (NOVRLP.GT.NRWBLK) GO TO 10 
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!          PARAMETERS ARE ACCEPTABLE - CONTINUE AT 100.                 
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
      GO TO 20 
   10 CONTINUE 
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!          PARAMETERS ARE INVALID.  SET IFLAG = 1, AND TERMINATE        
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
      IFLAG = 1 
      RETURN 
   20 CONTINUE 
!                                                                       
!***************************************************************        
!                                                                       
!               ****  FIRST, IN TOPBLK....                              
!                                                                       
!***************************************************************        
!                                                                       
!          ***  APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN            
!                 PIVOTING ....                                         
!                                                                       
!***************************************************************        
!                                                                       
      DO 110 I = 1, NRWTOP 
         IPLUS1 = I+1 
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               DETERMINE COLUMN PIVOT AND PIVOT INDEX                  
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
         IPVT = I 
         COLMAX = ABS(TOPBLK(I,I)) 
         DO 30 J = IPLUS1, NOVRLP 
            TEMPIV = ABS(TOPBLK(I,J)) 
            IF (TEMPIV.LE.COLMAX) GO TO 30 
            IPVT = J 
            COLMAX = TEMPIV 
   30    CONTINUE 
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               TEST FOR SINGULARITY:                                   
!                                                                       
!                       IF SINGULAR THEN TERMINATE AT 1000;             
!                                   ELSE CONTINUE.                      
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
         IF (PIVMAX+COLMAX.EQ.PIVMAX) GO TO 410 
         PIVMAX = MAX(COLMAX,PIVMAX) 
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               IF NECESSARY INTERCHANGE COLUMNS                        
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
         PIVOT(I) = IPVT 
         IF (IPVT.EQ.I) GO TO 60 
         DO 40 L = I, NRWTOP 
            SWAP = TOPBLK(L,IPVT) 
            TOPBLK(L,IPVT) = TOPBLK(L,I) 
            TOPBLK(L,I) = SWAP 
   40    CONTINUE 
         DO 50 L = 1, NRWBLK 
            SWAP = ARRAY(L,IPVT,1) 
            ARRAY(L,IPVT,1) = ARRAY(L,I,1) 
            ARRAY(L,I,1) = SWAP 
   50    CONTINUE 
   60    CONTINUE 
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               COMPUTE MULTIPLIERS AND PERFORM COLUMN                  
!                       ELIMINATION                                     
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
         COLPIV = TOPBLK(I,I) 
         DO 100 J = IPLUS1, NOVRLP 
            COLMLT = TOPBLK(I,J)/COLPIV 
            TOPBLK(I,J) = COLMLT 
            IF (IPLUS1.GT.NRWTOP) GO TO 80 
            DO 70 L = IPLUS1, NRWTOP 
               TOPBLK(L,J) = TOPBLK(L,J)-COLMLT*TOPBLK(L,I) 
   70       CONTINUE 
   80       CONTINUE 
            DO 90 L = 1, NRWBLK 
               ARRAY(L,J,1) = ARRAY(L,J,1)-COLMLT*ARRAY(L,I,1) 
   90       CONTINUE 
  100    CONTINUE 
  110 END DO 
!                                                                       
!***************************************************************        
!                                                                       
!          ****  IN EACH BLOCK ARRAY(,,K)....                           
!                                                                       
!***************************************************************        
!                                                                       
      INCR = 0 
      DO 320 K = 1, NBLOKS 
         KPLUS1 = K+1 
!                                                                       
!          *****************************************************        
!                                                                       
!          ***  FIRST APPLY NRWBLK-NRWTOP ROW ELIMINATIONS WITH         
!                       ROW PIVOTING....                                
!                                                                       
!          *****************************************************        
!                                                                       
         DO 180 J = NRWTP1, NRWBLK 
            JPLUS1 = J+1 
            JMINN = J-NRWTOP 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               DETERMINE ROW PIVOT AND PIVOT INDEX                     
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
            IPVT = JMINN 
            ROWMAX = ABS(ARRAY(JMINN,J,K)) 
            LOOP = JMINN+1 
            DO 120 I = LOOP, NRWBLK 
               TEMPIV = ABS(ARRAY(I,J,K)) 
               IF (TEMPIV.LE.ROWMAX) GO TO 120 
               IPVT = I 
               ROWMAX = TEMPIV 
  120       CONTINUE 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               TEST FOR SINGULARITY:                                   
!                                                                       
!                       IF SINGULAR THEN TERMINATE AT 1000;             
!                                   ELSE CONTINUE.                      
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
            IF (PIVMAX+ROWMAX.EQ.PIVMAX) GO TO 410 
            PIVMAX = MAX(ROWMAX,PIVMAX) 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               IF NECESSARY INTERCHANGE ROWS                           
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
            INCRJ = INCR+J 
            PIVOT(INCRJ) = INCR+IPVT+NRWTOP 
            IF (IPVT.EQ.JMINN) GO TO 140 
            DO 130 L = J, NCLBLK 
               SWAP = ARRAY(IPVT,L,K) 
               ARRAY(IPVT,L,K) = ARRAY(JMINN,L,K) 
               ARRAY(JMINN,L,K) = SWAP 
  130       CONTINUE 
  140       CONTINUE 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               COMPUTE MULTIPLERS                                      
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
            ROWPIV = ARRAY(JMINN,J,K) 
            DO 150 I = LOOP, NRWBLK 
               ARRAY(I,J,K) = ARRAY(I,J,K)/ROWPIV 
  150       CONTINUE 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               PERFORM ROW ELIMINATION WITH COLUMN INDEXING            
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
            DO 170 L = JPLUS1, NCLBLK 
               ROWMLT = ARRAY(JMINN,L,K) 
               DO 160 I = LOOP, NRWBLK 
                  ARRAY(I,L,K) = ARRAY(I,L,K)-ROWMLT*ARRAY(I,J,K) 
  160          CONTINUE 
  170       CONTINUE 
  180    CONTINUE 
!                                                                       
!          *****************************************************        
!                                                                       
!          ***  NOW APPLY NRWTOP COLUMN ELIMINATIONS WITH               
!                      COLUMN PIVOTING....                              
!                                                                       
!          *****************************************************        
!                                                                       
         DO 310 I = NRWEL1, NRWBLK 
            IPLUSN = I+NRWTOP 
            IPLUS1 = I+1 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               DETERMINE COLUMN PIVOT AND PIVOT INDEX                  
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
            IPVT = IPLUSN 
            COLMAX = ABS(ARRAY(I,IPVT,K)) 
            LOOP = IPLUSN+1 
            DO 190 J = LOOP, NCLBLK 
               TEMPIV = ABS(ARRAY(I,J,K)) 
               IF (TEMPIV.LE.COLMAX) GO TO 190 
               IPVT = J 
               COLMAX = TEMPIV 
  190       CONTINUE 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               TEST FOR SINGULARITY:                                   
!                                                                       
!                       IF SINGULAR THEN TERMINATE AT 1000;             
!                                   ELSE CONTINUE.                      
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
            IF (PIVMAX+COLMAX.EQ.PIVMAX) GO TO 410 
            PIVMAX = MAX(COLMAX,PIVMAX) 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               IF NECESSARY INTERCHANGE COLUMNS                        
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
            INCRN = INCR+IPLUSN 
            PIVOT(INCRN) = INCR+IPVT 
            IRWBLK = IPLUSN-NRWBLK 
            IF (IPVT.EQ.IPLUSN) GO TO 240 
            DO 200 L = I, NRWBLK 
               SWAP = ARRAY(L,IPVT,K) 
               ARRAY(L,IPVT,K) = ARRAY(L,IPLUSN,K) 
               ARRAY(L,IPLUSN,K) = SWAP 
  200       CONTINUE 
            IPVBLK = IPVT-NRWBLK 
            IF (K.EQ.NBLOKS) GO TO 220 
            DO 210 L = 1, NRWBLK 
               SWAP = ARRAY(L,IPVBLK,KPLUS1) 
               ARRAY(L,IPVBLK,KPLUS1) = ARRAY(L,IRWBLK,KPLUS1) 
               ARRAY(L,IRWBLK,KPLUS1) = SWAP 
  210       CONTINUE 
            GO TO 240 
  220       CONTINUE 
            DO 230 L = 1, NRWBOT 
               SWAP = BOTBLK(L,IPVBLK) 
               BOTBLK(L,IPVBLK) = BOTBLK(L,IRWBLK) 
               BOTBLK(L,IRWBLK) = SWAP 
  230       CONTINUE 
  240       CONTINUE 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               COMPUTE MULTIPLIERS AND PERFORM COLUMN                  
!                       ELIMINATION                                     
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
            COLPIV = ARRAY(I,IPLUSN,K) 
            DO 300 J = LOOP, NCLBLK 
               COLMLT = ARRAY(I,J,K)/COLPIV 
               ARRAY(I,J,K) = COLMLT 
               IF (I.EQ.NRWBLK) GO TO 260 
               DO 250 L = IPLUS1, NRWBLK 
                  ARRAY(L,J,K) = ARRAY(L,J,K)-COLMLT*ARRAY(L,IPLUSN,K) 
  250          CONTINUE 
  260          CONTINUE 
               JRWBLK = J-NRWBLK 
               IF (K.EQ.NBLOKS) GO TO 280 
               DO 270 L = 1, NRWBLK 
                  ARRAY(L,JRWBLK,KPLUS1) = ARRAY(L,JRWBLK,KPLUS1)-COLMLT&
     &               *ARRAY(L,IRWBLK,KPLUS1)                            
  270          CONTINUE 
               GO TO 300 
  280          CONTINUE 
               DO 290 L = 1, NRWBOT 
                  BOTBLK(L,JRWBLK) = BOTBLK(L,JRWBLK)-COLMLT*BOTBLK(L,  &
     &               IRWBLK)                                            
  290          CONTINUE 
  300       CONTINUE 
  310    CONTINUE 
         INCR = INCR+NRWBLK 
  320 END DO 
!                                                                       
!***************************************************************        
!                                                                       
!          ****  FINALLY, IN BOTBLK....                                 
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!          ***  APPLY NRWBOT ROW ELIMINATIONS WITH ROW                  
!                  PIVOTING....                                         
!                                                                       
!               IF BOT HAS JUST ONE ROW GO TO 500                       
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
      IF (NRWBOT.EQ.1) GO TO 400 
      DO 390 J = NRWTP1, NVRLP0 
         JPLUS1 = J+1 
         JMINN = J-NRWTOP 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               DETERMINE ROW PIVOT AND PIVOT INDEX                     
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
         IPVT = JMINN 
         ROWMAX = ABS(BOTBLK(JMINN,J)) 
         LOOP = JMINN+1 
         DO 330 I = LOOP, NRWBOT 
            TEMPIV = ABS(BOTBLK(I,J)) 
            IF (TEMPIV.LE.ROWMAX) GO TO 330 
            IPVT = I 
            ROWMAX = TEMPIV 
  330    CONTINUE 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               TEST FOR SINGULARITY:                                   
!                                                                       
!                       IF SINGULAR THEN TERMINATE AT 1000;             
!                                   ELSE CONTINUE.                      
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
         IF (PIVMAX+ROWMAX.EQ.PIVMAX) GO TO 410 
         PIVMAX = MAX(ROWMAX,PIVMAX) 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               IF NECESSARY INTERCHANGE ROWS                           
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
         INCRJ = INCR+J 
         PIVOT(INCRJ) = INCR+IPVT+NRWTOP 
         IF (IPVT.EQ.JMINN) GO TO 350 
         DO 340 L = J, NOVRLP 
            SWAP = BOTBLK(IPVT,L) 
            BOTBLK(IPVT,L) = BOTBLK(JMINN,L) 
            BOTBLK(JMINN,L) = SWAP 
  340    CONTINUE 
  350    CONTINUE 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               COMPUTE MULTIPLIERS                                     
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
         ROWPIV = BOTBLK(JMINN,J) 
         DO 360 I = LOOP, NRWBOT 
            BOTBLK(I,J) = BOTBLK(I,J)/ROWPIV 
  360    CONTINUE 
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               PERFORM ROW ELIMINATION WITH COLUMN INDEXING            
!                                                                       
!             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
         DO 380 L = JPLUS1, NOVRLP 
            ROWMLT = BOTBLK(JMINN,L) 
            DO 370 I = LOOP, NRWBOT 
               BOTBLK(I,L) = BOTBLK(I,L)-ROWMLT*BOTBLK(I,J) 
  370       CONTINUE 
  380    CONTINUE 
  390 END DO 
  400 CONTINUE 
!                                                                       
!***************************************************************        
!                                                                       
!          DONE PROVIDED THE LAST ELEMENT IS NOT ZERO                   
!                                                                       
!***************************************************************        
!                                                                       
      IF (PIVMAX+ABS(BOTBLK(NRWBOT,NOVRLP)).NE.PIVMAX) RETURN 
!                                                                       
!***************************************************************        
!                                                                       
!       ****  MATRIX IS SINGULAR - SET IFLAG = - 1.                     
!                                  TERMINATE AT 1000.                   
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
  410 CONTINUE 
      IFLAG = -1 
      RETURN 
      end subroutine                                        
                                                                        
        SUBROUTINE CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,            &
     &             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,JOB)             
!                                                                       
!***************************************************************        
!                                                                       
!  C R S L V E  SOLVES THE LINEAR SYSTEM                                
!                       A*X = B                                         
!  USING THE DECOMPOSITION ALREADY GENERATED IN  C R D C M P.           
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               *****  PARAMETERS  *****                                
!                                                                       
!       *** ON ENTRY  ...                                               
!                                                                       
!               TOPBLK - REAL(DP)(NRWTOP,NOVRLP)                
!                         OUTPUT FROM  C R D C M P                      
!                                                                       
!               NOVRLP - INTEGER                                        
!                         THE NUMBER OF COLUMNS IN WHICH SUCC-          
!                         ESSIVE BLOCKS OVERLAP, WHERE                  
!                                NOVRLP = NRWTOP + NRWBOT               
!                                                                       
!               NRWTOP - INTEGER                                        
!                         NUMBER OF ROWS IN THE BLOCK TOPBLK            
!                                                                       
!               ARRAY  - REAL(DP)(NRWBLK,NCLBLK,NBLOKS)         
!                         OUTPUT FROM  C R D C M P                      
!                                                                       
!               NRWBLK - INTEGER                                        
!                         NUMBER OF ROWS IN K-TH BLOCK                  
!                                                                       
!               NCLBLK - INTEGER                                        
!                         NUMBER OF COLUMNS IN K-TH BLOCK               
!                                                                       
!               NBLOKS - INTEGER                                        
!                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN          
!                         THE MATRIX A                                  
!                                                                       
!               BOTBLK - REAL(DP)(NRWBOT,NOVRLP)                
!                         OUTPUT FROM  C R D C M P                      
!                                                                       
!               NRWBOT - INTEGER                                        
!                         NUMBER OF ROWS IN THE BLOCK BOTBLK            
!                                                                       
!                PIVOT - INTEGER(N)                                     
!                         THE PIVOT VECTOR FROM  C R D C M P            
!                                                                       
!                    B - REAL(DP)(N)                            
!                         THE RIGHT HAND SIDE VECTOR                    
!                                                                       
!               JOB    - INTEGER, INDICATING:                           
!                      = 0: SOLVE A*X = B;                              
!                      NON-ZERO: SOLVE TRANSPOSE(A)*X = B.              
!                                                                       
!       *** ON RETURN  ...                                              
!                                                                       
!                    B - REAL(DP)(N)                            
!                         THE SOLUTION VECTOR                           
!                                                                       
!***************************************************************        
!                                                                       
        IMPLICIT NONE 
        REAL(DP) TOPBLK,ARRAY,BOTBLK,B 
        REAL(DP) DOTPRD,BJ,XINCRJ,BINCRJ,SWAP,BI 
        INTEGER NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*),     &
     &          JOB                                                     
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),              &
     &          BOTBLK(NRWBOT,*),B(*)                                   
        INTEGER NRWTP1,NRWBK1,NVRLP1,NRWBT1,NROWEL,NVRLP0,NBLKS1,       &
     &          NBKTOP,J,I,LOOP,INCR,INCRJ,INCRI,JPIVOT,JRWTOP,         &
     &          LL,L1,IPLUSN,INCRN,NRWTP0,NRWEL1,K,INCRTP,NRWBTL,       &
     &          IPVTN,NRWELL,IPVTI,L                                    
!                                                                       
!***************************************************************        
!                                                                       
!          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****             
!                                                                       
!***************************************************************        
!                                                                       
        NRWTP1 = NRWTOP+1 
        NRWBK1 = NRWBLK+1 
        NVRLP1 = NOVRLP+1 
        NRWTP0 = NRWTOP-1 
        NRWBT1 = NRWBOT+1 
        NROWEL = NRWBLK-NRWTOP 
        NRWEL1 = NROWEL+1 
        NVRLP0 = NOVRLP-1 
        NBLKS1 = NBLOKS+1 
        NBKTOP = NRWBLK+NRWTOP 
!                                                                       
!       IF JOB IS NON-ZERO, TRANSFER TO THE SECTION DEALING WITH        
!       TRANSPOSE(A)*X = B.                                             
!                                                                       
        IF ( JOB .NE. 0 ) GO TO 530 
!                                                                       
!***************************************************************        
!                                                                       
!               ****  FORWARD RECURSION  ****                           
!                                                                       
!***************************************************************        
!                                                                       
!          ***  FIRST, IN TOPBLK....                                    
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               FORWARD SOLUTION                                        
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
        DO 130 J = 1,NRWTOP 
           B(J) = B(J)/TOPBLK(J,J) 
           IF(J.EQ.NRWTOP)GO TO 120 
              BJ = -B(J) 
              LOOP = J+1 
              DO 110 I = LOOP,NRWTOP 
                 B(I) = B(I)+TOPBLK(I,J)*BJ 
  110                       CONTINUE 
  120                               CONTINUE 
  130                                    CONTINUE 
!                                                                       
!       ********************************************************        
!                                                                       
!          ***  IN EACH BLOCK ARRAY(,,K)....                            
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
        INCR = 0 
        DO 280 K = 1,NBLOKS 
           INCRTP = INCR+NRWTOP 
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               FORWARD MODIFICATION                                    
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
           DO 220 J = 1,NRWTOP 
              INCRJ = INCR+J 
              XINCRJ = -B(INCRJ) 
              DO 210 I = 1,NRWBLK 
                 INCRI = INCRTP+I 
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ 
  210                       CONTINUE 
  220                               CONTINUE 
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               FORWARD ELIMINATION                                     
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
           DO 240 J = NRWTP1,NRWBLK 
              INCRJ = INCR+J 
              JPIVOT = PIVOT(INCRJ) 
              IF(JPIVOT.EQ.INCRJ)GO TO 225 
                 SWAP = B(INCRJ) 
                 B(INCRJ) = B(JPIVOT) 
                 B(JPIVOT) = SWAP 
  225                       CONTINUE 
              BINCRJ = -B(INCRJ) 
              LOOP = J-NRWTP0 
              DO 230 I = LOOP,NRWBLK 
                 INCRI = INCRTP+I 
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BINCRJ 
  230                       CONTINUE 
  240                               CONTINUE 
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               FORWARD SOLUTION                                        
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
           DO 270 J = NRWBK1,NBKTOP 
              INCRJ = INCR+J 
              JRWTOP = J -NRWTOP 
              B(INCRJ) = B(INCRJ)/ARRAY(JRWTOP,J,K) 
              IF(J.EQ.NBKTOP)GO TO 260 
                 XINCRJ = -B(INCRJ) 
                 LOOP = J-NRWTP0 
                 DO 250 I = LOOP,NRWBLK 
                    INCRI = INCRTP+I 
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ 
  250                             CONTINUE 
  260                                        CONTINUE 
  270                                                CONTINUE 
           INCR = INCR+NRWBLK 
  280           CONTINUE 
!                                                                       
!       ********************************************************        
!                                                                       
!          ***  FINALLY, IN BOTBLK....                                  
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               FORWARD MODIFICATION                                    
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
        INCRTP = INCR+NRWTOP 
        DO 320 J = 1,NRWTOP 
           INCRJ = INCR+J 
           XINCRJ = -B(INCRJ) 
           DO 310 I = 1,NRWBOT 
              INCRI = INCRTP+I 
              B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ 
  310                 CONTINUE 
  320                      CONTINUE 
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               FORWARD ELIMINATION                                     
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
        IF(NRWBOT.EQ.1)GO TO 350 
           DO 340 J = NRWTP1,NVRLP0 
              INCRJ = INCR+J 
              JPIVOT = PIVOT(INCRJ) 
              IF(JPIVOT.EQ.INCRJ)GO TO 325 
                 SWAP = B(INCRJ) 
                 B(INCRJ) = B(JPIVOT) 
                 B(JPIVOT) = SWAP 
  325                       CONTINUE 
              BINCRJ = -B(INCRJ) 
              LOOP = J-NRWTP0 
              DO 330 I = LOOP,NRWBOT 
                 INCRI = INCRTP+I 
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*BINCRJ 
  330                       CONTINUE 
  340                               CONTINUE 
  350                                    CONTINUE 
!                                                                       
!***************************************************************        
!                                                                       
!               ****  BACKWARD RECURSION  ****                          
!                                                                       
!***************************************************************        
!                                                                       
!          ***  FIRST IN BOTBLK....                                     
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               BACKWARD SOLUTION                                       
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
        DO 430 LL = 1,NRWBOT 
           J = NVRLP1-LL 
           INCRJ = INCR+J 
           NRWBTL = NRWBT1-LL 
           B(INCRJ) = B(INCRJ)/BOTBLK(NRWBTL,J) 
           IF(LL.EQ.NRWBOT)GO TO 420 
              XINCRJ = -B(INCRJ) 
              LOOP = NRWBOT-LL 
              DO 410 I = 1,LOOP 
                 INCRI = INCRTP+I 
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ 
  410                       CONTINUE 
  420                               CONTINUE 
  430                                    CONTINUE 
!                                                                       
!       ********************************************************        
!                                                                       
!          ***  THEN IN EACH BLOCK ARRAY(,,K)....                       
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
        DO 490 L = 1,NBLOKS 
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               BACKWARD ELIMINATION                                    
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
           K = NBLKS1-L 
           INCR = INCR-NRWBLK 
           DO 450 L1 = NRWEL1,NRWBLK 
              I = NRWBLK+NRWEL1-L1 
              IPLUSN = I+NRWTOP 
              LOOP = IPLUSN+1 
              INCRN = INCR+IPLUSN 
              DOTPRD = B(INCRN) 
              DO 440 J = LOOP,NCLBLK 
                 INCRJ = INCR+J 
                 DOTPRD = DOTPRD-ARRAY(I,J,K)*B(INCRJ) 
  440                       CONTINUE 
              B(INCRN) = DOTPRD 
              IPVTN = PIVOT(INCRN) 
              IF(INCRN.EQ.IPVTN)GO TO 445 
                 SWAP = B(INCRN) 
                 B(INCRN) = B(IPVTN) 
                 B(IPVTN) = SWAP 
  445                       CONTINUE 
  450                               CONTINUE 
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               BACKWARD MODIFICATION                                   
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
           INCRTP = INCR+NRWTOP 
           DO 460 J = NRWBK1,NCLBLK 
              INCRJ = INCR+J 
              XINCRJ = -B(INCRJ) 
              DO 455 I = 1,NROWEL 
                 INCRI = INCRTP+I 
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ 
  455                       CONTINUE 
  460                               CONTINUE 
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               BACKWARD SOLUTION                                       
!                                                                       
!          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
           DO 480 LL = 1,NROWEL 
              J = NRWBK1-LL 
              INCRJ = INCR+J 
              NRWELL = NRWEL1-LL 
              B(INCRJ) = B(INCRJ)/ARRAY(NRWELL,J,K) 
              IF(LL.EQ.NROWEL)GO TO 470 
                 XINCRJ = -B(INCRJ) 
                 LOOP = NROWEL-LL 
                 DO 465 I = 1,LOOP 
                    INCRI = INCRTP+I 
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ 
  465                             CONTINUE 
  470                                        CONTINUE 
  480                                                CONTINUE 
  490                                                     CONTINUE 
!                                                                       
!       ********************************************************        
!                                                                       
!          ***  IN TOPBLK FINISH WITH....                               
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
!               BACKWARD ELIMINATION                                    
!                                                                       
!       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
!                                                                       
        DO 520 L = 1,NRWTOP 
           I = NRWTP1-L 
           LOOP = I+1 
           DOTPRD = B(I) 
           DO 510 J = LOOP,NOVRLP 
              DOTPRD = DOTPRD-TOPBLK(I,J)*B(J) 
  510                 CONTINUE 
           B(I) = DOTPRD 
           IPVTI = PIVOT(I) 
           IF(I.EQ.IPVTI)GO TO 515 
                 SWAP = B(I) 
                 B(I) = B(IPVTI) 
                 B(IPVTI) = SWAP 
  515                    CONTINUE 
  520                         CONTINUE 
!                                                                       
!       RETURN FROM THE SOLUTION OF A.X = B.                            
        RETURN 
!                                                                       
!       IF JOB IS NON-ZERO, SOLVE TRANSPOSE(A)*X = B:                   
!                                                                       
  530      CONTINUE 
                                                                        
!       FIRST, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U).           
                                                                        
        DO 540 I = 1,NRWTOP 
           IPVTI = PIVOT(I) 
           IF ( I .NE. IPVTI ) THEN 
              SWAP = B(I) 
              B(I) = B(IPVTI) 
              B(IPVTI) = SWAP 
           ENDIF 
           BI = -B(I) 
           LOOP = I+1 
           DO 535 J = LOOP,NOVRLP 
              B(J) = B(J) + BI*TOPBLK(I,J) 
  535               CONTINUE 
  540                  CONTINUE 
                                                                        
!       IN EACH BLOCK, K = 1,..,NBLOKS:                                 
                                                                        
        INCR = NRWTOP 
        DO 590 K = 1,NBLOKS 
                                                                        
!          FIRST, THE FORWARD SOLUTION.                                 
                                                                        
           DO 550 J = 1,NROWEL 
              INCRJ = INCR + J 
              DO 545 I = 1,J-1 
                 B(INCRJ) = B(INCRJ) - ARRAY(I,NRWTOP+J,K)*B(INCR+I) 
  545                     CONTINUE 
              B(INCRJ) = B(INCRJ)/ARRAY(J,NRWTOP+J,K) 
  550                CONTINUE 
                                                                        
!           FORWARD MODIFICATION.                                       
                                                                        
            DO 570 I = 1,NOVRLP 
               INCRI = INCR + NROWEL + I 
               LOOP = NRWBLK + I 
               DO 560 J = 1,NROWEL 
                  INCRJ = INCR + J 
                  B(INCRI) = B(INCRI) - ARRAY(J,LOOP,K)*B(INCRJ) 
  560                       CONTINUE 
  570                              CONTINUE 
                                                                        
!           NOW, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U). THIS    
!           CORRESPONDS TO THE LOOP 540 ABOVE.                          
                                                                        
            INCR = INCR + NROWEL 
            DO 580 I = 1,NRWTOP 
               INCRI = INCR + I 
               IPVTI = PIVOT(INCRI) 
               IF ( INCRI .NE. IPVTI ) THEN 
                  SWAP = B(INCRI) 
                  B(INCRI) = B(IPVTI) 
                  B(IPVTI) = SWAP 
               ENDIF 
               LOOP = NROWEL + I 
               BI = -B(INCRI) 
               DO 575 J = I+1,NOVRLP 
                  INCRJ = INCR+J 
                  L = NRWBLK + J 
                  B(INCRJ) = B(INCRJ) + BI*ARRAY(LOOP,L,K) 
  575                       CONTINUE 
  580                              CONTINUE 
            INCR = INCR + NRWTOP 
  590          CONTINUE 
                                                                        
!       FINALLY, FINISH WITH NRWBOT SOLUTIONS:                          
                                                                        
        DO 600 J = 1,NRWBOT 
           INCRJ = INCR + J 
           DO 595 I = 1,J-1 
              B(INCRJ) = B(INCRJ) - BOTBLK(I,J+NRWTOP)*B(INCR+I) 
  595               CONTINUE 
           B(INCRJ) = B(INCRJ)/BOTBLK(J,J+NRWTOP) 
  600         CONTINUE 
                                                                        
                                                                        
!       NOW, THE BACKWARD PASS:                                         
                                                                        
                                                                        
!       FIRST, BACKWARD SOLUTION IN BOTBLK:                             
                                                                        
        INCRJ = INCR + NRWBOT 
        DO 610 J = 1,NRWBOT-1 
           INCRJ = INCRJ - 1 
           DO 605 I = NRWBOT-J+1,NRWBOT 
              INCRI = INCR + I 
              B(INCRJ) = B(INCRJ) - BOTBLK(I,NOVRLP-J)*B(INCRI) 
  605               CONTINUE 
                                                                        
           IF ( INCRJ .NE. PIVOT(INCRJ) ) THEN 
              SWAP = B(INCRJ) 
              B(INCRJ) = B(PIVOT(INCRJ)) 
              B(PIVOT(INCRJ)) = SWAP 
           ENDIF 
  610         CONTINUE 
                                                                        
!       NOW DO THE DEFERRED OPERATIONS IN BOTBLOK:                      
                                                                        
        DO 620 J = 1,NRWTOP 
           INCRJ = INCR - J + 1 
           DO 615 I = 1,NRWBOT 
              INCRI = INCR + I 
              B(INCRJ) = B(INCRJ) - BOTBLK(I,NRWTP1-J)*B(INCRI) 
  615               CONTINUE 
  620                  CONTINUE 
                                                                        
                                                                        
!       NOW, IN EACH BLOCK, K = NBLOKS,..,1:                            
        DO 800 K = NBLOKS,1,-1 
                                                                        
!          FIRST, THE BACKSUBSTITUIONS:                                 
                                                                        
           DO 630 J = 1,NRWTOP 
              INCRJ = INCR - J + 1 
              LOOP = NBKTOP - J + 1 
              DO 625 I = 1,J-1 
                 INCRI = INCR - I + 1 
                 B(INCRJ) = B(INCRJ) - ARRAY(NRWBLK-I+1,LOOP,K)*B(INCRI) 
  625                     CONTINUE 
              B(INCRJ) = B(INCRJ)/ARRAY(NRWBLK-J+1,LOOP,K) 
  630               CONTINUE 
                                                                        
!          THEN THE BACKWARD SOLUTION IN THE KTH BLOCK:                 
                                                                        
           DO 650 J = 1,NROWEL 
              INCRJ = INCR - NRWTOP -J + 1 
              DO 645 I = 1,J+NRWTOP-1 
                 INCRI = INCRJ + I 
                 B(INCRJ) = B(INCRJ) -                                  &
     &           ARRAY(NRWBLK-NRWTOP-J+1+I,NRWBLK-J+1,K)*B(INCRI)       
  645                     CONTINUE 
              IF ( INCRJ .NE. PIVOT(INCRJ) ) THEN 
                 SWAP = B(INCRJ) 
                 B(INCRJ) = B(PIVOT(INCRJ)) 
                 B(PIVOT(INCRJ)) = SWAP 
              ENDIF 
  650               CONTINUE 
                                                                        
!          NOW, THE DEFERRED OPERATIONS ON B:                           
                                                                        
           INCR = INCR - NRWBLK 
           DO 660 J = 1,NRWTOP 
              INCRJ = INCR + J - NRWTOP 
              DO 655 I = 1,NRWBLK 
                 INCRI = INCR + I 
                 B(INCRJ) = B(INCRJ) - ARRAY(I,J,K)*B(INCRI) 
  655                    CONTINUE 
  660                          CONTINUE 
  800                             CONTINUE 
                                                                        
!       FINALLY, THE LAST SET OF BACK-SUBSTITUTIONS IN TOPBLK:          
                                                                        
        DO 900 J = 1,NRWTOP 
           INCRJ = NRWTOP -J + 1 
           DO 850 I = INCRJ+1,NRWTOP 
              B(INCRJ) = B(INCRJ) - TOPBLK(I,INCRJ)*B(I) 
  850               CONTINUE 
           B(INCRJ) = B(INCRJ)/TOPBLK(INCRJ,INCRJ) 
  900         CONTINUE 
!                                                                       
!       RETURN FROM THE SOLUTION OF A-TRANSPOSE.X = B                   
                                                                        
        RETURN 
      end subroutine                                        
                                                                        
                                                                        
      subroutine lufac(n, ndim, a, ip, ier) 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension a(ndim,n), ip(n) 
      intrinsic abs 
!  blas: daxpy, dscal, dswap. idamax                                    
                                                                        
      parameter ( zero = 0.0_dp, one = 1.0_dp ) 
                                                                        
!  The subroutine lufac is a very simple code to compute the            
!  LU decomposition (with partial pivoting) of the n by n matrix a.     
                                                                        
!  The LU factors are overwritten on a.  The integer array ip           
!  reflects the pairwise interchanges performed.  note that ip(k)       
!  therefore does not give the index in the original array of           
!  the k-th pivot.                                                      
                                                                        
!  On exit, the error flag ier is zero when no zero pivots are          
!  encountered.  Otherwise, ier is equal to the index of the            
!  step at which a zero pivot occurred.                                 
                                                                        
      ier = 0 
      ip(n) = 0 
                                                                        
!  Begin loop over columns 1 through n-1.  k is the current             
!  column index.                                                        
                                                                        
      do 100 k = 1, n-1 
                                                                        
!  Find the row index ipiv of the element of largest magnitude in       
!  column k.                                                            
                                                                        
         ipiv = k-1 + idamax(n-k+1, a(k,k), 1) 
         piv = a(ipiv,k) 
         if (piv .eq. zero) then 
            ier = k 
            return 
         endif 
         ip(k) = ipiv 
                                                                        
!  Perform interchanges if necessary.                                   
                                                                        
         if (ipiv .ne. k) then 
            call dswap(n-k+1, a(ipiv,k), ndim, a(k,k), ndim) 
         endif 
                                                                        
!  Save the (negative) multipliers in the subdiagonal elements of       
!  column k.                                                            
                                                                        
         call dscal(n-k, (-one/piv), a(k+1,k), 1) 
                                                                        
!  Update the remaining matrix.  Note that a(i,k) now contains          
!  the negative multipliers.                                            
                                                                        
         do 50 j = k+1, n 
            call daxpy(n-k, a(k,j), a(k+1,k), 1, a(k+1,j), 1) 
   50    continue 
                                                                        
!  End of loop over columns.                                            
                                                                        
  100 continue 
      if (a(n,n).eq.zero) ier = n 
      return 
      end subroutine                                        
                                                                        
      subroutine lusol (n, ndim, a, ip, b, x) 
      implicit integer (i-n)
      implicit real(dp) (a-h, o-z) 
      dimension a(ndim,n), ip(n), b(n), x(n) 
                                                                        
!  blas:  daxpy, dcopy                                                  
                                                                        
!  The subroutine lusol is a simple-minded routine to solve a           
!  linear system whose LU factors have been computed by lufac.          
!  On entry, the matrix a should contain the LU factors, and            
!  ip should contain the interchange array constructed by lufac.        
                                                                        
                                                                        
!  Copy the right-hand side b into x.                                   
                                                                        
      call dcopy (n, b, 1, x, 1) 
                                                                        
!  Forward solution with l (unit lower-triangular factor), which        
!  is stored in the strict lower triangle of a.                         
                                                                        
      do 20 k = 1, n-1 
         ipiv = ip(k) 
         if (ipiv .ne. k) then 
            tem = x(ipiv) 
            x(ipiv) = x(k) 
            x(k) = tem 
         endif 
         call daxpy ( n-k, x(k), a(k+1,k), 1, x(k+1), 1 ) 
   20 continue 
                                                                        
!  Backward solution with u (upper-triangular factor), which is stored  
!  in the upper triangle of a.                                          
                                                                        
      do 40 kb = n, 1, -1 
         x(kb) = x(kb)/a(kb,kb) 
         call daxpy(kb-1, (-x(kb)), a(1,kb), 1, x(1), 1) 
   40 continue 
                                                                        
      return 
      end subroutine                                        
                                                                        
      subroutine dcopy ( n, x, incx, y, incy ) 
      implicit real(dp) (a-h,o-z) 
      integer            n, incx, incy 
      dimension          x( * ), y( * ) 
                                                                        
!  dcopy  performs the operation                                        
!                                                                       
!     y := x                                                            
!                                                                       
!  nag fortran 77 version of the blas routine dcopy .                   
!  nag fortran 77 o( n ) basic linear algebra routine.                  
!                                                                       
!  -- written on 26-november-1982.                                      
!     sven hammarling, nag central office.                              
                                                                        
      integer            i     , ix    , iy 
                                                                        
      if( n.lt.1 )return 
                                                                        
      if( ( incx.eq.incy ).and.( incy.gt.0 ) )then 
         do 10, iy = 1, 1 + ( n - 1 )*incy, incy 
            y( iy ) = x( iy ) 
   10    continue 
      else 
         if( incx.ge.0 )then 
            ix = 1 
         else 
            ix = 1 - ( n - 1 )*incx 
         end if 
         if( incy.gt.0 )then 
            do 20, iy = 1, 1 + ( n - 1 )*incy, incy 
               y( iy ) = x( ix ) 
               ix      = ix + incx 
   20       continue 
         else 
            iy = 1 - ( n - 1 )*incy 
            do 30, i = 1, n 
               y( iy ) = x( ix ) 
               iy      = iy + incy 
               ix      = ix + incx 
   30       continue 
         end if 
      end if 
                                                                        
      return 
                                                                        
!     end of dcopy .                                                    
                                                                        
      end subroutine                                        
                                                                        
      subroutine daxpy ( n, alpha, x, incx, y, incy ) 
      implicit real(dp) (a-h,o-z) 
      integer            n, incx, incy 
      dimension   x( * ), y( * ) 
                                                                        
!  daxpy  performs the operation                                        
!                                                                       
!     y := alpha*x + y                                                  
!                                                                       
!                                                                       
!  modified nag fortran 77 version of the blas routine daxpy .          
!                                                                       
!  -- written on 3-september-1982.                                      
!     sven hammarling, nag central office.                              
                                                                        
      integer            i     , ix    , iy 
      parameter        ( zero  = 0.0+0 ) 
                                                                        
      if( n    .lt.1    )return 
      if( alpha.eq.zero )return 
                                                                        
      if( ( incx.eq.incy ).and.( incx.gt.0 ) )then 
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx 
            y( ix ) = alpha*x( ix ) + y( ix ) 
   10    continue 
      else 
         if( incy.ge.0 )then 
            iy = 1 
         else 
            iy = 1 - ( n - 1 )*incy 
         end if 
         if( incx.gt.0 )then 
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx 
               y( iy ) = alpha*x( ix ) + y( iy ) 
               iy      = iy + incy 
   20       continue 
         else 
            ix = 1 - ( n - 1 )*incx 
            do 30, i = 1, n 
               y( iy ) = alpha*x( ix ) + y( iy ) 
               ix      = ix + incx 
               iy      = iy + incy 
   30       continue 
         end if 
      end if 
      return 
                                                                        
!     end of daxpy .                                                    
                                                                        
      end subroutine                                        
                                                                        
      real(dp) function ddot  ( n, x, incx, y, incy ) 
      implicit real(dp) (a-h,o-z) 
      integer           n, incx, incy 
      dimension         x( * ), y( * ) 
                                                                        
!  ddot   returns the value                                             
!                                                                       
!     ddot   = x'y                                                      
!                                                                       
!                                                                       
!  modified nag fortran 77 version of the blas routine ddot  .          
!                                                                       
!  -- written on 21-september-1982.                                     
!     sven hammarling, nag central office.                              
                                                                        
      integer             i     , ix    , iy 
      parameter         ( zero  = 0.0_dp ) 
                                                                        
      sum = zero 
      if( n.ge.1 )then 
         if( ( incx.eq.incy ).and.( incx.gt.0 ) )then 
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx 
               sum = sum + x( ix )*y( ix ) 
   10       continue 
         else 
            if( incy.ge.0 )then 
               iy = 1 
            else 
               iy = 1 - ( n - 1 )*incy 
            end if 
            if( incx.gt.0 )then 
               do 20, ix = 1, 1 + ( n - 1 )*incx, incx 
                  sum = sum + x( ix )*y( iy ) 
                  iy  = iy  + incy 
   20          continue 
            else 
               ix = 1 - ( n - 1 )*incx 
               do 30, i = 1, n 
                  sum = sum + x( ix )*y( iy ) 
                  ix  = ix  + incx 
                  iy  = iy  + incy 
   30          continue 
            end if 
         end if 
      end if 
                                                                        
      ddot   = sum 
      return 
                                                                        
!     end of ddot  .                                                    
                                                                        
      end function
                                                                        
      subroutine dscal ( n, alpha, x, incx ) 
      implicit real(dp) (a-h,o-z) 
      integer          n, incx 
      dimension        x( * ) 
                                                                        
!  dscal  performs the operation                                        
!                                                                       
!     x := alpha*x                                                      
!                                                                       
!                                                                       
!  modified nag fortran 77 version of the blas routine dscal .          
!                                                                       
!  -- written on 26-november-1982.                                      
!     sven hammarling, nag central office.                              
                                                                        
      integer            ix 
      parameter        ( one   = 1.0_dp, zero  = 0.0_dp ) 
                                                                        
      if( n.ge.1 )then 
         if( alpha.eq.zero )then 
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx 
               x( ix ) = zero 
   10       continue 
         else if( alpha.eq.( -one ) )then 
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx 
               x( ix ) = -x( ix ) 
   20       continue 
         else if( alpha.ne.one )then 
            do 30, ix = 1, 1 + ( n - 1 )*incx, incx 
               x( ix ) = alpha*x( ix ) 
   30       continue 
         end if 
      end if 
                                                                        
      return 
                                                                        
!     end of dscal .                                                    
                                                                        
      end subroutine                                        
                                                                        
      subroutine dswap ( n, x, incx, y, incy ) 
      implicit real(dp) (a-h,o-z) 
      integer      n, incx, incy 
      dimension    x( * ), y( * ) 
                                                                        
!  dswap  performs the operations                                       
!                                                                       
!     temp := x,   x := y,   y := temp.                                 
!                                                                       
!                                                                       
!  modified nag fortran 77 version of the blas routine dswap .          
!                                                                       
!  -- written on 26-november-1982.                                      
!     sven hammarling, nag central office.                              
                                                                        
      integer            i     , ix    , iy 
                                                                        
      if( n.lt.1 )return 
                                                                        
      if( ( incx.eq.incy ).and.( incy.gt.0 ) )then 
         do 10, iy = 1, 1 + ( n - 1 )*incy, incy 
            temp    = x( iy ) 
            x( iy ) = y( iy ) 
            y( iy ) = temp 
   10    continue 
      else 
         if( incx.ge.0 )then 
            ix = 1 
         else 
            ix = 1 - ( n - 1 )*incx 
         end if 
         if( incy.gt.0 )then 
            do 20, iy = 1, 1 + ( n - 1 )*incy, incy 
               temp    = x( ix ) 
               x( ix ) = y( iy ) 
               y( iy ) = temp 
               ix      = ix + incx 
   20       continue 
         else 
            iy = 1 - ( n - 1 )*incy 
            do 30, i = 1, n 
               temp    = x( ix ) 
               x( ix ) = y( iy ) 
               y( iy ) = temp 
               iy      = iy + incy 
               ix      = ix + incx 
   30       continue 
         end if 
      end if 
                                                                        
      return 
                                                                        
!     end of dswap .                                                    
                                                                        
      end subroutine                                        
                                                                        
      function idamax( n, x, incx ) 
        implicit real(dp) (a-h,o-z) 
        integer         n, incx 
        dimension       x( * ) 
        integer :: idamax
!  idamax returns the smallest value of i such that                     
!                                                                       
!     abs( x( i ) ) = max( abs( x( j ) ) )                              
!                      j                                                
!                                                                       
!  nag fortran 77 version of the blas routine idamax.                   
!  nag fortran 77 o( n ) basic linear algebra routine.                  
!                                                                       
!  -- written on 31-may-1983.                                           
!     sven hammarling, nag central office.                              
                                                                        
      intrinsic           abs 
      integer             i     , imax  , ix 
                                                                        
      if( n.lt.1 )then 
         idamax = 0 
         return 
      end if 
                                                                        
      imax = 1 
      if( n.gt.1 )then 
         xmax = abs( x( 1 ) ) 
         ix   = 1 
         do i = 2, n 
            ix = ix + incx 
            if( xmax.lt.abs( x( ix ) ) )then 
               xmax = abs( x( ix ) ) 
               imax = i 
            end if 
         enddo
      end if
                                                                        
      idamax = imax 
      return 
                                                                        
!     end of idamax.                                                    
                                                                        
    end function idamax


      subroutine dload ( n, const, x, incx ) 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension  x(*) 
!                                                                       
!  dload  performs the operation                                        
!                                                                       
!     x = const*e,   e' = ( 1  1 ... 1 ).                               
!                                                                       
!                                                                       
!  nag fortran 77 o( n ) basic linear algebra routine.                  
!                                                                       
!  -- written on 22-september-1983.                                     
!     sven hammarling, nag central office.                              
!                                                                       
      parameter        ( zero = 0.0_dp ) 
                                                                        
      if( n.lt.1 )return 
                                                                        
      if( const.ne.zero )then 
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx 
            x( ix ) = const 
   10    continue 
      else 
         do 20, ix = 1, 1 + ( n - 1 )*incx, incx 
            x( ix ) = zero 
   20    continue 
      end if 
                                                                        
      return 
                                                                        
!     end of dload .                                                    
      end subroutine                                        
                                                                        
                                                                        
                                                                        
      subroutine maxpy ( nrow, ncol, alpha, xmat, nrowy, ymat ) 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension xmat(nrow, ncol), ymat(nrowy, ncol) 
                                                                        
!  Subroutine maxpy takes as input the scalar alpha and two matrices,   
!  xmat and ymat.  xmat has declared row dimension nrow, and            
!  ymat has declared row dimension nrowy, but both are                  
!  conceptually nrow by ncol.                                           
!  On output, (new ymat) is alpha*xmat+ (old ymat), by analogy          
!  with the vector blas routine saxpy.                                  
                                                                        
      do 100 j = 1, ncol 
      do 100 i = 1, nrow 
         ymat(i,j) = ymat(i,j) + alpha*xmat(i,j) 
  100 continue 
      return 
      end subroutine                                        
                                                                        
                                                                        
      subroutine matcop( nrow1, nrow2, nrow, ncol, xmat1, xmat2 ) 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension xmat1(nrow1, ncol), xmat2(nrow2, ncol) 
                                                                        
!  Given 2 matrices xmat1 and xmat2, where xmat1 has declared           
!  row dimension nrow1, xmat2 has declared row dimension nrow2,         
!  and both have column dimension ncol, the routine matcop copies       
!  rows 1 through nrow, and columns 1 through ncol from xmat1 into      
!  xmat2.                                                               
                                                                        
                                                                        
      if (nrow .le. 0 .or. ncol .le. 0) return 
      do 100 j = 1, ncol 
      do 100 i = 1, nrow 
           xmat2(i,j) = xmat1(i,j) 
  100 continue 
      return 
      end subroutine                                        
                                                                        
      subroutine mtload( nrow, ncol, const, nrowx, xmat ) 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension xmat(nrowx, ncol) 
                                                                        
!  mtload sets elements 1 through nrow, 1 through ncol, of the          
!  matrix xmat (whose declared row dimension is nrowx) to the           
!  scalar value const.                                                  
                                                                        
      if (nrow .le. 0 .or. ncol .le. 0)  return 
      do 100 j = 1, ncol 
      do 100 i = 1, nrow 
         xmat(i,j) = const 
  100 continue 
      return 
      end subroutine                                        
                                                                        
                                                                        
                                                                        
      subroutine mssq  ( nrow, ncol, xmat, scale, sumsq ) 
      implicit integer (i-n)
      implicit real(dp) (a-h,o-z) 
      dimension   xmat(nrow, *) 
                                                                        
!  Given the nrow by ncol matrix xmat, mssq returns values              
!  scale and sumsq such that                                            
!     (scale**2) * sumsq = sum of squares of xmat(i,j),                 
!  where  scale = max  abs(xmat(i,j)).                                  
                                                                        
!  mssq is a stripped-down matrix version of the blas routine sssq.     
                                                                        
      intrinsic    abs 
      parameter   ( one   = 1.0_dp, zero  = 0.0_dp ) 
                                                                        
      scale = zero 
      sumsq = one 
      if( nrow.ge.1 .and. ncol .ge. 1) then 
         do 10 i = 1, nrow 
         do 10 j = 1, ncol 
            if( xmat(i,j) .ne. zero )then 
               absxij = abs(xmat(i,j)) 
               if( scale .lt. absxij ) then 
                  sumsq = one + sumsq* (scale/absxij)**2 
                  scale = absxij 
               else 
                  sumsq = sumsq + (absxij/scale)**2 
               end if 
            end if 
   10    continue 
      end if 
      return 
      end subroutine                                        
                                                                        
      subroutine dssq  ( n, x, incx, scale, sumsq )
      implicit real(dp) (a-h,o-z) 
      integer            n, incx,ix
      dimension   x( * ) 
                                                                        
!  Given the n-vector x, dssq returns values scale and sumsq such that  
!     (scale**2) * sumsq = sum of squares of x(i),                      
!  where  scale = max  abs(x(i)).                                       
                                                                        
      intrinsic          abs 
      parameter        ( one   = 1.0_dp, zero  = 0.0_dp ) 
                                                                        
      scale = zero 
      sumsq = one 
      if( n.ge.1 )then 
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx 
            if( x( ix ).ne.zero )then 
               absxi = abs( x( ix ) ) 
               if( scale.lt.absxi )then 
                  sumsq = one   + sumsq*( scale/absxi )**2 
                  scale = absxi 
               else 
                  sumsq = sumsq + ( absxi/scale )**2 
               end if 
            end if 
   10    continue 
      end if 
                                                                        
      return 
                                                                        
!     end of dssq  .                                                    
                                                                        
      end subroutine                                        


    end module sepbvp

