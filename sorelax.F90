module sorelax
  use eomrelax, only : kp
  use eomrelax, only : sysrelax, action, cvgtest, constraint
#ifdef WITHIO
  use iotools, only : allwrite
#endif


  implicit none

  private
  

  public kp, relax

contains
  
  subroutine relax(itmax,conv,sorfacin,y,r,sortarg,dumpFiles)
    implicit none
               

    integer ,intent(in) :: itmax
    real(kp), intent(in) :: conv  
    real(kp), dimension(:,:), intent(inout) :: y
    real(kp), dimension(:), intent(in) :: r
    real(kp), dimension(:), intent(in) :: sorfacin
    real(kp), dimension(:), intent(in), optional :: sortarg
    real(kp), dimension(size(sorfacin,1)) :: sorfac
    logical, intent(in), optional :: dumpFiles

    real(kp), dimension(size(y,1)) :: dS,dy
    real(kp), dimension(size(y,1),size(y,1)) :: nabladS
      
    real(kp), dimension(size(y,1),size(y,2)) :: dvecy

    real(kp) :: Snow, Sold, DeltaS,SigmaS,varS
            
    real(kp) :: erreur,erreur2
    integer :: iter,itmean
    
    integer :: i,j
    integer :: n,ne,kmin,kmax,k
    
    logical, save :: dump = .false.
    logical, parameter :: display=.true.
    logical, parameter :: use_cvgtest=.false.

    logical, save :: dynsor = .false.
    integer, parameter :: itdynamic = 100

    integer :: itaction,itdisplay,iterdump,itconstraint

    integer :: nmode
    integer, parameter :: ngrid=2

    if (present(dumpFiles)) then
       dump = dumpFiles
    endif

    if (present(sortarg)) then
       dynsor = .true.
    endif
    
    sorfac = sorfacin

    itdisplay=250
    itaction=25
    iterdump=itdisplay
    itmean=0

    itconstraint=itmax

    n=size(r,1)
    ne=size(y,1)

    if (mod(n,2).ne.0) then
       stop 'Number of points has to be even'
    endif

    kmin=lbound(r,1)
    kmax=ubound(r,1)

    print *,'kmin kmax',kmin,kmax

    if (n.ne.size(y,2)) then
       write(*,*)'__relax: y and r sizes differ !'
       stop
    endif

    iter=0    
    erreur=conv+1._kp
    erreur2=1._kp

    write(*,*)
    write(*,*)'Relaxation in progress...'
    write(*,*)
    
    Snow=action(y,r)
    Sold=Snow
    DeltaS = huge(1._kp)

    do while ((abs(erreur).ge.conv).and.(iter.le.itmax))

       iter=iter+1
  
! first side boundary conditions
       k=kmin
       call sysrelax(k,kmin,kmax,y,r,dS,nabladS)
       
! relaxation on two grids         
! the fields at iter np1=n+1 are corrected from those at iter n such as:
! dS(Fnp1) = dS(Fn+FCorr) = dS(Fn) + NabladS()*FCorr = 0
! FCorr = -NabladS()^(-1) . dS(Fn)

!1grid
!       do k=kmin+1,kmax-1,1            
!          call sysrelax2(k,kmin,kmax,y,r,dS,nabladS,ind)
!          dy=dysolve(dS,nabladS)          
!          y(ind(:),k)=y(ind(:),k) - sorfac(ind(:))*dy(ind(:))
!       enddo
            
!2grids
!       do k=kmin+1,kmax-2,2            
!             call sysrelax2(k,kmin,kmax,y,r,dS,nabladS,ind)
!             dy=dysolve(dS,nabladS)          
!             y(ind(:),k)=y(ind(:),k) - sorfac(ind(:))*dy(ind(:))
!          enddo
          
!          do k=kmin+2,kmax-1,2  
!             call sysrelax2(k,kmin,kmax,y,r,dS,nabladS,ind)          
!             dy=dysolve(dS,nabladS)          
!             y(ind(:),k)=y(ind(:),k) - sorfac(ind(:))*dy(ind(:))
!          enddo


       if (dynsor) then
          if (iter.gt.itdynamic) then
             sorfac(:) = sortarg(:) + (sorfacin(:) - sortarg(:))/(real(iter,kp)/real(itdynamic,kp))
          endif
       endif
                           
       do nmode=1,ngrid
          do k=kmin+nmode, kmax-(ngrid-nmode+1), ngrid
             call sysrelax(k,kmin,kmax,y,r,dS,nabladS)
             dy=dysolve(dS,nabladS)          
             y(:,k)=y(:,k) - sorfac(:)*dy(:)
          enddo                   
       enddo

       
! second side boundary conditions
       k=kmax
       call sysrelax(k,kmin,kmax,y,r,dS,nabladS)   
       
      
       SigmaS = action(y,r)

       itmean = itmean + 1

       Snow = (Snow*real(itmean-1,kp) + SigmaS)/real(itmean,kp)
       varS = (varS*(real(itmean-1,kp)) + (SigmaS-Sold)**2)/real(itmean,kp)
       DeltaS = sqrt(varS)

       if ((mod(iter,itaction).eq.0).or.(iter.eq.1)) then         
          itmean=0
          erreur2=cvgtest(y,r)
          Sold=Snow
       endif
      

       erreur=DeltaS/abs(Snow)

          
       if ((display).and.(mod(iter,itdisplay).eq.0)) then
         
          
          write(*,*)
          write(*,*)
          write(*,*)
          write(*,*)
          write(*,*)
          write(*,*)
          write(*,*)
          write(*,*)'---------------------------------------------------------'
          write(*,*)'iter itmax = ',iter,itmax
          write(*,*)'S = ',SigmaS
          write(*,*)'DeltaS/S = ',erreur
          write(*,*)'constraint eq = ',erreur2
          write(*,*)
          write(*,*)'npoints  =',n
          write(*,*)'rmin rmax = ',r(kmin),r(kmax)
          write(*,*)
          write(*,*)'S.O.R factors =',sorfac
          write(*,*)'stationary precision required ',conv                 
          write(*,*)'---------------------------------------------------------'
          write(*,*)
          write(*,*)
          write(*,*)
          write(*,*)
          if (dump.and.(mod(iter,iterdump).eq.0)) then
          else
             write(*,*)
          endif
          if (mod(iter,itconstraint).eq.0) then
          else
             write(*,*)
          endif
          write(*,*)
       endif
      
       if (mod(iter,itconstraint).eq.0) then
          call constraint(y,r)
       endif
       
#ifdef WITHIO
       if ((dump).and.(mod(iter,iterdump).eq.0)) then     
          if (ne.eq.5) then
             call allwrite('liverelax.dat',r,y(1,:),y(2,:) & 
                  ,y(3,:),y(4,:),y(5,:))
          elseif (ne.eq.4) then
             call allwrite('liverelax.dat',r,y(1,:),y(2,:) & 
                  ,y(3,:),y(4,:))
          elseif (ne.eq.3) then
             call allwrite('liverelax.dat',r,y(1,:),y(2,:) &
                  ,y(3,:))
          elseif (ne.eq.2) then
             call allwrite('liverelax.dat',r,y(1,:),y(2,:))
          elseif (ne.eq.1) then
             call allwrite('liverelax.dat',r,y(1,:))
          else
             stop 'livedump not implemented'
          endif

#endif
       endif


       
    end do
    
    if (erreur.le.conv) then
       write(*,*)
       write(*,*)'stationary precision reached',erreur
       write(*,*)'number of iterations',iter   
    else
       write(*,*)'number of iterations reached',iter
       write(*,*)'stationaty precision = ',erreur
    endif

  end subroutine relax


  function dysolve(dS,nabladS)

    implicit none
    real(kp), dimension(:),intent(in) :: dS
    real(kp), dimension(:,:),intent(in) :: nabladS
    integer, dimension(size(dS,1)) :: ipiv
    integer :: ne,info
    real(kp), dimension(size(dS,1)) :: dysolve

!no choice for lapack
    integer, parameter :: dp = kind(1._8)
    real(dp), dimension(size(dS,1)) :: dStemp
    real(dp), dimension(size(dS,1),size(dS,1)) :: nabladS_LU

    integer :: i,j
    real(kp) :: buffer
    logical, parameter :: do_test = .false.


    ne=size(dS,1)

    if (ne.eq.1) then

       dysolve = dS(1)/nabladS(1,1)

    else

       nabladS_LU=nabladS
       dStemp=dS 
    
       call dgesv(ne,1,nabladS_LU,ne,ipiv,dStemp,ne,info)

       if (info.eq.0) then
          dysolve=dStemp
       else
          write(*,*)'info= ',info
          dysolve=0._kp
          read(*,*)
       endif
    
    endif

    if (do_test) then
       write(*,*)
       write(*,*)'ne= info= ',ne,info              
       write(*,*)'dysolve= ',dysolve
       write(*,*)
       do i=1,ne                    
          buffer=0._kp
          do j=1,ne
             buffer=nabladS(i,j)*dysolve(j)+buffer
          enddo
          write(*,*)'residual i: dS-nabladS*dy= ',i,dS(i)-buffer
       enddo
       read(*,*)
    endif
    
  end function dysolve


end module sorelax
