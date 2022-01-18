module eomrelax
  use precision, only : kp
  use chamber, only : stress_trace, rChamb, Mpl
  use gravity, only : potential, deriv_potential, deriv_second_potential
  use gravity, only : conformal, deriv_conformal, deriv_second_conformal

  implicit none

  private




  public kp
  public sysrelax, action, cvgtest, constraint


contains

  subroutine sysrelax(k,kmin,kmax,y,x,dS,nabladS) 

    !use with sorelax to provide the equations of motion, their linearisation
    !the action, additional constraints and additional convergence test

    implicit none 
    integer,intent(in) :: k,kmin,kmax
    real(kp), dimension(:,:),intent(inout) :: y
    real(kp), dimension(:),intent(in) :: x
    real(kp), dimension(:),intent(out) :: dS
    real(kp), dimension(:,:),intent(out) :: nabladS

    real(kp) :: buffer

    integer :: ne

    real(kp) :: h
    real(kp) :: xk,xkm1,xkp1

    real(kp) :: phik
    real(kp) :: phikm1, phikp1, phikm12, phikp12
    real(kp) :: phi1,phi2
    real(kp) :: phinm1,phin

    real(kp) :: Uphik, dUphik, d2Uphik
    real(kp) :: Omphik, dOmphik, d2Omphik
    real(kp) :: tracexk,rk

    ne=size(y,1)

    if (k.eq.kmin) then
       ! boundary conditions in rmin

       phi2=y(1,2)
       phi1=phi2
       y(1,1)=phi1

    else if (k.ge.kmax) then          
       ! boundary conditions in rmax     


       phinm1=y(1,kmax-1)
       phin=phinm1
       y(1,kmax)=phin

    else

       !discretization

       xk=x(k)
       rk = 1./xk
       xkm1=x(k-1)
       xkp1=x(k+1)

       h=(xkp1-xkm1)/2.

       phik=y(1,k)
       phikm1=y(1,k-1)
       phikp1=y(1,k+1)

       phikm12=(phik+phikm1)/2._kp
       phikp12=(phik+phikp1)/2._kp

       !the dilaton potential and its derivatives wrt phi
       !the 2nd derivative is needed for the quasilinearisation

       !potential for a massive dilaton
       Uphik = potential(phik)
       dUphik = deriv_potential(phik)
       d2Uphik = deriv_second_potential(phik)

       Omphik = conformal(phik)
       dOmphik = deriv_conformal(phik)
       d2Omphik = deriv_second_conformal(phik)

       tracexk = stress_trace(rk)

       !dS/dphi=0

!       dS(1) = (phik-phikm1)/h + (phik-phikp1)/h &
!            + h*(dUphik - 0.5_kp*dOmphik*tracexk) / xk**4

       dS(1) = (phik-phikm1) + (phik-phikp1) &
            + h*h*(dUphik - 0.5_kp*dOmphik*tracexk) / xk**4

       !jacobian of dS

!       nabladS(1,1) = 2._kp/h &
!            + h*(d2Uphik - 0.5_kp* d2Omphik*tracexk)/ xk**4

       nabladS(1,1) = 2._kp &
            + h*h*(d2Uphik - 0.5_kp* d2Omphik*tracexk)/ xk**4

    endif
       
  end subroutine sysrelax



  function action(y,x)
    implicit none

    real(kp) :: action
    real(kp), dimension(:),intent(in) :: x
    real(kp), dimension(:,:),intent(in) :: y

    integer :: i,imin,imax,ne,n
    real(kp) :: h, phi, Uphi, Omphi, strace
    real(kp) :: phidot,ri

    n=size(x,1)
    ne=size(y,1)

    imin=lbound(x,1)
    imax=ubound(x,1)

    action=0._kp

    do i=imin+1,imax-1

       h=(x(i+1)-x(i-1))/2._kp

       phi = y(1,i)
       Uphi = potential(phi)
       Omphi = conformal(phi)
       ri = 1._kp/x(i)
       strace = stress_trace(ri)*rChamb

       phidot = (y(1,i+1)-y(1,i-1))/2._kp

!       action = action - (0.5_kp*phidot**2/h + h*(Uphi-0.5*Omphi*strace)/x(i)**4)

       action = action - (0.5_kp*phidot**2 + h*h*(Uphi-0.5_kp*Omphi*strace)/x(i)**4)

    enddo


  end function action


  function cvgtest(y,x)
    implicit none
    real(kp) :: cvgtest
    real(kp), dimension(:),intent(in) :: x
    real(kp), dimension(:,:),intent(in) :: y

    cvgtest=0.

  end function cvgtest



  subroutine constraint(y,x)
    implicit none
    real(kp), dimension(:),intent(in) :: x
    real(kp), dimension(:,:),intent(inout) :: y 

  end subroutine constraint


end module eomrelax
