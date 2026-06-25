!   This file is part of futilities
!
!   Copyright (C) 2026 C. Ringeval
!   
!   futilities is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   futilities is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with futilities.  If not, see <https://www.gnu.org/licenses/>.

module fflint
  use, intrinsic :: iso_c_binding
  use precision, only : fdp, pidp
  implicit none

  private

  
  interface

     subroutine allocate_acb_t(n, ptr) bind(C)
       import C_INT, C_PTR
       integer(C_INT), value :: n
       type(C_PTR) :: ptr
     end subroutine allocate_acb_t
     
     subroutine free_acb_t(n, ptr) bind(C)
       import C_INT, C_PTR
       integer(C_INT), value :: n
       type(C_PTR) :: ptr
     end subroutine free_acb_t     

     function allocated_bytes_acb_t(n,ptr) bind (C)
       import C_INT, C_PTR
       integer(C_INT) :: allocated_bytes_acb_t
       integer(C_INT), value :: n
       type(C_PTR) :: ptr
     end function allocated_bytes_acb_t
     
     subroutine initialize_acb_t_real(n, x, ptr) bind(C)
       import C_INT, C_DOUBLE, C_PTR
       integer(C_INT), value :: n
       real(C_DOUBLE), dimension(n) :: x
       type(C_PTR) :: ptr
     end subroutine initialize_acb_t_real

     subroutine initialize_acb_t_cmpx(n, x, y, ptr) bind(C)
       import C_INT, C_DOUBLE, C_PTR
       integer(C_INT), value :: n
       real(C_DOUBLE), dimension(n) :: x, y
       type(C_PTR) :: ptr
     end subroutine initialize_acb_t_cmpx

     subroutine get_elliptic_thetas(theta1,theta2,theta3,theta4,z,tau,prec) bind(C)
       import C_DOUBLE, C_INT
       real(C_DOUBLE), dimension(2), intent(out) :: theta1,theta2,theta3,theta4
       real(C_DOUBLE), dimension(2), intent(in) :: z, tau
       integer(C_INT), value :: prec
     end subroutine get_elliptic_thetas
     
     function arb_fpwrap_double_exp_integral_ei(res,x,flags) bind(C)
       import C_INT, C_DOUBLE
       integer(C_INT) :: arb_fpwrap_double_exp_integral_ei
       real(C_DOUBLE) :: res
       real(C_DOUBLE), value :: x
       integer(C_INT), value :: flags
     end function arb_fpwrap_double_exp_integral_ei

     function arb_fpwrap_cdouble_exp_integral_ei(res,x,flags) bind(C)
       import C_INT, C_DOUBLE
       integer(C_INT) :: arb_fpwrap_cdouble_exp_integral_ei
       complex(C_DOUBLE) :: res
       complex(C_DOUBLE), value :: x
       integer(C_INT), value :: flags
     end function arb_fpwrap_cdouble_exp_integral_ei     

     function arb_fpwrap_double_exp_integral_e(res,s,x,flags) bind(C)
       import C_INT, C_DOUBLE
       integer(C_INT) :: arb_fpwrap_double_exp_integral_e
       real(C_DOUBLE) :: res
       real(C_DOUBLE), value :: s,x
       integer(C_INT), value :: flags
     end function arb_fpwrap_double_exp_integral_e

     function arb_fpwrap_cdouble_exp_integral_e(res,s,x,flags) bind(C)
       import C_INT, C_DOUBLE
       integer(C_INT) :: arb_fpwrap_cdouble_exp_integral_e
       complex(C_DOUBLE) :: res
       complex(C_DOUBLE), value :: s,x
       integer(C_INT), value :: flags
     end function arb_fpwrap_cdouble_exp_integral_e    


     
  end interface

  interface exp_integral_ei
     procedure arb_fpwrap_cdouble_exp_integral_ei, arb_fpwrap_double_exp_integral_ei
  end interface exp_integral_ei

    interface exp_integral_en
     procedure arb_fpwrap_cdouble_exp_integral_e, arb_fpwrap_double_exp_integral_e
  end interface exp_integral_en

  
  integer, parameter :: ntheta = 4

  logical, parameter :: display = .false.
  
  
  public free_acb_t, allocate_acb_t, allocated_bytes_acb_t
  public initialize_acb_t_real, initialize_acb_t_cmpx

!from flint  
  public :: ntheta, exp_integral_ei, exp_integral_en
  public :: elliptic_thetas, deriv_elliptic_thetas

!brute force series  
  public :: lambert_sine_series, deriv_elliptic_lnthetas
  public :: elliptic_thetas_fourier_series, deriv_elliptic_thetas_fourier_series
  
contains


  
  function elliptic_thetas(u,lnq)
    implicit none
    complex(fdp), dimension(ntheta) :: elliptic_thetas
    complex(fdp), intent(in) :: u,lnq
    
    real(C_DOUBLE), dimension(2), save :: z,tau
    real(C_DOUBLE), dimension(2), save :: t1, t2, t3, t4
!$omp threadprivate(z,tau,t1,t2,t3,t4)

!64bits precision 2^(-53)
    integer(C_INT) :: prec = 53
    
    z(1) = real(u/pidp,C_DOUBLE)
    z(2) = real(aimag(u/pidp),C_DOUBLE)

    tau(1) = real(aimag(lnq/pidp),fdp)
    tau(2) = -real(lnq/pidp,fdp)
    
    call get_elliptic_thetas(t1,t2,t3,t4,z,tau,prec)

    elliptic_thetas(1) = cmplx(t1(1),t1(2),fdp)
    elliptic_thetas(2) = cmplx(t2(1),t2(2),fdp)
    elliptic_thetas(3) = cmplx(t3(1),t3(2),fdp)
    elliptic_thetas(4) = cmplx(t4(1),t4(2),fdp)
    
  end function elliptic_thetas

  
! Without and with the (-1)^n, we compute:
!  
!   sum(n=1 to infty) of (-1)^n * q^2n/(1-q^2n) * sin(2nu)
!
!  
  function lambert_sine_series(u,q,tol)
    implicit none
    complex(fdp), dimension(4) :: lambert_sine_series
    complex(fdp), intent(in) :: u,q
    real(fdp), intent(in), optional :: tol

    complex(fdp) :: qn, twonu
    complex(fdp), dimension(4) :: term
    real(fdp) :: abserr, maxerr
    real(fdp) :: monen
    
    integer :: counter
    logical, parameter :: debug = .false.

    if (present(tol)) then
       maxerr = tol
    else
       maxerr = epsilon(1._fdp)
    endif
    
    lambert_sine_series = cmplx(0._fdp,0._fdp,fdp)
    term = cmplx(0._fdp,0._fdp,fdp)
    abserr = 1._fdp
    
    qn = q
    twonu = 2._fdp*u
    monen = -1._fdp
    counter = 0
    
    do while (abserr.gt.maxerr)
       term(4) = qn/(1._fdp - qn*qn) * sin(twonu)
       term(1) = term(4)*qn
       term(2) = term(1) * monen
       term(3) = term(4) * monen
       
       lambert_sine_series = lambert_sine_series + term
       
       qn = q*qn
       twonu = twonu + 2._fdp*u
       monen = -monen
       if (debug) counter = counter + 1
       abserr = maxval(abs(term))
    end do

    if (debug) then
       write(*,*)'lambert_sine_series:'
       write(*,*)'counter = ',counter
       write(*,*)'u =   q = ',u,q
       write(*,*)'errors  = ',term(1),term(2),term(3),term(4)
    endif
       
  end function lambert_sine_series

    
!NIST Handbook of Mathematical Functions page 529
  recursive function deriv_elliptic_lnthetas(u,lnq,tol) result(dlnthetas)
    complex(fdp), dimension(ntheta) :: dlnthetas, dlntasthe
    complex(fdp), intent(in) :: u,lnq
    real(fdp), intent(in), optional :: tol

    complex(fdp), parameter :: ipi = cmplx(0._fdp,pidp,fdp)
    complex(fdp), parameter :: pi2 = pidp*pidp
    
    complex(fdp), dimension(4) :: series
    complex(fdp) :: uolnq,tanu,q

    real(fdp), parameter :: toosmall = 0.25_fdp


    if (abs(lnq).le.toosmall) then
       uolnq = u/lnq
       if (abs(aimag(uolnq*ipi)).lt.real(-pidp/lnq,fdp)) then
          dlntasthe = -ipi/lnq * deriv_elliptic_lnthetas(-uolnq*ipi,pi2/lnq,tol) &
               + 2._fdp*uolnq
          dlnthetas(1) = dlntasthe(1)
          dlnthetas(2) = dlntasthe(4)
          dlnthetas(3) = dlntasthe(3)
          dlnthetas(4) = dlntasthe(2)
          return
       endif
    endif
    
    
    tanu = tan(u)
    series = lambert_sine_series(u,exp(lnq),tol)
    
    dlnthetas(1) = 4._fdp*series(1) + 1._fdp/tanu
    dlnthetas(2) = 4._fdp*series(2) - tanu
    dlnthetas(3) = 4._fdp*series(3)
    dlnthetas(4) = 4._fdp*series(4)
        
  end function deriv_elliptic_lnthetas


    
  function deriv_elliptic_thetas(u,lnq,tol)    
    implicit none
    complex(fdp), dimension(ntheta) :: deriv_elliptic_thetas
    complex(fdp), intent(in) :: u,lnq
    real(fdp), intent(in), optional :: tol    

    deriv_elliptic_thetas = elliptic_thetas(u,lnq) * deriv_elliptic_lnthetas(u,lnq,tol)
    
  end function deriv_elliptic_thetas

  


!direct calculation of thetas by fourier series  
  recursive function elliptic_thetas_fourier_series(u,lnq,tol,lnw) result(thetas)
    implicit none
    complex(fdp), intent(in) :: u,lnq
    real(fdp), intent(in), optional :: tol
    complex(fdp), intent(in), optional :: lnw
    
    complex(fdp), dimension(ntheta) :: thetas, tasthe
    
    real(fdp) :: twon, twonm1, n2, nmhalf2, monen
    real(fdp) :: abserr, maxerr
    
    complex(fdp) :: uolnq, sqrpiolnq
    complex(fdp) :: i2nm1u,i2nu
    complex(fdp) :: exp2nm1plus,exp2nm1minus,exp2nplus,exp2nminus

    complex(fdp), dimension(4) :: term

    real(fdp), parameter :: pi2 = pidp*pidp
    complex(fdp), parameter :: i = cmplx(0._fdp,1._fdp)
    complex(fdp), parameter :: ipi = cmplx(0._fdp,pidp)
    
    real(fdp), parameter :: logeps = log(epsilon(0._fdp))
    real(fdp), parameter :: toosmall = 0.25_fdp

    complex(fdp) :: lnnorm
    
    logical, parameter :: debug = .false.
    
    integer :: n

    if (present(lnw)) then
       lnnorm = lnw
    else
       lnnorm = cmplx(0._fdp,0._fdp,fdp)
    endif
    

    if (abs(lnq).le.toosmall) then
       uolnq = u/lnq
       sqrpiolnq = sqrt(-pidp/lnq)

       lnnorm = u*uolnq
!the NaN way
!       tasthe = exp(u*uolnq) * elliptic_thetas_fourier_series(-uolnq*ipi,pi2/lnq,tol)
       tasthe = elliptic_thetas_fourier_series(-uolnq*ipi,pi2/lnq,tol,lnnorm)
       
       thetas(1) = -i*sqrpiolnq*tasthe(1)
       thetas(2) = sqrpiolnq*tasthe(4)
       thetas(3) = sqrpiolnq*tasthe(3)
       thetas(4) = sqrpiolnq*tasthe(2)       

       return
    endif
            
    if (present(tol)) then
       maxerr = tol
    else
       maxerr = epsilon(1._fdp)
    endif

    thetas = cmplx(0._fdp,0._fdp,fdp)
    abserr = 1._fdp
    
    n = 0

    do while (abserr.gt.maxerr)
       n = n+1
       twon = n + n
       twonm1 = twon - 1
       n2 = n*n
       nmhalf2 = n2 + 0.25_fdp - n
       monen = (-1)**n

       i2nm1u = twonm1*u*i
       i2nu = twon*u*i
       exp2nm1plus = exp(nmhalf2*lnq + i2nm1u + lnnorm)
       exp2nm1minus = exp(nmhalf2*lnq - i2nm1u + lnnorm)
       exp2nplus = exp(n2*lnq + i2nu + lnnorm)
       exp2nminus = exp(n2*lnq - i2nu + lnnorm)

       term(1) = monen * (exp2nm1plus - exp2nm1minus)
       term(2) = exp2nm1plus + exp2nm1minus
       term(3) = exp2nplus + exp2nminus
       term(4) = monen * (exp2nplus + exp2nminus)

       thetas = thetas + term

       abserr = maxval(abs(term))
       
    end do
    
    thetas(1) = -thetas(1)/i
    thetas(3) = exp(lnnorm) + thetas(3)
    thetas(4) = exp(lnnorm) + thetas(4)
       
    if (debug) then
       write(*,*)'elliptic_thetas_fourier_series'
       write(*,*)'u = lnq = ',u,lnq
       write(*,*)'errors  = ',term
       write(*,*)'thetas  = ',thetas
    endif
    
  end function elliptic_thetas_fourier_series

!direct calculation of d(thetas)/du by fourier series  
  recursive function deriv_elliptic_thetas_fourier_series(u,lnq,tol,lnw) result(dthetas)
    implicit none
    complex(fdp), intent(in) :: u,lnq
    real(fdp), intent(in), optional :: tol
    complex(fdp), intent(in), optional :: lnw
    
    complex(fdp), dimension(ntheta) :: dthetas
    
    real(fdp) :: twon, twonm1, n2, nmhalf2, monen
    real(fdp) :: abserr, maxerr
    
    complex(fdp) :: uolnq, twouoipi
    complex(fdp) :: i2nm1u,i2nu
    complex(fdp) :: exp2nm1plus,exp2nm1minus,exp2nplus,exp2nminus
    
    complex(fdp), dimension(4) :: dterm

    real(fdp), parameter :: pi2 = pidp*pidp
    complex(fdp), parameter :: i = cmplx(0._fdp,1._fdp)
    complex(fdp), parameter :: ipi = cmplx(0._fdp,pidp)
    
    real(fdp), parameter :: logeps = log(epsilon(0._fdp))
    real(fdp), parameter :: toosmall = 0.25_fdp
    complex(fdp), dimension(ntheta) :: thetas, dtasthe

    complex(fdp) :: lnnorm
    
    logical, parameter :: debug = .false.
    
    integer :: n

    if (present(lnw)) then
       lnnorm = lnw
    else
       lnnorm = cmplx(0._fdp,0._fdp,fdp)
    endif
    

    if (abs(lnq).le.toosmall) then
       uolnq = u/lnq
       twouoipi = 2._fdp*u/ipi

!the Nan way       
!       piolnqthreehalfexpu2olnq = sqrt(-pidp/lnq)*(-pidp/lnq) * exp(u*uolnq)
!       thetas = elliptic_thetas_fourier_series(-uolnq*ipi,pi2/lnq,tol)
!       dtasthe = deriv_elliptic_thetas_fourier_series(-uolnq*ipi,pi2/lnq,tol)
!       dthetas(1) = piolnqthreehalfexpu2olnq * ( dtasthe(1) - thetas(1)*twouoipi )
!       dthetas(2) = piolnqthreehalfexpu2olnq * ( dtasthe(4) - thetas(4)*twouoipi ) * i
!       dthetas(3) = piolnqthreehalfexpu2olnq * ( dtasthe(3) - thetas(3)*twouoipi ) * i
!       dthetas(4) = piolnqthreehalfexpu2olnq * ( dtasthe(2) - thetas(2)*twouoipi ) * i
       
       lnnorm = u*uolnq + 1.5_fdp*log(-pidp/lnq)
       thetas = elliptic_thetas_fourier_series(-uolnq*ipi,pi2/lnq,tol,lnnorm)
       dtasthe = deriv_elliptic_thetas_fourier_series(-uolnq*ipi,pi2/lnq,tol,lnnorm)

       dthetas(1) = ( dtasthe(1) - thetas(1)*twouoipi )
       dthetas(2) = ( dtasthe(4) - thetas(4)*twouoipi ) * i
       dthetas(3) = ( dtasthe(3) - thetas(3)*twouoipi ) * i
       dthetas(4) = ( dtasthe(2) - thetas(2)*twouoipi ) * i
       
       return
    endif

    if (present(tol)) then
       maxerr = tol
    else
       maxerr = epsilon(1._fdp)
    endif
    

    dthetas = cmplx(0._fdp,0._fdp,fdp)
    abserr = 1._fdp
    
    n = 0
    
    do while (abserr.gt.maxerr)
       n = n + 1
       twon = n + n
       twonm1 = twon - 1
       n2 = n*n
       nmhalf2 = n2 + 0.25_fdp - n
       monen = (-1)**n

       i2nm1u = twonm1*u*i
       i2nu = twon*u*i
       exp2nm1plus = exp(nmhalf2*lnq + i2nm1u + lnnorm)
       exp2nm1minus = exp(nmhalf2*lnq - i2nm1u + lnnorm)
       exp2nplus = exp(n2*lnq + i2nu + lnnorm)
       exp2nminus = exp(n2*lnq - i2nu + lnnorm)
       
       dterm(1) = monen * twonm1 * (exp2nm1plus + exp2nm1minus)
       dterm(2) = twonm1 * (exp2nm1plus - exp2nm1minus)
       dterm(3) = twon * (exp2nplus - exp2nminus)
       dterm(4) = monen * twon * (exp2nplus - exp2nminus)

       dthetas = dthetas + dterm

       abserr = maxval(abs(dterm))
       
    end do

    
    dthetas(1) = -dthetas(1)
    dthetas(2) = -dthetas(2)/i
    dthetas(3) = -dthetas(3)/i
    dthetas(4) = -dthetas(4)/i
    
    if (debug) then
       write(*,*)'deriv_elliptic_thetas_fourier_series'
       write(*,*)'u = lnq = ',u,lnq
       write(*,*)'errors  = ',dterm
       write(*,*)'dthetas = ',dthetas
    endif
    
  end function deriv_elliptic_thetas_fourier_series


  
      
end module fflint
