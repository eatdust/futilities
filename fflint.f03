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

  public :: ntheta, exp_integral_ei, exp_integral_en
  public :: elliptic_thetas, deriv_elliptic_thetas
  public :: lambert_sine_series, deriv_elliptic_lnthetas
  
contains


  
!fortran to c string converter
  function f_c_string(fname)
    implicit none
    character(len=*), intent(in) :: fname
    character(kind=C_CHAR, len=len(fname)+1) :: f_c_string
    integer :: i,n

    n = len(fname)
    do i=1,n
       f_c_string(i:i) = fname(i:i)
    enddo
    f_c_string(n+1:n+1)=C_NULL_CHAR

  end function f_c_string

  
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
       write(*,*)'errors  = ',term(1),term(2),term(3),term(4)
    endif
       
  end function lambert_sine_series

    
!NIST Handbook of Mathematical Functions page 529  
  function deriv_elliptic_lnthetas(u,q,tol)
    complex(fdp), dimension(ntheta) :: deriv_elliptic_lnthetas
    complex(fdp), intent(in) :: u,q
    real(fdp), intent(in), optional :: tol
    
    complex(fdp), dimension(4) :: series
    complex(fdp) :: tanu

    tanu = tan(u)
    series = lambert_sine_series(u,q,tol)

    deriv_elliptic_lnthetas(1) = 4._fdp*series(1) + 1._fdp/tanu
    deriv_elliptic_lnthetas(2) = 4._fdp*series(2) - tanu
    deriv_elliptic_lnthetas(3) = 4._fdp*series(3)
    deriv_elliptic_lnthetas(4) = 4._fdp*series(4)
        
  end function deriv_elliptic_lnthetas

  

  function deriv_elliptic_thetas(u,lnq,tol)    
    implicit none
    complex(fdp), dimension(ntheta) :: deriv_elliptic_thetas
    complex(fdp), intent(in) :: u,lnq
    real(fdp), intent(in), optional :: tol
    
    complex(fdp), dimension(ntheta) :: dlnthetas, thetas

    deriv_elliptic_thetas = elliptic_thetas(u,lnq) * deriv_elliptic_lnthetas(u,exp(lnq),tol)
    
  end function deriv_elliptic_thetas

      
end module fflint
