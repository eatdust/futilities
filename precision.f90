module precision
  implicit none
  
  public
  integer, parameter :: fsp=kind(1.0_4)
  integer, parameter :: fdp=kind(1.0_8)
  integer, parameter :: fqp=kind(1.0_16)
  
  integer, parameter :: isp=4
  integer, parameter :: idp=8

  real(fsp), parameter :: pisp = 3.141592653589793238462643383279502884197169399375105820974944592_fsp
  real(fdp), parameter :: pidp = 3.141592653589793238462643383279502884197169399375105820974944592_fdp

end module precision
