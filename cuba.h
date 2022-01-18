  interface
     integer function integrand(ndim,x,ncomp,f,udata,nvec,ncore)
       implicit none
       integer, parameter :: dp = kind(1._8)
       integer :: ndim, ncomp
       integer :: udata, nvec, ncore
       real(dp), dimension(ndim,nvec) :: x
       real(dp), dimension(ncomp,nvec) :: f
     end function integrand
  end interface
