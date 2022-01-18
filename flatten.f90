module flatten
  implicit none

contains

  function flatten_indices(ndim,isize,ivec)
    implicit none
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: isize(ndim)
    integer(kind=8), dimension(ndim), intent(in) :: ivec

    integer(kind=8) :: flatten_indices

    integer(kind=8) :: q,j

    q = ivec(1)
    do j=2,ndim
       q = q + product(isize(1:j-1))*(ivec(j)-1)
    enddo

    flatten_indices = q

  end function flatten_indices


  function unflatten_indices(ndim,isize,q)
    implicit none
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), dimension(ndim) :: unflatten_indices
    integer(kind=8), intent(in) :: isize(ndim)
    integer(kind=8), intent(in) :: q

    integer(kind=8), dimension(ndim) :: ivec
    integer(kind=8) :: sum
    integer(kind=8) :: k

    if (ndim.eq.1) then
       ivec(1) = q
       return
    endif

    ivec(ndim) = (q-1)/product(isize(1:ndim-1)) + 1
    sum = 0

    do k=ndim-1,1,-1
       sum = sum + (ivec(k+1)-1) * product(isize(1:k))
       if (k.gt.1) then
          ivec(k) = (q-1 - sum)/product(isize(1:k-1)) + 1
       else
          ivec(k) = (q-1 - sum) + 1
       endif
    enddo

    unflatten_indices = ivec

  end function unflatten_indices

end module flatten
