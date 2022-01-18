  interface

     function cvoderhs(n,t,y)
       use precision, only : fdp
       implicit none
       integer, intent(in) :: n
       real(fdp), intent(in) :: t
       real(fdp), dimension(n), intent(in) :: y
       real(fdp), dimension(n) :: cvoderhs
     end function cvoderhs


     function cvodejac(n,t,y,f)
       use precision, only : fdp
       implicit none
       integer, intent(in) :: n
       real(fdp), intent(in) :: t
       real(fdp), dimension(n), intent(in) :: y,f
       real(fdp), dimension(n*n) :: cvodejac
     end function cvodejac

  end interface
