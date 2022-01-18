  interface

     subroutine initu(ncomp, nmsh, xx, nudim,u,rpar,ipar)
       use bvpprec
       implicit none
       integer :: ncomp, nmsh,nudim
       real(dp), dimension(*) :: rpar
       integer, dimension(*) :: ipar
       real(dp), dimension(*) :: xx
       real(dp), dimension(nudim,*) :: u
     end subroutine initu

     subroutine fsub(ncomp,x,u,f,rpar,ipar)
       use bvpprec
       implicit none          
       integer, intent(in) :: ncomp
       real(dp), intent(in) :: x
       real(dp), dimension(*) :: u
       real(dp), dimension(*) :: f
       real(dp), dimension(*) :: rpar
       integer, dimension(*) :: ipar
     end subroutine fsub

     subroutine dfsub(ncomp,x,u,df,rpar,ipar)
       use bvpprec
       implicit none          
       integer, intent(in) :: ncomp
       real(dp), intent(in) :: x
       real(dp), dimension(*) :: u
       real(dp), dimension(ncomp,*) :: df
       real(dp), dimension(*) :: rpar
       integer, dimension(*) :: ipar
     end subroutine dfsub

     subroutine gsub(i,ncomp,u,g,rpar,ipar)
       use bvpprec
       implicit none
       integer, intent(in) :: i,ncomp
       
       real(dp) :: g
       real(dp), dimension(*) :: u
       real(dp), dimension(*) :: rpar
       integer, dimension(*) :: ipar
     end subroutine gsub

     subroutine dgsub(i,ncomp,u,dg,rpar,ipar)
       use bvpprec
       implicit none
       integer, intent(in) :: i,ncomp
       real(dp), dimension(*) :: u
       real(dp), dimension(*) :: dg
       real(dp), dimension(*) :: rpar
       integer, dimension(*) :: ipar
     end subroutine dgsub

  end interface
