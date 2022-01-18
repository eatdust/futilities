module fifo
  implicit none

  private

  type qval
     integer :: id
  end type qval

  type node
     type(qval) :: val
     type(node), pointer :: next
     type(node), pointer :: prev
  end type node
 
  type ends
     integer :: size
     type(node), pointer :: front, back
  end type ends
 
  
  public ends, qval, node
  public nullify_queue_ends, free_queue
  public push_in_front, pop_out_front
  public push_in_back, pop_out_back
  public is_queue_empty, queue_size 

contains
 

  subroutine nullify_queue_ends(bounds)
    implicit none
    type(ends), intent(inout) :: bounds
    bounds%front => null()
    bounds%back => null()
    bounds%size = 0
  end subroutine nullify_queue_ends
 


  subroutine push_in_back(bounds,setval)
    implicit none
    type(ends), intent(inout) :: bounds
    type(qval), intent(in) :: setval
    type(node), pointer :: bnode

    allocate(bnode)
    bnode%val = setval
   
    if ( associated(bounds%back) ) then
       bounds%back%next => bnode
       bnode%prev => bounds%back
       bounds%back => bnode
    else
       bounds%back => bnode
       bounds%front => bnode
       bnode%prev => null()
    end if

    bnode%next => null()
    bounds%size = bounds%size + 1

    bnode => null()

  end subroutine push_in_back


  
  subroutine pop_out_back(bounds,getval)
    implicit none
    type(ends), intent(inout) :: bounds
    type(qval), intent(out) :: getval
    
    type(node) :: rbnode

    if (.not.associated(bounds%back)) then
       write(*,*)'pop_out_back: nothing to remove!'
       stop
    endif

    rbnode = bounds%back
    getval = rbnode%val

    deallocate(bounds%back)
    
    if (associated(rbnode%prev)) then
       bounds%back => rbnode%prev
       bounds%back%next => null()
    else
       bounds%front => null()
       bounds%back => null()
    end if
    
    rbnode%next => null()
    rbnode%prev => null()
    bounds%size = bounds%size - 1

  end subroutine pop_out_back




  subroutine push_in_front(bounds,setval)
    implicit none
    type(ends), intent(inout) :: bounds
    type(qval), intent(in) :: setval
    type(node), pointer :: fnode

    allocate(fnode)
    fnode%val = setval

    if ( associated(bounds%front) ) then
       bounds%front%prev => fnode
       fnode%next => bounds%front
       bounds%front => fnode
    else
       bounds%front => fnode
       bounds%back => fnode       
       fnode%next => null()
    end if

    fnode%prev => null()
    bounds%size = bounds%size + 1

    fnode => null()

  end subroutine push_in_front




 
  subroutine pop_out_front(bounds,getval)
    implicit none
    type(ends), intent(inout) :: bounds
    type(qval), intent(out) :: getval

    type(node) :: rfnode
 
    if (.not.associated(bounds%front)) then
       write(*,*)'pop_out_front: nothing to remove!'
       stop
    endif

    rfnode = bounds%front
    getval = rfnode%val

    deallocate(bounds%front)
    
    if (associated(rfnode%next)) then
       bounds%front => rfnode%next
       bounds%front%prev => null()
    else
       bounds%front => null()
       bounds%back => null()
    end if
    
    rfnode%next => null()
    rfnode%prev => null()
    bounds%size = bounds%size - 1

  end subroutine pop_out_front
 



  function is_queue_empty(bounds)
    implicit none
    logical :: is_queue_empty
    type(ends), intent(inout) :: bounds

    logical :: unsane = .false.

    unsane = xor(associated(bounds%front),associated(bounds%back))

    if (unsane) then
       stop 'is_queue_empty: queue screwed!'
    endif

    
    if ( associated(bounds%front) ) then
       is_queue_empty = .false.
    else
       is_queue_empty = .true.
    end if
    
  end function is_queue_empty
 


  function queue_size(qbounds)
    implicit none
    integer :: queue_size
    type(ends), intent(inout) :: qbounds
    
    type(node), pointer :: ptrnode
    integer :: sizecheck

    if (is_queue_empty(qbounds)) then
       queue_size = 0
       return
    endif
    
    ptrnode => qbounds%front

    queue_size = fwcount_queue(ptrnode) &
         + bwcount_queue(ptrnode) - 1

    ptrnode => qbounds%back
    sizecheck = fwcount_queue(ptrnode) &
         + bwcount_queue(ptrnode) - 1

    if (queue_size.ne.sizecheck) then
       write(*,*)'fwcount= bwcount= ',queue_size,sizecheck
       stop 'queue_size: queue is screwed!'
    endif

    ptrnode => null()

  end function queue_size


 
  recursive function fwcount_queue(anynode) result(count)
    implicit none
    integer :: count    
    type(node), intent(in) :: anynode
    
    count = 1
    
    if (associated(anynode%next)) then
       count = fwcount_queue(anynode%next) + 1
    endif

  end function fwcount_queue



  recursive function bwcount_queue(anynode) result(count)
    implicit none
    integer :: count    
    type(node), intent(in) :: anynode
  
    count = 1
    
    if (associated(anynode%prev)) then
       count = bwcount_queue(anynode%prev) + 1
    endif

  end function bwcount_queue



  
  subroutine fwprint_queue(qnode)
    implicit none
     type(node), intent(in), target :: qnode
     type(node), pointer :: ptrnode

     ptrnode => qnode

     do
        write(*,*)'fwprint_queue: ',ptrnode%val
        if (.not.associated(ptrnode%next)) exit
        ptrnode => ptrnode%next
     end do

     ptrnode => null()

  end subroutine fwprint_queue



  subroutine bwprint_queue(qnode)
    implicit none
    type(node), intent(in), target :: qnode
    type(node), pointer :: ptrnode

    ptrnode => qnode

    do
       write(*,*)'bwprint_queue: ',ptrnode%val
       if (.not.associated(ptrnode%prev)) exit
        ptrnode => ptrnode%prev
     end do

     ptrnode => null()

   end subroutine bwprint_queue




   subroutine free_queue(qbounds)
     implicit none
     type(ends), intent(inout) :: qbounds

     logical :: goforward=.false.

     if (is_queue_empty(qbounds)) then
        write(*,*)'free_queue: nothing to free!'
        return
     endif
    
     if (goforward) then
        if (associated(qbounds%front)) call fwfree_queue(qbounds%front)
     else
        if (associated(qbounds%back)) call bwfree_queue(qbounds%back)
     endif         
    
     call nullify_queue_ends(qbounds)
     
   end subroutine free_queue



   recursive subroutine fwfree_queue(ptrnode)
     implicit none
     type(node), pointer, intent(inout) :: ptrnode
    
     if (.not.associated(ptrnode)) return

     if (associated(ptrnode%next)) call fwfree_queue(ptrnode%next)
     ptrnode%prev => null()
     ptrnode%next => null()
     deallocate(ptrnode)

     
   end subroutine fwfree_queue



   recursive subroutine bwfree_queue(ptrnode)
     implicit none
     type(node), pointer, intent(inout) :: ptrnode

     if (.not.associated(ptrnode)) return

     if (associated(ptrnode%prev)) call bwfree_queue(ptrnode%prev)
     ptrnode%next => null()
     ptrnode%prev => null()
     deallocate(ptrnode)
     
   end subroutine bwfree_queue


end module fifo
