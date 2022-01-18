module scheduler
#ifdef MPISCHED
  use mpi
#endif
  use fifo, only : ends,qval
  use fifo, only : nullify_queue_ends, free_queue
  use fifo, only : push_in_back, pop_out_back
  use fifo, only : push_in_front, pop_out_front
  use fifo, only : is_queue_empty, queue_size
  
  implicit none


  private

  
  integer, parameter :: PollingTime = 1

  integer, parameter :: TagNode = 11
  integer, parameter :: CountOne = 1
  integer, parameter :: NullAssert = 0
  integer, parameter :: QLockFlag = -1
  integer, parameter :: QStarvFlag = 0
  integer, parameter :: QInitFlag = 10


#ifdef MPISCHED
  integer :: WinOnQ
  integer :: IntSize, DispSize
  integer(MPI_ADDRESS_KIND) :: QrdmaAddress
  integer(MPI_ADDRESS_KIND) :: WinSize
  integer(MPI_ADDRESS_KIND), parameter :: ZeroDisplace = 0
#else
  integer :: QrdmaAddress
#endif
 
  integer, parameter :: MinDebug = 1
  integer, parameter :: MediumDebug = 2
  integer, parameter :: MaxDebug = 3
  integer, parameter :: debugLevel = MinDebug


  type(ends) :: qbounds  

  logical, parameter :: ReDistribOnRestart = .true.

    
  public initialize_scheduler, free_scheduler, scheduler_save_queue
  public start_scheduling, irq_scheduler, stop_scheduling
  public restore_scheduler, scheduled_size
      
contains


  function scheduled_size()
    implicit none
    integer :: scheduled_size

    scheduled_size = qbounds%size

  end function scheduled_size


  subroutine initialize_scheduler(ntodo)    
    implicit none
    integer, intent(in) :: ntodo
    integer :: rank
       
    call create_new_queue(ntodo)

    call initialize_rdma()
    
  end subroutine initialize_scheduler



  subroutine restore_scheduler(nq)
    implicit none
    integer, intent(in) :: nq

    if (ReDistribOnRestart) then
       call redistribute_saved_queues(nq)
    else
       call restore_saved_queues()
    endif

    call initialize_rdma()

  end subroutine restore_scheduler




  subroutine initialize_rdma()
    implicit none
    

    logical :: flag
    integer :: rank, code
#ifdef MPISCHED
    integer(MPI_ADDRESS_KIND) :: taille, base, dispunit
#endif

    rank = get_mpi_rank()

#ifdef MPISCHED

    call MPI_TYPE_SIZE(MPI_INTEGER,IntSize,code)
         
    WinSize = IntSize !one element
    DispSize = IntSize


    call MPI_WIN_ALLOCATE(WinSize,DispSize,MPI_INFO_NULL &
         ,MPI_COMM_WORLD,QrdmaAddress,WinOnQ,code)

    call MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,rank,NullAssert,WinOnQ,code)

    call MPI_PUT(QInitFlag,CountOne,MPI_INTEGER,rank,ZeroDisplace &
         ,CountOne,MPI_INTEGER,WinOnQ,code)

    call MPI_WIN_UNLOCK(rank,WinOnQ,code)

    if (debugLevel.ge.minDebug) then
       write(*,*)
       write(*,*)'initialize_scheduler: Opening WIN on rank:',rank
   
       call MPI_WIN_GET_ATTR(winOnQ,MPI_WIN_SIZE,taille,flag,code)
       call MPI_WIN_GET_ATTR(winOnQ,MPI_WIN_DISP_UNIT,dispunit,flag,code)
       call MPI_WIN_GET_ATTR(winOnQ,MPI_WIN_BASE,base,flag,code)

       write(*,*)'window opened with:',taille,dispunit,base
       write(*,*)
    endif

#endif

  end subroutine initialize_rdma


 


  subroutine create_new_queue(ntodo)
    implicit none
    integer, intent(in) :: ntodo    
    type(qval) :: qbuffer
    integer :: nwork, rank, qsize
    integer :: i,j

    nwork = get_mpi_size()
    rank = get_mpi_rank()
    
    if (mod(ntodo,nwork).ne.0) then
       write(*,*)'ntodo= nwork= ',ntodo,nwork,mod(ntodo,nwork)
       write(*,*)'create_new_queue: non-equal queue sizes!'
    endif

    
!!!create a qsize queue in each of the mpi processeses
!!!clear the boundaries
!!    qsize = ntodo/nwork
!!    call nullify_queue_ends(qbounds)

!interleaves 0,1,2..nqsize-1 into the queue
!!    do i = 0,qsize-1
!!       qbuffer%id = rank + i*nwork
!!       call push_in_back(qbounds, qbuffer)
!!    end do

!allow for non-equals queue size
    call nullify_queue_ends(qbounds)
    
    j=rank
    qsize = 0
    do while (j.le.ntodo-1)
       qbuffer%id = j
       call push_in_back(qbounds, qbuffer)
       qsize = qsize + 1
       j = j + nwork
    enddo

!for debugging only
!    qbuffer%id=125
!    call push_in_back(qbounds, qbuffer)
!    qbuffer%id=466
!    call push_in_back(qbounds, qbuffer)
!    qbuffer%id=807
!    call push_in_back(qbounds, qbuffer)
!    qbuffer%id=1148
!    call push_in_back(qbounds, qbuffer)
!    call scheduler_save_queue(rank)

  end subroutine create_new_queue





  subroutine restore_saved_queues()
    use iofifo, only : read_queue_vals
    implicit none    
    type(qval), dimension(:), pointer :: ptrqvals => null()

    integer :: rank, nwork
    integer :: i, qsize

    rank = get_mpi_rank()
    nwork = get_mpi_size()

    call nullify_queue_ends(qbounds)

    if (nwork.gt.1) then
       call read_queue_vals(ptrqvals, rank)
    else
       call read_queue_vals(ptrqvals)
    endif

    qsize = size(ptrqvals,1)
    
    do i=1,qsize
       call push_in_back(qbounds,ptrqvals(i))
    enddo

    deallocate(ptrqvals)
    ptrqvals => null()

  end subroutine restore_saved_queues





  subroutine redistribute_saved_queues(nq)
    use fifo, only : node
    use iofifo, only : read_queue_vals
    implicit none
    integer, intent(in) :: nq

    type(qval), dimension(:), pointer :: ptrqvals => null()
    type(qval), dimension(:), allocatable :: valbuffer

    type(node), pointer :: ptrnode
    type(ends) :: qtmp

    integer :: rank, nwork, ntodo
    integer :: i,j,qsize

    rank = get_mpi_rank()
    nwork = get_mpi_size()


!aggregate the whole old queue in qtmp

    call nullify_queue_ends(qtmp)

    do i=0,nq-1

       if (nq.gt.1)  then
          call read_queue_vals(ptrqvals,i)
       else
          call read_queue_vals(ptrqvals)
       endif

       qsize = size(ptrqvals,1)
       
       do j=1,qsize
          call push_in_back(qtmp,ptrqvals(j))
       enddo
        
       deallocate(ptrqvals)
       ptrqvals => null()

    enddo

    ntodo = qtmp%size


!using a temporary array
    allocate(valbuffer(0:ntodo-1))

    i=0
    ptrnode => qtmp%front  
    do 
       valbuffer(i) = ptrnode%val
       if (.not.associated(ptrnode%next)) exit
       i=i+1
       ptrnode => ptrnode%next
    enddo
    ptrnode => null()
    if (i.ne.ntodo-1) stop 'redistribute_saved_queue: screwed!'

    call free_queue(qtmp)


    
    if (mod(ntodo,nwork).ne.0) then
       write(*,*)'ntodo= nwork= ',ntodo,nwork,mod(ntodo,nwork)
       write(*,*)'redistribute_saved_queues: queue of unequal sizes!'
    endif

!redistribute withing the new processes           
       
    call nullify_queue_ends(qbounds)
    
!    do i = 0,qsize-1
!       j = rank + i*nwork       
!       call push_in_back(qbounds, valbuffer(j))       
!    end do

    j=rank
    qsize = 0
    do while (j.le.ntodo-1)
       call push_in_back(qbounds, valbuffer(j))
       qsize = qsize + 1
       j = j + nwork
    enddo

    print *,'qsize',qsize,qbounds%size

    deallocate(valbuffer)

  end subroutine redistribute_saved_queues




  subroutine scheduler_save_queue(rank)
    use iofifo, only : dump_queue_vals
    implicit none

    integer :: rank, nsize
        
    nsize = get_mpi_size()

    if (is_queue_empty(qbounds)) return

    if (nsize.gt.1) then
       call dump_queue_vals(qbounds, rank)
    else
       call dump_queue_vals(qbounds)
    endif

  end subroutine scheduler_save_queue



  subroutine start_scheduling(outdex)
    implicit none
    integer, intent(out) :: outdex
    type(qval) :: frontval

    integer :: rank,code    

    call pop_out_front(qbounds,frontval)

    outdex = frontval%id
    rank = get_mpi_rank()


    if (get_mpi_size().eq.1) return

    if (debugLevel.ge.MaxDebug) then
       write(*,*)
       write(*,*)'start_scheduling: local store at RDMA'
       write(*,*)'rank= QrdmaAddress= ',rank,QrdmaAddress
       write(*,*)
    endif
          
  end subroutine start_scheduling




  subroutine irq_scheduler()
    implicit none

    integer, save :: chunk = 1
    integer :: nsize,rank
    
!initialization
    rank = get_mpi_rank()
    nsize = get_mpi_size()

    
    if (debugLevel.ge.MediumDebug) then
       write(*,*)'------------------------------------------------'
       write(*,*)'irq_scheduler:    rank= qsize= chunk=     ',rank,qbounds%size,chunk
       write(*,*)
    endif
    
    
    if (nsize.eq.1) return


            
    if (qbounds%size.gt.2*chunk) then

       call give_nodes(chunk)
     
    endif


    if (is_queue_empty(qbounds)) then
       
       call steal_nodes(chunk)
       
    endif
             
  end subroutine irq_scheduler


  

  function stop_scheduling()
    implicit none
    logical :: stop_scheduling

    integer :: rank

    stop_scheduling = .false.

    if (is_queue_empty(qbounds)) stop_scheduling=.true.


    rank = get_mpi_rank()


    if (debugLevel.ge.MediumDebug) then       
       write(*,*)'stop_scheduling:  rank= qsize= stop =     ',rank &
            ,qbounds%size,stop_scheduling
    endif


  end function stop_scheduling






  subroutine free_scheduler()
    implicit none
    integer :: code

#ifdef MPISCHED    
    call MPI_WIN_FREE(WinOnQ,code)
#endif


  end subroutine free_scheduler






  subroutine give_nodes(chunk)
    implicit none
    integer, intent(in) :: chunk
    type(qval) :: backval

    integer :: rank,code
    integer, dimension(:), allocatable, volatile :: msg

    integer :: targrank,nsize
    integer :: i,GetQrdma,ibuffer

    logical :: flag

#ifdef MPISCHED
    integer(MPI_ADDRESS_KIND) :: taille, base,shift
#endif
    
#ifdef MPISCHED

    rank = get_mpi_rank()
    nsize = get_mpi_size()

    if (debugLevel.ge.MediumDebug) then
       write(*,*)'give_nodes: entering       || rank= qsize=',rank,qbounds%size
       write(*,*)
    endif

    
    allocate(msg(chunk))

    do targrank = 0,nsize-1

       if (rank.eq.targrank) cycle


       if (debugLevel.ge.MaxDebug) then
          write(*,*)'give_nodes: lock/get/put   || orig= targ= ',rank,targrank
       endif


       call MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,targrank,NullAssert,WinOnQ,code)

       
!       call MPI_GET(GetQrdma,CountOne,MPI_INTEGER,targrank,ZeroDisplace &
!            ,CountOne,MPI_INTEGER,WinOnQ,code)
!       call MPI_PUT(QLockFlag,CountOne,MPI_INTEGER,targrank,ZeroDisplace &
!            ,CountOne,MPI_INTEGER,WinOnQ,code)

       call MPI_FETCH_AND_OP(QLockFlag,GetQrdma,MPI_INTEGER,targrank,ZeroDisplace &
          ,MPI_REPLACE,WinOnQ,code)


       call MPI_WIN_UNLOCK(targrank,WinOnQ,code)
      

       if (debugLevel.ge.MaxDebug) then
          write(*,*)'give_nodes: unlocked       || orig= targ= ',rank,targrank
          write(*,*)'give_nodes: orig= targ= GetQrdma=          ',rank,targrank,GetQrdma
       endif

!TEST
!       call MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,targrank,NullAssert,WinOnQ,code)
!       call MPI_GET(ibuffer,CountOne,MPI_INTEGER,targrank,ZeroDisplace &
!            ,CountOne,MPI_INTEGER,WinOnQ,code)
!       call MPI_WIN_UNLOCK(targrank,WinOnQ,code)
!       print *,'sjhould be -1',ibuffer
!ENDTEST


       select case (GetQrdma)

       case (QLockFlag)

          cycle

       case (QStarvFlag)

          do i=1,chunk       
             call pop_out_back(qbounds,backval)
             msg(i) = backval%id
          enddo

          if (debugLevel.ge.MinDebug) then
             write(*,*)
             write(*,*)'give_nodes:  SENDING INIT  || orig= targ= msg= ',rank,targrank,msg
          endif
          
          call MPI_SSEND(msg,chunk,MPI_INTEGER,targrank,TagNode,MPI_COMM_WORLD,code)

          if (debugLevel.ge.MinDebug) then
             write(*,*)'give_nodes:  SENDING DONE  || orig= targ= ',rank,targrank
             write(*,*)
          endif

          exit
          
       end select

    enddo

    deallocate(msg)

#endif
    
  end subroutine give_nodes





   
  subroutine steal_nodes(chunk)
    implicit none
    integer, intent(in) :: chunk

    type(qval) :: backval
 
    logical, volatile :: WorkAvailable
    integer, volatile :: request
    integer, dimension(:), allocatable, volatile :: msg

    integer :: code, rank, fromrank, nsize
    logical :: RunCompleted

    integer :: i

#ifdef MPISCHED
    integer, dimension(MPI_STATUS_SIZE) :: status    
#endif


    nsize = get_mpi_size()
    rank = get_mpi_rank()
    
   
#ifdef MPISCHED

    call MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,rank,NullAssert,WinOnQ,code)

    call MPI_PUT(QStarvFlag,CountOne,MPI_INTEGER,rank,ZeroDisplace &
         ,CountOne,MPI_INTEGER,WinOnQ,code)

    call MPI_WIN_UNLOCK(rank,WinOnQ,code)


    WorkAvailable = .false.
    RunCompleted = .false.    

    allocate(msg(chunk))

    if (debugLevel.ge.MediumDebug) then
       write(*,*)'steal_nodes: WAITING MSG   || rank= Qrdma= ',rank,QrdmaAddress
       write(*,*)
    endif
       

    call MPI_IRECV(msg,chunk,MPI_INTEGER,MPI_ANY_SOURCE,TagNode &
         ,MPI_COMM_WORLD,request,code)

    do
       
       call MPI_TEST(request,WorkAvailable,status,code)
       if (WorkAvailable) exit

       call check_all_starving(RunCompleted)
       if (RunCompleted) exit

       call sleep(PollingTime)

    end do


    if (WorkAvailable) then
    
       fromrank = status(MPI_SOURCE)

       do i=1,chunk
          backval%id = msg(i)          
          call push_in_back(qbounds,backval)          
       enddo

       if (debugLevel.ge.MinDebug) then
          write(*,*)'steal_nodes: RECEIVED DONE || rank= from= msg= ',rank,fromrank,msg
          write(*,*)
       endif
    
    else
       
       if (.not.RunCompleted) stop 'steal_nodes: polling bugged!'

    endif

    deallocate(msg)



#endif

  end subroutine steal_nodes





  subroutine check_all_starving(flag)
    implicit none
    logical, intent(out) :: flag

    integer :: targrank,rank,nsize,code
    integer :: nstarv

    integer :: GetMyFlag
    integer, volatile :: GetQFlag


    flag = .true.


#ifdef MPISCHED

    rank = get_mpi_rank()
    nsize = get_mpi_size()
    
    nstarv = 0

    do targrank = 0,nsize-1

       if (targrank.eq.rank) cycle

       call MPI_WIN_LOCK(MPI_LOCK_SHARED,targrank,NullAssert,WinOnQ,code)
       
       call MPI_GET(GetQFlag,CountOne,MPI_INTEGER,targrank,ZeroDisplace &
            ,CountOne,MPI_INTEGER,WinOnQ,code)

       call MPI_WIN_UNLOCK(targrank,WinOnQ,code)
          
                 
       if (debugLevel.ge.MaxDebug) then
          write(*,*)'check_all_starving: orig= targ= GetQFlag= ',rank,targrank,GetQFlag
       endif
     

       select case (GetQFlag)

       case (QLockFlag)
          flag = .false.
          exit

       case (QStarvFlag)
          nstarv = nstarv + 1

       case default
          flag = .false.
          exit

       end select
    enddo

!IMPORTANT:
!test myself only when I am sure than all the others are starving;
!otherwise they could still being accessing my RDMA while I am checking
!their

    if (flag) then

       if (nstarv.ne.nsize-1) stop 'check_all_starving: bug!!'

       call MPI_WIN_LOCK(MPI_LOCK_SHARED,rank,NullAssert,WinOnQ,code)
       
       call MPI_GET(GetMyFlag,CountOne,MPI_INTEGER,rank,ZeroDisplace &
            ,CountOne,MPI_INTEGER,WinOnQ,code)
       
       call MPI_WIN_UNLOCK(rank,WinOnQ,code)

       if (GetMyflag.eq.QLockFlag) flag=.false.
    endif

#endif
           
  end subroutine check_all_starving




  function get_mpi_rank()
    implicit none
    integer :: get_mpi_rank
    integer :: code

    get_mpi_rank = 0

#ifdef MPISCHED
    call mpi_comm_rank(MPI_COMM_WORLD,get_mpi_rank,code)
#endif
  end function get_mpi_rank



  function get_mpi_size()
    implicit none
    integer :: get_mpi_size
    integer :: code

    get_mpi_size = 1

#ifdef MPISCHED
    call mpi_comm_size(MPI_COMM_WORLD,get_mpi_size,code)
#endif
  end function get_mpi_size



 
!unused
  subroutine find_longest_queue(longqsize,longrank)
    implicit none
    integer, intent(out) :: longqsize
    integer, intent(out) :: longrank

    integer :: rank,qsize
   
    integer :: code,msize
    integer, dimension(2) :: inpair, outpair


    longrank = 0

#ifdef MPISCHED

    qsize = queue_size(qbounds)
    rank = get_mpi_rank()

    msize = 1

    inpair(1) = qsize
    inpair(2) = rank
    
    call MPI_ALLREDUCE(inpair,outpair,msize,MPI_2INTEGER,MPI_MAXLOC &
         ,MPI_COMM_WORLD,code)
   
    longqsize = outpair(1)
    longrank = outpair(2)
   
#endif    

  end subroutine find_longest_queue




end module scheduler
