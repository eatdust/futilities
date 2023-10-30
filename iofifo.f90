!   This file is part of futilities
!
!   Copyright (C) 2013-2021 C. Ringeval
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

module iofifo
  use fifo, only : ends, qval, node
  implicit none

  private

  public dump_queue_vals, read_queue_vals

  public write_queue_vals, check_queue_file

  integer, parameter :: lenIoRank = 4
  integer, parameter :: lenIoMax = 64
  
  logical, parameter :: display = .true.


contains

  subroutine int2char_zero(icount,clen,ccount)
    implicit none
    integer, intent(in) :: icount
    integer, intent(in) :: clen
    character(len=clen), intent(out) :: ccount

    character(len=lenIoMax) :: numToStrg
    character(len=clen) :: strg

    if (clen.gt.lenIoMax) then
       stop 'step2char: clen > lenIoMax!'
    endif

    write(numToStrg,*) icount
    strg = trim(adjustl(numToStrg)) 
    strg = adjustr(strg)
    call replace_char(strg,' ','0')

    ccount(1:clen) = strg(1:clen)

  end subroutine int2char_zero



  subroutine replace_char(strg,charold,charnew)
    implicit none
    character(len=*), intent(inout) :: strg
    character, intent(in) :: charold, charnew

    integer :: position
    position = 1

    do 
       position = scan(strg,charold)
       if (position.ne.0) then
          strg(position:position) = charnew
       else
          exit
       endif
    enddo
        
  end subroutine replace_char


  
  function queue_filename(rank)    
    implicit none
    character(len=lenIoMax) :: queue_filename
    integer, intent(in), optional :: rank

    character(len=lenIoRank) :: cRank


    if (present(rank)) then
       call int2char_zero(rank,lenIoRank,cRank)
       queue_filename = 'queue'//cRank//'.bin'
    else
       queue_filename = 'queue.bin'
    endif

  end function queue_filename



  subroutine dump_queue_vals(qbounds, rank)
    implicit none
    type(ends), intent(in) :: qbounds
    integer, intent(in), optional :: rank
    
    type(node), pointer :: ptrnode
    character(len=lenIoMax) :: filename

    integer, parameter :: nunit = 300
    logical, parameter :: goforward = .true.


    if (qbounds%size.eq.0) stop 'write_queue_vals: queue empty!'

    if (present(rank)) then
       filename = queue_filename(rank)
    else
       filename = queue_filename()
    endif
   
  
    open(file=filename,unit=nunit,form='unformatted',status='unknown')

    write(nunit) qbounds%size

    if (goforward) then

       ptrnode => qbounds%front

       do
          write(nunit) ptrnode%val
          if (.not.associated(ptrnode%next)) exit
          ptrnode => ptrnode%next
       end do

    else

       ptrnode => qbounds%back

       do
          write(nunit)ptrnode%val
          if (.not.associated(ptrnode%prev)) exit
          ptrnode => ptrnode%prev
       end do

    endif

    close(nunit)

    ptrnode => null()

  end subroutine dump_queue_vals


  function check_queue_file(rank)
    implicit none
    integer, intent(in), optional :: rank
    logical :: check_queue_file

    integer, parameter :: nunit = 301
    character(len=lenIoMax) :: filename

    logical :: ishere
    integer :: nsize,i

    if (present(rank)) then
       filename = queue_filename(rank)
    else
       filename = queue_filename()
    endif

    inquire(file=filename,exist=ishere) 
    
    check_queue_file = ishere
        
  end function check_queue_file


  
  subroutine read_queue_vals(ptrval, rank)
    implicit none
    type(qval), dimension(:), intent(inout), pointer :: ptrval
    integer, intent(in), optional :: rank

    integer, parameter :: nunit = 301
    character(len=lenIoMax) :: filename

    logical :: ishere
    integer :: nsize,i

    if (associated(ptrval)) stop 'read_queue: ptrval already full!'

    if (present(rank)) then
       filename = queue_filename(rank)
    else
       filename = queue_filename()
    endif

    inquire(file=filename,exist=ishere) 
    if (.not.ishere) then
       write(*,*)'searching file ',filename
       stop 'read_queue_vals: file not found!'
    endif

    open(file=filename,unit=nunit,form='unformatted',status='old')
    read(nunit) nsize

    if (display) then
       write(*,*)'read_queue_vals: reading queue... ',nsize
    endif

    allocate(ptrval(nsize))

    do i=1,nsize
       read(nunit) ptrval(i)
    enddo

    close(nunit)

  end subroutine read_queue_vals



  subroutine write_queue_vals(qbounds, rank)
    implicit none
    type(ends), intent(in) :: qbounds
    integer, intent(in), optional :: rank
    
    type(node), pointer :: ptrnode
    character(len=lenIoMax) :: filename

    integer, parameter :: nunit = 302
    logical, parameter :: goforward = .true.

    if (qbounds%size.eq.0) stop 'write_queue_vals: queue empty!'
    
    if (present(rank)) then
       filename = 'asc'//queue_filename(rank)
    else
       filename = 'asc'//queue_filename()
    endif

    print *,'filename ', filename, qbounds%size
  

    open(file=filename,unit=nunit,form='formatted',status='new')

    write(nunit,*) qbounds%size

    if (goforward) then
       print *,'val',associated(qbounds%front)
       ptrnode => qbounds%front


       do
          print *,'val',ptrnode%val
          write(nunit,*) ptrnode%val
          if (.not.associated(ptrnode%next)) exit
          ptrnode => ptrnode%next
       end do

    else

       ptrnode => qbounds%back

       do
          write(nunit,*) ptrnode%val
          if (.not.associated(ptrnode%prev)) exit
          ptrnode => ptrnode%prev
       end do

    endif

    close(nunit)

    ptrnode => null()

  end subroutine write_queue_vals




end module iofifo
