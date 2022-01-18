module iotools

  implicit none
  
  integer, parameter :: lenIoRank = 4
  integer, parameter :: lenIoMax = 64
  integer, parameter :: sp=kind(1._4)
  integer, parameter :: dp=kind(1._8)
  integer, parameter :: qp=kind(1._16)
  
 
  interface livewrite
     module procedure sp_livewrite, dp_livewrite, qp_livewrite
  end interface

  interface allwrite
     module procedure sp_allwrite, dp_allwrite, qp_allwrite
  end interface

  interface binallwrite
     module procedure sp_binallwrite, dp_binallwrite
  end interface binallwrite

  integer, parameter :: reclUnit = 4
  
  private

  public lenIoRank, lenIoMax
  public replace_char,int2char_zero,int2char_adjtrim
  public count_column
  public delete_file,livewrite,allwrite,binallwrite

contains

  

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


  

  subroutine int2char_zero(icount,clen,ccount)
    implicit none
    integer, intent(in) :: icount
    integer, intent(in) :: clen
    character(len=clen), intent(out) :: ccount

    character(len=lenIoMax) :: numToStrg
    character(len=clen) :: strg

    if (clen.gt.lenIoMax) then
       stop 'int2char_zero: clen > lenIoMax!'
    endif

    write(numToStrg,*) icount
    strg = trim(adjustl(numToStrg)) 
    strg = adjustr(strg)
    call replace_char(strg,' ','0')

    ccount(1:clen) = strg(1:clen)

  end subroutine int2char_zero




  subroutine int2char_adjtrim(icount,clen,ccount)
    implicit none
    integer, intent(in) :: icount
    integer, intent(in) :: clen
    character(len=clen), intent(out) :: ccount

    character(len=lenIoMax) :: numToStrg
    character(len=clen) :: strg

    if (clen.gt.lenIoMax) then
       stop 'int2char_adjtrim: clen > lenIoMax!'
    endif

    write(numToStrg,*) icount
    strg = trim(adjustl(numToStrg)) 
!    strg = adjustr(strg)
!    call replace_char(strg,' ','0')

    ccount(1:clen) = strg(1:clen)

  end subroutine int2char_adjtrim




  subroutine delete_file(name)
    implicit none
    character(len=*) :: name
    logical :: isthere

    inquire(file=name,exist=isthere)

    if (isthere) then
       open(unit=10,file=name)
       close(unit=10,status='delete')
    endif

  end subroutine delete_file





  function count_column(filename,delimiter)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: delimiter

    logical, parameter :: display = .false.

    integer :: count_column

    integer :: stat, j, num_delim
    character :: single_byte, CR, LF, column_delimiter

    integer, parameter :: nunit = 200

    LF = char(10) ! Line Feed
    CR = char(13) ! Carriage Return

    column_delimiter = delimiter    
    
    open (nunit, file=filename, form='unformatted', &
         access='direct', status='old', recl = 1, iostat=stat)
    if (stat .ne. 0) stop 'read_hearder: missing file!'
    
    ! process header line of the file
    j = 0
    num_delim = 0
    single_byte='o'

    do while ((single_byte .ne. CR) .and. (single_byte .ne. LF))
       j = j + 1
       read(nunit, rec = j) single_byte
       if (single_byte .eq. column_delimiter) then
          num_delim = num_delim + 1
       end if
       !write (*,'(I3,5x,I3,5x,A)') j, ichar(single_byte), single_byte
    end do
    close(nunit)

    if (display) then
       write (*,*)'delimiter ',delimiter,' found ',num_delim,'times'
    endif

    count_column = num_delim    
  end function count_column



  subroutine sp_livewrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(len=*) :: name    
    real(sp) :: x,a
    real(sp), optional :: b,c,d,e,f,g
      
    open(10,file=name,position='append',status='unknown')
    
    if (.not.present(b)) then
       write(10,100) x,a         
    elseif (.not.present(c)) then           
       write(10,100) x,a,b         
    elseif (.not.present(d)) then             
       write(10,100) x,a,b,c                    
    elseif (.not.present(e)) then         
       write(10,100) x,a,b,c,d                    
    elseif (.not.present(f)) then
       write(10,100) x,a,b,c,d,e            
    elseif (.not.present(g)) then
       write(10,100) x,a,b,c,d,e,f            
    else
       write(10,100) x,a,b,c,d,e,f,g            
    endif
    
    close(10)

100 format(8(ES25.16E3))

  end subroutine sp_livewrite


  subroutine dp_livewrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(len=*) :: name    
    real(dp) :: x,a
    real(dp), optional :: b,c,d,e,f,g
    
    open(10,file=name,position='append',status='unknown')
    
    if (.not.present(b)) then
       write(10,100) x,a         
    elseif (.not.present(c)) then           
       write(10,100) x,a,b         
    elseif (.not.present(d)) then             
       write(10,100) x,a,b,c                    
    elseif (.not.present(e)) then         
       write(10,100) x,a,b,c,d                    
    elseif (.not.present(f)) then
       write(10,100) x,a,b,c,d,e            
    elseif (.not.present(g)) then
       write(10,100) x,a,b,c,d,e,f            
    else
       write(10,100) x,a,b,c,d,e,f,g            
    endif
    
    close(10)
    
100 format(8(ES25.16E3))      
    
  end subroutine dp_livewrite



  subroutine qp_livewrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(len=*) :: name    
    real(qp) :: x,a
    real(qp), optional :: b,c,d,e,f,g
    
    open(10,file=name,position='append',status='unknown')
    
    if (.not.present(b)) then
       write(10,100) x,a         
    elseif (.not.present(c)) then           
       write(10,100) x,a,b         
    elseif (.not.present(d)) then             
       write(10,100) x,a,b,c                    
    elseif (.not.present(e)) then         
       write(10,100) x,a,b,c,d                    
    elseif (.not.present(f)) then
       write(10,100) x,a,b,c,d,e            
    elseif (.not.present(g)) then
       write(10,100) x,a,b,c,d,e,f            
    else
       write(10,100) x,a,b,c,d,e,f,g            
    endif
    
    close(10)
    
100 format(8(ES41.32E3))      
    
  end subroutine qp_livewrite


  subroutine sp_allwrite(name,x,a,b,c,d,e,f,g)
    implicit none
      character(*) :: name
      integer :: j,npts
      real(sp) :: x(:),a(:)
      real(sp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)

      npts=ubound(x,1)
      
      if (ubound(a,1).ne.npts) then
         write(*,*)'WARNING: vectors length differ'
      endif

!      write(*,*)'__write: save in ',name
      open(10,file=name,position='append',status='unknown')
      
      if (.not.present(b)) then
         do j=1,npts      
            write(10,100) x(j),a(j)
         enddo
      elseif (.not.present(c)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j)
         enddo
      elseif (.not.present(d)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j)            
         enddo
      elseif (.not.present(e)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j)            
         enddo
      elseif (.not.present(f)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j)            
         enddo
      elseif (.not.present(g)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
         enddo
      else
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
         enddo
      endif
      
      close(10)

100   format(8(ES25.16))      

    end subroutine sp_allwrite


    subroutine dp_allwrite(name,x,a,b,c,d,e,f,g)
      implicit none
      character(*) :: name
      integer :: j,npts
      real(dp) :: x(:),a(:)
      real(dp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)

      npts=ubound(x,1)
      
      if (ubound(a,1).ne.npts) then
         write(*,*)'WARNING: vectors length differ'
      endif

!      write(*,*)'__write: save in ',name
      open(10,file=name,position='append',status='unknown')
      
      if (.not.present(b)) then
         do j=1,npts      
            write(10,100) x(j),a(j)
         enddo
      elseif (.not.present(c)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j)
         enddo
      elseif (.not.present(d)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j)            
         enddo
      elseif (.not.present(e)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j)            
         enddo
      elseif (.not.present(f)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j)            
         enddo
      elseif (.not.present(g)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
         enddo
      else
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
         enddo
      endif
      
      close(10)

100   format(8(ES25.16))      

    end subroutine dp_allwrite

    subroutine qp_allwrite(name,x,a,b,c,d,e,f,g)
      implicit none
      character(*) :: name
      integer :: j,npts
      real(qp) :: x(:),a(:)
      real(qp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)

      npts=ubound(x,1)
      
      if (ubound(a,1).ne.npts) then
         write(*,*)'WARNING: vectors length differ'
      endif

!      write(*,*)'__write: save in ',name
      open(10,file=name,position='append',status='unknown')
      
      if (.not.present(b)) then
         do j=1,npts      
            write(10,100) x(j),a(j)
         enddo
      elseif (.not.present(c)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j)
         enddo
      elseif (.not.present(d)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j)            
         enddo
      elseif (.not.present(e)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j)            
         enddo
      elseif (.not.present(f)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j)            
         enddo
      elseif (.not.present(g)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
         enddo
      else
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
         enddo
      endif
      
      close(10)

100   format(8(ES41.32))      

    end subroutine qp_allwrite


    


  subroutine sp_binallwrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(*) :: name
    integer :: j,npts
    real :: x(:)
    real, optional :: a(:), b(:),c(:),d(:),e(:),f(:),g(:)
    
    integer :: datarecl
    integer :: recnum


    recnum = 0
    npts=ubound(x,1)

    if (present(a)) then
       if (ubound(a,1).ne.npts) then
          write(*,*)'WARNING: vectors length differ'
       endif
    end if
       
    write(*,*)'__write: save in ',name

    if (.not.present(a)) then
       datarecl=2*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j)
       enddo
    elseif (.not.present(b)) then
       datarecl=2*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j)
       enddo
    elseif (.not.present(c)) then
       datarecl=3*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j)
       enddo
    elseif (.not.present(d)) then
       datarecl=4*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j)            
       enddo
    elseif (.not.present(e)) then
       datarecl=5*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j)            
       enddo
    elseif (.not.present(f)) then
       datarecl=6*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j)            
       enddo
    elseif (.not.present(g)) then
       datarecl=7*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
       enddo
    else
       datarecl=8*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
       enddo
    endif
    
    close(10)
    
  end subroutine sp_binallwrite
  


  subroutine dp_binallwrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(*) :: name
    integer :: j,npts
    real(dp) :: x(:)
    real(dp), optional :: a(:), b(:),c(:),d(:),e(:),f(:),g(:)

    integer :: recnum
    integer :: datarecl

    npts=ubound(x,1)
    recnum=0

    if (present(a)) then
       if (ubound(a,1).ne.npts) then
          write(*,*)'WARNING: vectors length differ'
       endif
    end if
       
    write(*,*)'__write: save in ',name

    if (.not.present(a)) then
       datarecl=2*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j)
       enddo
    elseif (.not.present(b)) then
       datarecl = 4*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j)
       enddo
    elseif (.not.present(c)) then
       datarecl = 6*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum =recnum+1
          write(10,rec=recnum) x(j),a(j),b(j)
       enddo
    elseif (.not.present(d)) then
       datarecl = 8*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum =recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j)            
       enddo
    elseif (.not.present(e)) then
       datarecl = 10*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum = recnum + 1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j)            
       enddo
    elseif (.not.present(f)) then
       datarecl = 12*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum = recnum + 1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j)            
       enddo
    elseif (.not.present(g)) then
       datarecl = 14*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum = recnum + 1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
       enddo
    else
       datarecl = 16*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum = recnum + 1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
       enddo
    endif
    
    close(10)
        
  end subroutine dp_binallwrite



    
end module iotools

