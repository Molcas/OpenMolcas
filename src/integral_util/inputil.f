!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Put_ln(In_Line)
!***********************************************************************
! This function replaces function getln                                *
!                                                                      *
! It reads, broadcasts, and parses an input line                       *
!                                                                      *
! Blank lines or lines containing star (*) in column 1 are skipped     *
! Lines staring with exclaimation (!) in column 1 are skipped          *
!                                                                      *
! After this routine has been called, data can be retrieved using      *
! the subroutines                                                      *
!                                                                      *
!   Get_F(icol,array,n)  (for floating point values)                   *
!   Get_I(icol,iarry,n)  (for integer values)                          *
!   Get_S(icol,strgs,n)  (for character strings)                       *
!                                                                      *
! where icol is the first non-blank work to be taken, and n is the     *
! number of data.                                                      *
!                                                                      *
!***********************************************************************
      use getline_mod, only: Line, nCol, iStrt, iEnd
      implicit None
      Character(LEN=*) In_line

      Integer i, j, l, icom

      Line=In_Line
      l=len(line)
      Do i = 1, l
         if(ichar(line(i:i)).eq.9) line(i:i)=' '
         If (line(i:i).eq.';') line(i:l)=' '
      End Do
      ncol=0
      j=1
      Do
         icom=0
         do i=j,l
           if(line(i:i).eq.',') then
             icom=icom+1
             if(icom.ne.1) goto 20
           else
             if(line(i:i).ne.' ') goto 20
           end if
         end do
         Exit

 20      do j=i,l
           if(line(j:j).eq.' '.or.line(j:j).eq.',') goto 30
         end do
         j=l+1
 30      ncol=ncol+1
         istrt(ncol)=i
         iend(ncol)=j-1
      End Do

      End Subroutine Put_ln

      character(LEN=180) function get_ln(lunit)
      use getline_mod, only: Quit_On_Error
      Implicit None
      Integer lunit
      character(LEN=180), External::  get_ln_quit
      get_ln=get_ln_quit(lunit,1)
      if(Quit_On_Error) Then
        Call WarningMessage(2,'Error in Get_Ln')
        Call Quit_OnUserError()
      End If
      Return
      End function get_ln

      character(LEN=180) function get_ln_EOF(lunit)
      use getline_mod, only: Quit_On_Error
      Implicit None
      Integer lunit
      character(LEN=180), External::  get_ln_quit
      get_ln_EOF=get_ln_quit(lunit,0)
      if(Quit_On_Error) get_ln_EOF='EOF'
      End function get_ln_EOF

      character(LEN=180) function get_ln_quit(lunit,icritical)
!***********************************************************************
! This function replaces function getln                                *
!                                                                      *
! It reads, broadcasts, and parses an input line                       *
!                                                                      *
! Blank lines or lines containing star (*) in column 1 are skipped     *
! Lines staring with exclaimation (!) in column 1 are skipped          *
!                                                                      *
! After this routine has been called, data can be retrieved using      *
! the subroutines                                                      *
!                                                                      *
!   Get_F(icol,array,n)  (for floating point values)                   *
!   Get_I(icol,iarry,n)  (for integer values)                          *
!   Get_S(icol,strgs,n)  (for character strings)                       *
!                                                                      *
! where icol is the first non-blank work to be taken, and n is the     *
! number of data.                                                      *
!                                                                      *
!***********************************************************************
      use getline_mod, only: Line, iGetLine, nCol, iStrt, iEnd,
     &                       MyUnit, Quit_On_Error
      implicit None
      Integer lUnit, iCritical

      Integer i, j, l, icom
      Character(LEN=256) filename
      Quit_On_Error=.false.
      myunit=lunit
1     read(lunit,'(A)',err=100,end=200) line
        igetline=igetline+1
      if(line.eq.' '.or.line(1:1).eq.'*'.or.line(1:1).eq.'!') goto 1
      l=len(line)
      Do i = 1, l
         If (ichar(line(i:i)).eq.9) line(i:i)=' '
         If (line(i:i).eq.';') line(i:l)=' '
      End Do
      ncol=0
      j=1
 10   icom=0
      do i=j,l
        if(line(i:i).eq.',') then
          icom=icom+1
          if(icom.ne.1) goto 20
        else
          if(line(i:i).ne.' ') goto 20
        end if
      end do
      get_ln_quit=line
      return
 20   do j=i,l
        if(line(j:j).eq.' ' .or.  line(j:j).eq.',') goto 30
      end do
      j=l+1
 30   ncol=ncol+1
      istrt(ncol)=i
      iend(ncol)=j-1
      goto 10
!
100   filename=' '
      inquire(unit=lunit,name=filename)
      if(filename.ne.' ') then
        write(6,'(a,a)') 'Error reading file=',filename
      else
        write(6,'(a,i8)') 'Error reading unit=',lunit
      endif
        write(6,'(a)') 'Line: ',line(1:80)
      Quit_On_Error=.true.
200   if(icritical.eq.0) then
       Quit_On_Error=.true.
      return
      endif
!
      filename=' '
      inquire(unit=lunit,name=filename)
      if(filename.ne.' ') then
        write(6,'(a,a)') 'EOF reached for file=',filename
      else
        write(6,'(a,i8)') 'EOF reached for unit=',lunit
      endif

      Quit_On_Error=.true.
      end function get_ln_quit
!
      subroutine Get_F(icol,val,n)
      use getline_mod, only: Line, nCol, iStrt, iEnd
      implicit None
      integer iCol, n
      Real*8 val(n)

      Integer ic, i, i1, i2
      Character(LEN=80) string

      ic=icol
      do i=1,n
        if(ic.le.ncol) then
          i1=istrt(ic)
          i2=iend(ic)
          if(i1.le.i2) then
            string=' '
            string(LEN(string)+i1-i2:)=line(i1:i2)
            read(string,'(F80.0)',err=600,end=600) val(i)
          else
            val(i)=0.0D0
          end if
          ic=ic+1
        else
          write(6,110) icol+n-1,line
110       format(/' ERROR IN GET_F: TRYING TO READ',i4,' VALUES'/1x,a)
          Call FindErrorLine()
          Call WarningMessage(2,'Error in Get_F')
          Call Quit_OnUserError()
        end if
      end do
      return
!
600       Call FindErrorLine()
          Call WarningMessage(2,'Error in Get_F')
          Call Quit_OnUserError()
      end subroutine Get_F
!
      subroutine Get_F1(icol,val)
      implicit None
      Integer icol
      Real*8 val

      Real*8 dum(1)
      call Get_F(icol,dum,1)
      val=dum(1)
      end subroutine Get_F1

      subroutine Get_I(icol,ival,n)
      use getline_mod, only: Line, nCol, iStrt, iEnd
      implicit None
      Integer iCol, n
      integer ival(n)

      Integer ic, i, i1, i2
      Character(LEN=80) string
      ic=icol
      do i=1,n
        if(ic.le.ncol) then
          i1=istrt(ic)
          i2=iend(ic)
          if(i1.le.i2) then
            string=' '
            string(LEN(string)+i1-i2:)=line(i1:i2)
            read(string,'(i80)',err=600,end=600) ival(i)
          else
            ival(i)=0
          end if
          ic=ic+1
        else
          write(6,210) icol+n-1,line
210       format(/' ERROR IN GET_I: TRYING TO READ',i4,' VALUES'/1x,a)
          Call FindErrorLine()
          Call WarningMessage(2,'Error in Get_I')
          Call Quit_OnUserError()
        end if
      end do
      return
600       Call FindErrorLine()
          Call WarningMessage(2,'Error in Get_I')
          Call Quit_OnUserError()
      end subroutine Get_I
!
      subroutine Get_I1(icol,ival)
      implicit None
      integer iCol, ival

      integer idum(1)
      call Get_I(icol,idum,1)
      ival=idum(1)
      end subroutine Get_I1
!
      Subroutine Get_S(icol,str,n)
      use getline_mod
      Implicit None
      Integer iCol, n
      character(LEN=*) str(n)

      Integer ic, i
      ic=icol
      do i=1,n
        if(ic.le.ncol) then
          if(iend(ic).ge.istrt(ic)) then
            str(i)=line(istrt(ic):iend(ic))
          else
            str(i)=' '
          end if
          ic=ic+1
        else
          write(6,210) icol+n-1,line
210       format(/' ERROR IN GET_S: TRYING TO READ',i4,' STRINGS'/1x,a)
          Call FindErrorLine()
          Call WarningMessage(2,'Error in Get_S')
          Call Quit_OnUserError()
        end if
      end do
      end Subroutine Get_S

      subroutine Read_v(lunit,work,istrt,iend,inc,ierr)
      implicit None
      Integer lUnit, iStrt, iEnd, Inc, iErr
      real*8 work(iend)

      Integer i
      ierr=0
      read(lunit,*,err=100) (work(istrt+i),i=0,iend-istrt,inc)
      goto 110
100   ierr=1
110   continue
      return
      end subroutine Read_v

      subroutine Read_iv(lunit,iwork,istrt,iend,inc,ierr)
      implicit None
      Integer lUnit, iStrt, iEnd, Inc, iErr
      integer iwork(iend)

      Integer i
      ierr=0
      read(lunit,*,err=100) (iwork(istrt+i),i=0,iend-istrt,inc)
      goto 110
100   ierr=1
110   continue
      return
      end subroutine Read_iv

      subroutine FindErrorLine()
      use getline_mod, only: MyUnit, iGetLine
      Implicit None

      character(LEN=180) line
      Integer lunit, isave
      lunit=myunit
      isave=igetline
      rewind (lunit)
 2    read(lunit,'(a)', end=300) Line
      Call UpCase(Line)
      Line = adjustl(Line)
      if(Line(1:1).eq.'&') then
         line=line(2:)
          goto 3
      endif
      goto 2
 3    igetline=0
      write(6,'(a,a,a)') ' >>>>> Input file for module ',
     *  line(1:index(line,' ')),' <<<<<'
 1    read(lunit,'(A)',err=100,end=200) line
       igetline=igetline+1
       if(igetline.eq.isave) then
        write (6,*) '******   Error  *******'
        write (6,'(a)') line
        write (6,'(a)')
        Call WarningMessage(2,'Error in FindErrorLine')
        call Quit_OnUserError()
       endif
       if(isave-igetline.le.50) then
        write (6,'(a)') line
       endif
       goto 1
!       write(6,'(a)') ' >>>>> Input error <<<<<'
!       rewind(lunit)
!       igetline=0
!       goto 1
100   continue
200   continue
300   continue
      Call WarningMessage(1,'FindErrorLine:'//
     & ' Error in input was not located;'//
     & '  Please, check it manually!')
      return
      end subroutine FindErrorLine

      subroutine ResetErrorLine
      use getline_mod, only: igetline
      Implicit None
      igetline=0
      end subroutine ResetErrorLine
