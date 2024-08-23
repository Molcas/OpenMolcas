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
      use getline_mod, only: Line, iGetLine, nCol, iStrt, iEnd,         &
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
