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
      write(6,'(a,a,a)') ' >>>>> Input file for module ',               &
     &  line(1:index(line,' ')),' <<<<<'
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
      Call WarningMessage(1,'FindErrorLine:'//                          &
     & ' Error in input was not located;'//                             &
     & '  Please, check it manually!')
      return
      end subroutine FindErrorLine
