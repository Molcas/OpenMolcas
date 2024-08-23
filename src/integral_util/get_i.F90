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
