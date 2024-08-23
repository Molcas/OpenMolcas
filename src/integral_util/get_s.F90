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
