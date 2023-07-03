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

subroutine FixEqualSign2(Line,LuRd,Lu_UDIC,iRow,New_Line)

character*(*) Line
character*180 Temp_Line
character*180 Get_Ln
external Get_Ln

nLine = len(Line)
if (nLine > len(Temp_Line)) then
  call WarningMessage(2,'Error in FixEqualSign!')
  call Abend()
end if

Temp_Line = adjustl(Line)
ix = iCLast(Temp_Line,nLine)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the next line and determine if the lines should be merged.

Line = Get_Ln(LuRd)
Line = adjustl(Line)
iy = iCLast(Line,nLine)
call UpCase(Line)
!                                                                      *
!***********************************************************************
!                                                                      *
if (index(Line(1:iy),'END ') == 1) then
  iRow = iRow+1
  write(Lu_UDIC,'(A)') Temp_Line
  New_Line = 2

else if (index(Line(1:iy),' ') == 0) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! If the line contains two or more items we should merge the lines.

  ! Just one item

  iRow = iRow+1
  write(Lu_UDIC,'(A)') Temp_Line
  New_Line = 1

else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  Temp_Line(ix+2:ix+2) = '='
  ix = ix+2
  if (ix+2+iy > nLine) then
    call WarningMessage(2,'Problems merging lines!')
    call Abend()
  end if
  Temp_Line(ix+2:nLine) = Line(1:nLine-ix-1)
  Line = Temp_Line
  call UpCase(Line)
  New_Line = 0
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine FixEqualSign2
