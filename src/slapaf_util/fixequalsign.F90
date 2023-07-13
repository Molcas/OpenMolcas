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

subroutine FixEqualSign(Line,LuRd)

use Definitions, only: iwp

implicit none
character(len=*), intent(inout) :: Line
integer(kind=iwp), intent(in) :: LuRd
integer(kind=iwp) :: ix, iy, nLine
character(len=180) :: Temp_Line
integer(kind=iwp), external :: iCLast
character(len=180), external :: Get_Ln

nLine = len(Line)
if (nLine > len(Temp_Line)) then
  call WarningMessage(2,'Error in FixEqualSign!')
  call Abend()
end if

Temp_Line = adjustl(Line)
ix = iCLast(Temp_Line,nLine)
Temp_Line(ix+2:ix+2) = '='
ix = ix+2

Line = Get_Ln(LuRd)
Line = adjustl(Line)
iy = iCLast(Line,nLine)
if (ix+2+iy > nLine) then
  call WarningMessage(2,'Problems merging lines!')
  call Abend()
end if
Temp_Line(ix+2:nLine) = Line(1:nLine-ix-1)
Line = Temp_Line
call UpCase(Line)

return

end subroutine FixEqualSign
