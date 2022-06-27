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

subroutine ColChg(q,QMAX,QMIN,C1,C2,C3)

use Constants, only: Zero, Half
use Definitions, only: wp, u6

implicit none
real(kind=wp), intent(in) :: q, QMAX, QMIN
real(kind=wp), intent(out) :: C1, C2, C3
real(kind=wp) :: DNeg, DPos
character(len=20) :: Colour

! Assign tesserae colours for GeomView:

DPos = QMax*Half
DNeg = QMin*Half
! Total charge density < 0
if (q < DNeg) then
  Colour = 'Dark Blue'
else if (q < Zero) then
  Colour = 'Light Blue'
! Total charge density > 0
else if (q < DPos) then
  Colour = 'Pink'
else
  Colour = 'Red'
end if
call ColTss(u6,Colour,C1,C2,C3)

return

end subroutine ColChg
