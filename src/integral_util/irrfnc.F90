!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

function IrrFnc(iFnc)

use Symmetry_Info, only: nIrrep, iOper

implicit none
integer IrrFnc
integer iFnc
integer iTest(8), i, ix, iy, iz, jx, jy, jz, iCh
integer, external :: iNew

ix = iand(iFnc,1)
iy = iand(iFnc,2)/2
iz = iand(iFnc,4)/4
do i=1,nIrrep
  jx = iand(iOper(i-1),1)
  jy = iand(iOper(i-1),2)/2
  jz = iand(iOper(i-1),4)/4
  iCh = 1
  if ((ix /= 0) .and. (jx /= 0)) iCh = -iCh
  if ((iy /= 0) .and. (jy /= 0)) iCh = -iCh
  if ((iz /= 0) .and. (jz /= 0)) iCh = -iCh
  iTest(i) = iCh
end do
IrrFnc = iNew(iTest,nIrrep)-1

return

end function IrrFnc
