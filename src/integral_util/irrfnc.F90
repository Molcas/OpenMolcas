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

use Symmetry_Info, only: iOper, nIrrep
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IrrFnc
integer(kind=iwp), intent(in) :: iFnc
integer(kind=iwp) :: i, iCh, iTest(8)
logical(kind=iwp) :: ix, iy, iz, jx, jy, jz
integer(kind=iwp), external :: iNew

ix = btest(iFnc,0)
iy = btest(iFnc,1)
iz = btest(iFnc,2)
do i=1,nIrrep
  jx = btest(iOper(i-1),0)
  jy = btest(iOper(i-1),1)
  jz = btest(iOper(i-1),2)
  iCh = 1
  if (ix .and. jx) iCh = -iCh
  if (iy .and. jy) iCh = -iCh
  if (iz .and. jz) iCh = -iCh
  iTest(i) = iCh
end do
IrrFnc = iNew(iTest,nIrrep)-1

return

end function IrrFnc
