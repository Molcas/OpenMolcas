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

subroutine PPMmG(nHer,MmPPG,la,lb,lr)

nElem(i) = (i+1)*(i+2)/2

nHer = 0
MmPPG = 0

laplb = max(nElem(la+1),nElem(lb))**2
MmPPG = MmPPG+2*laplb

if (la > 0) then
  lamlb = max(nElem(la-1),nElem(lb))**2
else
  lamlb = 0
end if
MmPPG = MmPPG+2*lamlb

lalbp = max(nElem(la),nElem(lb+1))**2
MmPPG = MmPPG+2*lalbp

if (lb > 0) then
  lalbm = max(nElem(la),nElem(lb-1))**2
else
  lalbm = 0
end if
MmPPG = MmPPG+2*lalbm

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(lr)

end subroutine PPMmG
