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
! Copyright (C) 1998, Roland Lindh                                     *
!***********************************************************************

function iPD(iSO_,jSO_,iSOSym,nSOs)

use Basis_Info, only: nBas

implicit none
integer iPD
integer iSO_, jSO_, nSOs
integer iSOSym(2,nSOs)
integer iSO, jSO, iSym, iSOr, jSym, jSOr, ij

iPD = -999999

iSO = max(iSO_,jSO_)
jSO = min(iSO_,jSO_)
iSym = iSOSym(1,iSO)
iSOr = iSOSym(2,iSO)
jSym = iSOSym(1,jSO)
jSOr = iSOSym(2,jSO)
if (iSym == jSym) then
  ij = iSOr*(iSOr-1)/2+jSOr
else
  ij = (iSOr-1)*nBas(jSym)+jSOr
end if

iPD = ij

return

end function iPD
