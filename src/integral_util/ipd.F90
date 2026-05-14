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

use Index_Functions, only: nTri_Elem
use Basis_Info, only: nBas
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iPD
integer(kind=iwp), intent(in) :: iSO_, jSO_, nSOs, iSOSym(2,nSOs)
integer(kind=iwp) :: ij, iSO, iSOr, iSym, jSO, jSOr, jSym

iPD = -999999

iSO = max(iSO_,jSO_)
jSO = min(iSO_,jSO_)
iSym = iSOSym(1,iSO)
iSOr = iSOSym(2,iSO)
jSym = iSOSym(1,jSO)
jSOr = iSOSym(2,jSO)
if (iSym == jSym) then
  ij = nTri_Elem(iSOr-1)+jSOr
else
  ij = (iSOr-1)*nBas(jSym)+jSOr
end if

iPD = ij

return

end function iPD
