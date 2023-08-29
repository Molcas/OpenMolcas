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

function Cho_P_IndxParentDiag(iQ,iSym)
!
! Purpose: return index in global 1st reduced set of qualified iQ,
!          sym. iSym.

use Cholesky, only: Cho_Real_Par, IndRed, IndRed_G, iQuAB
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: Cho_P_IndxParentDiag
integer(kind=iwp), intent(in) :: iQ, iSym

if (Cho_Real_Par) then
  Cho_P_IndxParentDiag = IndRed_G(iQuAB(iQ,iSym),2)
else
  Cho_P_IndxParentDiag = IndRed(iQuAB(iQ,iSym),2)
end if

end function Cho_P_IndxParentDiag
