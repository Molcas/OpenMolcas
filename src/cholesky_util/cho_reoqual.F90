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

subroutine Cho_ReoQual(iQuAB,MxQ,nSym,iQScr,IDK,nK,nQ)
!
! Purpose: reorder iQuAB array to IDK ordering.

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: MxQ, nSym, IDK(*), nK(nSym), nQ(nSym)
integer(kind=iwp), intent(inout) :: iQuAB(MxQ,nSym)
integer(kind=iwp), intent(out) :: iQScr(MxQ)
integer(kind=iwp) :: iK, iSym, kID, lK

kID = 0
do iSym=1,nSym
  if (nQ(iSym) > 0) then
    iQScr(1:nQ(iSym)) = iQuAB(1:nQ(iSym),iSym)
    do iK=1,nK(iSym)
      lK = IDK(kID+iK)
      iQuAB(iK,iSym) = iQScr(lK)
    end do
    kID = kID+nQ(iSym)
  else
    iQuAB(1:nK(iSym),iSym) = 0
  end if
end do

end subroutine Cho_ReoQual
