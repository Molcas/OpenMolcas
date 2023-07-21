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

implicit none
integer MxQ, nSym
integer iQuAB(MxQ,nSym)
integer iQScr(MxQ)
integer IDK(*)
integer nK(nSym)
integer nQ(nSym)
integer iSym, kID
integer iK, lK

kID = 0
do iSym=1,nSym
  if (nQ(iSym) > 0) then
    call iCopy(nQ(iSym),iQuAB(1,iSym),1,iQScr,1)
    do iK=1,nK(iSym)
      lK = IDK(kID+iK)
      iQuAB(iK,iSym) = iQScr(lK)
    end do
    kID = kID+nQ(iSym)
  else
    do iK=1,nK(iSym)
      iQuAB(iK,iSym) = 0
    end do
  end if
end do

end subroutine Cho_ReoQual
