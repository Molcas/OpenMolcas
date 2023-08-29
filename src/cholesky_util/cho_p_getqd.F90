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

subroutine Cho_P_GetQD(QD)
!
! Purpose: copy qualified diagonal elements from global diagonal to
!          array QD.

use Cholesky, only: Diag_G, IndRed_G, iQuAB, nQual, nSym
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: QD(*)
integer(kind=iwp) :: iAB, iQ, iSym, kQD

kQD = 0
do iSym=1,nSym
  do iQ=1,nQual(iSym)
    iAB = IndRed_G(iQuAB(iQ,iSym),2)
    QD(kQD+iQ) = Diag_G(iAB)
  end do
  kQD = kQD+nQual(iSym)
end do

end subroutine Cho_P_GetQD
