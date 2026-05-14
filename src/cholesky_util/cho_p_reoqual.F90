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

subroutine Cho_P_ReoQual(iQScr,IDK,nK)

use Cholesky, only: Cho_Real_Par, iQuAB, MaxQual, nQual, nSym
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: iQScr(*)
integer(kind=iwp), intent(in) :: IDK(*), nK(*)

call Cho_ReoQual(iQuAB,MaxQual,nSym,iQScr,IDK,nK,nQual)
if (Cho_Real_Par) then
  call Cho_P_QualSwp()
  call Cho_ReoQual(iQuAB,MaxQual,nSym,iQScr,IDK,nK,nQual)
  call Cho_P_QualSwp()
end if

end subroutine Cho_P_ReoQual
