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

function Cho_IndxParentDiag_S(iQ,iSym)

use Cholesky, only: IndRed, iQuAB
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: Cho_IndxParentDiag_S
integer(kind=iwp) :: iQ, iSym
#include "cholesky.fh"

Cho_IndxParentDiag_S = IndRed(iQuAB(iQ,iSym),2)

end function Cho_IndxParentDiag_S
