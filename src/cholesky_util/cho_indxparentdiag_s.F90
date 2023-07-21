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

integer function Cho_IndxParentDiag_S(iQ,iSym)

use ChoSwp, only: iQuAB, IndRed

implicit none
integer iQ, iSym
#include "cholesky.fh"

Cho_IndxParentDiag_S = IndRed(iQuAB(iQ,iSym),2)

end function Cho_IndxParentDiag_S
