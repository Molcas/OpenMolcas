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

integer function Cho_P_IndxParentDiag(iQ,iSym)
!
! Purpose: return index in global 1st reduced set of qualified iQ,
!          sym. iSym.

implicit none
integer iQ, iSym
#include "cho_para_info.fh"
integer Cho_IndxParentDiag_S
external Cho_IndxParentDiag_S
integer Cho_IndxParentDiag_P
external Cho_IndxParentDiag_P

if (Cho_Real_Par) then
  Cho_P_IndxParentDiag = Cho_IndxParentDiag_P(iQ,iSym)
else
  Cho_P_IndxParentDiag = Cho_IndxParentDiag_S(iQ,iSym)
end if

end function Cho_P_IndxParentDiag
