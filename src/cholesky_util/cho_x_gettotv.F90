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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_X_GetTotV(nV,n)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: return total number of vectors (over all nodes).

use Cholesky, only: NumCho_Bak
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: n, nV(n)
#include "cholesky.fh"
#include "cho_para_info.fh"
integer(kind=iwp) :: iSym

#ifdef _DEBUGPRINT_
if (n < nSym) call Cho_Quit('Illegal input variable in Cho_X_GetTotV',104)
#endif

if (Cho_Real_Par) then
  do iSym=1,nSym
    nV(iSym) = NumCho_Bak(iSym)
  end do
else
  do iSym=1,nSym
    nV(iSym) = NumCho(iSym)
  end do
end if

end subroutine Cho_X_GetTotV
