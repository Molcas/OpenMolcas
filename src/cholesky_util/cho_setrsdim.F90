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

subroutine CHO_SETRSDIM(NDIMRS,MSYM,MRED,IRED,ILOC)
!
! Purpose: set reduced set dimension.

use Cholesky, only: MaxRed, nnBstR, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: MSYM, MRED, IRED, ILOC
integer(kind=iwp), intent(inout) :: NDIMRS(MSYM,MRED)

if (IRED <= MAXRED) NDIMRS(:,IRED) = NNBSTR(1:NSYM,ILOC)

end subroutine CHO_SETRSDIM
