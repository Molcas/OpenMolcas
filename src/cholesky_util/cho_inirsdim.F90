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

subroutine CHO_INIRSDIM()
!
! Purpose: initialize reduced set dimension.

use Cholesky, only: MaxRed, nDimRS, nnBstR, nSym, RstCho, XnPass
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ILOC, IRS, NSET

if (RSTCHO) then
  ILOC = 3
  do IRS=1,XNPASS
    call CHO_GETRED(IRS,ILOC,.false.)
    call CHO_SETREDIND(ILOC)
    nDimRS(1:NSYM,iRS) = NNBSTR(1:NSYM,ILOC)
  end do
  NSET = XNPASS
else
  nDimRS(1:NSYM,1) = NNBSTR(1:NSYM,1)
  NSET = 1
end if

nDimRS(:,NSET+1:MAXRED) = 0

end subroutine CHO_INIRSDIM
