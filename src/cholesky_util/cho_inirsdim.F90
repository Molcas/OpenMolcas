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

use ChoArr, only: nDimRS

implicit real*8(a-h,o-z)
#include "cholesky.fh"

if (RSTCHO) then
  ILOC = 3
  do IRS=1,XNPASS
    call CHO_GETRED(IRS,ILOC,.false.)
    call CHO_SETREDIND(ILOC)
    call ICOPY(NSYM,NNBSTR(:,ILOC),1,nDimRS(:,iRS),1)
  end do
  NSET = XNPASS
else
  call ICOPY(NSYM,NNBSTR(1,1),1,nDimRS,1)
  NSET = 1
end if

NUM = NSYM*(MAXRED-NSET)
call IZERO(nDimRS(1,NSET+1),NUM)

end subroutine CHO_INIRSDIM
