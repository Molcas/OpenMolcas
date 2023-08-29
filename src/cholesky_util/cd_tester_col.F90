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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine CD_Tester_Col(Col,nDim,iCol,nCol,Buf,lBuf)

use CD_Tester_mod, only: Mat
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, nCol, iCol(nCol), lBuf
real(kind=wp), intent(out) :: Col(nDim,nCol)
real(kind=wp), intent(in) :: Buf(lBuf)
integer(kind=iwp) :: i, kOff

#include "macros.fh"
unused_var(Buf)

do i=1,nCol
  kOff = nDim*(iCol(i)-1)
  Col(:,i) = Mat(kOff+1:kOff+nDim)
end do

return

end subroutine CD_Tester_Col
