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

use CDTHLP, only: Mat
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nDim, nCol, iCol(nCol), lBuf
real(kind=wp) :: Col(nDim,nCol), Buf(lBuf)
integer(kind=iwp) :: i, kOff

do i=1,nCol
  kOff = 1+nDim*(iCol(i)-1)
  call dCopy_(nDim,Mat(kOff),1,Col(1,i),1)
end do

return

#ifdef _WARNING_WORKAROUND_
! Avoid unused argument warnings
if (.false.) call Unused_real_array(Buf)
#endif

end subroutine CD_Tester_Col
