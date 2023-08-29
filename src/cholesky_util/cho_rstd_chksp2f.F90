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

subroutine Cho_RstD_ChkSP2F(iSP2F,l_iSP2F,nErr)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: l_iSP2F, iSP2F(l_iSP2F)
integer(kind=iwp), intent(out) :: nErr
integer(kind=iwp) :: i
integer(kind=iwp), allocatable :: iChk(:)

call mma_allocate(iChk,l_iSP2F,Label='iChk')

call Cho_RstD_GetInd3(iChk,l_iSP2F)

nErr = 0
do i=1,l_iSP2F
  if (iChk(i) /= iSP2F(i)) nErr = nErr+1
end do

call mma_deallocate(iChk)

end subroutine Cho_RstD_ChkSP2F
