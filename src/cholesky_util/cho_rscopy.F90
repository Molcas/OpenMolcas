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

subroutine CHO_RSCOPY(IRS1,IRS2)
!
! Purpose: copy red. set info from location IRS1 to IRS2.
!          Special action is taken with INDRED if IRS1=1 so that it
!          will point as expected for the "current" reduced set.

use Cholesky, only: iiBstR, iiBstRSh, IndRed, nnBstR, nnBstRSh, nnBstRT
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IRS1, IRS2
integer(kind=iwp) :: IAB
integer(kind=iwp) :: MSYM

MSYM = size(iiBstRSh,1)
nnBstRSh(:,:,IRS2) = nnBstRSh(:,:,IRS1)
iiBstRSh(:,:,IRS2) = iiBstRSh(:,:,IRS1)
iiBstR(1:MSYM,IRS2) = iiBstR(1:MSYM,IRS1)
nnBstR(1:MSYM,IRS2) = nnBstR(1:MSYM,IRS1)
if (IRS1 == 1) then
  do IAB=1,size(INDRED,1)
    INDRED(IAB,IRS2) = IAB
  end do
else
  IndRed(:,iRS2) = IndRed(:,iRS1)
end if
NNBSTRT(IRS2) = NNBSTRT(IRS1)

end subroutine CHO_RSCOPY
