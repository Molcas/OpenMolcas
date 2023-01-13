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
! Copyright (C) 1995, Roland Lindh                                     *
!               2000, Valera Veryazov                                  *
!               2014, Thomas Dresselhaus                               *
!***********************************************************************

subroutine MOEvalDer(MOValue,iDir,nMOs,nCoor,CCoor,CMOs,nCMO,DoIt)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: nDrv = 1, mAO = 4
integer(kind=iwp), intent(in) :: iDir, nMOs, nCoor, nCMO, DoIt(nMOs)
real(kind=wp), intent(out) :: MOValue(nCoor*nMOs)
real(kind=wp), intent(in) :: CCoor(3,nCoor), CMOs(nCMO)
integer(kind=iwp) :: I
real(kind=wp), allocatable :: MOTmp(:,:)

call mma_allocate(MOTmp,4,nCoor*nMOs)

call MOEval(MOTmp,nMOs,nCoor,CCoor,CMOs,nCMO,DoIt,nDrv,mAO)

! iDir = 1 then do dX
! iDir = 2 then do dY
! iDir = 3 then do dZ
!write(u6,*) 'iDir:',iDir
if ((iDir > 0) .and. (iDir < 4)) then
  MOValue(:) = MOTmp(iDir,:)
else ! do gradient
  do I=1,nCoor*nMOs
    MOValue(I) = MOTmp(2,I)+MOTmp(3,I)+MOTmp(4,I)
  end do
end if
call mma_deallocate(MOTmp)

return

end subroutine MOEvalDer
