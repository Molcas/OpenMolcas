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

use Definitions, only: wp, iwp

implicit none
integer, parameter :: nDrv = 1, mAO = 4
integer(kind=iwp), intent(in) :: iDir, nMOs, nCoor, nCMO, DoIt(nMOs)
real(kind=wp), intent(out) :: MOValue(nCoor*nMOs)
real(kind=wp), intent(in) :: CCoor(3,nCoor), CMOs(nCMO)
#include "WrkSpc.fh"
integer(kind=iwp) :: I, IJ, IJX, IJY, IJZ, iMoTmp

call GetMem('MOTMP','Allo','Real',iMoTmp,4*nCoor*nMOs)

call MOEval(Work(iMoTmp),nMOs,nCoor,CCoor,CMOs,nCMO,DoIt,nDrv,mAO)

! iDir = 1 then do dX
! iDir = 2 then do dY
! iDir = 3 then do dZ
!write(u6,*) 'iDir:',iDir
if (iDir > 0 .and. iDir < 4) then
  do I=1,nCoor*nMOs
    IJ = iDir+1+(I-1)*4
    MOValue(I) = Work(iMoTmp-1+IJ)
  end do
else ! do gradient
  do I=1,nCoor*nMOs
    IJX = 2+(I-1)*4
    IJY = 3+(I-1)*4
    IJZ = 4+(I-1)*4
    MOValue(I) = Work(iMoTmp-1+IJX)+Work(iMoTmp-1+IJY)+Work(iMoTmp-1+IJZ)
  end do
end if
call GetMem('MOTMP','Free','Real',iMoTmp,4*nCoor*nMOs)

return

end subroutine MOEvalDer
