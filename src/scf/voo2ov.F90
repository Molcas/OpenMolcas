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
! Copyright (C) 1992, Martin Schuetz                                   *
!               2017, Roland Lindh                                     *
!***********************************************************************

subroutine vOO2OV(v1,nOO,v2,mOV,nD,kOV)

use Interfaces_SCF, only: vOO2OV_inner
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOO, mOV, nD, kOV(nD)
real(kind=wp), intent(in) :: v1(nOO,nD)
real(kind=wp), intent(out) :: v2(mOV)
integer(kind=iwp) :: iD, iEnd, iSt

iEnd = 0
v2(:) = Zero
do iD=1,nD
  if (kOV(iD) < 1) cycle
  iSt = iEnd+1
  iEnd = iEnd+kOV(iD)
  call vOO2OV_inner(v1(:,iD),nOO,v2(iSt:iEnd),kOV(iD),iD)
end do

end subroutine vOO2OV
