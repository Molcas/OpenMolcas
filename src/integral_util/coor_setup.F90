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

subroutine Coor_Setup(iSD4,nSD,Coor)

use Basis_Info, only: dbsc
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSD, iSD4(0:nSD,4)
real(kind=wp), intent(out) :: Coor(3,4)
integer(kind=iwp) :: iCnt, iCnttp, jCnt, jCnttp, kCnt, kCnttp, lCnt, lCnttp

iCnttp = iSD4(13,1)
iCnt = iSD4(14,1)
jCnttp = iSD4(13,2)
jCnt = iSD4(14,2)
kCnttp = iSD4(13,3)
kCnt = iSD4(14,3)
lCnttp = iSD4(13,4)
lCnt = iSD4(14,4)

if (dbsc(iCnttp)%Aux) then
  Coor(:,1) = dbsc(jCnttp)%Coor(:,jCnt)
else
  Coor(:,1) = dbsc(iCnttp)%Coor(:,iCnt)
end if
Coor(:,2) = dbsc(jCnttp)%Coor(:,jCnt)

if (dbsc(kCnttp)%Aux) then
  Coor(:,3) = dbsc(lCnttp)%Coor(:,lCnt)
else
  Coor(:,3) = dbsc(kCnttp)%Coor(:,kCnt)
end if
Coor(:,4) = dbsc(lCnttp)%Coor(:,lCnt)

end subroutine Coor_Setup
