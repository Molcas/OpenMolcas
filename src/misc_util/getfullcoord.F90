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

subroutine GetFullCoord(Coor,FMass,FAtLbl,nFAtoms,lSlapaf)

use Constants, only: One, uToau
use Definitions, only: wp, iwp

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nFAtoms
real(kind=wp), intent(out) :: Coor(3,nFAtoms), FMass(nFAtoms)
character(len=LenIn), intent(out) :: FAtLbl(nFAtoms)
logical(kind=iwp), intent(in) :: lSlapaf
integer(kind=iwp) :: i, iAt, iOper(8), jAt, mCenter, nAtoms, nCenter, nOper, nSym
real(kind=wp) :: AMass, Cnew(3), Cold(3), RotVec(3)
character(len=LenIn) :: Byte4
integer(kind=iwp) :: jOper(3) = [2,3,5]
logical(kind=iwp), external :: EQ

call Get_iScalar('nSym',nSym)
call Get_iArray('Symmetry operations',iOper,nSym)
call Get_iScalar('Unique atoms',nAtoms)
if (nAtoms > nFAtoms) call SysAbendMsg('GetFullCoord','nAtoms > nFAtoms','')
call Get_cArray('Unique Atom Names',FAtLbl,LenIn*nAtoms)
if (lSlapaf) then
  call Get_dArray('Initial Coordinates',Coor,3*nAtoms)
else
  call Get_dArray('Unique Coordinates',Coor,3*nAtoms)
end if

call Get_Mass(FMass,nAtoms)
FMass(1:nAtoms) = FMass(1:nAtoms)/uToau

nCenter = nAtoms

if (nSym /= 1) then

  nOper = 0
  if (nSym == 2) nOper = 1
  if (nSym == 4) nOper = 2
  if (nSym == 8) nOper = 3
  nCenter = nAtoms
  do i=1,nOper
    RotVec(1) = merge(-One,One,btest(iOper(jOper(i)),0))
    RotVec(2) = merge(-One,One,btest(iOper(jOper(i)),1))
    RotVec(3) = merge(-One,One,btest(iOper(jOper(i)),2))
    mCenter = nCenter
    outer: do iAt=1,mCenter
      Cold(:) = Coor(:,iAt)
      Byte4 = FAtLbl(iAt)
      AMass = FMass(iAt)
      Cnew(:) = RotVec*Cold
      do jAt=1,nCenter
        if (Byte4 /= FAtLbl(jAt)) cycle
        if (EQ(Cnew,Coor(:,jAt))) cycle outer
      end do
      nCenter = nCenter+1
      if (nCenter > nFAtoms) call SysAbendMsg('GetFullCoord','nCenter > nFAtoms','')
      Coor(:,nCenter) = Cnew
      FAtLbl(nCenter) = Byte4
      FMass(nCenter) = AMass
    end do outer
  end do

end if

if (nCenter /= nFAtoms) call SysAbendMsg('GetFullCoord','nCenter /= nFAtoms','')

return

end subroutine GetFullCoord
