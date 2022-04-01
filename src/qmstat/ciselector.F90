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

subroutine CiSelector(nEqState,nState,STC,nCIRef,iCIInd,dCIRef)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: nEqState
integer(kind=iwp), intent(in) :: nState, nCIRef, iCIInd(nCIRef)
real(kind=wp), intent(in) :: STC(nState,nState), dCIRef(nCIRef)
integer(kind=iwp) :: indMAX, iRef, iState
real(kind=wp) :: dScal, dScalMAX

! Initial stuff

dScalMAX = Zero
indMAX = 1

! Compute relevant scalar products

do iState=1,nState
  dScal = Zero
  do iRef=1,nCIRef
    dScal = dScal+STC(iCIInd(iRef),iState)*dCIRef(iRef)
  end do
  dScal = abs(dScal)

  ! Test if largest

  if (dScal > dScalMAX) then
    dScalMAX = dScal
    indMAX = iState
  end if
end do

! If maximum overlap is small, scream!

if (dScalMAX < sqrt(Half)) then
  write(6,*)
  write(6,*) '   WARNING! Less than 50% of CISElect reference found. Consider to redefine reference!'
end if

! Now set nEqState

nEqState = indMAX

! Auf Wiedersehen

return

end subroutine CiSelector
