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

subroutine CiSelector(nEqState,nState,iSTC,nCIRef,iCIInd,dCIRef)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
#include "maxi.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: nEqState, nState, iSTC, nCIRef, iCIInd(MxState)
real(kind=wp) :: dCIRef(MxState)
integer(kind=iwp) :: indBase, indMAX, indx, iRef, iState
real(kind=wp) :: dScal, dScalMAX

! Initial stuff

dScalMAX = Zero
indMAX = 1

! Compute relevant scalar products

do iState=1,nState
  dScal = Zero
  do iRef=1,nCIRef
    indBase = nState*(iState-1)
    indx = indBase+iCIInd(iRef)-1
    dScal = dScal+Work(iSTC+indx)*dCIRef(iRef)
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
