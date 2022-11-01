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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, uToau
use Definitions, only: wp, iwp

implicit none
#include "Molcas.fh"
real(kind=wp) :: Coor(3,mxAtom), FMass(mxAtom) !IFG
character(len=LenIn) :: FAtLbl(mxAtom) !IFG
integer(kind=iwp) :: nFAtoms
logical(kind=iwp) :: lSlapaf
integer(kind=iwp) :: i, iAt, iOper(8), jAt, jOper, lw2, mCenter, nAtoms, nCenter, nOper, nSym
real(kind=wp) :: AMass, RotVec(3), Xnew, Xold, Xold2, Ynew, Yold, Yold2, Znew, Zold, Zold2
character(len=LenIn) :: Byte4
real(kind=wp), allocatable :: Mass(:), w1(:,:)
character(len=LenIn), allocatable :: AtomLbl(:)

call Get_iScalar('nSym',nSym)
call Get_iArray('Symmetry operations',iOper,nSym)
call Get_iScalar('Unique atoms',nAtoms)
call mma_allocate(AtomLbl,8*nAtoms,label='AtomLbl')
call Get_cArray('Unique Atom Names',AtomLbl,LenIn*nAtoms)
call mma_allocate(w1,3,8*nAtoms,label='w1')
if (lSlapaf) then
  call Get_dArray('Initial Coordinates',w1,3*nAtoms)
else
  call Get_dArray('Unique Coordinates',w1,3*nAtoms)
end if

call mma_allocate(Mass,8*nAtoms,label='Mass')
call Get_Mass(Mass,nAtoms)
call dScal_(nAtoms,One/uToau,Mass,1)

if (nSym == 1) then

  nFAtoms = nAtoms
  do i=1,nFAtoms
    FMass(i) = Mass(i)
    FAtLbl(i) = AtomLbl(i)
    Coor(:,i) = w1(:,i)
  end do

else

  lw2 = 0
  nOper = 0
  if (nSym == 2) nOper = 1
  if (nSym == 4) nOper = 2
  if (nSym == 8) nOper = 3
  nCenter = nAtoms
  do i=1,nOper
    jOper = i+1
    if (i == 3) jOper = 5
    RotVec(1) = One
    if (btest(iOper(jOper),0)) RotVec(1) = -One
    RotVec(2) = One
    if (btest(iOper(jOper),1)) RotVec(2) = -One
    RotVec(3) = One
    if (btest(iOper(jOper),2)) RotVec(3) = -One
    mCenter = nCenter
    outer: do iAt=1,mCenter
      Xold = w1(1,iAt)
      Yold = w1(2,iAt)
      Zold = w1(3,iAt)
      Byte4 = AtomLbl(lw2+iAt)
      AMass = Mass(lw2+iAt)
      FMass(lw2+iAt) = AMass
      Xnew = RotVec(1)*Xold
      Ynew = RotVec(2)*Yold
      Znew = RotVec(3)*Zold
      do jAt=1,nCenter
        if (Byte4 == AtomLbl(lw2+jAt)) then
          Xold2 = w1(1,jAt)
          Yold2 = w1(2,jAt)
          Zold2 = w1(3,jAt)
          if ((Xnew == Xold2) .and. (Ynew == Yold2) .and. (Znew == Zold2)) cycle outer
        end if
      end do
      nCenter = nCenter+1
      w1(1,nCenter) = Xnew
      w1(2,nCenter) = Ynew
      w1(3,nCenter) = Znew
      AtomLbl(lw2+nCenter) = Byte4
      Mass(lw2+nCenter) = AMass
    end do outer
  end do
  nFAtoms = nCenter

  do iAt=1,nCenter
    FAtLbl(iAt) = AtomLbl(lw2+iAt)
    FMass(iAt) = Mass(lw2+iAt)
    Coor(:,iAt) = w1(:,iAt)
  end do

end if

call mma_deallocate(AtomLbl)
call mma_deallocate(w1)
call mma_deallocate(Mass)

return

end subroutine GetFullCoord
