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
real(kind=wp) :: Coor(3,mxAtom), FMass(mxAtom) !IFG
character(len=LenIn) :: FAtLbl(mxAtom) !IFG
integer(kind=iwp) :: nFAtoms
logical(kind=iwp) :: lSlapaf
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iAt, iOper(8), jAt, jOper, lw1, lw2, mCenter, nAtoms, nCenter, nOper, nSym
real(kind=wp) :: Mass(MxAtom), AMass, RotVec(3), Xnew, Xold, Xold2, Ynew, Yold, Yold2, Znew, Zold, Zold2 !IFG
character(len=LenIn) :: AtomLbl(mxAtom), Byte4 !IFG

call Get_iScalar('nSym',nSym)
call Get_iArray('Symmetry operations',iOper,nSym)
call Get_iScalar('Unique atoms',nAtoms)
call Get_cArray('Unique Atom Names',AtomLbl,LenIn*nAtoms)
call GetMem('Coor','ALLO','REAL',lw1,3*8*nAtoms)
if (lSlapaf) then
  call Get_dArray('Initial Coordinates',Work(lw1),3*nAtoms)
else
  call Get_dArray('Unique Coordinates',Work(lw1),3*nAtoms)
end if

call Get_Mass(Mass,nAtoms)
call dScal_(nAtoms,One/uToau,Mass,1)

if (nSym == 1) then

  nFAtoms = nAtoms
  do i=0,nFAtoms-1
    FMass(i+1) = Mass(i+1)
    FAtLbl(i+1) = AtomLbl(i+1)
    Coor(1,i+1) = Work(lw1+3*i)
    Coor(2,i+1) = Work(lw1+3*i+1)
    Coor(3,i+1) = Work(lw1+3*i+2)
  end do

else

  lw2 = 1
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
    outer: do iAt=0,mCenter-1
      Xold = Work(lw1+iAt*3+0)
      Yold = Work(lw1+iAt*3+1)
      Zold = Work(lw1+iAt*3+2)
      Byte4 = AtomLbl(lw2+iAt)
      AMass = Mass(lw2+iAt)
      FMass(lw2+iAt) = AMass
      Xnew = RotVec(1)*Xold
      Ynew = RotVec(2)*Yold
      Znew = RotVec(3)*Zold
      do jAt=0,nCenter-1
        if (Byte4 == AtomLbl(Lw2+jAt)) then
          Xold2 = Work(lw1+jAt*3+0)
          Yold2 = Work(lw1+jAt*3+1)
          Zold2 = Work(lw1+jAt*3+2)
          if ((Xnew == Xold2) .and. (Ynew == Yold2) .and. (Znew == Zold2)) cycle outer
        end if
      end do
      nCenter = nCenter+1
      Work(lw1+(nCenter-1)*3+0) = Xnew
      Work(lw1+(nCenter-1)*3+1) = Ynew
      Work(lw1+(nCenter-1)*3+2) = Znew
      AtomLbl(lw2+nCenter-1) = Byte4
      Mass(lw2+nCenter-1) = AMass
    end do outer
  end do
  nFAtoms = nCenter

  do iAt=0,nCenter-1
    FAtLbl(iAt+1) = AtomLbl(lw2+iAt)
    FMass(iAt+1) = Mass(lw2+iAt)
    Coor(1,iAt+1) = Work(lw1+3*iAt)
    Coor(2,iAt+1) = Work(lw1+3*iAt+1)
    Coor(3,iAt+1) = Work(lw1+3*iAt+2)
  end do

end if

call GetMem('Coor','FREE','REAL',lw1,3*8*nAtoms)

return

end subroutine GetFullCoord
