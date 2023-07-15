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
! Copyright (C) 1991,1997, Roland Lindh                                *
!***********************************************************************

subroutine Cllct2(Strng,Vector,dVector,Val,nAtom,nCntr,mCntr,xyz,Grad,Ind,Typ,qMss,Lbl,lWrite,Deg,Hess,lIter)
!***********************************************************************
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             May '91                                                  *
!                                                                      *
!             Modified to be used in optimizations with constraints,   *
!             June '97 (R. Lindh)                                      *
!***********************************************************************

use Symmetry_Info, only: iOper, nIrrep
use Slapaf_Info, only: AtomLbl, Cx, dMass
use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: Strng
integer(kind=iwp), intent(in) :: nAtom, nCntr, mCntr, lIter
real(kind=wp), intent(out) :: Vector(3,nAtom), dVector(3,nAtom,3,nAtom), Val, xyz(3,nCntr+mCntr), Grad(3,nCntr+mCntr), &
                              qMss(nCntr+mCntr), Deg, Hess(3,nCntr+mCntr,3,nCntr+mCntr)
integer(kind=iwp), intent(out) :: Ind(nCntr+mCntr,2)
character(len=6), intent(in) :: Typ
character(len=8), intent(in) :: Lbl
logical(kind=iwp), intent(inout) :: lWrite
#include "print.fh"
#include "Molcas.fh"
integer(kind=iwp) :: i, iEnd, iFrst, iPhase, iPrint, iRout, isAtom, ixyz, j, jsAtom, lStrng, nCent, nPar1, nPar2
real(kind=wp) :: Axis(3), Perp_Axis(3,2)
character(len=LenIn5) :: Label
character(len=LenIn) :: AtName
character(len=3) :: Oper
logical(kind=iwp) :: ldB, lWarn
integer(kind=iwp), allocatable :: iDCR(:)
real(kind=wp), external :: D_Bend, D_Bond, D_Cart, D_Trsn

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 50
iPrint = nPrint(iRout)
ldB = .true.
lWarn = lWrite
if (iPrint > 20) lWrite = .true.
#ifdef _DEBUGPRINT_
call RecPrt(' In Cllct2: Coor',' ',Cx(:,:,lIter),3,nAtom)
call RecPrt('dMass',' ',dMass,1,nAtom)
#endif

iFrst = 1
iEnd = 1
lStrng = len(Strng)

! Pick up cartesian coordinates associated with the internal coordinate

call mma_allocate(iDCR,MxAtom,Label='iDCR')

nCent = nCntr+mCntr
do ixyz=1,nCent
  call NxtWrd(Strng,iFrst,iEnd)
  if (iEnd >= iFrst) then
    Label = Strng(iFrst:iEnd)
    nPar1 = index(Label,'(')
    nPar2 = index(Label,')')
  else
    Label = ' '
    nPar1 = 0
    nPar2 = 0
  end if
  iPhase = 0
  if ((nPar1 /= 0) .and. (nPar2 /= 0)) then
    AtName = '    '
    AtName = Label(1:nPar1-1)
    Oper = Label(nPar1+1:nPar2-1)
    call UpCase(Oper)
    if (index(Oper,'X') /= 0) iPhase = ieor(iPhase,1)
    if (index(Oper,'Y') /= 0) iPhase = ieor(iPhase,2)
    if (index(Oper,'Z') /= 0) iPhase = ieor(iPhase,4)

    ! Check if operator belongs to the current point group

    i = 0
    do j=1,nIrrep-1
      if (iPhase == iOper(j)) i = j
    end do
    iDCR(ixyz) = iOper(i)
    if (i == 0) then
      call WarningMessage(2,' Undefined symmetry operator')
      write(u6,'(A)') Oper
      call Quit_OnUserError()
    end if
    iFrst = iEnd+1
  else if ((nPar1 == 0) .and. (nPar2 == 0)) then
    AtName = '    '
    Oper = ' '
    if (iEnd >= iFrst) then
      AtName = Strng(iFrst:iEnd)
      if ((iEnd >= 1) .and. (iEnd < lStrng)) iFrst = iEnd+1
    end if
    iDCR(ixyz) = iOper(0)
  else
    call WarningMessage(2,' Syntax error in:'//Label)
    call Quit_OnUserError()
  end if

  ! Find corresponding coordinate

  if ((Typ(1:5) /= 'EDIFF') .and. (Typ(1:3) /= 'NAC') .and. (Typ(1:6) /= 'SPHERE') .and. (Typ(1:6) /= 'TRANSV')) then
    jsAtom = 0
    do isAtom=1,nAtom
      if (AtName == AtomLbl(isAtom)) jsAtom = isAtom
    end do
    if (jsAtom == 0) then
      call WarningMessage(2,' Unrecognizable atom label '//trim(AtName))
      call Quit_OnUserError()
    end if
  else
    jsAtom = ixyz
  end if

  ! Store away the unique center index and the operator

  Ind(ixyz,1) = jsAtom
  Ind(ixyz,2) = iPhase
  xyz(:,ixyz) = Cx(:,jsAtom,lIter)
  ! Generate actual coordinate
  if (btest(iPhase,0)) xyz(1,ixyz) = -xyz(1,ixyz)
  if (btest(iPhase,1)) xyz(2,ixyz) = -xyz(2,ixyz)
  if (btest(iPhase,2)) xyz(3,ixyz) = -xyz(3,ixyz)
  if (Typ == 'DISSOC') qMss(ixyz) = dMass(jsAtom)

end do  ! do ixyz=1,nCntr+mCntr

if (iPrint >= 99) then
  call RecPrt(' Coordinates',' ',xyz,3,nCntr+mCntr)
  call RecPrt('qMss',' ',qMss,1,nCntr+mCntr)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the internal coordinate

if (Typ == 'X     ') then
  Val = xyz(1,1)
  Grad(:,:) = Zero
  Hess(:,:,:,:) = Zero
  Grad(1,1) = One
  if (lWrite) write(u6,'(1X,A,A,2X,F10.4,A)') Lbl,' : x-component=',Val,'/ bohr'
  Deg = D_Cart(Ind(1,1),nIrrep)
else if (Typ == 'Y     ') then
  Val = xyz(2,1)
  Grad(:,:) = Zero
  Hess(:,:,:,:) = Zero
  Grad(2,1) = One
  if (lWrite) write(u6,'(1X,A,A,2X,F10.4,A)') Lbl,' : y-component=',Val,'/ bohr'
  Deg = D_Cart(Ind(1,1),nIrrep)
else if (Typ == 'Z     ') then
  Val = xyz(3,1)
  Grad(:,:) = Zero
  Hess(:,:,:,:) = Zero
  Grad(3,1) = One
  if (lWrite) write(u6,'(1X,A,A,2X,F10.4,A)') Lbl,' : z-component=',Val,'/ bohr'
  Deg = D_Cart(Ind(1,1),nIrrep)
else if (Typ == 'STRTCH') then
  call Strtch(xyz,nCntr,Val,Grad,lWrite,Lbl,Hess,ldB)
  Deg = D_Bond(Ind(1:2,1),Ind(1:2,2),nIrrep)
else if (Typ == 'LBEND1') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCntr,Val,Grad,lWrite,Lbl,Hess,ldB,Axis,Perp_Axis(1,1),.false.)
  Deg = D_Bend(Ind(1:3,1),Ind(1:3,2),nIrrep)
else if (Typ == 'LBEND2') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCntr,Val,Grad,lWrite,Lbl,Hess,ldB,Axis,Perp_Axis(1,2),.true.)
  Deg = D_Bend(Ind(1:3,1),Ind(1:3,2),nIrrep)
else if (Typ == 'BEND  ') then
  call Bend(xyz,nCntr,Val,Grad,lWrite,lWarn,Lbl,Hess,ldB)
  Deg = D_Bend(Ind(1:3,1),Ind(1:3,2),nIrrep)
else if (Typ == 'TRSN  ') then
  call Trsn(xyz,nCntr,Val,Grad,lWrite,lWarn,Lbl,Hess,ldB)
  Deg = D_Trsn(Ind(1:4,1),Ind(1:4,2),nIrrep)
else if (Typ == 'OUTOFP') then
  call OutOfP(xyz,nCntr,Val,Grad,lWrite,lWarn,Lbl,Hess,ldB)
  Deg = D_Trsn(Ind(1:4,1),Ind(1:4,2),nIrrep)
else if (Typ(1:3) == 'NAC') then
  call NACInt(xyz,nCntr,Val,Grad,lWrite,Lbl,Hess,ldB,lIter)
  Deg = One
else if (Typ(1:5) == 'EDIFF') then
  call ConInt(xyz,nCntr,Val,Grad,lWrite,Lbl,Hess,ldB,lIter)
  Deg = One
else if (Typ(1:6) == 'SPHERE') then
  call SphInt(xyz,nCntr,xyz,.false.,Val,Grad,lWrite,Lbl,Hess,ldB)
  Deg = One
else if (Typ(1:6) == 'TRANSV') then
  call Transverse(xyz,nCntr,Val,Grad,lWrite,Lbl,Hess,ldB)
  Deg = One
else if (Typ == 'DISSOC') then
  call Dissoc(xyz,nCntr,mCntr,qMss,Val,Grad,lWrite,Lbl,Hess,ldB)
  Deg = One
else
  call WarningMessage(2,' Type declaration is not supported: '//trim(Typ))
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
Deg = sqrt(Deg)

call ProjSym2(nAtom,nCent,Ind,xyz,iDCR,Grad,Vector,Hess,dVector)
if (iPrint >= 99) then
  call RecPrt(' symmetry adapted vector',' ',Vector,3,nAtom)
  call RecPrt(' symmetry adapted dvector',' ',dVector,3*nAtom,3*nAtom)
end if

call mma_deallocate(iDCR)

return

end subroutine Cllct2
