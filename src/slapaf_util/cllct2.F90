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

subroutine Cllct2(Strng,Vector,dVector,value,nAtom,nCntr,mCntr,xyz,Grad,Ind,type,qMss,Lbl,lWrite,Deg,Hess,lIter)
!***********************************************************************
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             May '91                                                  *
!                                                                      *
!             Modified to be used in optimizations with constraints,   *
!             June '97 (R. Lindh)                                      *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iOper
use Slapaf_Info, only: Cx, dMass, AtomLbl

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "Molcas.fh"
character(len=*) Strng
character(len=LenIn5) Label
character(len=LenIn) Name
character(len=8) Lbl
character Oper*3, type*6
real*8 Vector(3,nAtom), xyz(3,nCntr+mCntr), Grad(3,nCntr+mCntr), dVector(3,nAtom,3,nAtom), Axis(3), Perp_Axis(3,2), &
       qMss(nCntr+mCntr), Hess(3,nCntr+mCntr,3,nCntr+mCntr)
integer Ind(nCntr+mCntr,2), iDCR(MxAtom)
logical lWrite, ldB, lWarn
real*8, allocatable :: Not_Allocated(:,:)
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine SphInt(xyz,nCent,OfRef,RR0,Bf,l_Write,Label,dBf,ldB)
    integer nCent
    real*8 xyz(3,nCent)
    real*8, allocatable, target :: OfRef(:,:)
    real*8 RR0
    real*8 Bf(3,nCent)
    logical l_Write
    character(len=8) Label
    real*8 dBf(3,nCent,3,nCent)
    logical ldB
  end subroutine SphInt
end interface

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

! Pick up cartesian coordinates associated with the
! internal coordinate

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
    Name = '    '
    Name = Label(1:nPar1-1)
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
      write(6,'(A)') Oper
      call Quit_OnUserError()
    end if
    iFrst = iEnd+1
  else if ((nPar1 == 0) .and. (nPar2 == 0)) then
    Name = '    '
    Oper = ' '
    if (iEnd >= iFrst) then
      Name = Strng(iFrst:iEnd)
      if ((iEnd >= 1) .and. (iEnd < lStrng)) iFrst = iEnd+1
    end if
    iDCR(ixyz) = iOper(0)
  else
    call WarningMessage(2,' Syntax error in:'//Label)
    call Quit_OnUserError()
  end if

  ! Find corresponding coordinate

  if ((type(1:5) /= 'EDIFF') .and. (type(1:3) /= 'NAC') .and. (type(1:6) /= 'SPHERE') .and. (type(1:6) /= 'TRANSV')) then
    jsAtom = 0
    do isAtom=1,nAtom
      if (Name == AtomLbl(isAtom)) jsAtom = isAtom
    end do
    if (jsAtom == 0) then
      call WarningMessage(2,' Unrecognizable atom label '//Name)
      call Quit_OnUserError()
    end if
  else
    jsAtom = ixyz
  end if

  ! Store away the unique center index and the operator

  Ind(ixyz,1) = jsAtom
  Ind(ixyz,2) = iPhase
  call dcopy_(3,Cx(:,jsAtom,lIter),1,xyz(1,ixyz),1)
  ! Generate actual coordinate
  if (iand(iPhase,1) /= 0) xyz(1,ixyz) = -xyz(1,ixyz)
  if (iand(iPhase,2) /= 0) xyz(2,ixyz) = -xyz(2,ixyz)
  if (iand(iPhase,4) /= 0) xyz(3,ixyz) = -xyz(3,ixyz)
  if (type == 'DISSOC') qMss(ixyz) = dMass(jsAtom)

end do  ! do ixyz=1,nCntr+mCntr

if (iPrint >= 99) then
  call RecPrt(' Coordinates',' ',xyz,3,nCntr+mCntr)
  call RecPrt('qMss',' ',qMss,1,nCntr+mCntr)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the internal coordinate

if (type == 'X     ') then
  value = xyz(1,1)
  call dcopy_(3,[Zero],0,Grad,1)
  call dcopy_(9,[Zero],0,Hess,1)
  Grad(1,1) = One
  if (lWrite) write(6,'(1X,A,A,2X,F10.4,A)') Lbl,' : x-component=',value,'/ bohr'
  Deg = D_Cart(Ind(1,1),nIrrep)
else if (type == 'Y     ') then
  value = xyz(2,1)
  call dcopy_(3,[Zero],0,Grad,1)
  call dcopy_(9,[Zero],0,Hess,1)
  Grad(2,1) = One
  if (lWrite) write(6,'(1X,A,A,2X,F10.4,A)') Lbl,' : y-component=',value,'/ bohr'
  Deg = D_Cart(Ind(1,1),nIrrep)
else if (type == 'Z     ') then
  value = xyz(3,1)
  call dcopy_(3,[Zero],0,Grad,1)
  call dcopy_(9,[Zero],0,Hess,1)
  Grad(3,1) = One
  if (lWrite) write(6,'(1X,A,A,2X,F10.4,A)') Lbl,' : z-component=',value,'/ bohr'
  Deg = D_Cart(Ind(1,1),nIrrep)
else if (type == 'STRTCH') then
  call Strtch(xyz,nCntr,value,Grad,lWrite,Lbl,Hess,ldB)
  Deg = D_Bond(Ind,Ind(1,2),nIrrep)
else if (type == 'LBEND1') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCntr,value,Grad,lWrite,Lbl,Hess,ldB,Axis,Perp_Axis(1,1),.false.)
  Deg = D_Bend(Ind,Ind(1,2),nIrrep)
else if (type == 'LBEND2') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCntr,value,Grad,lWrite,Lbl,Hess,ldB,Axis,Perp_Axis(1,2),.true.)
  Deg = D_Bend(Ind,Ind(1,2),nIrrep)
else if (type == 'BEND  ') then
  call Bend(xyz,nCntr,value,Grad,lWrite,lWarn,Lbl,Hess,ldB)
  Deg = D_Bend(Ind,Ind(1,2),nIrrep)
else if (type == 'TRSN  ') then
  call Trsn(xyz,nCntr,value,Grad,lWrite,lWarn,Lbl,Hess,ldB)
  Deg = D_Trsn(Ind,Ind(1,2),nIrrep)
else if (type == 'OUTOFP') then
  call OutOfP(xyz,nCntr,value,Grad,lWrite,lWarn,Lbl,Hess,ldB)
  Deg = D_Trsn(Ind,Ind(1,2),nIrrep)
else if (type(1:3) == 'NAC') then
  call NACInt(xyz,nCntr,value,Grad,lWrite,Lbl,Hess,ldB,lIter)
  Deg = One
else if (type(1:5) == 'EDIFF') then
  call ConInt(xyz,nCntr,value,Grad,lWrite,Lbl,Hess,ldB,lIter)
  Deg = One
else if (type(1:6) == 'SPHERE') then
  call SphInt(xyz,nCntr,Not_Allocated,value,Grad,lWrite,Lbl,Hess,ldB)
  Deg = One
else if (type(1:6) == 'TRANSV') then
  call Transverse(xyz,nCntr,value,Grad,lWrite,Lbl,Hess,ldB)
  Deg = One
else if (type == 'DISSOC') then
  call Dissoc(xyz,nCntr,mCntr,qMss,value,Grad,lWrite,Lbl,Hess,ldB)
  Deg = One
else
  call WarningMessage(2,' Type declaration is not supported: '//type)
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

return

end subroutine Cllct2
