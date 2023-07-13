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
! Copyright (C) 1991, Roland Lindh                                     *
!               2008, Giovanni Ghigo                                   *
!***********************************************************************

subroutine CllCtoF(Strng,nCntr,mCntr,xyz,Temp,Ind,Typ,qMss,Lbl)
!***********************************************************************
!                                                                      *
!     Author: Giovanni Ghigo, Dep. of General and Organic Chemistry    *
!             University of Torino, ITALY                              *
!             July 2008                                                *
!     Adapted from  Cllct  by                                          *
!             Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             May '91                                                  *
!***********************************************************************

use Symmetry_Info, only: iOper, nIrrep
use Slapaf_Info, only: AtomLbl, Coor, dMass
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: Strng
integer(kind=iwp), intent(in) :: nCntr, mCntr
real(kind=wp), intent(out) :: xyz(3,nCntr+mCntr), Temp(3,nCntr+mCntr), qMss(nCntr+mCntr)
integer(kind=iwp), intent(out) :: Ind(nCntr+mCntr,2)
character(len=6), intent(in) :: Typ
character(len=8), intent(in) :: Lbl
#include "print.fh"
#include "Molcas.fh"
integer(kind=iwp) :: i, iEnd, iFrst, iPhase, iPrint, iRout, isAtom, ixyz, j, jsAtom, nAtom, nCent, nPar1, nPar2
real(kind=wp) :: Axis(3), Dummy(1), Perp_Axis(3,2), Val
logical(kind=iwp) :: ldB, lWarn, lWrite
character(len=LenIn1) :: Label
character(len=LenIn) :: AtName
character(len=3) :: Oper

nAtom = size(Coor,2)

iRout = 50
iPrint = nPrint(iRout)
ldB = .false.
lWrite = .true.
lWarn = .true.
if (iPrint >= 99) call RecPrt(' In CllCtoF: Coor',' ',Coor,3,nAtom)

iFrst = 1
nCent = nCntr+mCntr

! Pick up cartesian coordinates associated with the internal coordinate

do ixyz=1,nCent
  call NxtWrd(Strng,iFrst,iEnd)
  Label = Strng(iFrst:iEnd)
  nPar1 = index(Label,'(')
  nPar2 = index(Label,')')
  iPhase = 0
  if ((nPar1 /= 0) .and. (nPar2 /= 0)) then
    AtName = '    '
    AtName = Label(1:nPar1-1)
    Oper = Label(nPar1+1:nPar2-1)
    call UpCase(Oper)
    if (index(Oper,'X') /= 0) iPhase = ieor(iPhase,1)
    if (index(Oper,'Y') /= 0) iPhase = ieor(iPhase,2)
    if (index(Oper,'Z') /= 0) iPhase = ieor(iPhase,4)

    ! Check if operator belong to the current point group

    i = 0
    do j=1,nIrrep-1
      if (iPhase == iOper(j)) i = j
    end do
    if (i == 0) then
      call WarningMessage(2,'Error in CllctOF')
      write(u6,*) '*********** ERROR ***********'
      write(u6,*) ' Undefined symmetry operator '
      write(u6,*) '*****************************'
      write(u6,*) ' ',Oper
      write(u6,*)
      call Quit_OnUserError()
    end if
  else if ((nPar1 == 0) .and. (nPar2 == 0)) then
    AtName = '    '
    AtName = Strng(iFrst:iEnd)
    Oper = ' '
  else
    call WarningMessage(2,'Error in CllctOF')
    write(u6,*) '********** ERROR **********'
    write(u6,*) ' Syntax error in:'
    write(u6,*) ' ',Label
    write(u6,*) '***************************'
    write(u6,*)
    call Quit_OnUserError()
  end if

  ! Find corresponding coordinate

  jsAtom = 0
  do isAtom=1,nAtom
    if (AtName == AtomLbl(isAtom)) jsAtom = isAtom
  end do
  if (jsAtom == 0) then
    call WarningMessage(2,'Error in CllctOF')
    write(u6,*) '********** ERROR **********'
    write(u6,*) ' Unrecognizable atom label '
    write(u6,*) ' ',trim(AtName)
    write(u6,*) '***************************'
    write(u6,*)
    call Quit_OnUserError()
  end if

  ! Store away the unique center index and the operator

  Ind(ixyz,1) = jsAtom
  Ind(ixyz,2) = iPhase
  xyz(:,ixyz) = Coor(:,jsAtom)
  ! Generate actual coordinate
  if (btest(iPhase,0)) xyz(1,ixyz) = -xyz(1,ixyz)
  if (btest(iPhase,1)) xyz(2,ixyz) = -xyz(2,ixyz)
  if (btest(iPhase,2)) xyz(3,ixyz) = -xyz(3,ixyz)
  if (Typ == 'DISSOC') qMss(ixyz) = dMass(jsAtom)
  iFrst = iEnd+1
end do

if (iPrint >= 99) call RecPrt(' Coordinates',' ',xyz,3,nCent)
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the internal coordinate

if (Typ == 'X     ') then
  Val = xyz(1,1)
  Temp(:,:) = Zero
  Temp(1,1) = One
  if (lWrite) write(u6,'(1X,A,A,2X,F10.4,A)') Lbl,' : x-component=',Val,'/ bohr'
else if (Typ == 'Y     ') then
  Val = xyz(2,1)
  Temp(:,:) = Zero
  Temp(2,1) = One
  if (lWrite) write(u6,'(1X,A,A,2X,F10.4,A)') Lbl,' : y-component=',Val,'/ bohr'
else if (Typ == 'Z     ') then
  Val = xyz(3,1)
  Temp(:,:) = Zero
  Temp(3,1) = One
  if (lWrite) write(u6,'(1X,A,A,2X,F10.4,A)') Lbl,' : z-component=',Val,'/ bohr'
else if (Typ == 'STRTCH') then
  call Strtch(xyz,nCent,Val,Temp,lWrite,Lbl,Dummy,ldB)
else if (Typ == 'LBEND1') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCent,Val,Temp,lWrite,Lbl,Dummy,ldB,Axis,Perp_Axis(1,1),.false.)
else if (Typ == 'LBEND2') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCent,Val,Temp,lWrite,Lbl,Dummy,ldB,Axis,Perp_Axis(1,2),.true.)
else if (Typ == 'BEND  ') then
  call Bend(xyz,nCent,Val,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
else if (Typ == 'TRSN  ') then
  call Trsn(xyz,nCent,Val,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
else if (Typ == 'OUTOFP') then
  call OutOfP(xyz,nCent,Val,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
else if (Typ == 'DISSOC') then
  call Dissoc(xyz,nCntr,mCntr,qMss,Val,Temp,lWrite,Lbl,Dummy,ldB)
else
  call WarningMessage(2,'Error in CllctOF')
  write(u6,'(A,A)') ' Type declaration is corrupted: ',trim(Typ)
  call Quit_OnUserError()
end if

return

end subroutine CllCtoF
