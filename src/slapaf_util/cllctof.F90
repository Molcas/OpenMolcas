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

subroutine CllCtoF(Strng,nCntr,mCntr,xyz,Temp,Ind,type,qMss,Lbl)
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

use Symmetry_Info, only: nIrrep, iOper
use Slapaf_Info, only: dMass, AtomLbl, Coor

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "Molcas.fh"
character*(*) Strng
character Label*(LenIn1), Name*(LenIn), Oper*3, type*6, Lbl*8
real*8 xyz(3,nCntr+mCntr), Temp(3,nCntr+mCntr), qMss(nCntr+mCntr), Axis(3), Perp_Axis(3,2)
integer Ind(nCntr+mCntr,2)
logical lWrite, ldB, lWarn
dimension Dummy(1)

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
    Name = '    '
    Name = Label(1:nPar1-1)
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
      write(6,*) '*********** ERROR ***********'
      write(6,*) ' Undefined symmetry operator '
      write(6,*) '*****************************'
      write(6,*) ' ',Oper
      write(6,*)
      call Quit_OnUserError()
    end if
  else if ((nPar1 == 0) .and. (nPar2 == 0)) then
    Name = '    '
    Name = Strng(iFrst:iEnd)
    Oper = ' '
  else
    call WarningMessage(2,'Error in CllctOF')
    write(6,*) '********** ERROR **********'
    write(6,*) ' Syntax error in:'
    write(6,*) ' ',Label
    write(6,*) '***************************'
    write(6,*)
    call Quit_OnUserError()
  end if

  ! Find corresponding coordinate

  jsAtom = 0
  do isAtom=1,nAtom
    if (Name == AtomLbl(isAtom)) jsAtom = isAtom
  end do
  if (jsAtom == 0) then
    call WarningMessage(2,'Error in CllctOF')
    write(6,*) '********** ERROR **********'
    write(6,*) ' Unrecognizable atom label '
    write(6,*) ' ',Name
    write(6,*) '***************************'
    write(6,*)
    call Quit_OnUserError()
  end if

  ! Store away the unique center index and the operator

  Ind(ixyz,1) = jsAtom
  Ind(ixyz,2) = iPhase
  call dcopy_(3,Coor(1,jsAtom),1,xyz(1,ixyz),1)
  ! Generate actual coordinate
  if (iand(iPhase,1) /= 0) xyz(1,ixyz) = -xyz(1,ixyz)
  if (iand(iPhase,2) /= 0) xyz(2,ixyz) = -xyz(2,ixyz)
  if (iand(iPhase,4) /= 0) xyz(3,ixyz) = -xyz(3,ixyz)
  if (type == 'DISSOC') qMss(ixyz) = dMass(jsAtom)
  iFrst = iEnd+1
end do

if (iPrint >= 99) call RecPrt(' Coordinates',' ',xyz,3,nCent)
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the internal coordinate

if (type == 'X     ') then
  value = xyz(1,1)
  call dcopy_(3,[Zero],0,Temp,1)
  Temp(1,1) = One
  if (lWrite) write(6,'(1X,A,A,2X,F10.4,A)') Lbl,' : x-component=',value,'/ bohr'
else if (type == 'Y     ') then
  value = xyz(2,1)
  call dcopy_(3,[Zero],0,Temp,1)
  Temp(2,1) = One
  if (lWrite) write(6,'(1X,A,A,2X,F10.4,A)') Lbl,' : y-component=',value,'/ bohr'
else if (type == 'Z     ') then
  value = xyz(3,1)
  call dcopy_(3,[Zero],0,Temp,1)
  Temp(3,1) = One
  if (lWrite) write(6,'(1X,A,A,2X,F10.4,A)') Lbl,' : z-component=',value,'/ bohr'
else if (type == 'STRTCH') then
  call Strtch(xyz,nCent,value,Temp,lWrite,Lbl,Dummy,ldB)
else if (type == 'LBEND1') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCent,value,Temp,lWrite,Lbl,Dummy,ldB,Axis,Perp_Axis(1,1),.false.)
else if (type == 'LBEND2') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCent,value,Temp,lWrite,Lbl,Dummy,ldB,Axis,Perp_Axis(1,2),.true.)
else if (type == 'BEND  ') then
  call Bend(xyz,nCent,value,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
else if (type == 'TRSN  ') then
  call Trsn(xyz,nCent,value,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
else if (type == 'OUTOFP') then
  call OutOfP(xyz,nCent,value,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
else if (type == 'DISSOC') then
  call Dissoc(xyz,nCntr,mCntr,qMss,value,Temp,lWrite,Lbl,Dummy,ldB)
else
  call WarningMessage(2,'Error in CllctOF')
  write(6,'(A,A)') ' Type declaration is corrupted:',type
  call Quit_OnUserError()
end if

return

end subroutine CllCtoF
