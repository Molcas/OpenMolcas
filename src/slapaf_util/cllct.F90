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
!***********************************************************************

subroutine Cllct(Strng,Vector,value,nAtom,Coor,nCntr,mCntr,xyz,Temp,Ind,type,qMss,TMtrx,Lbl,lWrite,Deg,lAtom)
!***********************************************************************
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             May '91                                                  *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iOper
use Slapaf_Info, only: dMass, AtomLbl

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "Molcas.fh"
character*(*) Strng
character Label*(LenIn5), Name*(LenIn)
character Oper*3, type*6, Lbl*8
real*8 Coor(3,nAtom), Vector(3,nAtom), xyz(3,nCntr+mCntr), Temp(3,nCntr+mCntr), qMss(nCntr+mCntr), TMtrx(3,nAtom,3,(nCntr+mCntr)), &
       Axis(3), Perp_Axis(3,2)
integer Ind(nCntr+mCntr,2)
logical lWrite, ldB, lWarn, lAtom(nAtom)
dimension Dummy(1)

iRout = 50
iPrint = nPrint(iRout)
ldB = .false.
lWarn = lWrite
if (iPrint > 20) lWrite = .true.
if (iPrint >= 99) call RecPrt(' In Cllct: Coor',' ',Coor,3,nAtom)

iFrst = 1
nCent = nCntr+mCntr

! Pick up cartesian coordinates associated with the
! internal coordinate

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

    ! Check if operator belongs to the current point group

    i = 0
    do j=1,nIrrep-1
      if (iPhase == iOper(j)) i = j
    end do
    if (i == 0) then
      call WarningMessage(2,'Error in Cllclt')
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
    call WarningMessage(2,'Error in Cllclt')
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
    call WarningMessage(2,'Error in Cllclt')
    write(6,*) '********** ERROR **********'
    write(6,*) ' Unrecognizable atom label '
    write(6,*) ' ',Name
    write(6,*) '***************************'
    write(6,*)
    call Quit_OnUserError()
  end if
  lAtom(jsAtom) = .false.

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
  Deg = D_Cart(Ind,nIrrep)
else if (type == 'Y     ') then
  value = xyz(2,1)
  call dcopy_(3,[Zero],0,Temp,1)
  Temp(2,1) = One
  if (lWrite) write(6,'(1X,A,A,2X,F10.4,A)') Lbl,' : y-component=',value,'/ bohr'
  Deg = D_Cart(Ind,nIrrep)
else if (type == 'Z     ') then
  value = xyz(3,1)
  call dcopy_(3,[Zero],0,Temp,1)
  Temp(3,1) = One
  if (lWrite) write(6,'(1X,A,A,2X,F10.4,A)') Lbl,' : z-component=',value,'/ bohr'
  Deg = D_Cart(Ind,nIrrep)
else if (type == 'STRTCH') then
  call Strtch(xyz,nCent,value,Temp,lWrite,Lbl,Dummy,ldB)
  Deg = D_Bond(Ind,Ind(1,2),nIrrep)
else if (type == 'LBEND1') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCent,value,Temp,lWrite,Lbl,Dummy,ldB,Axis,Perp_Axis(1,1),.false.)
  Deg = D_Bend(Ind,Ind(1,2),nIrrep)
else if (type == 'LBEND2') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCent,value,Temp,lWrite,Lbl,Dummy,ldB,Axis,Perp_Axis(1,2),.true.)
  Deg = D_Bend(Ind,Ind(1,2),nIrrep)
else if (type == 'BEND  ') then
  call Bend(xyz,nCent,value,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
  Deg = D_Bend(Ind,Ind(1,2),nIrrep)
else if (type == 'TRSN  ') then
  call Trsn(xyz,nCent,value,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
  Deg = D_Trsn(Ind,Ind(1,2),nIrrep)
else if (type == 'OUTOFP') then
  call OutOfP(xyz,nCent,value,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
  Deg = D_Trsn(Ind,Ind(1,2),nIrrep)
else if (type == 'DISSOC') then
  call Dissoc(xyz,nCntr,mCntr,qMss,value,Temp,lWrite,Lbl,Dummy,ldB)
  Deg = One
else
  call WarningMessage(2,'Error in Cllclt')
  write(6,'(A,A)') ' Type declaration is corrupted:',type
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
Deg = sqrt(Deg)

! Symmetrize and normalize
! Observe that the B matrix is defined in both symmetry adapted
! internal coordinates and symmetry adapted cartesian coordinates.

if (iPrint >= 99) call RecPrt(' Transformation vector',' ',Temp,3,nCent)

! Symmetrization of P: 1/g Sum over G of X(G)*GPG
!
! i) project away symmetry breaking displacements
! ii) divide the contribution for each center with
!     the multiplicity, g/u
! iii) sum over unique centers.

call dcopy_(3*nAtom,[Zero],0,Vector,1)
call dcopy_(3*nAtom*3*nCent,[Zero],0,TMtrx,1)
do ixyz=1,nCent
  jsAtom = Ind(ixyz,1)
  iPhase = Ind(ixyz,2)
  tx = One
  ty = One
  tz = One

  ! Step i
  ! Project away nonsymmetric displacements

  ! Restrict loop to the stabilizers of the center.
  do iIrrep=0,nIrrep-1
    if ((Coor(1,jsAtom) /= Zero) .and. (iand(iOper(iIrrep),1) /= 0)) Go To 350
    if ((Coor(2,jsAtom) /= Zero) .and. (iand(iOper(iIrrep),2) /= 0)) Go To 350
    if ((Coor(3,jsAtom) /= Zero) .and. (iand(iOper(iIrrep),4) /= 0)) Go To 350

    if (iand(iOper(iIrrep),1) /= 0) tx = Zero
    if (iand(iOper(iIrrep),2) /= 0) ty = Zero
    if (iand(iOper(iIrrep),4) /= 0) tz = Zero
350 continue
  end do

  ! Step ii

  ! Rotate vector back to the unique center

  if (iand(iPhase,1) /= 0) tx = -tx
  if (iand(iPhase,2) /= 0) ty = -ty
  if (iand(iPhase,4) /= 0) tz = -tz
  TMtrx(1,jsAtom,1,ixyz) = tx
  TMtrx(2,jsAtom,2,ixyz) = ty
  TMtrx(3,jsAtom,3,ixyz) = tz
end do
call dGeMV_('N',3*nAtom,3*nCent,One,TMtrx,3*nAtom,Temp,1,Zero,Vector,1)
if (iPrint >= 99) then
  call RecPrt('TMtrx',' ',TMtrx,3*nAtom,3*nCent)
  call RecPrt(' symmetry adapted vector',' ',Vector,3,nAtom)
end if

return

end subroutine Cllct
