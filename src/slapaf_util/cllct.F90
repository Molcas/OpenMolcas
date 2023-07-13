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

subroutine Cllct(Strng,Vector,Val,nAtom,Coor,nCntr,mCntr,xyz,Temp,Ind,Typ,qMss,TMtrx,Lbl,lWrite,Deg,lAtom)
!***********************************************************************
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             May '91                                                  *
!***********************************************************************

use Symmetry_Info, only: iOper, nIrrep
use Slapaf_Info, only: AtomLbl, dMass
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: Strng
integer(kind=iwp), intent(in) :: nAtom, nCntr, mCntr
real(kind=wp), intent(out) :: Vector(3,nAtom), Val, xyz(3,nCntr+mCntr), Temp(3,nCntr+mCntr), qMss(nCntr+mCntr), &
                              TMtrx(3,nAtom,3,(nCntr+mCntr)), Deg
real(kind=wp), intent(in) :: Coor(3,nAtom)
integer(kind=iwp), intent(out) :: Ind(nCntr+mCntr,2)
character(len=6), intent(in) :: Typ
character(len=8), intent(in) :: Lbl
logical(kind=iwp), intent(inout) :: lWrite
logical(kind=iwp), intent(out) :: lAtom(nAtom)
#include "print.fh"
#include "Molcas.fh"
integer(kind=iwp) :: i, iEnd, iFrst, iIrrep, iPhase, iPrint, iRout, isAtom, ixyz, j, jsAtom, nCent, nPar1, nPar2
real(kind=wp) :: Axis(3), Dummy(1), Perp_Axis(3,2), tx, ty, tz
logical(kind=iwp) :: ldB, lWarn
character(len=LenIn5) :: Label
character(len=LenIn) :: AtName
character(len=3) :: Oper
real(kind=wp), external :: D_Bend, D_Bond, D_Cart, D_Trsn

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
    if (i == 0) then
      call WarningMessage(2,'Error in Cllclt')
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
    call WarningMessage(2,'Error in Cllclt')
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
    call WarningMessage(2,'Error in Cllclt')
    write(u6,*) '********** ERROR **********'
    write(u6,*) ' Unrecognizable atom label '
    write(u6,*) ' ',trim(AtName)
    write(u6,*) '***************************'
    write(u6,*)
    call Quit_OnUserError()
  end if
  lAtom(jsAtom) = .false.

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
  Deg = D_Cart(Ind(1,1),nIrrep)
else if (Typ == 'Y     ') then
  Val = xyz(2,1)
  Temp(:,:) = Zero
  Temp(2,1) = One
  if (lWrite) write(u6,'(1X,A,A,2X,F10.4,A)') Lbl,' : y-component=',Val,'/ bohr'
  Deg = D_Cart(Ind(1,1),nIrrep)
else if (Typ == 'Z     ') then
  Val = xyz(3,1)
  Temp(:,:) = Zero
  Temp(3,1) = One
  if (lWrite) write(u6,'(1X,A,A,2X,F10.4,A)') Lbl,' : z-component=',Val,'/ bohr'
  Deg = D_Cart(Ind(1,1),nIrrep)
else if (Typ == 'STRTCH') then
  call Strtch(xyz,nCent,Val,Temp,lWrite,Lbl,Dummy,ldB)
  Deg = D_Bond(Ind(1:2,1),Ind(1:2,2),nIrrep)
else if (Typ == 'LBEND1') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCent,Val,Temp,lWrite,Lbl,Dummy,ldB,Axis,Perp_Axis(1,1),.false.)
  Deg = D_Bend(Ind(1:3,1),Ind(1:3,2),nIrrep)
else if (Typ == 'LBEND2') then
  call CoSys(xyz,Axis,Perp_Axis)
  call LBend(xyz,nCent,Val,Temp,lWrite,Lbl,Dummy,ldB,Axis,Perp_Axis(1,2),.true.)
  Deg = D_Bend(Ind(1:3,1),Ind(1:3,2),nIrrep)
else if (Typ == 'BEND  ') then
  call Bend(xyz,nCent,Val,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
  Deg = D_Bend(Ind(1:3,1),Ind(1:3,2),nIrrep)
else if (Typ == 'TRSN  ') then
  call Trsn(xyz,nCent,Val,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
  Deg = D_Trsn(Ind(1:4,1),Ind(1:4,2),nIrrep)
else if (Typ == 'OUTOFP') then
  call OutOfP(xyz,nCent,Val,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
  Deg = D_Trsn(Ind(1:4,1),Ind(1:4,2),nIrrep)
else if (Typ == 'DISSOC') then
  call Dissoc(xyz,nCntr,mCntr,qMss,Val,Temp,lWrite,Lbl,Dummy,ldB)
  Deg = One
else
  call WarningMessage(2,'Error in Cllclt')
  write(u6,'(A,A)') ' Type declaration is corrupted: ',trim(Typ)
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

Vector(:,:) = Zero
TMtrx(:,:,:,:) = Zero
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
    if ((Coor(1,jsAtom) /= Zero) .and. btest(iOper(iIrrep),0)) cycle
    if ((Coor(2,jsAtom) /= Zero) .and. btest(iOper(iIrrep),1)) cycle
    if ((Coor(3,jsAtom) /= Zero) .and. btest(iOper(iIrrep),2)) cycle

    if (btest(iOper(iIrrep),0)) tx = Zero
    if (btest(iOper(iIrrep),1)) ty = Zero
    if (btest(iOper(iIrrep),2)) tz = Zero
  end do

  ! Step ii

  ! Rotate vector back to the unique center

  if (btest(iPhase,0)) tx = -tx
  if (btest(iPhase,1)) ty = -ty
  if (btest(iPhase,2)) tz = -tz
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
