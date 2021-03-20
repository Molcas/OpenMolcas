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

subroutine MyCoor(isAuto,Ox,Oy,Oz,Rx,Ry,Rz,iGx,iGy,iGz,iMagic,iCustOrig)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
!                                                                      *
!   Read Coordinates and calculate a cub for grid                      *
!   isAuto=1 - real job, else only print                               *
!   Origin(3) fix cub in space                                         *
!   Rx,Ry,Rz - size of cub                                             *
!   iGx,iGy,iGz - net                                                  *
!   iMagic = magic guess for net                                       *
!***********************************************************************

use Constants, only: One, Two, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: isAuto, iMagic, iCustOrig
real(kind=wp), intent(inout) :: Ox, Oy, Oz, Rx, Ry, Rz
integer(kind=iwp), intent(inout) :: iGx, iGy, iGz
#include "Molcas.fh"
#include "grid.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iAt, lw2, nCenter, NoOrig
real(kind=wp) :: rrx, rry, rrz, x529
character(len=128) :: Line
character(len=LenIn) :: Byte4

!----------------------------------------------------------------------*
!     Prologue                                                         *
!----------------------------------------------------------------------*
if (isLUSCUS == 1) then
  x529 = Angstrom
else
  x529 = One
end if
lw2 = 1
call get_iScalar('nSym',nSym)
call Get_nAtoms_All(nAtoms)
call Get_LblCnt_All(AtomLbl)
call GetMem('Coor','ALLO','REAL',ipCoor,3*nAtoms)
call Get_Coord_All(Work(ipCoor),nAtoms)
nCenter = nAtoms
! Check : are any strange names?
NoOrig = 0
do iAt=0,nCenter-1
  write(line,'(a)') AtomLbl(lw2+iAt)
  if (index(line,'Ori') /= 0) NoOrig = NoOrig+1
end do

if (ISLUSCUS == 1 .and. ISBINARY == 3) then
  write(LINE,'(2X,I8)') NCENTER-NOORIG
  call PRINTLINE(LID,LINE,10,0)
  call PRINTLINE(LID,LINE,0,0)
  if (isUHF == 1) then
    call PRINTLINE(LID_ab,LINE,10,0)
    call PRINTLINE(LID_ab,LINE,0,0)
  end if
  do IAT=0,NCENTER-1
    write(LINE,'(A)') ATOMLBL(LW2+IAT)
    if (index(LINE,'ORI') == 0) then
      Byte4 = ATOMLBL(LW2+IAT)(1:2)
      if (index('0123456789',Byte4(2:2)) /= 0) Byte4(2:2) = ' '
      write(LINE,'(1X,A2,2X,3F15.8)') Byte4,WORK(IPCOOR+3*IAT)*x529,WORK(IPCOOR+3*IAT+1)*x529,WORK(IPCOOR+3*IAT+2)*x529
      call PRINTLINE(LID,LINE,50,0)
      if (isUHF == 1) then
        call PRINTLINE(LID_ab,LINE,50,0)
      end if
    end if
  end do
end if

write(Line,'(A,I8)') 'Natom= ',nCenter-NoOrig
if (isBinary == 1) then
  write(LuVal) trim(Line)
  if (isUHF == 1) write(LuVal_ab) Line(1:15)
else
  write(LuVal,'(A)') trim(Line)
  if (isUHF == 1) write(LuVal_ab,'(A)') Line(1:15)
end if
do iAt=0,nCenter-1
  if (index(Line,'Ori') == 0) then
    write(line,'(A,2X,3F15.8)') AtomLbl(lw2+iAt),Work(ipCoor+3*iAt),Work(ipCoor+3*iAt+1),Work(ipCoor+3*iAt+2)
  else
    write(Line,'(A)') AtomLbl(lw2+iAt)
  end if
  if (isBinary == 1) then
    write(LuVal) trim(Line)
    if (isUHF == 1) write(LuVal_ab) trim(Line)
  else
    write(LuVal,'(A)') trim(Line)
    if (isUHF == 1) write(LuVal_ab,'(A)') trim(Line)
  end if
end do

if (iCustOrig == 1) goto 500

if (isAuto == 1) then
  !--------------------------------------------------------------------*
  !   Find Cub parameters                                              *
  !               Ox->RxMin, Rx->RxMax                                 *
  !--------------------------------------------------------------------*
  Ox = huge(Ox)
  Oy = huge(Oy)
  Oz = huge(Oz)
  Rx = -huge(Rx)
  Ry = -huge(Ry)
  Rz = -huge(Rz)

  do iAt=0,nCenter-1
    rrx = Work(ipCoor+3*iAt)
    rry = Work(ipCoor+3*iAt+1)
    rrz = Work(ipCoor+3*iAt+2)
    if (rrx < Ox) Ox = rrx
    if (rrx > Rx) Rx = rrx
    if (rry < Oy) Oy = rry
    if (rry > Ry) Ry = rry
    if (rrz < Oz) Oz = rrz
    if (rrz > Rz) Rz = rrz
  end do
  Rx = Rx-Ox
  Ry = Ry-Oy
  Rz = Rz-Oz

  ! and now, expand this cub to place all atoms inside

  Ox = int(Ox-TheGap)
  Oy = int(Oy-TheGap)
  Oz = int(Oz-TheGap)
  Rx = int(Rx+Two*TheGap)
  Ry = int(Ry+Two*TheGap)
  Rz = int(Rz+Two*TheGap)
end if ! finish of iAuto.
!----------------------------------------------------------------------*
!     Calculate corrected coords                                       *
!----------------------------------------------------------------------*

! make a stupid Patch: Cerius2 works well only with even nets!

Rx = 2*(int(Rx)/2)
Ry = 2*(int(Ry)/2)
Rz = 2*(int(Rz)/2)

if (iMagic > 0) then
  iGx = int(abs(Rx))*iMagic
  iGy = int(abs(Ry))*iMagic
  iGz = int(abs(Rz))*iMagic
end if

! make a stupid Patch: Cerius2 works well only with even nets!

iGx = (iGx+1)/2*2
iGy = (iGy+1)/2*2
iGz = (iGz+1)/2*2
mynCenter = nCenter

!----------------------------------------------------------------------*
!     Print coordinates of the system                                  *
!----------------------------------------------------------------------*
!call bXML('Coord')
!call iXML('nCoord',nCenter)
!do iAt=0,nCenter-1
!  call cXML('Atom',AtomLbl(lw2+iAt))
!  call daXML('Atom coord',Work(ipCoor+3*iAt),3)
!end do
!call eXML('Coord')

500 continue
write(u6,*)
write(u6,'(6X,A)') 'Cartesian coordinates:'
write(u6,'(6X,A)') '-----------------------------------------'
write(u6,'(6X,A)') 'No.  Label     X         Y         Z     '
write(u6,'(6X,A)') '-----------------------------------------'
do iAt=0,nCenter-1
  write(u6,'(4X,I4,3X,A,2X,3F10.5)') iAt+1,AtomLbl(lw2+iAt),Work(ipCoor+3*iAt),Work(ipCoor+3*iAt+1),Work(ipCoor+3*iAt+2)
end do
write(u6,'(6X,A)') '-----------------------------------------'
write(u6,'(6X,A,3F12.6)') 'Grid Origin      = ',Ox,Oy,Oz
write(u6,'(6X,A,3F12.6)') 'Grid Axis Length = ',Rx,Ry,Rz
write(u6,'(6X,A)') '-----------------------------------------'
write(u6,*)
write(u6,*)

!----------------------------------------------------------------------*
!     Normal exit                                                      *
!----------------------------------------------------------------------*
return

end subroutine MyCoor
