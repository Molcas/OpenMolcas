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

subroutine ChkVec(OrbFileName,iVer,NSYM_L,NBAS_L,NORB_L,InfoLbl,iRc)
! Purpose: To check if OrbFileName is a valid orbital file, and which
! information it contains.

use InpOrbFmt, only: Magic, mxVer
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: OrbFileName
integer(kind=iwp), intent(out) :: iVer, NSYM_L, NBAS_L(8), NORB_L(8), iRc
character(len=8), intent(out) :: InfoLbl
integer(kind=iwp) :: I, iIND, iOCC, iONE, iORB, ip, istatus, IUHF, jVer, LU
logical(kind=iwp) :: lExists
character(len=80) :: line
integer(kind=iwp), external :: isFreeUnit
#include "warnings.h"

call f_Inquire(OrbFileName,lExists)
if (.not. lExists) then
  iRC = _RC_IO_ERROR_READ_
  return
end if
LU = 99
LU = IsFreeUnit(LU)
call Molcas_Open(LU,OrbFileName)
! Check version!
read(LU,'(A80)',iostat=istatus) Line
if (istatus /= 0) then
  call Error()
  return
end if
iVer = 0
do jVer=1,mxVer
  if (Magic(jVer) == Line(1:len(Magic(jVer)))) iVer = jVer
end do
iRc = _RC_IO_ERROR_READ_
if (iVer == 0) return

! Find and read information section:
do
  read(LU,'(A80)',iostat=istatus) LINE
  if (istatus > 0) then
    call Error()
    return
  end if
  if (LINE == '#INFO') exit
end do
do
  read(LU,'(A80)',iostat=istatus) LINE
  if (istatus > 0) then
    call Error()
    return
  end if
  if (LINE(1:1) /= '*') exit
end do
read(LINE,*) IUHF,NSYM_L
do
  read(LU,'(A80)',iostat=istatus) LINE
  if (istatus > 0) then
    call Error()
    return
  end if
  if (LINE(1:1) /= '*') exit
end do
read(LINE,*) (NBAS_L(I),I=1,NSYM_L)
do
  read(LU,'(A80)',iostat=istatus) LINE
  if (istatus > 0) then
    call Error()
    return
  end if
  if (LINE(1:1) /= '*') exit
end do
read(LINE,*) (NORB_L(I),I=1,NSYM_L)

! Create a InfoLbl telling what data sets are available:
iORB = 0
iOCC = 0
iONE = 0
iIND = 0
! Find section InfoLbls:
if (IUHF == 0) then
  do
    read(LU,'(A80)',iostat=istatus) LINE
    if (istatus /= 0) exit
    if (LINE(1:4) == '#ORB') iORB = 1
    if (LINE(1:4) == '#OCC') iOCC = 1
    if (LINE(1:4) == '#ONE') iONE = 1
    if (LINE(1:4) == '#IND') iIND = 1
  end do
else
  do
    read(LU,'(A80)',iostat=istatus) LINE
    if (istatus /= 0) exit
    if (LINE(1:5) == '#UORB') iORB = 1
    if (LINE(1:5) == '#UOCC') iOCC = 1
    if (LINE(1:5) == '#UONE') iONE = 1
    if (LINE(1:4) == '#IND') iIND = 1
  end do
end if

! Intentionally reached end of file:
InfoLbl = ' '
ip = 0
ip = ip+iorb
if (iorb == 1) InfoLbl(ip:ip) = 'C'
ip = ip+iocc
if (iocc == 1) InfoLbl(ip:ip) = 'O'
ip = ip+ione
if (ione == 1) InfoLbl(ip:ip) = 'E'
ip = ip+iind
if (iind == 1) InfoLbl(ip:ip) = 'I'

close(LU)
iRC = _RC_ALL_IS_WELL_

return

contains

subroutine Error()
  close(LU)
  iRC = _RC_IO_ERROR_READ_
end subroutine Error

end subroutine ChkVec
