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

character*(*) OrbFileName
character*8 InfoLbl
character*80 line
dimension NBAS_L(8)
dimension NORB_L(8)
#include "warnings.h"
#include "inporbfmt.fh"
logical lExists

call f_Inquire(OrbFileName,lExists)
if (.not. lExists) Go To 921
LU = 99
LU = IsFreeUnit(LU)
call Molcas_Open(LU,OrbFileName)
! Check version!
read(LU,'(A80)',ERR=920,end=920) Line
iVer = 0
do jVer=1,mxVer
  if (Magic(jVer) == Line(1:len(Magic(jVer)))) iVer = jVer
end do
iRc = _RC_IO_ERROR_READ_
if (iVer == 0) return

! Find and read information section:
10 continue
read(LU,'(A80)',ERR=920) LINE
if (LINE /= '#INFO') goto 10
11 continue
read(LU,'(A80)',ERR=920) LINE
if (LINE(1:1) == '*') goto 11
read(LINE,*) IUHF,NSYM_L
12 continue
read(LU,'(A80)',ERR=920) LINE
if (LINE(1:1) == '*') goto 12
read(LINE,*) (NBAS_L(I),I=1,NSYM_L)
13 continue
read(LU,'(A80)',ERR=920) LINE
if (LINE(1:1) == '*') goto 13
read(LINE,*) (NORB_L(I),I=1,NSYM_L)

! Create a InfoLbl telling what data sets are available:
iORB = 0
iOCC = 0
iONE = 0
iIND = 0
! Find section InfoLbls:
if (IUHF == 0) then
21 continue
  read(LU,'(A80)',end=900,ERR=900) LINE
  if (LINE(1:4) == '#ORB') iORB = 1
  if (LINE(1:4) == '#OCC') iOCC = 1
  if (LINE(1:4) == '#ONE') iONE = 1
  if (LINE(1:4) == '#IND') iIND = 1
  goto 21
else
22 continue
  read(LU,'(A80)',end=900,ERR=900) LINE
  if (LINE(1:5) == '#UORB') iORB = 1
  if (LINE(1:5) == '#UOCC') iOCC = 1
  if (LINE(1:5) == '#UONE') iONE = 1
  if (LINE(1:4) == '#IND') iIND = 1
  goto 22
end if

! Intentionally reached end of file:
900 continue
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

920 continue
close(LU)
921 continue
iRC = _RC_IO_ERROR_READ_

return

end subroutine ChkVec
