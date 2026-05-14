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

subroutine CHO_PUTRED1(INFRED,NNBSTRSH,INDRED,INDRSH,ISP2F,MRED,MSYM,MMSHL,LMMBSTRT,IPASS,ILOC)
!
! Purpose: write index arrays for current reduced set (reduced set
!          IPASS).

use Cholesky, only: LuPri, LuRed, MaxRed, nnBstRT, nnShl, nSym
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: MRED, INFRED(MRED), MSYM, MMSHL, LMMBSTRT, IPASS, ILOC
integer(kind=iwp), intent(_IN_) :: NNBSTRSH(MSYM,MMSHL), INDRED(LMMBSTRT), INDRSH(LMMBSTRT), ISP2F(MMSHL)
integer(kind=iwp) :: IADR, IADR1, IOPT, LTOT
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_PUTRED1'

! Test dimensions.
! ----------------

if ((ILOC < 1) .or. (ILOC > 3)) call CHO_QUIT('ILOC error in '//SECNAM,104)

if (MSYM /= NSYM) call CHO_QUIT('NSYM error in '//SECNAM,104)

if (MMSHL /= NNSHL) call CHO_QUIT('NNSHL error in '//SECNAM,104)

if (LMMBSTRT /= NNBSTRT(1)) call CHO_QUIT('NNBSTRT(1) error in '//SECNAM,104)

if (LMMBSTRT < NNBSTRT(ILOC)) call CHO_QUIT('NNBSTRT(ILOC) error in '//SECNAM,104)

if ((IPASS < 1) .or. (IPASS > MAXRED)) call CHO_QUIT('IPASS error in '//SECNAM,104)

! Get first address.
! ------------------

IADR1 = INFRED(IPASS)
if (IADR1 < 0) then
  write(LUPRI,*) SECNAM,': negative address for reduced set ',IPASS,': ',IADR1
  call CHO_QUIT('Error in '//SECNAM,104)
end if

if (LOCDBG) write(LUPRI,*) SECNAM,': putting reduced set ',IPASS,' at addr: ',IADR1

! Write index arrays.
! -------------------

IOPT = 1
IADR = IADR1
LTOT = NSYM*NNSHL
call IDAFILE(LURED,IOPT,NNBSTRSH,LTOT,IADR)
IOPT = 1
IADR = IADR1+NSYM*NNSHL
LTOT = NNBSTRT(ILOC)
call IDAFILE(LURED,IOPT,INDRED,LTOT,IADR)
if (IPASS == 1) then
  IOPT = 1
  IADR = IADR1+NSYM*NNSHL+NNBSTRT(1)
  LTOT = NNBSTRT(1)
  call IDAFILE(LURED,IOPT,INDRSH,LTOT,IADR)
  IOPT = 1
  IADR = IADR1+NSYM*NNSHL+2*NNBSTRT(1)
  LTOT = NNSHL
  call IDAFILE(LURED,IOPT,ISP2F,LTOT,IADR)
end if

end subroutine CHO_PUTRED1
