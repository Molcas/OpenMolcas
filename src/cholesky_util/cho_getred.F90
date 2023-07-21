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

subroutine CHO_GETRED(IPASS,ILOC,LRSH)
!
! Purpose: read index arrays for current reduced set (reduced set
!          IPASS).

use ChoArr, only: iSP2F
use ChoSwp, only: nnBstRsh, InfRed, IndRSh, IndRed

implicit real*8(a-h,o-z)
integer IPASS, ILOC
logical LRSH
#include "cholesky.fh"
character*10 SECNAM
parameter(SECNAM='CHO_GETRED')
logical LOCDBG
parameter(LOCDBG=.false.)
integer CHO_ISUMELM
external CHO_ISUMELM

#ifdef _DEBUGPRINT_
! Test dimensions.
! ----------------

if (size(nnBstRSh,1) /= NSYM) call CHO_QUIT('NSYM error in '//SECNAM,104)

if (size(nnBstRsh,2) /= NNSHL) call CHO_QUIT('NNSHL error in '//SECNAM,104)

if (size(IndRed,1) /= NNBSTRT(1)) call CHO_QUIT('NNBSTRT(1) error in '//SECNAM,104)

if ((IPASS < 1) .or. (IPASS > size(InfRed)) call CHO_QUIT('IPASS error in '//SECNAM,104)
#endif

! Get first address.
! ------------------

IADR1 = INFRED(IPASS)
#ifdef _DEBUGPRINT_
if (IADR1 < 0) then
  write(LUPRI,*) SECNAM,': negative address for reduced set ',IPASS,': ',IADR1
  call CHO_QUIT('Error in '//SECNAM,104)
end if
#endif

if (LOCDBG) write(LUPRI,*) SECNAM,': getting reduced set ',IPASS,' at addr: ',IADR1

! Read index arrays.
! ------------------

IOPT = 2
IADR = IADR1
LTOT = NSYM*NNSHL
call IDAFILE(LURED,IOPT,NNBSTRSH(:,:,ILOC),LTOT,IADR)
IOPT = 2
IADR = IADR1+NSYM*NNSHL
LSAV = CHO_ISUMELM(NNBSTRSH(:,:,ILOC),NSYM*NNSHL)
LTOT = LSAV
call IDAFILE(LURED,IOPT,INDRED(:,ILOC),LTOT,IADR)
if (LRSH .and. (IPASS == 1)) then
  IOPT = 2
  IADR = IADR1+NSYM*NNSHL+LSAV
  LTOT = LSAV
  call IDAFILE(LURED,IOPT,INDRSH,LTOT,IADR)
  IOPT = 2
  IADR = IADR1+NSYM*NNSHL+2*LSAV
  LTOT = NNSHL
  call IDAFILE(LURED,IOPT,ISP2F,LTOT,IADR)
end if

end subroutine CHO_GETRED
