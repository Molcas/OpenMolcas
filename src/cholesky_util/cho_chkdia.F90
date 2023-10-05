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

subroutine CHO_CHKDIA(DIAG,ISYM,XM,YM,ZM,NNEGT,NNEG,NCONV)
!
! Purpose: 1) find min. in updated diagonal, XM (sym. ISYM only)
!          2) find max. in updated diagonal, YM (sym. ISYM only)
!          3) find abs. max. in updated diagonal, ZM (sym. ISYM only)
!          4) count #diagonals < 0.0, NNEGT
!          5) count #diagonals < THRNEG, NNEGT
!          6) count #screenable diagonals, NCONV
!
! -- also: a) zero negative diagonals < THRNEG
!          b) screen diagonal if requested (flag SCDIAG)
!          c) Keep track of most negative zeroed diagonal.

use Cholesky, only: Damp, DIAMNZ, IABMNZ, iiBstR, IndRed, LuPri, nnBstR, ScDiag, ThrCom, ThrNeg, TooNeg, WarNeg
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Diag(*)
integer(kind=iwp), intent(in) :: ISYM
real(kind=wp), intent(out) :: XM, YM, ZM
integer(kind=iwp), intent(out) :: NNEGT, NNEG, NCONV
integer(kind=iwp) :: IAB, JAB, JAB1, JAB2
real(kind=wp) :: TST
character(len=*), parameter :: SECNAM = 'CHO_CHKDIA'

! Initialization.
! ---------------

NNEGT = 0
NNEG = 0
NCONV = 0

if (NNBSTR(ISYM,2) > 0) then
  JAB1 = IIBSTR(ISYM,2)+1
  JAB2 = JAB1+NNBSTR(ISYM,2)-1
  XM = DIAG(INDRED(JAB1,2))
  YM = DIAG(INDRED(JAB1,2))
  ZM = abs(YM)
else
  XM = Zero
  YM = Zero
  ZM = Zero
  return
end if

! Find min. and max. diagonal and zero too negative diagonals.
! ------------------------------------------------------------

do JAB=JAB1,JAB2
  IAB = INDRED(JAB,2)  ! get address in first red. set
  XM = min(XM,DIAG(IAB))
  YM = max(YM,DIAG(IAB))
  if (DIAG(IAB) < Zero) then
    NNEGT = NNEGT+1
    if (DIAG(IAB) < THRNEG) then
      NNEG = NNEG+1
      if (DIAG(IAB) < TOONEG) then
        write(LUPRI,'(A,A,I12,1X,1P,D16.8)') SECNAM,': diagonal too negative: ',IAB,DIAG(IAB)
        write(LUPRI,'(A,A)') SECNAM,': shutting down Cholesky decomposition!'
        call CHO_QUIT('Diagonal too negative in '//SECNAM,104)
      end if
      if (DIAG(IAB) < WARNEG) write(LUPRI,'(A,A,I12,1X,1P,D16.8,A)') SECNAM,': Negative diagonal: ',IAB,DIAG(IAB),' (zeroed)'
      if (DIAG(IAB) < DIAMNZ) then
        DIAMNZ = DIAG(IAB)
        IABMNZ = IAB
      end if
      DIAG(IAB) = Zero
    end if
  end if
end do
ZM = max(abs(XM),abs(YM))

! Screen diagonal (if requested) and count the screenables as converged.
! NOTE: the screening is actually identical to setting up
!       reduced sets. However, doing the screening here will
!       force entries in later vectors of this integral pass
!       to have zero entries at screened diagonals.
! -----------------------------------------------------------

if (SCDIAG) then
  do JAB=JAB1,JAB2
    IAB = INDRED(JAB,2)  ! get address in first red. set
    TST = sqrt(abs(DIAG(IAB))*ZM)*DAMP(2)
    if (TST <= THRCOM) then
      NCONV = NCONV+1
      DIAG(IAB) = Zero
    end if
  end do
else
  do JAB=JAB1,JAB2
    IAB = INDRED(JAB,2)  ! get address in first red. set
    TST = sqrt(abs(DIAG(IAB))*ZM)*DAMP(2)
    if (TST <= THRCOM) NCONV = NCONV+1
  end do
end if

end subroutine CHO_CHKDIA
