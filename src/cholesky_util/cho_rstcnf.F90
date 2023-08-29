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

subroutine CHO_RSTCNF(NERR)
!
! Purpose: check restart configuration info.

use Cholesky, only: Cho_AdrVec, Damp, LuPri, ScDiag, Span, ThrCom, ThrDiag, ThrNeg, TooNeg, WarNeg, XCho_AdrVec, XDamp, XScDiag, &
                    XSpan, XThrCom, XThrDiag, XThrNeg, XTooNeg, XWarNeg
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: NERR
integer(kind=iwp) :: I, J
real(kind=wp) :: ERR
real(kind=wp), parameter :: ERRTOL = 1.0e-14_wp
character(len=*), parameter :: SWITCH(2) = [' ON','OFF']

NERR = 0

if (CHO_ADRVEC /= XCHO_ADRVEC) then
  write(LUPRI,'(A,I9,/,A,I9)') 'RESTART: addressing mode for vectors from restart file:',XCHO_ADRVEC, &
                               '         addressing mode for vectors from input       :',CHO_ADRVEC
  write(LUPRI,'(A,A)') '         Restart will fail - please specify correct address mode.'
  call CHO_QUIT('Cholesky restart failure in CHO_RSTCNF',105)
end if

ERR = abs(THRCOM-XTHRCOM)
if (ERR > ERRTOL) then
  write(LUPRI,'(A,D16.8,/,A,D16.8)') 'RESTART: decomposition threshold from restart file: ',XTHRCOM, &
                                     '         decomposition threshold from input       : ',THRCOM
  NERR = NERR+1
end if

ERR = abs(THRDIAG-XTHRDIAG)
if (ERR > ERRTOL) then
  write(LUPRI,'(A,D16.8,/,A,D16.8)') 'RESTART: init. diag. screening from restart file: ',XTHRDIAG, &
                                     '         init. diag. screening from input       : ',THRDIAG
  NERR = NERR+1
end if

ERR = abs(DAMP(1)-XDAMP(1))
if (ERR > ERRTOL) then
  write(LUPRI,'(A,D16.8,/,A,D16.8)') 'RESTART: 1st screening damping from restart file: ',XDAMP(1), &
                                     '         1st screening damping from input       : ',DAMP(1)
  NERR = NERR+1
end if

ERR = abs(DAMP(2)-XDAMP(2))
if (ERR > ERRTOL) then
  write(LUPRI,'(A,D16.8,/,A,D16.8)') 'RESTART: 2nd screening damping from restart file: ',XDAMP(2), &
                                     '         2nd screening damping from input       : ',DAMP(2)
  NERR = NERR+1
end if

if (SCDIAG .neqv. XSCDIAG) then
  if (XSCDIAG) then
    I = 1
  else
    I = 2
  end if
  if (SCDIAG) then
    J = 1
  else
    J = 2
  end if
  write(LUPRI,'(A,A,/,A,A)') 'RESTART: diag. screening from restart file: ',SWITCH(I), &
                             '         diag. screening from input       : ',SWITCH(J)
  NERR = NERR+1
end if

ERR = abs(THRNEG-XTHRNEG)
if (ERR > ERRTOL) then
  write(LUPRI,'(A,D16.8,/,A,D16.8)') 'RESTART: neg. diag. threshold from restart file: ',XTHRNEG, &
                                     '         neg. diag. threshold from input       : ',THRNEG
  NERR = NERR+1
end if

ERR = abs(WARNEG-XWARNEG)
if (ERR > ERRTOL) then
  write(LUPRI,'(A,D16.8,/,A,D16.8)') 'RESTART: neg. diag. warn thr. from restart file: ',XWARNEG, &
                                     '         neg. diag. warn thr. from input       : ',WARNEG
  NERR = NERR+1
end if

ERR = abs(TOONEG-XTOONEG)
if (ERR > ERRTOL) then
  write(LUPRI,'(A,D16.8,/,A,D16.8)') 'RESTART: too neg. diag. thr. from restart file: ',XTOONEG, &
                                     '         too neg. diag. thr. from input       : ',TOONEG
  NERR = NERR+1
end if

ERR = abs(SPAN-XSPAN)
if (ERR > ERRTOL) then
  write(LUPRI,'(A,D16.8,/,A,D16.8)') 'RESTART: span factor from restart file: ',XSPAN, &
                                     '         span factor from input       : ',SPAN
  NERR = NERR+1
end if

end subroutine CHO_RSTCNF
