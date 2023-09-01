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

subroutine CHO_SUBTR(XINT,WRK,LWRK,ISYM)
!
! Purpose: driver for subtracting contributions from previous vectors
!          from the qualified integrals (in XINT).

use Cholesky, only: Cho_DiaChk, Cho_IOVec, LuPri, nnBstR, nnBstRT, nQual, NumCho, Tol_DiaChk
use Cholesky_procedures, only: Cho_VecBuf_Subtr
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: LWRK
real(kind=wp), intent(inout) :: XINT(*), WRK(LWRK)
integer(kind=iwp), intent(in) :: ISYM
integer(kind=iwp) :: KDIAG, KEND, NERR
real(kind=wp) :: TOL
logical(kind=iwp) :: FXDMEM
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_SUBTR'

! Return if nothing to do.
! ------------------------

if (NUMCHO(ISYM) < 1) then ! no prev. vectors.
  if (LOCDBG) write(LUPRI,*) SECNAM,': nothing done because NUMCHO = ',NUMCHO(ISYM),' (sym. ',ISYM,')'
  return
else if (NNBSTR(ISYM,2) < 1) then ! nothing to do (this sym.)
  if (LOCDBG) write(LUPRI,*) SECNAM,': nothing done because NNBSTR = ',NNBSTR(ISYM,2),' (sym. ',ISYM,')'
  return
else if (NQUAL(ISYM) < 1) then ! no qualifieds in this sym.
  if (LOCDBG) write(LUPRI,*) SECNAM,': nothing done because NQUAL  = ',NQUAL(ISYM),' (sym. ',ISYM,')'
  return
end if

! Debug: read original diagonal and check that these elements are
!        included in the integrals
! ---------------------------------------------------------------

if (CHO_DIACHK .or. LOCDBG) then
  KDIAG = 1
  KEND = KDIAG+NNBSTRT(1)
  LWRK = LWRK-KEND+1
  if (LWRK < 0) then
    write(LUPRI,*) SECNAM,': diagonal/integral check skipped due to insufficient memory'
  else
    TOL = TOL_DIACHK
    NERR = 0
    call CHO_CHKINTO(XINT,WRK(KDIAG),ISYM,NERR,TOL,.true.)
    if (NERR /= 0) then
      write(LUPRI,*) SECNAM,': ',NERR,' diagonal errors found!'
      write(LUPRI,*) '          #tests: ',NQUAL(ISYM)
      !write(LUPRI,*) '          Printing integrals:'
      !call CHO_OUTPUT(XINT,1,NNBSTR(ISYM,2),1,NQUAL(ISYM),NNBSTR(ISYM,2),NQUAL(ISYM),1,LUPRI)
      call CHO_QUIT('Diagonal errors in '//SECNAM,104)
    else
      write(LUPRI,*) SECNAM,': comparison of qual. integrals and original diagonal: no errors !'
    end if
  end if
end if

! Subtract contributions for vectors in buffer.
! (Returns immediately if nothing to do.)
! ---------------------------------------------

call CHO_VECBUF_SUBTR(XINT,WRK,LWRK,ISYM,.true.,.true.)

! Subtract contributions for vectors on disk.
! -------------------------------------------

if ((CHO_IOVEC == 3) .or. (CHO_IOVEC == 4)) then
  FXDMEM = CHO_IOVEC == 4
  call CHO_SUBTR1(XINT,WRK,LWRK,ISYM,FXDMEM)
else
  call CHO_SUBTR0(XINT,WRK,LWRK,ISYM)
end if

end subroutine CHO_SUBTR
