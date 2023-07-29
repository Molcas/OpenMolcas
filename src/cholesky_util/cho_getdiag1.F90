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

subroutine CHO_GETDIAG1(DIAG,BUF,IBUF,LENBUF,NDUMP)
!
! Purpose: read diagonal in first reduced set.

use Cholesky, only: IndRed, IndRSh, INF_DIAG, IPRINT, mmBstRT, nnBstRT, nSym, RstDia
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: LENBUF, IBUF(4,LENBUF), NDUMP
real(kind=wp) :: Diag(*), BUF(LENBUF)
integer(kind=iwp) :: IOPT, IRED, IRS, ISYLST(8), ISYM, NSYLST
logical(kind=iwp), parameter :: LOCDBG = .false.

! Read diagonal from file.
! ------------------------

if (RSTDIA) then
  IOPT = 2
  call CHO_IODIAG(DIAG,IOPT)
else
  call FZERO(DIAG,NNBSTRT(1))
  call IZERO(INDRSH,NNBSTRT(1))
  call IZERO(INDRED,NNBSTRT(1))
  call CHO_RDDBUF(DIAG,BUF,IBUF,INDRSH,INDRED,LENBUF,MMBSTRT,NDUMP)
  call CHO_GADGOP(DIAG,NNBSTRT(1),'+')
  call CHO_GAIGOP(INDRSH,NNBSTRT(1),'+')
  call CHO_GAIGOP(INDRED,NNBSTRT(1),'+')
end if

! Copy info to current reduced set (IRED=2).
! Also set up IRED=3 (although it should be redundant).
! -----------------------------------------------------

do IRS=2,3
  call CHO_RSCOPY(1,IRS)
end do

! Print.
! ------

if (LOCDBG .or. (IPRINT >= INF_DIAG)) then
  do ISYM=1,NSYM
    ISYLST(ISYM) = ISYM
  end do
  NSYLST = NSYM
  IRED = 1
  call CHO_PRTDIA(DIAG,ISYLST,NSYLST,IRED)
  if (LOCDBG) then
    IRED = 2
    call CHO_PRTDIA(DIAG,ISYLST,NSYLST,IRED)
  end if
end if

end subroutine CHO_GETDIAG1
