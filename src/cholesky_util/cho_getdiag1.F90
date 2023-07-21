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

use ChoSwp, only: IndRSh, IndRed

implicit real*8(a-h,o-z)
real*8 Diag(*), BUF(LENBUF)
integer IBUF(4,LENBUF)
#include "cholesky.fh"
#include "choprint.fh"
character*12 SECNAM
parameter(SECNAM='CHO_GETDIAG1')
logical LOCDBG
parameter(LOCDBG=.false.)
integer ISYLST(8)
parameter(INFOD=INF_DIAG)
parameter(TINY=1.0D-14)

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

if (LOCDBG .or. (IPRINT >= INFOD)) then
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
