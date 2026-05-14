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

subroutine CHO_VECRD(SCR,LSCR,JVEC1,IVEC2,ISYM,JNUM,IREDC,MUSED)
!
! Purpose: read as many vectors as fit into SCR array starting
!          at vector JVEC1 and reading at most until vector IVEC2.
!          On exit, JNUM is the number of vectors read.
!          On entry as well as exit, IREDC identifies the reduced
!          set stored in core (at position "3"; use -1 if none
!          or unkown). Vectors are taken from buffer, if possible.
!
! NOTE: if no vectors can be read, JNUM=0 and MUSED=0 are returned,
!       but execution is NOT stopped here!!!

use Cholesky, only: RUN_EXTERNAL, RUN_MODE
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LSCR, JVEC1, IVEC2, ISYM
real(kind=wp), intent(inout) :: SCR(LSCR)
integer(kind=iwp), intent(out) :: JNUM, MUSED
integer(kind=iwp), intent(inout) :: IREDC
integer(kind=iwp) :: JN, JV1, KS, LFT, MU
logical(kind=iwp) :: DOREAD

! Initialize.
! -----------

JNUM = 0
MUSED = 0
if (LSCR < 1) return

! Copy vectors from buffer (if possible).
! Only relevant for "external" runs (else the vectors are stored in
! current reduced set).
! -----------------------------------------------------------------

if (RUN_MODE == RUN_EXTERNAL) call CHO_VECBUF_RETRIEVE(SCR,LSCR,JVEC1,IVEC2,ISYM,JNUM,IREDC,MUSED)

! Read remaining vectors from disk.
! ---------------------------------

DOREAD = .true.
JV1 = JVEC1+JNUM
LFT = LSCR-MUSED
if ((IVEC2 >= JV1) .and. (LFT > 0)) then
  KS = MUSED+1
  JN = 0
  MU = 0
  call CHO_VECRD1(SCR(KS),LFT,JV1,IVEC2,ISYM,JN,IREDC,MU,DOREAD)
  JNUM = JNUM+JN
  MUSED = MUSED+MU
end if

end subroutine CHO_VECRD
