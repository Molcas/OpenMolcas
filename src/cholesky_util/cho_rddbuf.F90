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

subroutine CHO_RDDBUF(DIAG,BUF,IBUF,INDRSH,INDRED,LENBUF,LMMBSTRT,NDUMP)
!
! Purpose: read diagonal from disk and set first reduced set indices.

use Cholesky, only: iiBstR, iiBstRsh, iSP2F, lBuf, LuPri, LuScr
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: Diag(*)
integer(kind=iwp), intent(in) :: LENBUF, LMMBSTRT, NDUMP
real(kind=wp), intent(out) :: BUF(LENBUF)
integer(kind=iwp), intent(out) :: IBUF(4,LENBUF), INDRSH(LMMBSTRT), INDRED(LMMBSTRT,3)
integer(kind=iwp) :: IAB, IDUMP, ISHLAB, ISYMAB, IUNIT, L, LENGTH
character(len=*), parameter :: SECNAM = 'CHO_RDDBUF'

if (LENBUF < LBUF) then
  write(LUPRI,'(//,1X,A,A)') SECNAM,': LENBUF >= LBUF required!'
  write(LUPRI,'(1X,A,I10)') 'LENBUF = ',LENBUF
  write(LUPRI,'(1X,A,I10,/)') 'LBUF   = ',LBUF
  call CHO_QUIT('Buffer error in '//SECNAM,102)
end if

IUNIT = LUSCR
LUSCR = -1
rewind(IUNIT)

do IDUMP=1,NDUMP
  call CHO_RDBUF(LENGTH,BUF,IBUF,LBUF,IUNIT)
  if (IDUMP == NDUMP) call CHO_CLOSE(IUNIT,'DELETE')
  do L=1,LENGTH
    if (IBUF(2,L) > 0) then
      ISHLAB = IBUF(1,L)
      ISYMAB = IBUF(3,L)
      IAB = IIBSTR(ISYMAB,1)+IIBSTRSH(ISYMAB,ISHLAB,1)+IBUF(2,L)
      DIAG(IAB) = BUF(L)
      INDRSH(IAB) = ISP2F(ISHLAB)
      INDRED(IAB,1) = IBUF(4,L)
    end if
  end do
end do

end subroutine CHO_RDDBUF
