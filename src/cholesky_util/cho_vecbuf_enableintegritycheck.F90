!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_VecBuf_EnableIntegrityCheck(irc)
!
! Thomas Bondo Pedersen, September 2012.
!
! Enable integrity check of buffer: allocate and store norm and sum
! of each vector in the buffer.

use ChoArr, only: nDimRS
use ChoSwp, only: InfVec
use ChoVecBuf, only: CHVBFI, CHVBUF, ip_CHVBFI_SYM, ip_CHVBUF_SYM, l_CHVBFI_SYM, nVec_in_Buf
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: irc
#include "cholesky.fh"
#include "choprint.fh"
real(kind=wp), external :: Cho_dSumElm, dDot_
integer(kind=iwp) :: ip, ipV, iSym, jRed, jVec, l_ChVBfI

! Set return code
irc = 0

! Not implemented in internal run mode (since the buffered vectors
! may change length and hence norm and sum, causing too much book
! keeping activity for a debug feature).
if (RUN_MODE /= RUN_EXTERNAL) return

! Return if no buffer is allocated
if (.not. allocated(CHVBUF)) return

! Return if already enabled
if (allocated(CHVBFI)) return

! Check that nDimRS is allocated
if (.not. allocated(nDimRS)) then
  irc = 1
  return
end if

! Allocate and store norm and sum of each vector in the buffer
l_ChVBfI = 0
do iSym=1,nSym
  l_ChVBfI_Sym(iSym) = 2*nVec_in_Buf(iSym)
  l_ChVBfI = l_ChVBfI+l_ChVBfI_Sym(iSym)
end do
if (l_ChVBfI > 0) then
  call mma_allocate(CHVBFI,l_ChVBfI,Label='CHVBFI')
  ip = 1
  do iSym=1,nSym
    ip_ChVBfI_Sym(iSym) = ip
    ip = ip+l_ChVBfI_Sym(iSym)
  end do
  do iSym=1,nSym
    ipV = ip_ChVBuf_Sym(iSym)
    ip = ip_ChvBfI_Sym(iSym)
    do jVec=1,nVec_in_Buf(iSym)
      jRed = InfVec(jVec,2,iSym)
      CHVBFI(ip) = sqrt(dDot_(nDimRS(iSym,jRed),CHVBUF(ipV),1,CHVBUF(ipV),1))
      CHVBFI(ip+1) = Cho_dSumElm(CHVBUF(ipV),nDimRS(iSym,jRed))
      ipV = ipV+nDimRS(iSym,jRed)
      ip = ip+2
    end do
  end do
  if (iPrint > 2) call Cho_VecBuf_PrtRef('@NABLE')
  write(LuPri,'(A)') 'Cholesky vector buffer integrity checks enabled'
else
  call iZero(l_ChVBfI_Sym,nSym)
  call iZero(ip_ChVBfI_Sym,nSym)
end if

end subroutine Cho_VecBuf_EnableIntegrityCheck
